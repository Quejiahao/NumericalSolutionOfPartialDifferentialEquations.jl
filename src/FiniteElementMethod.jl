Base.@kwdef mutable struct TriangleGrid{T<:Number}
    Tes_index::Array{Int,2}
    A::Array{T,2}
    fval::Array{T,1}
    isnot_boundary::Array{Bool,1}
    f::Function
    # len::Int
    # wid::Int = len
end

function construct_grid_point_FEM(len::Int, wid::Int = len)
    return repeat([range(0.0, 1.0, wid + 2);], wid + 2),
    repeat([range(0.0, 1.0, len + 2);], inner = len + 2)
end

function construct_grid_FEM(f::Function, len::Int, wid::Int = len)
    grid_point_x, grid_point_y = construct_grid_point_FEM(len, wid)
    return f.(grid_point_x, grid_point_y)[_isnot_boundary(len, wid)]
end

function _isnot_boundary(len::Int, wid::Int = len)
    return [zeros(Bool, wid + 2) [
        zeros(Bool, 1, len)
        ones(Bool, wid, len)
        zeros(Bool, 1, len)
    ] zeros(Bool, wid + 2)][:]
end

function construct_triangle_grid(f::Function, len::Int, wid::Int = len)
    grid_point_x, grid_point_y = construct_grid_point_FEM(len, wid)
    fval = f.(grid_point_x, grid_point_y)
    # _grid_point_index = div.(4 : (2 * (wid + 1) * (len + 1) + 3), 2)
    _grid_point_index =
        transpose(deleteat!([1:((wid+2)*(len+1));], 1:(wid+2):((wid+2)*(len+1))))
    Tes_index = [
        _grid_point_index.-1 _grid_point_index
        _grid_point_index _grid_point_index.+(wid+1)
        _grid_point_index.+(wid+1) _grid_point_index.+(wid+2)
    ]
    isnot_boundary = _isnot_boundary(len, wid)
    return TriangleGrid(
        fval = fval,
        Tes_index = Tes_index,
        A = [
            transpose(grid_point_x)
            transpose(grid_point_y)
        ],
        isnot_boundary = isnot_boundary,
        f = f,
    )
end

@doc raw"""
$$A^e = \begin{pmatrix}
    x_1^1 & x_1^2 & x_1^3 \\
    x_2^1 & x_2^2 & x_2^3
\end{pmatrix}$$
"""
function construct_element_stiffness_matrix_with_abs_det_A_e_2(
    Ae::Array{T,2};
    isnot_boundarye = [true; true; true],
    A_e::Array{T,2} = Ae[:, 2:3] .- Ae[:, 1],
) where {T<:Number}
    nabla_lambdas_hat_x = [
        -1.0 -1.0
        1.0 0.0
        0.0 1.0
    ] .* isnot_boundarye

    abs_det_A_e_2 = abs(det(A_e)) * 2.0
    A_e_star = [
        A_e[2, 2] -A_e[1, 2]
        -A_e[2, 1] A_e[1, 1]
    ]
    nabla_lambdas_hat_x_mul_A_e_star = nabla_lambdas_hat_x * A_e_star
    Ke =
        nabla_lambdas_hat_x_mul_A_e_star * transpose(nabla_lambdas_hat_x_mul_A_e_star) ./
        abs_det_A_e_2
    return Ke, abs_det_A_e_2
end

function construct_element_stiffness_matrix(Ae::Array{T,2}; kw...) where {T<:Number}
    return construct_element_stiffness_matrix_with_abs_det_A_e_2(Ae; kw...)[1]
end

function _fval_element_load_vector(
    Ae::Array{T,2},
    f::Function,
    integral_method::Symbol = :piecewise_linear,
) where {T<:Number}
    if integral_method === :piecewise_linear
        return [f(Ae[:, i]...) for i = 1:size(Ae, 2)]
    else
        error("Unknown integral method: ", integral_method)
    end
end

function construct_element_load_vector(
    Ae::Array{T,2},
    f::Function;
    isnot_boundarye = [true; true; true],
    integral_method::Symbol = :piecewise_linear,
    abs_det_A_e_2::T = abs(det(Ae[:, 2:3] .- Ae[:, 1])) * 2.0,
    fvale::Array{T,1} = _fval_element_load_vector(Ae, f, integral_method),
) where {T<:Number}
    if integral_method === :piecewise_linear
        return abs_det_A_e_2 ./ 48 .* [
            2 1 1
            1 2 1
            1 1 2
        ] .* isnot_boundarye .* transpose(isnot_boundarye) * fvale
    elseif integral_method === :HCubature    #   I do not know what is this package algorithm
        return abs_det_A_e_2 ./ 2 .* hcubature(
            x ->
                f((1 - x[2]) * x[1], x[2]) * (1 - x[2]) .* [
                    1 - x[1] - x[2]
                    x[1]
                    x[2]
                ] .* isnot_boundarye,
            [0.0, 0.0],
            [1.0, 1.0],
        )[1]
    else
        error("Unknown integral method: ", integral_method)
    end
end

function construct_stiffness_matrix_and_load_vector(triangle_grid::TriangleGrid; kw...)
    _e_index_i = repeat(1:3, 3) #[1, 1, 2, 1, 2, 3]
    _e_index_j = repeat(1:3, inner = 3) #[1, 2, 2, 3, 3, 3]
    Tes_num = size(triangle_grid.Tes_index, 2)
    Vs_type = typeof(triangle_grid.fval[1])
    Is = zeros(Int, length(_e_index_i), Tes_num)
    Js = zeros(Int, length(_e_index_j), Tes_num)
    K_Vs = zeros(Vs_type, length(_e_index_i), Tes_num)
    f_Vs = zeros(Vs_type, 3, Tes_num)

    for e = 1:Tes_num
        Te_index = triangle_grid.Tes_index[:, e]
        Ae = triangle_grid.A[:, Te_index]
        fvale = triangle_grid.fval[Te_index]
        isnot_boundarye = triangle_grid.isnot_boundary[Te_index]
        Ke, abs_det_A_e_2 = construct_element_stiffness_matrix_with_abs_det_A_e_2(
            Ae;
            isnot_boundarye = isnot_boundarye,
        )
        fe = construct_element_load_vector(
            Ae,
            triangle_grid.f;
            isnot_boundarye = isnot_boundarye,
            abs_det_A_e_2 = abs_det_A_e_2,
            fvale = fvale,
            kw...,
        )

        Is[:, e] = Te_index[_e_index_i]
        Js[:, e] = Te_index[_e_index_j]
        K_Vs[:, e] = Ke[:] # [Ke[_e_index_i[k], _e_index_j[k]] for k = 1:length(_e_index_i)]
        f_Vs[:, e] = fe
    end

    grid_point_num = size(triangle_grid.A, 2)
    K = dropzeros(sparse(Is[:], Js[:], K_Vs[:], grid_point_num, grid_point_num))
    f = Array(sparsevec(triangle_grid.Tes_index[:], f_Vs[:], grid_point_num))

    return K, f
end

function delete_zero_boundary(K, f, triangle_grid::TriangleGrid)
    return K[triangle_grid.isnot_boundary, triangle_grid.isnot_boundary],
    f[triangle_grid.isnot_boundary]
end

function solve_poissions_equation_FEM(
    f::Function,
    # boundary::Vector{Vector{T2}}  zero boundary
    ;
    solver::Function = \,
    len::Int = 31,
    wid::Int = len,
    kw...,
)
    triangle_grid = construct_triangle_grid(f, len, wid)
    K, load_vector = construct_stiffness_matrix_and_load_vector(triangle_grid; kw...)
    K, load_vector = delete_zero_boundary(K, load_vector, triangle_grid)
    return solver(K, load_vector)
end

function solve_poissions_equation(
    f::Vector{T1},
    boundary::Vector{Vector{T2}},
    scheme::Symbol;
    fun::Function = identity,
    solver::Function = \,
    kw...,
) where {T1<:Number,T2<:Number}
    if scheme === :FDM
        @info "Use FDM scheme"
        return solve_poissions_equation(f, boundary; solver = solver, kw...)
    elseif scheme === :FEM
        # only zeros boundary for FEM, may add boundary later
        @info "Use FEM scheme"
        return solve_poissions_equation_FEM(fun; solver = solver, kw...)
    else
        error("Unknown scheme: ", scheme)
    end
end

function homework4(;
    max_log2_size = 11,
    scheme = :FEM,
    f = (x, y) -> -2.0 * (x^2 + y^2 - x - y),
    u = (x, y) -> x * (x - 1) * y * (y - 1),
    # f = (x, y) -> 2 * pi^2 * sin(pi * x) * sin(pi * y),
    # u = (x, y) -> sin(pi * x) * sin(pi * y),
    norm_p = Inf,
    norm_k = 0,
    kw...,
)
    err = zeros(max_log2_size)
    for i = 1:max_log2_size
        err[i] = test_solve_poissions_equation_known(;
            scheme = scheme,
            size = 2^i - 1,
            bound_func = [(x, y) -> 0, (x, y) -> 0, (x, y) -> 0, (x, y) -> 0],  # this argument is currently useless for FEM, only zero boundary is supported
            f = f,
            u = u,
            norm_p = norm_p,
            norm_k = norm_k,
            kw...,
        )
    end
    h = string.("2^{", string.([-1:-1:-11;]), "}")
    log2_err = log2.(err)
    order = [
        0
        log2_err[1:(end-1)] - log2_err[2:end]
    ]
    return [h err order]
end
