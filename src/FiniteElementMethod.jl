Base.@kwdef mutable struct TriangleGrid{T1,T2,T3<:Number}
    Tes_index::Array{Int,2}
    A::Array{T1,2}
    fval::Array{T2,1}
    isnot_boundary::Array{Bool,1}
    boundary::SparseArrays.SparseVector{T3,Int64}
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

function _construct_boundary(
    boundary::Vector{Vector{T}},
    isnot_boundary::Array{Bool,1},
    len::Int,
    wid::Int = len,
) where {T<:Number}
    bound_index = findall(x -> x == false, isnot_boundary)
    _boundary = zeros(T, length(bound_index))
    _boundary[(wid+1+2*len):(-2):(wid+3)] = boundary[1][2:(end-1)]
    _boundary[1:(wid+2)] = boundary[2]
    _boundary[(wid+4):2:(wid+2+2*len)] = boundary[3][2:(end-1)]
    _boundary[end:(-1):(wid+3+2*len)] = boundary[4]
    __boundary = sparsevec(bound_index, _boundary)
    return __boundary
end

function construct_triangle_grid(
    f::Function,
    boundary::Vector{Vector{T}},
    len::Int,
    wid::Int = len,
) where {T<:Number}
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
    __boundary = _construct_boundary(boundary, isnot_boundary, len, wid)
    return TriangleGrid(
        fval = fval,
        Tes_index = Tes_index,
        A = [
            transpose(grid_point_x)
            transpose(grid_point_y)
        ],
        isnot_boundary = isnot_boundary,
        boundary = __boundary,
        f = f,
    )
end

@doc raw"""
Boundary conditions integration.
"""
function bound_int()

end

@doc raw"""
$$A^e = \begin{pmatrix}
    x_1^1 & x_1^2 & x_1^3 \\
    x_2^1 & x_2^2 & x_2^3
\end{pmatrix}$$
"""
function construct_element_stiffness_matrix_with_abs_det_A_e_2_g_lambdae(
    Ae::Array{T1,2};
    isnot_boundarye = [true; true; true],
    A_e::Array{T1,2} = Ae[:, 2:3] .- Ae[:, 1],
    boundarye::SparseArrays.SparseVector{T2,Int64} = [0.0; 0.0; 0.0],
) where {T1,T2<:Number}
    nabla_lambdas_hat_x = [
        -1.0 -1.0
        1.0 0.0
        0.0 1.0
    ]

    abs_det_A_e_2 = abs(det(A_e)) * 2.0
    A_e_star = [
        A_e[2, 2] -A_e[1, 2]
        -A_e[2, 1] A_e[1, 1]
    ]
    nabla_lambdas_hat_x_mul_A_e_star = nabla_lambdas_hat_x * A_e_star
    Ke =
        nabla_lambdas_hat_x_mul_A_e_star * transpose(nabla_lambdas_hat_x_mul_A_e_star) ./
        abs_det_A_e_2
    g_lambdae = Ke * boundarye
    Ke .*= isnot_boundarye .* transpose(isnot_boundarye)
    return Ke, abs_det_A_e_2, g_lambdae
end

function construct_element_stiffness_matrix(Ae::Array{T,2}; kw...) where {T<:Number}
    return construct_element_stiffness_matrix_with_abs_det_A_e_2_g_lambdae(Ae; kw...)[1]
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
        f_lambdae =
            abs_det_A_e_2 ./ 48 .* [
                2 1 1
                1 2 1
                1 1 2
            ] .* isnot_boundarye .* transpose(isnot_boundarye) * fvale
        return f_lambdae
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
        boundarye = triangle_grid.boundary[Te_index]
        Ke, abs_det_A_e_2, g_lambdae =
            construct_element_stiffness_matrix_with_abs_det_A_e_2_g_lambdae(
                Ae;
                isnot_boundarye = isnot_boundarye,
                boundarye = boundarye,
            )
        fe = construct_element_load_vector(
            Ae,
            triangle_grid.f;
            isnot_boundarye = isnot_boundarye,
            abs_det_A_e_2 = abs_det_A_e_2,
            fvale = fvale,
            kw...,
        )
        fe -= g_lambdae

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
    boundary::Vector{Vector{T}};
    solver::Function = \,
    len::Int = 31,
    wid::Int = len,
    kw...,
) where {T<:Number}
    triangle_grid = construct_triangle_grid(f, boundary, len, wid)
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
        return solve_poissions_equation_FEM(fun, boundary; solver = solver, kw...)
    else
        error("Unknown scheme: ", scheme)
    end
end

function u_alpha(x, y; alpha = 0.05)
    r = sqrt(x .^ 2 + y .^ 2)
    alpha_r = alpha + r
    u = sin(1 ./ alpha_r)
    return u
end

function f_alpha(x, y; alpha = 0.05)
    r = sqrt(x .^ 2 + y .^ 2)
    alpha_r = alpha + r
    cos_ar = cos(1 ./ alpha_r)
    f =
        1 ./ (alpha_r) .* (
            -cos_ar ./ r + (2 .* cos_ar) ./ (alpha_r .^ 2) -
            sin(1 ./ alpha_r) ./ (alpha_r .^ 3)
        )
    return f
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
    is_norm = true,
    kw...,
)
    min_log2_size = 2
    err = zeros(max_log2_size - min_log2_size + 1)
    Ts = zeros(max_log2_size - min_log2_size + 1)
    for i = min_log2_size:max_log2_size
        t = time()
        err[i-min_log2_size+1] = test_solve_poissions_equation_known(;
            scheme = scheme,
            size = 2^i - 1,
            bound_func = [u, u, u, u],
            f = f,
            u = u,
            norm_p = norm_p,
            norm_k = norm_k,
            is_norm = is_norm,
            kw...,
        )
        Ts[i-min_log2_size+1] = time() - t
    end
    h = string.("2^{", string.([-min_log2_size:-1:-max_log2_size;]), "}")
    log2_err = log2.(err)
    order = [
        0
        log2_err[1:(end-1)] - log2_err[2:end]
    ]
    return [h err order Ts]
end

function homework5(; norm_p = 2, norm_k = 1, is_norm = false, kw...)
    return homework4(;
        f = f_alpha,
        u = u_alpha,
        norm_p = norm_p,
        norm_k = norm_k,
        is_norm = is_norm,
        kw...,
    )
end

function construct_1d_stiffness_matrix(N::Int)
    pi_3N = pi / (3 * N)
    K = con_tri_diag(N + 1, 2 * pi_3N, pi / (6 * N))
    K[0, 0] = K[N+1, N+1] = pi_3N
    return K
end

function construct_1d_load_vector(N)

end
