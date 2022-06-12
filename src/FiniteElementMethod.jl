function construct_triangle_grid(f::Function, len::Int, wid::Int = len)
    grid_point_x = repeat([range(0.0, 1.0, wid + 2);], wid + 2)
    grid_point_y = repeat([range(0.0, 1.0, len + 2);], inner = len + 2)
    fval = f.(grid_point_x, grid_point_y)
    # _grid_point_index = div.(4 : (2 * (wid + 1) * (len + 1) + 3), 2)
    _grid_point_index =
        repeat(deleteat!([1:((wid+2)*(len+1));], 1:(wid+2):((wid+2)*(len+1))), inner = 2)
    T = [_grid_point_index _grid_point_index .+ (wid + 1) _grid_point_index .- 1]
    T[2:2:end, 3] .+= (wid + 3)
    return Dict(
        :fval => fval,
        :T => T,
    )
end

@doc raw"""
$$A^e = \begin{pmatrix}
    x_1^1 & x_1^2 & x_1^3 \\
    x_2^1 & x_2^2 & x_2^3
\end{pmatrix}$$
"""
function construct_element_stiffness_matrix(
    Ae::Array{T,2};
    isnot_boundary = [true; true; true],
    A_e::Array{T,2} = Ae[:, 2:3] .- Ae[:, 1],
) where {T<:Number}
    nabla_lambdas_hat_x = [
        -1.0 -1.0
        1.0 0.0
        0.0 1.0
    ] .* isnot_boundary

    det_Ae_2 = det(A_e) * 2.0
    A_e_star = [
        A_e[2, 2] -A_e[1, 2]
        -A_e[2, 1] A_e[1, 1]
    ]
    nabla_lambdas_hat_x_mul_A_e_star = nabla_lambdas_hat_x * A_e_star
    return Symmetric(
        nabla_lambdas_hat_x_mul_A_e_star * transpose(nabla_lambdas_hat_x_mul_A_e_star) ./
        det_Ae_2,
    )
end

function construct_stiffness_matrix()

end

function _fval_element_load_vector(
    Ae::Array{T,2},
    f::Function,
    isnot_boundary = [true; true; true],
    integral_method::Symbol = :piecewise_linear,
) where {T<:Number}
    if integral_method == :piecewise_linear
        return [f(Ae[:, i]...) for i = 1:size(Ae, 2)]
    else
        error("Unknown integral method: ", integral_method)
    end
end

function construct_element_load_vector(
    Ae::Array{T,2},
    f::Function;
    isnot_boundary = [true; true; true],
    integral_method::Symbol = :piecewise_linear,
    A_e::Array{T,2} = Ae[:, 2:3] .- Ae[:, 1],
    fval::Array{T,2} = _fval_element_load_vector(Ae, f, isnot_boundary, integral_method),
) where {T<:Number}
    if integral_method == :piecewise_linear
        return det(Ae) ./ 24 * [
            2 1 1
            1 2 1
            1 1 2
        ] * fval
    else
        error("Unknown integral method: ", integral_method)
    end
end

function construct_load_vector()

end

function solve_poissions_equation_FEM(
    f::Function,
    # boundary::Vector{Vector{T2}}  zero boundary
    ;
    solver::Function = \,
    kw...,
) where {T1<:Number,T2<:Number}
    # calc A_e
    sti_mat = construct_stiffness_matrix()
    loa_mat = construct_load_vector()
    return solver(sti_mat, loa_mat)
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
        return solve_poissions_equation(f, boundary; solver = solver, kw...)
    elseif scheme === :FEM
        return solve_poissions_equation_FEM(fun; solver = solver, kw...)
    end
end
