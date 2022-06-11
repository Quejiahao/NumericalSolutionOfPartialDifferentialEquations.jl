@doc raw"""
$$A^e = \begin{pmatrix}
    x_1^1 & x_1^2 & x_1^3 \\
    x_2^1 & x_2^2 & x_2^3
\end{pmatrix}$$
"""
function construct_element_stiffness_matrix(
    Ae::Array{T,2},
    isnot_boundary = [true; true; true],
) where {T<:Number}
    nabla_lambdas_hat_x = [
        -1.0 -1.0
        1.0 0.0
        0.0 1.0
    ] .* isnot_boundary

    A_e = Ae[:, 2:3] .- Ae[:, 1]
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

function construct_element_load_vector()

end

function construct_load_vector()

end

function solve_poissions_equation_FEM(
    f::Vector{T1},
    boundary::Vector{Vector{T2}};
    solver::Function = \,
    kw...,
) where {T1<:Number,T2<:Number}
    sti_mat = construct_stiffness_matrix()
    loa_mat = construct_load_vector()
    return solver(sti_mat, loa_mat)
end

function solve_poissions_equation(
    f::Vector{T1},
    boundary::Vector{Vector{T2}},
    scheme::Symbol;
    solver::Function = \,
    kw...,
) where {T1<:Number,T2<:Number}
    if scheme === :FDM
        return solve_poissions_equation(f, boundary; solver = solver, kw...)
    elseif scheme === :FEM
        return solve_poissions_equation_FEM(f, boundary; solver = solver, kw...)
    end
end
