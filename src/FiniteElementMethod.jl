function construct_element_stiffness_matrix()
    
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
        return solve_poissions_equation(
            f,
            boundary;
            solver = solver,
            kw...,
        )
    elseif scheme === :FEM
        return solve_poissions_equation_FEM(
            f,
            boundary;
            solver = solver,
            kw...,
        )
    end
end
