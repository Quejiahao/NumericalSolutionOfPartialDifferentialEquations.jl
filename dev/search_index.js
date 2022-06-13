var documenterSearchIndex = {"docs":
[{"location":"LaplaceEquation/#Laplace-Equation","page":"Laplace Equation","title":"Laplace Equation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NumericalSolutionOfPartialDifferentialEquations","category":"page"},{"location":"#NumericalSolutionOfPartialDifferentialEquations","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NumericalSolutionOfPartialDifferentialEquations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NumericalSolutionOfPartialDifferentialEquations]","category":"page"},{"location":"#Base.:*-Tuple{LinearAlgebra.UniformScaling{Bool}, AbstractVecOrMat}","page":"Home","title":"Base.:*","text":"Default operation overloads for I (identity matrix) are too slow, so I write these overloads.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.add_boundary_condition-Union{Tuple{T2}, Tuple{T1}, Tuple{Vector{T1}, Array{Vector{T2}, 1}}} where {T1<:Number, T2<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.add_boundary_condition","text":"add_boundary_condition(\n    f::Vector{T1},\n    boundary::Vector{Vector{T2}},\n) where {T1<:Number,T2<:Number}\n\nAdd Dirichlet boundary conditions to the right of difference equation.\n\nLeft:   boundary[1]\nBottom: boundary[2]\nRight:  boundary[3]\nTop:    boundary[4]\n\nThey are column vectors and the orientation is counterclockwise.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_boundary-Union{Tuple{T}, Tuple{Vector{T}, Int64}, Tuple{Vector{T}, Int64, Int64}} where T<:Function","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_boundary","text":"construct_boundary(\n    bound_func::Vector{T},\n    len::Int,\n    wid::Int = len,\n) where {T<:Function}\n\nUse bound_func to construct boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_element_stiffness_matrix_with_abs_det_A_e_2-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_element_stiffness_matrix_with_abs_det_A_e_2","text":"A^e = beginpmatrix\n    x_1^1  x_1^2  x_1^3 \n    x_2^1  x_2^2  x_2^3\nendpmatrix\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_grid-Union{Tuple{T}, Tuple{T, Int64}, Tuple{T, Int64, Int64}} where T<:Function","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_grid","text":"construct_grid(f::T, len::Int, wid::Int = len) where {T<:Function}\n\nUse f to construct grid.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_laplacian-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_laplacian","text":"construct_laplacian(; size::Int = 32, dim::Int = 2)\n\nConstruct Laplacian by this method:\n\nboldsymbol L_mathrmdim = sum_i = 1^m left( bigotimes_j = 1^i - 1 boldsymbol I_n right) otimes boldsymbol L_1 otimes left(bigotimes_j = i + 1^m boldsymbol I_nright)\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.homework1-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.homework1","text":"Need ~500 GiB Memory!!!\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.plot_1d_heat_equ_solution-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{Vector{T1}, Int64, Int64}} where {T1<:Number, T2<:Number, T3<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.plot_1d_heat_equ_solution","text":"Also used by convection diffusion equation\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.solve_heat_equation-Union{Tuple{T4}, Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{T1, Vector{T2}}} where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.solve_heat_equation","text":"solve_heat_equation(\n    thermal_diffusivity::T1,\n    initial::Vector{T2};\n    solver::Function = \\,\n    total_time::T3 = 1.0,\n    time_step_num::Int = 1024,\n    theta::T4 = 0.5,\n    kw...,\n) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}\n\nSolve Heat equation with Dirichlet boundary conditions (zero boundary and given initial).\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.solve_heat_equation_explicit-Union{Tuple{T2}, Tuple{T1}, Tuple{T1, T2}} where {T1<:Number, T2<:Union{Function, Vector{<:Number}}}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.solve_heat_equation_explicit","text":"Three Scheme\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.solve_poissions_equation-Union{Tuple{T2}, Tuple{T1}, Tuple{Vector{T1}, Array{Vector{T2}, 1}}} where {T1<:Number, T2<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.solve_poissions_equation","text":"solve_poissions_equation(\n    f::Vector{T1},\n    boundary::Vector{Vector{T2}},\n) where {T1<:Number,T2<:Number}\n\nSolve Poisson's equation with Dirichlet boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.test_solve_heat_equation-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.test_solve_heat_equation","text":"(x t) in (0 1)^2\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_known-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_known","text":"# 已知精确解\n# 共轭梯度?LU分解?\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_unknown-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_unknown","text":"# 无精确解 log-log plot\n\n\n\n\n\n","category":"method"}]
}
