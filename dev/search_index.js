var documenterSearchIndex = {"docs":
[{"location":"LaplaceEquation/#Laplace-Equation","page":"Laplace Equation","title":"Laplace Equation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NumericalSolutionOfPartialDifferentialEquations","category":"page"},{"location":"#NumericalSolutionOfPartialDifferentialEquations","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NumericalSolutionOfPartialDifferentialEquations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NumericalSolutionOfPartialDifferentialEquations]","category":"page"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.add_boundary_condition-Union{Tuple{T2}, Tuple{T1}, Tuple{Vector{T1}, Array{Vector{T2}, 1}}} where {T1<:Number, T2<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.add_boundary_condition","text":"add_boundary_condition(\n    f::Vector{T1},\n    boundary::Vector{Vector{T2}},\n) where {T1<:Number,T2<:Number}\n\nAdd Dirichlet boundary conditions to the right of difference equation.\n\nLeft:   boundary[1]\nBottom: boundary[2]\nRight:  boundary[3]\nTop:    boundary[4]\n\nThey are column vectors and the orientation is counterclockwise.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_boundary-Union{Tuple{T}, Tuple{Vector{T}, Int64}, Tuple{Vector{T}, Int64, Int64}} where T<:Function","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_boundary","text":"construct_boundary(\n    bound_func::Vector{T},\n    len::Int,\n    wid::Int = len,\n) where {T<:Function}\n\nUse bound_func to construct boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_grid-Union{Tuple{T}, Tuple{T, Int64}, Tuple{T, Int64, Int64}} where T<:Function","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_grid","text":"construct_grid(f::T, len::Int, wid::Int = len) where {T<:Function}\n\nUse f to construct grid.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.construct_laplacian-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.construct_laplacian","text":"construct_laplacian(; size::Int = 32, dim::Int = 2)\n\nConstruct Laplacian by this method:\n\nboldsymbol L_mathrmdim = sum_i = 1^m left( bigotimes_j = 1^i - 1 boldsymbol I_n right) otimes boldsymbol L_1 otimes left(bigotimes_j = i + 1^m boldsymbol I_nright)\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.solve_poissions_equation-Union{Tuple{T2}, Tuple{T1}, Tuple{Vector{T1}, Array{Vector{T2}, 1}}} where {T1<:Number, T2<:Number}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.solve_poissions_equation","text":"solve_poissions_equation(\n    f::Vector{T1},\n    boundary::Vector{Vector{T2}},\n) where {T1<:Number,T2<:Number}\n\nSolve Poisson's equation with Dirichlet boundary conditions.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_known-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_known","text":"# 已知精确解\n# 共轭梯度?LU分解?\n\n\n\n\n\n","category":"method"},{"location":"#NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_unknown-Tuple{}","page":"Home","title":"NumericalSolutionOfPartialDifferentialEquations.test_solve_poissions_equation_unknown","text":"# 无精确解 log-log plot\n\n\n\n\n\n","category":"method"}]
}
