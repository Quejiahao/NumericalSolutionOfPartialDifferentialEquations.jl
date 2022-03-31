module NumericalSolutionOfPartialDifferentialEquations

export construct_laplacian,
    conjugate_gradient_inverse_divide,
    solve_poissions_equation,
    test_solve_poissions_equation_known,
    test_solve_poissions_equation_unknown

if isdefined(ENV, :NOMKL) && ENV["NOMKL"] === "1"
    @info "Do not use MKL, using default BLAS"
else
    try
        @eval Main begin
            using MKL
        end
    catch
        @info "MKL not found, using default BLAS"
    end
end

using LinearAlgebra
using SparseArrays
using Optim

import Base: *, \

"""
Default operation overloads for `I` (identity matrix) are too slow, so I write these overloads.
"""
*(J::UniformScaling{Bool}, A::AbstractVecOrMat) = J.λ ? A : zero(A)
*(A::AbstractVecOrMat, J::UniformScaling{Bool}) = J.λ ? A : zero(A)
\(J::UniformScaling{Bool}, A::AbstractVecOrMat) = J.λ ? A : throw(SingularException(1))

function con_tri_diag(
    size::Int,
    main_diag::T,
    first_above::T,
    first_below::T = first_above,
) where {T<:Number}
    Is = [1:size; 1:(size-1); 2:size]
    Js = [1:size; 2:size; 1:(size-1)]
    Vs = [
        fill(main_diag, size)
        fill(first_above, size - 1)
        fill(first_below, size - 1)
    ]
    return sparse(Is, Js, Vs)
end

@doc raw"""
    construct_laplacian(; size::Int = 32, dim::Int = 2)

Construct Laplacian by this method:

```math
\boldsymbol L_{\mathrm{dim}} = \sum_{i = 1}^{m} \left( \bigotimes_{j = 1}^{i - 1} \boldsymbol I_n \right) \otimes \boldsymbol L_1 \otimes \left(\bigotimes_{j = i + 1}^{m} \boldsymbol I_n\right)
```
"""
function construct_laplacian(; size::Int = 32, dim::Int = 2)
    lap_1d = con_tri_diag(size, -2.0, 1.0)

    if dim === 1
        return lap_1d
    end

    lap_dim = kron(lap_1d, sparse(I, size * (dim - 1), size * (dim - 1)))
    for i = 2:(dim-1)
        lap_dim += kron(
            sparse(I, size * (i - 1), size * (i - 1)),
            lap_1d,
            sparse(I, size * (dim - i), size * (dim - i)),
        )
    end
    lap_dim += kron(sparse(I, size * (dim - 1), size * (dim - 1)), lap_1d)

    return lap_dim
end

function plot_2d_solution(U::Vector{T}, len::Int, wid::Int = len; kw...) where {T<:Number}
    xs = range(0.0, 1.0, wid + 2)[2:end-1]
    ys = range(1.0, 0.0, len + 2)[2:end-1]
    x_grid = [x for x in xs for y in ys]
    y_grid = [y for x in xs for y in ys]
    @eval (@__MODULE__) begin
        using Plots
        theme(:orange)
        display(plot($x_grid, $y_grid, $U, st = :surface; $kw...))
    end
end

function conjugate_gradient_inverse_divide(A, b)
    f = x -> (x' * A / 2 - b') * x
    g! = (G, x) -> (G[:] = (A * x - b))
    return optimize(f, g!, zeros(length(b)), ConjugateGradient()).minimizer
end

include("PoissonsEquation.jl")

include("HeatEquation.jl")

end
