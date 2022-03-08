module NumericalSolutionOfPartialDifferentialEquations

export construct_laplacian, solve_poissions_equation, test_solve_poissions_equation

try
    @eval Main begin
        using MKL
    end
catch
    @info "MKL not found, using default BLAS"
end

using LinearAlgebra
using SparseArrays

@doc raw"""
    construct_laplacian(; size::Int = 32, dim::Int = 2)

Construct Laplacian by this method:

```math
\boldsymbol L_{\mathrm{dim}} = \sum_{i = 1}^{m} \left( \bigotimes_{j = 1}^{i - 1} \boldsymbol I_n \right) \otimes \boldsymbol L_1 \otimes \left(\bigotimes_{j = i + 1}^{m} \boldsymbol I_n\right)
```
"""
function construct_laplacian(; size::Int = 32, dim::Int = 2)
    Is = [1:size; 1:(size-1); 2:size]
    Js = [1:size; 2:size; 1:(size-1)]
    Vs = [fill(-2.0, size); fill(1.0, 2 * size - 2)]
    lap_1d = sparse(Is, Js, Vs)

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

include("PoissonsEquation.jl")

end
