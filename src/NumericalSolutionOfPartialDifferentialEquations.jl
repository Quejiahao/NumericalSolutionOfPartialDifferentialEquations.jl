module NumericalSolutionOfPartialDifferentialEquations

export construct_laplacian

try
    @eval Main begin
        using MKL
    end
catch
    @info "MKL not found, using default BLAS"
end

using LinearAlgebra
using SparseArrays

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

end
