module NumericalSolutionOfPartialDifferentialEquations

export construct_laplacian,
    conjugate_gradient_inverse_divide,
    solve_poissions_equation,
    solve_poissions_equation_FEM,
    test_solve_poissions_equation_known,
    test_solve_poissions_equation_unknown,
    test_solve_heat_equation,
    test_solve_poissions_equation_FEM,
    homework1,
    homework2,
    homework3,
    homework4,
    homework5,
    homework6

if haskey(ENV, "NOMKL") && ENV["NOMKL"] === "1"
    @info "Do not use MKL, using default BLAS"
else
    try
        @eval Main begin
            # using MKL
        end
    catch
        @info "MKL not found, using default BLAS"
    end
end

using LinearAlgebra
using SparseArrays
using Optim
using HCubature

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
    Is, Js, Vs = Int[], Int[], T[]
    if main_diag !== zero(main_diag)
        Is = [Is; 1:size]
        Js = [Js; 1:size]
        Vs = [Vs; fill(main_diag, size)]
    end
    if first_above !== zero(first_above)
        Is = [Is; 1:(size-1)]
        Js = [Js; 2:size]
        Vs = [Vs; fill(first_above, size - 1)]
    end
    if first_below !== zero(first_below)
        Is = [Is; 2:size]
        Js = [Js; 1:(size-1)]
        Vs = [Vs; fill(first_below, size - 1)]
    end
    if length(Is) === 0
        return spzeros(size, size)
    end
    return sparse(Is, Js, Vs, size, size)
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

function construct_initial(
    initial::T1;
    space_step_num::Int = 31,
    total_space::T2 = 1.0,
    boundary::Bool = false,
    kw...,
) where {T1<:Function,T2<:Number}
    xs = range(0.0, total_space, space_step_num + 2)
    if !boundary
        xs = xs[2:end-1]
    end
    return initial.(xs)
end

function plot_2d_solution(
    U::Vector{T},
    len::Int,
    wid::Int = len;
    scheme = :FDM,
    kw...,
) where {T<:Number}
    xs = range(0.0, 1.0, wid + 2)[2:end-1]
    if scheme === :FEM
        ys = range(0.0, 1.0, len + 2)[2:end-1]
    else
        ys = range(1.0, 0.0, len + 2)[2:end-1]
    end
        x_grid = [x for x in xs for y in ys]
        y_grid = [y for x in xs for y in ys]
    @eval (@__MODULE__) begin
        using Plots
        theme(:orange)
        display(plot($x_grid, $y_grid, $U, st = :surface; $kw...))
    end
end

function plot_1d_solution(
    x::Vector{T},
    y::Vector{T};
    kw...,
) where {T<:Number}
    @eval (@__MODULE__) begin
        using Plots
        theme(:orange)
        display(plot($x, $y; $kw...))
    end
end

function conjugate_gradient_inverse_divide(A, b)
    f = x -> (x' * A / 2 - b') * x
    g! = (G, x) -> (G[:] = (A * x - b))
    return optimize(f, g!, zeros(length(b)), ConjugateGradient()).minimizer
end

function norm_for_1d_grid(
    f::Vector{T},
    p::Real,
    int_size::Vector{T} = 1 / (len(f) + 1),
    k::Int = 0;
    is_norm::Bool = true,
) where {T<:Number}
    if k === 0
        f_all = f
    elseif k === 1
        f_all =
            [
                f[2] - f[1]
                f[3:end] - f[1:end-2]
                f[end] - f[end-1]
            ] ./ int_size
        if is_norm
            f_all = [f_all; f[:]]
        end
    end
    if p === Inf
        return norm(f_all, Inf)
    else
        return sum(abs.(f_all) .^ p .* int_size) .^ (1 / p)
    end
end

"""
I do not know how to write the version that `k` greater than 1.

is_norm: norm or seminorm
"""
function norm_for_grid(
    f::Vector{T},
    p::Real,
    len::Int,
    k::Int = 0,
    wid::Int = len;
    is_norm::Bool = true,
) where {T<:Number}
    if k === 0
        f_all = f
    elseif k === 1
        f = reshape(f, len, wid)
        f_1 =
            [
                f[2:2, :] - f[1:1, :]
                f[3:end, :] - f[1:end-2, :]
                f[end:end, :] - f[end-1:end-1, :]
            ] .* (1 + len)
        f_2 =
            [f[:, 2] - f[:, 1] f[:, 3:end] - f[:, 1:end-2] f[:, end] - f[:, end-1]] .*
            (1 + wid)
        f_all = [f_1[:]; f_2[:]]
        if is_norm
            f_all = [f_all; f[:]]
        end
    end
    if p === Inf
        return norm(f_all, Inf)
    else
        return sum(abs.(f_all) .^ p ./ ((1 + len) * (1 + wid))) .^ (1 / p)
    end
end

include("PoissonsEquation.jl")

include("HeatEquation.jl")

include("ConvectionDiffusionEquation.jl")

include("FiniteElementMethod.jl")

include("SpectralMethod.jl")

end
