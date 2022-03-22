"""
    add_boundary_condition(
        f::Vector{T1},
        boundary::Vector{Vector{T2}},
    ) where {T1<:Number,T2<:Number}

Add Dirichlet boundary conditions to the right of difference equation.

    Left:   boundary[1]
    Bottom: boundary[2]
    Right:  boundary[3]
    Top:    boundary[4]

They are column vectors and the orientation is counterclockwise.
"""
function add_boundary_condition(
    f::Vector{T1},
    boundary::Vector{Vector{T2}},
) where {T1<:Number,T2<:Number}
    len = length(boundary[1]) - 2
    wid = length(boundary[2]) - 2

    if length(f) !== len * wid
        error("The length of `f` must be compatible to boundary.")
    end

    f_with_boundary = f
    f_with_boundary[1:len] += boundary[1][2:(end-1)]
    f_with_boundary[len:len:end] += boundary[2][2:(end-1)]
    f_with_boundary[end:(-1):(end-len+1)] += boundary[3][2:(end-1)]
    f_with_boundary[(end-len+1):(-len):1] += boundary[4][2:(end-1)]

    return f_with_boundary
end

"""
    construct_boundary(
        bound_func::Vector{T},
        len::Int,
        wid::Int = len,
    ) where {T<:Function}

Use `bound_func` to construct boundary conditions.
"""
function construct_boundary(
    bound_func::Vector{T},
    len::Int,
    wid::Int = len,
) where {T<:Function}
    return [
        bound_func[1].(0.0, range(1.0, 0.0, len + 2)),
        bound_func[2].(range(0.0, 1.0, wid + 2), 0.0),
        bound_func[3].(1.0, range(0.0, 1.0, len + 2)),
        bound_func[4].(range(1.0, 0.0, wid + 2), 1.0),
    ]
end

"""
    construct_grid(f::T, len::Int, wid::Int = len) where {T<:Function}

Use `f` to construct grid.
"""
function construct_grid(f::T, len::Int, wid::Int = len) where {T<:Function}
    return f.(range(0.0, 1.0, wid + 2)[2:end-1]', range(1.0, 0.0, len + 2)[2:end-1])[:]
end

"""
    solve_poissions_equation(
        f::Vector{T1},
        boundary::Vector{Vector{T2}},
    ) where {T1<:Number,T2<:Number}

Solve Poisson's equation with Dirichlet boundary conditions.
"""
function solve_poissions_equation(
    f::Vector{T1},
    boundary::Vector{Vector{T2}},
) where {T1<:Number,T2<:Number}
    size = length(boundary[1]) - 2
    f_with_boundary = add_boundary_condition(f ./ ((size + 1)^2), boundary)
    lap = construct_laplacian(; size = size)
    return -lap \ f_with_boundary
end

function solve_poissions_equation(
    f::T1,
    bound_func::Vector{T2},
    len::Int,
    wid::Int = len,
) where {T1<:Function,T2<:Function}
    return solve_poissions_equation(
        construct_grid(f, len, wid),
        construct_boundary(bound_func, len, wid),
    )
end

function solve_poissions_equation(
    f::Vector{T1},
    bound_func::Vector{T2},
    len::Int,
    wid::Int = len,
) where {T1<:Number,T2<:Function}
    return solve_poissions_equation(f, construct_boundary(bound_func, len, wid))
end

function solve_poissions_equation(
    f::T1,
    boundary::Vector{Vector{T2}},
    len::Int,
    wid::Int = len,
) where {T1<:Function,T2<:Number}
    return solve_poissions_equation(construct_grid(f, len, wid), boundary)
end

"""
    # 已知精确解
    # 共轭梯度?LU分解?
"""
function test_solve_poissions_equation_known(;
    size = 31,
    plotgui = false,
    bound_func = [
        (x, y) -> sin(pi * y),
        (x, y) -> 0,
        (x, y) -> MathConstants.e * sin(pi * y),
        (x, y) -> 0,
    ],
    f = (x, y) -> (pi^2 - 1) * exp(x) * sin(pi * y),
    u = (x, y) -> exp(x) * sin(pi * y),
)
    U = solve_poissions_equation(f, bound_func, size)
    plotgui && plot_2d_solution(U, size)
    return norm(U - construct_grid(u, size), Inf)
end

"""
    # 无精确解 log-log plot
"""
function test_solve_poissions_equation_unknown(;
    max_size = 512,
    plotgui = false,
    bound_func = [(x, y) -> 0, (x, y) -> 0, (x, y) -> 0, (x, y) -> 0],
    f = (x, y) -> sinc(4 * x * y),
)
    log2_max_size = floor(Int, log2(max_size))
    U = [solve_poissions_equation(f, bound_func, 1)]

    if log2_max_size < 2
        return U
    end

    log2_err = zeros(log2_max_size - 1)
    for i = 2 : log2_max_size
        size = 2 ^ i - 1
        push!(U, solve_poissions_equation(f, bound_func, size))
        log2_err[i - 1] = log2(norm(reshape(U[i], size, size)[2:2:end,2:2:end][:] - U[i - 1], Inf))
    end

    if plotgui
        plot_index = log2_max_size > 5 ? 5 : log2_max_size
        plot_2d_solution(U[plot_index], 2 ^ plot_index - 1)
    end

    return U, log2_err
end
