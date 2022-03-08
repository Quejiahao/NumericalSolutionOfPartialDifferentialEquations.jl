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
    f_with_boundary = add_boundary_condition(f ./ ((size + 1) ^ 2), boundary)
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

function test_solve_poissions_equation(; size = 32, plotgui = false)
    bound_func = [
        (x, y) -> sin(pi * y),
        (x, y) -> 0,
        (x, y) -> MathConstants.e * sin(pi * y),
        (x, y) -> 0,
    ]
    f(x, y) = (pi^2 - 1) * exp(x) * sin(pi * y)
    u(x, y) = exp(x) * sin(pi * y)
    U = solve_poissions_equation(f, bound_func, size)
    err = norm(U - construct_grid(u, size), Inf)

    if plotgui
        wid = len = size
        xs = range(0.0, 1.0, wid + 2)[2:end-1]
        ys = range(1.0, 0.0, len + 2)[2:end-1]
        x_grid = [x for x in xs for y in ys]
        a_grid = [y for x in xs for y in ys]
        @eval (@__MODULE__) begin
            using Plots
            theme(:orange)
            display(plot($x_grid, $a_grid, $U, st = :surface))
        end
    end

    return err
end
