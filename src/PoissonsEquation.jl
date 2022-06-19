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
function construct_grid(
    f::T,
    len::Int,
    wid::Int = len;
    scheme = :FDM,
    kw...,
) where {T<:Function}
    if scheme === :FEM
        return construct_grid_FEM(f, len, wid)
    else
        return f.(range(0.0, 1.0, wid + 2)[2:end-1]', range(1.0, 0.0, len + 2)[2:end-1])[:]
    end
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
    boundary::Vector{Vector{T2}};
    solver::Function = \,
    kw...,
) where {T1<:Number,T2<:Number}
    size = length(boundary[1]) - 2
    f_with_boundary = add_boundary_condition(f ./ ((size + 1)^2), boundary)
    lap = construct_laplacian(; size = size)
    return solver(-lap, f_with_boundary)
end

function solve_poissions_equation(
    f::T1,
    bound_func::Vector{T2},
    len::Int,
    wid::Int = len;
    scheme::Symbol = :FDM,
    kw...,
) where {T1<:Function,T2<:Function}
    return solve_poissions_equation(
        construct_grid(f, len, wid),
        construct_boundary(bound_func, len, wid),
        scheme;
        fun = f,
        len = len,
        wid = wid,
        kw...,
    )
end

function solve_poissions_equation(
    f::Vector{T1},
    bound_func::Vector{T2},
    len::Int,
    wid::Int = len;
    scheme::Symbol = :FDM,
    kw...,
) where {T1<:Number,T2<:Function}
    return solve_poissions_equation(
        f,
        construct_boundary(bound_func, len, wid),
        scheme;
        kw...,
    )
end

function solve_poissions_equation(
    f::T1,
    boundary::Vector{Vector{T2}},
    len::Int,
    wid::Int = len;
    scheme::Symbol = :FDM,
    kw...,
) where {T1<:Function,T2<:Number}
    return solve_poissions_equation(
        construct_grid(f, len, wid),
        boundary,
        scheme;
        fun = f,
        len = len,
        wid = wid,
        kw...,
    )
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
    norm_p = Inf,
    kw...,
)
    U = solve_poissions_equation(f, bound_func, size; kw...)
    plotgui && plot_2d_solution(U, size; kw...)
    return norm_for_grid(U - construct_grid(u, size; kw...), norm_p, size)
end

"""
    # 无精确解 log-log plot
"""
function test_solve_poissions_equation_unknown(;
    max_size = 511,
    plotgui = false,
    bound_func = [(x, y) -> 0, (x, y) -> 0, (x, y) -> 0, (x, y) -> 0],
    f = (x, y) -> sinc(4 * x * y),
    kw...,
)
    log2_max_size = floor(Int, log2(max_size + 1))
    U = [solve_poissions_equation(f, bound_func, 1)]

    if log2_max_size < 2
        return U
    end

    log2_err = zeros(log2_max_size - 1)
    for i = 2:log2_max_size
        size = 2^i - 1
        push!(U, solve_poissions_equation(f, bound_func, size; kw...))
        log2_err[i-1] =
            log2(norm(reshape(U[i], size, size)[2:2:end, 2:2:end][:] - U[i-1], Inf))
    end

    if plotgui
        plot_index = log2_max_size > 5 ? 5 : log2_max_size
        plot_2d_solution(U[plot_index], 2^plot_index - 1; kw...)
    end

    return U, log2_err
end

"""
Need ~500 GiB Memory!!!
"""
function homework1()
    for i = 1:14
        print(test_solve_poissions_equation_known(; size = 2^i - 1))
        print("\n")
    end
    log2_err = test_solve_poissions_equation_unknown(; max_size = 2^14 - 1)[2]
    print(log2_err)

    test_solve_poissions_equation_known(;
        plotgui = true,
        xlabel = "x",
        ylabel = "y",
        zlabel = "U",
        reuse = false,
    )
    test_solve_poissions_equation_unknown(;
        plotgui = true,
        xlabel = "x",
        ylabel = "y",
        zlabel = "U",
        reuse = false,
    )
    theme(:default)
    b, k = [ones(13) [1:13;]] \ -log2_err
    scatter(
        [1:13;],
        -log2_err;
        markersize = 8,
        label = "-\\log_2\\lrvv{U_h - U_{h / 2}}_\\infty",
        xlabel = "-\\log_2 h",
        reuse = false,
    )
    plot!(x -> (k * x + b), [1:13;]; label = "fitting curve")
end
