@doc raw"""
隐式时间 $\approv 2~3$ 显式时间
$$\dfrac{\partial u}{\partial t} = a \dfrac{\partial^2 u}{\partial x^2}$$
第一类边界条件
三种
条件不满足时 0.5 1 0.5 2?
三对角解法 LU 分解
N L^\infty order L^2 order
sparse div???
"""

function construct_heat_grid(
    u::T1;
    time_step_num::Int = 1024,
    space_step_num::Int = 31,
    total_time::T2 = 1.0,
    total_space::T3 = 1.0,
    kw...,
) where {T1<:Function,T2<:Number,T3<:Number}
    return u.(
        range(0.0, total_space, space_step_num + 2)[2:end-1],
        range(0.0, total_time, time_step_num + 1)',
    )
end

function con_heat_equ_diff_opt(;
    space_step_num::Int = 31,
    mu::T1 = 0.5,
    theta::T2 = 0.5,
) where {T1<:Number,T2<:Number}
    return theta == 0.0 ? I :
           con_tri_diag(space_step_num, 1.0 + 2 * mu * theta, -mu * theta)
end

function con_pre_temp_factor(;
    space_step_num::Int = 31,
    mu::T1 = 0.5,
    theta::T2 = 0.5,
) where {T1<:Number,T2<:Number}
    return theta == 1.0 ? I :
           con_tri_diag(space_step_num, 1.0 + 2 * mu * (theta - 1.0), mu * (1.0 - theta))
end

@doc raw"""
    solve_heat_equation(
        thermal_diffusivity::T1,
        initial::Vector{T2};
        solver::Function = \,
        total_time::T3 = 1.0,
        time_step_num::Int = 1024,
        theta::T4 = 0.5,
        kw...,
    ) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}

Solve Heat equation with Dirichlet boundary conditions (zero boundary and given initial).
"""
function solve_heat_equation(
    thermal_diffusivity::T1,
    initial::Vector{T2};
    solver::Function = \,
    total_time::T3 = 1.0,
    time_step_num::Int = 1024,
    theta::T4 = 0.5,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}
    space_step_num = length(initial)
    time_step = total_time / time_step_num
    space_step = 1.0 / (space_step_num + 1)
    @show mu = thermal_diffusivity * time_step / space_step^2
    @show theta
    result = [initial zeros(space_step_num, time_step_num)]
    heat_equ_diff_opt =
        con_heat_equ_diff_opt(; space_step_num = space_step_num, mu = mu, theta = theta)
    pre_temp_factor =
        con_pre_temp_factor(; space_step_num = space_step_num, mu = mu, theta = theta)
    # pre_to_next_factor = solver(heat_equ_diff_opt, pre_temp_factor)
    for i = 1:time_step_num
        pre_temp_with_factor_and_boundary = pre_temp_factor * result[:, i] # + boundary # current zero boundary
        result[:, i+1] = solver(heat_equ_diff_opt, pre_temp_with_factor_and_boundary)
        # result[:, i+1] = pre_to_next_factor * result[:, i]
    end
    return result
end

function solve_heat_equation(
    thermal_diffusivity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Function}
    return solve_heat_equation(
        thermal_diffusivity,
        construct_initial(initial; kw...);
        kw...,
    )
end

"""
    Three Scheme
"""
function solve_heat_equation_explicit(
    thermal_diffusivity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Union{Function,Vector{<:Number}}}
    return solve_heat_equation(thermal_diffusivity, initial; theta = 0.0, kw...)
end

function solve_heat_equation_imexplicit(
    thermal_diffusivity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Union{Function,Vector{<:Number}}}
    return solve_heat_equation(thermal_diffusivity, initial; theta = 1.0, kw...)
end

function solve_heat_equation_crank_nicolson(
    thermal_diffusivity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Union{Function,Vector{<:Number}}}
    return solve_heat_equation(thermal_diffusivity, initial; theta = 0.5, kw...)
end

"""
    Also used by convection diffusion equation
"""
function plot_1d_heat_equ_solution(
    U::Vector{T1},
    space_step_num::Int,
    time_step_num::Int;
    total_space::T2 = 1.0,
    total_time::T3 = 1.0,
    boundary::Bool = false,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number}
    xs = range(0.0, total_space, space_step_num + 2)
    if !boundary
        xs = xs[2:end-1]
    end
    ts = range(0.0, total_time, time_step_num + 1)
    x_grid = [x for t in ts for x in xs]
    t_grid = [t for t in ts for x in xs]
    @eval (@__MODULE__) begin
        using Plots
        theme(:orange)
        display(plot($x_grid, $t_grid, $U, st = :surface; $kw...))
    end
end

@doc raw"""
$$(x, t) \in (0, 1)^2$$
"""
function test_solve_heat_equation(;
    thermal_diffusivity = 0.5,
    initial = x -> sin(pi * x),
    u = (x, t) -> sin(pi * x) * exp(-thermal_diffusivity * pi^2 * t),
    space_step_num::Int = 31,
    time_step_num = 1024,
    plotgui = false,
    kw...,
)
    U1 = solve_heat_equation_explicit(
        thermal_diffusivity,
        initial;
        space_step_num = space_step_num,
        time_step_num = time_step_num,
        kw...,
    )
    U2 = solve_heat_equation_imexplicit(
        thermal_diffusivity,
        initial;
        space_step_num = space_step_num,
        time_step_num = time_step_num,
        kw...,
    )
    U3 = solve_heat_equation_crank_nicolson(
        thermal_diffusivity,
        initial;
        space_step_num = space_step_num,
        time_step_num = time_step_num,
        kw...,
    )
    u_grid = construct_heat_grid(
        u;
        space_step_num = space_step_num,
        time_step_num = time_step_num,
    )

    if plotgui
        plot_1d_heat_equ_solution(
            U1[:],
            space_step_num,
            time_step_num;
            reuse = false,
            kw...,
        )
        plot_1d_heat_equ_solution(
            U2[:],
            space_step_num,
            time_step_num;
            reuse = false,
            kw...,
        )
        plot_1d_heat_equ_solution(
            U3[:],
            space_step_num,
            time_step_num;
            reuse = false,
            kw...,
        )
        plot_1d_heat_equ_solution(
            u_grid,
            space_step_num,
            time_step_num;
            reuse = false,
            kw...,
        )
    end

    err1 = U1 - u_grid
    err2 = U2 - u_grid
    err3 = U3 - u_grid
    return [
        U1,
        U2,
        U3,
        u_grid,
        norm(err1, Inf),
        norm(err2, Inf),
        norm(err3, Inf),
        norm(norm.(eachcol(err1)), Inf) / sqrt(space_step_num + 1),
        norm(norm.(eachcol(err2)), Inf) / sqrt(space_step_num + 1),
        norm(norm.(eachcol(err3)), Inf) / sqrt(space_step_num + 1),
    ]
end

function homework2(; mu = 10.0, i_num = 8)
    result = []
    thermal_diffusivity = 0.5
    for i = 3:i_num
        space_step_num = 2^i - 1
        time_step_num = ceil(Int, thermal_diffusivity * 4^i / mu)
        append!(
            result,
            [
                test_solve_heat_equation(;
                    thermal_diffusivity = thermal_diffusivity,
                    space_step_num = space_step_num,
                    time_step_num = time_step_num,
                ),
            ],
        )
    end
    err = [result[i][end-5:end] for i = eachindex(result)]
    log2_err = [log2.(x) for x in err]
    append!(result, [[err; log2_err]])
    return result
end
