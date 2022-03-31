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

function con_heat_equ_diff_opt(;
    size::Int = 31,
    mu::T1 = 0.5,
    theta::T2 = 0.5,
) where {T1<:Number,T2<:Number}
    return theta == 0.0 ? I : con_tri_diag(size, 1.0 + 2 * mu * theta, -mu * theta)
end

function con_pre_temp_factor(;
    size::Int = 31,
    mu::T1 = 0.5,
    theta::T2 = 0.5,
) where {T1<:Number,T2<:Number}
    return theta == 1.0 ? I :
           con_tri_diag(size, 1.0 + 2 * mu * (theta - 1.0), mu * (1.0 - theta))
end

@doc raw"""
    solve_heat_equation(
        f::Vector{T1},
        boundary::Vector{Vector{T2}};
        solver::Function = \,
        kw...,
    ) where {T1<:Number,T2<:Number}

Solve Heat equation with Dirichlet boundary conditions.
"""
function solve_heat_equation(
    thermal_diffusivity::T1,
    initial::Vector{T2};
    solver::Function = \,
    total_time::T3 = 1.0,
    time_step_num::Int = 2048,
    theta::T4 = 0.5,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}
    size = length(initial)
    time_step = total_time / time_step_num
    space_step = 1.0 / (size + 1)
    mu = thermal_diffusivity * time_step / space_step^2
    result = [initial zeros(size, time_step_num)]
    heat_equ_diff_opt = con_heat_equ_diff_opt(; size = size, mu = mu, theta = theta)
    pre_temp_factor = con_pre_temp_factor(; size = size, mu = mu, theta = theta)
    for i = 1:time_step_num
        pre_temp_with_factor_and_boundary = pre_temp_factor * result[:, i] # + boundary # current zero boundary
        result[:, i+1] = solver(heat_equ_diff_opt, pre_temp_with_factor_and_boundary)
    end
    return result
end

function solve_heat_equation(
    thermal_diffusivity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Function}
    return solve_heat_equation(
        construct_grid(f, len, wid),
        construct_boundary(bound_func, len, wid);
        kw...,
    )
end

function solve_heat_equation(
    thermal_diffusivity::T1,
    initial::T2;
    space_steps::Int = 31,
    space_step_num::Int = 32,
    mu::T3 = 0.5,
    theta::T4 = 0.5,
    kw...,
) where {T1<:Number,T2<:Function,T3<:Number,T4<:Number}
    return solve_heat_equation(
        construct_grid(f, len, wid),
        construct_boundary(bound_func, len, wid);
        kw...,
    )
end

function homework2() end
