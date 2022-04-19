construct_conv_diff_grid = construct_heat_grid

function con_pre_wave_factor_upwind(;
    space_step_num::Int = 31,
    nu::T = 0.5,
) where {T<:Number}
    return nu == 0.0 ? I : con_tri_diag(space_step_num, 1.0 - nu, nu, 0.0)
end

function con_pre_wave_factor_lax_wendroff(;
    space_step_num::Int = 31,
    nu::T = 0.5,
) where {T<:Number}
    return nu == 0.0 ? I :
           con_tri_diag(
        space_step_num,
        1.0 - nu^2,
        0.5 * nu * (nu + 1.0),
        0.5 * nu * (nu - 1.0),
    )
end

function con_pre_wave_factor_beam_warming(;
    space_step_num::Int = 31,
    nu::T = 0.5,
) where {T<:Number}
    size = space_step_num
    main_diag = 0.5 * (nu - 1.0) * (nu - 2.0)
    first_above = nu * (2 - nu)
    second_above = 0.5 * nu * (nu - 1.0)
    Is = [1:size; 1:(size-1); 1:(size-2)]
    Js = [1:size; 2:size; 3:size]
    Vs = [
        fill(main_diag, size)
        fill(first_above, size - 1)
        fill(second_above, size - 2)
    ]
    return sparse(Is, Js, Vs)
end

function solve_conv_diff_equation(
    velocity::T1,
    initial::Vector{T2};
    solver::Function = \,
    total_time::T3 = 1.0,
    time_step_num::Int = 1024,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number}
    space_step_num = length(initial)
    time_step = total_time / time_step_num
    space_step = 1.0 / (space_step_num + 1)
    @show nu = velocity * time_step / space_step
end
