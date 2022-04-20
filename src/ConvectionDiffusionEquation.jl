construct_conv_diff_grid = construct_heat_grid

function con_pre_wave_factor_upwind(;
    space_step_num::Int = 31,
    nu::T = 0.5,
    periodic::Bool = true,
    kw...,
) where {T<:Number}
    if nu > zero(nu)
        pre_wave_factor = con_tri_diag(space_step_num, 1.0 - nu, nu, 0.0)
        if periodic
            pre_wave_factor[end, 1] = nu
        end
    elseif nu < zero(nu)
        pre_wave_factor = con_tri_diag(space_step_num, 1.0 + nu, 0.0, -nu)
        if periodic
            pre_wave_factor[1, end] = -nu
        end
    else
        return I
    end
    return pre_wave_factor
end

function con_pre_wave_factor_lax_wendroff(;
    space_step_num::Int = 31,
    nu::T = 0.5,
    periodic::Bool = true,
    kw...,
) where {T<:Number}
    if nu > zero(nu)
        pre_wave_factor = con_tri_diag(
            space_step_num,
            1.0 - nu^2,
            0.5 * nu * (nu + 1.0),
            0.5 * nu * (nu - 1.0),
        )
        if periodic
            pre_wave_factor[end, 1] = 0.5 * nu * (nu + 1.0)
            pre_wave_factor[1, end] = 0.5 * nu * (nu - 1.0)
        end
    elseif nu < zero(nu)
        error(
            "Velocity must be non-negative for Lax-Wendroff scheme, negative version has not been complemented yet!",
        )
    else
        return I
    end
    return pre_wave_factor
end

function con_pre_wave_factor_beam_warming(;
    space_step_num::Int = 31,
    nu::T = 0.5,
    periodic::Bool = true,
    kw...,
) where {T<:Number}
    if nu < zero(nu)
        error(
            "Velocity must be non-negative for Lax-Wendroff scheme, negative version has not been complemented yet!",
        )
    elseif nu === zero(nu)
        return I
    end
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
    if periodic
        Is = [Is; size - 1; size; size]
        Js = [Js; 1; 2; 1]
        Vs = [Vs; second_above; second_above; first_above]
    end
    return sparse(Is, Js, Vs)
end

function con_pre_wave_factor(; scheme::Symbol = :lax_wendroff, kw...)
    if scheme === :upwind
        return con_pre_wave_factor_upwind(; kw...)
    elseif scheme == :lax_wendroff
        return con_pre_wave_factor_lax_wendroff(; kw...)
    elseif scheme == :beam_warming
        return con_pre_wave_factor_beam_warming(; kw...)
    else
        error("Unknown scheme: %s", scheme)
    end
end

function solve_conv_diff_equation(
    velocity::T1,
    initial::Vector{T2};
    total_time::T3 = 1.0,
    time_step_num::Int = 1024,
    total_space::T4 = 1.0,
    scheme::Symbol = :lax_wendroff,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}
    space_step_num = length(initial)
    time_step = total_time / time_step_num
    space_step = total_space / (space_step_num + 1)
    @show nu = velocity * time_step / space_step
    result = [initial zeros(space_step_num, time_step_num)]
    pre_wave_factor =
        con_pre_wave_factor(; scheme = scheme, space_step_num = space_step_num, nu = nu)
    for i = 1:time_step_num
        result[:, i+1] = pre_wave_factor * result[:, i]
    end
    return result
end

function solve_conv_diff_equation(
    velocity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Function}
    return solve_conv_diff_equation(velocity, construct_initial(initial; kw...); kw...)
end

plot_1d_conv_diff_equ_solution = plot_1d_heat_equ_solution

function test_solve_conv_diff_equation(;
    velocity = 1.0,
    initial = x -> sin(pi * x),
    u = (x, t) -> sin(pi * x - velocity * pi * t),
    total_time = 2.0,
    time_step_num = 1024,
    total_space = 2.0,
    space_step_num = 31,
    plotgui = false,
    kw...,
)
    U = Dict()
    err = Dict()
    U[:true_solution] = construct_conv_diff_grid(
        u;
        time_step_num = time_step_num,
        space_step_num = space_step_num,
        total_time = total_time,
        total_space = total_space,
    )
    if plotgui
        plot_1d_conv_diff_equ_solution(
            U[:true_solution][:],
            space_step_num,
            time_step_num;
            total_space = total_space,
            total_time = total_time,
            reuse = false,
            kw...,
        )
    end
    test_schemes = [:upwind, :lax_wendroff, :beam_warming]
    for scheme in test_schemes
        U[scheme] = solve_conv_diff_equation(
            velocity,
            initial;
            total_time = total_time,
            time_step_num = time_step_num,
            total_space = total_space,
            space_step_num = space_step_num,
            scheme = scheme,
            kw...,
        )
        err[scheme] = U[scheme] - U[:true_solution]
        if plotgui
            plot_1d_conv_diff_equ_solution(
                U[scheme][:],
                space_step_num,
                time_step_num;
                total_space = total_space,
                total_time = total_time,
                reuse = false,
                kw...,
            )
        end
    end
    return U, err
end

function homework3(; kw...,)
    test_solve_conv_diff_equation(; kw...,)
end
