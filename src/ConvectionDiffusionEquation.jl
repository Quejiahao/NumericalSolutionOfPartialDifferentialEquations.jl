function construct_conv_diff_grid(
    u::T1;
    time_step_num::Int = 1024,
    space_step_num::Int = 31,
    total_time::T2 = 1.0,
    total_space::T3 = 1.0,
    kw...,
) where {T1<:Function,T2<:Number,T3<:Number}
    return u.(
        range(0.0, total_space, space_step_num + 2),
        range(0.0, total_time, time_step_num + 1)',
    )
end

function con_pre_wave_factor_upwind(;
    space_step_num::Int = 31,
    nu::T = 0.5,
    periodic::Bool = true,
    kw...,
) where {T<:Number}
    if nu > zero(nu)
        pre_wave_factor = con_tri_diag(space_step_num, 1.0 - nu, 0.0, nu)
        if periodic
            pre_wave_factor[1, end-1] = nu
        end
    elseif nu < zero(nu)
        pre_wave_factor = con_tri_diag(space_step_num, 1.0 + nu, -nu, 0.0)
        if periodic
            pre_wave_factor[end, 2] = -nu
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
            0.5 * nu * (nu - 1.0),
            0.5 * nu * (nu + 1.0),
        )
        if periodic
            pre_wave_factor[1, end-1] = 0.5 * nu * (nu + 1.0)
            pre_wave_factor[end, 2] = 0.5 * nu * (nu - 1.0)
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
    first_below = nu * (2 - nu)
    second_below = 0.5 * nu * (nu - 1.0)
    Is = [1:size; 2:size; 3:size]
    Js = [1:size; 1:(size-1); 1:(size-2)]
    Vs = [
        fill(main_diag, size)
        fill(first_below, size - 1)
        fill(second_below, size - 2)
    ]
    if periodic
        Is = [Is; 1; 2; 1]
        Js = [Js; size - 2; size - 1; size - 1]
        Vs = [Vs; second_below; second_below; first_below]
    end
    return sparse(Is, Js, Vs)
end

function con_pre_wave_factor(; scheme::Symbol = :lax_wendroff, kw...)
    if scheme === :upwind
        return con_pre_wave_factor_upwind(; kw...)
    elseif scheme === :lax_wendroff
        return con_pre_wave_factor_lax_wendroff(; kw...)
    elseif scheme === :beam_warming
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
    flow_in_out::Bool = false,
    periodic::Bool = true,
    kw...,
) where {T1<:Number,T2<:Number,T3<:Number,T4<:Number}
    if flow_in_out
        periodic = false
    end
    space_step_num = length(initial)
    time_step = total_time / time_step_num
    space_step = total_space / (space_step_num - 1) # -1 for boundary
    @show nu = velocity * time_step / space_step
    result = [initial zeros(space_step_num, time_step_num)]
    pre_wave_factor = con_pre_wave_factor(;
        scheme = scheme,
        space_step_num = space_step_num,
        nu = nu,
        periodic = periodic,
    )
    for i = 1:time_step_num
        result[:, i+1] = pre_wave_factor * result[:, i]
        if flow_in_out
            if scheme === :beam_warming
                result[2, i+1] = result[3, i+1]
            end
            result[1, i+1] = result[2, i+1]
            if scheme === :lax_wendroff
                result[end, i+1] = result[end-1, i+1]
            end
        end
    end
    return result
end

function solve_conv_diff_equation(
    velocity::T1,
    initial::T2;
    kw...,
) where {T1<:Number,T2<:Function}
    return solve_conv_diff_equation(
        velocity,
        construct_initial(initial; boundary = true, kw...);
        kw...,
    )
end

plot_1d_conv_diff_equ_solution(U, space_step_num, time_step_num; kw...) =
    plot_1d_heat_equ_solution(U, space_step_num, time_step_num; boundary = true, kw...)

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

function solve_finite_volume_riemann_problem(;
    x_max = 4.0,
    t_max = 1.0,
    a1 = 2.0,
    a2 = 1.0,
    time_step = 2^(-6),
    space_step = 2^(-5),
)
    space_step_num = ceil(Int, 2 * x_max / space_step + 1)
    time_step_num = ceil(Int, t_max / time_step + 1)
    result = zeros(space_step_num, time_step_num)
    result[1:round(Int, space_step_num / 2), 1] .= a1
    result[round(Int, space_step_num / 2)+1:end, 1] .= a2
    v = time_step / space_step
    function _flux(ul, ur)
        if ul >= 0 && ur >= 0
            u2 = ul
        elseif ul <= 0 && ur <= 0
            u2 = ur
        elseif ul >= 0 && ur <= 0
            if ul + ur >= 0
                u2 = ul
            elseif ul + ur < 0
                u2 = ur
            end
        elseif ul <= 0 && ur > 0
            u2 = 0
        end
        return u2^2 / 2
    end
    for i = 1:time_step_num-1
        F = _flux.(result[1:end-1, i], result[2:end, i])
        result[2:end-1, i+1] = result[2:end-1, i] + v * (F[1:end-1] - F[2:end])
        result[1, i+1] = result[2, i+1]
        result[end, i+1] = result[end-1, i+1]
    end
    return result
end

function homework3(; nu = 0.5, i_max = 10, kw...)
    sine_wave = Dict(
        :velocity => 1.0,
        :initial => x -> sin(pi * x),
        :u => (x, t) -> sin(pi * (x - t)),
        :total_space => 2.0,
        :space_step_num => 31,
        :total_time => 1.0,
        :time_step_num => 512,
        :plotgui => false,
    )

    _total_space_square_wave = 3.0
    function _u_square_wave(x, t)
        # (mod(x - t, _total_space_square_wave) > 1.0) ? 1.0 : 3.0
        xt = x - t
        # return (xt <= 1.0) && (xt >= 0.0) ? 3.0 : 1.0
        return (xt <= 1.0) ? 3.0 : 1.0
    end
    square_wave = Dict(
        :velocity => 1.0,
        :initial => x -> ((x > 1.0) ? 1.0 : 3.0),
        :u => _u_square_wave,
        :total_space => _total_space_square_wave,
        :space_step_num => 31,
        :total_time => 0.6,
        :time_step_num => 512,
        :plotgui => false,
        :flow_in_out => true,
    )
    # triangle_wave = Dict(

    results = []
    for i = 4:i_max
        # @show sine_wave[:space_step_num] = 2^i - 1
        # sine_wave[:time_step_num] = round(
        #     Int,
        #     sine_wave[:total_time] * (sine_wave[:space_step_num] + 1) /
        #     sine_wave[:total_space] / nu,
        # )
        # push!(results, test_solve_conv_diff_equation(; sine_wave..., kw...))
        @show square_wave[:space_step_num] = 2^i - 1
        square_wave[:time_step_num] = round(
            Int,
            square_wave[:total_time] * (square_wave[:space_step_num] + 1) /
            square_wave[:total_space] / nu,
        )
        push!(results, test_solve_conv_diff_equation(; square_wave..., kw...))
    end

    results_err = []
    test_schemes = [:upwind, :lax_wendroff, :beam_warming]
    for scheme in test_schemes
        l_infty_err = [
            norm(results[i][1][scheme] - results[i][1][:true_solution], Inf) for
            i = 1:length(results)
        ]
        l_infty_err_log2 = log2.(l_infty_err)
        l_infty_order = [0; l_infty_err_log2[1:end-1] - l_infty_err_log2[2:end]]
        l_2_err = [
            norm(
                norm.(eachcol(results[i][1][scheme] - results[i][1][:true_solution]),),
                Inf,
            ) / sqrt(2^(i + 3)) for i = 1:length(results)
        ]
        l_2_err_log2 = log2.(l_2_err)
        l_2_order = [0; l_2_err_log2[1:end-1] - l_2_err_log2[2:end]]
        push!(results_err, [l_2_order l_2_err [-4:-1:-i_max;] l_infty_err l_infty_order])
    end

    if 0 === 1
        for scheme in test_schemes
            for i in [4, 7]
                xs = [range(0.0, 2 * pi, 2^(i + 3) + 1);]
                plot(
                    xs,
                    [results[i][1][scheme][:, end], results[i][1][:true_solution][:, end]],
                    xtickfontsize = 18,
                    ytickfontsize = 18,
                    xlabel = "\$x\$",
                    xguidefontsize = 18,
                    label = ["\$U\$" "\$u\$"],
                    legendfontsize = 18,
                )
            end
        end
        a1 = 1.0
        a2 = 0.0
        xs = [-4.0:2^(-5):4.0;]
        function _y(x)
            if x <= 0.5
                return 1.0
            elseif x <= 0.5
                return x
            else
                return 0.0
            end
        end
        function _z(x)
            if x <= 0.0
                return 1.0
            elseif x <= 0.0
                return x
            else
                return 0.0
            end
        end
        ys = _y.(xs)
        zs = _z.(xs)
        plot(
            xs,
            [solve_finite_volume_riemann_problem(; a1 = a1, a2 = a2)[:, end], zs, ys],
            xtickfontsize = 18,
            ytickfontsize = 18,
            xlabel = "\$x\$",
            xguidefontsize = 18,
            label = ["Conservation" "Upwind" "\$u\$"],
            legendfontsize = 18,
        )
    end
    return results, results_err
end
