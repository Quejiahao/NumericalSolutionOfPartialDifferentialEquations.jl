function construct_chebyshev_matrix(N::Int)
    if N === 0
        D = 0.0
        x = 1.0
    else
        x = -cos.(pi .* (0:N) ./ N)
        c = [2; ones(N - 1, 1); 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N + 1)
        dX = X - X'
        D = (c * (1 ./ c)') ./ (dX + I)
        diag_ind = diagind(D)
        D[diag_ind] = D[diag_ind] - sum(D; dims = 2)
    end
    return D, x
end

function homework6(;
    max_log2_size = 11,
    f = x -> 2 * sin(x),
    u = sin,
    norm_p = Inf,
    norm_k = 0,
    is_norm = true,
    plotgui = false,
)
    min_log2_size = 2
    err_FEM = zeros(max_log2_size - min_log2_size + 1)
    err_SM = zeros(max_log2_size - min_log2_size + 1)
    Ts_FEM = zeros(max_log2_size - min_log2_size + 1)
    Ts_SM = zeros(max_log2_size - min_log2_size + 1)
    for i = min_log2_size:max_log2_size
        N = 2^i - 1
        """
            Finite Element method
        """
        t = time()
        K = construct_1d_stiffness_matrix(N)[2:N, 2:N]
        F = construct_1d_load_vector(N, f)[2:N]
        U = K \ F
        xs = (1:N-1) .* pi ./ N
        U_true = u.(xs)
        if i == 3 && plotgui
            # plot_1d_solution(Array(xs), U)
        end
        err_FEM[i-min_log2_size+1] = norm_for_1d_grid(
            U - U_true,
            norm_p,
            repeat([pi / N], N - 1),
            norm_k;
            is_norm = is_norm,
        )
        Ts_FEM[i-min_log2_size+1] = time() - t
        """
            Spectral method
        """
        t = time()
        D, x = construct_chebyshev_matrix(N)
        D = D ./ (pi / 2)
        x = pi / 2 .+ x .* (pi / 2)
        K = -D^2 + I
        F = f.(x)
        K = K[2:N, 2:N]
        F = F[2:N]
        U = K \ F
        U_true = u.(x[2:N])
        if i == 3 && plotgui
            plot_1d_solution(x[2:N], U)
        end
        err_SM[i-min_log2_size+1] = norm_for_1d_grid(
            U - U_true,
            norm_p,
            x[2:N] - x[1:N-1],
            norm_k;
            is_norm = is_norm,
        )
        Ts_SM[i-min_log2_size+1] = time() - t
    end
    h = string.("2^{", string.([-min_log2_size:-1:-max_log2_size;]), "}")
    log2_err_FEM = log2.(err_FEM)
    log2_err_SM = log2.(err_SM)
    order_FEM = [
        0
        log2_err_FEM[1:(end-1)] - log2_err_FEM[2:end]
    ]
    order_SM = [
        0
        log2_err_SM[1:(end-1)] - log2_err_SM[2:end]
    ]
    return [h err_FEM order_FEM Ts_FEM err_SM order_SM Ts_SM]
end
