using NumericalSolutionOfPartialDifferentialEquations
import NumericalSolutionOfPartialDifferentialEquations.u_alpha
import NumericalSolutionOfPartialDifferentialEquations.f_alpha
using Test

@testset "NumericalSolutionOfPartialDifferentialEquations.jl" begin
    @test test_solve_poissions_equation_known() < 0.000835802511
    log2_err = test_solve_poissions_equation_unknown()[2]
    @test abs(2 + log2_err[end] - log2_err[end-1]) < 0.00014

    bound_func = [(x, y) -> 0, (x, y) -> 0, (x, y) -> 0, (x, y) -> 0]   # this argument is currently useless for FEM, only zero boundary is supported
    @test test_solve_poissions_equation_known(;
        scheme = :FEM,
        bound_func = bound_func,
        f = (x, y) -> 2 * pi^2 * sin(pi * x) * sin(pi * y),
        u = (x, y) -> sin(pi * x) * sin(pi * y),
    ) < 0.00240531659
    @test test_solve_poissions_equation_known(;
        scheme = :FEM,
        bound_func = bound_func,
        f = (x, y) -> -2.0 * (x^2 + y^2 - x - y),
        u = (x, y) -> x * (x - 1) * y * (y - 1),
    ) < 0.000162601471
    @test test_solve_poissions_equation_known(;
        scheme = :FEM,
        bound_func = [u_alpha, u_alpha, u_alpha, u_alpha],
        f = f_alpha,
        u = u_alpha,
        norm_p = 2,
        norm_k = 1,
        is_norm = false,
        integral_method = :HCubature
    ) < 0.0105706
end
