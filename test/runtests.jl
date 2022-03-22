using NumericalSolutionOfPartialDifferentialEquations
using Test

@testset "NumericalSolutionOfPartialDifferentialEquations.jl" begin
    @test test_solve_poissions_equation_known() < 0.0008358025108066159
    log2_err = test_solve_poissions_equation_unknown()[2]
    @test abs(2 + log2_err[end] - log2_err[end-1]) < 0.00014
end
