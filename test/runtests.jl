using NumericalSolutionOfPartialDifferentialEquations
using Test

@testset "NumericalSolutionOfPartialDifferentialEquations.jl" begin
    @test test_solve_poissions_equation() < 0.0007856741026288639
end
