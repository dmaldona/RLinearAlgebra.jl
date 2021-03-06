n = 10
d = 10
κ = 8.0

A = RLinearAlgebra.generate_matrix(n, d, d, κ)
b = randn(n)
@testset "Call basic RPM via API" begin
    sol = LinearSolver(TypeRPM())
    x = rsolve(sol, A, b)
    x_lu = A\b
end

A = [1 0 0; 0 2 0; 0 0 3]
@testset "Test row distribution" begin
    p = RLinearAlgebra.distribution(SVDistribution(), A)
    @test sum(p) ≈ 1.0
    @test p[3] ≈ 0.5
    p = RLinearAlgebra.distribution(UFDistribution(), A)
    @test sum(p) ≈ 1.0
    @test p[3] ≈ 1/3
end
