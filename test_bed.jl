using Krylov
using LinearAlgebra
using Random

include("random_matrix.jl")
include("blendenpik_gauss.jl")

n = 1000
d = 20
κ = 1e29

A = generate_matrix(n, d, d, κ)
b = randn(n)

x_lu = A \ b
r_lu = A*x_lu - b

x_bg = blendenpick_gauss(A, b);
r_bg = A*x_bg - b

printstyled("""
SOLVING: $n × $d full column rank least squares problem.

LU SOLVER:
\tResidual Norm: $(norm(r_lu))
\tNormal System Residual Norm: $(norm(A'*r_lu))

BLENDENPIK_GAUSS:
\tResidual Norm: $(norm(r_bg))
\tNormal System Residual Norm: $(norm(A'*r_bg))
"""
)