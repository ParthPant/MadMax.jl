using MadMax
using Test

# @testset "MadMax.jl" begin
#     # Write your tests here.
# end

f(x) = x^2-6
BisectionEval(f, 2, 3)
RegulaFalsiEval(f, 2, 3)