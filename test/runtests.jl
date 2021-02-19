using MadMax
using Test

# @testset "MadMax.jl" begin
#     # Write your tests here.
# end

X = [4 1 1 -2
     1 5 2 6
     1 2 3 4 ]
    
x, a = GaussJacobiEval(X, [0;0;0], N=50)
println(x)