using MadMax
using Test

# @testset "MadMax.jl" begin
#     # Write your tests here.
# end

X = [ 2 -1  0 -7
     -1  2 -1 -1 
      0 -1  2 -1]
    
x, a = GaussSeidelEval(X, [0;0;0], N=5)
x, a = GaussJacobiEval(X, [0;0;0], N=5)
x, a = SOREval(X, [5.5;4.5;3.5], 1.1716, N=2)
println(x)