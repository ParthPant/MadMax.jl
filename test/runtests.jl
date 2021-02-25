using MadMax
using Test

# @testset "MadMax.jl" begin
#     # Write your tests here.
# end

X = [-4 -1  0
     -1  4 -1]
    
x, a = GaussSeidelEval(X, [0.25;0.25], N=3)
x, a = GaussJacobiEval(X, [0.25;0.25], N=7)
x, a = SOREval(X, [0.25;0.25], 1.1716, N=2)

for i in a
      println(i)
end
