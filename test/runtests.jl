using MadMax
using Test

# @testset "MadMax.jl" begin
#     # Write your tests here.
# end

X = [ 3 -6  2 -23
     -4  1 -1  8
      1 -3  7 -17]
    
X2 = [3  2  0 -5
      2  3 -1 -4
      0 -1  2 -1]

x, a = GaussSeidelEval(X2, [0.9;-3.1;0.9], N=3)
x, a = GaussJacobiEval(X, [0.9;-3.1;0.9], N=7)
x, a = SOREval(X2, [0;0;0], 1.25, N=3)
x, a = SOREval2(X2, [0;0;0], 1.25, N=3)

for i in a
      println(i)
end
