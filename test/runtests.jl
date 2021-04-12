using MadMax
using Test

# BisectionEval
# RegulaFalsiEval
# NewtonEval
# GaussJacobiEval
# GaussSeidelEval
# SOREval
# SOREval2
# EulerODE
# ModifiedEulerODE

X = [5 5 15 -40
      5 15 35 -104
      15 35 99 -290]
    
X2 = [3  2  0 -5
      2  3 -1 -4
      0 -1  2 -1]

x, a = GaussSeidelEval(X2, [0.9;-3.1;0.9], N=3)
x, a = GaussJacobiEval(X, [0;0;0], N=7)
x, a = SOREval(X2, [0;0;0], 1.25, N=3)
x, a = SOREval2(X2, [0;0;0], 0.8, N=7)

f(x) = x^24-2.6x+1.6 
RegulaFalsiEval(f, 0, 1, tol=0.5*10^-3)
NewtonEval(f, 0.05, N=4, tol=0.5*10^-4)

f(x, y) = x^2+y
ModifiedEulerODE(f, 0, 1, 0.1, 0.05)