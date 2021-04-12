module MadMax

using PrettyTables
using Flux
using LinearAlgebra

export BisectionEval
export RegulaFalsiEval
export NewtonEval
export GaussJacobiEval
export GaussSeidelEval
export SOREval
export SOREval2
export EulerODE
export ModifiedEulerODE

function getLowerTriangular(X)
    R = zeros(size(X)...)
    for i in 1:size(X)[1]
        for j in 1:size(X)[2]
            if i>j
                R[i,j] = X[i, j]
            end
        end
    end
    return R
end

function getUpperTriangular(X)
    R = zeros(size(X)...)
    for i in 1:size(X)[1]
        for j in 1:size(X)[2]
            if i<j
                R[i,j] = X[i, j]
            end
        end
    end
    return R
end

function SOREval(X, x0, ω=1; tol=0.005, N=5, verbose=true)
    A = X[:, 1:end-1]
    b = -X[:, end]

    D = Diagonal(A)
    U = getUpperTriangular(A)
    L = getLowerTriangular(A)
    
    H = inv(D+ω*L) * ((1-ω)D - ω*U)
    c = ω*inv(D+ω*L)*b

    if verbose
        println("H=")
        show(stdout, "text/plain", H)
        println("")
        println("")
        println("c=")
        show(stdout, "text/plain", c)
    end

    x = x0 
    approximations = [x0]

    n = 0
    while n != N 
        x = H*x + c
        approximations = [approximations; [x]]
        n += 1
    end

    if verbose
        println("")
        println("")
        for (i,a) in enumerate(approximations)
            print("$(i-1):")
            println(a)
        end
    end

    return x, approximations
end

function SOREval2(X, x0, ω=1; tol=0.005, N=5, verbose=true)
    A = X[:, 1:end-1]
    b = -X[:, end]

    D = Diagonal(A)
    U = getUpperTriangular(A)
    L = getLowerTriangular(A)
    
    x = x0 
    approximations = [x0]

    H = inv(D+ω*L)*ω
    if verbose
        println("H=")
        show(stdout, "text/plain", H)
        println("")
    end

    n = 0
    while n != N 
        x_k = last(approximations)
        r_k = b - A*x_k
        e_k = H*r_k

        x = e_k  + x_k
    
        if verbose
            header = ["r_$(n+1)" "e_$(n+1)" "x_$(n+1)"]
            data = hcat(r_k, e_k, x)
            pretty_table(data,  header)
        end

        approximations = [approximations; [x]]
        n += 1
    end

    return x, approximations
end

function GaussJacobiEval(X, x0; tol=0.005, N=5, verbose=true)
    A = X[:, 1:end-1]
    b = -X[:, end]

    D = Diagonal(A)
    U = getUpperTriangular(A)
    L = getLowerTriangular(A)
    
    H = -inv(D) * (L+U)
    c = inv(D)*b

    if verbose
        println("H=")
        show(stdout, "text/plain", H)
        println("")
        println("")
        println("c=")
        show(stdout, "text/plain", c)
    end

    x = x0 
    approximations = [x0]

    n = 0
    while n != N 
        x = H*x + c
        approximations = [approximations; [x]]
        n += 1
    end

    if verbose
        println("")
        println("")
        for (i,a) in enumerate(approximations)
            print("$(i-1):")
            println(a)
        end
    end
    
    return x, approximations
end

function GaussSeidelEval(X, x0; tol=0.005, N=5, verbose=true)
    A = X[:, 1:end-1]
    b = -X[:, end]

    D = Diagonal(A)
    U = getUpperTriangular(A)
    L = getLowerTriangular(A)
    
    H = -inv(L+D)*U
    c = inv(L+D)*b

    if verbose
        println("H=")
        show(stdout, "text/plain", H)
        println("")
        println("")
        println("c=")
        show(stdout, "text/plain", c)
    end

    x = x0 
    approximations = [x0]

    n = 0
    while n != N 
        x = H*x + c
        approximations = [approximations; [x]]
        n += 1
    end

    if verbose
        println("")
        println("")
        for (i,a) in enumerate(approximations)
            print("$(i-1):")
            println(a)
        end
    end

    return x, approximations
end

function NewtonEval(f::Function, x0; tol = 0.005, N=100, verbose=true)
    error = Base.Inf64

    # find gradient 
    x = 0
    df(x) = gradient(f, x)[1]

    approximations = Array{Float64}(undef, 1)
    approximations = [x0]
    errors = [error]
    n = 0

    # iterate
    while error > tol && n != N
        x_last = last(approximations)
        delta = f(x_last) / df(x_last) 
        x_n = x_last - delta

        approximations = [approximations; x_n]

        # get error
        error = abs(x_last - x_n)
        errors = [errors; error]

        n += 1
    end

    if verbose
        headers = ["S.N" "x_n" "Error"]
        data = hcat(0:n, approximations, errors)
        
        pretty_table(data, headers)
        println("$n iterations")
    end

    return approximations, error
end

function BisectionEval(f::Function, a, b; tol = 0.005, N = 100, verbose=true)
    # make expr
    error = Base.Inf64 
    
    # find starting interval
    interval = (a, b)

    # iterate
    approximations = []
    errors = []
    n = 0

    if verbose
        as = []
        bs = []
        prod_s =[]
    end

    while error > tol && n!=N
        # get m_k
        m_k = sum(interval)/2

        # get error
        if length(approximations) > 0
            error = abs(last(approximations) - m_k)
        end
        push!(errors, error)

        if verbose
            push!(as, interval[1])
            push!(bs, interval[2])
        end

        # append the appproximation
        push!(approximations, m_k)

        # get new interval
        interval = (f(interval[1])*f(m_k)) < 0 ? (interval[1], m_k) : (m_k, interval[2])

        prod = (f(interval[1])*f(m_k)) < 0 ? "<0" : ">0" 
        verbose && push!(prod_s, prod) 
        
        n += 1
    end

    if verbose
        data = hcat(1:n, as, bs, approximations, prod_s, errors)
        header = ["S.No" "a_k-1" "b_k-1" "m_k" "f(a_k-1)*f(m_k)" "error"]
        pretty_table(data, header)
        println("$n iterations")
    end

    return approximations, errors
end

function RegulaFalsiEval(f::Function, a, b; tol = 0.005, N = 100, verbose=true)
    # make expr
    error = Base.Inf64 

    # find starting interval
    interval = (a, b)

    # iterate
    approximations = []
    errors = []
    n = 0

    if verbose
        as = []
        bs = []
        f_as = []
        f_bs = []
        prod_s =[]
    end

    while error > tol && n<=N
        # get m_k

        m_k = (interval[1]*f(interval[2]) - interval[2]*f(interval[1])) / (f(interval[2]) - f(interval[1])) 

        if verbose
            push!(as, interval[1])
            push!(bs, interval[2])
            push!(f_as, f(interval[1]))
            push!(f_bs, f(interval[2]))
        end

        # get error
        if length(approximations) > 0
            error = abs(last(approximations) - m_k)
        end
        push!(errors, error)

        # append the appproximation
        push!(approximations, m_k)

        # get new interval
        interval = (f(interval[1])*f(m_k)) < 0 ? (interval[1], m_k) : (m_k, interval[2])

        prod = (f(interval[1])*f(m_k)) < 0 ? "<0" : ">0" 
        verbose && push!(prod_s, prod) 

        n += 1
    end

    if verbose
        data = hcat(1:n, as, bs, f_as, f_bs, approximations, prod_s, errors)
        header = ["S.No" "a_k-1" "b_k-1" "f(a)" "f(b)" "m_k" "f(a_k-1)*f(m_k)" "error"]
        pretty_table(data, header)
        println("$n iterations")
    end

    return approximations, errors
end

function EulerODE(f::Function, x0, y0, x, h; verbose = true)
    x_n = x0
    y_n = y0
    xs = []
    ys = []
    n = 0
    push!(xs, x_n)
    push!(ys, y_n)
    while x_n != x
        y_n = y_n + h*f(x_n, y_n)
        x_n += h
        n += 1
        push!(xs, x_n)
        push!(ys, y_n)
    end

    if verbose
        data = hcat(0:n, xs, ys)
        header = ["S.No" "x" "y"]
        pretty_table(data, header)
        println("$n iterations")
    end
end

function ModifiedEulerODE(f::Function, x0, y0, x, h; verbose = true)
    x_n = x0
    y_n = y0
    xs = []
    ys = []
    y1s = []
    n = 0
    push!(xs, x_n)
    push!(ys, y_n)
    push!(y1s, 0)
    while x_n != x
        y_n_1 = y_n + h*f(x_n, y_n)
        y_n = y_n + h*0.5*(f(x_n, y_n) + f(x_n+h, y_n_1))
        x_n += h
        n += 1
        push!(y1s, y_n_1)
        push!(xs, x_n)
        push!(ys, y_n)
    end

    if verbose
        data = hcat(0:n, xs, y1s, ys)
        header = ["S.No" "x" "y1" "y"]
        pretty_table(data, header)
        println("$n iterations")
    end
end

end
