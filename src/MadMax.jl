module MadMax

using PrettyTables
using Flux

export BisectionEval
export RegulaFalsiEval
export NewtonEval

function NewtonEval(f::Function, x0; tol = 0.005, N=100, verbose=true)
    error = Base.Inf64
    # find gradient 
    df(x) = gradient(f, x)[1]

    approximations = [x0]
    errors = [error]
    n = 0

    # iterate
    while error > tol && n != N
        x_last = last(approximations)
        delta = f(x_last) / df(x_last)
        x_n = x_last - delta

        push!(approximations, x_n)

        # get error
        error = error(x_last - x_n)
        push!(errors, error)

        n += 1
    end

    if verbose
        headers = ["S.N" "x_n" "Error"]
        data = hcat(0:n, approximations, errors)
        
        pretty_table(data, header)
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

end
