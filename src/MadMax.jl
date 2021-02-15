module MadMax

using PrettyTables

export BisectionEval
export RegulaFalsiEval

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
