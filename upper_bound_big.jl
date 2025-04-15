"In this file, we compute a upper bound on the Lyapunov exponent lambda_u and lambda_f.
To do this, we use arithmetic intervals.
We compute the target quantity is estimated over small intervals.
In particular, we refine the study by reducing the size of the intervals
where the bound obtained is the worst."

### Files used 

include("structure_big.jl")

### Upper bound: general method

function split_interval(I ::str_interval) ::Tuple{str_interval, str_interval}
    "This function splits an interval of type str_interval into two new intervals of equal size."
    mid = (inf(I.inter) + sup(I.inter)) / 2
    I1 = str_interval(interval(inf(I.inter), mid), I.length+1, true, I.bound) 
    I2 = str_interval(interval(mid, sup(I.inter)), I.length+1, true, I.bound)
    return I1,I2
end

function split_interval(I ::Interval{Float64}) ::Tuple{Interval{Float64},Interval{Float64}}
    "We define the function split_interval for interval of type Interval"
    mid = (inf(I) + sup(I)) / 2
    I1 = interval(inf(I), mid) 
    I2 = interval(mid, sup(I))
    return I1,I2
end

function split_worst_intervals(l_input ::Vector{str_interval}, q ::Float64) ::Vector{str_interval}
    "This function identifies the intervals with the worst bound
    (a proportion of intervals equal to 1-q in the list of intervals l_input).
    It then divides these intervals in two and returns the list of updated intervals l_output."
    bound = quantile([I.bound for I in l_input], q)
    len_output = length(l_input) + count(I -> I.bound > bound, l_input)
    l_output = Vector{str_interval}(undef, len_output)
    j = 1
    for I in l_input
        if I.bound > bound
            I1, I2 = split_interval(I)
            l_output[j] = I1
            j += 1
            l_output[j] = I2
        else
            l_output[j] = I
        end
        j += 1
    end
    return l_output
end

function compute_of_bounds(T ::T, l ::Vector{str_interval}, g ::Function, K=nothing ::Union{Float64,Nothing}, phi=nothing ::Union{phi,Nothing}) ::Vector{str_interval}
    "Udpdate the upper bound by computing the Birkhoff sum of the function g
    on the new intervals of the list of intervals l"
    for I in filter(I -> I.new, l)
        I.new = false
        Birkhoff_sum = interval(0, 0)
        m = div(I.length, 2) + 1
        inter = I.inter
        for j in 1:m
            Birkhoff_sum += g(inter, T, m, K, phi)
            I.bound = min(sup(Birkhoff_sum) / j, I.bound)
            inter = T.value(inter)
        end
    end
    return l
end

### Upper bound: lambda_f

function phi_n(phi ::phi, T ::T, x ::Union{Float64,Interval{Float64}}, s ::Union{Float64,Interval{Float64}}, n ::Int64) ::Union{Float64,Interval{Float64}}
    "Compute the value of phi^{(n)}."
    if n == 0
        return s 
    end
    return phi.value(x, phi_n(phi, T, T.value(x), s, n-1))
end

function min_btw_phi_n_K(phi ::phi, T ::T, x ::Union{BigFloat,Interval{BigFloat}}, s ::Union{BigFloat,Interval{BigFloat}}, n ::Int64, K ::BigFloat) ::Union{BigFloat,Interval{BigFloat}}
    "Compute recursevely the value of the minimum between phi^{(n)} and K."
    if n == 0
        return s
    end
    return phi.value(x, min(min_btw_phi_n_K(phi, T, T.value(x), s, n-1, BigFloat(K)), BigFloat(K)))
end

function verify_condition_on_K(phi ::phi, T ::T, l_input ::Vector{Interval{Float64}}, K ::Float64, N_phi ::Int64, split ::Bool)
    "For all interval I in the list l_input, verify the condition on K. 
    If the condition is true for I: add the subdision in two egal part of I to l_ouput if split is true,
    and add I to l_output if split is false" 
    l_output = Vector{Interval{Float64}}()
    for I in l_input
        if all(x -> x > K, [sup(phi_n(phi, T, I, K, n)) for n in 1:N_phi])
            if split
                I1, I2 = split_interval(I)
                append!(l_output, [I1,I2])
            else
                push!(l_output, I)
            end
        end
    end
    return l_output
end

function find_K(phi ::phi, T ::T, N ::Int64, coeff ::Float64, K_max ::Float64, iteration ::Int64) ::Float64
    "We find a constant K (as small as possible) which dominates the function q.
    To do that, we estime the probability of extinction and therefore its maximum.
    Thus, we tested if K is suitable for increasingly large values of K"
    l_estime_K = [(i-1)/2^N for i in 1:2^N]
    bound = maximum([phi_n(phi, T, x, 0., N) for x in l_estime_K])
    K = bound + coeff * (1 - bound)
    l = [interval(0., 1.)]
    for i in 1:iteration 
        l = verify_condition_on_K(phi, T, l, K, div(i, 2) + 1, true)
        if isempty(l)
            return K
        end
    end
    while K < K_max
        K = K + (1 - K) * coeff
        l = verify_condition_on_K(phi, T, l, K, div(iteration, 2) + 1, false)
        if isempty(l)
            return K
        end
    end 
    error("K has not been found")
end

function fct_for_upper_bound_lambda_f(inter ::Interval{BigFloat}, T ::T, m ::Int64, K ::BigFloat, phi ::phi) ::Interval{BigFloat}
    "This function is the one whose Birkhoff sum can be used to bound lambda_f."
    return phi.log_sderivate(inter, min_btw_phi_n_K(phi, T, T.value(inter), interval(K, K), m, K))
end

function upper_bound_lambda_f(phi ::phi, T ::T, N ::Int64, iteration ::Int64, N_find_K ::Int64, K_max ::Float64, coeff ::Float64, iteration_find_K ::Int64) ::Float64
    "This function allow us to find a good upper bound of lambda_f"
    l = [str_interval(interval(BigFloat((i-1)/2^N), BigFloat(i/2^N)), N, true, BigFloat(1000.0)) for i in 1:2^N]
    K = BigFloat(find_K(phi, T, N_find_K, coeff, K_max, iteration_find_K))
    for _ in 1:iteration
        l = compute_of_bounds(T, l, fct_for_upper_bound_lambda_f, K, phi)
        l = split_worst_intervals(l, 0.95)
    end
    return maximum([I.bound for I in l])
end

### Upper bound: lambda_u

function fct_for_upper_bound_lambda_u(inter ::Interval{BigFloat}, T ::T, m, K, phi) ::Interval{BigFloat}
    "This function is the one whose Birkhoff sum can be used to bound lambda_f."
    return log(T.derivate(inter))
end

function upper_bound_lambda_u(T ::T, N ::Int64, iteration ::Int64) ::Float64
    "This function allow us to find a good upper bound of lambda_u"
    l = [str_interval(interval((i-1)/2^N,i/2^N), N, true, 1000.0) for i in 1:2^N]
    for _ in 1:iteration 
        l = compute_of_bounds(T, l, fct_for_upper_bound_lambda_u)
        l = split_worst_intervals(l, 0.95)
    end
    return maximum([I.bound for I in l])
end
