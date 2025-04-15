"In this file, we compute a lower bound on the Lyapunov exponents lambda_u and lambda_f.
To do this, we estimate the periodic orbits of the transformation T and
we compute a lower bound on the Lyapunov exponents thanks to the estimation of the 
periodic orbits."

### Files used 

include("structure_big.jl")

### Periodic orbits

"The next four fonctions allow us to find a representative of each sequence of size n
with integer values between 0 and d up to rotation and not containing a repeating pattern.
These sequences allow us to estimate the periodic orbits of the T transformation."

function is_rotation_of_minimal(seq ::Vector{Int64}, n ::Int64) ::Bool
    for i in 1:n-1
        rotated_seq = vcat(seq[i+1:end], seq[1:i])
        if rotated_seq < seq
            return false
        end
    end
    return true
end

function has_repeating_pattern(seq ::Vector{Int64}, n ::Int64) ::Bool
    for len in 1:div(n, 2)
        if n % len == 0
            pattern = seq[1:len]
            is_repeating = true
            for i in 1:div(n, len)
                if seq[(i-1)*len+1:i*len] != pattern
                    is_repeating = false
                    break
                end
            end
            if is_repeating
                return true
            end
        end
    end
    return false
end

function generate_minimal_non_repeating_sequences(d ::Int64, n ::Int64) ::Vector{Vector{Int64}}
    minimal_sequences = Vector{Vector{Int64}}()
    seq = fill(0, n)
    while seq != nothing
        if is_rotation_of_minimal(seq, n) && !has_repeating_pattern(seq, n)
            push!(minimal_sequences, copy(seq)) 
        end
        seq = next_sequence(seq, d)
    end
    return minimal_sequences
end

function next_sequence(seq ::Vector{Int64}, d ::Int64) ::Union{Nothing,Vector{Int64}}
    n = length(seq)
    for i in n:-1:1
        if seq[i] < d - 1
            seq[i] += 1
            for j in i+1:n
                seq[j] = 0
            end
            return seq
        end
    end
    return nothing
end

function dichotomy(g ::Function, z ::Float64, x_min ::Float64, x_max ::Float64, epsi ::Float64) ::Float64
    "This function allow us to an approximation of the equation g(x)=z with an error less than epsi
    in the interval [x_min,x_max] thanks to the method of dichotomy."
    t = (x_min + x_max) / 2
    while x_max - x_min > epsi
        if g(t) > z
            x_max = t
        else
            x_min = t
        end
        t = (x_min + x_max) / 2
    end
    return t
end

function d_adique_to_real(suite ::Vector{Int64}, d ::Int64, i ::Int64) ::Vector{Float64}
    "This function takes as input a sequence of size i of integers between 0 and d-1
    and returns the real number between 0 and 1 having this decimal developement in base d."
    value = [sum(suite[(j+i-k) % i+1] * d^(j-1) for j in 1:i) for k in 1:i]
    return value
end

function compose(f ::Function, x ::Union{Float64,Interval{Float64}}, n ::Int64) ::Union{Float64,Interval{Float64}}
    "This function compute the value of the n-th composition of a function f"
    if n == 0
        return x
    end
    return compose(f, f(x), n-1)
end

function periodic_orbit(T ::T, N ::Int64, M ::Int64) ::Tuple{Vector{Vector{Float64}},Float64}
    "This function allow us to find an approximation with an error smaller than 2^-N
    of all periodic orbit of the transformation T of lenght smaller than M."
    d = T.degree
    epsi = 1 / (2^N)
    seq_orbite = Vector{Vector{Float64}}()
    push!(seq_orbite, [0.0])
    for i in 1:M
        seq = Vector{Vector{Int64}}()
        g = x -> compose(T.value, x, i) - x
        if i > 1
            seq = generate_minimal_non_repeating_sequences(T.degree, i)
        end
        if i == 1 && d > 2
            seq = [[j] for j in 1:d-2]
        end
        for x in seq 
            image_orb_x = d_adique_to_real(x, T.degree, i)
            orb_x = Vector{Float64}()
            xm = 0.0
            xM = 1.0
            for f in image_orb_x
                a = dichotomy(g, Float64(f), xm, xM, epsi)
                push!(orb_x, a)
                Ta = T.value(a)
                err = epsi * T.max_derivate
                xm = Ta - err
                xM = Ta + err
            end
            push!(seq_orbite, orb_x) 
        end
    end
    return seq_orbite, epsi
end

### Lower bound: lambda_f

function phi_n_per(phi, orb ::Vector{Float64}, len ::Int64, s ::Float64, n ::Int64, i ::Int64, err ::Float64) ::Interval{Float64}
    "Compute the value of phi^{(n)} along an estimate of a periodic orbit orb
    of lenght len and with an error smaller than err"
    if n == 0
        return interval(s, s)
    end
    return phi.value(interval(orb[i+1] - err, orb[i+1] + err), phi_n_per(phi, orb, len, s, n-1, (i+1) % len, err))
end


function fct_for_lower_bound_lambda_f(phi ::phi, orbx ::Vector{Float64}, i ::Int64, n ::Int64, epsi ::Float64) ::Float64
    "This function is the one whose Birkhoff sum along a periodic orbit can be used to bound lambda_f."
    bound = 0.
    q_Tx = phi_n_per(phi, orbx, i, 0., n, 0, epsi)
    for j in i:-1:1
        bound += inf(phi.log_sderivate(interval(orbx[j] - epsi, orbx[j] + epsi), q_Tx))
        q_Tx = phi.value(interval(orbx[j] - epsi, orbx[j] + epsi), q_Tx)
    end
    return bound / i
end

### Lower bound: lambda_u

function fct_for_lower_bound_lambda_u(T ::T, orb ::Vector{Float64}, epsi ::Float64) ::Float64
    "This function is the one whose Birkhoff sum along a periodic orbit can be used to bound lambda_u."
    bound = interval(0, 0)
    i = length(orb)
    for j in 1:i
        bound += log(T.derivate(interval(orb[j] - epsi, orb[j] + epsi)))
    end
    return inf(bound) / i
end