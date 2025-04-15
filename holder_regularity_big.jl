"In this file, we compute and plot the Holder regularity of the probability of extinction.
We use the bound on the Lyapunov exponent lambda_u and lambda_f find thanks to the files
lower_bound.jl and upper_bound.jl."

### Files used 

include("structure_big.jl")
include("upper_bound_big.jl")
include("lower_bound_big.jl")

### Plot Holder regularity

function plot_bound_holder_reg(phi, T ::T, first_parameter ::Float64, size_step ::Float64, nb_step ::Int64, M_orbit ::Int64, iteration ::Int64, plt ::Bool)::Tuple{Vector{Float64}, Vector{Float64}}
    "This function plot the Holder regularity of the probability of exctinction of the law of reproduction
    phi and the transformation T for parameter l=first_paramter + size_step * i for i in 0::nb_step-1. 
    To do this, we compute the periodic orbit of T of size smaller than M with an error smaller than 2^-N."
    "step 1: periodic orbit"
    orbit, epsi = periodic_orbit(T, 40, M_orbit)
    "step 2: lambda_u"
    lower_lambda_u = maximum([fct_for_lower_bound_lambda_u(T, orbx, epsi) for orbx in orbit])
    upper_lambda_u = upper_bound_lambda_u(T, 6, iteration)
    "step 3: lambda_f"
    upper_lambda_f = zeros(Float64, nb_step)
    for p in 1:nb_step
        l = first_parameter + (p - 1) * size_step
        upper_lambda_f[p] = min(upper_bound_lambda_f(phi(l), T, 6, iteration, 12, 0.99, 0.01, 20), 0)
    end
    lower_lambda_f = zeros(Float64, nb_step)
    for p in 1:nb_step
        l = first_parameter + (p - 1) * size_step
        lower_lambda_f[p] = maximum([fct_for_lower_bound_lambda_f(phi(l), orbx, length(orbx), 100, epsi) for orbx in orbit])
    end
    "step 4: plot"
    lower_regularity = upper_lambda_f ./ (-upper_lambda_u)
    upper_regularity = lower_lambda_f ./ (-lower_lambda_u)
    if plt == true
        L = [first_parameter + i * size_step for i in 0:nb_step-1]
        p = plot(L .- size_step, upper_regularity, color = :blue, seriestype = :step)
        plot!(p, L, lower_regularity, color = :red, seriestype = :step)
        display(p)
    end
   return upper_regularity, lower_regularity
end 