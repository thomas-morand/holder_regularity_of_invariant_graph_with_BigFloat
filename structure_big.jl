"The four files: structure.jl, upper_bound.jl, lower_bound.jl, and holder_regularity.jl
allow us to compute the Hölder regularity of the probability of extinction of a Galton-Watson process
in dynamical environnments in the uniformly supercritical case.
To do this, we need to estimate two Lyapunov exponents lambda_u and lambda_f. 
We estimate a lower bound of these Lyapunov exponents in the file lower_bound.jl
and a upper bound in the file upper_bound.jl.
The file holder_regularity.jl is used to link the different files 
and plot the Holder regularity of the probability of extinction of a family of Galton-Watson process
in dynamical environnments parametrised by a real parameter.
Here is the article associated with these programmes: ##############################
"

"This file contains the structures used in the code.
The first structure contains the necessary information about the transformation,
the second about the laws of reproduction.
Finally, a mutable structure allows us to manage intervals for interval arithmetic."

### Packages

using IntervalArithmetic
using Plots
using Dates
using StatsBase
setprecision(100)

### Transformations

"The transformation T is a C^1 uniformly dilating transformation of the circle seen as [0,1].
We also ask that T(0)=0 (but this is not limiting: if T(0)!=0,
simply apply the corresponding rotation to phi)."

struct T
    "The transformation T"
    value ::Function
    "The derivate of T"
    derivate ::Function
    "The maximum of the derivate of T"
    max_derivate ::Float64
    "The topological degree of T"
    degree ::Int64
end

### Examples of transformations

T1(e ::Float64) = T(
    x -> 2*x + e*sin(2*π*x),
    x -> 2 + 2*π*e*cos(2*π*x) ,
    2 + 2*π*e,
    2
)
"The absolute value of e is less than 1/(π*2) which is approximately 0.159"

T2(n ::Int64) = T(
    x -> n*x,
    x -> n,
    Float64(n),
    n
)
"n is an integer bigger than 2"

### Laws of reproduction

"phi is a family of probability generating function.
The first variable x allow us to parametrize the family and take value in the circle.
So, s->phi(x,s) is the probability generating function of law mu_x"

struct phi
    "The law of reproduction"
    value ::Function
    "The logarithm of the derivate with respect to s of the law of reproduction"
    log_sderivate ::Function
end

### Examples of laws of reproduction

phi1(l ::Float64) = phi(
    (x, s) -> exp((s-1) * exp(l + cos(2*x*π))),
    (x, s) -> log(exp(l + cos(2*x*π)) * phi1(l).value(x, s))
)

phi2(l ::Float64) = phi(
    (x, s) -> exp((s-1) * exp(l - cos(2*x*π))),
    (x, s) -> log(exp(l - cos(2*x*π)) * phi2(l).value(x, s))
)

phi3(l ::Float64) = phi(
    (x, s) -> (atan(-cos(2 * x * π) - l) / π + 1 / 2) / 
              (1 - (1 - (atan(-cos(2 * x * π) - l) / π + 1 / 2)) * s),
    (x, s) -> log((atan(-cos(2 * x * π) - l) / π + 1 / 2) * 
                  (1 - (atan(-cos(2 * x * π) - l) / π + 1 / 2)) / 
                  (1 - (1 - (atan(-cos(2 * x * π) - l) / π + 1 / 2)) * s)^2)
)

phi4(l ::Float64, a ::Float64)= phi(
    (x, s) -> phi2(l).value(x + a, s),
    (x, s) -> phi2(l).log_sderivate(x + a, s)
)

### Intervals

"The structure str_interval allows us to manage intervals in the context of interval arithmetic."

mutable struct str_interval
    "The interval"
    inter ::Interval{BigFloat}
    "The integer N such that 2^(-N) is the length of the interval"
    length ::Int64
    "A boolean to know if the interval was obtained at the last stage"
    new ::Bool
    "The bound obtained on the interval (for the quantity we wish to bound on this interval)."
    bound ::BigFloat
end