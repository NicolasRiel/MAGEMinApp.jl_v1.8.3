#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMinApp
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Nerone, S., Dominguez, H., Moyen, J-F.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#

const MAGERES_AGE_ESTIMATORS = [("ysg",   "Youngest grain"),
                                ("yc2s",  "Youngest cluster YC2σ(3+)"),
                                ("wm",    "Weighted mean"),
                                ("bayes", "Bayesian (uniform)")]

const MAGERES_BAYES_STEPS  = 20_000
const MAGERES_BAYES_BURNIN = 5_000

"""
    mageres_weighted_mean(x::Vector{Float64}, s::Vector{Float64})

    Inverse-variance weighted mean of `x` with 1σ `s`, its 1σ uncertainty and the MSWD.
"""
function mageres_weighted_mean(x::Vector{Float64}, s::Vector{Float64})
    w  = 1 ./ s .^ 2
    mu = sum(w .* x) / sum(w)
    return mu, 1 / sqrt(sum(w)), length(x) > 1 ? sum(w .* (x .- mu) .^ 2) / (length(x) - 1) : NaN
end

"""
    mageres_erfc(x::Float64)

    Complementary error function with a fractional error below 1.2e-7 everywhere (Chebyshev fit, Press et al.,
    Numerical Recipes), accurate in the tails.
"""
function mageres_erfc(x::Float64)
    z = abs(x)
    t = 1.0 / (1.0 + 0.5z)
    r = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
            t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? r : 2.0 - r
end

"""
    mageres_normal_mass(lo::Float64, hi::Float64)

    Probability that a standard normal variable lies in [`lo`, `hi`], computed on the tail that avoids cancellation.
"""
function mageres_normal_mass(lo::Float64, hi::Float64)
    if lo > 0
        return 0.5 * (mageres_erfc(lo / sqrt(2.0)) - mageres_erfc(hi / sqrt(2.0)))
    else
        return 0.5 * (mageres_erfc(-hi / sqrt(2.0)) - mageres_erfc(-lo / sqrt(2.0)))
    end
end

"""
    mageres_youngest_grain(x, s)

    Youngest single-spot age of `x` [kyr] and its 1σ.
"""
function mageres_youngest_grain(x::Vector{Float64}, s::Vector{Float64})
    isempty(x) && return NaN, NaN
    i = argmin(x)
    return x[i], s[i]
end

"""
    mageres_youngest_cluster(x, s; n_min = 3)

    Weighted mean and 1σ of the youngest cluster of at least `n_min` ages of `x` that overlap at 2σ
    (YC2σ(3+), Dickinson & Gehrels 2009), or NaN when there is none.
"""
function mageres_youngest_cluster(x::Vector{Float64}, s::Vector{Float64}; n_min::Int = 3)
    p = sortperm(x)
    for i in eachindex(p)
        hi = x[p[i]] + 2s[p[i]]
        g  = [p[i]]
        for j in i+1:length(p)
            x[p[j]] - 2s[p[j]] <= hi || break
            hi = min(hi, x[p[j]] + 2s[p[j]])
            push!(g, p[j])
        end
        if length(g) >= n_min
            mu, se, _ = mageres_weighted_mean(x[g], s[g])
            return mu, se
        end
    end
    return NaN, NaN
end

"""
    mageres_uniform_loglik(a, b, x, s)

    Log-likelihood of ages `x` with Gaussian 1σ `s` drawn from a uniform crystallisation distribution between the
    ages `a` < `b` [kyr].
"""
function mageres_uniform_loglik(a::Float64, b::Float64, x::Vector{Float64}, s::Vector{Float64})
    L = -length(x) * log(b - a)
    @inbounds for i in eachindex(x)
        L += log(max(mageres_normal_mass((a - x[i]) / s[i], (b - x[i]) / s[i]), 1e-300))
    end
    return L
end

"""
    mageres_bayesian_eruption_age(x, s; steps, burnin, seed = 1)

    Bayesian eruption age (Keller, Schoene & Samperton 2018): Metropolis sampling of the young and old edges (a, b) of
    a uniform relative crystallisation distribution, observed with the Gaussian analytical uncertainties `s`, under a
    flat prior; the step size is tuned during burn-in. Returns the posterior mean and standard deviation of the
    young edge `a` [kyr], i.e. the end of crystallisation taken as the eruption age.
"""
function mageres_bayesian_eruption_age(x::Vector{Float64}, s::Vector{Float64}; steps::Int = MAGERES_BAYES_STEPS,
                                       burnin::Int = MAGERES_BAYES_BURNIN, seed::Integer = 1)
    n = length(x)
    n == 0 && return NaN, NaN
    n == 1 && return x[1], s[1]
    rng  = MersenneTwister(seed)
    sm   = sort(s)[(n + 1) ÷ 2]
    a, b = minimum(x), maximum(x)
    b    = max(b, a + sm)
    L    = mageres_uniform_loglik(a, b, x, s)
    h    = max((b - a) / n, sm / 2)
    acc  = 0
    sum1 = 0.0
    sum2 = 0.0
    m    = 0
    for k in 1:steps
        a1 = a + h * randn(rng)
        b1 = b + h * randn(rng)
        if b1 > a1
            L1 = mageres_uniform_loglik(a1, b1, x, s)
            if log(rand(rng)) < L1 - L
                a, b, L = a1, b1, L1
                acc += 1
            end
        end
        if k <= burnin && k % 100 == 0
            r    = acc / 100
            h   *= r < 0.2 ? 0.7 : r > 0.4 ? 1.4 : 1.0
            acc  = 0
        end
        if k > burnin
            sum1 += a
            sum2 += a * a
            m    += 1
        end
    end
    mu = sum1 / m
    return mu, sqrt(max(sum2 / m - mu^2, 0.0))
end

"""
    mageres_eruption_age_estimate(method::AbstractString, x, s; steps = MAGERES_BAYES_STEPS, seed = 1)

    Eruption age [kyr] and its 1σ from single-spot ages `x` with 1σ `s` by `method` (see `MAGERES_AGE_ESTIMATORS`).
"""
function mageres_eruption_age_estimate(method::AbstractString, x::Vector{Float64}, s::Vector{Float64};
                                       steps::Int = MAGERES_BAYES_STEPS, seed::Integer = 1)
    isempty(x) && return NaN, NaN
    method == "ysg"  && return mageres_youngest_grain(x, s)
    method == "yc2s" && return mageres_youngest_cluster(x, s)
    method == "wm"   && return mageres_weighted_mean(x, s)[1:2]
    method == "bayes" &&
        return mageres_bayesian_eruption_age(x, s; steps = steps, burnin = min(MAGERES_BAYES_BURNIN, steps ÷ 4),
                                             seed = seed)
    throw(ArgumentError("unknown eruption-age estimator $(method)"))
end

const MAGERES_AGE_CACHE      = Dict{UInt64,Tuple{Float64,Float64}}()
const MAGERES_AGE_CACHE_LOCK = ReentrantLock()
const MAGERES_AGE_CACHE_MAX  = 200_000

"""
    mageres_cached_estimate(key::UInt64, method, x, s; steps = MAGERES_BAYES_STEPS, seed = 1)

    `mageres_eruption_age_estimate` memoised by `key` (which must identify the ages, the method and the settings).
"""
function mageres_cached_estimate(key::UInt64, method::AbstractString, x::Vector{Float64}, s::Vector{Float64};
                                 steps::Int = MAGERES_BAYES_STEPS, seed::Integer = 1)
    r = lock(() -> get(MAGERES_AGE_CACHE, key, nothing), MAGERES_AGE_CACHE_LOCK)
    r === nothing || return r
    r = mageres_eruption_age_estimate(method, x, s; steps = steps, seed = seed)
    lock(MAGERES_AGE_CACHE_LOCK) do
        length(MAGERES_AGE_CACHE) >= MAGERES_AGE_CACHE_MAX && empty!(MAGERES_AGE_CACHE)
        MAGERES_AGE_CACHE[key] = r
    end
    return r
end
