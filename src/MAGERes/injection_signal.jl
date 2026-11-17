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

const MAGERES_T_INJ_MODES  = ("constant", "signal")
const MAGERES_NOISE_MIN_STEP_YR = 1.0

"""
    mageres_injection_signal_on(o::MAGEResOptions)

    Whether the injection temperature follows the time signal of `o`.
"""
mageres_injection_signal_on(o::MAGEResOptions) = o.T_injection_mode == "signal"

"""
    mageres_injection_peak(o::MAGEResOptions, t::Real)

    First-order term [°C] at time `t` [yr]: Gaussian heating-then-cooling pulses of amplitude `T_inj_peak_C` above the
    base temperature (negative for a cooling trough), with standard deviation `T_inj_peak_rise_yr` before each peak and
    `T_inj_peak_decay_yr` after it (equal widths give a symmetric Gaussian). There are `T_inj_peak_cycles` peaks, at
    `T_inj_peak_yr` + k `T_inj_peak_interval_yr` (k = 0, 1, …); where pulses overlap the strongest one applies, so
    every cycle reaches the same amplitude.
"""
function mageres_injection_peak(o::MAGEResOptions, t::Real)
    g = 0.0
    for k in 0:max(o.T_inj_peak_cycles, 1)-1
        tp = o.T_inj_peak_yr + k * o.T_inj_peak_interval_yr
        s  = t < tp ? o.T_inj_peak_rise_yr : o.T_inj_peak_decay_yr
        g  = max(g, exp(-0.5 * ((t - tp) / s)^2))
    end
    return o.T_inj_peak_C * g
end

"""
    mageres_injection_sine(o::MAGEResOptions, t::Real)

    Second-order term [°C] at time `t` [yr]: `T_inj_sine_C` sin(2π t / `T_inj_sine_period_yr` + `T_inj_sine_phase_deg`).
"""
mageres_injection_sine(o::MAGEResOptions, t::Real) =
    o.T_inj_sine_C * sin(2π * t / o.T_inj_sine_period_yr + deg2rad(o.T_inj_sine_phase_deg))

"""
    mageres_injection_noise(o::MAGEResOptions, ts::AbstractVector{<:Real})

    Third-order term [°C] at the times `ts` [yr]: an Ornstein–Uhlenbeck process of standard deviation `T_inj_noise_C`
    and correlation time `T_inj_noise_corr_yr`, drawn from `T_inj_noise_seed` on a grid starting at t = 0 with a step
    of a quarter of the correlation time (at least `MAGERES_NOISE_MIN_STEP_YR`), linearly interpolated. The values
    at a given time do not depend on the other requested times or on the run time.
"""
function mageres_injection_noise(o::MAGEResOptions, ts::AbstractVector{<:Real})
    out = zeros(length(ts))
    (o.T_inj_noise_C > 0 && !isempty(ts)) || return out
    h   = max(o.T_inj_noise_corr_yr / 4, MAGERES_NOISE_MIN_STEP_YR)
    n   = max(ceil(Int, max(maximum(ts), 0.0) / h) + 2, 2)
    rng = MersenneTwister(o.T_inj_noise_seed)
    r   = exp(-h / o.T_inj_noise_corr_yr)
    x   = Vector{Float64}(undef, n)
    x[1] = o.T_inj_noise_C * randn(rng)
    for k in 2:n
        x[k] = r * x[k-1] + o.T_inj_noise_C * sqrt(1 - r^2) * randn(rng)
    end
    for (i, t) in enumerate(ts)
        u = max(t, 0.0) / h
        k = floor(Int, u)
        w = u - k
        out[i] = (1 - w) * x[k+1] + w * x[k+2]
    end
    return out
end

"""
    mageres_injection_T_components(o::MAGEResOptions, ts::AbstractVector{<:Real})

    Injection temperature [°C] at the times `ts` [yr] and its parts: `base` (`T_injection_C`), `peak`, `sine` and
    `noise` terms, the `total` clamped to [`T_inj_min_C`, `T_inj_max_C`], and `clipped` where the clamp acted. In
    constant mode the total is `T_injection_C` and the other terms are zero.
"""
function mageres_injection_T_components(o::MAGEResOptions, ts::AbstractVector{<:Real})
    n = length(ts)
    if !mageres_injection_signal_on(o)
        z = zeros(n)
        return (base = fill(o.T_injection_C, n), peak = z, sine = copy(z), noise = copy(z),
                total = fill(o.T_injection_C, n), clipped = falses(n))
    end
    peak  = [mageres_injection_peak(o, t) for t in ts]
    sine  = [mageres_injection_sine(o, t) for t in ts]
    noise = mageres_injection_noise(o, ts)
    raw   = o.T_injection_C .+ peak .+ sine .+ noise
    total = clamp.(raw, o.T_inj_min_C, o.T_inj_max_C)
    return (base = fill(o.T_injection_C, n), peak = peak, sine = sine, noise = noise, total = total,
            clipped = BitVector(total .!= raw))
end

"""
    mageres_injection_T(o::MAGEResOptions, t::Real)

    Injection temperature [°C] at time `t` [yr].
"""
mageres_injection_T(o::MAGEResOptions, t::Real) =
    mageres_injection_signal_on(o) ? mageres_injection_T_components(o, [Float64(t)]).total[1] : o.T_injection_C

"""
    mageres_injection_times(o::MAGEResOptions)

    Times [yr] of the injections up to the run time.
"""
mageres_injection_times(o::MAGEResOptions) =
    o.injection_period_yr > 0 ? collect(o.first_injection_yr:o.injection_period_yr:o.run_time_yr) : Float64[]

"""
    mageres_injection_T_max(o::MAGEResOptions)

    Highest injection temperature [°C] over the injections of the run (`T_injection_C` without injection).
"""
function mageres_injection_T_max(o::MAGEResOptions)
    ts = mageres_injection_times(o)
    isempty(ts) && return o.T_injection_C
    return maximum(mageres_injection_T_components(o, ts).total)
end
