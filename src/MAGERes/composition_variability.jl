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

const MAGERES_COMP_VAR_LAWS   = ("scaled", "constant")
const MAGERES_COMP_VAR_REF_WF = 0.01

"""
    mageres_comp_sigma(o::MAGEResOptions, w::Real)

    Relative standard deviation (fraction) of a component of mass fraction `w`: `comp_var_ref_pct` (w / 1 wt%)^(−b)
    with b = `comp_var_exponent` for the "scaled" law, `comp_var_ref_pct` for the "constant" law, at most
    `comp_var_max_pct`.
"""
function mageres_comp_sigma(o::MAGEResOptions, w::Real)
    cap = o.comp_var_max_pct / 100
    o.comp_var_law == "constant" && return min(o.comp_var_ref_pct / 100, cap)
    w > 0 || return cap
    return min(o.comp_var_ref_pct / 100 * (w / MAGERES_COMP_VAR_REF_WF)^(-o.comp_var_exponent), cap)
end

"""
    mageres_oxide_mass_fractions(oxides::Vector{String}, x::AbstractVector{<:Real})

    Mass fractions of the molar bulk `x` over `oxides`.
"""
function mageres_oxide_mass_fractions(oxides::Vector{String}, x::AbstractVector{<:Real})
    m = [x[k] * MAGERES_OXIDE_MOLAR_MASS[oxides[k]] for k in eachindex(oxides)]
    s = sum(m)
    return s > 0 ? m ./ s : m
end

"""
    mageres_comp_sigmas(o::MAGEResOptions, oxides::Vector{String}, x, C)

    Relative standard deviations of the molar bulk `x` over `oxides` (from its mass fractions) and of the trace
    elements `C` [µg/g].
"""
mageres_comp_sigmas(o::MAGEResOptions, oxides::Vector{String}, x::AbstractVector{<:Real}, C::AbstractVector{<:Real}) =
    ([mageres_comp_sigma(o, w) for w in mageres_oxide_mass_fractions(oxides, x)],
     [mageres_comp_sigma(o, 1e-6 * c) for c in C])

"""
    mageres_lognormal_factor(σ::Real, ξ::Real)

    Multiplicative factor of mean 1 and relative standard deviation `σ` for the standard normal deviate `ξ`.
"""
function mageres_lognormal_factor(σ::Real, ξ::Real)
    s = sqrt(log1p(σ^2))
    return exp(s * ξ - s^2 / 2)
end

"""
    mageres_perturb(x, σ, ξ; closure)

    `x` with every entry multiplied by its lognormal factor (`σ`, `ξ`); zero entries stay zero. With `closure`, the
    result is rescaled to the sum of `x`.
"""
function mageres_perturb(x::AbstractVector{<:Real}, σ::AbstractVector{<:Real}, ξ::AbstractVector{<:Real};
                         closure::Bool)
    y = [x[k] * mageres_lognormal_factor(σ[k], ξ[k]) for k in eachindex(x)]
    (closure && sum(y) > 0) && (y .*= sum(x) / sum(y))
    return y
end

"""
    mageres_comp_rng(o::MAGEResOptions, stream::Integer, k::Integer, part::Integer)

    Random generator of draw `k` of `stream` (1 injection, 2 host) and `part` (0 major oxides, 1 trace elements),
    seeded from `comp_var_seed`, so that each draw is independent of the others and of the run length.
"""
mageres_comp_rng(o::MAGEResOptions, stream::Integer, k::Integer, part::Integer) =
    MersenneTwister(UInt32[o.comp_var_seed, stream, k, part])

"""
    mageres_injection_base_bulk(o::MAGEResOptions)

    Normalised molar bulk of the injected magma before variability.
"""
mageres_injection_base_bulk(o::MAGEResOptions) =
    mageres_normalized(isempty(o.injection_bulk) ? o.magma_bulk : o.injection_bulk)

"""
    mageres_injection_composition(o::MAGEResOptions, k::Integer, C::AbstractVector{<:Real})

    Normalised molar bulk and trace elements [µg/g] (perturbed from `C`) of injected batch `k`. Without
    `comp_var_injection` the base composition is returned.
"""
function mageres_injection_composition(o::MAGEResOptions, k::Integer, C::AbstractVector{<:Real})
    x = mageres_injection_base_bulk(o)
    o.comp_var_injection || return x, collect(Float64, C)
    σx, σc = mageres_comp_sigmas(o, o.oxides, x, C)
    xk = mageres_perturb(x, σx, randn(mageres_comp_rng(o, 1, k, 0), length(x)); closure = true)
    Ck = mageres_perturb(C, σc, randn(mageres_comp_rng(o, 1, k, 1), length(C)); closure = false)
    return xk, Ck
end

"""
    mageres_host_profile_y(o::MAGEResOptions)

    Heights [m] of the host-profile nodes: spaced for half the re-equilibration tolerance along the geotherm, and at
    most a quarter of `comp_var_host_corr_m` apart with host variability.
"""
function mageres_host_profile_y(o::MAGEResOptions)
    H  = 1e3 * o.domain_thickness_km
    dT = o.geotherm_gradient_C_km * o.domain_thickness_km
    nn = max(2, ceil(Int, dT / (0.5 * o.dT_tol)) + 1)
    o.comp_var_host && (nn = max(nn, ceil(Int, 4H / o.comp_var_host_corr_m) + 1))
    return collect(range(-H, 0.0; length = nn))
end

"""
    mageres_depth_deviates(rng, n::Int, y::Vector{Float64}, L::Float64)

    `n` × `length(y)` standard normal deviates, each row an Ornstein–Uhlenbeck sequence along the heights `y` with
    correlation length `L`.
"""
function mageres_depth_deviates(rng, n::Int, y::Vector{Float64}, L::Float64)
    Ξ = zeros(n, length(y))
    isempty(y) && return Ξ
    Ξ[:, 1] .= randn(rng, n)
    for j in 2:length(y)
        r = exp(-abs(y[j] - y[j-1]) / L)
        Ξ[:, j] .= r .* Ξ[:, j-1] .+ sqrt(1 - r^2) .* randn(rng, n)
    end
    return Ξ
end

"""
    mageres_host_compositions(o::MAGEResOptions, y::Vector{Float64}, C::AbstractVector{<:Real})

    Normalised molar host bulks and host trace elements [µg/g] (perturbed from `C`) at the heights `y`. With
    `comp_var_host` every component varies along `y` with correlation length `comp_var_host_corr_m`; otherwise the
    base composition is returned at every height.
"""
function mageres_host_compositions(o::MAGEResOptions, y::Vector{Float64}, C::AbstractVector{<:Real})
    x0 = mageres_normalized(o.host_bulk)
    C0 = collect(Float64, C)
    o.comp_var_host || return [copy(x0) for _ in y], [copy(C0) for _ in y]
    σx, σc = mageres_comp_sigmas(o, o.oxides, x0, C0)
    L  = o.comp_var_host_corr_m
    Ξx = mageres_depth_deviates(mageres_comp_rng(o, 2, 0, 0), length(x0), y, L)
    Ξc = mageres_depth_deviates(mageres_comp_rng(o, 2, 0, 1), length(C0), y, L)
    return [mageres_perturb(x0, σx, view(Ξx, :, j); closure = true) for j in eachindex(y)],
           [mageres_perturb(C0, σc, view(Ξc, :, j); closure = false) for j in eachindex(y)]
end

"""
    mageres_host_bulk_at(o::MAGEResOptions, y::Real)

    Normalised molar host bulk at height `y` [m], linearly interpolated between the host-profile nodes.
"""
function mageres_host_bulk_at(o::MAGEResOptions, y::Real)
    o.comp_var_host || return mageres_normalized(o.host_bulk)
    ys    = mageres_host_profile_y(o)
    xs, _ = mageres_host_compositions(o, ys, Float64[])
    yy    = clamp(Float64(y), ys[1], ys[end])
    j     = clamp(searchsortedlast(ys, yy), 1, length(ys) - 1)
    w     = (yy - ys[j]) / (ys[j+1] - ys[j])
    return (1 - w) .* xs[j] .+ w .* xs[j+1]
end
