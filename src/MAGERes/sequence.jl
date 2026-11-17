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

const MAGERES_SEQ_N_GRAINS   = [3, 5, 10, 20, 50, 100, 200, 500, 1000]
const MAGERES_SEQ_N_DRAWS    = 8
const MAGERES_SEQ_N_ERUPT    = 20
const MAGERES_SEQ_SUB_STEPS  = 4_000

"""
    MAGEResSeqEntry

    One sample of the eruption sequence: the zircon sample (one eruption, or the merged eruptions of one output), a key
    identifying its synthetic analyses, erupted area [m², full section], area-weighted eruption temperature [°C], time
    [yr] of the last eruption before it (-Inf for the first) and whether it belongs to the selected output.
"""
struct MAGEResSeqEntry
    sample   :: MAGEResZirconSample
    key      :: UInt64
    area     :: Float64
    T        :: Float64
    t_prev   :: Float64
    selected :: Bool
end

"""
    MAGEResSeqStats

    Statistics of one sequence entry. From the erupted zircon mass (no analytical noise): mass-weighted mean, 5, 50
    and 95 % ages, standard deviation, skewness and excess kurtosis, youngest zircon age [kyr before t_ref], mean
    residence time (mean age − eruption age) [kyr], antecryst fraction (zircon older than the previous eruption) and
    inherited fraction. From the synthetic single spots: number of spots and of inherited spots, weighted mean and 1σ,
    MSWD, 5 and 95 % spot ages.
"""
struct MAGEResSeqStats
    erupt_age  :: Float64
    mean       :: Float64
    q05        :: Float64
    q50        :: Float64
    q95        :: Float64
    sd         :: Float64
    skew       :: Float64
    kurt       :: Float64
    youngest   :: Float64
    residence  :: Float64
    antecryst  :: Float64
    inherited  :: Float64
    n_spots    :: Int
    n_inh      :: Int
    wm         :: Float64
    wm_sigma   :: Float64
    mswd       :: Float64
    spot_q05   :: Float64
    spot_q95   :: Float64
end

"""
    mageres_seq_key(l::MAGEResLayer, o::MAGEResOptions)

    Key of the synthetic analyses of eruption layer `l` with the dating settings of `o`.
"""
mageres_seq_key(l::MAGEResLayer, o::MAGEResOptions) =
    hash((l.tracers, l.t, o.zircon_method, o.zircon_sigma_kyr, o.zircon_eruption_age_Ma, o.zircon_synthetic_grains))

"""
    mageres_zircon_sequence(st::MAGEResState; t_ref = st.t, t_sel = -Inf, group = "eruption", bounds = Float64[])

    The zircon-bearing eruptions of `st` as `MAGEResSeqEntry`, ages before `t_ref` [yr], those after `t_sel` [yr]
    flagged as selected. With `group = "output"` the eruptions between consecutive output times `bounds` [yr] are
    merged into one entry.
"""
function mageres_zircon_sequence(st::MAGEResState; t_ref::Float64 = st.t, t_sel::Float64 = -Inf,
                                 group::AbstractString = "eruption", bounds::Vector{Float64} = Float64[])
    te = mageres_te(st)
    (te === nothing || !te.zircon) && return MAGEResSeqEntry[]
    layers  = st.pluton.layers
    n_tr    = mageres_zrc_N_index(te, te.n_fam) + 1 + te.n_bins
    idx     = [k for (k, l) in enumerate(layers) if length(l.tracers) >= n_tr]
    samples = mageres_zircon_samples(st; t_ref = t_ref, zoned = false)
    t_all   = [l.t for l in layers]
    prev(t) = (j = findlast(<(t), t_all); j === nothing ? -Inf : t_all[j])
    entries = MAGEResSeqEntry[]
    for (s, k) in zip(samples, idx)
        l = layers[k]
        s.inherited + sum(s.mass) > 0 || continue
        push!(entries, MAGEResSeqEntry(s, mageres_seq_key(l, st.opts), l.area, l.T, prev(l.t), l.t > t_sel))
    end
    (group != "output" || isempty(bounds) || isempty(entries)) && return entries
    merged = MAGEResSeqEntry[]
    lo     = -Inf
    for hi in bounds
        g = [e for e in entries if lo < e.sample.t <= hi]
        lo = hi
        isempty(g) && continue
        A  = sum(e.area for e in g)
        T  = A > 0 ? sum(e.area * e.T for e in g) / A : NaN
        s  = mageres_combine_zircon_samples([e.sample for e in g], st.opts.zircon_synthetic_grains)
        push!(merged, MAGEResSeqEntry(s, hash([e.key for e in g]), A, T, g[1].t_prev, any(e -> e.selected, g)))
    end
    return merged
end

"""
    mageres_quantile(v::Vector{Float64}, p::Float64)

    Linearly interpolated quantile `p` of `v` (NaN when empty).
"""
function mageres_quantile(v::Vector{Float64}, p::Float64)
    isempty(v) && return NaN
    w = sort(v)
    h = 1 + (length(w) - 1) * p
    i = floor(Int, h)
    return i >= length(w) ? w[end] : w[i] + (h - i) * (w[i+1] - w[i])
end

"""
    mageres_seq_stats(e::MAGEResSeqEntry, bin_yr::Float64, t_ref::Float64)

    `MAGEResSeqStats` of entry `e` with zircon age bins of `bin_yr` [yr] and ages before `t_ref` [yr]. The erupted mass
    of each bin is spread uniformly over the part of the bin before the (last) eruption, and the moments include that
    spread exactly.
"""
function mageres_seq_stats(e::MAGEResSeqEntry, bin_yr::Float64, t_ref::Float64)
    s    = e.sample
    ea   = (t_ref - s.t) / 1e3
    t0   = [(j - 1) * bin_yr for j in eachindex(s.mass)]
    t1   = [min(j * bin_yr, s.t) for j in eachindex(s.mass)]
    ok   = [j for j in eachindex(s.mass) if s.mass[j] > 0 && t1[j] > t0[j]]
    M    = sum(s.mass[ok]; init = 0.0)
    nan  = NaN
    mu = q05 = q50 = q95 = sd = sk = ku = yo = ant = nan
    if M > 0
        m  = s.mass[ok]
        a  = [(t_ref - (t0[j] + t1[j]) / 2) / 1e3 for j in ok]
        w  = [(t1[j] - t0[j]) / 1e3 for j in ok]
        mu = sum(m .* a) / M
        d  = a .- mu
        m2 = sum(m .* (d .^ 2 .+ w .^ 2 ./ 12)) / M
        m3 = sum(m .* (d .^ 3 .+ d .* w .^ 2 ./ 4)) / M
        m4 = sum(m .* (d .^ 4 .+ d .^ 2 .* w .^ 2 ./ 2 .+ w .^ 4 ./ 80)) / M
        sd = sqrt(m2)
        sk = m2 > 0 ? m3 / m2^1.5 : nan
        ku = m2 > 0 ? m4 / m2^2 - 3 : nan
        p  = sortperm(a)
        c  = cumsum(m[p]) ./ M
        function q(f)
            i = searchsortedfirst(c, f)
            i = min(i, length(p))
            c0 = i == 1 ? 0.0 : c[i-1]
            return (a[p[i]] - w[p[i]] / 2) + w[p[i]] * (f - c0) / max(c[i] - c0, eps())
        end
        q05, q50, q95 = q(0.05), q(0.5), q(0.95)
        yo  = a[p[1]] - w[p[1]] / 2
        ant = isfinite(e.t_prev) ?
              sum(s.mass[j] * clamp((e.t_prev - t0[j]) / (t1[j] - t0[j]), 0.0, 1.0) for j in ok) / M : 0.0
    end
    tot = M + s.inherited
    g   = s.grains
    wm, ws, ms = length(g) >= 1 ? mageres_weighted_mean(g, s.sigma) : (nan, nan, nan)
    return MAGEResSeqStats(ea, mu, q05, q50, q95, sd, sk, ku, yo, mu - ea, ant, tot > 0 ? s.inherited / tot : nan,
                           length(g), s.n_inherited, wm, ws, ms, mageres_quantile(g, 0.05), mageres_quantile(g, 0.95))
end

"""
    mageres_task_foreach(f, r)

    Call `f(i)` for every `i` of `r` in spawned tasks, at most `Threads.nthreads()` at a time, and wait for all of
    them. Unlike `Threads.@threads`, it does not mark a threaded region, so it can run while a MAGEMin
    `@threads :static` minimisation of a live run is in progress.
"""
function mageres_task_foreach(f, r)
    idx = collect(r)
    isempty(idx) && return nothing
    nt  = max(1, min(Threads.nthreads(), length(idx)))
    nxt = Threads.Atomic{Int}(1)
    tasks = [Threads.@spawn begin
                 while true
                     k = Threads.atomic_add!(nxt, 1)
                     k > length(idx) && break
                     f(idx[k])
                 end
             end for _ in 1:nt]
    foreach(wait, tasks)
    return nothing
end

"""
    mageres_seq_estimates(entries::Vector{MAGEResSeqEntry}, t_ref::Float64)

    For every entry and every estimator of `MAGERES_AGE_ESTIMATORS`: error of the eruption-age estimate from all the
    synthetic single spots relative to the true (last) eruption age, and its 1σ [kyr], as a matrix of tuples
    (entries × estimators). Memoised per entry; the entries run in parallel tasks.
"""
function mageres_seq_estimates(entries::Vector{MAGEResSeqEntry}, t_ref::Float64)
    R = fill((NaN, NaN), length(entries), length(MAGERES_AGE_ESTIMATORS))
    mageres_task_foreach(eachindex(entries)) do i
        e = entries[i]
        y = e.sample.grains .- (t_ref - e.sample.t) / 1e3
        for (j, (m, _)) in enumerate(MAGERES_AGE_ESTIMATORS)
            R[i, j] = mageres_cached_estimate(hash((e.key, m, "all")), m, y, e.sample.sigma)
        end
    end
    return R
end

"""
    mageres_seq_grain_test(entries::Vector{MAGEResSeqEntry}, t_ref::Float64)

    Eruption-age estimators against the number of dated grains N (`MAGERES_SEQ_N_GRAINS`): for up to
    `MAGERES_SEQ_N_ERUPT` entries evenly spread over the sequence, `MAGERES_SEQ_N_DRAWS` random subsets of N of their
    single spots. Returns N, and per estimator the median absolute error [kyr] and the fraction of estimates within 2σ
    of the true eruption age (NaN where no entry has N spots).
"""
function mageres_seq_grain_test(entries::Vector{MAGEResSeqEntry}, t_ref::Float64)
    pool = [e for e in entries if length(e.sample.grains) >= 3]
    isempty(pool) && return Int[], zeros(0, 0), zeros(0, 0)
    pick = length(pool) <= MAGERES_SEQ_N_ERUPT ? pool :
           pool[unique(round.(Int, range(1, length(pool); length = MAGERES_SEQ_N_ERUPT)))]
    n_max = maximum(length(e.sample.grains) for e in pick)
    Ns    = [n for n in MAGERES_SEQ_N_GRAINS if n <= n_max]
    nm    = length(MAGERES_AGE_ESTIMATORS)
    err   = [Float64[] for _ in Ns, _ in 1:nm]
    cov   = [Bool[] for _ in Ns, _ in 1:nm]
    jobs  = [(e, a, d) for e in pick for a in eachindex(Ns) for d in 1:MAGERES_SEQ_N_DRAWS
             if length(e.sample.grains) >= Ns[a]]
    out   = Vector{Vector{Tuple{Float64,Float64}}}(undef, length(jobs))
    mageres_task_foreach(eachindex(jobs)) do k
        e, a, d = jobs[k]
        N   = Ns[a]
        sub = randperm(MersenneTwister(hash((e.key, N, d))), length(e.sample.grains))[1:N]
        y   = e.sample.grains[sub] .- (t_ref - e.sample.t) / 1e3
        sg  = e.sample.sigma[sub]
        out[k] = [mageres_cached_estimate(hash((e.key, m, N, d)), m, y, sg; steps = MAGERES_SEQ_SUB_STEPS, seed = d)
                  for (m, _) in MAGERES_AGE_ESTIMATORS]
    end
    for (k, (_, a, _)) in enumerate(jobs), j in 1:nm
        v, s = out[k][j]
        isfinite(v) || continue
        push!(err[a, j], abs(v))
        push!(cov[a, j], abs(v) <= 2s)
    end
    med = [isempty(err[a, j]) ? NaN : mageres_quantile(err[a, j], 0.5) for a in eachindex(Ns), j in 1:nm]
    cv  = [isempty(cov[a, j]) ? NaN : count(cov[a, j]) / length(cov[a, j]) for a in eachindex(Ns), j in 1:nm]
    return Ns, med, cv
end
