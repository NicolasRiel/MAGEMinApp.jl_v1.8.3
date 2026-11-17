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

const MAGERES_SEQ_GROUP_OPTIONS = [Dict("label" => "One sample per eruption",            "value" => "eruption"),
                                   Dict("label" => "One sample per output (merged)",     "value" => "output")]

const MAGERES_SEQ_MAX_VIOLINS   = 40
const MAGERES_SEQ_COLOR         = "rgb(0,114,189)"
const MAGERES_SEQ_SEL_COLOR     = "rgb(200,40,40)"
const MAGERES_INJ_T_COLOR      = "rgb(217,95,2)"
const MAGERES_ESTIMATOR_COLORS  = ["rgb(0,114,189)", "rgb(230,140,30)", "rgb(60,160,60)", "rgb(126,47,142)"]

"""
    MAGEResSequence

    Everything the sequence figures and table need: the entries, their statistics, the eruption-age estimates
    (entries × estimators, error and 1σ [kyr]), the reference time [yr] and the dating options.
"""
struct MAGEResSequence
    entries :: Vector{MAGEResSeqEntry}
    stats   :: Vector{MAGEResSeqStats}
    est     :: Matrix{Tuple{Float64,Float64}}
    t_ref   :: Float64
    opts    :: MAGEResOptions
end

"""
    mageres_sequence(st::MAGEResState; t_ref = st.t, t_sel = -Inf, group = "eruption", bounds = Float64[])

    `MAGEResSequence` of the eruptions of `st` (see `mageres_zircon_sequence`), or `nothing` without zircon.
"""
function mageres_sequence(st::MAGEResState; t_ref::Float64 = st.t, t_sel::Float64 = -Inf,
                          group::AbstractString = "eruption", bounds::Vector{Float64} = Float64[])
    te = mageres_te(st)
    (te === nothing || !te.zircon) && return nothing
    entries = mageres_zircon_sequence(st; t_ref = t_ref, t_sel = t_sel, group = group, bounds = bounds)
    stats   = [mageres_seq_stats(e, te.bin_yr, t_ref) for e in entries]
    return MAGEResSequence(entries, stats, mageres_seq_estimates(entries, t_ref), t_ref, st.opts)
end

"""
    mageres_sequence_empty(seq)

    Placeholder figure when `seq` has no zircon or no zircon-bearing eruption, otherwise `nothing`.
"""
function mageres_sequence_empty(seq)
    seq === nothing &&
        return mageres_placeholder_figure("Zircon requires trace elements with a Zr saturation model"; height = 360)
    isempty(seq.entries) && return mageres_placeholder_figure("No eruption carrying zircon yet"; height = 360)
    return nothing
end

"""
    mageres_seq_axis(title::AbstractString; kw...)

    Framed axis with outside ticks, extra attributes `kw`.
"""
mageres_seq_axis(title::AbstractString; kw...) =
    attr(; title = title, showline = true, linecolor = "black", mirror = true, ticks = "outside", zeroline = false,
           kw...)

"""
    mageres_seq_eruption_axis(seq::MAGEResSequence; kw...)

    Eruption-age axis [kyr before t_ref], reversed so that time runs from left to right.
"""
mageres_seq_eruption_axis(seq::MAGEResSequence; kw...) =
    mageres_seq_axis(@sprintf("Eruption age before t = %.1f kyr (kyr)", seq.t_ref / 1e3); autorange = "reversed", kw...)

"""
    mageres_seq_label(e::MAGEResSeqEntry)

    Hover label of entry `e`.
"""
mageres_seq_label(e::MAGEResSeqEntry) =
    mageres_eruption_label(e.sample.t_first, e.sample.t, e.sample.n_eruptions) *
    @sprintf(" — %d spots, %d inherited", length(e.sample.grains), e.sample.n_inherited)

"""
    mageres_seq_band!(traces, x, lo, hi, name, fill; legend = true)

    Append to `traces` a filled band between the curves `lo` and `hi` over `x`.
"""
function mageres_seq_band!(traces, x, lo, hi, name::AbstractString, fill::AbstractString; legend::Bool = true)
    push!(traces, scatter(; x = x, y = hi, mode = "lines", line = attr(width = 0), legendgroup = name,
                            showlegend = false, hoverinfo = "skip"))
    push!(traces, scatter(; x = x, y = lo, mode = "lines", line = attr(width = 0), fill = "tonexty",
                            fillcolor = fill, name = name, legendgroup = name, showlegend = legend,
                            hoverinfo = "skip"))
    return traces
end

"""
    mageres_injection_T_trace(T_injection, t_ref, x_lo, x_hi; yaxis, legend = nothing)

    Line of the injection temperature `T_injection = (t [yr], T [°C])` against age before `t_ref` [yr], kept within the
    age range [`x_lo`, `x_hi`] [kyr], on axis `yaxis`; `nothing` without data in that range.
"""
function mageres_injection_T_trace(T_injection, t_ref::Float64, x_lo::Float64, x_hi::Float64; yaxis::AbstractString,
                              legend = nothing)
    T_injection === nothing && return nothing
    a = (t_ref .- T_injection[1]) ./ 1e3
    k = [i for i in eachindex(a) if x_lo <= a[i] <= x_hi && isfinite(T_injection[2][i])]
    length(k) >= 2 || return nothing
    kw = legend === nothing ? (;) : (; legend = legend)
    return scatter(; x = a[k], y = T_injection[2][k], mode = "lines+markers", name = "injection T (°C, right axis)",
                     xaxis = "x", yaxis = yaxis, line = attr(color = MAGERES_INJ_T_COLOR, width = 1), marker = attr(color = MAGERES_INJ_T_COLOR, size = 3),
                     hovertemplate = "%{y:.0f} °C<extra>injection T</extra>", kw...)
end

"""
    mageres_injection_T_axis(; kw...)

    Right-hand axis of the injection temperature.
"""
mageres_injection_T_axis(; kw...) =
    attr(; title = attr(text = "Injection T (°C)", font = attr(color = MAGERES_INJ_T_COLOR)), side = "right",
           showgrid = false, zeroline = false, ticks = "outside", tickfont = attr(color = MAGERES_INJ_T_COLOR),
           linecolor = MAGERES_INJ_T_COLOR, showline = true, kw...)

"""
    mageres_zircon_sequence_figure(seq; T_injection = nothing)

    Zircon age against eruption age, ages before t_ref, x axis limited to the eruption ages. Up to
    `MAGERES_SEQ_MAX_VIOLINS` samples: one violin of the synthetic single-spot ages per sample, with the spots, all of
    the same width. Above: quantile bands of the spot ages (5–95 % and 25–75 %) and their median over the sequence,
    with the spots of the selected output. Both show the median and 5–95 % ages of the erupted zircon mass (no
    analytical noise) and the line where zircon age equals eruption age; the distance above that line is the time
    between crystallisation and eruption. Samples of the selected output are drawn in red. With `T_injection`
    (times [yr], temperatures [°C]), the temperature of every injection is drawn on the right axis.
"""
function mageres_zircon_sequence_figure(seq; T_injection = nothing)
    p = mageres_sequence_empty(seq)
    p === nothing || return p
    o      = sortperm([s.erupt_age for s in seq.stats]; rev = true)
    E      = seq.entries[o]
    S      = seq.stats[o]
    xe     = [s.erupt_age for s in S]
    span   = max(maximum(xe) - minimum(xe), 1.0)
    x_lo   = minimum(xe) - 0.03 * span
    x_hi   = maximum(xe) + 0.03 * span
    gaps   = [g for g in -diff(xe) if g > 0]
    width  = length(E) == 1 ? 0.1 * span :
             max(0.7 * (isempty(gaps) ? span : mageres_quantile(gaps, 0.5)), 0.01 * span)
    bands  = length(E) > MAGERES_SEQ_MAX_VIOLINS
    hov    = [mageres_seq_label(e) for e in E]
    traces = GenericTrace[]
    rgba(c, a) = replace(c, "rgb" => "rgba", ")" => ",$(a))")
    if bands
        k  = [i for i in eachindex(E) if length(E[i].sample.grains) >= 2]
        q(f) = [mageres_quantile(E[i].sample.grains, f) for i in k]
        mageres_seq_band!(traces, xe[k], q(0.05), q(0.95), "spots 5–95 %", rgba(MAGERES_SEQ_COLOR, 0.18))
        mageres_seq_band!(traces, xe[k], q(0.25), q(0.75), "spots 25–75 %", rgba(MAGERES_SEQ_COLOR, 0.35))
        push!(traces, scatter(; x = xe[k], y = q(0.5), mode = "lines", name = "spots median", text = hov[k],
                                line = attr(color = MAGERES_SEQ_COLOR, width = 1.5),
                                hovertemplate = "median %{y:.2f} kyr<extra>%{text}</extra>"))
        shown = false
        for (i, e) in enumerate(E)
            e.selected && !isempty(e.sample.grains) || continue
            g = e.sample.grains
            push!(traces, scatter(; x = xe[i] .+ 0.004 * span .* (rand(MersenneTwister(i), length(g)) .- 0.5), y = g,
                                    mode = "markers", name = "spots of the selected output", legendgroup = "sel",
                                    showlegend = !shown,
                                    marker = attr(color = MAGERES_SEQ_SEL_COLOR, size = 3, opacity = 0.6),
                                    hovertemplate = "%{y:.2f} kyr<extra>$(hov[i])</extra>"))
            shown = true
        end
    else
        for (i, e) in enumerate(E)
            g = e.sample.grains
            isempty(g) && continue
            col = e.selected ? MAGERES_SEQ_SEL_COLOR : MAGERES_SEQ_COLOR
            push!(traces, PlotlyJS.violin(; y = g, x0 = xe[i], width = width, name = hov[i], showlegend = false,
                                            points = "all", jitter = 0.5, pointpos = 0, spanmode = "hard",
                                            line = attr(color = col, width = 1), fillcolor = rgba(col, 0.3),
                                            marker = attr(color = col, size = 3, opacity = 0.7),
                                            hoveron = "violins+points",
                                            hovertemplate = "%{y:.2f} kyr<extra>$(hov[i])</extra>"))
        end
    end
    m = [i for i in eachindex(S) if isfinite(S[i].q50)]
    if !isempty(m)
        push!(traces, scatter(; x = xe[m], y = [S[i].q50 for i in m], mode = "lines", name = "zircon mass median",
                                text = hov[m], line = attr(color = "black", width = 1.5),
                                hovertemplate = "mass median %{y:.2f} kyr<extra>%{text}</extra>"))
        push!(traces, scatter(; x = xe[m], y = [S[i].q95 for i in m], mode = "lines", name = "zircon mass 5–95 %",
                                legendgroup = "mass", line = attr(color = "black", width = 1, dash = "dot"),
                                text = hov[m], hovertemplate = "mass 95 %: %{y:.2f} kyr<extra>%{text}</extra>"))
        push!(traces, scatter(; x = xe[m], y = [S[i].q05 for i in m], mode = "lines", name = "zircon mass 5–95 %",
                                legendgroup = "mass", showlegend = false,
                                line = attr(color = "black", width = 1, dash = "dot"),
                                text = hov[m], hovertemplate = "mass 5 %: %{y:.2f} kyr<extra>%{text}</extra>"))
    end
    push!(traces, scatter(; x = [x_lo, x_hi], y = [x_lo, x_hi], mode = "lines", name = "zircon age = eruption age",
                            line = attr(color = "grey", dash = "dash", width = 1), hoverinfo = "skip"))
    tb = mageres_injection_T_trace(T_injection, seq.t_ref, x_lo, x_hi; yaxis = "y2")
    tb === nothing || push!(traces, tb)
    note = mageres_method_label(seq.opts) *
           @sprintf("<br>%d samples; red: selected output", length(E)) *
           (bands ? "<br>many samples: spot quantile bands (try one sample per output)" : "")
    layout = Layout(; title        = "Zircon age vs eruption age (synthetic single spots, autocrysts and antecrysts)",
                      height       = 560,
                      plot_bgcolor = "white",
                      showlegend   = true,
                      legend       = attr(orientation = "h", y = -0.18),
                      xaxis        = mageres_seq_axis(@sprintf("Eruption age before t = %.1f kyr (kyr)", seq.t_ref / 1e3);
                                                      range = [x_hi, x_lo]),
                      yaxis        = mageres_seq_axis(@sprintf("Zircon age before t = %.1f kyr (kyr)", seq.t_ref / 1e3);
                                                      mirror = tb === nothing),
                      yaxis2       = mageres_injection_T_axis(; overlaying = "y", anchor = "x"),
                      margin       = attr(l = 80, r = tb === nothing ? 30 : 80, t = 50, b = 60),
                      violinmode   = "overlay",
                      annotations  = [attr(x = 0.01, xref = "paper", xanchor = "left", y = 0.99, yref = "paper",
                                           yanchor = "top", showarrow = false, align = "left", text = note,
                                           bgcolor = "rgba(255,255,255,0.85)", bordercolor = "#ccc", borderwidth = 1,
                                           font = attr(size = 11))])
    return PlotlyJS.plot(traces, layout)
end

"""
    mageres_panel_legend(y::Float64)

    Horizontal legend whose top is at `y` (paper coordinates), left-aligned under a panel.
"""
mageres_panel_legend(y::Float64) =
    attr(orientation = "h", x = 0.0, xanchor = "left", y = y, yanchor = "top", bgcolor = "rgba(255,255,255,0)",
         font = attr(size = 11), tracegroupgap = 2, traceorder = "normal")

"""
    mageres_zircon_sequence_stats_figure(seq; T_injection = nothing)

    Statistics of the zircon of each sample against eruption age, in four panels, each with its legend below it:
    (a) age before eruption of the erupted zircon mass (5–95 % band, mean = mean residence time, youngest zircon)
    with the weighted mean of the synthetic spots ± 2σ; (b) MSWD of the spots; (c) skewness and excess kurtosis of
    the zircon mass; (d) the erupted zircon mass split into autocrysts (crystallised since the previous eruption),
    antecrysts (crystallised before it) and inherited zircon, stacked to 100 %.
    With `T_injection` (times [yr], temperatures [°C]), the temperature of every injection is drawn on the right axis of (a).
"""
function mageres_zircon_sequence_stats_figure(seq; T_injection = nothing)
    p = mageres_sequence_empty(seq)
    p === nothing || return p
    S   = seq.stats
    o   = sortperm([s.erupt_age for s in S]; rev = true)
    S   = S[o]
    E   = seq.entries[o]
    x   = [s.erupt_age for s in S]
    hov = [mageres_seq_label(e) for e in E]
    dom = [[0.0, 0.17], [0.29, 0.44], [0.54, 0.68], [0.79, 1.0]]
    span = max(maximum(x) - minimum(x), 1.0)
    x_lo = minimum(x) - 0.03 * span
    x_hi = maximum(x) + 0.03 * span
    tr  = GenericTrace[]
    ln(y, name, col, ax, lg; dash = "solid", kw...) =
        scatter(; x = x, y = y, mode = "lines+markers", name = name, xaxis = "x", yaxis = ax, legend = lg, text = hov,
                  line = attr(color = col, dash = dash, width = 1.5), marker = attr(color = col, size = 5),
                  hovertemplate = "%{y:.3g}<extra>%{text}</extra>", kw...)
    push!(tr, scatter(; x = x, y = [s.q95 - s.erupt_age for s in S], mode = "lines", yaxis = "y4", xaxis = "x",
                        legend = "legend", line = attr(width = 0), showlegend = false, hoverinfo = "skip"))
    push!(tr, scatter(; x = x, y = [s.q05 - s.erupt_age for s in S], mode = "lines", yaxis = "y4", xaxis = "x",
                        legend = "legend", fill = "tonexty", fillcolor = "rgba(0,114,189,0.18)", line = attr(width = 0),
                        name = "zircon mass 5–95 %", hoverinfo = "skip"))
    push!(tr, ln([s.residence for s in S], "mean (residence time)", MAGERES_SEQ_COLOR, "y4", "legend"))
    push!(tr, ln([s.youngest - s.erupt_age for s in S], "youngest zircon", "rgb(60,160,60)", "y4", "legend";
                 dash = "dot"))
    push!(tr, scatter(; x = x, y = [s.wm - s.erupt_age for s in S], mode = "markers", name = "spots: weighted mean ± 2σ",
                        xaxis = "x", yaxis = "y4", legend = "legend", text = hov,
                        marker = attr(color = "black", size = 6, symbol = "diamond"),
                        error_y = attr(type = "data", array = [2s.wm_sigma for s in S], visible = true,
                                       color = "rgb(90,90,90)", thickness = 1, width = 0),
                        hovertemplate = "%{y:.2f} kyr<extra>%{text}</extra>"))
    tb = mageres_injection_T_trace(T_injection, seq.t_ref, x_lo, x_hi; yaxis = "y5", legend = "legend")
    tb === nothing || push!(tr, tb)
    push!(tr, ln([s.mswd for s in S], "MSWD of spots", "black", "y3", "legend2"))
    push!(tr, scatter(; x = [minimum(x), maximum(x)], y = [1.0, 1.0], mode = "lines", xaxis = "x", yaxis = "y3",
                        legend = "legend2", name = "MSWD = 1", line = attr(color = "grey", dash = "dash", width = 1),
                        hoverinfo = "skip"))
    push!(tr, ln([s.skew for s in S], "skewness", "rgb(230,140,30)", "y2", "legend3"))
    push!(tr, ln([s.kurt for s in S], "excess kurtosis", "rgb(126,47,142)", "y2", "legend3"; dash = "dot"))
    inh  = [isfinite(s.inherited) ? s.inherited : 0.0 for s in S]
    ante = [isfinite(s.antecryst) ? s.antecryst : 0.0 for s in S]
    for (name, y, col) in (("autocrysts (since previous eruption)", 100 .* (1 .- inh) .* (1 .- ante), "rgb(0,114,189)"),
                           ("antecrysts (before previous eruption)", 100 .* (1 .- inh) .* ante, "rgb(200,40,40)"),
                           ("inherited", 100 .* inh, "rgb(150,150,150)"))
        push!(tr, scatter(; x = x, y = y, mode = "lines", stackgroup = "fractions", name = name, xaxis = "x",
                            yaxis = "y", legend = "legend4", text = hov, line = attr(color = col, width = 0.5),
                            fillcolor = replace(col, "rgb" => "rgba", ")" => ",0.6)"),
                            hovertemplate = "%{y:.1f} %<extra>%{text}</extra>"))
    end
    layout = Layout(; title        = "Zircon statistics through the eruption sequence",
                      height       = 1150,
                      plot_bgcolor = "white",
                      showlegend   = true,
                      margin       = attr(l = 80, r = tb === nothing ? 30 : 80, t = 50, b = 110),
                      xaxis        = mageres_seq_axis(@sprintf("Eruption age before t = %.1f kyr (kyr)", seq.t_ref / 1e3);
                                                      anchor = "y", range = [x_hi, x_lo]),
                      yaxis        = mageres_seq_axis("Erupted zircon mass (%)"; domain = dom[1], anchor = "x",
                                                      range = [0, 100]),
                      yaxis2       = mageres_seq_axis("Skewness, kurtosis"; domain = dom[2], anchor = "x"),
                      yaxis3       = mageres_seq_axis("MSWD"; domain = dom[3], anchor = "x", type = "log"),
                      yaxis4       = mageres_seq_axis("Age before eruption (kyr)"; domain = dom[4], anchor = "x",
                                                      rangemode = "tozero", mirror = tb === nothing),
                      yaxis5       = mageres_injection_T_axis(; overlaying = "y4", anchor = "x"),
                      legend       = mageres_panel_legend(dom[4][1] - 0.012),
                      legend2      = mageres_panel_legend(dom[3][1] - 0.012),
                      legend3      = mageres_panel_legend(dom[2][1] - 0.012),
                      legend4      = mageres_panel_legend(-0.075))
    return PlotlyJS.plot(tr, layout)
end

"""
    mageres_estimator_summary(err::Vector{Float64}, sig::Vector{Float64})

    Mean bias, RMS error [kyr] and fraction within 2σ of the finite estimates `err` with 1σ `sig`.
"""
function mageres_estimator_summary(err::Vector{Float64}, sig::Vector{Float64})
    k = [i for i in eachindex(err) if isfinite(err[i])]
    isempty(k) && return NaN, NaN, NaN
    e = err[k]
    return sum(e) / length(e), sqrt(sum(e .^ 2) / length(e)), count(abs.(e) .<= 2 .* sig[k]) / length(k)
end

"""
    mageres_eruption_age_test_figure(seq)

    Eruption-age estimators tested against the known eruption ages: (a) estimate − true eruption age ± 2σ per sample
    for every estimator of `MAGERES_AGE_ESTIMATORS`, with bias, RMS error and 2σ coverage over the sequence;
    (b) median absolute error and (c) 2σ coverage against the number of dated grains (`mageres_seq_grain_test`).
"""
function mageres_eruption_age_test_figure(seq)
    p = mageres_sequence_empty(seq)
    p === nothing || return p
    S    = seq.stats
    E    = seq.entries
    x    = [s.erupt_age for s in S]
    span = max(maximum(x) - minimum(x), 1.0)
    nm   = length(MAGERES_AGE_ESTIMATORS)
    hov  = [mageres_seq_label(e) for e in E]
    tr   = GenericTrace[]
    for (j, (_, name)) in enumerate(MAGERES_AGE_ESTIMATORS)
        err = [seq.est[i, j][1] for i in eachindex(E)]
        sig = [seq.est[i, j][2] for i in eachindex(E)]
        col = MAGERES_ESTIMATOR_COLORS[j]
        b, r, c = mageres_estimator_summary(err, sig)
        label   = @sprintf("%s: bias %+.2f, RMS %.2f kyr, %.0f %% within 2σ", name, b, r, 100c)
        push!(tr, scatter(; x = x .+ 0.006 * span * (j - (nm + 1) / 2), y = err, mode = "markers", name = label,
                            legendgroup = name, legend = "legend", xaxis = "x", yaxis = "y", text = hov,
                            marker = attr(color = col, size = 6, line = attr(color = "black", width = 0.5)),
                            error_y = attr(type = "data", array = 2 .* sig, visible = true, color = col,
                                           thickness = 1, width = 0),
                            hovertemplate = "%{y:.2f} kyr<extra>%{text}</extra>"))
    end
    push!(tr, scatter(; x = [minimum(x) - 0.05span, maximum(x) + 0.05span], y = [0.0, 0.0], mode = "lines",
                        xaxis = "x", yaxis = "y", line = attr(color = "black", dash = "dash", width = 1),
                        showlegend = false, hoverinfo = "skip"))
    Ns, med, cv = mageres_seq_grain_test(E, seq.t_ref)
    if !isempty(Ns)
        for (j, (_, name)) in enumerate(MAGERES_AGE_ESTIMATORS)
            col = MAGERES_ESTIMATOR_COLORS[j]
            push!(tr, scatter(; x = Ns, y = med[:, j], mode = "lines+markers", name = name, legendgroup = name,
                                legend = "legend2", xaxis = "x2", yaxis = "y2", line = attr(color = col),
                                marker = attr(color = col, size = 6),
                                hovertemplate = "N = %{x}: %{y:.2f} kyr<extra>$(name)</extra>"))
            push!(tr, scatter(; x = Ns, y = 100 .* cv[:, j], mode = "lines+markers", name = name, legendgroup = name,
                                legend = "legend2", showlegend = false, xaxis = "x2", yaxis = "y3",
                                line = attr(color = col),
                                marker = attr(color = col, size = 6),
                                hovertemplate = "N = %{x}: %{y:.0f} %<extra>$(name)</extra>"))
        end
        push!(tr, scatter(; x = [Ns[1], Ns[end]], y = [95.0, 95.0], mode = "lines", xaxis = "x2", yaxis = "y3",
                            legend = "legend2", name = "95 % (ideal 2σ coverage)",
                            line = attr(color = "grey", dash = "dash", width = 1), hoverinfo = "skip"))
    end
    layout = Layout(; title        = "Eruption-age estimators tested against the model eruption ages (" *
                                     mageres_method_label(seq.opts) * ")",
                      height       = 1100,
                      plot_bgcolor = "white",
                      showlegend   = true,
                      margin       = attr(l = 80, r = 30, t = 50, b = 110),
                      xaxis        = mageres_seq_eruption_axis(seq; anchor = "y"),
                      yaxis        = mageres_seq_axis("Estimate − eruption age (kyr)"; domain = [0.67, 1.0],
                                                      anchor = "x"),
                      xaxis2       = mageres_seq_axis("Number of dated grains N"; type = "log", anchor = "y3",
                                                      tickvals = Ns, ticktext = string.(Ns)),
                      yaxis2       = mageres_seq_axis("Median |error| (kyr)"; domain = [0.27, 0.45], anchor = "x2",
                                                      rangemode = "tozero"),
                      yaxis3       = mageres_seq_axis("Within 2σ (%)"; domain = [0.0, 0.19], anchor = "x2",
                                                      range = [0, 105]),
                      legend       = mageres_panel_legend(0.585),
                      legend2      = mageres_panel_legend(-0.075),
                      annotations  = [attr(x = 0.5, xref = "paper", xanchor = "center", y = 0.47, yref = "paper",
                                           yanchor = "bottom", showarrow = false, font = attr(size = 12),
                                           text = @sprintf("Random subsets of N spots: %d draws for each of up to %d samples",
                                                           MAGERES_SEQ_N_DRAWS, MAGERES_SEQ_N_ERUPT))])
    return PlotlyJS.plot(tr, layout)
end

"""
    mageres_sequence_csv(seq)

    CSV table of the sequence, one row per sample: eruption times and age, erupted area and temperature, spot counts,
    zircon-mass and spot statistics, and every eruption-age estimate with its 1σ and error [kyr].
"""
function mageres_sequence_csv(seq::MAGEResSequence)
    f(v) = isfinite(v) ? @sprintf("%.6g", v) : "NaN"
    est  = [m for (m, _) in MAGERES_AGE_ESTIMATORS]
    head = vcat(["first_eruption_time_yr", "last_eruption_time_yr", "n_eruptions",
                 "eruption_age_kyr_before_t_ref", "erupted_area_m2", "erupted_T_C", "selected_output",
                 "n_spots", "n_inherited_spots", "mass_mean_age_kyr", "mass_q05_age_kyr", "mass_q50_age_kyr",
                 "mass_q95_age_kyr", "mass_sd_kyr", "mass_skewness", "mass_excess_kurtosis",
                 "youngest_zircon_age_kyr", "mean_residence_kyr", "antecryst_fraction", "inherited_fraction",
                 "spot_weighted_mean_kyr", "spot_weighted_mean_1sigma_kyr", "spot_mswd", "spot_q05_kyr",
                 "spot_q95_kyr"],
                reduce(vcat, (["$(m)_age_kyr", "$(m)_1sigma_kyr", "$(m)_error_kyr"] for m in est)))
    lines = [@sprintf("# MAGERes eruption sequence; ages in kyr before t_ref = %.1f yr; %s", seq.t_ref,
                      mageres_method_label(seq.opts)), join(head, ",")]
    for (i, (e, s)) in enumerate(zip(seq.entries, seq.stats))
        z = e.sample
        row = vcat([f(z.t_first), f(z.t), string(z.n_eruptions), f(s.erupt_age), f(e.area), f(e.T),
                    string(e.selected), string(s.n_spots), string(s.n_inh)],
                   f.([s.mean, s.q05, s.q50, s.q95, s.sd, s.skew, s.kurt, s.youngest, s.residence, s.antecryst,
                       s.inherited, s.wm, s.wm_sigma, s.mswd, s.spot_q05, s.spot_q95]),
                   reduce(vcat, ([f(seq.est[i, j][1] + s.erupt_age), f(seq.est[i, j][2]), f(seq.est[i, j][1])]
                                 for j in eachindex(est))))
        push!(lines, join(row, ","))
    end
    return join(lines, "\n") * "\n"
end
