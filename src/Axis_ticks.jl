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

const TICK_QUANTITIES      = ("T", "P", "X", "mu", "other")
const TICK_QUANTITY_LABELS = Dict("T" => "Temperature", "P" => "Pressure", "X" => "Composition",
                                  "mu" => "Chemical potential μ", "other" => "Other axes")
const TICK_MAX_MAJOR       = 200
const TICK_MAX_MINOR       = 2000

default_axis_tick_axis()  = Dict{String,Any}("mode" => "divisions", "n_div" => 5, "step" => nothing, "start" => nothing,
                                             "n_minor" => 0, "label_every" => 1, "decimals" => nothing)
default_axis_tick_style() = Dict{String,Any}("direction" => "outside", "major_len" => 4.0, "minor_len" => 2.0, "mirror" => true)
default_axis_tick_spec()  = Dict{String,Any}("version" => 1, "style" => default_axis_tick_style(),
                                             "axes" => Dict{String,Any}(q => default_axis_tick_axis() for q in TICK_QUANTITIES))

const AXIS_TICKS = Ref{Dict{String,Any}}(default_axis_tick_spec())

function tick_get(d, k::AbstractString)
    d isa AbstractDict || return nothing
    haskey(d, k)         && return d[k]
    haskey(d, Symbol(k)) && return d[Symbol(k)]
    return nothing
end

function tick_float(v)
    v isa Bool && return nothing
    v isa Real && return isfinite(v) ? Float64(v) : nothing
    if v isa AbstractString
        x = tryparse(Float64, strip(v))
        return x !== nothing && isfinite(x) ? x : nothing
    end
    return nothing
end

function tick_int(v, default::Int, lo::Int, hi::Int)
    x = tick_float(v)
    return x === nothing ? default : clamp(round(Int, x), lo, hi)
end

function tick_bool(v, default::Bool)
    v isa Bool && return v
    v isa AbstractString && lowercase(v) in ("true", "false") && return lowercase(v) == "true"
    v isa AbstractVector && return !isempty(v)
    return default
end

function parse_axis_tick_axis(d)
    out  = default_axis_tick_axis()
    mode = tick_get(d, "mode")
    out["mode"]        = mode isa AbstractString && mode in ("divisions", "step") ? String(mode) : "divisions"
    out["n_div"]       = tick_int(tick_get(d, "n_div"), 5, 1, 50)
    step               = tick_float(tick_get(d, "step"))
    out["step"]        = step !== nothing && step > 0 ? step : nothing
    out["start"]       = tick_float(tick_get(d, "start"))
    out["n_minor"]     = tick_int(tick_get(d, "n_minor"), 0, 0, 10)
    out["label_every"] = tick_int(tick_get(d, "label_every"), 1, 1, 50)
    dec                = tick_float(tick_get(d, "decimals"))
    out["decimals"]    = dec === nothing ? nothing : clamp(round(Int, dec), 0, 6)
    return out
end

function parse_axis_tick_style(d)
    out = default_axis_tick_style()
    dir = tick_get(d, "direction")
    out["direction"] = dir isa AbstractString && dir in ("outside", "inside") ? String(dir) : "outside"
    for (k, def) in (("major_len", 4.0), ("minor_len", 2.0))
        x = tick_float(tick_get(d, k))
        out[k] = x === nothing ? def : clamp(x, 0.0, 20.0)
    end
    out["mirror"] = tick_bool(tick_get(d, "mirror"), true)
    return out
end

function parse_axis_tick_spec(d)
    out  = default_axis_tick_spec()
    out["style"] = parse_axis_tick_style(tick_get(d, "style"))
    axes = tick_get(d, "axes")
    for q in TICK_QUANTITIES
        out["axes"][q] = parse_axis_tick_axis(tick_get(axes, q))
    end
    return out
end

function tick_quantity(title)
    t = lowercase(strip(string(something(title, ""))))
    occursin("path", t)                          && return "other"
    startswith(t, "μ")                           && return "mu"
    startswith(t, "composition")                 && return "X"
    occursin(r"^(pressure|p\s*[\(\[])", t)      && return "P"
    (occursin(r"(^|\s)t\s*[\(\[]", t) || occursin("temperature", t)) && return "T"
    return "other"
end

tick_scale(title) = tick_quantity(title) == "P" && occursin("gpa", lowercase(string(something(title, "")))) ? 0.1 : 1.0

function tick_snap(v::Float64)
    r = round(v, digits = 10)
    return r == 0.0 ? 0.0 : r
end

function tick_label(v::Real, decimals)
    decimals === nothing && return svg_tick_label(v)
    s = Printf.format(Printf.Format("%.$(decimals)f"), Float64(v))
    return occursin(r"^-0(\.0*)?$", s) ? s[2:end] : s
end

function axis_tick_values(ax::AbstractDict, lo::Real, hi::Real; scale::Real = 1.0)
    a, b = Float64(min(lo, hi)), Float64(max(lo, hi))
    span = b - a
    if !(isfinite(span) && span > 0)
        return (major = [a], minor = Float64[], labels = [tick_label(a, ax["decimals"])])
    end
    tol   = 1e-9 * span
    major = Float64[]
    step  = 0.0
    if ax["mode"] == "step" && ax["step"] !== nothing
        st = Float64(ax["step"]) * scale
        s0 = ax["start"] === nothing ? 0.0 : Float64(ax["start"]) * scale
        k0 = ceil((a - tol - s0) / st)
        k1 = floor((b + tol - s0) / st)
        if 2 <= k1 - k0 + 1 <= TICK_MAX_MAJOR
            major = [tick_snap(s0 + k * st) for k in k0:k1]
            step  = st
        end
    end
    if isempty(major)
        n     = ax["n_div"]
        step  = span / n
        major = [k == n ? b : tick_snap(a + k * step) for k in 0:n]
    end

    minor = Float64[]
    m     = ax["n_minor"]
    if m > 0
        dm   = step / (m + 1)
        base = major[1]
        k0   = ceil(Int, (a - tol - base) / dm)
        k1   = floor(Int, (b + tol - base) / dm)
        if k1 - k0 + 1 <= TICK_MAX_MINOR
            minor = [tick_snap(base + k * dm) for k in k0:k1 if mod(k, m + 1) != 0]
        end
    end

    e      = ax["label_every"]
    labels = [mod(i - 1, e) == 0 ? tick_label(v, ax["decimals"]) : "" for (i, v) in enumerate(major)]
    return (major = major, minor = minor, labels = labels)
end

axis_ticks_for(spec, title, lo, hi) =
    axis_tick_values(spec["axes"][tick_quantity(title)], lo, hi; scale = tick_scale(title))

function plotly_axis_tick_attr!(ax::AbstractDict, t, style)
    dir = style["direction"]
    ax[:tickmode]  = "array"
    ax[:tickvals]  = t.major
    ax[:ticktext]  = t.labels
    ax[:ticks]     = dir
    ax[:ticklen]   = style["major_len"]
    ax[:tickwidth] = 1
    ax[:tickcolor] = "black"
    ax[:showline]  = true
    ax[:linecolor] = "black"
    ax[:linewidth] = 1
    ax[:mirror]    = style["mirror"] ? "ticks" : true
    ax[:showgrid]  = false
    ax[:zeroline]  = false
    ax[:minor]     = Dict{Any,Any}(:tickmode => "array", :tickvals => t.minor, :ticks => isempty(t.minor) ? "" : dir,
                                   :ticklen => style["minor_len"], :tickwidth => 0.75, :tickcolor => "black", :showgrid => false)
    for k in (:tick0, :dtick, :nticks)
        delete!(ax, k)
    end
    return ax
end

function layout_axis!(layout, key::Symbol)
    old = get(layout.fields, key, nothing)
    ax  = Dict{Any,Any}()
    if old isa AbstractDict
        for (k, v) in (old isa PlotlyBase.PlotlyAttribute ? old.fields : old)
            ax[Symbol(k)] = v
        end
    end
    layout.fields[key] = ax
    return ax
end

function layout_axis_title(ax::AbstractDict)
    t = get(ax, :title, "")
    t isa AbstractDict && return string(something(get(t, :text, nothing), ""))
    return string(something(t, ""))
end

function apply_axis_ticks!(layout, spec; xr = nothing, yr = nothing, xtitle = nothing, ytitle = nothing)
    style = spec["style"]
    for (key, r, t) in ((:xaxis, xr, xtitle), (:yaxis, yr, ytitle))
        ax = layout_axis!(layout, key)
        rr = r === nothing ? get(ax, :range, nothing) : r
        (rr === nothing || length(rr) < 2) && continue
        tt = t === nothing ? layout_axis_title(ax) : t
        plotly_axis_tick_attr!(ax, axis_ticks_for(spec, tt, Float64(rr[1]), Float64(rr[end])), style)
    end
    return layout
end

function ticked_layout(layout_in, spec; logo::Bool = true)
    l = Layout(copy(layout_in.fields))
    logo && (l.fields[:images] = diagram_logo_images())
    return apply_axis_ticks!(l, spec)
end

function svg_axis_tick_kwargs(spec, xt, yt)
    style = spec["style"]
    return (xticks = xt.major, yticks = yt.major, xticklabels = xt.labels, yticklabels = yt.labels,
            xminor = xt.minor, yminor = yt.minor, tick_len = Float64(style["major_len"]),
            minor_len = Float64(style["minor_len"]), tick_dir = String(style["direction"]), mirror = style["mirror"] == true)
end

svg_axis_tick_kwargs(spec, xtitle, ytitle, xr, yr) =
    svg_axis_tick_kwargs(spec, axis_ticks_for(spec, xtitle, xr[1], xr[end]), axis_ticks_for(spec, ytitle, yr[1], yr[end]))

current_axis_tick_spec(store_data) = store_data === nothing ? AXIS_TICKS[] : parse_axis_tick_spec(store_data)

axis_ticks_path() = joinpath(mageminapp_user_dir(), "axis_ticks.json")

function load_axis_tick_prefs(; path::AbstractString = axis_ticks_path())
    isfile(path) || return default_axis_tick_spec()
    try
        return parse_axis_tick_spec(JSON.parsefile(path))
    catch e
        @warn "Failed to read axis tick preferences, using defaults" path=path exception=e
        return default_axis_tick_spec()
    end
end

function save_axis_tick_prefs(spec; path::AbstractString = axis_ticks_path())
    try
        mkpath(dirname(path))
        tmp = path * ".tmp"
        open(tmp, "w") do io
            JSON.print(io, spec, 2)
        end
        mv(tmp, path; force = true)
        return true
    catch e
        @warn "Failed to save axis tick preferences" path=path exception=e
        return false
    end
end

tick_unit_hint(q::AbstractString, gpa::Bool) =
    q == "T"  ? "Step and start in °C." :
    q == "P"  ? "Step and start in $(gpa ? "GPa" : "kbar")." :
    q == "X"  ? "Step and start as composition fraction (0-1)." :
    q == "mu" ? "Step and start in kJ/mol." : "Step and start in axis units."

tick_p_factor(q::AbstractString, gpa::Bool) = q == "P" && gpa ? 0.1 : 1.0

tick_disp(v, f) = v === nothing ? nothing : round(v * f, sigdigits = 12)

function axis_tick_fields(spec, q::AbstractString, gpa::Bool)
    ax, st = spec["axes"][q], spec["style"]
    f = tick_p_factor(q, gpa)
    return (ax["mode"], ax["n_div"], tick_disp(ax["step"], f), tick_disp(ax["start"], f), ax["n_minor"], ax["label_every"],
            ax["decimals"], st["direction"], st["major_len"], st["minor_len"], st["mirror"])
end

function axis_tick_spec_from_fields(spec, q::AbstractString, gpa::Bool, mode, n_div, step, start, n_minor, every, decimals,
                                    dir, majlen, minlen, mirror)
    out = deepcopy(spec)
    f   = tick_p_factor(q, gpa)
    s, s0 = tick_float(step), tick_float(start)
    out["axes"][q] = Dict{String,Any}("mode" => mode, "n_div" => n_div, "step" => tick_disp(s, 1 / f), "start" => tick_disp(s0, 1 / f),
                                      "n_minor" => n_minor, "label_every" => every, "decimals" => decimals)
    out["style"] = Dict{String,Any}("direction" => dir, "major_len" => majlen, "minor_len" => minlen, "mirror" => mirror)
    return parse_axis_tick_spec(out)
end

function axis_ticks_callbacks!(app)
    callback!(
        app,
        Output("axis-ticks-store",    "data"),
        Output("axis-ticks-status",   "children"),
        Output("axis-ticks-unit",     "children"),
        Output("axis-ticks-mode",     "value"),
        Output("axis-ticks-ndiv",     "value"),
        Output("axis-ticks-step",     "value"),
        Output("axis-ticks-start",    "value"),
        Output("axis-ticks-nminor",   "value"),
        Output("axis-ticks-every",    "value"),
        Output("axis-ticks-decimals", "value"),
        Output("axis-ticks-dir",      "value"),
        Output("axis-ticks-majlen",   "value"),
        Output("axis-ticks-minlen",   "value"),
        Output("axis-ticks-mirror",   "value"),

        Input("axis-ticks-apply",       "n_clicks"),
        Input("axis-ticks-reset",       "n_clicks"),
        Input("axis-ticks-quantity",    "value"),
        Input("pressure-unit-dropdown", "value"),

        State("axis-ticks-store",    "data"),
        State("axis-ticks-mode",     "value"),
        State("axis-ticks-ndiv",     "value"),
        State("axis-ticks-step",     "value"),
        State("axis-ticks-start",    "value"),
        State("axis-ticks-nminor",   "value"),
        State("axis-ticks-every",    "value"),
        State("axis-ticks-decimals", "value"),
        State("axis-ticks-dir",      "value"),
        State("axis-ticks-majlen",   "value"),
        State("axis-ticks-minlen",   "value"),
        State("axis-ticks-mirror",   "value"),
    ) do _apply, _reset, q, pressure_unit, store, fields...

        bid = pushed_button( callback_context() )
        q   = q in TICK_QUANTITIES ? String(q) : "T"
        gpa = pressure_unit == "gpa"

        status = ""
        if bid == "axis-ticks-apply"
            AXIS_TICKS[] = axis_tick_spec_from_fields(AXIS_TICKS[], q, gpa, fields...)
            status = save_axis_tick_prefs(AXIS_TICKS[]) ? "$(TICK_QUANTITY_LABELS[q]) ticks applied and saved to $(axis_ticks_path())." :
                                                          "$(TICK_QUANTITY_LABELS[q]) ticks applied (preferences could not be saved)."
        elseif bid == "axis-ticks-reset"
            AXIS_TICKS[] = default_axis_tick_spec()
            status = save_axis_tick_prefs(AXIS_TICKS[]) ? "All axis ticks reset to defaults." : "Reset to defaults (preferences could not be saved)."
        end

        spec     = AXIS_TICKS[]
        store_up = (store !== nothing && parse_axis_tick_spec(store) == spec) ? no_update() : spec
        return store_up, status, tick_unit_hint(q, gpa), axis_tick_fields(spec, q, gpa)...
    end
    return app
end
