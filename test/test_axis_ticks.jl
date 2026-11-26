using Test
using MAGEMinApp
using PlotlyJS

const M = MAGEMinApp

spec_with(q; style = Dict{String,Any}(), kw...) = begin
    s = M.default_axis_tick_spec()
    merge!(s["axes"][q], Dict{String,Any}(String(k) => v for (k, v) in kw))
    merge!(s["style"], style)
    M.parse_axis_tick_spec(s)
end

render_layout(; kw...) = begin
    io = IOBuffer()
    c  = M.SVGCanvas(758, 700, 80, 124, 60, 50, 400.0, 1200.0, 1.0, 25.0)
    M.svg_layout_layers(io, c, Set{String}(); kw...)
    String(take!(io))
end

@testset "tick quantity" begin
    @test M.tick_quantity("Temperature [Celsius]") == "T"
    @test M.tick_quantity("Temperature [°C]") == "T"
    @test M.tick_quantity("Event 2: temperature [Celsius]") == "T"
    @test M.tick_quantity("Event 2 T [°C]") == "T"
    @test M.tick_quantity("T(C)") == "T"
    @test M.tick_quantity("Pressure [kbar]") == "P"
    @test M.tick_quantity("Pressure [GPa]") == "P"
    @test M.tick_quantity("P [GPa]") == "P"
    @test M.tick_quantity("P(GPa)") == "P"
    @test M.tick_quantity("Composition [X0 -> X1]") == "X"
    @test M.tick_quantity("μ(SiO2) [kJ/mol]") == "mu"
    @test M.tick_quantity("Pressure-Temperature path") == "other"
    @test M.tick_quantity("") == "other"
    @test M.tick_quantity(nothing) == "other"
    @test M.tick_scale("Pressure [GPa]") == 0.1
    @test M.tick_scale("Pressure [kbar]") == 1.0
    @test M.tick_scale("Temperature [°C]") == 1.0
end

@testset "tick values" begin
    d = M.default_axis_tick_axis()
    t = M.axis_tick_values(d, 400.0, 1200.0)
    @test t.major == [400.0 + k * (1200.0 - 400.0) / 5 for k in 0:5]
    @test isempty(t.minor)
    @test t.labels == ["400", "560", "720", "880", "1040", "1200"]

    t = M.axis_tick_values(spec_with("P"; mode = "step", step = 5.0)["axes"]["P"], 1.0, 25.0)
    @test t.major == [5.0, 10.0, 15.0, 20.0, 25.0]

    t = M.axis_tick_values(spec_with("P"; mode = "step", step = 5.0, start = 1.0)["axes"]["P"], 1.0, 25.0)
    @test t.major == [1.0, 6.0, 11.0, 16.0, 21.0]

    t = M.axis_tick_values(spec_with("P"; mode = "step", step = 5.0)["axes"]["P"], 0.1, 2.5; scale = 0.1)
    @test t.major == [0.5, 1.0, 1.5, 2.0, 2.5]
    @test t.labels == ["0.5", "1", "1.5", "2", "2.5"]

    t = M.axis_tick_values(spec_with("X"; mode = "step", step = 0.1)["axes"]["X"], 0.3, 1.0)
    @test t.major == [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    @test t.labels == ["0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"]

    t = M.axis_tick_values(spec_with("T"; mode = "step", step = 1000.0)["axes"]["T"], 400.0, 1200.0)
    @test length(t.major) == 6

    t = M.axis_tick_values(spec_with("T"; mode = "step", step = 1e-4)["axes"]["T"], 400.0, 1200.0)
    @test length(t.major) == 6

    t = M.axis_tick_values(spec_with("T"; mode = "step", step = 200.0, n_minor = 1)["axes"]["T"], 450.0, 1250.0)
    @test t.major == [600.0, 800.0, 1000.0, 1200.0]
    @test t.minor == [500.0, 700.0, 900.0, 1100.0]

    t = M.axis_tick_values(spec_with("T"; n_div = 4, n_minor = 3)["axes"]["T"], 0.0, 100.0)
    @test length(t.minor) == 12
    @test isempty(intersect(t.minor, t.major))
    @test all(0.0 .<= t.minor .<= 100.0)

    t = M.axis_tick_values(spec_with("T"; label_every = 2, decimals = 1)["axes"]["T"], 0.0, 100.0)
    @test t.labels == ["0.0", "", "40.0", "", "80.0", ""]

    t = M.axis_tick_values(d, 5.0, 5.0)
    @test t.major == [5.0] && isempty(t.minor)

    t = M.axis_tick_values(d, 25.0, 1.0)
    @test first(t.major) == 1.0 && last(t.major) == 25.0
end

@testset "spec parsing" begin
    s = M.parse_axis_tick_spec(nothing)
    @test s == M.default_axis_tick_spec()
    s = M.parse_axis_tick_spec(Dict("style" => Dict("direction" => "sideways", "major_len" => "abc", "minor_len" => 99, "mirror" => "false"),
                                    "axes" => Dict("T" => Dict("mode" => "step", "n_div" => -3, "step" => -2, "n_minor" => 50,
                                                               "label_every" => 0, "decimals" => "2"))))
    @test s["style"]["direction"] == "outside"
    @test s["style"]["major_len"] == 4.0
    @test s["style"]["minor_len"] == 20.0
    @test s["style"]["mirror"] == false
    @test s["axes"]["T"]["n_div"] == 1
    @test s["axes"]["T"]["step"] === nothing
    @test s["axes"]["T"]["n_minor"] == 10
    @test s["axes"]["T"]["label_every"] == 1
    @test s["axes"]["T"]["decimals"] == 2
    @test s["axes"]["P"] == M.default_axis_tick_axis()
    sym = M.parse_axis_tick_spec(Dict(:axes => Dict(:X => Dict(:mode => "step", :step => 0.2))))
    @test sym["axes"]["X"]["mode"] == "step" && sym["axes"]["X"]["step"] == 0.2
end

@testset "fields round trip" begin
    s = spec_with("P"; mode = "step", step = 5.0, start = 2.0)
    f = M.axis_tick_fields(s, "P", true)
    @test f[3] == 0.5 && f[4] == 0.2
    s2 = M.axis_tick_spec_from_fields(s, "P", true, f...)
    @test s2 == s
    s3 = M.axis_tick_spec_from_fields(s, "P", true, "step", 5, 0.3, nothing, 1, 1, nothing, "inside", 5, 2, false)
    @test s3["axes"]["P"]["step"] == 3.0
    @test s3["style"]["direction"] == "inside" && s3["style"]["mirror"] == false
end

@testset "preferences file" begin
    dir  = mktempdir()
    path = joinpath(dir, "sub", "axis_ticks.json")
    @test M.load_axis_tick_prefs(path = path) == M.default_axis_tick_spec()
    s = spec_with("T"; mode = "step", step = 100.0, n_minor = 1)
    @test M.save_axis_tick_prefs(s; path = path)
    @test M.load_axis_tick_prefs(path = path) == s
    write(path, "{not json")
    @test M.load_axis_tick_prefs(path = path) == M.default_axis_tick_spec()
end

@testset "plotly layout" begin
    l = Layout(xaxis_title = "Temperature [Celsius]", yaxis_title = "Pressure [GPa]", xaxis_range = [400.0, 1200.0],
               yaxis_range = [0.1, 2.5], xaxis = attr(fixedrange = true, tick0 = 400.0, dtick = 160.0),
               yaxis = attr(fixedrange = true), images = [attr(source = "a"), attr(source = "b")])
    s  = spec_with("P"; mode = "step", step = 5.0, n_minor = 4)
    l2 = M.ticked_layout(l, s)
    @test length(l.fields[:images]) == 2
    @test length(l2.fields[:images]) == 1
    @test !haskey(l.fields[:xaxis], :tickvals)
    xa, ya = l2.fields[:xaxis], l2.fields[:yaxis]
    @test xa[:tickvals] == [400.0 + k * 160.0 for k in 0:5]
    @test !haskey(xa, :dtick) && !haskey(xa, :tick0)
    @test xa[:fixedrange] == true
    @test ya[:tickvals] == [0.5, 1.0, 1.5, 2.0, 2.5]
    @test length(ya[:minor][:tickvals]) == 4 * 4 + 4
    @test xa[:minor][:ticks] == ""
    @test ya[:mirror] == "ticks" && ya[:ticks] == "outside"
end

@testset "svg layout" begin
    golden_a = read(joinpath(@__DIR__, "golden_layout_a.svgfrag"), String)
    golden_b = read(joinpath(@__DIR__, "golden_layout_b.svgfrag"), String)
    @test render_layout(xticks = [400 + k * 160.0 for k in 0:5], yticks = [1 + k * 4.8 for k in 0:5],
                        xtitle = "Temperature [°C]", ytitle = "Pressure [kbar]", title = "Test") == golden_a
    @test render_layout(xticks = [500.0, 800.0], yticks = [5.0, 10.0], xtitle = "", ytitle = "Y", title = "",
                        xticklabels = ["a", "b"], xtick_angle = -45) == golden_b

    s  = spec_with("P"; style = Dict("direction" => "inside", "mirror" => false), mode = "step", step = 5.0, n_minor = 4)
    s["axes"]["T"] = M.parse_axis_tick_axis(Dict("mode" => "step", "step" => 200.0, "n_minor" => 1, "label_every" => 2))
    kw = M.svg_axis_tick_kwargs(s, "Temperature [°C]", "Pressure [kbar]", [400.0, 1200.0], [1.0, 25.0])
    svg = render_layout(; kw..., xtitle = "Temperature [°C]", ytitle = "Pressure [kbar]", title = "")

    ids = [m.captures[1] for m in eachmatch(r"id=\"([^\"]+)\"", svg)]
    @test allunique(ids)
    @test count("<g ", svg) == count("</g>", svg)
    @test !occursin("<image", svg)
    @test !occursin(r"<text id=\"Tick_[^>]*><tspan[^>]*>\u00a0</tspan>", svg)
    minor_block = match(r"<g id=\"Ticks_minor\"[^>]*>(.*?)</g>"s, svg)
    @test minor_block !== nothing
    n_minor = length(kw.xminor) + length(kw.yminor)
    @test count("<line", minor_block.captures[1]) == n_minor
    @test count("<line", svg) == n_minor + length(kw.xticks) + length(kw.yticks)
    labels = [m.captures[1] for m in eachmatch(r"<text id=\"Tick_[xy]_\d+\"[^>]*><tspan[^>]*>([^<]*)</tspan>", svg)]
    @test all(l -> !occursin(r"\d\.\d{6,}", l), labels)
    @test "400" in labels && !("600" in labels) && "800" in labels
    xs = [parse(Float64, m.captures[1]) for m in eachmatch(r"<line x1=\"([-0-9.]+)\"", svg)]
    @test all(80 - 4 .<= xs .<= 758 - 124 + 4)

    @test !occursin("Ticks_minor", render_layout(; M.svg_axis_tick_kwargs(M.default_axis_tick_spec(), "Temperature [°C]", "Pressure [kbar]",
                                                                           [400.0, 1200.0], [1.0, 25.0])...,
                                                   xtitle = "T", ytitle = "P", title = ""))
end
