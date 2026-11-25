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
using Test
using Statistics
using MAGEMinApp

const MV = MAGEMinApp

if !isdefined(MV, :mageres_comp_sigma)
    @info "MAGERes disabled in this build: composition-variability tests skipped"
else

cv(; kw...) = (o = MV.MAGEResOptions(; kw...); MV.mageres_validate_options(o); o)

@testset "variability law" begin
    o = cv()
    @test MV.mageres_comp_sigma(o, 0.01) ≈ 0.04
    @test MV.mageres_comp_sigma(o, 0.66) ≈ 0.04 * 66^-0.15
    @test MV.mageres_comp_sigma(o, 1e-6) ≈ 0.04 * 1e-4^-0.15
    @test MV.mageres_comp_sigma(o, 1e-12) == 0.5
    @test MV.mageres_comp_sigma(o, 0.0) == 0.5
    ws = 10 .^ range(-9, 0; length = 50)
    @test all(diff([MV.mageres_comp_sigma(o, w) for w in ws]) .<= 0)
    c = cv(comp_var_law = "constant", comp_var_ref_pct = 7.0)
    @test all(MV.mageres_comp_sigma(c, w) == 0.07 for w in ws)
    @test_throws ArgumentError cv(comp_var_law = "other")
    @test_throws ArgumentError cv(comp_var_ref_pct = -1.0)
    @test_throws ArgumentError cv(comp_var_host_corr_m = 0.0)
    @test_throws ArgumentError cv(comp_var_seed = -1)
end

@testset "lognormal draws" begin
    o  = cv(comp_var_injection = true)
    x0 = MV.mageres_injection_base_bulk(o)
    C0 = collect(o.te_injection_ppm)
    xs = Vector{Vector{Float64}}(undef, 20_000); Cs = similar(xs)
    for k in eachindex(xs)
        xs[k], Cs[k] = MV.mageres_injection_composition(o, k, C0)
    end
    X = reduce(hcat, xs); C = reduce(hcat, Cs)
    @test all(sum(X; dims = 1) .≈ sum(x0))
    @test all(X[x0 .== 0, :] .== 0)
    σx, σc = MV.mageres_comp_sigmas(o, o.oxides, x0, C0)
    for e in (findfirst(==("Lu"), o.te_elements), findfirst(==("La"), o.te_elements), findfirst(==("Zr"), o.te_elements))
        @test mean(C[e, :]) ≈ C0[e] rtol = 0.01
        @test std(C[e, :]) / C0[e] ≈ σc[e] rtol = 0.05
    end
    i_ti = findfirst(==("TiO2"), o.oxides)
    @test mean(X[i_ti, :]) ≈ x0[i_ti] rtol = 0.01
    @test std(X[i_ti, :]) / x0[i_ti] ≈ σx[i_ti] rtol = 0.1
    @test MV.mageres_injection_composition(o, 7, C0) == MV.mageres_injection_composition(o, 7, C0)
    @test MV.mageres_injection_composition(o, 7, C0) != MV.mageres_injection_composition(o, 8, C0)
    @test first(MV.mageres_injection_composition(o, 7, C0)) == first(MV.mageres_injection_composition(o, 7, Float64[]))
    @test MV.mageres_injection_composition(cv(comp_var_injection = true, comp_var_seed = 2), 7, C0) !=
          MV.mageres_injection_composition(o, 7, C0)
    off = cv()
    @test MV.mageres_injection_composition(off, 7, C0) == (x0, C0)
end

@testset "layered host" begin
    o  = cv(comp_var_host = true, comp_var_host_corr_m = 1_000.0, domain_thickness_km = 60.0, grid_h_max_m = 1000.0)
    y  = MV.mageres_host_profile_y(o)
    @test maximum(diff(y)) <= o.comp_var_host_corr_m / 4 + 1e-9
    @test length(MV.mageres_host_profile_y(cv(domain_thickness_km = 60.0))) <= length(y)
    xs, Cs = MV.mageres_host_compositions(o, y, o.te_host_ppm)
    x0 = MV.mageres_normalized(o.host_bulk)
    @test length(xs) == length(y) && all(sum(x) ≈ sum(x0) for x in xs)
    i_si = findfirst(==("SiO2"), o.oxides)
    s  = [x[i_si] / x0[i_si] - 1 for x in xs]
    r1 = cor(s[1:end-1], s[2:end])
    @test r1 ≈ exp(-(y[2] - y[1]) / o.comp_var_host_corr_m) atol = 0.08
    @test xs[1] != xs[end]
    @test MV.mageres_host_bulk_at(o, y[3]) ≈ xs[3]
    @test MV.mageres_host_bulk_at(o, (y[3] + y[4]) / 2) ≈ (xs[3] .+ xs[4]) ./ 2
    off = cv()
    xo, Co = MV.mageres_host_compositions(off, MV.mageres_host_profile_y(off), off.te_host_ppm)
    @test all(x == x0 for x in xo) && all(C == off.te_host_ppm for C in Co)
    @test MV.mageres_host_material(off, -5e3).bulk == MV.mageres_host_material(off).bulk
end

@testset "preview figure" begin
    @test MV.mageres_comp_var_figure(cv()) !== nothing
    @test MV.mageres_comp_var_figure(cv(comp_var_host = true, comp_var_injection = true, comp_var_law = "constant")) !== nothing
end

@testset "constant-property run with variability" begin
    o = cv(thermodynamics = "constant", te_enabled = false, eruption_trigger = "melt", comp_var_host = true,
           comp_var_injection = true, run_time_yr = 5_000.0, injection_period_yr = 1_000.0, output_every_yr = 5_000.0)
    st = MV.mageres_init_state(o)
    MV.mageres_run!(st, o.run_time_yr)
    L  = st.source.layers
    @test length(L) >= 4
    b  = [l.moles ./ sum(l.moles) for l in L]
    @test all(b[k] != b[k+1] for k in 1:length(b)-1)
    @test b[2] ≈ first(MV.mageres_injection_composition(o, 2, Float64[])) rtol = 1e-6
    @test abs(MV.mageres_budget(st).mass_error) < 1e-8
end

@testset "MAGEMin run with variability" begin
    o  = cv(comp_var_host = true, comp_var_injection = true, run_time_yr = 3_000.0, output_every_yr = 1_000.0)
    st = MV.mageres_init_state(o; progressbar = false)
    pr = st.thermo.profile
    ox = MV.mageres_ox_rows(st)
    @test pr.density[ox, 1] ./ sum(pr.density[ox, 1]) != pr.density[ox, end] ./ sum(pr.density[ox, end])
    MV.mageres_run!(st, o.run_time_yr)
    L  = st.source.layers
    @test length(L) >= 3
    b  = [l.moles ./ sum(l.moles) for l in L]
    @test all(b[k] != b[k+1] for k in 1:length(b)-1)
    B  = MV.mageres_budget(st)
    @test abs(B.mass_error) < 1e-6
    @test abs(B.budget_error) < 1e-3
end

end
