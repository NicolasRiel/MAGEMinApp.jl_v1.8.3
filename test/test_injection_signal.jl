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
using MAGEMinApp

const MI = MAGEMinApp

if !isdefined(MI, :mageres_injection_T)
    @info "MAGERes disabled in this build: injection-signal tests skipped"
else

signal(; kw...) = (o = MI.MAGEResOptions(; T_injection_mode = "signal", kw...); MI.mageres_validate_options(o); o)

@testset "constant mode" begin
    o = MI.MAGEResOptions(T_injection_C = 1012.5)
    @test MI.mageres_injection_T(o, 0.0) == 1012.5
    @test MI.mageres_injection_T(o, 7.3e4) == 1012.5
    @test MI.mageres_injection_T_max(o) == 1012.5
    signal_zero = signal(T_injection_C = 1012.5)
    @test all(MI.mageres_injection_T_components(signal_zero, collect(0.0:1000.0:1e5)).total .== 1012.5)
end

@testset "first order: Gaussian peak" begin
    o = signal(T_inj_peak_C = 80.0, T_inj_peak_yr = 30_000.0, T_inj_peak_rise_yr = 2_000.0,
               T_inj_peak_decay_yr = 8_000.0)
    @test MI.mageres_injection_T(o, 30_000.0) ≈ 1030.0
    @test MI.mageres_injection_T(o, 28_000.0) - 950.0 ≈ 80.0 * exp(-0.5)
    @test MI.mageres_injection_T(o, 38_000.0) - 950.0 ≈ 80.0 * exp(-0.5)
    @test MI.mageres_injection_T(o, 0.0) - 950.0 < 1e-6
    T = MI.mageres_injection_T_components(o, collect(0.0:100.0:1e5)).total
    k = argmax(T)
    @test all(diff(T[1:k]) .>= 0) && all(diff(T[k:end]) .<= 0)
    s = signal(T_inj_peak_C = 50.0, T_inj_peak_yr = 40_000.0)
    @test MI.mageres_injection_T(s, 35_000.0) ≈ MI.mageres_injection_T(s, 45_000.0)
    c = signal(T_inj_peak_C = -60.0, T_inj_peak_yr = 30_000.0)
    @test MI.mageres_injection_T(c, 30_000.0) ≈ 890.0
end

@testset "first order: repeated cycles" begin
    one = signal(T_inj_peak_C = 80.0, T_inj_peak_yr = 20_000.0, T_inj_peak_rise_yr = 3_000.0, T_inj_peak_decay_yr = 6_000.0)
    rep = signal(T_inj_peak_C = 80.0, T_inj_peak_yr = 20_000.0, T_inj_peak_rise_yr = 3_000.0, T_inj_peak_decay_yr = 6_000.0,
                 T_inj_peak_cycles = 3, T_inj_peak_interval_yr = 30_000.0)
    ts = collect(0.0:500.0:1.2e5)
    @test MI.mageres_injection_T_components(one, ts).total ==
          MI.mageres_injection_T_components(signal(T_inj_peak_C = 80.0, T_inj_peak_yr = 20_000.0,
              T_inj_peak_rise_yr = 3_000.0, T_inj_peak_decay_yr = 6_000.0, T_inj_peak_interval_yr = 5_000.0), ts).total
    for k in 0:2
        @test MI.mageres_injection_T(rep, 20_000.0 + 30_000.0k) ≈ 1030.0
        @test MI.mageres_injection_peak(rep, 20_000.0 + 30_000.0k + 4_000.0) ≈ MI.mageres_injection_peak(one, 24_000.0)
    end
    @test MI.mageres_injection_peak(rep, 110_000.0 + 20_000.0) < 1e-6
    tight = signal(T_inj_peak_C = 80.0, T_inj_peak_yr = 20_000.0, T_inj_peak_cycles = 4, T_inj_peak_interval_yr = 5_000.0)
    @test maximum(MI.mageres_injection_T_components(tight, ts).peak) ≈ 80.0
    @test_throws ArgumentError signal(T_inj_peak_cycles = 0)
    @test_throws ArgumentError signal(T_inj_peak_cycles = 2, T_inj_peak_interval_yr = 0.0)
    d = MI.mageres_options_to_dict(rep)
    d["T_inj_peak_cycles"] = 3.0
    @test MI.mageres_options_from_dict(d).T_inj_peak_cycles == 3
end

@testset "second order: sinusoid" begin
    o = signal(T_inj_sine_C = 15.0, T_inj_sine_period_yr = 4_000.0)
    @test MI.mageres_injection_T(o, 0.0) ≈ 950.0
    @test MI.mageres_injection_T(o, 1_000.0) ≈ 965.0
    @test MI.mageres_injection_T(o, 3_000.0) ≈ 935.0
    @test MI.mageres_injection_T(o, 41_000.0) ≈ 965.0
    p = signal(T_inj_sine_C = 15.0, T_inj_sine_period_yr = 4_000.0, T_inj_sine_phase_deg = 90.0)
    @test MI.mageres_injection_T(p, 0.0) ≈ 965.0
    both = signal(T_inj_sine_C = 15.0, T_inj_sine_period_yr = 4_000.0, T_inj_peak_C = 80.0, T_inj_peak_yr = 20_000.0)
    c = MI.mageres_injection_T_components(both, [21_000.0])
    @test c.total[1] ≈ 950.0 + c.peak[1] + 15.0
end

@testset "third order: perturbation" begin
    o  = signal(T_inj_noise_C = 20.0, T_inj_noise_corr_yr = 1_000.0, T_inj_noise_seed = 7)
    ts = collect(0.0:250.0:2e6)
    n1 = MI.mageres_injection_noise(o, ts)
    @test n1 == MI.mageres_injection_noise(o, ts)
    @test n1 != MI.mageres_injection_noise(signal(T_inj_noise_C = 20.0, T_inj_noise_corr_yr = 1_000.0,
                                                  T_inj_noise_seed = 8), ts)
    sd = sqrt(sum(n1 .^ 2) / length(n1))
    @test 0.85 * 20.0 < sd < 1.15 * 20.0
    r1 = sum(n1[1:end-1] .* n1[2:end]) / sum(n1 .^ 2)
    @test abs(r1 - exp(-250.0 / 1000.0)) < 0.05
    @test MI.mageres_injection_noise(o, [12_345.0])[1] ≈ MI.mageres_injection_noise(o, [0.0, 12_345.0, 9e5])[2]
end

@testset "clamp and validation" begin
    o = signal(T_inj_peak_C = 600.0, T_inj_peak_yr = 10_000.0, T_inj_max_C = 1300.0)
    c = MI.mageres_injection_T_components(o, [10_000.0])
    @test c.total[1] == 1300.0
    @test c.clipped[1]
    @test_throws ArgumentError signal(T_inj_min_C = 1000.0, T_inj_max_C = 900.0)
    @test_throws ArgumentError signal(T_inj_peak_rise_yr = 0.0)
    @test_throws ArgumentError signal(T_inj_sine_period_yr = 0.0)
    @test_throws ArgumentError signal(T_inj_noise_corr_yr = -1.0)
    @test_throws ArgumentError MI.mageres_validate_options(MI.MAGEResOptions(T_injection_mode = "ramp"))
end

@testset "options round trip and old configs" begin
    o = signal(T_inj_peak_C = 50.0, T_inj_sine_C = 10.0, T_inj_noise_C = 5.0, T_inj_noise_seed = 3)
    r = MI.mageres_options_from_dict(MI.mageres_options_to_dict(o))
    @test MI.mageres_injection_T_components(r, collect(0.0:1000.0:1e5)).total ==
          MI.mageres_injection_T_components(o, collect(0.0:1000.0:1e5)).total
    d = MI.mageres_options_to_dict(MI.MAGEResOptions())
    for k in collect(keys(d))
        (startswith(k, "T_inj_") || k == "T_injection_mode") && delete!(d, k)
    end
    @test MI.mageres_options_from_dict(d).T_injection_mode == "constant"
    d["T_inj_trend_C"] = 10.0
    @test MI.mageres_options_from_dict(d).T_injection_mode == "constant"
    @test MI.mageres_injection_T_figure(o) !== nothing
    @test MI.mageres_injection_T_figure(MI.MAGEResOptions()) !== nothing
    old = (T_injection_C = 900.0, T_injection_mode = "signal", T_inj_peak_C = 40.0, T_inj_trend_C = 5.0,
           oxides = Any["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "H2O"])
    r = MI.mageres_options_from_stored(old)
    @test r.T_inj_peak_C == 40.0 && r.T_inj_peak_cycles == 1 && r.T_injection_mode == "signal"
    @test r.oxides isa Vector{String}
end

@testset "constant-property run follows the signal" begin
    o = signal(thermodynamics = "constant", te_enabled = false, eruption_trigger = "melt",
               run_time_yr = 6_000.0, injection_period_yr = 1_000.0, output_every_yr = 6_000.0,
               T_inj_peak_C = 120.0, T_inj_peak_yr = 6_000.0, T_inj_peak_rise_yr = 2_000.0,
               T_inj_sine_C = 10.0, T_inj_sine_period_yr = 2_000.0, T_inj_noise_C = 10.0)
    st = MI.mageres_init_state(o)
    MI.mageres_run!(st, o.run_time_yr)
    @test length(st.source.layers) >= 5
    @test all(l.T ≈ MI.mageres_injection_T(o, l.t) for l in st.source.layers)
    @test st.source.layers[end].T > st.source.layers[1].T
end

@testset "equilibrated-melt injection: trace elements" begin
    o  = MI.MAGEResOptions(injection_melt_only = true, T_injection_mode = "signal")
    mm = MI.mageres_init_magemin(o)
    th = MI.mageres_build_thermo(o, mm; progressbar = false)
    te = th.te
    P  = MI.mageres_lithostatic_kbar(o, -1e3 * o.lens_depth_km)
    x  = MI.mageres_normalized(isempty(o.injection_bulk) ? o.magma_bulk : o.injection_bulk)
    for T in (880.0, 1000.0)
        pb = MI.mageres_equilibrate!(th, o, [x], [P], [T]; stage = "test")[1][1]
        @test 0 < pb.te.w_M < 1
        Cm = MI.mageres_injection_melt_te(pb, te, te.C_injection)
        for e in eachindex(Cm)
            e == te.i_zr && te.zircon && continue
            @test pb.te.w_M * Cm[e] + pb.te.w_S * pb.te.D[e] * Cm[e] ≈ te.C_injection[e] rtol = 1e-8
        end
        te.zircon && isfinite(pb.te.C_sat) && @test Cm[te.i_zr] <= pb.te.C_sat + 1e-9
        MI.mageres_injection_equilibrium!(th, o, T)
        rho = MI.mageres_density_mass(th.injection_density, th.Mox)
        r   = MI.MAGERES_IN0 + length(o.oxides) .+ (0:length(te.elements)-1)
        @test th.injection_density[r] ≈ rho .* Cm .* 1e-6 rtol = 1e-6
        @test th.injection.frac_M > 0.8
    end
    phi880  = MI.mageres_injection_equilibrium!(th, o, 880.0)
    phi1000 = MI.mageres_injection_equilibrium!(th, o, 1000.0)
    @test 0 < phi880 < phi1000 <= 1
    @test MI.mageres_injection_equilibrium!(th, o, 880.0) == phi880
    @test MI.mageres_injection_volume(o, phi880) ≈ phi880 * MI.mageres_injection_volume(o)
    bulk = MI.MAGEResOptions(injection_melt_only = false)
    @test MI.mageres_injection_volume(bulk, phi880) == MI.mageres_injection_volume(bulk)
    i_rb = findfirst(==("Rb"), te.elements)
    d880 = (MI.mageres_injection_equilibrium!(th, o, 880.0); copy(th.injection_density))
    d1000 = (MI.mageres_injection_equilibrium!(th, o, 1000.0); copy(th.injection_density))
    r = MI.MAGERES_IN0 + length(o.oxides) .+ (0:length(te.elements)-1)
    @test d880[r][i_rb] / MI.mageres_density_mass(d880, th.Mox) > d1000[r][i_rb] / MI.mageres_density_mass(d1000, th.Mox)
    MI.mageres_injection_equilibrium!(th, o, 880.0)
    @test th.injection_density == d880
end

end
