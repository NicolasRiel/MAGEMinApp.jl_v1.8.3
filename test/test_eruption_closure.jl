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

const MC = MAGEMinApp

function closure_setup(; opening = "split")
    o  = MC.MAGEResOptions(thermodynamics = "constant", te_enabled = false, opening_mode = opening,
                           eruption_trigger = "melt", run_time_yr = 10_000.0, injection_period_yr = 1_000.0,
                           grid_h_min_m = 62.5, lens_thickness_m = 100.0)
    st = MC.mageres_init_state(o; progressbar = false)
    for _ in 1:8
        MC.mageres_step!(st)
    end
    return o, st
end

function indicator_remap(o, st, kin)
    old = st.grid
    n   = MC.mageres_ncells(old)
    Qi  = zeros(n + 1, n)
    for m in 1:n
        Qi[m, m] = old.area[m]
    end
    nc = MC.mageres_kin_profile(st.chamber, kin)
    g  = MC.mageres_build_grid(o, nc)
    Qn, led = MC.mageres_remap(old, Qi, st.chamber, kin, g; mode = :compact,
                               injected_density = vcat(zeros(n), 1.0), inflow_density = y -> zeros(n + 1))
    w = vec(sum(Qn[1:n, :]; dims = 2)) ./ old.area
    return w, Qn, g, nc, led
end

wavy(n, a, b) = [a + b * sin(1.7 * i) for i in 1:n]
swavy(c, a, b) = MC.mageres_mirror_stations!(wavy(length(c.S), a, b), c.S)

@testset "event kinematics carry every chamber cell exactly once" begin
    for opening in ("split", "up", "down")
        o, st = closure_setup(; opening = opening)
        c     = st.chamber
        n     = length(c.S)
        f, r  = MC.mageres_station_heights(c)
        h     = r .- f
        ch    = [m for m in 1:MC.mageres_ncells(st.grid) if st.grid.frac[m] > 0]
        q1    = swavy(c, 0.3, 0.1)
        q2    = swavy(c, 0.4, 0.1)
        kins  = Dict(
            "identity"      => MC.mageres_kin_identity(c),
            "uniform close" => MC.mageres_kin_close(c, h .* swavy(c, 0.12, 0.08)),
            "zoned close"   => MC.mageres_kin_close(c, h .* swavy(c, 0.10, 0.05);
                                                    b1 = f .+ h .* swavy(c, 0.3, 0.1), b2 = f .+ h .* swavy(c, 0.6, 0.1),
                                                    q = hcat(q1, q2, 1 .- q1 .- q2)),
            "opening"       => MC.mageres_kin_close(c, -h .* swavy(c, 0.10, 0.05)),
            "injection"     => MC.mageres_kin_inject(c, MC.mageres_lens_delta(o, c, 4e5), f .+ h .* swavy(c, 0.5, 0.3)))
        for (label, kin) in kins
            w, Qn, g, nc, led = indicator_remap(o, st, kin)
            @test maximum(abs(w[m] - 1) for m in ch) < 1e-9
            @test MC.mageres_chamber_area(nc) - MC.mageres_chamber_area(c) ≈ MC.mageres_kin_opened_area(c, kin) rtol = 1e-9
            @test nc.T .* nc.shape ≈ kin.new[:, end] .- kin.new[:, 1] atol = 1e-9
        end
        _, Qn, _, _, _ = indicator_remap(o, st, kins["injection"])
        @test sum(Qn[end, :]) ≈ 4e5 / MC.MAGERES_SYMMETRY rtol = 1e-9
    end
end

@testset "lens-shaped closure compresses uniformly" begin
    o, st = closure_setup()
    c    = st.chamber
    f, r = MC.mageres_station_heights(c)
    kin  = MC.mageres_kin_close(c, 0.2 .* (r .- f))
    nc   = MC.mageres_kin_profile(c, kin)
    g    = MC.mageres_build_grid(o, nc)
    Qn, _ = MC.mageres_remap(st.grid, permutedims(copy(st.grid.area)), c, kin, g; mode = :compact,
                             injected_density = [1.0], inflow_density = y -> [1.0])
    inner = [m for m in 1:MC.mageres_ncells(g) if g.frac[m] >= 1 - 1e-12]
    @test !isempty(inner)
    @test all(Qn[1, m] / g.area[m] ≈ 1 / 0.8 for m in inner)
end

@testset "sill opens at the requested heights" begin
    o, st = closure_setup()
    c    = st.chamber
    n    = length(c.S)
    f, r = MC.mageres_station_heights(c)
    kin  = MC.mageres_kin_inject(c, MC.mageres_lens_delta(o, c, 6e5), f .+ (r .- f) .* wavy(n, 0.5, 0.3))
    _, Qn, g, nc, _ = indicator_remap(o, st, kin)
    hits = 0
    for j in 1:MC.mageres_ncells(g)
        Qn[end, j] / g.area[j] > 1 - 1e-9 || continue
        s, ξ = MC.mageres_frame_coords(nc.F, MC.MAGEResPoint(g.xc[j], g.yc[j]))
        k  = clamp(searchsortedlast(c.S, s), 1, n - 1)
        t  = (s - c.S[k]) / (c.S[k+1] - c.S[k])
        lo = (1 - t) * kin.new[k, 2] + t * kin.new[k+1, 2]
        hi = (1 - t) * kin.new[k, 3] + t * kin.new[k+1, 3]
        @test lo - g.h[j] / 2 - 1e-6 <= ξ <= hi + g.h[j] / 2 + 1e-6
        hits += 1
    end
    @test hits > 0
end

@testset "cut cells keep host rock out of the chamber" begin
    o  = MC.MAGEResOptions(thermodynamics = "constant", te_enabled = false, eruption_trigger = "melt",
                           run_time_yr = 10_000.0, injection_period_yr = 1_000.0)
    st = MC.mageres_init_state(o; progressbar = false)
    for _ in 1:8
        MC.mageres_step!(st)
    end
    g    = st.grid
    rows = MC.mageres_ox_rows(st)
    xm   = MC.mageres_normalized(isempty(o.injection_bulk) ? o.magma_bulk : o.injection_bulk)
    xh   = MC.mageres_normalized(o.host_bulk)
    inner = [n for n in 1:MC.mageres_ncells(g) if g.frac[n] >= 1 - 1e-9]
    outer = [n for n in 1:MC.mageres_ncells(g) if g.frac[n] <= 1e-9]
    @test !isempty(inner)
    dev(n, x) = (v = st.Q[rows, n]; maximum(abs.(v ./ sum(v) .- x ./ sum(x))))
    @test maximum(dev(n, xm) for n in inner) < 1e-10
    @test maximum(dev(n, xh) for n in outer) < 1e-10
    @test vec(sum(st.Q; dims = 2)) ≈ MC.mageres_expected_content(st.ledger) rtol = 1e-10
end

@testset "station closure is exact for every pair of segments" begin
    S = collect(range(-2500.0, 2500.0; length = 41))
    D = [100.0 * exp(-((0.5 * (S[k] + S[k+1])) / 400)^2) + (isodd(k) ? 3.0 : -2.0) for k in 1:40]
    D = 0.5 .* (D .+ reverse(D))
    x = MC.mageres_station_closure(S, D)
    @test x[1] == 0 && x[end] == 0
    for b in 1:20
        k = 2b - 1
        @test (S[k+1] - S[k]) * (x[k] + x[k+1]) / 2 + (S[k+2] - S[k+1]) * (x[k+1] + x[k+2]) / 2 ≈ D[k] + D[k+1] rtol = 1e-12
    end
    @test x ≈ reverse(x)
end

@testset "chamber profile stays symmetric" begin
    o, st = closure_setup()
    c    = st.chamber
    n    = length(c.S)
    f, r = MC.mageres_station_heights(c)
    kin  = MC.mageres_kin_close(c, (r .- f) .* wavy(n, 0.12, 0.08))
    nc   = MC.mageres_kin_profile(c, kin)
    @test nc.shape ≈ reverse(nc.shape)
    left = [i for i in 1:n if c.S[i] <= 0]
    @test MC.mageres_chamber_area(nc) ≈ 2 * MC.mageres_trapz(c.S[left], nc.T .* nc.shape[left]) rtol = 1e-12
end
