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

const MZ = MAGEMinApp

const ZA_OPTS = MZ.MAGEResOptions(run_time_yr = 3_000.0, output_every_yr = 1_000.0)
const ZA_ST   = (st = MZ.mageres_init_state(ZA_OPTS; progressbar = false); MZ.mageres_run!(st, ZA_OPTS.run_time_yr); st)

zr_total(st, c) = (te = MZ.mageres_te(st);
                   st.Q[MZ.mageres_te_rows(st)[te.i_zr], c] +
                   MZ.MAGERES_ZRC_ZR_FRAC * sum(sum(st.Q[MZ.mageres_zrc_mass_rows(st, f), c]) for f in 1:te.n_fam))

@testset "host-origin row" begin
    st = ZA_ST
    r  = MZ.mageres_host_row(st)
    e  = MZ.mageres_expected_content(st.ledger)
    @test sum(st.Q[r, :]) ≈ e[r] rtol = 1e-10
    @test st.ledger.injected[r] == 0.0
    @test all(st.Q[r, c] <= MZ.mageres_cell_mass(st, c) * (1 + 1e-9) for c in 1:MZ.mageres_ncells(st.grid))
    g   = st.grid
    tip = 1e3 * ZA_OPTS.lens_length_km / 2 + ZA_OPTS.grid_h_max_m + 300.0
    far = [c for c in 1:MZ.mageres_ncells(g)
           if abs(MZ.mageres_frame_coords(st.chamber.F, MZ.MAGEResPoint(g.xc[c], g.yc[c]))[1]) > tip]
    @test !isempty(far)
    @test all(st.Q[r, c] ≈ MZ.mageres_cell_mass(st, c) for c in far)
    b = MZ.mageres_budget(st)
    @test b.mass_error < 1e-10
    @test b.te_error < 1e-10
    @test all(isnan(MZ.mageres_layer_host_mass(l, MZ.mageres_te(st))) == false for l in st.source.layers)
    @test all(MZ.mageres_layer_host_mass(l, MZ.mageres_te(st)) == 0.0 for l in st.source.layers)
    d = MZ.mageres_volume_diagnostics(st)
    @test 0 <= d.host_body < 1
    tf, nf = MZ.mageres_feeding_series([MZ.MAGEResOutput(st.t, "end", st)])
    @test tf == [st.t]
    @test nf == [count(==(6.0), MZ.mageres_material_class(st))]
end

@testset "host zircon release law" begin
    o = ZA_OPTS
    @test MZ.mageres_disaggregation_release(o, o.phi_perc) == 0.0
    @test MZ.mageres_disaggregation_release(o, 0.5 * (o.phi_perc + o.phi_disaggregation)) ≈ 0.5
    @test MZ.mageres_disaggregation_release(o, o.phi_disaggregation) == 1.0
    @test MZ.mageres_disaggregation_release(o, 0.9) == 1.0
end

@testset "disaggregation exchange" begin
    st   = deepcopy(ZA_ST)
    g    = st.grid
    b    = findall(st.body)
    r    = MZ.mageres_host_row(st)
    host = [c for c in 1:MZ.mageres_ncells(g) if g.frac[c] == 0 && !st.body[c]]
    c    = host[argmin([minimum(hypot(g.xc[h] - g.xc[k], g.yc[h] - g.yc[k]) for k in b) for h in host])]
    conn = falses(MZ.mageres_ncells(g)); conn[c] = true; conn[b] .= true
    melt = zeros(MZ.mageres_ncells(g)); melt[c] = ZA_OPTS.phi_disaggregation + 0.05
    fill!(st.melt_out, 0.0)
    tot0 = vec(sum(st.Q; dims = 2))
    Xc   = copy(st.Q[:, c])
    hb0  = sum(st.Q[r, b])
    zr   = MZ.mageres_zircon_rows(st)
    zb0  = sum(st.Q[zr, b])
    done = MZ.mageres_disaggregate!(st, conn, melt)
    @test done == [c]
    @test vec(sum(st.Q; dims = 2)) ≈ tot0 rtol = 1e-12
    @test sum(st.Q[r, b]) > hb0
    @test sum(st.Q[zr, b]) - zb0 ≈ sum(Xc[zr]) - sum(st.Q[zr, c]) rtol = 1e-9
    @test st.Q[r, c] < 0.5 * MZ.mageres_cell_mass(st, c)
    @test MZ.mageres_disaggregate!(st, conn, melt) == Int[]
    cold = deepcopy(ZA_ST)
    @test MZ.mageres_disaggregate!(cold, conn, zeros(MZ.mageres_ncells(g))) == Int[]
end
