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

const MN = MAGEMinApp

const ZN_OPTS = MN.MAGEResOptions(run_time_yr = 3_000.0, output_every_yr = 1_000.0)
const ZN_ST   = (st = MN.mageres_init_state(ZN_OPTS; progressbar = false); MN.mageres_run!(st, ZN_OPTS.run_time_yr); st)

zn_total(st, c) = (te = MN.mageres_te(st);
                   st.Q[MN.mageres_te_rows(st)[te.i_zr], c] +
                   MN.MAGERES_ZRC_ZR_FRAC * sum(sum(st.Q[MN.mageres_zrc_mass_rows(st, f), c]) for f in 1:te.n_fam))

@testset "zircon nucleation with existing grains" begin
    st0 = ZN_ST
    te  = MN.mageres_te(st0)
    c   = something(findfirst(st0.body), findfirst(k -> st0.grid.frac[k] >= 0.5 && st0.thermo.points[k].te.w_M > 0,
                                                 1:MN.mageres_ncells(st0.grid)))
    T    = MN.mageres_temperature(st0)[c]
    p    = st0.thermo.points[c]
    q    = p.te
    mc   = MN.mageres_cell_mass(st0, c)
    meff = mc * (q.w_M + q.D[te.i_zr] * q.w_S)
    Cs   = 1e-6 * MN.mageres_zr_saturation(p, T)
    iz   = MN.mageres_te_rows(st0)[te.i_zr]
    fn   = te.n_fam
    for (S, expect) in ((2.0, true), (1.05, false))
        st = deepcopy(st0)
        for f in 1:te.n_fam
            st.Q[MN.mageres_zrc_N_row(st, f), c] = 0.0
            st.Q[MN.mageres_zrc_mass_rows(st, f), c] .= 0.0
        end
        st.Q[MN.mageres_zrc_N_row(st, 1), c]        = 1e3
        st.Q[MN.mageres_zrc_mass_rows(st, 1)[1], c] = 1e3 * MN.MAGERES_ZRC_RHO * 4 / 3 * π * (100e-6)^3
        st.Q[iz, c] = S * Cs * meff
        z0 = zn_total(st, c)
        MN.mageres_zircon_cell!(st, c, T, 1e3 * MN.MAGERES_SECONDS_PER_YEAR, 2, fn, false)
        @test (st.Q[MN.mageres_zrc_N_row(st, fn), c] > 0) == expect
        @test zn_total(st, c) ≈ z0 rtol = 1e-12
        @test st.Q[MN.mageres_zrc_N_row(st, 1), c] == 1e3
    end
end
