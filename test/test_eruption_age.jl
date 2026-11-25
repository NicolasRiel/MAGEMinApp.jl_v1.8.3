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
using Random
using MAGEMinApp

const MA = MAGEMinApp

if !isdefined(MA, :mageres_bayesian_eruption_age)
    @info "MAGERes disabled in this build: eruption-age tests skipped"
else

@testset "error function" begin
    @test MA.mageres_erfc(0.0) ≈ 1.0 rtol = 1e-6
    @test MA.mageres_erfc(1.0) ≈ 0.157299207050285 rtol = 1e-6
    @test MA.mageres_erfc(3.0) ≈ 2.209049699858544e-5 rtol = 1e-6
    @test MA.mageres_erfc(-1.0) ≈ 2 - 0.157299207050285 rtol = 1e-6
    @test MA.mageres_normal_mass(-10.0, 10.0) ≈ 1.0 rtol = 1e-6
    @test MA.mageres_normal_mass(5.0, 6.0) ≈ 2.866515718791939e-7 - 9.865876450377e-10 rtol = 1e-5
end

@testset "youngest cluster" begin
    @test all(isnan, MA.mageres_youngest_cluster([1.0, 2.0], [0.1, 0.1]))
    mu, se = MA.mageres_youngest_cluster([5.0, 5.0, 5.0, 50.0], [1.0, 1.0, 1.0, 1.0])
    @test mu ≈ 5.0
    @test se ≈ 1 / sqrt(3)
    @test MA.mageres_youngest_cluster([1.0, 10.0, 10.1, 10.2], [0.1, 0.1, 0.1, 0.1])[1] ≈ 10.1
end

@testset "eruption-age estimators" begin
    rng  = MersenneTwister(11)
    bias = Dict(m => Float64[] for (m, _) in MA.MAGERES_AGE_ESTIMATORS)
    inside = 0
    nrep   = 30
    for r in 1:nrep
        x = 10 .+ 40 .* rand(rng, 100)
        s = fill(3.0, 100)
        y = x .+ s .* randn(rng, 100)
        for (m, _) in MA.MAGERES_AGE_ESTIMATORS
            e, es = MA.mageres_eruption_age_estimate(m, y, s; seed = r)
            push!(bias[m], e - 10)
            m == "bayes" && abs(e - 10) <= 2es && (inside += 1)
        end
    end
    mean(v) = sum(v) / length(v)
    @test abs(mean(bias["bayes"])) < 1.0
    @test inside / nrep >= 0.8
    @test mean(bias["ysg"]) < -1.0
    @test mean(bias["wm"]) > 15.0
    @test MA.mageres_eruption_age_estimate("bayes", [4.0], [1.0]) == (4.0, 1.0)
    @test all(isnan, MA.mageres_eruption_age_estimate("bayes", Float64[], Float64[]))
    @test_throws ArgumentError MA.mageres_eruption_age_estimate("none", [1.0], [1.0])
end

seq_sample(t, mass, grains; inh = 0.0) =
    MA.MAGEResZirconSample(t, mass, zeros(length(mass)), inh, 1.0, grains, fill(1.0, length(grains)), 0,
                           [1.0], [vcat(inh, mass)], MA.MAGEResZirconGrain[], t, 1)

@testset "sequence statistics" begin
    s = seq_sample(10_000.0, [0.0, 0.0, 1.0, 0.0, 0.0], [7.2, 7.5, 7.8]; inh = 1.0)
    e = MA.MAGEResSeqEntry(s, UInt64(1), 1.0, 800.0, 2_500.0, false)
    q = MA.mageres_seq_stats(e, 1000.0, 10_000.0)
    @test q.erupt_age ≈ 0.0
    @test q.mean ≈ 7.5
    @test q.q05 ≈ 7.05
    @test q.q95 ≈ 7.95
    @test q.youngest ≈ 7.0
    @test q.residence ≈ 7.5
    @test q.sd ≈ 1 / sqrt(12)
    @test abs(q.skew) < 1e-9
    @test q.kurt ≈ -1.2
    @test q.antecryst ≈ 0.5
    @test q.inherited ≈ 0.5
    @test q.n_spots == 3
    first = MA.MAGEResSeqEntry(s, UInt64(1), 1.0, 800.0, -Inf, false)
    @test MA.mageres_seq_stats(first, 1000.0, 10_000.0).antecryst == 0.0
end

@testset "sequence figures and table" begin
    rng = MersenneTwister(3)
    E   = [MA.MAGEResSeqEntry(seq_sample(t, [1.0, 2.0, 3.0, 1.0, 0.0], 1 .+ 3 .* rand(rng, 30) .+ (10 - t / 1e3)),
                              hash(t), 1e6, 850.0, t - 2000, t > 6000) for t in (2_500.0, 4_000.0, 8_000.0)]
    te  = 1000.0
    seq = MA.MAGEResSequence(E, [MA.mageres_seq_stats(e, te, 10_000.0) for e in E],
                             MA.mageres_seq_estimates(E, 10_000.0), 10_000.0, MA.MAGEResOptions())
    @test size(seq.est) == (3, length(MA.MAGERES_AGE_ESTIMATORS))
    for f in (MA.mageres_zircon_sequence_figure, MA.mageres_zircon_sequence_stats_figure,
              MA.mageres_eruption_age_test_figure)
        @test f(seq) !== nothing
        @test f(nothing) !== nothing
    end
    Ns, med, cv = MA.mageres_seq_grain_test(E, 10_000.0)
    @test Ns == [3, 5, 10, 20]
    @test size(med) == (4, length(MA.MAGERES_AGE_ESTIMATORS))
    csv = MA.mageres_sequence_csv(seq)
    lines = split(strip(csv), "\n")
    @test length(lines) == 2 + length(E)
    @test length(split(lines[2], ",")) == length(split(lines[3], ","))
end

@testset "host cells feeding the chamber" begin
    t  = collect(0.0:1_000.0:10_000.0)
    n  = [0, 0, 2, 5, 9, 12, 10, 7, 7, 4, 1]
    tr = MA.mageres_feeding_trace((t, n), 10_000.0, 0.0, 10.0; yaxis = "y7")
    @test tr.fields[:y] == reverse(n)
    @test tr.fields[:x] == reverse(10.0 .- t ./ 1e3)
    @test MA.mageres_feeding_trace((t, n), 10_000.0, 20.0, 30.0; yaxis = "y7") === nothing
    @test MA.mageres_feeding_trace(nothing, 10_000.0, 0.0, 10.0; yaxis = "y7") === nothing
    rng  = MersenneTwister(3)
    E    = [MA.MAGEResSeqEntry(seq_sample(t, [1.0, 2.0, 3.0, 1.0, 0.0], 1 .+ 3 .* rand(rng, 30) .+ (10 - t / 1e3)),
                               hash(t), 1e6, 850.0, t - 2000, t > 6000) for t in (2_500.0, 4_000.0, 8_000.0)]
    seq2 = MA.MAGEResSequence(E, [MA.mageres_seq_stats(e, 1000.0, 10_000.0) for e in E],
                              MA.mageres_seq_estimates(E, 10_000.0), 10_000.0, MA.MAGEResOptions())
    @test MA.mageres_zircon_sequence_stats_figure(seq2; feeding = (t, n)) !== nothing
end

end
