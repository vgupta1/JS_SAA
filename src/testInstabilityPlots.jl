### Instability Plots
# Generates the code for the intability/suboptimality in the paper
# Presentation uses similar plots but with different paramters, see presentation plots.jl

using Distributions, Random, DelimitedFiles
include("../src/JS_SAA_main.jl")


function genSinglePathCurves(subcase; seed=8675409, usePoisson=true)
	K = 10000
	N = 10

	if subcase == :TradeoffBad
		p0 = .3
		s = .5
		path = string("../Results/TradeoffBad", usePoisson, ".csv")
	elseif subcase == :TradeoffGoodP0
		p0 = .75
		s = .5
		path = string("../Results/TradeoffGoodP0", usePoisson, ".csv")
	elseif subcase == :TradeoffGoodS
		p0 = .3
		s = .2
		path = string("../Results/TradeoffGoodS", usePoisson, ".csv")
	end
	Random.seed!(seed)

	#Gen Pks uniformly over interval for now
	ps = rand(K) 
	ps = .6 .+ .3 * rand(K)

	f = open(path, "w")

	if usePoisson 
		Nhats = rand(Poisson(N), K)
		Nhats[Nhats .== 0] .= 1
	else
		Nhats = fill(N, K)
	end
	mhats = JS.nv_sim_path(ps, Nhats);

	alpha_grid = range(0, stop=2N, length=100)

	#note, this actually computes N * ZLoo
	#loo_perf = map(a->JS.nv_loo(mhats, p0, a, Nhats, s), alpha_grid)
	alphaOR, jstar, outOR = JS.nv_oracle_alpha(mhats, ps, p0, alpha_grid, Nhats, s)
	alphaLOO, jstar, loo_perf = JS.nv_loo_alpha(mhats, p0, alpha_grid, Nhats, s)
	loo_perf /= N

	saa_subopt = map(a-> JS.nv_saa(mhats, p0, a, Nhats, s), alpha_grid) ./ N
	saa_zero = saa_subopt[1]
	saa_subopt = saa_subopt .- saa_zero
	instab = loo_perf - saa_subopt .- saa_zero;

	#Header
	writedlm(f, ["Alpha" "OR" "LOO" "SAA-SubOpt" "Instability"], ',')
	writedlm(f, [vec(alpha_grid) outOR loo_perf saa_subopt instab], ',')
	close(f)
end

genSinglePathCurves(:TradeoffBad)
genSinglePathCurves(:TradeoffGoodP0)
genSinglePathCurves(:TradeoffGoodS)

genSinglePathCurves(:TradeoffBad, usePoisson=false)
genSinglePathCurves(:TradeoffGoodP0, usePoisson=false)
genSinglePathCurves(:TradeoffGoodS, usePoisson=false)
