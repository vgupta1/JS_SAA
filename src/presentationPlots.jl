###
# Makes the data for series of plots used in the presentation
###
using Distributions, Random, DelimitedFiles
include("../src/JS_SAA_main.jl")
  
#Then update the bad example to the one from the paper.  
function genSinglePathCurves(seed=8675309)
	K = 10000
	s = .91
	d = 75
	N = 10
	f = open("../Results/single_path_curves_presentation.csv", "w")

	Random.seed!(seed)

	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	ps = [ps rand(Dirichlet(.025 * ones(d)), K - floor(Int, K/2)) ]
	@assert size(ps, 2) == K
	@assert size(ps, 1) == d

	#gen problems
	supp = collect(1:d)
	cs = JS.getNewsVendorCosts(repeat(supp, inner=(1, K)), s, K)
	xs = JS.genSSAAtrainers(repeat(supp, inner=(1, K)), s, K)
	Nhats = rand(Poisson(N), K)
	Nhats[ Nhats .== 0] .= 1
	lams = ones(K)

	full_info = JS.zstar(xs, cs, ps, lams)

	#gen data
	mhats = JS.sim_path(ps, Nhats)
	p0 = JS.get_GM_anchor(mhats)
	alpha_grid = range(0, stop=3N, length=100)

	#note, this actually computes N * ZLoo
	alphaOR, jstar, outOR =  JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
	alphaLOO, jloo, loo_perf = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
	loo_perf /= N

	saa_subopt = map(a-> JS.zbar(xs, cs, p0, a, mhats, mhats, lams), alpha_grid) ./ N
	saa_zero = saa_subopt[1]
	saa_subopt = saa_subopt .- saa_zero
	instab = loo_perf - saa_subopt .- saa_zero;

	#Header
	writedlm(f, ["alpha" "OR" "LOO" "SAASubopt" "Instab"], ',')
	writedlm(f, [vec(alpha_grid) outOR loo_perf saa_subopt instab], ',')
	close(f)
end


#K = 10,0000 Same set up as above
function genHistogramPlots_multiple(numRuns = 1000)
	Random.seed!(8675309)
	outPath = "../Results/UPDATED_singleKUnifNewsvendor"
	usePoisson=true

	#Gen the data
	K = 10000
	s = .91
	d = 75
	N = 10

	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	ps = [ps rand(Dirichlet(.025 * ones(d)), K - floor(Int, K/2)) ]
	@assert size(ps, 2) == K
	@assert size(ps, 1) == d
	p0 = ones(d) ./ d

	alpha_grid = range(0, stop=3N, length=100)
	lams = ones(K)

	cs = JS.getNewsVendorCosts(repeat(collect(1:d), inner=(1, K)), s, K)
	xs = JS.genSSAAtrainers(repeat(collect(1:d), inner=(1, K)), s, K)

	f = open("$(outPath)_$(K)_$(numRuns).csv", "w")

	#Header
	writedlm(f, ["Run" "Method" "TruePerf" "time" "alpha"], ',')

	#Compute the full-info value once for reference
	t =
	  @elapsed full_info = JS.zstar(xs, cs, ps, lams)
	println(full_info)
	writedlm(f, [1 "FullInfo" full_info t 0.], ',')

	Nhats = fill(N, K)
	for iRun = 1:numRuns
		println("iRun \t:", iRun)
		if usePoisson
			Nhats = rand(Poisson(N), K)
			Nhats[ Nhats.== 0 ] .= 1
		end
		mhats = JS.sim_path(ps, Nhats);

		#for data-driven shrinkage anchor
		phat_avg = JS.get_GM_anchor(mhats)

		#Gen the SAA const
		# t = 
		#   @elapsed perf_SAA = JS.zbar(xs, cs, p0, 0., mhats, ps, lams)
		# writedlm(f, [iRun "SAA" perf_SAA t 0.0], ',')

		#Gen the Oracle const
		t = 
		  @elapsed alphaOR, min_indx, or_alpha_curve = JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
		writedlm(f, [iRun "Oracle" or_alpha_curve[min_indx] t alphaOR], ',')

		#Gen the LOO cost with 1/d shrinkage
		# t = 
		#   @elapsed alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
		# writedlm(f, [iRun "LOO_unif" or_alpha_curve[min_indx] t alphaLOO], ',')

		#Gen the LOO cost with the phatAvg shrinkage
		t = 
		  @elapsed alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, phat_avg, alpha_grid)
		writedlm(f, [iRun "LOO_avg" or_alpha_curve[min_indx] t alphaLOO], ',')

		#The SAA version
		writedlm(f, [iRun "SAA" or_alpha_curve[1] -1 0.0], ',')

		flush(f)
	end #endRun
	close(f)	
end

function genBadSinglePathCurves(seed=8675409)
	K = 10000
	N = 10
	s = .5
	p0 = .3

	#Gen Pks uniformly over interval for now
	ps = rand(K) 
	ps = .6 .+ .3 * rand(K)

	f = open("../Results/single_path_curves_bad_presentation.csv", "w")

	Nhats = rand(Poisson(N), K)
	Nhats[ Nhats .== 0] .= 1
	mhats = JS.nv_sim_path(ps, Nhats);

	Random.seed!(seed)
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
	writedlm(f, ["alpha" "OR" "LOO" "SAASubopt" "Instab"], ',')
	writedlm(f, [vec(alpha_grid) outOR loo_perf saa_subopt instab], ',')
	close(f)
end

##
#genHistogramPlots_multiple()
genBadSinglePathCurves()
#genSinglePathCurves()
