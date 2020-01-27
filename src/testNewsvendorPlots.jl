###
# Generates a histogram for SAA, Alpha LOO, and Alpha OR
# Used for illustrative picture in presentation
###
using Distributions, Random, DelimitedFiles
include("../src/JS_SAA_main.jl")

# Generates data that can be used to create a histogram of costs
# ps are half uniform, and half slightly concentrated at 1/d
# all data is synthetic
function runTest(numRuns, K, outPath; usePoisson=true)
	s = .95
	d = 10
	N = 20
	Random.seed!(8675309)
	p0 = ones(d)/d
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	ps = [ps rand(Dirichlet(3 * ones(d)), K - floor(Int, K/2)) ]
	@assert size(ps, 2) == K
	@assert size(ps, 1) == d
	alpha_grid = range(0, stop=50, length=75)
	lams = ones(K)

	cs, xs = JS.genNewsvendorsDiffSupp(repeat(1:d, inner=(1, K)), s, K)

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
		if usePoisson
			Nhats = rand(Poisson(N), K)
		end
		mhats = JS.sim_path(ps, Nhats);

		#for data-driven shrinkage anchor
		phat_avg = vec(mean(mhats ./ Nhats', dims=2))
	
		#Gen the SAA const
		t = 
		  @elapsed perf_SAA = JS.zbar(xs, cs, p0, 0., mhats, ps, lams)
		writedlm(f, [iRun "SAA" perf_SAA t 0.0], ',')

		#Gen the Oracle const
		t = 
		  @elapsed alphaOR, min_indx, or_alpha_curve = JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
		writedlm(f, [iRun "Oracle" or_alpha_curve[min_indx] t alphaOR], ',')

		#Gen the LOO cost with 1/d shrinkage
		t = 
		  @elapsed alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
		writedlm(f, [iRun "LOO_unif" or_alpha_curve[min_indx] t alphaLOO], ',')

		#Gen the LOO cost with the phatAvg shrinkage
		t = 
		  @elapsed alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, phat_avg, alpha_grid)
		writedlm(f, [iRun "LOO_avg" or_alpha_curve[min_indx] t alphaLOO], ',')

		flush(f)
	end #endRun
	close(f)	
end
runTest(10, 10, "../Results/tempFile", usePoisson=false)

# Sept 2019:  This is the run that is used the paper to generate histogram plots for intro.
#Run again with K = 1000 for comparison graph in presentation
#runTest(1000, 1000, "../Results/singleKUnifNewsvendor", usePoisson=false)