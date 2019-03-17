###
# Generates a histogram for SAA, Alpha LOO, and Alpha OR
# Used for illustrative picture in presentation
###
using Distributions
include("../src/JS_SAA_main.jl")

# Generates data that can be used to create a histogram of costs
# ps are half uniform, and half slightly concentrated at 1/d
# all data is synthetic
function runTest(numRuns, K, outPath)
	const s = .95
	const d = 10
	const N = 20
	srand(8675309)
	p0 = ones(d)/d
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	ps = [ps rand(Dirichlet(3 * ones(d)), K - floor(Int, K/2)) ]
	@assert size(ps, 2) == K
	@assert size(ps, 1) == d
	alpha_grid = linspace(0, 50, 75)
	lams = ones(K)

	cs, xs = JS.genNewsvendors(collect(1:d), s * ones(K), K)
	f = open("$(outPath)_$(K)_$(numRuns).csv", "w")

	#Header
	writecsv(f, ["Run" "Method" "TruePerf" "time" "alpha"])

	#Compute the full-info value once for reference
	tic()
	full_info = JS.zstar(xs, cs, ps, lams)
	t = toc()
	println(full_info)
	writecsv(f, [1 "FullInfo" full_info t 0.])

	for iRun = 1:numRuns
		Nhats = rand(Poisson(N), K)
		mhats = JS.sim_path(ps, Nhats);

		#for data-driven shrinkage anchor
		phat_avg = vec(mean(mhats ./ Nhats', 2))
	
		#Gen the SAA const
		tic()
		perf_SAA = JS.zbar(xs, cs, p0, 0., mhats, ps, lams)
		t = toc()
		writecsv(f, [iRun "SAA" perf_SAA t 0.0])

		#Gen the Oracle const
		tic()
		alphaOR, min_indx, or_alpha_curve = JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
		t = toc()
		writecsv(f, [iRun "Oracle" or_alpha_curve[min_indx] t alphaOR])

		#Gen the LOO cost with 1/d shrinkage
		tic()
		alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
		t = toc()
		writecsv(f, [iRun "LOO_unif" or_alpha_curve[min_indx] t alphaLOO])

		#Gen the LOO cost with the phatAvg shrinkage
		tic()
		alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, phat_avg, alpha_grid)
		t = toc()
		writecsv(f, [iRun "LOO_avg" or_alpha_curve[min_indx] t alphaLOO])

		flush(f)
	end #endRun
	close(f)	
end
runTest(10, 10, "../Results/tempFile")

#Run once with K = 1 for initial graph in presentation
runTest(1000, 1, "../Results/singleKUnifNewsvendor")

#Run again with K = 1000 for comparison graph in presentation
runTest(1000, 1000, "../Results/singleKUnifNewsvendor")