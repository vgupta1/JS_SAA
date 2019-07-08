###
# Pred/Pres Synthetic Experiment
##
#### True distributions drawn as dirichlet
#Anchors depend on features X (features drawn gaussian)
#By varying alpha0 in generation, we can vary importance of features

##VG Improved Model would add some noise (unobserved) to the generation of Anchors
##By varying amount of unobserved noise you can change how informative the features are.  

using Distributions, PyPlot, Random, LinearAlgebra, DelimitedFiles
include("../JS_SAA_main.jl")
include("PP_helpers.jl")

const seed = 8675309 
const K = 1000  #VG up this to 10000
const d = 10
const alpha0 = 10
const N = 10
const s = .85
const numRuns = 20  #VG up this 100 or 200


##Some algorithmic tuning stuff
alpha_grid = range(0., stop=100, length=100);  #VG up the range too
h_grid = range(.1, stop=5, length=50)

Random.seed!(seed)

## Generating True distributions
Xs = randn(d, K) #Features are Gaussian

p0s = zeros(d, K)
for ix = 1:K
    p0s[:, ix] .= exp.(Xs[:, ix]) / sum( exp.(Xs[:, ix]) )  #deterministic MNL spec.
end

ps = zeros(d, K)
for ix = 1:K
    ps[:, ix] = rand(Dirichlet(alpha0 * p0s[:, ix]))
end

#Gen Subproblems
cs, xs = JS.genNewsvendorsDiffSupp(repeat(1:d, 1, K), s, K);


#setup output
f = open("PredPresResults_$(alpha0)_$(numRuns)_$(seed).csv", "w")
writedlm(f, ["K" "d" "N" "Method" "TruePerf" "alpha"])

for iRun = 1:numRuns
	Nhats = rand(Poisson(N), K)
	Nhats[Nhats .== 0] .= 1  #cludge to deal with the zero values
	mhats = JS.sim_path(ps, Nhats)
	phats = mhats ./ Nhats'
	phat_avg = JS.get_GM_anchor(mhats)

	#Full-Info Performance for reference
	if iRun == 1
		full_info = JS.zstar(xs, cs, ps, ones(K))
		writedlm(f, [K d N "FullInfo" full_info 0.])
	end

	#compute an SAA solution and true perf
	perf_SAA = JS.zbar(xs, cs, phat_avg, 0., mhats, ps, ones(K))
	writedlm(f, [K d N "SAA" perf_SAA 0.])

	#compute a Shrunken SAA and true perf (GM anchor)
	alphaOR_GM, min_indx, or_alpha_curve_GM = JS.oracle_alpha(xs, cs, mhats, ps, ones(K), phat_avg, alpha_grid)
	perf_SSAA_GM = or_alpha_curve_GM[min_indx]
	writedlm(f, [K d N "S-SAA" perf_SSAA_GM alphaOR_GM])

	#First fit the weights in an oracle fashion.
	hstar, errors, anchors = tune_h_mse(phats, ps, Xs, h_grid)
	println("hstar:\t", hstar)

	#compute predictive/prescriptive sol
	perf_PP = zbar_vec(xs, cs, anchors, 1e10, mhats, ps)
	writedlm(f, [K d N "PP" perf_PP hstar])

	#compute a Shrunken SAA with pred/pres anchors
	alphaOR_PP, min_indx, or_alpha_curve_PP = oracle_alpha_vec(xs, cs, mhats, ps, anchors, alpha_grid)
	perf_SSAA_PP = or_alpha_curve_PP[min_indx]
	writedlm(f, [K d N "S-PP" perf_SSAA_PP alphaOR_PP])

end
close(f)
