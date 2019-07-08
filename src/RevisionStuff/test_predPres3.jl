###
# Pred/Pres Synthetic Experiment 3
##
#### True distributions drawn as dirichlet
#Anchors depend on a UNIVARIATE Feature X + a noise term (MNL)
#By varying noise term in generation, we can vary importance of features
#By varying alpha0 in dirichlet, we vary importance of shrinkage.
#Feature is univariate to help the NW methods 

using Distributions, PyPlot, Random, LinearAlgebra, DelimitedFiles
include("../JS_SAA_main.jl")
include("PP_helpers.jl")

const seed = 8675309 
const K = 100  #VG up this to 10000
const d = 10
const alpha0 = 10
const N = 10
const s = .85
const numRuns = 50  #VG up this 100 or 200
const sigma = 3
const X_std = 10.

##Some algorithmic tuning stuff
alpha_grid = range(0., stop=100, length=100);  #VG up the range too
h_grid = range(.01, stop=5, length=100)

Random.seed!(seed)

#####
## Generating True distributions
beta = randn(d, 1)

Xs = X_std * randn(1, K) #Features are Gaussian
noise = sigma * randn(d, K)


p0s = zeros(d, K)
for ix = 1:K
    p0s[:, ix] .= exp.(beta * Xs[:, ix]  .+ noise[:, ix]) / sum( exp.(beta * Xs[:, ix] .+ noise[:, ix]) )
end

ps = zeros(d, K)
for ix = 1:K
    ps[:, ix] = rand(Dirichlet(alpha0 * p0s[:, ix]))
end
########

#Gen Subproblems
cs, xs = JS.genNewsvendorsDiffSupp(repeat(1:d, 1, K), s, K);

#setup output
f = open("temp3Results_$(alpha0)_$(sigma)_$(X_std)_$(numRuns)_$(seed).tsv", "w")
writedlm(f, ["K" "d" "N" "Method" "TruePerf" "alpha" "bandwidth"])

for iRun = 1:numRuns
	Nhats = rand(Poisson(N), K)
	Nhats[Nhats .== 0] .= 1  #cludge to deal with the zero values
	mhats = JS.sim_path(ps, Nhats)
	phats = mhats ./ Nhats'
	phat_avg = JS.get_GM_anchor(mhats)

	#Full-Info Performance for reference
	if iRun == 1
		full_info = JS.zstar(xs, cs, ps, ones(K))
		writedlm(f, [K d N "FullInfo" full_info 0. -1.])
	end

	#compute an SAA solution and true perf
	perf_SAA = JS.zbar(xs, cs, phat_avg, 0., mhats, ps, ones(K))
	writedlm(f, [K d N "SAA" perf_SAA 0. -1.])

	#compute a Shrunken SAA and true perf (GM anchor)
	#VG For now use oracle to just get a sense
	alphaOR_GM, min_indx, or_alpha_curve_GM = JS.oracle_alpha(xs, cs, mhats, ps, ones(K), phat_avg, alpha_grid)
	perf_SSAA_GM = or_alpha_curve_GM[min_indx]
	writedlm(f, [K d N "S-SAA" perf_SSAA_GM alphaOR_GM -1.])

	##NW values 
	#First fit the NW weights by minimizing MSE
	hstar, errors, anchors = tune_h_mse(phats, ps, Xs, h_grid)
	println("h-NW:\t", hstar)

	#compute the NW solution 
	perf_PP = zbar_vec(xs, cs, anchors, 1e10, mhats, ps)
	writedlm(f, [K d N "NW" perf_PP 0. hstar])

	#compute a Shrunken SAA with NW anchors just for fun.
	alphaOR_PP, min_indx, or_alpha_curve_PP = oracle_alpha_vec(xs, cs, mhats, ps, anchors, alpha_grid)
	perf_SSAA_PP = or_alpha_curve_PP[min_indx]
	writedlm(f, [K d N "HS-NW" perf_SSAA_PP alphaOR_PP hstar])

	##NW featbar values
	#Fit new NW weights by minimizing MSE
	hstar, errors, anchors = tune_h_mse_bar(phats, ps, Xs, h_grid)
	println("h-NW-F:\t", hstar)

	#Compute the NW with featbar full vector
	perf_PP = zbar_vec(xs, cs, anchors, 1e10, mhats, ps)
	writedlm(f, [K d N "NWbar" perf_PP 0. hstar])

	#Compute the S-NW solution
    alphaOR, jstar, kstar, or_alpha_h_curve = oracle_alpha_h(xs, cs, mhats, Xs, ps, get_anchors, alpha_grid, h_grid)
    perf_SNW = or_alpha_h_curve[jstar, kstar]
    println("h-SNW:\t", h_grid[kstar])
    writedlm(f, [K d N "S-NW" perf_SNW alphaOR h_grid[kstar]])
end
close(f)
