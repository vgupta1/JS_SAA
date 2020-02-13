##  a bare bones Test-Suite for some critical functions
using Test, Distributions, Random, LinearAlgebra
using DelimitedFiles

include("../src/JS_SAA_main.jl")
include("consts.jl")

@testset "allTests" begin 

#Tests shrink/crossval other simpler helpers
@testset "Vector Functions Bernoulli Newsvendor" begin
	Random.seed!(8675309)
	K = 100
	ps = rand(K)
	alpha = .1
	s = .9
	p0 = .5
	Nhat = 20

	@test isapprox(JS.exp_nv_obj(ps, p0, alpha, Nhat, s), 
					0.4670605847641166) 
	@test isapprox(JS.exp_nv_fullInfo(ps, s), 
					0.4472649361852286)
end

@testset "bestAlpha Bernoulli Newsvendor" begin
	Random.seed!(8675309)
	K = 100
	ps = rand(K)
	Nhats = rand(Poisson(20), K)
	alpha_grid = range(0, stop=80, length=80)
	ss = rand(K)/4 .+ .75
	p0 = .5
	mhats = JS.nv_sim_path(ps, Nhats)
	alphaOR, jstar, out = JS.nv_oracle_alpha(mhats, ps, p0, alpha_grid, Nhats, ss)
	@test isapprox(out[jstar], 0.4323864805448839)
	@test isapprox(alphaOR, 2.0253164556962027)

	alphaAP, jstar, out = JS.nv_apriori_alpha(ps, p0, alpha_grid, Nhats, ss)
	@test isapprox(out[jstar], 0.4338668678652777)
	@test isapprox(alphaAP, 2.025316455696202)

	alphaLOO, jstar, out = JS.nv_loo_alpha(mhats, p0, alpha_grid, Nhats, ss)
	@test isapprox(out[jstar], 7.926847617453273)
	@test isapprox(alphaLOO, 1.0126582278481013)

	p0_grid = range(0, stop=1, length=20)
	p0OR, alphaOR, jstar = JS.nv_oracle_both(mhats, ps, p0_grid, alpha_grid, Nhats, ss)
	@test isapprox(p0OR, 0.21052631578947367)
	@test isapprox(alphaOR, 14.177215189873417)
	@test jstar == 15

	p0LOO, alphaLOO, jstar = JS.nv_loo_both(mhats, p0_grid, alpha_grid, Nhats, ss)
	@test isapprox(p0LOO, 0.15789473684210525)
	@test isapprox(alphaLOO, 11.139240506329115)
	@test jstar == 12

end

@testset "simpleHelpers" begin
	K = 3
	d = 10
	N = 10
	Random.seed!(8675309)

	p0 = ones(d)/d
	ps = rand(Dirichlet(ones(d)), K)
	shrunk_p = JS.shrink(ps[:, 1], p0, .5, N)
		
	for i = 1:d
		@test isapprox(cShrunkP1[i], shrunk_p[i])
	end

	#generate proper data and calc mhat
	Nhats = rand(Poisson(N), K)
	mhats = JS.sim_path(ps, Nhats)
	
	pGM = JS.get_GM_anchor(mhats)
	for i = 1:d
		@test isapprox(cGenAnchor[i], pGM[i])
	end

	Nhats[1] = 0
	mhats = JS.sim_path(ps, Nhats)	
	pGM = JS.get_GM_anchor(mhats)
	for i = 1:d
		@test isapprox(cGenAnchorZero[i], pGM[i])
	end

end

@testset "shrunkenSAA" begin
	#discrete newsvendors supported 1:d
	K = 10
	s = .95
	d = 10
	N = 20
	Random.seed!(8675309)

	#gen an "interesting" distribution of ps still centered at 1/d
	p0 = ones(d)/d
	anchor = vcat(1., zeros(d-1))
	p0 = .5 * p0 + .5 * anchor
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	qs = rand(Dirichlet(5 * ones(d)), K-floor(Int, K/2))
	ps = [ps qs]


	alpha_grid = range(0, stop=50, length=75)
	lams = ones(K)

	Nhats = rand(Poisson(N), K)
	mhats = JS.sim_path(ps, Nhats);

	supps =  repeat(collect(1:d), outer=(1, K))
	cs = JS.getNewsVendorCosts(supps, s, K)
	xs = JS.genSSAAtrainers(supps, s)

	#First some pointwise tests
	@test isapprox(JS.z_k(xs[1], cs[:, 1], mhats[:, 1], ps[:, 1], 1, 1, (p0, 1)), 3.2530462468808903)
	@test isapprox(JS.zbar(xs, cs, mhats, ps, lams, (p0, 1)), 4.428156384118435)
	@test isapprox(JS.zstar(xs, cs, ps, lams), 4.361044771184344)

	@test isapprox(JS.zLOO_k_unsc(xs[1], cs[:, 1], mhats[:, 1], (p0, 1)), 49.0 )

	#Now check whole curve
	alphaOR, jstar, oracle_curve = JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
	@test isapprox(2.027027027027027, alphaOR)
	@test isapprox(jstar, 4)

	for i = 1:length(cOracleCurve)
		@test isapprox(cOracleCurve[i], oracle_curve[i])
	end

	alphaLOO, jstar, loo_curve = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
	@test isapprox(alphaLOO, 0.)
	@test isapprox(jstar, 1)
	for i = 1:length(cLooCurve)
		@test isapprox(cLooCurve[i], loo_curve[i])
	end

	#Check the entire cross-val curve for K-fold
	alphaCV, jstar, cv_curve = JS.cv_alpha(xs, cs, mhats, p0, alpha_grid, 5)
	@test isapprox(alphaCV, 8.783783783783784)
	@test isapprox(jstar, 14)

	for i = 1:length(cCVCurve)
		@test isapprox(cCVCurve[i], cv_curve[i])
	end

end #shrunkenSAA



@testset "KS Robust tests" begin
	#discrete newsvendors supported 1:d
	K = 10
	s = .95
	d = 10
	N = 20
	Random.seed!(8675309)

	#gen an "interesting" distribution of ps still centered at 1/d
	p0 = ones(d)/d
	anchor = vcat(1., zeros(d-1))
	p0 = .5 * p0 + .5 * anchor
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	qs = rand(Dirichlet(5 * ones(d)), K-floor(Int, K/2))
	ps = [ps qs]

	alpha_grid = range(0, stop=50, length=75)
	lams = ones(K)
	Nhats = rand(Poisson(N), K)
	mhats = JS.sim_path(ps, Nhats);

	supps =  repeat(collect(1:d), outer=(1, K))

	#static sets
	cs = JS.getNewsVendorCosts(supps, s, K)
	xsKS = JS.genKSTrainers(supps, s, K, :aPriori)
	xsKS2 = JS.genKSTrainers(supps, s, K, :crossVal)

	#First some pointwise tests for different gamma
	@test isapprox(JS.zbar(xsKS, cs, mhats, ps, lams, 0),    4.427097278114614)
	@test isapprox(JS.zbar(xsKS, cs, mhats, ps, lams, .025), 4.38864696562911)
	@test isapprox(JS.zbar(xsKS, cs, mhats, ps, lams, .05),   4.6573883655136825)

	#one test for the cross-val version
	Gamma_grid = range(0, stop=1-s, length=51)
	@test isapprox(JS.zbar(xsKS2, cs, mhats, ps, lams, (Gamma_grid, 5)), 4.559215177976772)
	
end #KS Tests

@testset "optimizingAnchors" begin
	#discrete newsvendors supported 1:d
	K = 10
	s = .95
	d = 10
	N = 20
	Random.seed!(8675309)

	#gen an "interesting" distribution of ps still centered at 1/d
	p0 = ones(d)/d
	anchor = vcat(1., zeros(d-1))
	p0 = .5 * p0 + .5 * anchor
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	qs = rand(Dirichlet(5 * ones(d)), K-floor(Int, K/2))
	ps = [ps qs]

	alpha_grid = range(0, stop=50, length=75)
	lams = ones(K)

	Nhats = rand(Poisson(N), K)
	mhats = JS.sim_path(ps, Nhats);

	supps =  repeat(collect(1:d), outer=(1, K))
	cs = JS.getNewsVendorCosts(supps, s, K)
	xs = JS.genSSAAtrainers(supps, s)

	Random.seed!(123456789)
	optp0, optAlpha, zstar = JS.loo_anchor(xs, cs, mhats, numClusters=5)
	@test isapprox(optAlpha, 17.100025815457034)
	@test isapprox(zstar, 87.0)
	for i = 1:length(cOptAnchor1)
		@test isapprox(cOptAnchor1[i], optp0[i])
	end

	optp0, optAlpha, zstar = JS.loo_anchor(xs, cs, mhats, numClusters=-1)
	@test isapprox(optAlpha, 4.174826098336609)
	@test isapprox(zstar, 86.6)
	for i = 1:length(cOptAnchor2)
		@test isapprox(cOptAnchor2[i], optp0[i])
	end

	optp0, optAlpha, zstar = JS.opt_oracle_anchor(xs, cs, ps, mhats; numClusters = 5)
	@test isapprox(optAlpha, 4.466926141845736)
	@test isapprox(zstar, 4.382226891260787)
	for i = 1:length(cOptAnchor3)
		@test isapprox(cOptAnchor3[i], optp0[i])
	end

	optp0, optAlpha, zstar = JS.opt_oracle_anchor(xs, cs, ps, mhats; numClusters = -1)
	@test isapprox(optAlpha, 3.542168883777246)
	@test isapprox(zstar, 4.361044771184344)
	for i = 1:length(cOptAnchor4)
		@test isapprox(cOptAnchor4[i], optp0[i])
	end


end #end optimizingAnchors

@testset "betaAnchors" begin
	#discrete newsvendors supported 1:d
	K = 10
	s = .95
	d = 10
	N = 20
	Random.seed!(8675309)

	#gen an "interesting" distribution of ps still centered at 1/d
	p0 = ones(d)/d
	anchor = vcat(1., zeros(d-1))
	p0 = .5 * p0 + .5 * anchor
	ps = rand(Dirichlet(ones(d)), floor(Int, K/2))
	qs = rand(Dirichlet(5 * ones(d)), K-floor(Int, K/2))
	ps = [ps qs]

	alpha_grid = range(0, stop=50, length=75)
	lams = ones(K)

	Nhats = rand(Poisson(N), K)
	mhats = JS.sim_path(ps, Nhats);

	supps =  repeat(collect(1:d), outer=(1, K))
	cs = JS.getNewsVendorCosts(supps, s, K)
	xs = JS.genSSAAtrainers(supps, s)

	theta2_grid = range(1e-6, stop=3, length=5)
	mu_grid = range(1e-6, stop=1, length=5)
	alphaLOO, p0LOO, zLOO = JS.loo_betaAnchor(xs, cs, mhats, alpha_grid, theta2_grid, mu_grid)
	@test isapprox(alphaLOO, 10.81081081081081)
	@test isapprox(zLOO, 87.0)

	for i = 1:length(p0LOO)
		@test isapprox(p0LOO[i], cBetaAnchor[i])
	end

	alphaOR, p0OR, zOR = JS.oracle_betaAnchor(xs, cs, mhats, ps, lams, alpha_grid, theta2_grid, mu_grid)
	@test isapprox(alphaOR, 1.3513513513513513)
	@test isapprox(zOR, 4.361044771184344)
	for i = 1:length(p0OR)
		@test isapprox(p0OR[i], cBetaAnchorOR[i])
	end
end #end betaAnchors



end #AllTEsts

