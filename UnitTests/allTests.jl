##  a bare bones Test-Suite for some critical functions
using Test, Distributions, Random, LinearAlgebra
#using DelimitedFiles

include("../src/JS_SAA_main.jl")
include("consts.jl")

@testset "allTests" begin 

#Tests shrink/crossval other simpler helpers
@testset "Vector Functions Bernolli Newsvendor" begin
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
	cs, xs = JS.genNewsvendorsDiffSupp(supps, s, K)

	#First some pointwise tests
	@test isapprox(JS.z_k(xs[1], cs[:, 1], p0, 1, mhats[:, 1], ps[:, 1], 1, 1), 3.2530462468808903)
	@test isapprox(JS.zbar(xs, cs, p0, 1, mhats, ps, lams), 4.428156384118435)
	@test isapprox(JS.zstar(xs, cs, ps, lams), 4.361044771184344)

	@test isapprox(JS.zLOO_k_unsc(xs[1], cs[:, 1], p0, 1, mhats[:, 1]), 49.0 )

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

end #shrunkenSAA

end #AllTEsts

