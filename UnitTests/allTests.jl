##  a bare bones Test-Suite for some critical functions
using Base.Test
using Distributions

include("../src/JS_SAA_main.jl")
include("consts.jl")

@testset "All Tests" begin 

#Tests shrink/crossval other simpler helpers
@testset "Vector Functions Baby Newsvendor" begin
	srand(8675309)
	K = 100
	ps = rand(K);
	alpha = .1
	s = .9
	p0 = .5
	Nhat = 20

	@test isapprox(JS.exp_nv_obj(ps, p0, alpha, Nhat, s), 
					0.4670605847641166) 
	@test isapprox(JS.exp_nv_fullInfo(ps, s), 
					0.4472649361852286)
end

@testset "bestAlpha Functions" begin
	srand(8675309)
	K = 100
	ps = rand(K);
	Nhats = rand(Poisson(20), K)
	alpha_grid = linspace(0, 80, 80)
	ss = rand(K)/4 + .75
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

	p0_grid = linspace(0, 1, 20)
	p0OR, alphaOR, jstar = JS.nv_oracle_both(mhats, ps, p0_grid, alpha_grid, Nhats, ss)
	@test isapprox(p0OR, 0.21052631578947367)
	@test isapprox(alphaOR, 14.177215189873417)
	@test jstar == 15

	p0LOO, alphaLOO, jstar = JS.nv_loo_both(mhats, p0_grid, alpha_grid, Nhats, ss)
	@test isapprox(p0LOO, 0.15789473684210525)
	@test isapprox(alphaLOO, 11.139240506329115)
	@test jstar == 12

end

@testset "genNewsvendorFunctions" begin
 	srand(8675309)
	const K = 20
	supp = collect(1:10)
	p0 = rand(10)
	p0 /= sum(p0)  
	alpha = .12
	pk = ones(10)/10  #unif distribution
	mhat = JS.sim_path(pk, 8)

	cs, xs = JS.genNewsvendors(supp, rand(K), K)
	sols = [x(p0, alpha, mhat) for x in xs]
	for ix = 1:length(sols)
		@test isapprox(sols[ix], cSols[ix])
	end

	#evaluate all the costs at x = 5
	costs = zeros(K)
	for k = 1:K
		costs[k] = dot(mhat/8, [cs[ix, k](5) for ix = 1:10])
		@test isapprox(costs[k], cCosts[k])
	end

end #end genNewsNvendor Functions

end #all tests

