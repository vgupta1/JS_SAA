##  a bare bones Test-Suite for some critical functions
using Base.Test
include("../src/babyNewsVendor.jl")

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
	Nhats = 20 * floor.(Int, rand(K))
	alpha_grid = linspace(0, 80, 80)
	ss = rand(K)/4 + .75
	p0 = .5
	mhats = JS.sim_path(ps, Nhats)
	alphaOR, jstar, out = JS.oracle_alpha(mhats, ps, p0, alpha_grid, Nhats, ss)
	@test isapprox(out[jstar], 7.72698286347512)
	@test isapprox(alphaOR, 0)
	
	alphaAP, jstar, out = JS.apriori_alpha(ps, p0, alpha_grid, Nhats, ss)
	@test isapprox(out[jstar], 0.5096223445262258)
	@test isapprox(alphaAP, 1.0126582278481013)

end

end #all tests

