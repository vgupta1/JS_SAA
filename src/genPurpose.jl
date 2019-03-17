## General purpose functions 

#User is expected to pass functions
#x_k(p0, alpha, mhatk)  : provides a solution with data mhatk
#c_ik(x)      		
#These are often collected into arrays of functions
# xs			: [x_k for k = 1:K ]
# c_k  			: [c_ik for i = 1:d ]
# cs			: (d, K) Matrix with elements c_ik

#General purpose functions
###
function shrink(phat_k, p0, alpha, Nhat_k)
	if Nhat_k <= 0
		return p0
	end
	@. alpha * p0 / (Nhat_k + alpha) + Nhat_k * phat_k / (Nhat_k + alpha)
end

#need to be careful because some problem smay have zero data
function get_GM_anchor(mhats)
	Nhats = vec(sum(mhats, 1))
	non_zero_indx = Nhats .> 0
	vec(mean(mhats[:, non_zero_indx] ./ Nhats[non_zero_indx]', 2))
end

function sim_path(p_k, N::Int)
	@assert N > 0 "Something weird"
    mhats_k = zeros(length(p_k))
    for j in rand(Categorical(p_k), N)
        mhats_k[j] += 1
    end
    mhats_k
end

#ps is a matrix of pk
function sim_path(ps, Nhats)
    K = size(ps, 2)
    out = zeros(size(ps))
    for k = 1:K
        out[:, k] = sim_path(ps[:, k], Nhats[k])
    end
    out    
end

#empirical bayes moment-matched estimates. 
#returns p0, alpha0
#VG this method does not seem to consistenty yield alpha in [0, 1]
function eb_mm_estimates(mhats)
	Nhats = sum(mhats, 1)
	p0 = mean(mhats ./ Nhats, 2)
	K = size(mhats, 2)

	C0 = sum(mhats.^2)/K
	sq_norm_p0 = norm(p0)^2
	Nbar = mean(Nhats)
	var_n = mean( @.(Nhats^2 - Nhats) )

	alpha0 = C0 - sq_norm_p0 * var_n - Nbar
	@assert alpha0 > 0 "Moment matching fails"
	alpha0 /= (1 - sq_norm_p0) / var_n
	alpha0 = 1 / alpha0 - 1

	return p0, alpha0
end

#solves an approximation to the Apriori MSE equation for alpha
#for large K this is like AlphaOR for MSE
#for now, super lazy and just search a grid.
#returns AlphaEstimate, minimizingAlphaIndex 
function mse_estimates(mhats, supps, p0, alpha_grid)
	K = size(mhats, 2)
	Nhats = sum(mhats, 1)
	phats = mhats ./ Nhats
	mu0s = vec(p0' * supps)

	#summary statistics
	muhats = Vector{Float64}(K)
	sigmas = Vector{Float64}(K)
	mse0 = Vector{Float64}(K)
	for k = 1:K
		muhats[k] = dot(phats[:, k], supps[:, k])
		sigmas[k] = dot(phats[:, k], supps[:, k].^2) - muhats[k]^2
		mse0[k] = dot(phats[:, k], @.((supps[:, k] - mu0s[k])^2) )
	end

	#worker function
	function mse(a)
		out = 0.
		for k = 1:K
			out += (Nhats[k] - a^2) / (Nhats[k] + a)^2 * sigmas[k]
			out += ( a / (a + Nhats[k]) )^2 * mse0[k]
		end
		out/K
	end
	out = map(mse, alpha_grid)
	alpha_grid[argmin(out)], argmin(out)
end

function z_k(x_k, c_k, p0, alpha, mhat_k, ps_k, lam_k, lamavg)
	x = x_k(p0, alpha, mhat_k)
	cs = map(c_ik -> c_ik(x), c_k)
	lam_k/lamavg * dot(ps_k, cs) 
end

#actually returns zLOO_k * N * lamavg
function zLOO_k_unsc(x_k, c_k, p0, alpha, mhat_k)
	out = 0.
	mhatloo = mhat_k[:] #copy data
	#only compute for terms where mhat_k[i] > 0
	#compute a base solution to resuse as necessary
	x_base = x_k(p0, alpha, mhat_k)
	for i = 1:length(mhat_k)
		#correct previous toggle
		if i > 1 && mhat_k[i - 1] > 0
			mhatloo[i - 1] += 1
		end

		if mhat_k[i] > 0
			mhatloo[i] -= 1	
		    x = x_k(p0, alpha, mhatloo)
		else 
			x = x_base
		end
		out += mhat_k[i] * c_k[i](x)
	end
	out
end

#average true performance
function zbar(xs, cs, p0, alpha, mhats, ps, lams)
	lamavg = mean(lams)
	K = size(cs, 2)
	out = 0.
	#VG maybe change this to a more numerically stable way to do avg?
	for k = 1:K
		out += z_k(xs[k], cs[:, k], p0, alpha, mhats[:, k], ps[:, k], lams[k], lamavg)
	end
	out/K
end

#full-info for scaling/comparison
#uses the fact that data scale doesn't matter.
zstar(xs, cs, ps, lams) = zbar(xs, cs, ps[:, 1], 0., ps, ps, lams)

function zLOObar_unsc(xs, cs, p0, alpha, mhats)
	K = size(cs, 2)
	out = 0.
	for k = 1:K
		out += zLOO_k_unsc(xs[k], cs[:, k], p0, alpha, mhats[:, k])
	end
	out /K
end

#returns alphaOR, minimizingAlphaIndex, curveInAlpha
function oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
	alphaOR = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zbar(xs, cs, p0, alpha, mhats, ps, lams)
		if out[j] < best_val
			jstar = j
			alphaOR = alpha
			best_val = out[j]
		end
	end
	return alphaOR, jstar, out
end	

#returns alphaLOO, jstar, and cuveInAlpha_unsc
function loo_alpha(xs, cs, mhats, p0, alpha_grid)
	alphaLOO = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zLOObar_unsc(xs, cs, p0, alpha, mhats)
		if out[j] < best_val
			jstar = j
			alphaLOO = alpha
			best_val = out[j]
		end
	end
	return alphaLOO, jstar, out
end	

###  
#### Deprecated
#Convenience functions to generate a sequence of newsvendor problems
#  with common support, potentially different service levels    
#  assumes supp is ordered vector of Reals
# returns cs, xs
# function genNewsvendors(supp, ss, K)
# 	@assert issorted(supp) "Supp must be an ordered vector"

# 	#Generic computation of the sth quantile 
# 	function x_k(p0, alpha, mhat_k, s) 
# 	    const Nhat_k = sum(mhat_k)
# 	    palpha = JS.shrink(mhat_k./Nhat_k, p0, alpha, Nhat_k)
# 	    indx = quantile(Categorical(palpha), s)
# 	    supp[indx]
# 	end

# 	function c_ik(i, x, s)
# 	    supp[i] > x ? s/(1-s) * (supp[i] - x) : (x - supp[i])
# 	end

# 	xs  = [(p0, alpha, mhat_k)-> x_k(p0, alpha, mhat_k, ss[k]) for k = 1:K]
# 	cs = [x->c_ik(i, x, ss[k]) for i = 1:length(supp), k = 1:K]
# 	cs, xs
# end

###  
# Ideally should have a broadcast implementation that combines
# previous two...
function genNewsvendorsDiffSupp(supps, s, K)
	#Generic computation of the sth quantile 
	function x_k(p0, alpha, mhat_k, k, s) 
	    Nhat_k = sum(mhat_k)
	    palpha = JS.shrink(mhat_k./Nhat_k, p0, alpha, Nhat_k)
	    if !isapprox(sum(palpha), 1) || sum(mhat_k .< 0) > 0 
	    	println("Nhat_k:\t", Nhat_k)
	    	println("alpha:\t", alpha)
	    	println("mhat_k:\n", mhat_k)
	    	println("phat_k:\n", palpha)
	    end

	    indx = quantile(Categorical(palpha), s)
	    supps[indx, k]
	end

	function c_ik(i, k, x, s)
	    supps[i, k] > x ? s/(1-s) * (supps[i, k] - x) : (x - supps[i, k])
	end

	xs  = [(p0, alpha, mhat_k)-> x_k(p0, alpha, mhat_k, k, s) for k = 1:K]
	cs = [x->c_ik(i, k, x, s) for i = 1:size(supps, 1), k = 1:K]
	cs, xs
end

#### Helper function for binning data
#dat_vec is vector of actual realizations
#computes a grid of size d
#returns 
	#d vector of lefthand binpoints 
	#d vector of counts
	#raw histogram object for safety
#assumes NA's already filtered from calcualtion
function bin_data(dat_vec, d; TOL = .001)
	#Explicitly make bins to prevent histogram from rounding
	l = minimum(dat_vec)
	u = maximum(dat_vec)
	bin_width = (u - l)/d
    bin_edges = collect(range(l, stop=u, length=d + 1))
    bin_edges[end] += TOL * bin_width  #inflate the last one a little...
    bin_edges[1] -= TOL * bin_width
    dat_hist = fit(Histogram, dat_vec, bin_edges, closed =:left)
    bin_edges[1:d], dat_hist.weights, dat_hist    
end
