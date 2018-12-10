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
shrink(phat_k, p0, alpha, Nhat_k) = @. alpha * p0 / (Nhat_k + alpha) + Nhat_k * phat_k / (Nhat_k + alpha)

function sim_path(p_k, N::Int)
    mhats_k = zeros(length(p_k))
    for j in rand(Categorical(p_k), N)
        mhats_k[j] += 1
    end
    mhats_k
end

#empirical bayes moment-matched estimates. 
#returns p0, alpha0
#VG this method does not seem to consistenty yield alpha in [0, 1]
function eb_mm_estimates(mhats)
	Nhats = sum(mhats, 1)
	p0 = mean(mhats ./ Nhats, 2)
	const K = size(mhats, 2)

	const C0 = sum(mhats.^2)/K
	const sq_norm_p0 = norm(p0)^2
	const Nbar = mean(Nhats)
	const var_n = mean( @.(Nhats^2 - Nhats) )

	println(C0)
	println(var_n)

	alpha0 = C0 - sq_norm_p0 * var_n - Nbar
	@assert alpha0 > 0 "Moment matching fails"
	alpha0 /= (1 - sq_norm_p0) / var_n
	alpha0 = 1 / alpha0 - 1

	return p0, alpha0
end

#solves an approximation to the Apriori MSE equation for alpha
#for large K this is like AlphaOR for MSE
#for now, super lazy and just search a grid.
#returns minimizingAlphaIndex, and AlphaEstimate
function mse_estimates(mhats, supp, p0, alpha_grid)
	const K = size(mhats, 2)
	Nhats = sum(mhats, 1)
	phats = mhats ./ Nhats
	const mu0 = dot(p0, supp)

	#summary statistics
	muhats = Vector{Float64}(K)
	sigmas = Vector{Float64}(K)
	mse0 = Vector{Float64}(K)
	for k = 1:K
		muhats[k] = dot(phats[:, k], supp)
		sigmas[k] = dot(phats[:, k], supp.^2) - muhats[k]^2
		mse0[k] = dot(phats[:, k], @.((supp - mu0)^2) )
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
	indmin(out), alpha_grid[indmin(out)]
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

function z_k(x_k, c_k, p0, alpha, mhat_k, ps_k, lam_k, lamavg)
	x = x_k(p0, alpha, mhat_k)
	cs = map(c_ik -> c_ik(x), c_k)
	lam_k/lamavg * dot(ps_k, cs) 
end

#actually returns zLOO_k * N * lamavg
function zLOO_k_unsc(x_k, c_k, p0, alpha, mhat_k)
	out = 0.
	mhatloo = mhat_k[:]
	#only compute for terms where mhat_k[i] > 0
	for i = 1:length(mhat_k)
		if mhat_k[i] > 0
		   mhatloo[i] -= 1	
		end

		#correct previous toggle
		if i > 1 && mhat_k[i - 1] > 0
			mhatloo[i - 1] += 1
		end

		x = x_k(p0, alpha, mhatloo)
		out += mhat_k[i] * c_k[i](x)
	end
	out
end

#average true performance
function zbar(xs, cs, p0, alpha, mhats, ps, lams)
	const lamavg = mean(lams)
	const K = size(cs, 2)
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
	const K = size(cs, 2)
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
#Convenience functions to generate a sequence of newsvendor problems
#  with common support, potentially different service levels    
#  assumes supp is ordered vector of Reals
# returns cs, xs
function genNewsvendors(supp, ss, K)
	@assert issorted(supp) "Supp must be an ordered vector"

	#Generic computation of the sth quantile 
	function x_k(p0, alpha, mhat_k, s) 
	    const Nhat_k = sum(mhat_k)
	    palpha = JS.shrink(mhat_k./Nhat_k, p0, alpha, Nhat_k)
	    indx = quantile(Categorical(palpha), s)
	    supp[indx]
	end

	function c_ik(i, x, s)
	    supp[i] > x ? s/(1-s) * (supp[i] - x) : (x - supp[i])
	end

	xs  = [(p0, alpha, mhat_k)-> x_k(p0, alpha, mhat_k, ss[k]) for k = 1:K]
	cs = [x->c_ik(i, x, ss[k]) for i = 1:length(supp), k = 1:K]
	cs, xs
end

###  
# Ideally should have a broadcast implementaiton that combines
# previous two...
function genNewsvendorsDiffSupp(supps, s, K)
	#Generic computation of the sth quantile 
	function x_k(p0, alpha, mhat_k, k, s) 
	    const Nhat_k = sum(mhat_k)
	    palpha = JS.shrink(mhat_k./Nhat_k, p0, alpha, Nhat_k)
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
