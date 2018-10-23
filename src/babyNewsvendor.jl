### 
#Workhorse functions for Stein Shrinkage SAA paper
###
module JS
using Distributions

###
#General purpose functions
###
#VG consider also write an inplace version. 
shrink(phat_k, p0, alpha, Nhat_k) = @. alpha * p0 / (Nhat_k + alpha) + Nhat_k * phat_k / (Nhat_k + alpha)


###
# Specialized functions for Newsvendor with Ber(p) demand
# unit holding, backorder cost: b = s/(1-s)
# Because of Bernoulli demand, represent data as mhat = count zeta ==1
# xk = 1 iff p(alpha) > 1-s

# prob xk = 1
nv_prob_k(p, p0, alpha, Nhat_k, s) = ccdf(Binomial(Nhat_k, p), (Nhat_k + alpha) * (1 - s) - alpha * p0)
    
#vectorized way to simulate counts
#returns vector mhat_k for 1s only.  
sim_path(ps, Ns) = @.sum(rand(Bernoulli(ps), Ns))

#E[Z_k(alpha)], i.e. expected true obj at x(alpha)
function exp_nv_obj_k(p, p0, alpha, Nhat_k, s)
    const b = s / (1 - s)
    (1 - p - b * p) * nv_prob_k(p, p0, alpha, Nhat_k, s) + b * p
end

#expected value of fullInfo solution
function exp_nv_fullInfo_k(p, s) 
	const b = s / (1 - s)
	const xstar = p > 1 - s ? 1. : 0.
	b * p + (1 - p - b * p) * xstar
end

#Z_k(alpha)..  mhat is # zeta-k == 1
function nv_obj_k(mhat, p, p0, alpha, Nhat_k, s)
	const p_alpha = shrink(mhat / Nhat_k, p0, alpha, Nhat_k)
	if p_alpha > 1 - s #xk = 1
		return 1 - p
	else #xk = 0
		return s / (1 - s) * p
	end
	return 0.  #never reached
end

#Computes Nint * ZLoo_k(alpha)
#Nint is the overall intensity for observations
#Assumes lambdas = 1
function nv_loo_k(mhat_k, p0, alpha, Nhat_k, s)
	out = 0.	
	#i corresponding to demand = 1
	if mhat_k > 0
		p_alpha = shrink((mhat_k - 1)/(Nhat_k - 1), p0, alpha, Nhat_k - 1) 
		#only non-zero if x neq 1
		out += p_alpha > 1 - s ? 0. : mhat_k * s / (1 - s) 
	end

	#i corresponding to demand = 0
	if mhat_k < Nhat_k
		p_alpha = shrink(mhat_k / (Nhat_k - 1), p0, alpha, Nhat_k - 1)
		#only non-zero if x  = 1
		out += p_alpha > 1 - s ? (Nhat_k - mhat_k) : 0.
	end
	return out
end

#this is a funny saa-ish term in the Loo-Stability decomposition
#Returns Nint * saa-ish term.  
function nv_saa_k(mhat_k, p0, alpha, Nhat_k, s)
	const p_alpha = shrink(mhat_k/Nhat_k, p0, alpha, Nhat_k)
	if p_alpha > 1 - s #xk = 1
		return Nhat_k - mhat_k
	else #xk = 0
		return s / (1 - s) * mhat_k
	end
	return 0.  #never reached
end


#Optimzes Z(\alpha) over alpha_grid
#returns alphaOR, minimizingAlphaIndex, curveInAlpha
function oracle_alpha(mhats, ps, p0, alpha_grid, Nhats, ss)
	out = map(a-> nv_obj(mhats, ps, p0, a, Nhats, ss), alpha_grid)
	jstar = indmin(out)
	return alpha_grid[jstar], jstar, out
end

#returns p0OR, alphaOR, jOR for alpha
function oracle_both(mhats, ps, p0_grid, alpha_grid, Nhats, ss)
	out = zeros(length(alpha_grid))
	or_perf, p0OR, alphaOR, jOR = Inf, 0., 0., 0
	for p0 in p0_grid
		out[:] = map(a-> nv_obj(mhats, ps, p0, a, Nhats, ss), alpha_grid)
		jstar = indmin(out)
		if out[jstar] < or_perf
			jOR = jstar
			or_perf = out[jstar]
			p0OR = p0
			alphaOR = alpha_grid[jstar]
		end
	end
	return p0OR, alphaOR, jOR
end

#optimizes ZLoo(\alpha)over alpha_grid
#Note curve in alpha = Nintensity * lamavg * ZLoo(alpha)
#returns alpha^LOO, maximizingAlphaIndex, curveInAlpha
function loo_alpha(mhats, p0, alpha_grid, Nhats, ss)
	out = map(a-> nv_loo(mhats, p0, a, Nhats, ss), alpha_grid)
	jstar = indmin(out)
	return alpha_grid[jstar], jstar, out
end

#loo estimate optimizing boht p0 and alpha0
function loo_both(mhats, p0_grid, alpha_grid, Nhats, ss)
	out = zeros(length(alpha_grid))
	loo_perf = Inf
	p0LOO = 0.
	alphaLOO = 0.
	jLOO = -1
	for p0 in p0_grid
		out[:] = map(a-> nv_loo(mhats, p0, a, Nhats, ss), alpha_grid)
		jstar = indmin(out)
		if out[jstar] < loo_perf
			jLOO = jstar
			loo_perf = out[jstar]
			p0LOO = p0
			alphaLOO = alpha_grid[jstar]
		end
	end
	return p0LOO, alphaLOO, jLOO
end



#Optimzes E[Z(\alpha)] over alpha_grid
#returns alphaAP, minimizingAlphaIndex, curveInAlpha
function apriori_alpha(ps, p0, alpha_grid, Nhats, ss)
	out = map(a-> exp_nv_obj(ps, p0, a, Nhats, ss), alpha_grid)
	jstar = indmin(out)
	return alpha_grid[jstar], jstar, out
end


#vector versions assume lam_k = 1 for all k
nv_loo(mhats, p0, alpha, Nhats, ss) =  mean( nv_loo_k.(mhats, p0, alpha, Nhats, ss) )
nv_obj(mhats, ps, p0, alpha, Nhats, ss) = mean( nv_obj_k.(mhats, ps, p0, alpha, Nhats, ss) )
nv_saa(mhats, p0, alpha, Nhats, ss) = mean( nv_saa_k.(mhats, p0, alpha, Nhats, ss) )

exp_nv_obj(ps, p0, alpha, Nhats, ss) = mean( exp_nv_obj_k.(ps, p0, alpha, Nhats, ss) )
nv_var_obj(ps, p0, alpha, Nhats, ss) = mean( nv_var_obj_k.(ps, p0, alpha, Nhats, ss) )
exp_nv_fullInfo(ps, s) = mean(exp_nv_fullInfo_k.(ps, s))

end #closes module