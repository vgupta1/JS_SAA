## General purpose functions 

#User is expected to pass functions
#x_k(alpha, mhatk)  : provides a solution with data mhatk
#c_ki(x)      		
#These are often collected into arrays of functions
# xs			: [x_k for k = 1:K ]
# c_k  			: [c_ki for i = 1:d ]
# cs			: (d, K) Matrix with elements c_ki

#General purpose functions
###
#VG consider also write an inplace version. 
shrink(phat_k, p0, alpha, Nhat_k) = @. alpha * p0 / (Nhat_k + alpha) + Nhat_k * phat_k / (Nhat_k + alpha)

function z_k(x_k, c_k, alpha, mhat_k, ps_k, lam_k, lamavg)
	x = x_k(alpha, mhat_k)
	cs = map(c_ki -> c_ki(x), c_k)
	lam_k/lamavg * dot(ps_k, cs) 
end

#actually returns zLOO_k * N * lamavg
function zLOO_k_unsc(x_k, c_k, alpha, mhat_k)
	out = 0.
	mhatloo = mhat_k[:]
	for i = 1:length(mhat_k)
		mhatloo[i] -= 1	
		if i > 1 
			mhatloo[i - 1] += 1
		end
		x = x_k(alpha, mhatloo)
		out += mhat_k[i] * c_k[i](x)
	end
	out
end

#average true performance
function zbar(xs, cs, alpha, mhats, ps, lams)
	const lamavg = mean(lams)
	const K = size(cs, 2)
	out = 0.
	#maybe change this to a more numerically stable way to do avg?
	for k = 1:k
		out += z_k(xs[k], cs[:, k], alpha, mhats[:, k], ps[:, k], lams[k], lamavg)
	end
	out/K
end

function zLOObar_unsc(xs, cs, alpha, mhats)
	const K = size(cs, 2)
	out = 0.
	for k = 1:K
		out += zLOO_k_unsc(xs[k], cs[:, k], alpha, mhats[:, k])
	end
	out /K
end

#returns alphaOR, minimizingAlphaIndex, curveInAlpha
function oracle_alpha(xs, cs, mhats, ps, lams, alpha_grid)
	alphaOR = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zbar(xs, cs, alpha, mhats, ps, lams)
		if out[j] < best_val
			jstar = j
			alphaOR = alpha
		end
	end
	return alphaOR, jstar, out
end	

#returns alphaLOO, jstar, and cuveInAlpha_unsc
function loo_alpha(xs, cs, mhats, alpha_grid)
	alphaLOO = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zLOObar_unsc(xs, cs, alpha, mhats)
		if out[j] < best_val
			jstar = j
			alphaLOO = alpha
		end
	end
	return alphaLOO, jstar, out
end	
