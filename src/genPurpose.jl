## General purpose functions 

#User is expected to pass functions
#x_k(mhatk, hyperparam)  : provides a solution with data mhatk using hyperparam (e.g., p0, alpha)
#c_ik(x)      		
#These are often collected into arrays of functions
# xs			: [x_k for k = 1:K]
# c_k  			: [c_ik for i = 1:d]
# cs			: (d, K) Matrix with elements c_ik

#General purpose functions
###
function shrink(phat_k, p0, alpha, Nhat_k; alphaMax = 1e10)
	if Nhat_k == 0 || alpha > alphaMax
		return p0
	end
	@. alpha * p0 / (Nhat_k + alpha) + Nhat_k * phat_k / (Nhat_k + alpha)
end

#need to be careful because some problems may have zero data
function get_GM_anchor(mhats)
	Nhats = vec(sum(mhats, dims=1))
	non_zero_indx = Nhats .> 0
	vec(mean(mhats[:, non_zero_indx] ./ Nhats[non_zero_indx]', dims=2))
end

function sim_path(p_k, N::Int)
    mhats_k = zeros(length(p_k))
    
    if N == 0 #allow simulating zero data points
    	return mhats_k
    end

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
	Nhats = sum(mhats, dims=1)
	p0 = mean(mhats ./ Nhats, dims=2)
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
	Nhats = sum(mhats, dims=1)
	phats = mhats ./ Nhats
	mu0s = vec(p0' * supps)

	#summary statistics
	muhats = Vector{Float64}(undef, K)
	sigmas = Vector{Float64}(undef, K)
	mse0 = Vector{Float64}(undef, K)
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

#kth element of true perf
function z_k(x_k, c_k, mhat_k, ps_k, lam_k, lamavg, hyper_param)
	x = x_k(mhat_k, hyper_param)
	cs = map(c_ik -> c_ik(x), c_k)
	lam_k/lamavg * dot(ps_k, cs) 
end

#actually returns zLOO_k * N * lamavg
#hyper_param = p0, alpha
function zLOO_k_unsc(x_k, c_k, mhat_k, hyper_param)
	out = 0.
	mhatloo = mhat_k[:] #copy data
	#only compute for terms where mhat_k[i] > 0
	#compute a base solution to resuse as necessary
	x_base = x_k(mhat_k, hyper_param)
	for i = 1:length(mhat_k)
		#correct previous toggle
		if i > 1 && mhat_k[i - 1] > 0
			mhatloo[i - 1] += 1
		end

		if mhat_k[i] > 0
			mhatloo[i] -= 1	
		    x = x_k(mhatloo, hyper_param)
		else 
			x = x_base
		end
		out += mhat_k[i] * c_k[i](x)
	end
	out
end


#average true performance
#hyper_param = p0, alpha
function zbar(xs, cs, mhats, ps, lams, hyper_param)
	lamavg = mean(lams)
	K = size(cs, 2)
	out = 0.
	#VG maybe change this to a more numerically stable way to do avg?
	for k = 1:K
		out += z_k(xs[k], cs[:, k], mhats[:, k], ps[:, k], lams[k], lamavg, hyper_param)
	end
	out/K
end

#full-info for scaling/comparison
#uses the fact that data scale doesn't matter.
#Only works for the regular xs (which is lame)
zstar(xs, cs, ps, lams) = zbar(xs, cs, ps, ps, lams, (ps[:, 1], 0.))

#hyper_param = p0, alpha,
function zLOObar_unsc(xs, cs, mhats, hyper_param)
	K = size(cs, 2)
	out = 0.
	for k = 1:K
		out += zLOO_k_unsc(xs[k], cs[:, k], mhats[:, k], hyper_param)
	end
	out /K
end

#actually returns N lam_avg * ZCVbar
#hyper_param = p0, alpha 
function zCVbar_unsc(xs, cs, cv_data, hyper_param)
	out = 0.
	numFolds = length(cv_data)
	lams = ones(length(xs))

	for fold = 1:numFolds
		#form the training set and testing distribution
		#be careful about empty problems in testing
		train_data = sum(cv_data[setdiff(1:numFolds, fold)])

		#dispatch to the zbar method, again using scaling not important
		out += zbar(xs, cs, train_data, cv_data[fold], lams, hyper_param)
	end
	return out
end

#returns alphaOR, minimizingAlphaIndex, curveInAlpha
function oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
	alphaOR = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zbar(xs, cs, mhats, ps, lams, (p0, alpha))
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
		out[j] = zLOObar_unsc(xs, cs, mhats, (p0, alpha))
		if out[j] < best_val
			jstar = j
			alphaLOO = alpha
			best_val = out[j]
		end
	end
	return alphaLOO, jstar, out
end	

function cv_alpha(xs, cs, mhats, p0, alpha_grid, numFolds)
	#divy up the data
	cv_data = split_cv(mhats, numFolds)

	alphaCV = 0.
	jstar = -1
	best_val = Inf
	out = zeros(length(alpha_grid))
	for (j, alpha) in enumerate(alpha_grid)
		out[j] = zCVbar_unsc(xs, cs, cv_data, (p0, alpha))
		if out[j] < best_val
			jstar = j
			alphaCV = alpha
			best_val = out[j]
		end
	end
	return alphaCV, jstar, out
end	

#maps ys to simplex somewhat safely
function soft_max(ys)
	ymax = maximum(ys)
	weights = exp.(ys .- ymax)
	weights /= sum(weights)
	@assert isprobvec(weights) "soft_max fail: \n $ys \n $weights"

	weights
end

#optimizes choice of anchor and alpha by approx minimizing LOO
#Heuristic usese particle swarm a 2 starts. 
#Passing numClusters = -1 makes anchor a linear comb of all phats
#return anchor, alpha, loo val
function loo_anchor(xs, cs, mhats; numClusters = 20, init_sqrt_alpha = 1.,
                    time_limit = 60., iterations=1000, store_trace=false, info=false)
    #First cluster to find points.
    phats = mhats ./ sum(mhats, dims=1)
    if numClusters < 0
    	mix_comp = phats
    	numClusters = size(phats, 2)
    else
	    cluster_res = kmeans(phats, numClusters; maxiter=50)
    	mix_comp = cluster_res.centers       
	end    

    #add the grandmean for ease
    mix_comp = hcat(mix_comp, JS.get_GM_anchor(mhats))

    # #### DEBUG
    # #confirm that everyone in mix_comp isn't dumb
    # for ix = 1:size(mix_comp, 2)
    # 	@assert isprobvec(mix_comp[:, ix]) "Not a probability component"  
    # end
    # if sum( mix_comp .== NaN ) > 0
    # 	throw("One of the mixture components has a NaN")
    # end
    # if sum( mix_comp .== Inf ) > 0
    # 	throw("One of the mix components as an Inf")
    # end
    # println("Mix Comp Passed Test")
    # #END DEBUG


    #write aux function 
    function f(ys)
        #use softmax to ensure simplex
        weights = soft_max(ys[1:end-1])
        p0 = vec(mix_comp * weights)        
        alpha = (ys[end])^2
        @assert isprobvec(p0) "P0 Failed here: \n $weights \n $ys[1:end-1]"
        JS.zLOObar_unsc(xs, cs, mhats, (p0, alpha))            
    end
    
    #optimize with two starting points and take best one
    x01 = [ones(numClusters + 1); init_sqrt_alpha]
    init_val1 = f(x01)
    @time res1 = optimize(f, x01, 
                    ParticleSwarm(n_particles=10), 
                    Optim.Options(time_limit=time_limit, iterations=iterations, store_trace=store_trace))
    z1 = Optim.minimum(res1)
    
    x02 = [zeros(numClusters); 1.; init_sqrt_alpha]
    init_val2 = f(x02)
    @time res2 = optimize(f, x02, 
                    ParticleSwarm(n_particles=10), 
                    Optim.Options(time_limit=time_limit, iterations=iterations, store_trace=store_trace))
    z2 = Optim.minimum(res2)
    
    if z1 < z2    
        xstar = Optim.minimizer(res1)
        zstar = z1
    else
        xstar = Optim.minimizer(res2)
        zstar = z2
    end
    info && println("Perc: Improv over LOO GM:\t", 1- zstar / init_val2)
    
    weights = soft_max(xstar[1:end-1])
    p0 = vec(mix_comp * weights)       
    alpha = (xstar[end])^2
    p0, alpha, zstar
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

###Helper function for splitting data for k-fold cross-validation
#This is not a memory efficient implementation
#Just lazy split/iterations
function split_cv(mhats, numFolds)
	d, K = size(mhats)
    Ntot = convert(Int, sum(mhats))
    splits = mod.(1:Ntot, numFolds) .+ 1
    shuffle!(splits)
    cv_data = [zero(mhats) for fold = 1:numFolds]

    #iterate through each k in order
    ix_split = 1

    for k = 1:K
        for i = 1:d
            for j = 1:mhats[i, k]
                cv_data[splits[ix_split]][i, k] += 1
                ix_split += 1            
            end
        end
    end
    return cv_data
end



