### Deprecated functions to be excised eventually

#maps ys to simplex somewhat safely
function soft_max(ys)
	ymax = maximum(ys)
	weights = exp.(ys .- ymax)
	weights /= sum(weights)
	@assert isprobvec(weights) "soft_max fail: \n $ys \n $weights"

	weights
end

###########

###########
#  Some functions to optimize the anchors
function genMixComp(mhats, numClusters; maxiter=50)
    #First cluster to find points.
    phats = mhats ./ sum(mhats, dims=1)
    if numClusters < 0
    	mix_comp = phats
    	numClusters = size(phats, 2)
    else
	    cluster_res = kmeans(phats, numClusters; maxiter=maxiter)
    	mix_comp = cluster_res.centers       
	end    

    #add the grandmean for ease
    hcat(mix_comp, JS.get_GM_anchor(mhats))
end

#optimizes choice of anchor and alpha by approx minimizing LOO
#Heuristic usese particle swarm a 2 starts. 
#Passing numClusters = -1 makes anchor a linear comb of all phats
#return anchor, alpha, loo val
function loo_anchor(xs, cs, mhats; numClusters = 20, init_sqrt_alpha = 1.,
                    time_limit = 60., iterations=1000, store_trace=false, info=false)
	mix_comp = genMixComp(mhats, numClusters)

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
    x01 = [ones(size(mix_comp, 2)); init_sqrt_alpha]
    init_val1 = f(x01)
    @time res1 = optimize(f, x01, 
                    ParticleSwarm(n_particles=10), 
                    Optim.Options(time_limit=time_limit, iterations=iterations, store_trace=store_trace))
    z1 = Optim.minimum(res1)
    
    x02 = [zeros(size(mix_comp, 2) - 1); 1.; init_sqrt_alpha]
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

#optimizes choice of anchor and alpha by approx minimizing oracle criteria
#Heuristic usese particle swarm a 2 starts. 
#Passing numClusters = -1 makes anchor a linear comb of all phats
#return anchor, alpha, loo val
function opt_oracle_anchor(xs, cs, ps, mhats; numClusters = 20, init_sqrt_alpha = 1.,
                    time_limit = 60., iterations=1000, store_trace=false, info=false)
	mix_comp = genMixComp(mhats, numClusters)
	lams = ones(length(xs))
    #write aux function 
    function f(ys)
        #use softmax to ensure simplex
        weights = soft_max(ys[1:end-1])
        p0 = vec(mix_comp * weights)        
        alpha = (ys[end])^2
        @assert isprobvec(p0) "P0 Failed here: \n $weights \n $ys[1:end-1]"
        JS.zbar(xs, cs, mhats, ps, lams, (p0, alpha))
    end
    
    #optimize with two starting points and take best one
    x01 = [ones(size(mix_comp, 2)); init_sqrt_alpha]
    init_val1 = f(x01)
    @time res1 = optimize(f, x01, 
                    ParticleSwarm(n_particles=10), 
                    Optim.Options(time_limit=time_limit, iterations=iterations, store_trace=store_trace))
    z1 = Optim.minimum(res1)
    
    x02 = [zeros(size(mix_comp, 2) - 1); 1.; init_sqrt_alpha]
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
