####Newsvendor Factories

#simple discrete quantile
function nv_quantile(pk, supp_k, s) 
	if s < 0
		return supp_k[1]
	end

	cdf_ = 0.
    for i = 1:length(pk)
        cdf_ += pk[i]
        if cdf_ >= s
        	return supp_k[i]
        end
    end
    return supp_k[end]
end

#returns a d x K matrix of univariate functions 
function getNewsVendorCosts(supps, s, K)
	function c_ik(i, k, x, s)
	    supps[i, k] > x ? s/(1-s) * (supps[i, k] - x) : (x - supps[i, k])
	end
	return [x->c_ik(i, k, x, s) for i = 1:size(supps, 1), k = 1:K]
end

#returns an K array of functions which train S-SAA
#hyperparsms for trainers:  shrinkage anchor and shrinakge amt
function genSSAAtrainers(supps, s, K)
	#Generic computation of the sth quantile 
	function x_k(mhat_k, k, s, p0, alpha) 
	    Nhat_k = sum(mhat_k)
	    palpha = JS.shrink(mhat_k./Nhat_k, p0, alpha, Nhat_k)
		nv_quantile(palpha, supps[:, k], s)	    
	end

	return [(mhat_k, (p0, alpha))-> x_k(mhat_k, k, s, p0, alpha) for k = 1:K]
end


function genKSTrainers(supps, s, K, method)
	#Computes KS Order quantity assuming Gamma small enough
	function x_k(mhat_k, k, s, Gamma::Number) 
	    Nhat_k = sum(mhat_k)
	    if Nhat_k == 0  #when there's no data, return maximum
	    	return supps[end, k]
	    end
	    @assert Gamma <= 1-s "KS Newsvendor not support Gamma > 1-s: $Gamma > $(1-s)"
	    phat = mhat_k./Nhat_k
	    qlo = nv_quantile(phat, supps[:, k], ceil(Nhat_k * (s - Gamma))/Nhat_k)
	    qhigh = nv_quantile(phat, supps[:, k], floor(Nhat_k * (s + Gamma) + 1)/Nhat_k)
	    return s * qhigh + (1 - s) * qlo
	end

	#does K-fold cross-validation to determine Gamma
	cs = getNewsVendorCosts(supps, s, K)
	function x_k(mhat_k, k, s, Gamma_grid, numFolds) 
	    Nhat_k = round(Int, sum(mhat_k))
	    #divy the data into cross-val sets
    	splits = mod.(1:Nhat_k, numFolds) .+ 1
    	shuffle!(splits)
    	cv_data = [zero(mhat_k) for fold = 1:numFolds]

    	ix_split = 1
        for i = 1:length(mhat_k)
            for j = 1:round(Int, mhat_k[i])
                cv_data[splits[ix_split]][i] += 1
                ix_split += 1            
            end
        end

	    #for each set train and evaluate
	    out = zero(Gamma_grid)
	    for (ix, Gamma) in enumerate(Gamma_grid)
		    for fold = 1:numFolds
		    	train_data = sum(cv_data[setdiff(1:numFolds, fold)])
		    	Nhat_train = sum(train_data)
		    	x = x_k(train_data, k, s, Gamma) 
	    		costs = [c(x) for c in cs[:, k]]
	    		out[ix] += dot(cv_data[fold], costs) / (Nhat_k - Nhat_train)
	    	end
	    	out[ix] /= numFolds
	    end

	    #identify who is best and retrain
	    Gammastar = Gamma_grid[argmin(out)]
	    #println(k, "\t", Gammastar)
	    x_k(mhat_k, k, s, Gammastar) 
	end

	if method == :aPriori
		xs  = [(mhat_k, Gamma)-> x_k(mhat_k, k, s, Gamma) for k = 1:K]
	elseif method == :crossVal
		xs = [(mhat_k, (Gamma_grid, numFolds)) -> x_k(mhat_k, k, s, Gamma_grid, numFolds) for k = 1:K]
	else
		throw("Method must be one of :aPriori or :crossVal")
	end

	xs
end
