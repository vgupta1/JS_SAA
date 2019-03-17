###
# Backtest using Rossman data
#
# For increasing K, 
#  Look at N consecutive days of data
#  Compute various methods
#  Record output by looking at performance on next numTestDays
###
using Distributions
include("../src/JS_SAA_main.jl")

#supp_full, ps_full are d x K matrices with true info per problem
#binned_data is just the binned_data, may have "NA"
#dates is the associated dates
#adds new subproblems in order appear in files
#s is the service level, N is the average amount of data per problem
function rollingTest(K_grid, supp_full, ps_full, binned_data_full, dates, outPath, N_grid, s;
					onlySAA = false, numTestDays=10)
	const Kmax = maximum(K_grid)
	@assert Kmax <= size(supp_full, 2) "K_grid exceeds available subproblems"
	@assert size(supp_full) == size(ps_full) "supp_full and ps_full have incompatible dimensions"
	const d = size(supp_full, 1)
	const numDataPoints = size(binned_data_full, 1)

	#For safety, trim inputs to size Kmax
	supp_full = view(supp_full, 1:d, 1:Kmax)
	ps_full = view(ps_full, 1:d, 1:Kmax)
	binned_data = view(binned_data_full, 1:numDataPoints, 1:Kmax) 
	dates = view(dates, 1:numDataPoints)

	p0 = ones(d) / d
	alpha_grid = range(0, stop=180, length=120)

	#set up output file
	f = open("$(outPath).csv", "w")
	writecsv(f, ["StartDate" "K" "d" "N" "Method" "TruePerf" "time" "alpha"])

	#generate all Kmax subproblems upfront and store in memory
	cs_full, xs_full = JS.genNewsvendorsDiffSupp(supp_full, s, Kmax)
	lam_full = ones(Kmax)

	for N in N_grid
		for ix_start = 1:N:numDataPoints
			#skip if partial window for training or test
			if ix_start + N + numTestDays - 1 > numDataPoints
				continue
			end

			#compute Training Data
			#within the window
			#Throw away NA's 
			mhats_full = zeros(Int, d, Kmax)
			for k = 1:Kmax
				for i = 1:N
					if binned_data[ix_start + i - 1, k] == "NA"
						continue
					end
					mhats_full[binned_data[ix_start + i - 1, k], k] += 1
				end
			end

			Nhats_full = vec(sum(mhats_full, 1))
			@assert length(Nhats_full) == Kmax	

			#build the out of sample training set
			mhats_out_full = zeros(Int, d, Kmax)
			for k = 1:Kmax			
				for i = 1:numTestDays
					if binned_data[ix_start + N + i - 1, k] == "NA"
						continue
					else
						mhats_out_full[binned_data[ix_start + N + i - 1, k], k] += 1
					end
				end
			end

			for K in K_grid
				println("($(N), $(ix_start), $(K))")
				#Take views on everything for simplicity
				lams = view(lam_full, 1:K)
				supp = view(supp_full, 1:d, 1:K)
				Nhats = view(Nhats_full, 1:K)
				mhats = view(mhats_full, 1:d, 1:K)
				cs = view(cs_full, 1:d, 1:K)
				xs = view(xs_full, 1:K)
				mhats_out = view(mhats_out_full, 1:d, 1:K)

				N_out = sum(mhats_out, 1)
				ps = mhats_out ./ N_out
				ps[:, vec(N_out) .== 0] = 0.

				#for data-driven shrinkage anchor
				#problems with no data do not contribute to the anchor
				phat_avg = JS.get_GM_anchor(mhats)
				if !isapprox(sum(phat_avg), 1) || minimum(phat_avg) < 0
					throw("Something went terriby wrong.  Grand avg phat is not a probability distribution")
				end

				#Compute the full-info value once for reference
				tic()
				full_info = JS.zstar(xs, cs, ps, lams)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "FullInfo" full_info t 0.])

				#SAA
				tic()
				perf_SAA = JS.zbar(xs, cs, p0, 0., mhats, ps, lams)

				t = toq()
				writecsv(f, [dates[ix_start] K d N "SAA" perf_SAA t 0.0])

				#Gen the Oracle cost with 1/d anchor
				tic()
				alphaOR, min_indx, or_alpha_curve = JS.oracle_alpha(xs, cs, mhats, ps, lams, p0, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "Oracle" or_alpha_curve[min_indx] t alphaOR])

				#Gen the LOO cost with 1/d anchor
				tic()
				alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, p0, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "LOO_unif" or_alpha_curve[min_indx] t alphaLOO])

				##MSE version of alpha
				tic()
				alphaMSE, min_indx = JS.mse_estimates(mhats, supp, p0, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "MSE" or_alpha_curve[min_indx] t alphaMSE])

				#Gen the Oracle cost with GM Anchor
				tic()
				alphaOR_GM, min_indx, or_alpha_curve_GM = JS.oracle_alpha(xs, cs, mhats, ps, lams, phat_avg, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "OraclePhat" or_alpha_curve_GM[min_indx] t alphaOR_GM])

				#Gen the LOO cost with the GM Anchor
				tic()
				alphaLOO, min_indx, looUnsc_curve = JS.loo_alpha(xs, cs, mhats, phat_avg, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "LOO_avg" or_alpha_curve_GM[min_indx] t alphaLOO])

				##MSE version of alpha with GM
				tic()
				alphaMSe, min_indx = JS.mse_estimates(mhats, supp, phat_avg, alpha_grid)
				t = toq()
				writecsv(f, [dates[ix_start] K d N "MSE" or_alpha_curve_GM[min_indx] t alphaMSE])

			end  #end K Loop
			flush(f)
		end #end window loop
	end #end N loop
	close(f)
	"$(outPath).csv"
end #end function 