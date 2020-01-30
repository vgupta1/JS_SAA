## A small driver file to run the 4 
# synthetic data tests in parallel
# data calibrated to the RossmanKaggleData
#Run like this
#julia -p 4 -L syntheticDataHarness.jl test_SyntheticData.jl numRuns outPathStub d usePoisson
#passed arguments 
	#ARGS[1] numRuns
	#ARGs[2] partial path for output.  
	#ARGS[3] d  
	#ARGS[4] bool for using the poisson amount of data assumption.

const numRuns = parse(Int, ARGS[1])
const spath = ARGS[2]
const param_path = "../RossmanKaggleData/CleanedData/"
const d = parse(Int, ARGS[3])
const usePoisson = parse(Bool, ARGS[4])

const s = .95
K_grid = vcat(round.(Int, 2 .^(3:.25:10)), 1115)
N_grid = [10]

outPath = "$(spath)_syntheticRossman_$(s)_$(usePoisson)_$(4*numRuns)"
#old naming specification.  keep for a bit.  
#outPath = "$(spath)_Ross_$(maximum(K_grid))_$(d)_$(N)_$(s)_$(usePoisson)_$(4*numRuns)"

#First read in the data and parse it appropriately
#Do this on a single processor bc it should be fast.
ps_full = readdlm("../RossmanKaggleData/CleanedData/ps_full$(d).csv", ',')
supp_full = readdlm("../RossmanKaggleData/CleanedData/support$(d).csv", ',')

start_time = time_ns()
file_a = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_a_", N_grid, s, seed=8675309, usePoisson=usePoisson)
file_b = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_b_", N_grid, s, seed=5167462266, usePoisson=usePoisson)
file_c = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_c_", N_grid, s, seed=5164174290, usePoisson=usePoisson)
file_d = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_d_", N_grid, s, seed=112456059, usePoisson=usePoisson)

# ######
file_a = fetch(file_a)
file_b = fetch(file_b)
file_c = fetch(file_c)
file_d = fetch(file_d)

time_stamp = (time_ns() - start_time) * 1e-9

##read everyone in, throw away a line
data, header = readdlm(file_a, ',', header=true)

data_t = readdlm(file_b, ',', skipstart=1)
data_t[:, 1] .+= numRuns
data = vcat(data, data_t)

data_t = readdlm(file_c, ',', skipstart=1)
data_t[:, 1] .+= 2numRuns
data = vcat(data, data_t)

data_t = readdlm(file_d, ',', skipstart=1)
data_t[:, 1] .+= 3numRuns
data = vcat(data, data_t)

f = open("$(outPath).csv", "w")
writedlm(f, header, ',')
writedlm(f, data, ',')
close(f)

println("Num Paths: \t $(4*numRuns) \t Time:", time_stamp )