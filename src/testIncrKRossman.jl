## A small driver file to run the 4 increasing K test using synthetic rossman data 
#and combine
#run like this
#julia -p 4 -L increasingKHarness.jl testIncrKRossman.jl numRuns outPathStub d N usePoisson
#passed arguments 
	#ARGS[1] is numRuns
	#ARGs[2] is partial path for output.  
	#ARGS[3] is d  
	#ARGS[4] is N
	#Args[5] is bool for using the poisson amount of data assumption.

const numRuns = parse(Int, ARGS[1])
const spath = ARGS[2]
const param_path = "../RossmanKaggleData/Results/"
const d = parse(Int, ARGS[3])
const N = parse(Int, ARGS[4])
const usePoisson = parse(Bool, ARGS[5])

const s = .95
K_grid = collect(100:100:1000)
outPath = "$(spath)_Ross_$(maximum(K_grid))_$(d)_$(N)_$(s)_$(usePoisson)_$(4*numRuns)"

#First read in the data and parse it appropriately
#Do this on a single processor bc it should be fast.
counts = readcsv("../RossmanKaggleData/Results/counts$(d).csv")
counts[ counts .== "NA" ] = 0  #R still dumps empties...
counts = convert(Array{Int64, 2}, counts[2:end, 2:end]')
ps_full = counts ./ sum(counts, 1)

supp_full = readcsv("../RossmanKaggleData/Results/support$(d).csv")
supp_full = convert(Array{Float64, 2}, supp_full[2:end, 2:end]')

tic()
file_a = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_a_", N, s, seed=8675309, usePoisson=usePoisson)
file_b = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_b_", N, s, seed=5167462266, usePoisson=usePoisson)
file_c = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_c_", N, s, seed=5164174290, usePoisson=usePoisson)
file_d = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_d_", N, s, seed=112456059, usePoisson=usePoisson)

# ######
file_a = fetch(file_a)
file_b = fetch(file_b)
file_c = fetch(file_c)
file_d = fetch(file_d)

time_stamp = toc()

##read everyone in, throw away a line
data, header = readcsv(file_a, header=true)

data_t = readcsv(file_b, skipstart=1)
data_t[:, 1] += numRuns
data = vcat(data, data_t)

data_t = readcsv(file_c, skipstart=1)
data_t[:, 1] += 2numRuns
data = vcat(data, data_t)

data_t = readcsv(file_d, skipstart=1)
data_t[:, 1] += 3numRuns
data = vcat(data, data_t)

f = open("$(outPath).csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)

println("Num Paths: \t $(4*numRuns) \t Time:", time_stamp )