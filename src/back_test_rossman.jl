#julia -L back_test_harness.jl back_test_rossman.jl outPathStub d
#passed arguments 
	#ARGS[1] is partial path for output.  
	#ARGS[2] is d must be one of 20, 50, 1000
	#ARGS[3] specifies the data stub:  one of AdjSales_NoWeekends, AdjSales_NoWeekends_RowShuffle, etc.

include("back_test_harness2.jl")

const spath = ARGS[1]
const param_path = ARGS[3] #"../RossmanKaggleData/Results/"
const d = parse(Int, ARGS[2])
const s = .95

K_grid = vcat(1, collect(10:10:90), collect(100:100:1000), 1115)
N_grid = [10, 20, 40]
#K_grid = [1115]
# K_grid = vcat(collect(100:100:1000), 1115)
N_grid = [10]

outPath = "$(spath)_RossBacktest_$(maximum(K_grid))_$(d)_$(s)_$(maximum(N_grid))"

#First read in the data and parse it appropriately
#Do this on a single processor bc it should be fast.
ps_full = readcsv("../RossmanKaggleData/Results/ps_full$(d).csv")
supp_full = readcsv("../RossmanKaggleData/Results/support$(d).csv")

@assert minimum(ps_full) >=0 "ps_full has negative entries"

#Load shuffled data
tdata = readcsv("../RossmanKaggleData/Results/$(param_path)_Binned$(d).csv")

binned_data = tdata[2:end, 2:end]  #column header = sotres, row header = dates
dates= tdata[2:end, 1] 

tic()
back_test2(K_grid, supp_full, ps_full, binned_data, dates, outPath, N_grid, s, onlySAA=false, numTestDays=10)
toc()

# tic()
# file_a = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_a_", N_grid, s, seed=8675309, usePoisson=usePoisson)
# file_b = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_b_", N_grid, s, seed=5167462266, usePoisson=usePoisson)
# file_c = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_c_", N_grid, s, seed=5164174290, usePoisson=usePoisson)
# file_d = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_d_", N_grid, s, seed=112456059, usePoisson=usePoisson)

# # ######
# file_a = fetch(file_a)
# file_b = fetch(file_b)
# file_c = fetch(file_c)
# file_d = fetch(file_d)

# time_stamp = toc()

# ##read everyone in, throw away a line
# data, header = readcsv(file_a, header=true)

# data_t = readcsv(file_b, skipstart=1)
# data_t[:, 1] += numRuns
# data = vcat(data, data_t)

# data_t = readcsv(file_c, skipstart=1)
# data_t[:, 1] += 2numRuns
# data = vcat(data, data_t)

# data_t = readcsv(file_d, skipstart=1)
# data_t[:, 1] += 3numRuns
# data = vcat(data, data_t)

# f = open("$(outPath).csv", "w")
# writecsv(f, header)
# writecsv(f, data)
# close(f)

# println("Num Paths: \t $(4*numRuns) \t Time:", time_stamp )