#julia -L rollingHarness.jl test_Rolling.jl outPathStub d datastub
#passed arguments 
	#ARGS[1] is partial path for output.  
	#ARGS[2] is d must be one of 20, 50, 1000
	#ARGS[3] specifies the data stub:  one of AdjSales_NoWeekends, AdjSales_NoWeekends_RowShuffle, etc.

using Distributed

const spath = ARGS[1]
const param_path = ARGS[3] # e.g. "AdjSales_NoWeekends"
const d = parse(Int, ARGS[2])
const s = .95

@assert d in [20, 50, 1000] "Only bins of size 20, 50, 1000 currently supported"

K_grid = vcat(round.(Int, 2 .^(3:.25:10)), 1115)
N_grid = [10, 20, 40]

outPath = "$(spath)_RossRolling_$(s)__$(d)_$(param_path)"

#First read in the data and parse it appropriately
ps_full = readdlm("../RossmanKaggleData/CleanedData/ps_full$(d).csv", ',')
supp_full = readdlm("../RossmanKaggleData/CleanedData/support$(d).csv", ',')

@assert minimum(ps_full) >=0 "ps_full has negative entries"

#Load shuffled data
tdata = readdlm("../RossmanKaggleData/CleanedData/$(param_path)_Binned$(d).csv", ',')

binned_data = tdata[2:end, 2:end]  #drop column header = stores, drop row header = dates
dates= tdata[2:end, 1] #keep track of dates for fun

#currently run single-threaded for ease
@elapsed rollingTest(K_grid, supp_full, ps_full, binned_data, dates, outPath, N_grid, s, onlySAA=false, numTestDays=10)

