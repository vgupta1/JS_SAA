#julia -L back_test_harness.jl back_test_rossman.jl outPathStub d
#passed arguments 
	#ARGS[1] is partial path for output.  
	#ARGS[2] is d must be one of 20, 50, 1000
	#ARGS[3] specifies the data stub:  one of AdjSales_NoWeekends, AdjSales_NoWeekends_RowShuffle, etc.

include("back_test_harness2.jl")

const spath = ARGS[1]
const param_path = ARGS[3] # e.g. "AdjSales_NoWeekends"
const d = parse(Int, ARGS[2])
const s = .95

@assert d in [20, 50, 1000] "Only bins of size 20, 50, 1000 currently supported"

K_grid = vcat(1, collect(10:10:90), collect(100:100:1000), 1115)
N_grid = [10, 20, 40]

outPath = "$(spath)_RossBacktest_$(maximum(K_grid))_$(d)_$(s)_$(maximum(N_grid))"

#First read in the data and parse it appropriately
ps_full = readcsv("../RossmanKaggleData/Results/ps_full$(d).csv")
supp_full = readcsv("../RossmanKaggleData/Results/support$(d).csv")

@assert minimum(ps_full) >=0 "ps_full has negative entries"

#Load shuffled data
tdata = readcsv("../RossmanKaggleData/Results/$(param_path)_Binned$(d).csv")

binned_data = tdata[2:end, 2:end]  #drop column header = stores, drop row header = dates
dates= tdata[2:end, 1] #keep track of dates for fun

#currently run single-threaded for ease
tic()
back_test2(K_grid, supp_full, ps_full, binned_data, dates, outPath, N_grid, s, onlySAA=false, numTestDays=10)
toc()

