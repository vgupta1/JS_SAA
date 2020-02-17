#julia -p 10 -L rollingHarness.jl test_Rolling.jl outPathStub d datastub
#passed arguments 
	#ARGS[1] is partial path for output.  
	#ARGS[2] is d must be one of 20, 50, 1000
#Use this grid of N
#N_grid = [10, 20, 40]

using Distributed

const spath = ARGS[1]
const d = parse(Int, ARGS[2])
@assert d in [20, 50, 1000] "Only bins of size 20, 50, 1000 currently supported"

const param_path = "AdjSales_NoWeekends_ShuffleWithinCol"
const s = .95
K_grid = vcat(round.(Int, 2 .^(3:.25:10)), 1115)
outPath = "$(spath)_Ross_$(s)__$(d)_$(param_path)"

#First read in the data and parse it appropriately
ps_full = readdlm("../RossmanKaggleData/CleanedData/ps_full$(d).csv", ',')
supp_full = readdlm("../RossmanKaggleData/CleanedData/support$(d).csv", ',')

#Load shuffled data
tdata = readdlm("../RossmanKaggleData/CleanedData/$(param_path)_Binned$(d).csv", ',')

binned_data = tdata[2:end, 2:end]  #drop column header = stores, drop row header = dates
dates= tdata[2:end, 1] #keep track of dates for fun

#Separate into 10 processors.  
#Divy them
half_pts = floor(Int, size(binned_data, 1)/2)
start_time = time_ns()

#For N=10, 20 split K_grid into 2 and data into 2 (sacrifice 1 point)
file10_aa = @spawn rollingTest(K_grid[1:2:end], supp_full, ps_full, 
					 binned_data[1:half_pts, :], dates[1:half_pts], 
					 "$(outPath)_thread_10_a_a", 
					 [10], s)
file10_ab = @spawn rollingTest(K_grid[1:2:end], supp_full, ps_full, 
					 binned_data[(half_pts + 1):end, :], dates[(half_pts + 1):end], 
					 "$(outPath)_thread_10_a_b", 
					 [10], s)
file10_ba = @spawn rollingTest(K_grid[2:2:end], supp_full, ps_full, 
					 binned_data[1:half_pts, :], dates[1:half_pts], 
					 "$(outPath)_thread_10_b_a", 
					 [10], s)
file10_bb = @spawn rollingTest(K_grid[2:2:end], supp_full, ps_full, 
					 binned_data[(half_pts + 1):end, :], dates[(half_pts + 1):end], 
					 "$(outPath)_thread_10_b_b", 
					 [10], s)

file20_aa = @spawn rollingTest(K_grid[1:2:end], supp_full, ps_full, 
					 binned_data[1:half_pts, :], dates[1:half_pts], 
					 "$(outPath)_thread_20_a_a", 
					 [20], s)
file20_ab = @spawn rollingTest(K_grid[1:2:end], supp_full, ps_full, 
					 binned_data[(half_pts + 1):end, :], dates[(half_pts + 1):end], 
					 "$(outPath)_thread_20_a_b", 
					 [20], s)
file20_ba = @spawn rollingTest(K_grid[2:2:end], supp_full, ps_full, 
					 binned_data[1:half_pts, :], dates[1:half_pts], 
					 "$(outPath)_thread_20_b_a", 
					 [20], s)
file20_bb = @spawn rollingTest(K_grid[2:2:end], supp_full, ps_full, 
					 binned_data[(half_pts + 1):end, :], dates[(half_pts + 1):end], 
					 "$(outPath)_thread_20_b_b", 
					[20], s)

#N = 40, just split K = 2 processors
file40_a = @spawn rollingTest(K_grid[1:2:end], supp_full, ps_full, 
					 binned_data, dates, 
					 "$(outPath)_thread_40_a", 
					[40], s)
file40_b = @spawn rollingTest(K_grid[2:2:end], supp_full, ps_full, 
					 binned_data, dates, 
					 "$(outPath)_thread_40_b", 
					[40], s)

# ######
file10_aa = fetch(file10_aa)
file10_ab = fetch(file10_ab)
file10_ba = fetch(file10_ba)
file10_bb = fetch(file10_bb)

file20_aa = fetch(file20_aa)
file20_ab = fetch(file20_ab)
file20_ba = fetch(file20_ba)
file20_bb = fetch(file20_bb)

file40_a = fetch(file40_a)
file40_b = fetch(file40_b)

time_stamp = (time_ns() - start_time) * 1e-9

##read everyone in, throw away a line
data_10_aa, header = readdlm(file10_aa, ',', header=true)
data_10_ab, header = readdlm(file10_ab, ',', header=true)
data_10_ba, header = readdlm(file10_ba, ',', header=true)
data_10_bb, header = readdlm(file10_bb, ',', header=true)

data_20_aa, header = readdlm(file20_aa, ',', header=true)
data_20_ab, header = readdlm(file20_ab, ',', header=true)
data_20_ba, header = readdlm(file20_ba, ',', header=true)
data_20_bb, header = readdlm(file20_bb, ',', header=true)

data_40_a, header = readdlm(file40_a, ',', header=true)
data_40_b, header = readdlm(file40_b, ',', header=true)

data = reduce(vcat, (data_10_aa, data_10_ab, data_10_ba, data_10_bb, 
					 data_20_aa, data_20_ab, data_20_ba, data_20_bb, 
					 data_40_a, data_40_b) )

f = open("$(outPath).csv", "w")
writedlm(f, header, ',')
writedlm(f, data, ',')
close(f)

println("Time:", time_stamp )
println("$(outPath).csv")