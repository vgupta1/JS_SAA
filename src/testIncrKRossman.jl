## A small driver file to run the 4 increasing K test using synthetic rossman data 
#and combine
#run like this
#julia -p 4 -L increasingKHarness.jl testIncrKRossman.jl numRuns outPathStub d 
#passed arguments 
	#ARGS[1] is numRuns
	#ARGs[2] is partial path for output.  
	#ARGS[3] is 

const spath = ARGS[2]
const param_path = "../RossmanKaggleData/Results/"
const numRuns = parse(Int, ARGS[1])
const d = parse(Int, ARGS[3])

K_grid = collect(100:100:1000)
outPath = "$(spath)_Rossman_$(maximum(K_grid))_$(d)_$(4*numRuns)"

#First read in the data and parse it appropriately
#Do this on a single processor bc it should be fast.
counts = readcsv("../RossmanKaggleData/Results/counts$(d).csv")
counts[ counts .== "NA" ] = 0  #R still dumps empties...
counts = convert(Array{Int64, 2}, counts[2:end, 2:end]')
ps_full = counts ./ sum(counts, 1)

supp_full = readcsv("../RossmanKaggleData/Results/support$(d).csv")
supp_full = convert(Array{Float64, 2}, supp_full[2:end, 2:end]')

tic()
file_a = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_a_", seed=8675309)
file_b = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_b_", seed=5167462266)
file_c = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_c_", seed=5164174290)
file_d = @spawn convInKtest(numRuns, K_grid, supp_full, ps_full, "$(outPath)_d_", seed=112456059)

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

#strip the name of file_a to make the numbers better
# println("Filea \t", file_a)
# println("spath \t", spath)

# indx = search(file_a, spath)[end] + 1

f = open("$(outPath).csv", "w")
writecsv(f, header)
writecsv(f, data)
close(f)

println("Num Paths: \t $(4*numRuns) \t Time:", time_stamp )