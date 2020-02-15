### JS_SAA main file
module JS

using Distributions, StatsBase, LinearAlgebra, Random

#These used by the methods that optimize the anchor
using Optim, Clustering

include("genPurpose.jl")
include("deprecated.jl")
include("newsvendorfactory.jl")
include("babyNewsvendor.jl")

end #ends module 