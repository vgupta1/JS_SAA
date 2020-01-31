### JS_SAA main file
module JS

using Distributions, StatsBase, LinearAlgebra, Random

include("genPurpose.jl")
include("newsvendorfactory.jl")
include("babyNewsvendor.jl")

end #ends module 