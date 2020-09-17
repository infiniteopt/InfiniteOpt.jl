module Collections

# Import the necessary packages.
using ..InfiniteOpt
import JuMP
const JuMPC = JuMP.Containers

include("VectorTuple.jl")
include("DualDict.jl")

end # end module
