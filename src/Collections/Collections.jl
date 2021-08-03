module Collections

# Import the necessary packages.
using ..InfiniteOpt
import JuMP
const JuMPC = JuMP.Containers

include("vectorize.jl")
include("VectorTuple.jl")

end # end module
