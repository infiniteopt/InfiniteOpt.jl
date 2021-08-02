module Collections

# Import the necessary packages.
using ..InfiniteOpt
import JuMP, AbstractTrees
const JuMPC = JuMP.Containers

include("vectorize.jl")
include("VectorTuple.jl")
include("binarytrees.jl")

end # end module
