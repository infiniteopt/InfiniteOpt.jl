module InfOpt

# Import the necessary packages.
import JuMP
using MathOptInterface
const MOI = MathOptInterface

# Export user methods and datatypes.
export InfiniteModel

# Import all of the datatpyes, methods, and definitions.
include("model.jl")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("show.jl")

end # module
