module InfOpt

# Import the necessary packages.
import JuMP
using MathOptInterface
const MOI = MathOptInterface

# Export user methods and datatypes.
export InfiniteModel, @infinite_variable, @point_variable, @global_variable,
GeneralVariableRef, FiniteVariableRef, InfiniteVariableRef, PointVariableRef,
GlobalVariableRef, InfOptVariable

# Import all of the datatpyes, methods, and definitions.
include("model.jl")
include("variables.jl")
include("expressions.jl")
include("operators.jl")
include("constraints.jl")
include("objective.jl")
include("macros.jl")
include("show.jl")

end # module
