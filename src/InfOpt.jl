module InfOpt

# Import the necessary packages.
import JuMP
import MathOptInterface
const MOI = MathOptInterface
import Distributions

# Export model object datatype
export InfiniteModel

# Export macros and constants
export @infinite_variable, @point_variable, @global_variable, @infinite_parameter,
Infinite, Global, Point, Parameter

# Export infinite parameter datatypes
export InfOptParameter, ParameterRef

# Export variable datatypes
export InfOptVariable, InfiniteVariable, PointVariable, GlobalVariable,
GeneralVariableRef, InfiniteVariableRef, MeasureFiniteVariableRef,
FiniteVariableRef, PointVariableRef, GlobalVariableRef, InfOptVariableRef

# Export infinite parameter set types
export AbstractInfiniteSet, IntervalSet, DistributionSet, DiscreteSet

# Export infinite parameter functions
export get_parameter_ref, set_parameter_refs, add_parameter_ref

# Export constraint datatypes
export GeneralConstraintRef, InfiniteConstraintRef, MeasureConstraintRef,
FiniteConstraintRef

# Export measure datatypes
export AbstractMeasureData, DiscretizeData, Measure, MeasureRef

# Export measure methods
export add_measure, measure

# Import all of the datatpyes, methods, and definitions.
include("datatypes.jl")
include("parameters.jl")
include("variables.jl")
include("measure.jl")
include("expressions.jl")
include("operators.jl")
include("constraints.jl")
include("objective.jl")
include("macros.jl")
include("show.jl")

end # module
