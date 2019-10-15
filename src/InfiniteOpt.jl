module InfiniteOpt

# Import the necessary packages.
import JuMP
import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const JuMPC = JuMP.Containers
import Distributions
import Random

# Import all of the datatpyes, methods, macros, and definitions.
include("datatypes.jl")
include("parameters.jl")
include("variables.jl")
include("reduced_variables.jl")
include("expressions.jl")
include("measures.jl")
include("operators.jl")
include("constraints.jl")
include("objective.jl")
include("measure_expansions.jl")
include("macros.jl")
include("optimize.jl")
include("results.jl")
include("show.jl")
include("utilities.jl")

include("TranscriptionOpt/TranscriptionOpt.jl")
using .TranscriptionOpt

# Export model object datatype
export InfiniteModel

# Export macros and constants
export @infinite_variable, @point_variable, @global_variable, @infinite_parameter,
@BDconstraint, @finite_parameter, Infinite, Global, Point, Parameter

# Export infinite parameter datatypes
export InfOptParameter, ParameterRef

# Export variable datatypes
export InfOptVariable, InfiniteVariable, PointVariable, GlobalVariable,
GeneralVariableRef, InfiniteVariableRef, MeasureFiniteVariableRef,
FiniteVariableRef, PointVariableRef, GlobalVariableRef, InfOptVariableRef,
ReducedInfiniteVariableRef, AbstractReducedInfo, ReducedInfiniteInfo

# Export infinite parameter set types
export AbstractInfiniteSet, IntervalSet, DistributionSet

# Export infinite parameter methods
export build_parameter, add_parameter, infinite_set, set_infinite_set,
num_parameters, parameter_by_name, all_parameters, num_supports, has_supports,
set_supports, add_supports, delete_supports, supports, used_by_constraint,
used_by_measure, used_by_variable, is_used, group_id, is_independent,
set_independent, unset_independent, is_finite_parameter

# Export variable methods
export used_by_objective, infinite_variable_ref, parameter_refs,
set_parameter_refs, add_parameter_ref, used_by_point_variable, parameter_values,
eval_supports, used_by_reduced_variable

# Export expression datatypes
export InfiniteExpr, ParameterExpr, MeasureExpr

# Export constraint datatypes
export GeneralConstraintRef, InfiniteConstraintRef, MeasureConstraintRef,
FiniteConstraintRef, BoundedScalarConstraint

# Export measure datatypes
export AbstractMeasureData, DiscreteMeasureData, MultiDiscreteMeasureData,
Measure, MeasureRef

# Export measure methods
export add_measure, measure, measure_function, measure_data, expand,
expand_all_measures!

# Export transcription datatypes
export TranscriptionData, TranscriptionModel

# Export transcription methods
export is_transcription_model, transcription_data, transcription_variable,
transcription_constraint, transcription_model

# Export optimize methods
export optimizer_model, set_optimizer_model, optimizer_model_ready,
set_optimizer_model_ready, build_optimizer_model!, optimizer_model_key

# Export result query methods
export map_value, map_optimizer_index, map_dual, map_shadow_price

# Export functions that will be included in future JuMP releases


# Export printing methods
export bound_string

# Export support generation and fill-in methods
export fill_in_supports!, generate_supports

end # module
