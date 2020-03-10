module InfiniteOpt

# Import the necessary packages.
import JuMP
import MathOptInterface
import Distributions
import Random
import MutableArithmetics
const _MA = MutableArithmetics
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const JuMPC = JuMP.Containers

# Import all of the datatpyes, methods, macros, and definitions.
include("datatypes.jl")
include("infinite_sets.jl")
include("parameters.jl")
include("variables.jl")
include("reduced_variables.jl")
include("expressions.jl")
include("mutable_arithmetics.jl")
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

include("MeasureEvalMethods/MeasureEvalMethods.jl")
using .MeasureEvalMethods: generate_measure_data, mc_sampling, gauss_legendre,
gauss_hermite, gauss_laguerre, infinite_transform, generate_supports_and_coeffs,
default_set_types, default_methods, eval_method_registry, register_eval_method

# Export the measure eval functions
export generate_measure_data, mc_sampling, gauss_legendre, gauss_hermite,
gauss_laguerre, infinite_transform, measure_dispatch,
generate_supports_and_coeffs, eval_method_registry, register_eval_method

# Export model object datatype
export InfiniteModel

# Export macros
export @infinite_variable, @point_variable, @hold_variable, @infinite_parameter,
@BDconstraint, @finite_parameter, @set_parameter_bounds, @add_parameter_bounds

# Export constants
export Infinite, Hold, Point, Parameter, Sampling, Quad, default_set_types,
default_methods

# Export infinite parameter datatypes
export InfOptParameter, ParameterRef

# Export infinite parameter set types
export AbstractInfiniteSet, IntervalSet, DistributionSet

# Export infinite parameter methods
export build_parameter, add_parameter, infinite_set, set_infinite_set,
num_parameters, parameter_by_name, all_parameters, num_supports, has_supports,
set_supports, add_supports, delete_supports, supports, used_by_constraint,
used_by_measure, used_by_variable, is_used, group_id, is_independent,
set_independent, unset_independent, is_finite_parameter, supports_in_set,
generate_and_add_supports!, generate_support_values

# Export variable datatypes
export InfOptVariable, InfiniteVariable, PointVariable, HoldVariable,
GeneralVariableRef, InfiniteVariableRef, MeasureFiniteVariableRef,
FiniteVariableRef, PointVariableRef, HoldVariableRef, InfOptVariableRef,
ReducedInfiniteVariableRef, AbstractReducedInfo, ReducedInfiniteInfo,
ParameterBounds

# Export variable methods
export used_by_objective, infinite_variable_ref, parameter_refs,
set_parameter_refs, add_parameter_ref, used_by_point_variable, parameter_values,
eval_supports, used_by_reduced_variable, has_parameter_bounds, parameter_bounds,
set_parameter_bounds, add_parameter_bound, delete_parameter_bound,
delete_parameter_bounds

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
expand_all_measures!, expect, support_sum

# Export measure evaluation functions
export generate_measure_data, MC_sampling, Gauss_Legendre, Gauss_Hermite,
Gauss_Laguerre, infinite_transform, support_formatting, measure_dispatch,
set_measure_defaults, measure_defaults

# Export transcription datatypes
export TranscriptionData, TranscriptionModel

# Export transcription methods
export is_transcription_model, transcription_data, transcription_variable,
transcription_constraint, transcription_model

# Export optimize methods
export optimizer_model, set_optimizer_model, optimizer_model_ready,
set_optimizer_model_ready, build_optimizer_model!, optimizer_model_key,
optimizer_model_variable, optimizer_model_constraint

# Export result query methods
export map_value, map_optimizer_index, map_dual, map_shadow_price,
map_lp_rhs_perturbation_range, map_lp_objective_perturbation_range

# Export functions that will be included in future JuMP releases


# Export printing methods
export bound_string

# Export support generation and fill-in methods
export fill_in_supports!, generate_supports

end # module
