module InfiniteOpt

# Import the necessary packages.
import JuMP
import MathOptInterface
import Distributions
import DataStructures

# Make useful aliases
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const JuMPC = JuMP.Containers
const MOIUC = MOIU.CleverDicts

# Import the Collections module
include("Collections/Collections.jl")
using .Collections: VectorTuple, same_structure, DualDict

# Export the collections and methods  from Collections
export VectorTuple, same_structure, DualDict

# Import all of the datatpyes, methods, macros, and definitions.
include("datatypes.jl")
include("infinite_sets.jl")
include("general_variables.jl")
# include("parameters.jl")
# include("variable_basics.jl")
# include("infinite_variables.jl")
# include("point_variables.jl")
# include("hold_variables.jl")
# include("reduced_variables.jl")
# include("expressions.jl")
# include("measures.jl")
# include("constraints.jl")
# include("objective.jl")
# include("measure_expansions.jl")
# include("macros.jl")
# include("optimize.jl")
# include("results.jl")
# include("show.jl")
include("utilities.jl")

# Import the other submodules
include("TranscriptionOpt/TranscriptionOpt.jl")
using .TranscriptionOpt

# include("MeasureEvalMethods/MeasureEvalMethods.jl")
# using .MeasureEvalMethods: generate_measure_data, mc_sampling, gauss_legendre,
# gauss_hermite, gauss_laguerre, infinite_transform, generate_supports_and_coeffs,
# trapezoid

# Export the measure eval functions and constants
export generate_measure_data, mc_sampling, gauss_legendre, gauss_hermite,
gauss_laguerre, trapezoid, infinite_transform, measure_dispatch,
generate_supports_and_coeffs

# Export core datatypes
export AbstractDataObject, InfiniteModel

# Export macros
export @infinite_variable, @point_variable, @hold_variable, @infinite_parameter,
@BDconstraint, @finite_parameter, @set_parameter_bounds, @add_parameter_bounds

# Export constants
export Infinite, Hold, Point, Parameter, sampling, quadrature

# Export index types
export AbstractInfOptIndex, ObjectIndex, IndependentParameterIndex,
DependentParametersIndex, DependentParameterIndex, FiniteParameterIndex,
InfiniteVariableIndex, PointVariableIndex, ReducedInfiniteVariableIndex,
HoldVariableIndex, MeasureIndex, ConstraintIndex, InfiniteParameterIndex

# Export infinite set types
export AbstractInfiniteSet, InfiniteScalarSet, InfiniteArraySet, IntervalSet,
UniDistributionSet, MultiDistributionSet, CollectionSet

# Export infinite set methods
export collection_sets, supports_in_set, generate_support_values

# Export parameter datatypes
export InfOptParameter, ScalarParameter, IndependentParameter, FiniteParameter,
DependentParameters, ScalarParameterData, MultiParameterData,
IndependentParameterRef, DependentParameterRef, FiniteParameterRef

# Export parameter methods
export build_parameter, add_parameter, infinite_set, set_infinite_set,
num_parameters, parameter_by_name, all_parameters, num_supports, has_supports,
set_supports, add_supports, delete_supports, supports, used_by_constraint,
used_by_measure, used_by_infinite_variable, is_used, generate_and_add_supports!

# Export variable datatypes
export InfOptVariable, InfiniteVariable, ReducedInfiniteVariable,
PointVariable, HoldVariable, ParameterBounds, VariableData, GeneralVariableRef,
DispatchVariableRef, InfiniteVariableRef, ReducedInfiniteVariableRef,
MeasureFiniteVariableRef, FiniteVariableRef, PointVariableRef, HoldVariableRef,
DecisionVariableRef, UserDecisionVariableRef

# Export variable methods
export used_by_objective, infinite_variable_ref, parameter_refs,
set_parameter_refs, add_parameter_ref, used_by_point_variable, parameter_values,
eval_supports, used_by_reduced_variable, has_parameter_bounds, parameter_bounds,
set_parameter_bounds, add_parameter_bound, delete_parameter_bound,
delete_parameter_bounds, parameter_list, raw_parameter_refs, raw_parameter_values,
intervals

# Export measure datatypes
export AbstractMeasureData, DiscreteMeasureData, MultiDiscreteMeasureData,
Measure, MeasureData, MeasureRef

# Export measure methods
export add_measure, measure, measure_function, measure_data, expand,
expand_all_measures!, expect, support_sum, measure_name, measure_data_in_hold_bounds,
make_point_variable_ref, make_reduced_variable_ref, expand_measure, integral,
set_integral_defaults, integral_defaults, coefficients, weight_function,
default_weight, add_measure_variable, delete_reduced_variable,
delete_internal_reduced_variable, expand_measures, reduction_info

# Export constraint datatypes
export BoundedScalarConstraint, ConstraintData, InfOptConstraintRef,
InfiniteConstraintRef, FiniteConstraintRef

# Export transcription datatypes
export TranscriptionData, TranscriptionModel

# Export transcription methods
export is_transcription_model, transcription_data, transcription_variable,
transcription_constraint, transcription_model

# Export optimize methods
export optimizer_model, set_optimizer_model, optimizer_model_ready,
set_optimizer_model_ready, build_optimizer_model!, optimizer_model_key,
optimizer_model_variable, optimizer_model_constraint, variable_supports,
constraint_supports, constraint_parameter_refs, add_infinite_model_optimizer,
clear_optimizer_model_build!

# Export result query methods
export map_value, map_optimizer_index, map_dual, map_shadow_price,
map_lp_rhs_perturbation_range, map_lp_objective_perturbation_range

# Export printing methods
export bound_string

# Export support generation and fill-in methods
export fill_in_supports!, generate_supports

end # module
