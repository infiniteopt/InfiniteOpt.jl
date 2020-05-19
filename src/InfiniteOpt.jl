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

# Import methods
include("datatypes.jl")
include("infinite_sets.jl")
include("general_variables.jl")
include("scalar_parameters.jl")
include("array_parameters.jl")
include("variable_basics.jl")
include("infinite_variables.jl")
include("reduced_variables.jl")
include("point_variables.jl")
include("hold_variables.jl")
include("expressions.jl")
include("measures.jl")
include("constraints.jl")
include("macros.jl")
include("objective.jl")
include("measure_expansions.jl")
include("optimize.jl")
# include("results.jl")
include("show.jl")
include("utilities.jl")

# Import TranscriptionOpt
include("TranscriptionOpt/TranscriptionOpt.jl")
using .TranscriptionOpt

# Import MeasureToolbox
include("MeasureToolbox/MeasureToolbox.jl")
using .MeasureToolbox: Automatic, UniTrapezoid, UniMCSampling,
UniIndepMCSampling, Quadrature, GaussHermite, GaussLegendre, GaussLaguerre,
MultiMCSampling, MultiIndepMCSampling, generate_integral_data,
uni_integral_defaults, set_uni_integral_defaults, integral,
multi_integral_defaults, set_multi_integral_defaults, expect, support_sum,
@integral, @expect, @support_sum

# Export the measure toolbox functions and constants
export Automatic, UniTrapezoid, UniMCSampling, UniIndepMCSampling, Quadrature,
GaussHermite, GaussLegendre, GaussLaguerre, MultiMCSampling,
MultiIndepMCSampling, generate_integral_data, uni_integral_defaults,
set_uni_integral_defaults, integral, multi_integral_defaults,
set_multi_integral_defaults, expect, support_sum, @integral, @expect,
@support_sum

# Export core datatypes
export AbstractDataObject, InfiniteModel

# Export macros
export @independent_parameter, @dependent_parameters, @infinite_parameter,
@finite_parameter, @infinite_variable, @point_variable, @hold_variable,
@infinite_parameter, @BDconstraint, @set_parameter_bounds, @add_parameter_bounds,
@measure

# Export constants
export Infinite, Hold, Point, sampling, quadrature

# Export index types
export AbstractInfOptIndex, ObjectIndex, IndependentParameterIndex,
DependentParametersIndex, DependentParameterIndex, FiniteParameterIndex,
InfiniteVariableIndex, PointVariableIndex, ReducedVariableIndex,
HoldVariableIndex, MeasureIndex, ConstraintIndex, InfiniteParameterIndex,
FiniteVariableIndex

# Export infinite set types
export AbstractInfiniteSet, InfiniteScalarSet, InfiniteArraySet, IntervalSet,
UniDistributionSet, MultiDistributionSet, CollectionSet, UserDefined,
MCSample, UniformGrid, Mixture, All, WeightedSample

# Export infinite set methods
export collection_sets, supports_in_set, generate_support_values

# Export parameter datatypes
export InfOptParameter, ScalarParameter, IndependentParameter, FiniteParameter,
DependentParameters, ScalarParameterData, MultiParameterData,
IndependentParameterRef, DependentParameterRef, FiniteParameterRef,
ScalarParameterRef, InfiniteParameter

# Export parameter methods
export build_parameter, add_parameter, add_parameters,
infinite_set, set_infinite_set, parameter_by_name, num_parameters,
all_parameters, num_supports, has_supports,
set_supports, add_supports, delete_supports, supports, used_by_constraint,
used_by_measure, used_by_infinite_variable, is_used, generate_and_add_supports!,
significant_digits

# Export variable datatypes
export InfOptVariable, InfiniteVariable, ReducedVariable,
PointVariable, HoldVariable, ParameterBounds, VariableData, GeneralVariableRef,
DispatchVariableRef, InfiniteVariableRef, ReducedVariableRef,
MeasureFiniteVariableRef, FiniteVariableRef, PointVariableRef, HoldVariableRef,
DecisionVariableRef, UserDecisionVariableRef

# Export variable methods
export used_by_objective, infinite_variable_ref, parameter_refs,
set_parameter_refs, add_parameter_ref, used_by_point_variable, parameter_values,
eval_supports, used_by_reduced_variable, has_parameter_bounds, parameter_bounds,
set_parameter_bounds, add_parameter_bounds, delete_parameter_bounds,
parameter_list, raw_parameter_refs, raw_parameter_values,
intervals, dispatch_variable_ref, internal_reduced_variable

# Export measure datatypes
export AbstractMeasureData, DiscreteMeasureData, FunctionalDiscreteMeasureData,
Measure, MeasureData, MeasureRef

# Export measure methods
export default_weight, support_label, coefficient_function, coefficients,
weight_function, measure_data_in_hold_bounds, build_measure, min_num_supports,
add_supports_to_parameters, add_measure, measure_function, measure_data,
is_analytic, measure, num_measures, all_measures

# Export objective methods
export objective_has_measures

# Export constraint datatypes
export BoundedScalarConstraint, ConstraintData, InfOptConstraintRef,
InfiniteConstraintRef, FiniteConstraintRef

# Export constraint methods
export original_parameter_bounds

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
export bound_string, set_string, measure_data_string

# Export support generation and fill-in methods
export fill_in_supports!, generate_supports

end # module
