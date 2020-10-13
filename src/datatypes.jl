################################################################################
#                                INDEX TYPES
################################################################################
"""
    AbstractInfOptIndex

An abstract type for all index objects used in `InfiniteOpt`.
"""
abstract type AbstractInfOptIndex end

"""
    ObjectIndex <: AbstractInfOptIndex

An abstract type for indices of objects stored in `MOI.Utilities.CleverDicts`.
"""
abstract type ObjectIndex <: AbstractInfOptIndex end

"""
    IndependentParameterIndex <: ObjectIndex

A `DataType` for storing the index of a [`IndependentParameter`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct IndependentParameterIndex <: ObjectIndex
    value::Int64
end

"""
    DependentParametersIndex <: ObjectIndex

A `DataType` for storing the index of a [`DependentParameters`](@ref) object.

**Fields**
- `value::Int64`: The index value.
"""
struct DependentParametersIndex <: ObjectIndex
    value::Int64
end

"""
    DependentParameterIndex <: AbstractInfOptIndex

A `DataType` for storing the index of an indiviudal parameter in a
[`DependentParameters`](@ref) object.

**Fields**
- `object_index::DependentParametersIndex`: The index of the parameter collection.
- `param_index::Int`: The index of the individual parameter in the above object.
"""
struct DependentParameterIndex <: AbstractInfOptIndex
    object_index::DependentParametersIndex
    param_index::Int
end

# Define convenient alias for infinite parameters
const InfiniteParameterIndex = Union{IndependentParameterIndex, DependentParameterIndex}

"""
    FiniteParameterIndex <: ObjectIndex

A `DataType` for storing the index of a [`FiniteParameter`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct FiniteParameterIndex <: ObjectIndex
    value::Int64
end

"""
    InfiniteVariableIndex <: ObjectIndex

A `DataType` for storing the index of a [`InfiniteVariable`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct InfiniteVariableIndex <: ObjectIndex
    value::Int64
end

"""
    ReducedVariableIndex <: ObjectIndex

A `DataType` for storing the index of a [`ReducedVariable`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct ReducedVariableIndex <: ObjectIndex
    value::Int64
end

"""
    PointVariableIndex <: ObjectIndex

A `DataType` for storing the index of a [`PointVariable`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct PointVariableIndex <: ObjectIndex
    value::Int64
end

"""
    HoldVariableIndex <: ObjectIndex

A `DataType` for storing the index of a [`HoldVariable`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct HoldVariableIndex <: ObjectIndex
    value::Int64
end

"""
    DerivativeIndex <: ObjectIndex

A `DataType` for storing the index of a [`Derivative`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct DerivativeIndex <: ObjectIndex
    value::Int64
end

# Define convenient aliases
const FiniteVariableIndex = Union{PointVariableIndex, HoldVariableIndex,
                                  FiniteParameterIndex}

"""
    MeasureIndex <: ObjectIndex

A `DataType` for storing the index of a [`Measure`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct MeasureIndex <: ObjectIndex
    value::Int64
end

"""
    ConstraintIndex <: ObjectIndex

A `DataType` for storing the index of constraint objects.

**Fields**
- `value::Int64`: The index value.
"""
struct ConstraintIndex <: ObjectIndex
    value::Int64
end

## Extend the CleverDicts key access methods
# index_to_key
function MOIUC.index_to_key(::Type{C},
                            index::Int64) where {C <: ObjectIndex}
    return C(index)
end

# key_to_index
function MOIUC.key_to_index(key::ObjectIndex)::Int64
    return key.value
end

# Extend Base functions
Base.length(index::AbstractInfOptIndex) = 1
Base.broadcastable(index::AbstractInfOptIndex) = Ref(index)
Base.iterate(index::AbstractInfOptIndex) = (index, true)
Base.iterate(index::AbstractInfOptIndex, state) = nothing
Base.hash(v::ObjectIndex, h::UInt) = hash(v.value, h)

################################################################################
#                            INFINITE SET TYPES
################################################################################
"""
    AbstractInfiniteSet

An abstract type for sets that characterize infinite parameters.
"""
abstract type AbstractInfiniteSet end

"""
    InfiniteScalarSet <: AbstractInfiniteSet

An abstract type for infinite sets that are one-dimensional.
"""
abstract type InfiniteScalarSet <: AbstractInfiniteSet end

"""
    IntervalSet <: InfiniteScalarSet

A `DataType` that stores the lower and upper interval bounds for infinite
parameters that are continuous over a certain that interval. This is for use
with a [`IndependentParameter`](@ref).

**Fields**
- `lower_bound::Float64` Lower bound of the infinite parameter.
- `upper_bound::Float64` Upper bound of the infinite parameter.
"""
struct IntervalSet <: InfiniteScalarSet
    lower_bound::Float64
    upper_bound::Float64
    function IntervalSet(lower::Real, upper::Real)
        if lower > upper
            error("Invalid interval set bounds, lower bound is greater than " *
                  "upper bound.")
        end
        return new(lower, upper)
    end
end

"""
    UniDistributionSet{T <: Distributions.UnivariateDistribution} <: InfiniteScalarSet

A `DataType` that stores the distribution characterizing an infinite parameter that
is random. This is for use with a [`IndependentParameter`](@ref).

**Fields**
- `distribution::T` Distribution of the random parameter.
"""
struct UniDistributionSet{T <: Distributions.UnivariateDistribution} <: InfiniteScalarSet
    distribution::T
end

"""
    InfiniteArraySet <: AbstractInfiniteSet

An abstract type for multi-dimensional infinite sets.
"""
abstract type InfiniteArraySet <: AbstractInfiniteSet end

# Make convenient Union for below
const NonUnivariateDistribution = Union{Distributions.MultivariateDistribution,
                                        Distributions.MatrixDistribution}

"""
    MultiDistributionSet{T <: NonUnivariateDistribution} <: InfiniteArraySet

A `DataType` that stores the distribution characterizing a collection of
infinite parameters that follows its form. This is for use with
[`DependentParameters`](@ref).

**Fields**
- `distribution::T` Distribution of the random parameters.
"""
struct MultiDistributionSet{T <: NonUnivariateDistribution} <: InfiniteArraySet
    distribution::T
end

# make convenient alias for distribution sets
const DistributionSet = Union{UniDistributionSet, MultiDistributionSet}

"""
    CollectionSet{T <: InfiniteScalarSet} <: InfiniteArraySet

A `DataType` that stores a collection of `InfiniteScalarSet`s characterizing a
collection of infinite parameters that follows its form. This is for use with
[`DependentParameters`](@ref).

**Fields**
- `sets::Array{T}` The collection of scalar sets.
"""
struct CollectionSet{T <: InfiniteScalarSet} <: InfiniteArraySet
    sets::Vector{T}
end

################################################################################
#                         DERIVATIVE EVALUATION METHODS
################################################################################
"""
    AbstractDerivativeMethod

An abstract type for storing derivative evaluation data that is pertinent to its 
reformation/transcription. 
"""
abstract type AbstractDerivativeMethod end 

"""
    GenerativeDerivativeMethod <: AbstractDerivativeMethod

An abstract type for derivative evaluation method types that will require support 
generation when employed (e.g., internal node points associated with orthogonal 
collocation). Such methods can be used with derivatives that on independent 
infinite parameters, but cannot be used for ones that depend on dependent parameters.
"""
abstract type GenerativeDerivativeMethod <: AbstractDerivativeMethod end 

"""
    OCTechnique

An abstract type for the method used to carry out orthogonal collocation.
"""
abstract type OCTechnique end

"""
    Lobatto <: OCTechnique

A quadrature method label for orthogonal collocation method that generates
internal nodes between public supports using Lobatto quadrature method.
"""
struct Lobatto <: OCTechnique end

"""
    OrthogonalCollocation <: GenerativeDerivativeMethod 

A `DataType` for storing information about orthogonal collocation method
for derivative evaluation. Note that the constructor for this method is of the
form: 
```julia 
    OrthogonalCollocation(num_nodes::Int, [technique::Type{<:OCTechnique} = Labatto])
```
where `num_nodes` is total number of nodes for each collocation interval. In 
practice, this corresponds to `num_nodes = num_internal_nodes + 2`. 

**Fields**
- `num_internal_nodes::Int`: The number of internal collocation points (nodes) 
  between the each support pair.
- `technique::Type{<:OCTechnique}`: The method used to produce the points.
"""
struct OrthogonalCollocation <: GenerativeDerivativeMethod 
    num_internal_nodes::Int
    technique::DataType
    # make the constructor 
    function OrthogonalCollocation(num_nodes::Int, 
        technique::Type{<:OCTechnique} = Lobatto
        )::OrthogonalCollocation
        num_nodes >= 2 || error("Must specify at least 2 collocation points (i.e., " *
                                "the bounds of each support interval with no internal " * 
                                "support nodes).")
        return new(num_nodes - 2, technique)
    end
end

"""
    NonGenerativeDerivativeMethod <: AbstractDerivativeMethod

An abstract type for derivative evaluation method types that do not require the 
definition of additional support points. Such methods are amendable to any 
derivative in InfiniteOpt including those with dependent infinite parameter 
dependencies.
"""
abstract type NonGenerativeDerivativeMethod <: AbstractDerivativeMethod end

"""
    FDTechnique

An abstract data type for labels of specific techniques applied in the finite 
difference method in derivative evaluation.
"""
abstract type FDTechnique end

"""
    Forward <: FDTechnique

A technique label for finite difference method that implements a forward 
difference approximation.
"""
struct Forward <: FDTechnique end

"""
    Central <: FDTechnique

A technique label for finite difference method that implements a central 
difference approximation.
"""
struct Central <: FDTechnique end

"""
    Backward <: FDTechnique

A technique label for finite difference method that implements a backward 
difference approximation.
"""
struct Backward <: FDTechnique end

"""
    FiniteDifference <: NonGenerativeDerivativeMethod

A `DataType` for information about finite difference method applied to 
a derivative evaluation. Note that the constructor is of the form:
```julia 
    FiniteDifference([technique::Type{<:FDTechnique} = Backward],
                     [add_boundary_constr::Bool = true])
```
where `technique` is the indicated finite difference method to be applied and 
`add_boundary_constr` indicates if the finite difference equation corresponding to 
a boundary support should be included. Thus, for backward difference since
corresponds to the terminal point and for forward difference this corresponds to 
the initial point. We recommend using `add_boundary_constr = false` when an final 
condition is given with a backward method or when an initial condition is given 
with a forward method. Note that this argument is ignored for central finite 
difference which cannot include any boundary points.

**Fields** 
- `technique::Type{<:FDTechnique}`: Mathematical technqiue behind finite difference
- `add_boundary_constraint::Bool`: Indicate if the boundary constraint should be 
  included in the transcription (e.g., the terminal boundary backward equation for 
  backward difference)
"""
struct FiniteDifference <: NonGenerativeDerivativeMethod 
    technique::DataType
    add_boundary_constraint::Bool
    # set the constructor 
    function FiniteDifference(technique::Type{<:FDTechnique} = Backward, 
                              add_boundary_constr::Bool = true)
        return new(technique, add_boundary_constr)
    end
end

"""
    support_label(method::GenerativeDerivativeMethod)

Return the support label associated with `method` if there is one, errors otherwise. 
This should be extended for any `GenerativeDerivativeMethod`.
"""
function support_label(method::AbstractDerivativeMethod)::DataType
    error("`support_label` not defined for derivative methods of type `$(typeof(method))`.")
end

# Extend support_label for OrthogonalCollocation
function support_label(method::OrthogonalCollocation)::DataType
    return OrthogonalCollocationNode
end

################################################################################
#                              PARAMETER TYPES
################################################################################
"""
    InfOptParameter <: JuMP.AbstractVariable

An abstract type for all parameters used in InfiniteOpt.
"""
abstract type InfOptParameter <: JuMP.AbstractVariable end

"""
    ScalarParameter <: InfOptParameter

An abstract type for scalar parameters used in InfiniteOpt.
"""
abstract type ScalarParameter <: InfOptParameter end

"""
    IndependentParameter{T <: InfiniteScalarSet,
                         M <: AbstractDerivativeMethod} <: ScalarParameter

A `DataType` for storing independent scalar infinite parameters.

**Fields**
- `set::T`: The infinite set that characterizes the parameter.
- `supports::DataStructures.SortedDict{Float64, Set{DataType}}`: The support points
   used to discretize the parameter and their associated type labels stored as
   `DataTypes`s which should be a subtype of [`AbstractSupportLabel`](@ref).
- `sig_digits::Int`: The number of significant digits used to round the support values.
- `derivative_method::M`: The derivative evaluation method used for derivatives that
   are conducted with respect to this parameter.
"""
struct IndependentParameter{T <: InfiniteScalarSet, 
                            M <: AbstractDerivativeMethod} <: ScalarParameter
    set::T
    supports::DataStructures.SortedDict{Float64, Set{DataType}} # Support to label set
    sig_digits::Int
    derivative_method::M
end

"""
    FiniteParameter <: ScalarParameter

A `DataType` for storing finite parameters meant to be nested in expressions
and replaced with their values at runtime.

**Fields**
- `value::Float64`: The parameter value.
"""
struct FiniteParameter <: ScalarParameter
    value::Float64
end

"""
    DependentParameters{T <: InfiniteArraySet, 
                        M <: NonGenerativeDerivativeMethod} <: InfOptParameter

A `DataType` for storing a collection of dependent infinite parameters.

**Fields**
- `set::T`: The infinite set that characterizes the parameters.
- `supports::Dict{Vector{Float64}, Set{DataType}}`: Support dictionary where keys
              are supports and the values are the set of labels for each support.
- `sig_digits::Int`: The number of significant digits used to round the support values.
- `derivative_methods::Vector{M}`: The derivative evaluation methods associated with 
  each parameter.
"""
struct DependentParameters{T <: InfiniteArraySet, 
                           M <: NonGenerativeDerivativeMethod} <: InfOptParameter
    set::T
    supports::Dict{Vector{Float64}, Set{DataType}} # Support to label set
    sig_digits::Int
    derivative_methods::Vector{M}
end

# Define convenient alias for infinite types
const InfiniteParameter = Union{IndependentParameter, DependentParameters}

"""
    AbstractDataObject

An abstract type for `DataType`s that store core variable `DataType`s and their
model specific information (e.g., dependency mappings). These are what are
stored in the `InfiniteModel` `CleverDict`s.
"""
abstract type AbstractDataObject end

"""
    ScalarParameterData{P <: ScalarParameter} <: AbstractDataObject

A mutable `DataType` for storing `ScalarParameter`s and their data.

**Fields**
- `parameter::P`: The scalar parameter.
- `object_num::Int`: The location of the corresponding `ObjectIndex` in
    `InfiniteModel.param_object_indices` (given by `InfiniteModel.last_object_num`).
- `parameter_num::Int`: Given by `InfiniteModel.last_param_num` (updated when
                        prior parameters are deleted)
- `name::String`: The name used for printing.
- `infinite_var_indices::Vector{InfiniteVariableIndex}`: Indices of dependent
   infinite variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{ConstraintIndex}`: Indices of dependent constraints.
- `in_objective::Bool`: Is this used in objective? This should be true only for finite parameters.
- `has_internal_supports::Bool`: Does this parameter have internal supports?
- `has_derivative_supports::Bool`: Have any derivative specfic supports been added?
- `has_deriv_constrs::Bool`: Have any derivative evaluation constraints been added 
                             to the infinite model associated with this parameter?
"""
mutable struct ScalarParameterData{P <: ScalarParameter} <: AbstractDataObject
    parameter::P
    object_num::Int
    parameter_num::Int
    name::String
    infinite_var_indices::Vector{InfiniteVariableIndex}
    derivative_indices::Vector{DerivativeIndex}
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{ConstraintIndex}
    in_objective::Bool
    has_internal_supports::Bool
    has_derivative_supports::Bool
    has_deriv_constrs::Bool
end

# Convenient constructor
function ScalarParameterData(param::P,
                             object_num::Int,
                             parameter_num::Int,
                             name::String = ""
                             ) where {P <: ScalarParameter}
    return ScalarParameterData{P}(param, object_num, parameter_num, name,
                                  InfiniteVariableIndex[], DerivativeIndex[], 
                                  MeasureIndex[], ConstraintIndex[], false, false, 
                                  false, false)
end

"""
    MultiParameterData{P <: DependentParameters} <: AbstractDataObject

A mutable `DataType` for storing [`DependentParameters`](@ref) and their data.

**Fields**
- `parameters::P`: The parameter collection.
- `object_num::Int`: The location of the corresponding `ObjectIndex` in
   `InfiniteModel.param_object_indices` (given by `InfiniteModel.last_object_num`).
- `parameter_nums::UnitRange{Int}`: Given by `InfiniteModel.last_param_num`
                                    (updated when prior parameters are deleted)
- `names::Vector{String}`: The names used for printing each parameter.
- `infinite_var_indices::Vector{InfiniteVariableIndex}`: Indices of
   dependent infinite variables.
- `derivative_indices::Vector{Vector{DerivativeIndex}} `: Indices of dependent derivatives.
- `measure_indices::Vector{Vector{MeasureIndex}}`: Indices of dependent measures.
- `constraint_indices::Vector{Vector{ConstraintIndex}}`: Indices of dependent
  constraints.
- `has_internal_supports::Bool`: Does this parameter have internal supports?
- `has_deriv_constrs::Bool`: Have any derivative evaluation constraints been added 
                             to the infinite model associated with this parameter?
"""
mutable struct MultiParameterData{P <: DependentParameters} <: AbstractDataObject
    parameters::P
    object_num::Int
    parameter_nums::UnitRange{Int}
    names::Vector{String}
    infinite_var_indices::Vector{InfiniteVariableIndex}
    derivative_indices::Vector{Vector{DerivativeIndex}} 
    measure_indices::Vector{Vector{MeasureIndex}}
    constraint_indices::Vector{Vector{ConstraintIndex}}
    has_internal_supports::Bool
    has_deriv_constrs::Bool
end

# Convenient constructor 
function MultiParameterData(params::P,
                            object_num::Int,
                            parameter_nums::UnitRange{Int},
                            names::Vector{String},
                            ) where {P <: DependentParameters}
    return MultiParameterData{P}(params, object_num, parameter_nums, names,
                                 InfiniteVariableIndex[],
                                 [DerivativeIndex[] for i in eachindex(names)],
                                 [MeasureIndex[] for i in eachindex(names)],
                                 [ConstraintIndex[] for i in eachindex(names)],
                                 false, false)
end

################################################################################
#                             PARAMETER BOUNDS
################################################################################
"""
    ParameterBounds{P <: GeneralVariableRef}

A `DataType` for storing intervaled bounds of parameters. This is used to define
subdomains of [`HoldVariable`](@ref)s and [`BoundedScalarConstraint`](@ref)s.
Note that the GeneralVariableRef must pertain to infinite parameters.

**Fields**
- `intervals::Dict{GeneralVariableRef, IntervalSet}`: A dictionary
  of interval bounds on infinite parameters.
"""
struct ParameterBounds{P <: JuMP.AbstractVariableRef}
    intervals::Dict{P, IntervalSet}
end

################################################################################
#                               VARIABLE TYPES
################################################################################
"""
    InfOptVariable <: JuMP.AbstractVariable

An abstract type for infinite, reduced, point, and hold variables.
"""
abstract type InfOptVariable <: JuMP.AbstractVariable end

"""
    InfiniteVariable{P <: GeneralVariableRef} <: InfOptVariable

A `DataType` for storing core infinite variable information. Note that indices
that refer to the same dependent parameter group must be in the same tuple element.
It is important to note that `info.start` should contain a start value function
that generates the start value for a given infinite parameter support. This
function should map a support to a start value using user-formatting if
`is_vector_start = false`, otherwise it should do the mapping using a single
support vector as input. Also, the variable reference type `P` must pertain to
infinite parameters.

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, Function}`: JuMP variable information.
- `parameter_refs::VectorTuple{P}`: The infinite parameter references that
                                    parameterize the variable.
- `parameter_nums::Vector{Int}`: The parameter numbers of `parameter_refs`.
- `object_nums::Vector{Int}`: The parameter object numbers associated with `parameter_refs`.
- `is_vector_start::Bool`: Does the start function take support values formatted as vectors?
"""
struct InfiniteVariable{P <: JuMP.AbstractVariableRef} <: InfOptVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, Function}
    parameter_refs::VectorTuple{P}
    parameter_nums::Vector{Int}
    object_nums::Vector{Int}
    is_vector_start::Bool
end

"""
    ReducedVariable{I <: GeneralVariableRef} <: InfOptVariable

A `DataType` for storing reduced infinite variables which partially support an
infinite variable.

**Fields**
- `infinite_variable_ref::I`: The original infinite/derivvative variable.
- `eval_supports::Dict{Int, Float64}`: The original parameter tuple linear indices
                                     to the evaluation supports.
- `parameter_nums::Vector{Int}`: The parameter numbers associated with the reduced
                                 `parameter_refs`.
- `object_nums::Vector{Int}`: The parameter object numbers associated with the
                              reduced `parameter_refs`.
"""
struct ReducedVariable{I <: JuMP.AbstractVariableRef} <: InfOptVariable
    infinite_variable_ref::I
    eval_supports::Dict{Int, Float64}
    parameter_nums::Vector{Int}
    object_nums::Vector{Int}
end

"""
    PointVariable{I <: GeneralVariableRef} <: InfOptVariable

A `DataType` for storing point variable information. Note that the elements
`parameter_values` field must match the format of the parameter reference tuple
defined in [`InfiniteVariable`](@ref)

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}` JuMP Variable information.
- `infinite_variable_ref::I` The infinite variable/derivative reference
    associated with the point variable.
- `parameter_values::Vector{Float64}` The infinite parameter values
    defining the point.
"""
struct PointVariable{I <: JuMP.AbstractVariableRef} <: InfOptVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    infinite_variable_ref::I
    parameter_values::Vector{Float64}
end

"""
    HoldVariable{P <: GeneralVariableRef} <: InfOptVariable

A `DataType` for storing hold variable information.

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}` JuMP variable information.
- `parameter_bounds::ParameterBounds{P}` Valid parameter sub-domains
"""
struct HoldVariable{P <: JuMP.AbstractVariableRef} <: InfOptVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    parameter_bounds::ParameterBounds{P}
end

"""
    VariableData{V <: InfOptVariable} <: AbstractDataObject

A mutable `DataType` for storing `InfOptVariable`s and their data.

**Fields**
- `variable::V`: The scalar variable.
- `name::String`: The name used for printing.
- `lower_bound_index::Union{ConstraintIndex, Nothing}`: Index of lower bound constraint.
- `upper_bound_index::Union{ConstraintIndex, Nothing}`: Index of upper bound constraint.
- `fix_index::Union{ConstraintIndex, Nothing}`: Index on fixing constraint.
- `zero_one_index::Union{ConstraintIndex, Nothing}`: Index of binary constraint.
- `integrality_index::Union{ConstraintIndex, Nothing}`: Index of integer constraint.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{ConstraintIndex}`: Indices of dependent constraints.
- `in_objective::Bool`: Is this used in objective?
- `point_var_indices::Vector{PointVariableIndex}`: Indices of dependent point variables.
- `reduced_var_indices::Vector{ReducedVariableIndex}`: Indices of dependent reduced variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
"""
mutable struct VariableData{V <: InfOptVariable} <: AbstractDataObject
    variable::V
    name::String
    lower_bound_index::Union{ConstraintIndex, Nothing}
    upper_bound_index::Union{ConstraintIndex, Nothing}
    fix_index::Union{ConstraintIndex, Nothing}
    zero_one_index::Union{ConstraintIndex, Nothing}
    integrality_index::Union{ConstraintIndex, Nothing}
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{ConstraintIndex}
    in_objective::Bool
    point_var_indices::Vector{PointVariableIndex} # InfiniteVariables only
    reduced_var_indices::Vector{ReducedVariableIndex} # InfiniteVariables only
    derivative_indices::Vector{DerivativeIndex} # infinite and reduced only
    function VariableData(var::V, name::String = "") where {V <: InfOptVariable}
        return new{V}(var, name, nothing, nothing, nothing, nothing, nothing,
                   MeasureIndex[], ConstraintIndex[], false, PointVariableIndex[],
                   ReducedVariableIndex[], DerivativeIndex[])
    end
end

################################################################################
#                              DERIVATIVE TYPES
################################################################################
"""
    Derivative{V <: GeneralVariableRef} <: InfOptVariable

A `DataType` for storing core infinite derivative information. This follows a 
derivative of the form: ``\\frac{\\partial x(\\alpha, \\hdots)}{\\partial \\alpha}`` 
where ``x(\\alpha, \\hdots)`` is an infinite variable and ``\\alpha`` is an infinite 
parameter. Here, both ``x`` and ``\\alpha`` must be scalars. 

It is important to note that `info.start` should contain a start value function
that generates the start value for a given infinite parameter support. This
function should map a support to a start value using user-formatting if
`is_vector_start = false`, otherwise it should do the mapping using a single
support vector as input. Also, the variable reference type `V` must pertain to
infinite variables and parameters.

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, Function}`: JuMP variable information.
- `is_vector_start::Bool`: Does the start function take support values formatted as vectors?
- `variable_ref::V`: The variable reference of the infinite variable argument.
- `parameter_ref::V`: The variable reference of the infinite parameter the defines the
   differential operator.
"""
struct Derivative{V <: JuMP.AbstractVariableRef} <: InfOptVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, Function}
    is_vector_start::Bool
    variable_ref::V # could be ref of infinite/reduced variable/derivative or measure (top of derivative)
    parameter_ref::V # a scalar infinite parameter ref (bottom of derivative)
end

################################################################################
#                               MEASURE TYPES
################################################################################
"""
    AbstractMeasureData

An abstract type to define data for measures to define the behavior of
[`Measure`](@ref).
"""
abstract type AbstractMeasureData end

"""
    DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                        Vector{<:JuMP.AbstractVariableRef}},
                        N, B <: Union{Float64, Vector{Float64}}
                        } <: AbstractMeasureData

A DataType for immutable measure abstraction data where the
abstraction is of the form:
``measure = \\int_{\\tau \\in T} f(\\tau) w(\\tau) d\\tau \\approx \\sum_{i = 1}^N \\alpha_i f(\\tau_i) w(\\tau_i)``.
The supports and coefficients are immutable (i.e., they will not change
even if supports are changed for the underlying infinite parameter.) This
type can be used for both 1-dimensional and multi-dimensional measures.

**Fields**
- `parameter_refs::P`: The infinite parameter(s) over which the integration occurs.
                       These can be comprised of multiple independent parameters,
                       but dependent parameters cannot be mixed with other types.
- `coefficients::Vector{Float64}`: Coefficients ``\\alpha_i`` for the above
                                   measure abstraction.
- `supports::Array{Float64, N}`: Supports points ``\\tau_i``. This is a `Vector`
                                 if only one parameter is given, otherwise it is
                                 a `Matrix` where the supports are stored column-wise.
- `label::DataType`: Label for the support points ``\\tau_i`` when stored in the
                   infinite parameter(s), stemming from [`AbstractSupportLabel`](@ref).
- `weight_function::Function`: Weighting function ``w`` must map an individual
                               support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bound in accordance with ``T``, this denotes the
                    intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bound.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                           Vector{<:JuMP.AbstractVariableRef}},
                           N, B <: Union{Float64, Vector{Float64}}
                           } <: AbstractMeasureData
    parameter_refs::P
    coefficients::Vector{Float64}
    supports::Array{Float64, N} # supports are stored column-wise
    label::DataType # label that will used when the supports are added to the model
    weight_function::Function # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function DiscreteMeasureData(param_ref::V, coeffs::Vector{<:Real},
                                 supps::Vector{<:Real}, label::DataType,
                                 weight_func::Function,
                                 lower_bound::Real,
                                 upper_bound::Real,
                                 expect::Bool
                                 ) where {V <: JuMP.AbstractVariableRef}
        return new{V, 1, Float64}(param_ref, coeffs, supps, label, weight_func,
                                  lower_bound, upper_bound, expect)
    end
    # multi constructor
    function DiscreteMeasureData(param_refs::Vector{V}, coeffs::Vector{<:Real},
                                 supps::Matrix{<:Real}, label::DataType,
                                 weight_func::Function,
                                 lower_bound::Vector{<:Real},
                                 upper_bound::Vector{<:Real},
                                 expect::Bool
                                 ) where {V <: JuMP.AbstractVariableRef}
        return new{Vector{V}, 2, Vector{Float64}}(param_refs, coeffs, supps,
                                                  label, weight_func, lower_bound,
                                                  upper_bound, expect)
    end
end

"""
    FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                  Vector{<:JuMP.AbstractVariableRef}},
                                  B <: Union{Float64, Vector{Float64}}
                                  } <: AbstractMeasureData

A DataType for mutable measure abstraction data where the
abstraction is of the form:
``measure = \\int_{\\tau \\in T} f(\\tau) w(\\tau) d\\tau \\approx \\sum_{i = 1}^N \\alpha_i f(\\tau_i) w(\\tau_i)``.
This abstraction is equivalent to that of [`DiscreteMeasureData`](@ref), but
the difference is that the supports are not fully known at the time of measure
creation. Thus, functions are stored that will be used to generate the
concrete support points ``\\tau_i`` and their coefficients ``\\alpha_i`` when
the measure is evaluated (expanded). These supports are identified/generated
in accordance with the `label` with a gaurantee that at least `num_supports` are
generated. For example, if `label = MCSample` and `num_supports = 100` then
the measure will use all of the supports stored in the `parameter_refs` with the
label `MCSample` and will ensure there are at least 100 are generated. This
type can be used for both 1-dimensional and multi-dimensional measures.

**Fields**
- `parameter_refs::P`: The infinite parameter(s) over which the integration occurs.
                     These can be comprised of multiple independent parameters,
                     but dependent parameters cannot be mixed with other types.
- `coeff_function::Function`: Coefficient generation function making ``\\alpha_i``
                              for the above measure abstraction. It should take
                              all the supports as input (formatted as an Array)
                              and return the corresponding vector of coefficients.
- `min_num_supports::Int`: Specifies the minimum number of supports ``\\tau_i``
                       desired in association with `parameter_refs` and `label`.
- `label::DataType`: Label for the support points ``\\tau_i`` which are/will be
                   stored in the infinite parameter(s), stemming from [`AbstractSupportLabel`](@ref).
- `weight_function::Function`: Weighting function ``w`` must map an individual
                              support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bounds in accordance with ``T``, this denotes the
                  intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bounds.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                     Vector{<:JuMP.AbstractVariableRef}},
                                     B <: Union{Float64, Vector{Float64}}
                                     } <: AbstractMeasureData
    parameter_refs::P
    coeff_function::Function # supports --> coefficient vector
    min_num_supports::Int # minimum number of supports
    label::DataType # support label of included supports
    weight_function::Function # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function FunctionalDiscreteMeasureData(param_ref::V, coeff_func::Function,
                                           num_supps::Int, label::DataType,
                                           weight_func::Function,
                                           lower_bound::Real,
                                           upper_bound::Real,
                                           expect::Bool
                                           ) where {V <: JuMP.AbstractVariableRef}
        return new{V, Float64}(param_ref, coeff_func, num_supps, label, weight_func,
                               lower_bound, upper_bound, expect)
    end
    # multi constructor
    function FunctionalDiscreteMeasureData(param_refs::Vector{V},
                                           coeff_func::Function,
                                           num_supps::Int, label::DataType,
                                           weight_func::Function,
                                           lower_bound::Vector{<:Real},
                                           upper_bound::Vector{<:Real},
                                           expect::Bool
                                           ) where {V <: JuMP.AbstractVariableRef}
        return new{Vector{V}, Vector{Float64}}(param_refs, coeff_func, num_supps,
                                               label, weight_func, lower_bound,
                                               upper_bound, expect)
    end
end

"""
    Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}

A `DataType` for measure abstractions. The abstraction is determined by `data`
and is enacted on `func` when the measure is evaluated (expended).

**Fields**
- `func::T` The `InfiniteOpt` expression to be measured.
- `data::V` Data of the abstraction as described in a `AbstractMeasureData`
            concrete subtype.
- `object_nums::Vector{Int}`: The parameter object numbers of the evaluated
                              measure expression (i.e., the object numbers of
                              `func` excluding those that belong to `data`).
- `parameter_nums::Vector{Int}`: The parameter numbers that parameterize the
                                 evaluated measure expression. (i.e., the
                                 parameter numbers of `func` excluding those
                                 that belong to `data`).
- `constant_func::Bool`: Indicates if `func` is not parameterized by the infinite
                         parameters in `data`. (i.e., do the object numbers of
                         `func` and `data` have no intersection?) This is useful
                         to enable analytic evaluations if possible.
"""
struct Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}
    func::T
    data::V
    object_nums::Vector{Int}
    parameter_nums::Vector{Int}
    constant_func::Bool
end

"""
    MeasureData <: AbstractDataObject

A mutable `DataType` for storing [`Measure`](@ref)s and their data.

**Fields**
- `measure::Measure`: The measure structure.
- `name::String`: The base name used for printing `name(meas_expr d(par))`.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{ConstraintIndex}`: Indices of dependent constraints.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `in_objective::Bool`: Is this used in objective?
"""
mutable struct MeasureData <: AbstractDataObject
    measure::Measure
    name::String
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{ConstraintIndex}
    derivative_indices::Vector{DerivativeIndex}
    in_objective::Bool
    function MeasureData(measure::Measure, name::String = "measure")
        return new(measure, name, MeasureIndex[], ConstraintIndex[], 
                   DerivativeIndex[], false)
    end
end

################################################################################
#                              CONSTRAINT TYPES
################################################################################
"""
    BoundedScalarConstraint{F <: JuMP.AbstractJuMPScalar,
                            S <: MOI.AbstractScalarSet,
                            P <: GeneralVariableRef
                            } <: JuMP.AbstractConstraint

A `DataType` that stores scalar constraints that are defined over a sub-domain
of infinite parameters.

**Fields**
- `func::F` The JuMP object.
- `set::S` The MOI set.
- `bounds::ParameterBounds{P}` Set of valid parameter
    sub-domains that further boundconstraint.
- `orig_bounds::ParameterBounds{P}` Set of the constraint's
    original parameter sub-domains (not considering hold variables)
"""
struct BoundedScalarConstraint{F <: JuMP.AbstractJuMPScalar,
                               S <: MOI.AbstractScalarSet,
                               P <: JuMP.AbstractVariableRef
                               } <: JuMP.AbstractConstraint
    func::F
    set::S
    bounds::ParameterBounds{P}
    orig_bounds::ParameterBounds{P}
end

"""
    ConstraintData <: AbstractDataObject

A mutable `DataType` for storing constraints and their data.

**Fields**
- `constraint::JuMP.AbstractConstraint`: The scalar constraint.
- `object_nums::Vector{Int}`: The object numbers of the parameter objects that the
                              constraint depends on.
- `name::String`: The name used for printing.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `is_info_constraint::Bool`: Is this is constraint based on variable info (e.g., lower bound)
"""
mutable struct ConstraintData <: AbstractDataObject
    constraint::JuMP.AbstractConstraint
    object_nums::Vector{Int}
    name::String
    measure_indices::Vector{MeasureIndex}
    is_info_constraint::Bool
end

################################################################################
#                                INFINITE MODEL
################################################################################
"""
    InfiniteModel <: JuMP.AbstractModel

A `DataType` for storing all of the mathematical modeling information needed to
model an optmization problem with an infinite-dimensional decision space.

**Fields**
- `independent_params::MOIUC.CleverDict{IndependentParameterIndex, ScalarParameterData{IndependentParameter}}`:
   The independent parameters and their mapping information.
- `dependent_params::MOIUC.CleverDict{DependentParametersIndex, MultiParameterData}`:
   The dependent parameters and their mapping information.
- `finite_params::MOIUC.CleverDict{FiniteParameterIndex, ScalarParameterData{FiniteParameter}}`:
   The finite parameters and their mapping information.
- `name_to_param::Union{Dict{String, AbstractInfOptIndex}, Nothing}`:
   Field to help find a parameter given the name.
- `last_param_num::Int`: The last parameter number to be used.
- `param_object_indices::Vector{Union{IndependentParameterIndex, DependentParametersIndex}}`:
  The collection of parameter object indices in creation order.
- `infinite_vars::MOIUC.CleverDict{InfiniteVariableIndex, <:VariableData{<:InfiniteVariable}}`:
   The infinite variables and their mapping information.
- `reduced_vars::MOIUC.CleverDict{ReducedVariableIndex, <:VariableData{<:ReducedVariable}}`:
   The reduced infinite variables and their mapping information.
- `point_vars::MOIUC.CleverDict{PointVariableIndex, <:VariableData{<:PointVariable}}`:
   The point variables and their mapping information.
- `hold_vars::MOIUC.CleverDict{HoldVariableIndex, <:VariableData{<:HoldVariable}}`:
   The hold variables and their mapping information.
- `name_to_var::Union{Dict{String, AbstractInfOptIndex}, Nothing}`:
   Field to help find a variable given the name.
- `has_hold_bounds::Bool`:
   Does any variable have parameter bounds?
- `derivatives::MOIUC.CleverDict{DerivativeIndex, <:VariableData{<:Derivative}}`:
  The derivatives and their mapping information.
- `deriv_lookup::Dict{<:Tuple, DerivativeIndex}`: Map derivative variable-parameter 
  pairs to a derivative index to prevent duplicates.
- `measures::MOIUC.CleverDict{MeasureIndex, MeasureData}`:
   The measures and their mapping information.
- `integral_defaults::Dict{Symbol}`:
   The default keyword arguments for [`integral`](@ref).
- `constraints::MOIUC.CleverDict{ConstraintIndex, ConstraintData}`:
   The constraints and their mapping information.
- `name_to_constr::Union{Dict{String, ConstraintIndex}, Nothing}`:
   Field to help find a constraint given the name.
- `objective_sense::MOI.OptimizationSense`: Objective sense.
- `objective_function::JuMP.AbstractJuMPScalar`: Finite scalar function.
- `objective_has_measures::Bool`: Does the objective contain measures?
- `obj_dict::Dict{Symbol, Any}`: Store Julia symbols used with `InfiniteModel`
- `optimizer_constructor`: MOI optimizer constructor (e.g., Gurobi.Optimizer).
- `optimizer_model::JuMP.Model`: Model used to solve `InfiniteModel`
- `ready_to_optimize::Bool`: Is the optimizer_model up to date.
- `ext::Dict{Symbol, Any}`: Store arbitrary extension information.
"""
mutable struct InfiniteModel <: JuMP.AbstractModel
    # Parameter Data
    independent_params::MOIUC.CleverDict{IndependentParameterIndex, ScalarParameterData{<:IndependentParameter}}
    dependent_params::MOIUC.CleverDict{DependentParametersIndex, MultiParameterData}
    finite_params::MOIUC.CleverDict{FiniteParameterIndex, ScalarParameterData{FiniteParameter}}
    name_to_param::Union{Dict{String, AbstractInfOptIndex}, Nothing}
    last_param_num::Int
    param_object_indices::Vector{Union{IndependentParameterIndex, DependentParametersIndex}}

    # Variable Data
    infinite_vars::MOIUC.CleverDict{InfiniteVariableIndex, <:VariableData{<:InfiniteVariable}}
    reduced_vars::MOIUC.CleverDict{ReducedVariableIndex, <:VariableData{<:ReducedVariable}}
    point_vars::MOIUC.CleverDict{PointVariableIndex, <:VariableData{<:PointVariable}}
    hold_vars::MOIUC.CleverDict{HoldVariableIndex, <:VariableData{<:HoldVariable}}
    name_to_var::Union{Dict{String, AbstractInfOptIndex}, Nothing}
    has_hold_bounds::Bool

    # Derivative Data 
    derivatives::MOIUC.CleverDict{DerivativeIndex, <:VariableData{<:Derivative}}
    deriv_lookup::Dict{<:Tuple, DerivativeIndex}

    # Measure Data
    measures::MOIUC.CleverDict{MeasureIndex, MeasureData}

    # Constraint Data
    constraints::MOIUC.CleverDict{ConstraintIndex, ConstraintData}
    name_to_constr::Union{Dict{String, ConstraintIndex}, Nothing}

    # Objective Data
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    objective_has_measures::Bool

    # Objects
    obj_dict::Dict{Symbol, Any}

    # Optimize Data
    optimizer_constructor::Any
    optimizer_model::JuMP.Model
    ready_to_optimize::Bool

    # Extensions
    ext::Dict{Symbol, Any}
end

"""
    InfiniteModel([optimizer_constructor];
                  [OptimizerModel::Function = TranscriptionModel,
                  caching_mode::MOIU.CachingOptimizerMode = MOIU.AUTOMATIC,
                  bridge_constraints::Bool = true, optimizer_model_kwargs...])

Return a new infinite model where an optimizer is specified if an
`optimizer_constructor` is given. The optimizer
can also later be set with the [`JuMP.set_optimizer`](@ref) call. By default
the `optimizer_model` data field is initialized with a
[`TranscriptionModel`](@ref), but a different type of model can be assigned via
[`set_optimizer_model`](@ref) as can be required by extensions.

**Example**
```jldoctest
julia> using InfiniteOpt, JuMP, Ipopt;

julia> model = InfiniteModel()
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Measures: 0
Derivatives: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> model = InfiniteModel(Ipopt.Optimizer)
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Measures: 0
Derivatives: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```
"""
function InfiniteModel(; OptimizerModel::Function = TranscriptionModel,
                       kwargs...)::InfiniteModel
    return InfiniteModel(# Parameters
                         MOIUC.CleverDict{IndependentParameterIndex, ScalarParameterData{<:IndependentParameter}}(),
                         MOIUC.CleverDict{DependentParametersIndex, MultiParameterData}(),
                         MOIUC.CleverDict{FiniteParameterIndex, ScalarParameterData{FiniteParameter}}(),
                         nothing, 0,
                         Union{IndependentParameterIndex, DependentParametersIndex}[],
                         # Variables
                         MOIUC.CleverDict{InfiniteVariableIndex, VariableData{InfiniteVariable{GeneralVariableRef}}}(),
                         MOIUC.CleverDict{ReducedVariableIndex, VariableData{ReducedVariable{GeneralVariableRef}}}(),
                         MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}(),
                         MOIUC.CleverDict{HoldVariableIndex, VariableData{HoldVariable{GeneralVariableRef}}}(),
                         nothing, false,
                         # Derivatives
                         MOIUC.CleverDict{DerivativeIndex, VariableData{Derivative{GeneralVariableRef}}}(),
                         Dict{Tuple{GeneralVariableRef, GeneralVariableRef}, DerivativeIndex}(),
                         # Measures
                         MOIUC.CleverDict{MeasureIndex, MeasureData}(),
                         # Constraints
                         MOIUC.CleverDict{ConstraintIndex, ConstraintData}(),
                         nothing,
                         # Objective
                         MOI.FEASIBILITY_SENSE,
                         zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}),
                         false,
                         # Object dictionary
                         Dict{Symbol, Any}(),
                         # Optimize data
                         nothing, OptimizerModel(;kwargs...), false,
                         # Extensions
                         Dict{Symbol, Any}()
                         )
end

## Set the optimizer_constructor depending on what it is
# MOI.OptimizerWithAttributes
function _set_optimizer_constructor(model::InfiniteModel,
                                    constructor::MOI.OptimizerWithAttributes)::Nothing
    model.optimizer_constructor = constructor.optimizer_constructor
    return
end

# No attributes
function _set_optimizer_constructor(model::InfiniteModel, constructor)::Nothing
    model.optimizer_constructor = constructor
    return
end

# Dispatch for InfiniteModel call with optimizer constructor
function InfiniteModel(optimizer_constructor;
                       OptimizerModel::Function = TranscriptionModel,
                       kwargs...)::InfiniteModel
    model = InfiniteModel()
    model.optimizer_model = OptimizerModel(optimizer_constructor; kwargs...)
    _set_optimizer_constructor(model, optimizer_constructor)
    return model
end

# Define basic InfiniteModel extensions
Base.broadcastable(model::InfiniteModel) = Ref(model)
JuMP.object_dictionary(model::InfiniteModel)::Dict{Symbol, Any} = model.obj_dict

# Define basic accessors
_last_param_num(model::InfiniteModel)::Int = model.last_param_num
_param_object_indices(model::InfiniteModel) = model.param_object_indices

################################################################################
#                             OBJECT REFERENCES
################################################################################
"""
    GeneralVariableRef <: JuMP.AbstractVariableRef

A `DataType` that serves as the principal variable reference in `InfiniteOpt`
for building variable expressions. It contains the needed information to
create a variable type specifc reference (e.g., [`InfiniteVariableRef`](@ref))
via [`dispatch_variable_ref`](@ref) to obtain the correct subtype of
[`DispatchVariableRef`](@ref) based off of `index_type`. This allows us to
construct expressions using concrete containers unlike previous versions of
`InfiniteOpt` which provides us a significant performance boost.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `raw_index::Int64`: The raw index to be used in the `index_type` constructor.
- `index_type::DataType`: The concrete [`AbstractInfOptIndex`](@ref) type/constructor.
- `param_index::Int`: The index of a parameter in [`DependentParameters`](@ref).
  This is ignored for other variable types.
"""
struct GeneralVariableRef <: JuMP.AbstractVariableRef
    model::InfiniteModel
    raw_index::Int64
    index_type::DataType
    param_index::Int # for DependentParameterRefs
    function GeneralVariableRef(model::InfiniteModel, raw_index,
                               index_type::DataType, param_index::Int = -1)
       return new(model, Int64(raw_index), index_type, param_index)
    end
end

"""
    DispatchVariableRef <: JuMP.AbstractVariableRef

An abstract type for variable references that are created from
[`GeneralVariableRef`](@ref)s and are used to dispatch to the appropriate
methods for that particular variable/parameter/measure type.
"""
abstract type DispatchVariableRef <: JuMP.AbstractVariableRef end

"""
    IndependentParameterRef <: DispatchVariableRef

A `DataType` for independent infinite parameters references that parameterize
infinite variables.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::IndependentParameterIndex`: Index of the parameter in model.
"""
struct IndependentParameterRef <: DispatchVariableRef
    model::InfiniteModel
    index::IndependentParameterIndex
end

"""
    DependentParameterRef <: DispatchVariableRef

A `DataType` for dependent infinite parameter references that parameterize
infinite variables.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::DependentParameterIndex`: Index of the dependent parameter.
"""
struct DependentParameterRef <: DispatchVariableRef
    model::InfiniteModel
    index::DependentParameterIndex
end

"""
    InfiniteVariableRef <: DispatchVariableRef

A `DataType` for untranscripted infinite dimensional variable references (e.g.,
second stage variables, time dependent variables).

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::InfiniteVariableIndex`: Index of the variable in model.
"""
struct InfiniteVariableRef <: DispatchVariableRef
    model::InfiniteModel
    index::InfiniteVariableIndex
end

"""
    ReducedVariableRef <: DispatchVariableRef

A `DataTyp`e for partially transcripted infinite dimensional variable references.
This is used to expand measures that contain infinite variables that are not
fully transcripted by the measure.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::ReducedVariableIndex`: Index of the variable in model.
"""
struct ReducedVariableRef <: DispatchVariableRef
    model::InfiniteModel
    index::ReducedVariableIndex
end

"""
    DerivativeRef <: DispatchVariableRef

A `DataType` for untranscripted derivative references.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::DerivativeIndex`: Index of the derivative in model.
"""
struct DerivativeRef <: DispatchVariableRef
    model::InfiniteModel
    index::DerivativeIndex
end

"""
    MeasureFiniteVariableRef <: DispatchVariableRef

An abstract type to define finite variable and measure references.
"""
abstract type MeasureFiniteVariableRef <: DispatchVariableRef end

"""
    MeasureRef <: FiniteVariableRef

A `DataType` for referring to measure abstractions.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::MeasureIndex`: Index of the measure in model.
"""
struct MeasureRef <: MeasureFiniteVariableRef
    model::InfiniteModel
    index::MeasureIndex
end

"""
    FiniteVariableRef <: MeasureFiniteVariableRef

An abstract type to define new finite variable references.
"""
abstract type FiniteVariableRef <: MeasureFiniteVariableRef end

"""
    PointVariableRef <: FiniteVariableRef

A `DataType` for variables defined at a transcipted point (e.g., second stage
variable at a particular scenario, dynamic variable at a discretized time point).

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::PointVariableIndex`: Index of the variable in model.
"""
struct PointVariableRef <: FiniteVariableRef
    model::InfiniteModel
    index::PointVariableIndex
end

"""
    HoldVariableRef <: FiniteVariableRef

A `DataType` for finite fixed variable references (e.g., first stage variables,
steady-state variables).

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::HoldVariableIndex`: Index of the variable in model.
"""
struct HoldVariableRef <: FiniteVariableRef
    model::InfiniteModel
    index::HoldVariableIndex
end

"""
    FiniteParameterRef <: FiniteVariableRef

A `DataType` for finite parameters references who are replaced with their values
at the transcription step.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::FiniteParameterIndex`: Index of the parameter in model.
"""
struct FiniteParameterRef <: FiniteVariableRef
    model::InfiniteModel
    index::FiniteParameterIndex
end

## Define convenient aliases
const DecisionVariableRef = Union{InfiniteVariableRef, ReducedVariableRef,
                                  PointVariableRef, HoldVariableRef, 
                                  DerivativeRef}

const UserDecisionVariableRef = Union{InfiniteVariableRef, PointVariableRef,
                                      HoldVariableRef, DerivativeRef}

const ScalarParameterRef = Union{IndependentParameterRef, FiniteParameterRef}

"""
    InfOptConstraintRef{S <: JuMP.AbstractShape}

A `DataType` for constraints that are in `InfiniteModel`s

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::ConstraintIndex`: Index of the constraint in model.
- `shape::JuMP.AbstractShape`: Shape of the constraint
"""
struct InfOptConstraintRef{S <: JuMP.AbstractShape}
    model::InfiniteModel
    index::ConstraintIndex
    shape::S
end

# Make dumby model type for calling @expression 
struct _DumbyModel <: JuMP.AbstractModel end
const _Model = _DumbyModel()

################################################################################
#                            PARAMETER BOUND METHODS
################################################################################
## Modify parameter dictionary to expand any multidimensional parameter keys
# Case where dictionary is already in correct form
function _expand_parameter_tuple(
    param_bounds::NTuple{N, Pair{GeneralVariableRef, IntervalSet}}
    )::Dict{GeneralVariableRef, IntervalSet} where {N}
    return Dict(param_bounds...)
end

# Case where dictionary contains vectors
function _expand_parameter_dict(
    param_bounds::NTuple{N, Pair{<:Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}, IntervalSet}}
    )::Dict{GeneralVariableRef, IntervalSet} where {N}
    # Initialize new dictionary
    new_dict = Dict{GeneralVariableRef, IntervalSet}()
    # Find vector keys and expand
    for (key, set) in param_bounds
        # expand over the array of parameters if this is
        if isa(key, AbstractArray)
            for param in key
                new_dict[param] = set
            end
        # otherwise we have parameter reference
        else
            new_dict[key] = set
        end
    end
    return new_dict
end

# Case where dictionary contains vectors
function _expand_parameter_dict(param_bounds::Tuple)
    error("Invalid parameter bound dictionary format.")
end

# Constructor for expanding array parameters
function ParameterBounds(intervals::NTuple{N, Pair}
    )::ParameterBounds{GeneralVariableRef} where {N}
    return ParameterBounds(_expand_parameter_dict(intervals))
end

# Default method
function ParameterBounds()::ParameterBounds{GeneralVariableRef}
    return ParameterBounds(Dict{GeneralVariableRef, IntervalSet}())
end

# Make dictionary accessor
function intervals(pb::ParameterBounds{P})::Dict{P, IntervalSet} where {P}
    return pb.intervals
end

# Extend simple 1 argument Base dispatches
for op = (:length, :isempty, :keys, :iterate)
    @eval Base.$op(bounds::ParameterBounds) = $op(intervals(bounds))
end

# Extend simple 2 argument Base dispatches where the second is arbitrary
for op = (:getindex, :haskey, :iterate)
    @eval Base.$op(bounds::ParameterBounds, arg) = $op(intervals(bounds), arg)
end

# Extend Base.:(==)
function Base.:(==)(bounds1::ParameterBounds, bounds2::ParameterBounds)::Bool
    return intervals(bounds1) == intervals(bounds2)
end

# Extend Base.copy
Base.copy(bounds::ParameterBounds) = ParameterBounds(copy(intervals(bounds)))

# Extend Base.setindex!
function Base.setindex!(pb::ParameterBounds{P}, value::IntervalSet,
                        index::P)::IntervalSet where {P}
    return intervals(pb)[index] = value
end

# Extend Base.delete!
function Base.delete!(pb::ParameterBounds{P}, key::P)::ParameterBounds{P} where {P}
    delete!(intervals(pb), key)
    return pb
end

# Extend Base.merge
function Base.merge(pb1::ParameterBounds{P},
                    pb2::ParameterBounds{P})::ParameterBounds{P} where {P}
    new_dict = merge(intervals(pb1), intervals(pb2))
    return ParameterBounds(new_dict)
end

# Extend Base.merge!
function Base.merge!(pb1::ParameterBounds{P},
                     pb2::ParameterBounds{P})::ParameterBounds{P} where {P}
    merge!(intervals(pb1), intervals(pb2))
    return pb1
end

# Extend Base.filter
function Base.filter(f::Function,
                     pb::ParameterBounds{P})::ParameterBounds{P} where {P}
    new_dict = filter(f, intervals(pb))
    return ParameterBounds(new_dict)
end
