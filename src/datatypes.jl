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
    ParameterFunctionIndex <: ObjectIndex

A `DataType` for storing the index of a [`ParameterFunction`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct ParameterFunctionIndex <: ObjectIndex
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
    SemiInfiniteVariableIndex <: ObjectIndex

A `DataType` for storing the index of a [`SemiInfiniteVariable`](@ref).

**Fields**
- `value::Int64`: The index value.
"""
struct SemiInfiniteVariableIndex <: ObjectIndex
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
    FiniteVariableIndex <: ObjectIndex

A `DataType` for storing the index of a `JuMP.ScalarVariable`.

**Fields**
- `value::Int64`: The index value.
"""
struct FiniteVariableIndex <: ObjectIndex
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
const FiniteIndex = Union{PointVariableIndex, FiniteVariableIndex,
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
InOptConstraintIndex <: ObjectIndex

A `DataType` for storing the index of a constraint.

**Fields**
- `value::Int64`: The index value.
"""
struct InfOptConstraintIndex <: ObjectIndex
    value::Int64
end

## Extend the CleverDicts key access methods
# index_to_key
function MOIUC.index_to_key(
    ::Type{C},
    index::Int64
    ) where {C <: ObjectIndex}
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
    AbstractInfiniteDomain

An abstract type for domains that characterize infinite parameters.
"""
abstract type AbstractInfiniteDomain end

"""
    InfiniteScalarDomain <: AbstractInfiniteDomain

An abstract type for infinite domains that are one-dimensional.
"""
abstract type InfiniteScalarDomain <: AbstractInfiniteDomain end

"""
    IntervalDomain <: InfiniteScalarDomain

A `DataType` that stores the lower and upper interval bounds for infinite
parameters that are continuous over a certain that interval. This is for use
with a [`IndependentParameter`](@ref).

**Fields**
- `lower_bound::Float64` Lower bound of the infinite parameter.
- `upper_bound::Float64` Upper bound of the infinite parameter.
"""
struct IntervalDomain <: InfiniteScalarDomain
    lower_bound::Float64
    upper_bound::Float64
    function IntervalDomain(lower::Real, upper::Real)
        if lower > upper
            error("Invalid interval domain bounds, lower bound is greater than " *
                  "upper bound.")
        end
        return new(lower, upper)
    end
end

"""
    UniDistributionDomain{T <: Distributions.UnivariateDistribution} <: InfiniteScalarDomain

A `DataType` that stores the distribution characterizing an infinite parameter that
is random. This is for use with a [`IndependentParameter`](@ref).

**Fields**
- `distribution::T` Distribution of the random parameter.
"""
struct UniDistributionDomain{T <: Distributions.UnivariateDistribution} <: InfiniteScalarDomain
    distribution::T
end

"""
    InfiniteArrayDomain <: AbstractInfiniteDomain

An abstract type for multi-dimensional infinite domains.
"""
abstract type InfiniteArrayDomain <: AbstractInfiniteDomain end

# Make convenient Union for below
const NonUnivariateDistribution = Union{Distributions.MultivariateDistribution,
                                        Distributions.MatrixDistribution}

"""
    MultiDistributionDomain{T <: NonUnivariateDistribution} <: InfiniteArrayDomain

A `DataType` that stores the distribution characterizing a collection of
infinite parameters that follows its form. This is for use with
[`DependentParameters`](@ref).

**Fields**
- `distribution::T` Distribution of the random parameters.
"""
struct MultiDistributionDomain{T <: NonUnivariateDistribution} <: InfiniteArrayDomain
    distribution::T
end

# Extend Base.:(==) for relevant cases
function Base.:(==)(d1::D, d2::D)::Bool where {D <: MultiDistributionDomain}
    return d1.distribution == d2.distribution
end

# make convenient alias for distribution domains
const DistributionDomain = Union{UniDistributionDomain, MultiDistributionDomain}

"""
    CollectionDomain{T <: InfiniteScalarDomain} <: InfiniteArrayDomain

A `DataType` that stores a collection of `InfiniteScalarDomain`s characterizing a
collection of infinite parameters that follows its form. This is for use with
[`DependentParameters`](@ref).

**Fields**
- `domains::Array{T}` The collection of scalar domains.
"""
struct CollectionDomain{T <: InfiniteScalarDomain} <: InfiniteArrayDomain
    domains::Vector{T}
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
    IndependentParameter{T <: InfiniteScalarDomain} <: ScalarParameter

A `DataType` for storing independent scalar infinite parameters.

**Fields**
- `domain::T`: The infinite domain that characterizes the parameter.
"""
struct IndependentParameter{T <: InfiniteScalarDomain} <: ScalarParameter
    domain::T
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
    DependentParameters{T <: InfiniteArrayDomain, 
                        M <: NonGenerativeDerivativeMethod} <: InfOptParameter

A `DataType` for storing a collection of dependent infinite parameters.

**Fields**
- `domain::T`: The infinite domain that characterizes the parameters.
"""
struct DependentParameters{T <: InfiniteArrayDomain} <: InfOptParameter
    domain::T
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
- `parameter_func_indices::Vector{ParameterFunctionIndex}`: Indices of dependent
   infinite parameter functions.
- `infinite_var_indices::Vector{InfiniteVariableIndex}`: Indices of dependent
   infinite variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `in_objective::Bool`: Is this used in objective? This should be true only for finite parameters.
- `has_deriv_constrs::Bool`: Have any derivative evaluation constraints been added 
                             to the infinite model associated with this parameter?
"""
mutable struct ScalarParameterData{P <: ScalarParameter} <: AbstractDataObject
    parameter::P
    object_num::Int
    parameter_num::Int
    name::String
    parameter_func_indices::Vector{ParameterFunctionIndex}
    infinite_var_indices::Vector{InfiniteVariableIndex}
    derivative_indices::Vector{DerivativeIndex}
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{InfOptConstraintIndex}
    in_objective::Bool
    has_deriv_constrs::Bool # TODO maybe remove this?
end

# Convenient constructor
function ScalarParameterData(param::P,
                             object_num::Int,
                             parameter_num::Int,
                             name::String = ""
                             ) where {P <: ScalarParameter}
    return ScalarParameterData{P}(param, object_num, parameter_num, name,
                                  ParameterFunctionIndex[], InfiniteVariableIndex[], 
                                  DerivativeIndex[], MeasureIndex[], 
                                  InfOptConstraintIndex[], false, false)
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
- `parameter_func_indices::Vector{ParameterFunctionIndex}`: Indices of
   dependent infinite parameter functions.
- `infinite_var_indices::Vector{InfiniteVariableIndex}`: Indices of
   dependent infinite variables.
- `derivative_indices::Vector{Vector{DerivativeIndex}} `: Indices of dependent derivatives.
- `measure_indices::Vector{Vector{MeasureIndex}}`: Indices of dependent measures.
- `constraint_indices::Vector{Vector{InfOptConstraintIndex}}`: Indices of dependent
  constraints.
- `has_deriv_constrs::Bool`: Have any derivative evaluation constraints been added 
                             to the infinite model associated with this parameter?
"""
mutable struct MultiParameterData{P <: DependentParameters} <: AbstractDataObject
    parameters::P
    object_num::Int
    parameter_nums::UnitRange{Int}
    names::Vector{String}
    parameter_func_indices::Vector{ParameterFunctionIndex}
    infinite_var_indices::Vector{InfiniteVariableIndex}
    derivative_indices::Vector{Vector{DerivativeIndex}} 
    measure_indices::Vector{Vector{MeasureIndex}}
    constraint_indices::Vector{Vector{InfOptConstraintIndex}}
    has_deriv_constrs::Vector{Bool} # TODO maybe remove?
end

# Convenient constructor 
function MultiParameterData(params::P,
                            object_num::Int,
                            parameter_nums::UnitRange{Int},
                            names::Vector{String},
                            ) where {P <: DependentParameters}
    return MultiParameterData{P}(params, object_num, parameter_nums, names,
                                 ParameterFunctionIndex[], InfiniteVariableIndex[],
                                 [DerivativeIndex[] for i in eachindex(names)],
                                 [MeasureIndex[] for i in eachindex(names)],
                                 [InfOptConstraintIndex[] for i in eachindex(names)],
                                 zeros(Bool, length(names)))
end

################################################################################
#                             DOMAIN RESTRICTIONS
################################################################################
"""
    DomainRestrictions{P <: GeneralVariableRef}

A `DataType` for storing interval domains that constrain particular infinite 
parameters to a subdomain relative to their full domain. This is used to define
subdomains of [`DomainRestrictedConstraint`](@ref)s.
Note that the GeneralVariableRef must pertain to infinite parameters.

The constructor syntax is
```julia 
DomainRestrictions(restrictions...)
```
where each argument of `restrictions` is one of the following forms:
- `pref => value`
- `pref => [lb, ub]`
- `pref => IntervalDomain(lb, ub)`
- `prefs => value`
- `prefs => [lb, ub]`
- `prefs => IntervalDomain(lb, ub)`.
Note that `pref` and `prefs` must correspond to infinite parameters. 

**Fields**
- `intervals::Dict{GeneralVariableRef, IntervalDomain}`: A dictionary
  of interval bounds on infinite parameters.
"""
struct DomainRestrictions{P <: JuMP.AbstractVariableRef}
    intervals::Dict{P, IntervalDomain}
end

################################################################################
#                        PARAMETER FUNCTION OBJECTS
################################################################################
"""
    ParameterFunction{F <: Function, VT <: VectorTuple}

A `DataType` for storing known functions of infinite parameters. These equate to arbitrary 
functions that take support instances of infinite parameters `parameter_refs` in 
as input and compute a scalar value as output via `func`. These can then can 
incorporated in expressions via [`ParameterFunctionRef`](@ref)s.

**Fields**
- `func::F`: The function the takes infinite parameters as input and provide a 
            scalar number as output.
- `parameter_refs::VT`: The infinite parameter references that serve as 
                                   inputs to `func`. Their formatting is analagous 
                                   to those of infinite variables. 
- `parameter_nums::Vector{Int}`: The parameter numbers of `parameter_refs`.
- `object_nums::Vector{Int}`: The parameter object numbers associated with `parameter_refs`.
"""
struct ParameterFunction{F <: Function, VT <: Collections.VectorTuple}
    func::F
    parameter_refs::VT
    object_nums::Vector{Int}
    parameter_nums::Vector{Int}
end

"""
    ParameterFunctionData{F <: ParameterFunction} <: AbstractDataObject

A mutable `DataType` for storing `ParameterFunction`s and their data.

**Fields**
- `func::F`: The parameter function.
- `name::String`: The name used for printing.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `semi_infinite_var_indices::Vector{SemiInfiniteVariableIndex}`: Indices of dependent semi-infinite variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
"""
mutable struct ParameterFunctionData{F <: ParameterFunction} <: AbstractDataObject
    func::F
    name::String
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{InfOptConstraintIndex}
    semi_infinite_var_indices::Vector{SemiInfiniteVariableIndex}
    derivative_indices::Vector{DerivativeIndex}
    function ParameterFunctionData(func::F, name::String = "") where {F <: ParameterFunction}
        return new{F}(func, name, MeasureIndex[], InfOptConstraintIndex[], 
                      SemiInfiniteVariableIndex[], DerivativeIndex[])
    end
end

################################################################################
#                               VARIABLE TYPES
################################################################################
"""
    InfiniteVariable{F <: Function, VT <: VectorTuple} <: JuMP.AbstractVariable

A `DataType` for storing core infinite variable information. Note that indices
that refer to the same dependent parameter group must be in the same tuple element.
It is important to note that `info.start` should contain a start value function
that generates the start value for a given infinite parameter support. This
function should map a support to a start value using user-formatting if
`is_vector_start = false`, otherwise it should do the mapping using a single
support vector as input.

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, F}`: JuMP variable information.
  Here the start value is a function that maps the parameter values to a start value.
- `parameter_refs::VT`: The infinite parameter references that parameterize the 
  variable.
- `parameter_nums::Vector{Int}`: The parameter numbers of `parameter_refs`.
- `object_nums::Vector{Int}`: The parameter object numbers associated with `parameter_refs`.
- `is_vector_start::Bool`: Does the start function take support values formatted as vectors?
"""
struct InfiniteVariable{F <: Function, VT <: Collections.VectorTuple} <: JuMP.AbstractVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, F}
    parameter_refs::VT
    parameter_nums::Vector{Int}
    object_nums::Vector{Int}
    is_vector_start::Bool
end

"""
    SemiInfiniteVariable{I <: GeneralVariableRef} <: JuMP.AbstractVariable

A `DataType` for storing semi-infinite variables which partially support an
infinite variable.

**Fields**
- `infinite_variable_ref::I`: The original infinite/derivvative variable.
- `eval_supports::Dict{Int, Float64}`: The original parameter tuple linear indices
                                     to the evaluation supports.
- `parameter_nums::Vector{Int}`: The parameter numbers associated with the evaluated
                                 `parameter_refs`.
- `object_nums::Vector{Int}`: The parameter object numbers associated with the
                              evaluated `parameter_refs`.
"""
struct SemiInfiniteVariable{I <: JuMP.AbstractVariableRef} <: JuMP.AbstractVariable
    infinite_variable_ref::I
    eval_supports::Dict{Int, Float64}
    parameter_nums::Vector{Int}
    object_nums::Vector{Int}
end

"""
    PointVariable{I <: GeneralVariableRef} <: JuMP.AbstractVariable

A `DataType` for storing point variable information. Note that the elements
`parameter_values` field must match the format of the parameter reference tuple
defined in [`InfiniteVariable`](@ref)

**Fields**
- `info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}`: JuMP Variable information.
- `infinite_variable_ref::I`: The infinite variable/derivative reference
    associated with the point variable.
- `parameter_values::Vector{Float64}`: The infinite parameter values
    defining the point.
"""
struct PointVariable{I <: JuMP.AbstractVariableRef} <: JuMP.AbstractVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    infinite_variable_ref::I
    parameter_values::Vector{Float64}
end

"""
    VariableData{V <: JuMP.AbstractVariable} <: AbstractDataObject

A mutable `DataType` for storing variables and their data.

**Fields**
- `variable::V`: The scalar variable.
- `name::String`: The name used for printing.
- `lower_bound_index::Union{InfOptConstraintIndex, Nothing}`: Index of lower bound constraint.
- `upper_bound_index::Union{InfOptConstraintIndex, Nothing}`: Index of upper bound constraint.
- `fix_index::Union{InfOptConstraintIndex, Nothing}`: Index on fixing constraint.
- `zero_one_index::Union{InfOptConstraintIndex, Nothing}`: Index of binary constraint.
- `integrality_index::Union{InfOptConstraintIndex, Nothing}`: Index of integer constraint.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `in_objective::Bool`: Is this used in objective?
- `point_var_indices::Vector{PointVariableIndex}`: Indices of dependent point variables.
- `semi_infinite_var_indices::Vector{SemiInfiniteVariableIndex}`: Indices of dependent semi-infinite variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `deriv_constr_indices::Vector{InfOptConstraintIndex}`: Indices of dependent derivative evaluation constraints.
"""
mutable struct VariableData{V <: JuMP.AbstractVariable} <: AbstractDataObject
    variable::V
    name::String
    lower_bound_index::Union{InfOptConstraintIndex, Nothing}
    upper_bound_index::Union{InfOptConstraintIndex, Nothing}
    fix_index::Union{InfOptConstraintIndex, Nothing}
    zero_one_index::Union{InfOptConstraintIndex, Nothing}
    integrality_index::Union{InfOptConstraintIndex, Nothing}
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{InfOptConstraintIndex}
    in_objective::Bool
    point_var_indices::Vector{PointVariableIndex} # InfiniteVariables only
    semi_infinite_var_indices::Vector{SemiInfiniteVariableIndex} # InfiniteVariables only
    derivative_indices::Vector{DerivativeIndex} # infinite and semi-infinite only
    deriv_constr_indices::Vector{InfOptConstraintIndex} # Derivatives only
end

# Define constructor
function VariableData(
    var::V, 
    name::String = ""
    )::VariableData{V} where {V <: JuMP.AbstractVariable}
    return VariableData{V}(var, name, nothing, nothing, nothing, nothing, nothing,
                           MeasureIndex[], InfOptConstraintIndex[], false, 
                           PointVariableIndex[], SemiInfiniteVariableIndex[], 
                           DerivativeIndex[], InfOptConstraintIndex[])
end

################################################################################
#                              DERIVATIVE TYPES
################################################################################
# TODO modify to store derivatives of arbitrary order
"""
    Derivative{F <: Function, V <: GeneralVariableRef} <: JuMP.AbstractVariable

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
- `info::JuMP.VariableInfo{Float64, Float64, Float64, F}`: JuMP variable information.
- `is_vector_start::Bool`: Does the start function take support values formatted as vectors?
- `variable_ref::V`: The variable reference of the infinite variable argument.
- `parameter_ref::V`: The variable reference of the infinite parameter the defines the
   differential operator.
"""
struct Derivative{F <: Function, V <: JuMP.AbstractVariableRef} <: JuMP.AbstractVariable
    info::JuMP.VariableInfo{Float64, Float64, Float64, F}
    is_vector_start::Bool
    variable_ref::V # could be ref of infinite/semi-infinite variable/derivative or measure (top of derivative)
    parameter_ref::V # a scalar infinite parameter ref (bottom of derivative)
end

################################################################################
#                               MEASURE TYPES
################################################################################
"""
    AbstractMeasureData

An abstract type for measures to define the behavior of a
[`Measure`](@ref).
"""
abstract type AbstractMeasureData end

"""
    IntegralData{P <: Union{JuMP.AbstractVariableRef, Vector{<:JuMP.AbstractVariableRef}},
                 D <: AbstractInfiniteDomain,
                 F <: Union{Nothing, Function}} <: AbstractMeasureData

A `DataType` for defining integral measures and storing their needed canonical
information; namely, the infinite parameter(s) that act as the independent 
variable, the domain of the integral, and a weighting function (if there is 
one).

**Fields**
- `parameter_refs::P`: The infinite parameter(s) that act as the independent 
variable(s).
- `domain::D`: The domain of the integral (must be a sub-domain of the parameter 
domain).
- `weight_func::F`: A function of the form `w(d)::Float64`, where `d` are the 
infinite parameters, that is multiplied against the integrant function.   
"""
struct IntegralData{P <: Union{JuMP.AbstractVariableRef,
                               Vector{<:JuMP.AbstractVariableRef}},
                    D <: AbstractInfiniteDomain,
                    F <: Union{Nothing, Function}} <: AbstractMeasureData
    parameter_refs::P
    domain::D # needs to be a subset or equal to the parameter domain(s)
    weight_func::F

    # Univariate constructor
    function IntegralData(
        pref::P, 
        domain::D, 
        func::F = nothing
        ) where {P <: JuMP.AbstractVariableRef, 
                 D <: InfiniteScalarDomain, 
                 F <: Union{Nothing, Function}}
        return new{P, D, F}(pref, domain, func)
    end

    # Multivariate constructor
    function IntegralData(
        prefs::P, 
        domain::D, 
        func::F = nothing
        ) where {P <: Vector{<:JuMP.AbstractVariableRef}, 
                 D <: InfiniteArrayDomain, 
                 F <: Union{Nothing, Function}}
        return new{P, D, F}(prefs, domain, func)
    end
end

"""
    ExpectationData{P <: Union{JuMP.AbstractVariableRef, Vector{<:JuMP.AbstractVariableRef}},
                    F <: Union{Nothing, Function}} <: AbstractMeasureData 

A `DataType` for storing expectation operators. Principally, this includes 
the infinite parameter(s) the expectation is with respect to and the underlying 
probability density function. 

**Fields**
- `parameter_refs::P`: The infinite parameter(s) the expectation is over.
- `pdf::F`: The pdf function used by the expectation.
"""
struct ExpectationData{P <: Union{JuMP.AbstractVariableRef,
                                  Vector{<:JuMP.AbstractVariableRef}},
                       F <: Union{Nothing, Function}} <: AbstractMeasureData 
    parameter_refs::P # these also define the domain implicitly
    pdf::F

    # Expectation constructor
    function ExpectationData(
        prefs::P,
        func::F = nothing
        ) where {P <: Union{JuMP.AbstractVariableRef, Vector{<:JuMP.AbstractVariableRef}}, 
                 F <: Union{Nothing, Function}}
        return new{P, F}(prefs, func)
    end
end

"""
    Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}

A `DataType` for measure objects. The type of measure is determined by `data`
and is enacted on `func` when the measure is evaluated (expanded).

**Fields**
- `func::T` The `InfiniteOpt` expression to be measured.
- `data::V` Data needed to describe the measure.
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
    MeasureData{M <: Measure} <: AbstractDataObject

A mutable `DataType` for storing [`Measure`](@ref)s and their data.

**Fields**
- `measure::M`: The measure structure.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `in_objective::Bool`: Is this used in objective?
"""
mutable struct MeasureData{M <: Measure} <: AbstractDataObject
    measure::M
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{InfOptConstraintIndex}
    derivative_indices::Vector{DerivativeIndex}
    in_objective::Bool
end

function MeasureData(measure::M) where {M <: Measure}
    return MeasureData{M}(measure, MeasureIndex[], InfOptConstraintIndex[], 
                          DerivativeIndex[], false)
end

################################################################################
#                              CONSTRAINT TYPES
################################################################################
"""
    DomainRestrictedConstraint{C <: JuMP.AbstractConstraint, 
                               P <: GeneralVariableRef
                               } <: JuMP.AbstractConstraint

A `DataType` for creating a constraint with enforced `DomainRestrictions`. For 
example this may pertain to a boundary condition.

**Fields**
- `constraint::C`: The optimization constraint.
- `restrictions::DomainRestrictions{P}`: The restrictions that determine the 
   sub-domain of the constraint.
"""
struct DomainRestrictedConstraint{C <: JuMP.AbstractConstraint, 
                                  P <: JuMP.AbstractVariableRef
                                  } <: JuMP.AbstractConstraint
    constraint::C
    restrictions::DomainRestrictions{P}
end

"""
    ConstraintData{C <: JuMP.AbstractConstraint} <: AbstractDataObject

A mutable `DataType` for storing constraints and their data.

**Fields**
- `constraint::C`: The constraint.
- `object_nums::Vector{Int}`: The object numbers of the parameter objects that the
                              constraint depends on.
- `name::String`: The name used for printing.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `is_info_constraint::Bool`: Is this is constraint based on variable info 
   (e.g., lower bound)
"""
mutable struct ConstraintData{C <: JuMP.AbstractConstraint} <: AbstractDataObject
    constraint::C
    object_nums::Vector{Int}
    name::String
    measure_indices::Vector{MeasureIndex}
    is_info_constraint::Bool
end

################################################################################
#                                TRANSFORM API
################################################################################
"""
    AbstractTransformAttr

An abstract type for attributes used by a transformation backend that are stored 
in an `InfiniteModel`'s cache as the model is created. 
"""
abstract type AbstractTransformAttr end

"""
    FiniteParameterAttr <: AbstractTransformAttr

A finite parameter attribute that is used by a transformation backend. This 
is intended to be used by transformation backends that require/allow additional 
information about finite parameters to be specified by the user.
"""
abstract type FiniteParameterAttr <: AbstractTransformAttr end

"""
    InfiniteParameterAttr <: AbstractTransformAttr

An infinite parameter attribute that is used by a transformation backend. This 
is intended to be used by transformation backends that require/allow additional 
information about infinite parameters to be specified by the user.
"""
abstract type InfiniteParameterAttr <: AbstractTransformAttr end

"""
    VariableAttr <: AbstractTransformAttr

A variable attribute that is used by a transformation backend. This is 
intended to be used by transformation backends that require/allow additional 
information about variables to be specified by the user.
"""
abstract type VariableAttr <: AbstractTransformAttr end

"""
    DerivativeAttr <: AbstractTransformAttr

A derivative attribute that is used by a transformation backend. This is 
intended to be used by transformation backends that require/allow additional 
information about derivatives to be specified by the user.
"""
abstract type DerivativeAttr <: AbstractTransformAttr end

"""
    MeasureAttr <: AbstractTransformAttr

A measure attribute that is used by a transformation backend. This is 
intended to be used by transformation backends that require/allow additional 
information about measure to be specified by the user.
"""
abstract type MeasureAttr <: AbstractTransformAttr end

"""
    ConstraintAttr <: AbstractTransformAttr

A constraint attribute that is used by a transformation backend. This is 
intended to be used by transformation backends that require/allow additional 
information about constraints to be specified by the user.
"""
abstract type ConstraintAttr <: AbstractTransformAttr end

"""
    ModelAttr <: AbstractTransformAttr

A model attribute that is used by a transformation backend. This is 
intended to be used by transformation backends that require/allow additional 
information about the model to be specified by the user.
"""
abstract type ModelAttr <: AbstractTransformAttr end

"""
    TransformAttrCache

A convenient container for storing all the transformation attributes stored in 
an `InfiniteModel` that can be used by the transformation backend.
"""
struct TransformAttrCache
    finite_parameters::Dict{Tuple{FiniteParameterIndex, FiniteParameterAttr}, Any}
    independent_parameters::Dict{Tuple{IndependentParameterIndex, InfiniteParameterAttr}, Any}
    dependent_params::Dict{Tuple{DependentParametersIndex, InfiniteParameterAttr}, Any} # TODO maybe allow for parameter-wsie attributes
    infinite_variables::Dict{Tuple{InfiniteVariableIndex, VariableAttr}, Any}
    semi_ifninite_variables::Dict{Tuple{SemiInfiniteVariableIndex, VariableAttr}, Any}
    point_variables::Dict{Tuple{PointVariableIndex, VariableAttr}, Any}
    derivatives::Dict{Tuple{DerivativeIndex, DerivativeAttr}, Any}
    measures::Dict{Tuple{MeasureIndex, MeasureAttr}, Any}
    constraints::Dict{Tuple{InfOptConstraintIndex, ConstraintAttr}, Any}
    model::Dict{ModelAttr, Any}

    # Constructor 
    function TransformAttrCache()
        return new(
            Dict{Tuple{FiniteParameterIndex, FiniteParameterAttr}, Any}(),
            Dict{Tuple{IndependentParameterIndex, InfiniteParameterAttr}, Any}(),
            Dict{Tuple{DependentParametersIndex, InfiniteParameterAttr}, Any}(),
            Dict{Tuple{InfiniteVariableIndex, VariableAttr}, Any}(),
            Dict{Tuple{SemiInfiniteVariableIndex, VariableAttr}, Any}(),
            Dict{Tuple{PointVariableIndex, VariableAttr}, Any}(),
            Dict{Tuple{DerivativeIndex, DerivativeAttr}, Any}(),
            Dict{Tuple{MeasureIndex, MeasureAttr}, Any}(),
            Dict{Tuple{InfOptConstraintIndex, ConstraintAttr}, Any}(),
            Dict{ModelAttr, Any}()
            )
    end
end

"""
    ObjectWithAttributes{O, D <: Dict}

This serves as a wrapper type to store a modeling object (e.g., an infinite 
parameter) and [`AbstractTransformAttr`](@ref)s that should be added when the 
object is added to the model.

**Fields**
- `object:O`: The modeling object to be added.
- `attributes::D`: The dictionary of attributes to be added (`attr` => `value`). 
"""
struct ObjectWithAttributes{O, D <: Dict}
    object::O
    attributes::D
end

"""
    AbstractTransformBackend

An abstract type for transformation interfaces/models that act as a backend for 
`InfiniteModel`s.
"""
abstract type AbstractTransformBackend end

# TODO maybe add more types if needed

################################################################################
#                                INFINITE MODEL
################################################################################
const DefaultSigDigits = 12

"""
    InfiniteModel <: JuMP.AbstractModel

A `DataType` for storing all of the mathematical modeling information needed to
model an optmization problem with an infinite-dimensional decision space.
"""
mutable struct InfiniteModel <: JuMP.AbstractModel
    # Parameter Data
    independent_params::MOIUC.CleverDict{IndependentParameterIndex, ScalarParameterData{<:IndependentParameter}}
    dependent_params::MOIUC.CleverDict{DependentParametersIndex, MultiParameterData}
    finite_params::MOIUC.CleverDict{FiniteParameterIndex, ScalarParameterData{FiniteParameter}}
    name_to_param::Union{Dict{String, AbstractInfOptIndex}, Nothing}
    last_param_num::Int
    param_object_indices::Vector{Union{IndependentParameterIndex, DependentParametersIndex}}
    param_functions::MOIUC.CleverDict{ParameterFunctionIndex, <:ParameterFunctionData}

    # Variable Data
    infinite_vars::MOIUC.CleverDict{InfiniteVariableIndex, <:VariableData{<:InfiniteVariable}}
    semi_infinite_vars::MOIUC.CleverDict{SemiInfiniteVariableIndex, <:VariableData{<:SemiInfiniteVariable}}
    semi_lookup::Dict{<:Tuple, SemiInfiniteVariableIndex}
    point_vars::MOIUC.CleverDict{PointVariableIndex, <:VariableData{<:PointVariable}}
    point_lookup::Dict{<:Tuple, PointVariableIndex}
    finite_vars::MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}
    name_to_var::Union{Dict{String, AbstractInfOptIndex}, Nothing}
    significant_digits::Int

    # Derivative Data 
    derivatives::MOIUC.CleverDict{DerivativeIndex, <:VariableData{<:Derivative}}
    deriv_lookup::Dict{<:Tuple, DerivativeIndex}

    # Measure Data
    measures::MOIUC.CleverDict{MeasureIndex, <:MeasureData}

    # Constraint Data
    constraints::MOIUC.CleverDict{InfOptConstraintIndex, <:ConstraintData}
    constraint_restrictions::Dict{InfOptConstraintIndex, <:DomainRestrictions}
    name_to_constr::Union{Dict{String, InfOptConstraintIndex}, Nothing}

    # Objective Data
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    objective_has_measures::Bool # TODO maybe remove this?

    # Function Registration
    registrations::Vector{Any}
    func_lookup::Dict{Tuple{Symbol, Int}, Function}

    # Objects
    obj_dict::Dict{Symbol, Any}

    # Transformation/Optimize Data
    transform_attrs::TransformAttrCache
    transform_backend::Union{AbstractTransformBackend, Nothing}
    backend_ready::Bool

    # Extensions
    ext::Dict{Symbol, Any}
end

"""
    InfiniteModel([backend::AbstractTransformBackend]; [sig_digits::Int = DefaultSigDigits, kwargs..])

The core modeling object used to store infinite-dimensional optimization formulations. 
Optionally, a transformation backend `backend` can be specified that will altimately 
be used to transform and solve a model. Can specify the number of significant digits 
that should be used to process supports (e.g., variable points). Certain backends 
and/or extension packages might allow the specification of keyword arguments to 
provide [`ModelAttr`](@ref)s. 

**Example**
```julia-repl
julia> model = InfiniteModel()

TODO finish
```
"""
function InfiniteModel(; sig_digits::Int = DefaultSigDigits)
    return InfiniteModel(
        # Parameters
        MOIUC.CleverDict{IndependentParameterIndex, ScalarParameterData{<:IndependentParameter}}(),
        MOIUC.CleverDict{DependentParametersIndex, MultiParameterData}(),
        MOIUC.CleverDict{FiniteParameterIndex, ScalarParameterData{FiniteParameter}}(),
        nothing, 0,
        Union{IndependentParameterIndex, DependentParametersIndex}[],
        MOIUC.CleverDict{ParameterFunctionIndex, ParameterFunctionData{<:ParameterFunction}}(),
        # Variables
        MOIUC.CleverDict{InfiniteVariableIndex, VariableData{<:InfiniteVariable}}(),
        MOIUC.CleverDict{SemiInfiniteVariableIndex, VariableData{SemiInfiniteVariable{GeneralVariableRef}}}(),
        Dict{Tuple{GeneralVariableRef, Dict{Int, Float64}}, SemiInfiniteVariableIndex}(),
        MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}(),
        Dict{Tuple{GeneralVariableRef, Vector{Float64}}, PointVariableIndex}(),
        MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}(),
        nothing, sig_digits,
        # Derivatives
        MOIUC.CleverDict{DerivativeIndex, VariableData{<:Derivative}}(),
        Dict{Tuple{GeneralVariableRef, GeneralVariableRef}, DerivativeIndex}(),
        # Measures
        MOIUC.CleverDict{MeasureIndex, MeasureData{<:Measure}}(),
        # Constraints
        MOIUC.CleverDict{InfOptConstraintIndex, ConstraintData{<:JuMP.AbstractConstraint}}(),
        Dict{InfOptConstraintIndex, DomainRestrictions{GeneralVariableRef}}(),
        nothing,
        # Objective
        MOI.FEASIBILITY_SENSE,
        zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}),
        false,
        # registration
        RegisteredFunction[],
        Dict{Tuple{Symbol, Int}, Function}(),
        # Object dictionary
        Dict{Symbol, Any}(),
        # Transform data
        TransformAttrCache(), nothing, false,
        # Extensions
        Dict{Symbol, Any}()
    )
end

# TODO make constructor with transform method


# Define basic InfiniteModel extension functions
Base.broadcastable(model::InfiniteModel) = Ref(model)

"""
    JuMP.object_dictionary(model::InfiniteModel)::Dict{Symbol, Any}

Return the dictionary that maps the symbol name of a macro defined object (e.g., 
a parameter, variable, or constraint) to the corresponding object. Objects are 
registered to a specific symbol in the macros. For example, 
`@variable(model, x[1:2, 1:2])` registers the array of variables
`x` to the symbol `:x`.
"""
JuMP.object_dictionary(model::InfiniteModel) = model.obj_dict

"""
    Base.empty!(model::InfiniteModel)::InfiniteModel

Clear out `model` of everything. 
"""
function Base.empty!(model::InfiniteModel)
    # parameters
    empty!(model.independent_params)
    empty!(model.dependent_params)
    empty!(model.finite_params)
    model.name_to_param = nothing
    model.last_param_num = 0
    empty!(model.param_object_indices)
    empty!(model.param_functions)
    # variables
    empty!(model.infinite_vars)
    empty!(model.semi_infinite_vars)
    empty!(model.semi_lookup)
    empty!(model.point_vars)
    empty!(model.point_lookup)
    empty!(model.finite_vars)
    model.name_to_var = nothing
    model.significant_digits = DefaultSigDigits
    # derivatives and measures
    empty!(model.derivatives)
    empty!(model.deriv_lookup)
    empty!(model.measures)
    # constraints
    empty!(model.constraints)
    empty!(model.constraint_restrictions)
    model.name_to_constr = nothing
    # objective
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective_function = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    model.objective_has_measures = false
    # other stuff
    empty!(model.registrations)
    empty!(model.obj_dict)
    empty!(model.func_lookup)
    model.transform_attrs = TransformAttrCache()
    model.transform_backend = nothing
    model.backend_ready = false
    empty!(model.ext)
    return model
end

# Define basic accessors
_last_param_num(model::InfiniteModel) = model.last_param_num
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
    ParameterFunctionRef <: DispatchVariableRef

A `DataType` for infinite parameter function references.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::ParameterFunctionIndex`: Index of the infinite parameter function.
"""
struct ParameterFunctionRef <: DispatchVariableRef
    model::InfiniteModel
    index::ParameterFunctionIndex
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
    SemiInfiniteVariableRef <: DispatchVariableRef

A `DataType` for partially transcripted infinite dimensional variable references.
This is used to expand measures that contain infinite variables that are not
fully transcripted by the measure.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::SemiInfiniteVariableIndex`: Index of the variable in model.
"""
struct SemiInfiniteVariableRef <: DispatchVariableRef
    model::InfiniteModel
    index::SemiInfiniteVariableIndex
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
    MeasureRef <: DispatchVariableRef

A `DataType` for referring to measure abstractions.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::MeasureIndex`: Index of the measure in model.
"""
struct MeasureRef <: DispatchVariableRef
    model::InfiniteModel
    index::MeasureIndex
end

"""
    FiniteRef <: DispatchVariableRef

An abstract type for variable references that are finite.
"""
abstract type FiniteRef <: DispatchVariableRef end

"""
    PointVariableRef <: FiniteRef

A `DataType` for variables defined at a transcipted point (e.g., second stage
variable at a particular scenario, dynamic variable at a discretized time point).

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::PointVariableIndex`: Index of the variable in model.
"""
struct PointVariableRef <: FiniteRef
    model::InfiniteModel
    index::PointVariableIndex
end

"""
    FiniteVariableRef <: FiniteRef

A `DataType` for finite fixed variable references (e.g., first stage variables,
steady-state variables).

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::FiniteVariableIndex`: Index of the variable in model.
"""
struct FiniteVariableRef <: FiniteRef
    model::InfiniteModel
    index::FiniteVariableIndex
end

"""
    FiniteParameterRef <: FiniteRef

A `DataType` for finite parameters references who are replaced with their values
at the transcription step.

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::FiniteParameterIndex`: Index of the parameter in model.
"""
struct FiniteParameterRef <: FiniteRef
    model::InfiniteModel
    index::FiniteParameterIndex
end

## Define convenient aliases
const DecisionVariableRef = Union{InfiniteVariableRef, SemiInfiniteVariableRef,
                                  PointVariableRef, FiniteVariableRef, 
                                  DerivativeRef}

const UserDecisionVariableRef = Union{InfiniteVariableRef, PointVariableRef,
                                      FiniteVariableRef, DerivativeRef}

const ScalarParameterRef = Union{IndependentParameterRef, FiniteParameterRef}

"""
    InfOptConstraintRef

A `DataType` for constraints that are in `InfiniteModel`s

**Fields**
- `model::InfiniteModel`: Infinite model.
- `index::InfOptConstraintIndex`: Index of the constraint in model.
"""
struct InfOptConstraintRef
    model::InfiniteModel
    index::InfOptConstraintIndex
end

################################################################################
#                            PARAMETER BOUND METHODS
################################################################################
## Modify parameter dictionary to expand any multidimensional parameter keys
# Case where tuple is already in correct form
function _expand_parameter_tuple(
    param_restrictions::NTuple{N, Pair{GeneralVariableRef, IntervalDomain}}
    )::Dict{GeneralVariableRef, IntervalDomain} where {N}
    return Dict{GeneralVariableRef, IntervalDomain}(param_restrictions...)
end

# Case where tuple contains vectors
function _expand_parameter_tuple(
    param_restrictions::NTuple{N, Pair{<:Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}, IntervalDomain}}
    )::Dict{GeneralVariableRef, IntervalDomain} where {N}
    # Initialize new dictionary
    new_dict = Dict{GeneralVariableRef, IntervalDomain}()
    # Find vector keys and expand
    for (key, set) in param_restrictions
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

# Handle symbolic domain input 
_process_domain(d::Real) = IntervalDomain(d, d)
function _process_domain(d::Vector{<:Real})
    if length(d) != 2
        error("Invalid domain restriction syntax. Expected restrictions of the ",
              "format `DomainRestrictions(restrictions...)` where each ",
              "argument in `restrictions` is one of the following forms:",
              "\n- `pref => value`\n- `pref => [lb, ub]`\n- `prefs => value`",
              "\n- `prefs => [lb, ub]`")
    end
    return IntervalDomain(d[1], d[2])
end

# Case where tuple has symbolic domains
function _expand_parameter_tuple(
    param_restrictions::NTuple{N, Pair{<:Any, <:Union{Real, Vector{<:Real}}}}
    )::Dict{GeneralVariableRef, IntervalDomain} where {N}
    new_tuple = Tuple(p[1] => _process_domain(p[2]) for p in param_restrictions)
    return _expand_parameter_tuple(new_tuple)
end

# Case where tuple contains other stuff
function _expand_parameter_tuple(param_restrictions)
    error("Invalid domain restriction syntax. Expected restrictions of the ",
          "format `DomainRestrictions(restrictions...)` where each ",
          "argument in `restrictions` is one of the following forms:",
          "\n- `pref => value`\n- `pref => [lb, ub]`\n- `prefs => value`",
          "\n- `prefs => [lb, ub]`")
end

# Constructor for expanding array parameters
function DomainRestrictions(
    intervals::NTuple{N, Pair}
    )::DomainRestrictions{GeneralVariableRef} where {N}
    return DomainRestrictions(_expand_parameter_tuple(intervals))
end

# Convenient constructor
function DomainRestrictions(args...)::DomainRestrictions{GeneralVariableRef}
    return DomainRestrictions(args)
end

# Default method
function DomainRestrictions()::DomainRestrictions{GeneralVariableRef}
    return DomainRestrictions(Dict{GeneralVariableRef, IntervalDomain}())
end

# Make dictionary accessor
function intervals(dr::DomainRestrictions{P})::Dict{P, IntervalDomain} where {P}
    return dr.intervals
end

# Extend simple 1 argument Base dispatches
for op = (:length, :isempty, :keys, :iterate)
    @eval Base.$op(restrictions::DomainRestrictions) = $op(intervals(restrictions))
end

# Extend simple 2 argument Base dispatches where the second is arbitrary
for op = (:getindex, :haskey, :iterate)
    @eval Base.$op(restrictions::DomainRestrictions, arg) = $op(intervals(restrictions), arg)
end

# Extend Base.:(==)
function Base.:(==)(
    restrictions1::DomainRestrictions, 
    restrictions2::DomainRestrictions
    )::Bool
    return intervals(restrictions1) == intervals(restrictions2)
end

# Extend Base.copy
function Base.copy(restrictions::DomainRestrictions) 
    return DomainRestrictions(copy(intervals(restrictions)))
end

# Extend Base.setindex!
function Base.setindex!(
    dr::DomainRestrictions{P}, 
    value::IntervalDomain,
    index::P
    )::IntervalDomain where {P}
    return intervals(dr)[index] = value
end

# Extend Base.delete!
function Base.delete!(
    dr::DomainRestrictions{P}, 
    key::P
    )::DomainRestrictions{P} where {P}
    delete!(intervals(dr), key)
    return dr
end

# Extend Base.merge
function Base.merge(
    dr1::DomainRestrictions{P},
    dr2::DomainRestrictions{P}
    )::DomainRestrictions{P} where {P}
    new_dict = merge(intervals(dr1), intervals(dr2))
    return DomainRestrictions(new_dict)
end

# Extend Base.merge!
function Base.merge!(
    dr1::DomainRestrictions{P},
    dr2::DomainRestrictions{P}
    )::DomainRestrictions{P} where {P}
    merge!(intervals(dr1), intervals(dr2))
    return dr1
end

# Extend Base.filter
function Base.filter(
    f::Function,
    dr::DomainRestrictions{P}
    )::DomainRestrictions{P} where {P}
    new_dict = filter(f, intervals(dr))
    return DomainRestrictions(new_dict)
end
