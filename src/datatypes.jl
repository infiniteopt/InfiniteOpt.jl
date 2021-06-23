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
#                      GENERATIVE SUPPORT INFORMATION TYPES
################################################################################
"""
    AbstractGenerativeInfo

An abstract type for storing information about generating supports that are made 
based on existing supports as required by certain measures and/or derivatives 
that depend on a certain independent infinite parameter. Such as the case with 
internal collocation supports.
"""
abstract type AbstractGenerativeInfo end 

"""
    NoGenerativeSupports <: AbstractGenerativeInfo

A `DataType` to signify that no generative supports will be generated for the 
measures and/or the derivatives. Has no fields.
"""
struct NoGenerativeSupports <: AbstractGenerativeInfo end

"""
    UniformGenerativeInfo <: AbstractGenerativeInfo

A `DataType` for generative supports that will be generated in a uniform manner 
over finite elements (i.e., in between the existing supports). These generative 
supports are described by the `support_basis` which lie in a nominal domain [0, 1]. 
The constructor is of the form:
```
    UniformGenerativeInfo(support_basis::Vector{<:Real}, label::DataType, 
                          [lb::Real = 0, ub::Real = 1])
```
where the `support_basis` is defined over [`lb`, `ub`].

**Fields**
- `support_basis::Vector{Float64}`: The basis of generative supports defined in 
   [0, 1] that will be transformed for each finite element.
- `label::DataType`: The unique label to be given to each generative support.
"""
struct UniformGenerativeInfo <: AbstractGenerativeInfo
    support_basis::Vector{Float64}
    label::DataType
    function UniformGenerativeInfo(basis::Vector{<:Real}, label::DataType, 
                                   lb::Real = 0, ub::Real = 1)
        if minimum(basis) < lb || maximum(basis) > ub
            error("Support basis violate the given lower and upper bounds. " * 
                  "Please specify the appropriate lower bound and upper bounds.")
        end
        return new((basis .- lb) ./ (ub - lb),  label)
    end
end

# Extend Base.:(==)
function Base.:(==)(info1::UniformGenerativeInfo, info2::UniformGenerativeInfo)::Bool 
    return info1.support_basis == info2.support_basis && info1.label == info2.label
end

################################################################################
#                       BASIC DERIVATIVE EVALUATION TYPES
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
collocation). Such methods can be used with derivatives that depend on independent 
infinite parameters, but cannot be used for ones that depend on dependent parameters.
"""
abstract type GenerativeDerivativeMethod <: AbstractDerivativeMethod end 

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
    FiniteDifference{T <: FDTechnique} <: NonGenerativeDerivativeMethod

A `DataType` for information about finite difference method applied to 
a derivative evaluation. Note that the constructor is of the form:
```julia 
    FiniteDifference([technique::FDTechnique = Backward()],
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
- `technique::T`: Mathematical technqiue behind finite difference
- `add_boundary_constraint::Bool`: Indicate if the boundary constraint should be 
  included in the transcription (e.g., the terminal boundary backward equation for 
  backward difference)
"""
struct FiniteDifference{T <: FDTechnique} <: NonGenerativeDerivativeMethod 
    technique::T
    add_boundary_constraint::Bool
    # set the constructor 
    function FiniteDifference(technique::T = Backward(), 
        add_boundary_constr::Bool = true) where {T <: FDTechnique}
        return new{T}(technique, add_boundary_constr)
    end
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
    IndependentParameter{T <: InfiniteScalarDomain,
                         M <: AbstractDerivativeMethod,
                         I <: AbstractGenerativeInfo} <: ScalarParameter

A `DataType` for storing independent scalar infinite parameters.

**Fields**
- `domain::T`: The infinite domain that characterizes the parameter.
- `supports::DataStructures.SortedDict{Float64, Set{DataType}}`: The support points
   used to discretize the parameter and their associated type labels stored as
   `DataTypes`s which should be a subtype of [`AbstractSupportLabel`](@ref).
- `sig_digits::Int`: The number of significant digits used to round the support values.
- `derivative_method::M`: The derivative evaluation method used for derivatives that
   are conducted with respect to this parameter.
- `gnerative_supp_info::I`: The info associated with any generative supports that will 
   need to be generated for measures and/or derivatives based on existing supports. 
"""
struct IndependentParameter{T <: InfiniteScalarDomain, 
                            M <: AbstractDerivativeMethod,
                            I <: AbstractGenerativeInfo} <: ScalarParameter
    domain::T
    supports::DataStructures.SortedDict{Float64, Set{DataType}} # Support to label set
    sig_digits::Int
    derivative_method::M
    generative_supp_info::I
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
- `supports::Dict{Vector{Float64}, Set{DataType}}`: Support dictionary where keys
              are supports and the values are the set of labels for each support.
- `sig_digits::Int`: The number of significant digits used to round the support values.
- `derivative_methods::Vector{M}`: The derivative evaluation methods associated with 
  each parameter.
"""
struct DependentParameters{T <: InfiniteArrayDomain, 
                           M <: NonGenerativeDerivativeMethod} <: InfOptParameter
    domain::T
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
- `parameter_func_indices::Vector{ParameterFunctionIndex}`: Indices of dependent
   infinite parameter functions.
- `infinite_var_indices::Vector{InfiniteVariableIndex}`: Indices of dependent
   infinite variables.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `in_objective::Bool`: Is this used in objective? This should be true only for finite parameters.
- `generative_measures::Vector{MeasureIndex}`: Indices of measures that use `parameter.generative_supp_info`.
- `has_internal_supports::Bool`: Does this parameter have internal supports?
- `has_generative_supports::Bool`: Have any generative supports been added?
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
    generative_measures::Vector{MeasureIndex}
    has_internal_supports::Bool
    has_generative_supports::Bool
    has_deriv_constrs::Bool
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
                                  InfOptConstraintIndex[], false, MeasureIndex[], 
                                  false, false, false)
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
- `has_internal_supports::Bool`: Does this parameter have internal supports?
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
    has_internal_supports::Bool
    has_deriv_constrs::Vector{Bool}
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
                                 false, zeros(Bool, length(names)))
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

An abstract type to define data for measures to define the behavior of
[`Measure`](@ref).
"""
abstract type AbstractMeasureData end

"""
    DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                        Vector{<:JuMP.AbstractVariableRef}},
                        N, B <: Union{Float64, Vector{Float64}},
                        F <: Function
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
- `weight_function::F`: Weighting function ``w`` must map an individual
                               support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bound in accordance with ``T``, this denotes the
                    intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bound.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                           Vector{<:JuMP.AbstractVariableRef}},
                           N, B <: Union{Float64, Vector{Float64}},
                           F <: Function
                           } <: AbstractMeasureData
    parameter_refs::P
    coefficients::Vector{Float64}
    supports::Array{Float64, N} # supports are stored column-wise
    label::DataType # label that will used when the supports are added to the model
    weight_function::F # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function DiscreteMeasureData(
        param_ref::V, coeffs::Vector{<:Real},
        supps::Vector{<:Real}, 
        label::DataType,
        weight_func::F,
        lower_bound::Real,
        upper_bound::Real,
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F <: Function}
        return new{V, 1, Float64, F}(param_ref, coeffs, supps, label, weight_func,
                                     lower_bound, upper_bound, expect)
    end
    # multi constructor
    function DiscreteMeasureData(
        param_refs::Vector{V}, 
        coeffs::Vector{<:Real},
        supps::Matrix{<:Real}, 
        label::DataType,
        weight_func::F,
        lower_bound::Vector{<:Real},
        upper_bound::Vector{<:Real},
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F <: Function}
        return new{Vector{V}, 2, Vector{Float64}, F}(param_refs, coeffs, supps,
                                                     label, weight_func, lower_bound,
                                                     upper_bound, expect)
    end
end

"""
    FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                  Vector{<:JuMP.AbstractVariableRef}},
                                  B <: Union{Float64, Vector{Float64}},
                                  I <: AbstractGenerativeInfo,
                                  F1 <: Function,
                                  F2 <: Function
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

For 1-dimensional measures over independent infinite parameters, the 
`generative_supp_info` specifies the info needed to make generative supports based 
on those with that exist with `label`. Note that only 1 kind of generative 
supports are allowed for each infinite parameter.

**Fields**
- `parameter_refs::P`: The infinite parameter(s) over which the integration occurs.
                     These can be comprised of multiple independent parameters,
                     but dependent parameters cannot be mixed with other types.
- `coeff_function::F1`: Coefficient generation function making ``\\alpha_i``
                              for the above measure abstraction. It should take
                              all the supports as input (formatted as an Array)
                              and return the corresponding vector of coefficients.
- `min_num_supports::Int`: Specifies the minimum number of supports ``\\tau_i``
                       desired in association with `parameter_refs` and `label`.
- `label::DataType`: Label for the support points ``\\tau_i`` which are/will be
                   stored in the infinite parameter(s), stemming from [`AbstractSupportLabel`](@ref).
- `generative_supp_info::I`: Information needed to generate supports based on other 
   existing ones.
- `weight_function::F2`: Weighting function ``w`` must map an individual
                              support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bounds in accordance with ``T``, this denotes the
                  intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bounds.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                     Vector{<:JuMP.AbstractVariableRef}},
                                     B <: Union{Float64, Vector{Float64}},
                                     I <: AbstractGenerativeInfo,
                                     F1 <: Function,
                                     F2 <: Function
                                     } <: AbstractMeasureData
    parameter_refs::P
    coeff_function::F1 # supports (excluding generative)--> coefficient vector (includes generative)
    min_num_supports::Int # minimum number of supports
    label::DataType # support label of included supports
    generative_supp_info::I
    weight_function::F2 # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function FunctionalDiscreteMeasureData(
        param_ref::V, 
        coeff_func::F1,
        num_supps::Int, 
        label::DataType,
        gen_info::I,
        weight_func::F2,
        lower_bound::Real,
        upper_bound::Real,
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, I <: AbstractGenerativeInfo,
                 F1 <: Function, F2 <: Function}
        return new{V, Float64, I, F1, F2}(param_ref, coeff_func, num_supps, label, 
                                          gen_info, weight_func, lower_bound, 
                                          upper_bound, expect)
    end
    # multi constructor
    function FunctionalDiscreteMeasureData(
        param_refs::Vector{V},
        coeff_func::F1,
        num_supps::Int, 
        label::DataType,
        weight_func::F2,
        lower_bound::Vector{<:Real},
        upper_bound::Vector{<:Real},
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F1 <: Function, F2 <: Function}
        return new{Vector{V}, Vector{Float64}, NoGenerativeSupports, F1, F2}(
            param_refs, coeff_func, num_supps, label, NoGenerativeSupports(), 
            weight_func, lower_bound, upper_bound, expect)
    end
end

# Convenient Dispatch constructor 
function FunctionalDiscreteMeasureData(
    param_refs::Vector{V},
    coeff_func::Function,
    num_supps::Int, 
    label::DataType,
    info::NoGenerativeSupports,
    weight_func::Function,
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real},
    expect::Bool
    ) where {V <: JuMP.AbstractVariableRef}
    return FunctionalDiscreteMeasureData(param_refs, coeff_func, num_supps,
                                         label, weight_func, lower_bound,
                                         upper_bound, expect)
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
    MeasureData{M <: Measure} <: AbstractDataObject

A mutable `DataType` for storing [`Measure`](@ref)s and their data.

**Fields**
- `measure::M`: The measure structure.
- `name::String`: The base name used for printing `name(meas_expr d(par))`.
- `measure_indices::Vector{MeasureIndex}`: Indices of dependent measures.
- `constraint_indices::Vector{InfOptConstraintIndex}`: Indices of dependent constraints.
- `derivative_indices::Vector{DerivativeIndex}`: Indices of dependent derivatives.
- `in_objective::Bool`: Is this used in objective?
"""
mutable struct MeasureData{M <: Measure} <: AbstractDataObject
    measure::M
    name::String
    measure_indices::Vector{MeasureIndex}
    constraint_indices::Vector{InfOptConstraintIndex}
    derivative_indices::Vector{DerivativeIndex}
    in_objective::Bool
end

function MeasureData(measure::M, name::String = "measure") where {M <: Measure}
    return MeasureData{M}(measure, name, MeasureIndex[], InfOptConstraintIndex[], 
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
- `param_functions::MOIUC.CleverDict{ParameterFunctionIndex, ParameterFunctionData{ParameterFunction}}`: 
  The infinite parameter functions and their mapping information.
- `infinite_vars::MOIUC.CleverDict{InfiniteVariableIndex, <:VariableData{<:InfiniteVariable}}`:
   The infinite variables and their mapping information.
- `semi_infinite_vars::MOIUC.CleverDict{SemiInfiniteVariableIndex, <:VariableData{<:SemiInfiniteVariable}}`:
   The semi-infinite variables and their mapping information.
- `semi_lookup::Dict{<:Tuple, SemiInfiniteVariableIndex}`: Look-up if a variable already already exists.
- `point_vars::MOIUC.CleverDict{PointVariableIndex, <:VariableData{<:PointVariable}}`:
   The point variables and their mapping information.
- `point_lookup::Dict{<:Tuple, PointVariableIndex}`: Look-up if a variable already exists.
- `finite_vars::MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}`:
   The finite variables and their mapping information.
- `name_to_var::Union{Dict{String, AbstractInfOptIndex}, Nothing}`:
   Field to help find a variable given the name.
- `derivatives::MOIUC.CleverDict{DerivativeIndex, <:VariableData{<:Derivative}}`:
  The derivatives and their mapping information.
- `deriv_lookup::Dict{<:Tuple, DerivativeIndex}`: Map derivative variable-parameter 
  pairs to a derivative index to prevent duplicates.
- `measures::MOIUC.CleverDict{MeasureIndex, <:MeasureData}`:
   The measures and their mapping information.
- `integral_defaults::Dict{Symbol}`:
   The default keyword arguments for [`integral`](@ref).
- `constraints::MOIUC.CleverDict{InfOptConstraintIndex, <:ConstraintData}`:
   The constraints and their mapping information.
- `constraint_restrictions::Dict{InfOptConstraintIndex, <:DomainRestrictions}` Map constraints 
  to their domain restrictions if they have any.
- `name_to_constr::Union{Dict{String, InfOptConstraintIndex}, Nothing}`:
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
    param_functions::MOIUC.CleverDict{ParameterFunctionIndex, <:ParameterFunctionData}

    # Variable Data
    infinite_vars::MOIUC.CleverDict{InfiniteVariableIndex, <:VariableData{<:InfiniteVariable}}
    semi_infinite_vars::MOIUC.CleverDict{SemiInfiniteVariableIndex, <:VariableData{<:SemiInfiniteVariable}}
    semi_lookup::Dict{<:Tuple, SemiInfiniteVariableIndex}
    point_vars::MOIUC.CleverDict{PointVariableIndex, <:VariableData{<:PointVariable}}
    point_lookup::Dict{<:Tuple, PointVariableIndex}
    finite_vars::MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}
    name_to_var::Union{Dict{String, AbstractInfOptIndex}, Nothing}

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
function InfiniteModel(; 
    OptimizerModel::Function = TranscriptionModel,
    kwargs...
    )::InfiniteModel
    return InfiniteModel(# Parameters
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
                         nothing,
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
                         # Object dictionary
                         Dict{Symbol, Any}(),
                         # Optimize data
                         nothing, OptimizerModel(; kwargs...), false,
                         # Extensions
                         Dict{Symbol, Any}()
                         )
end

## Set the optimizer_constructor depending on what it is
# MOI.OptimizerWithAttributes
function _set_optimizer_constructor(
    model::InfiniteModel,
    constructor::MOI.OptimizerWithAttributes
    )::Nothing
    model.optimizer_constructor = constructor.optimizer_constructor
    return
end

# No attributes
function _set_optimizer_constructor(model::InfiniteModel, constructor)::Nothing
    model.optimizer_constructor = constructor
    return
end

# Dispatch for InfiniteModel call with optimizer constructor
function InfiniteModel(
    optimizer_constructor;
    OptimizerModel::Function = TranscriptionModel,
    kwargs...
    )::InfiniteModel
    model = InfiniteModel()
    model.optimizer_model = OptimizerModel(optimizer_constructor; kwargs...)
    _set_optimizer_constructor(model, optimizer_constructor)
    return model
end

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
JuMP.object_dictionary(model::InfiniteModel)::Dict{Symbol, Any} = model.obj_dict

"""
    Base.empty!(model::InfiniteModel)::InfiniteModel

Clear out `model` of everything except the optimizer information and return the 
cleared model. 
"""
function Base.empty!(model::InfiniteModel)::InfiniteModel 
    # Clear everything except the solver information
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
    empty!(model.obj_dict)
    empty!(model.optimizer_model)
    model.ready_to_optimize = false
    empty!(model.ext)
    return model
end

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

# Make dumby model type for calling @expression 
struct _DumbyModel <: JuMP.AbstractModel end
const _Model = _DumbyModel()

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
