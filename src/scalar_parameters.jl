################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::IndependentParameterIndex
                               )::IndependentParameterRef
    return IndependentParameterRef(model, index)
end

function dispatch_variable_ref(model::InfiniteModel,
                               index::FiniteParameterIndex
                               )::FiniteParameterRef
    return FiniteParameterRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::ScalarParameterData{<:IndependentParameter}
                          )::IndependentParameterIndex
    return MOIUC.add_item(model.independent_params, object)
end

function _add_data_object(model::InfiniteModel,
                          object::ScalarParameterData{<:FiniteParameter}
                          )::FiniteParameterIndex
    return MOIUC.add_item(model.finite_params, object)
end

# Extend _data_dictionary
function _data_dictionary(pref::IndependentParameterRef)::MOIUC.CleverDict
    return JuMP.owner_model(pref).independent_params
end

function _data_dictionary(pref::FiniteParameterRef)::MOIUC.CleverDict
    return JuMP.owner_model(pref).finite_params
end

# Extend _data_object
function _data_object(pref::ScalarParameterRef
                     )::ScalarParameterData
    return _data_dictionary(pref)[JuMP.index(pref)]
end

################################################################################
#                             CORE OBJECT METHODS
################################################################################
# Extend _core_variable_object for IndependentParameterRefs
function _core_variable_object(pref::IndependentParameterRef)::IndependentParameter
    return _data_object(pref).parameter
end

# Extend _core_variable_object for FiniteParameterRefs
function _core_variable_object(pref::FiniteParameterRef)::FiniteParameter
    return _data_object(pref).parameter
end

# Extend _set_core_variable_object for IndependentParameterRefs
function _set_core_variable_object(pref::IndependentParameterRef,
                                   param::IndependentParameter)::Nothing
    _data_object(pref).parameter = param
    return
end

# Extend _set_core_variable_object for FiniteParameterRefs
function _set_core_variable_object(pref::FiniteParameterRef,
                                   param::FiniteParameter)::Nothing
    _data_object(pref).parameter = param
    return
end

################################################################################
#                             PARAMETER DEFINITION
################################################################################
# Internal structure for building InfOptParameters
mutable struct _ParameterInfoExpr
    has_lb::Bool
    lower_bound::Any
    has_ub::Bool
    upper_bound::Any
    has_dist::Bool
    distribution::Any
    has_set::Bool
    set::Any
end

# Default constructor
function _ParameterInfoExpr(; lower_bound = NaN, upper_bound = NaN,
                            distribution = NaN, set = NaN)
    # isnan(::Expr) is not defined so we need to do !== NaN
    return _ParameterInfoExpr(lower_bound !== NaN, lower_bound,
                              upper_bound !== NaN, upper_bound,
                              distribution !== NaN, distribution,
                              set !== NaN, set)
end

# Internal function for use in processing valid key word arguments
function _is_set_keyword(kw::Expr)
    return kw.args[1] in [:set, :lower_bound, :upper_bound, :distribution]
end

# Extend to assist in building InfOptParameters
function JuMP._set_lower_bound_or_error(_error::Function,
                                        info::_ParameterInfoExpr, lower)
    info.has_lb && _error("Cannot specify parameter lower_bound twice")
    info.has_dist && _error("Cannot specify parameter lower_bound and " *
                            "distribution")
    info.has_set && _error("Cannot specify parameter lower_bound and set")
    info.has_lb = true
    info.lower_bound = lower
    return
end

# Extend to assist in building InfOptParameters
function JuMP._set_upper_bound_or_error(_error::Function,
                                        info::_ParameterInfoExpr, upper)
    info.has_ub && _error("Cannot specify parameter upper_bound twice")
    info.has_dist && _error("Cannot specify parameter upper_bound and " *
                            "distribution")
    info.has_set && _error("Cannot specify parameter upper_bound and set")
    info.has_ub = true
    info.upper_bound = upper
    return
end

# Extend to assist in building InfOptParameters
function _dist_or_error(_error::Function, info::_ParameterInfoExpr, dist)
    info.has_dist && _error("Cannot specify parameter distribution twice")
    (info.has_lb || info.has_ub) && _error("Cannot specify parameter " *
                                           "distribution and upper/lower bounds")
    info.has_set && _error("Cannot specify parameter distribution and set")
    info.has_dist = true
    info.distribution = dist
    return
end

# Extend to assist in building InfOptParameters
function _set_or_error(_error::Function, info::_ParameterInfoExpr, set)
    info.has_set && _error("Cannot specify variable fixed value twice")
    (info.has_lb || info.has_ub) && _error("Cannot specify parameter set and " *
                                           "upper/lower bounds")
    info.has_dist && _error("Cannot specify parameter set and distribution")
    info.has_set = true
    info.set = set
    return
end

# Construct an expression to build an infinite set (use with @infinite_macro)
function _constructor_set(_error::Function, info::_ParameterInfoExpr)
    if (info.has_lb || info.has_ub) && !(info.has_lb && info.has_ub)
        _error("Must specify both an upper bound and a lower bound")
    elseif info.has_lb
        check = :(isa($(info.lower_bound), Number))
        return :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a number."))
    elseif info.has_dist
        check = :(isa($(info.distribution), Distributions.UnivariateDistribution))
        return :($(check) ? UniDistributionSet($(info.distribution)) : error("Distribution must be a Distributions.UnivariateDistribution."))
    elseif info.has_set
        check1 = :(isa($(info.set), InfiniteScalarSet))
        check2 = :(isa($(info.set), Distributions.UnivariateDistribution))
        return :($(check1) ? $(info.set) : ($(check2) ? UniDistributionSet($(info.set)) : error("Set must be a subtype of InfiniteScalarSet.")))
    else
        _error("Must specify upper/lower bounds, a distribution, or a set")
    end
end

# Check that supports don't violate the set bounds
function _check_supports_in_bounds(_error::Function,
                                   supports::Union{Real, Vector{<:Real}},
                                   set::AbstractInfiniteSet)
    if !supports_in_set(supports, set)
        _error("Supports violate the set domain bounds.")
    end
    return
end

"""
    build_parameter(_error::Function, set::InfiniteScalarSet,
                           [num_params::Int = 1; num_supports::Int = 0,
                           supports::Union{Real, Vector{<:Real}} = Number[],
                           independent::Bool = false,
                           sig_fig::Int = 5])::IndependentParameter

Returns a [`IndependentParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. Errors if supports violate the
bounds associated with `set`. This is meant to primarily serve as a
helper method for [`@independent_parameter`](@ref).

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> build_parameter(error, IntervalSet(0, 3), supports = Vector(0:3))
IndependentParameter{IntervalSet}(IntervalSet(0.0, 3.0), DataStructures.SortedDict(0.0 => Set([]),1.0 => Set([]),2.0 => Set([]),3.0 => Set([])))
```
"""
function build_parameter(_error::Function,
    set::S;
    num_supports::Int = 0,
    supports::Union{Real, Vector{<:Real}} = Real[],
    sig_fig::Int = 5,
    extra_kw_args...
    )::IndependentParameter{S} where {S <: InfiniteScalarSet}
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    label = Set{Symbol}([UserDefined])
    length_supports = length(supports)
    if length_supports != 0
        _check_supports_in_bounds(_error, supports, set)
        num_supports == 0 || @warn("Ignoring num_supports since supports is not empty.")
    elseif num_supports != 0
        supports, label = generate_support_values(set, num_supports = num_supports,
                                                  sig_figs = sig_fig)
    end
    supports_dict = DataStructures.SortedDict{Float64, Set{Symbol}}(
                                                   i => label for i in supports)
    if length(supports_dict) != length_supports
        @warn("Support points are not unique, eliminating redundant points.")
    end
    return IndependentParameter(set, supports_dict)
end

"""
    build_parameter(_error::Function, value::Real)::FiniteParameter

Returns a [`FiniteParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. This is meant to primarily serve as
a helper method for [`@finite_parameter`](@ref).

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> build_finite_parameter(error, 1)
FiniteParameter(1.0)
```
"""

function build_parameter(_error::Function, value::Real;
                         extra_kw_args...)::FiniteParameter
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    return FiniteParameter(value)
end

"""
    add_parameter(model::InfiniteModel,
                  p::Union{IndependentParameterRef, FiniteParameterRef},
                  [name::String = ""])::ParameterRef

Returns a [`ParameterRef`](@ref) associated with the parameter `p` that is added
to `model`. This adds a parameter to the model in a manner similar to
`JuMP.add_variable`. This can be used to add parameters with the use of
[`@infinite_parameter`](@ref). [`build_parameter`](@ref) should be used to
construct `p`.

**Example**
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> p = build_independent_parameter(error, IntervalSet(0, 3), supports = Vector(0:3))
IndependentParameter{IntervalSet}(IntervalSet(0.0, 3.0), DataStructures.SortedDict(0.0 => Set([]),1.0 => Set([]),2.0 => Set([]),3.0 => Set([])))

julia> param_ref = add_parameter(model, p, "name")
name
```
"""
function add_parameter(model::InfiniteModel, p::IndependentParameter,
                       name::String="")::GeneralVariableRef
    obj_num = model.last_object_num += 1
    param_num = model.last_param_num += 1
    data_object = ScalarParameterData(p, obj_num, param_num, name)
    obj_index = _add_data_object(model, data_object)
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

function add_parameter(model::InfiniteModel, p::FiniteParameter,
                       name::String="")::GeneralVariableRef
    data_object = ScalarParameterData(p, -1, -1, name)
    obj_index = _add_data_object(model, data_object)
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

"""
    used_by_constraint(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used by a constraint or false otherwise.

**Example**
```julia-repl
julia> used_by_constraint(t)
true
```
"""
function used_by_constraint(pref::ScalarParameterRef)::Bool
    return length(_data_object.constraint_indices)
end

"""
    used_by_measure(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used by a measure or false otherwise.

**Example**
```julia-repl
julia> used_by_measure(t)
false
```
"""
function used_by_measure(pref::ScalarParameterRef)::Bool
    return length(_data_object.measure_indices)
end

"""
    used_by_variable(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used by an infinite variable or false otherwise.

**Example**
```julia-repl
julia> used_by_variable(t)
true
```
"""
function used_by_variable(pref::ScalarParameterRef)::Bool
    return length(_data_object.infinite_var_indices)
end

"""
    is_used(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used in the model or false otherwise.

**Example**
```julia-repl
julia> is_used(t)
true
```
"""
function is_used(pref::ScalarParameterRef)::Bool
    return used_by_measure(pref) || used_by_constraint(pref) || used_by_variable(pref)
end


"""
    JuMP.name(pref::Union{IndependentParameterRef, FiniteParameterRef})::String

Extend the [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) function to
accomodate infinite parameters. Returns the name string associated with `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @independent_parameter(model, t in [0, 1]))
julia> name(t)
"t"
```
"""
function JuMP.name(pref::ScalarParameterRef)::String
    return _data_object(pref).name
end

"""
    JuMP.set_name(pref::ScalarParameterRef, name::String)

Extend the [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String))
function to accomodate infinite parameters. Set a new base name to be associated
with `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @independent_parameter(model, t in [0, 1]))
julia> set_name(t, "time")

julia> name(t)
"time"
```
"""
function JuMP.set_name(pref::ScalarParameterRef, name::String)
    _data_object(pref).name = name
    return
end

"""
    num_parameters(model::InfiniteModel)::Int

Return the number of scalar parameters currently present in `model`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> num_parameters(model)
1
```
"""
function num_parameters(model::InfiniteModel)::Int
    return num_independent_parameters(model) + num_finite_parameters(model)
end

"""
    num_independent_parameters(model::InfiniteModel)::Int

Return the number of independent parameters currently present in `model`.
"""
function num_independent_parameters(model::InfiniteModel)::Int
    return length(model.independent_params)
end

"""
    num_finite_parameters(model::InfiniteModel)::Int

Return the number of finite parameters currently present in `model`.
"""
function num_finite_parameters(model::InfiniteModel)::Int
    return length(model.finite_params)
end

# Internal functions
_parameter_set(pref::IndependentParameterRef) = _core_variable_object(pref).set
_parameter_supports(pref::IndependentParameterRef) = _core_variable_object(pref).supports
function _update_parameter_set(pref::IndependentParameterRef, set::AbstractInfiniteSet)
    param = _core_variable_object(pref)
    supports = param.supports
    new_param = IndependentParameter(set, supports)
    _set_core_variable_object(pref, new_param)
    if is_used(pref)
        set_optimizer_model_ready(JuMP.owner_model(pref), false)
    end
    return
end
function _update_parameter_supports(pref::IndependentParameterRef,
                                    supports::DataStructures.SortedDict{Float64, Set{Symbol}})
    param = _core_variable_object(pref)
    set = param.set
    new_param = IndependentParameter(set, supports)
    _set_core_variable_object(pref, new_param)
    if is_used(pref)
        set_optimizer_model_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    infinite_set(pref::IndependentParameterRef)::AbstractInfiniteSet

Return the infinite set associated with `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> infinite_set(t)
[0, 1]
```
"""
function infinite_set(pref::IndependentParameterRef)::AbstractInfiniteSet
    return _parameter_set(pref)
end

"""
    set_infinite_set(pref::IndependentParameterRef, set::AbstractInfiniteSet)

Specify the infinite set of `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_infinite_set(t, IntervalSet(0, 2))

julia> infinite_set(t)
[0, 2]
```
"""
function set_infinite_set(pref::IndependentParameterRef, set::AbstractInfiniteSet)
    _update_parameter_set(pref, set)
    return
end

"""
    JuMP.has_lower_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_lower_bound` function to accomodate infinite parameters.
Return true if the set associated with `pref` has a defined lower bound or if a
lower bound can be found. Extensions with user-defined infinite set types
should extend `JuMP.has_lower_bound(set::NewType)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> has_lower_bound(t)
true
```
"""
function JuMP.has_lower_bound(pref::IndependentParameterRef)::Bool
    set = _parameter_set(pref)
    return JuMP.has_lower_bound(set)
end

"""
    JuMP.lower_bound(pref::IndependentParameterRef)::Number

Extend the `JuMP.lower_bound` function to accomodate infinite parameters.
Returns the lower bound associated with the infinite set. Errors if such a bound
is not well-defined.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> lower_bound(t)
0.0
```
"""
function JuMP.lower_bound(pref::IndependentParameterRef)::Number
    set = _parameter_set(pref)
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    return JuMP.lower_bound(set)
end

"""
    JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Number)

Extend the `JuMP.set_lower_bound` function to accomodate infinite parameters.
Updates the infinite set lower bound if such an operation is supported. Set
extensions that seek to employ this should extend
`JuMP.set_lower_bound(set::NewType, lower::Number)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_lower_bound(t, -1)

julia> lower_bound(t)
-1.0
```
"""
function JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Number)
    set = _parameter_set(pref)
    new_set = JuMP.set_lower_bound(set, lower)
    _update_parameter_set(pref, new_set)
    return
end

"""
    JuMP.has_upper_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_upper_bound` function to accomodate infinite parameters.
Return true if the set associated with `pref` has a defined upper bound or if a
upper bound can be found. Extensions with user-defined sets should extend
`JuMP.has_upper_bound(set::NewType)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> has_upper_bound(t)
true
```
"""
function JuMP.has_upper_bound(pref::IndependentParameterRef)::Bool
    set = _parameter_set(pref)
    return JuMP.has_upper_bound(set)
end

"""
    JuMP.upper_bound(pref::IndependentParameterRef)::Number

Extend the `JuMP.upper_bound` function to accomodate infinite parameters.
Returns the upper bound associated with the infinite set. Errors if such a bound
is not well-defined. Extensions with user-defined set types should extend
`JuMP.has_upper_bound(set::NewType)` and `JuMP.upper_bound(set::NewType)` if
appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> upper_bound(t)
1.0
```
"""
function JuMP.upper_bound(pref::IndependentParameterRef)::Number
    set = _parameter_set(pref)
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    return JuMP.upper_bound(set)
end

"""
    JuMP.set_upper_bound(pref::IndependentParameterRef, lower::Number)

Extend the `JuMP.set_upper_bound` function to accomodate infinite parameters.
Updates the infinite set upper bound if and only if it is an IntervalSet. Errors
otherwise. Extensions with user-defined infinite sets should extend
`JuMP.set_upper_bound(set::NewType, upper::Number)` if appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_upper_bound(t, 2)

julia> upper_bound(t)
2.0
```
"""
function JuMP.set_upper_bound(pref::IndependentParameterRef, upper::Number)
    set = _parameter_set(pref)
    new_set = JuMP.set_upper_bound(set, upper)
    _update_parameter_set(pref, new_set)
    return
end

"""
    num_supports(pref::IndependentParameterRef)::Int

Return the number of support points associated with `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> num_supports(t)
2
```
"""
function num_supports(pref::IndependentParameterRef)::Int
    return length(_parameter_supports(pref))
end

"""
    has_supports(pref::IndependentParameterRef)::Bool

Return true if `pref` has supports or false otherwise.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> has_supports(t)
true
```
"""
has_supports(pref::IndependentParameterRef)::Bool = num_supports(pref) > 0

"""
    supports(pref::IndependentParameterRef)::Vector

Return the support points associated with `pref`. Errors if there are no
supports.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> supports(t)
2-element Array{Int64,1}:
 0
 1
```
"""
function supports(pref::IndependentParameterRef)::Vector
    has_supports(pref) || error("Parameter $pref does not have supports.")
    return _parameter_supports(pref)
end

"""
    set_supports(pref::IndependentParameterRef, supports::Vector{<:Number}; [force = false])

Specify the support points for `pref`. Errors if the supports violate the bounds
associated with the infinite set. Warns if the points are not unique. If `force`
this will overwrite exisiting supports otherwise it will error if there are
existing supports.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_supports(t, [0, 1])

julia> supports(t)
2-element Array{Int64,1}:
 0
 1
```
"""
function set_supports(pref::IndependentParameterRef, supports::Vector{<:Number};
                      force = false)
    set = _parameter_set(pref)
    _check_supports_in_bounds(error, supports, set)
    if has_supports(pref) && !force
        error("Unable set supports for $pref since it already has supports." *
              " Consider using `add_supports` or use set `force = true` to " *
              "overwrite the existing supports.")
    end
    if !(is_independent(pref)) &&
       sum(values(pref.model.param_to_group_id) .== group_id(pref)) > 1
        _update_parameter_supports(pref, supports)
    else
        unique_supports = unique(supports)
        if length(unique_supports) != length(supports)
            @warn("Support points are not unique, eliminating redundant points.")
        end
        _update_parameter_supports(pref, unique_supports)
    end
    return
end

"""
    add_supports(pref::IndependentParameterRef, supports::Union{Number, Vector{<:Number}})

Add additional support points for `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> add_supports(t, 0.5)

julia> supports(t)
3-element Array{Float64,1}:
 0.0
 1.0
 0.5

julia> add_supports(t, [0.25, 1])

julia> supports(t)
4-element Array{Float64,1}:
 0.0
 1.0
 0.5
 0.25
```
"""
function add_supports(pref::IndependentParameterRef, supports::Union{Number,
                                                          Vector{<:Number}})
    set = _parameter_set(pref)
    _check_supports_in_bounds(error, supports, set)
    current_supports = _parameter_supports(pref)
    if !(is_independent(pref)) &&
       sum(values(pref.model.param_to_group_id) .== group_id(pref)) > 1
        new_supports = [current_supports; supports]
    else
        new_supports = unique([current_supports; supports])
    end
    _update_parameter_supports(pref, new_supports)
    return
end

"""
    delete_supports(pref::IndependentParameterRef)

Delete the support points for `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> delete_supports(t)

julia> supports(t)
ERROR: Parameter t does not have supports.
```
"""
function delete_supports(pref::IndependentParameterRef)
    _update_parameter_supports(pref, Int[])
    return
end

"""
    is_finite_parameter(pref::IndependentParameterRef)::Bool

Return a `Bool` indicating if `pref` is a finite parameter.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @finite_parameter(model, cost, 42))
julia> is_finite_parameter(cost)
true
```
"""
function is_finite_parameter(pref::IndependentParameterRef)::Bool
    set = infinite_set(pref)
    if isa(set, IntervalSet) && set.lower_bound == set.upper_bound
        return true
    end
    return false
end

"""
    JuMP.value(pref::IndependentParameterRef)::Number

Return the value of `pref` so long as it is a finite parameter. Errors if it is
an infinite parameter.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @finite_parameter(model, cost, 42))
julia> value(cost)
42
```
"""
function JuMP.value(pref::IndependentParameterRef)::Number
    is_finite_parameter(pref) || error("$pref is an infinite parameter.")
    return supports(pref)[1]
end

"""
    JuMP.set_value(pref::IndependentParameterRef, value::Number)

Set the value of `pref` so long as it is a finite parameter. Errors if it is
an infinite parameter.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @finite_parameter(model, cost, 42))
julia> set_value(cost, 27)

julia> value(cost)
27
```
"""
function JuMP.set_value(pref::IndependentParameterRef, value::Number)
    is_finite_parameter(pref) || error("$pref is an infinite parameter.")
    set_infinite_set(pref, IntervalSet(value, value))
    set_supports(pref, [value], force = true)
    return
end

"""
    fill_in_supports!(model::InfiniteModel; [num_supports::Int = 10,
                      sig_fig::Int = 5])

Automatically generate support points for all infinite parameters in model. User
can specify the number of significant figures kept after decimal point for the
auto-generated supports wtih `sig_fig`. This calls
[`fill_in_supports!`](@ref fill_in_supports!(::ParameterRef)) for each parameter
in the model. See [`fill_in_supports!`](@ref fill_in_supports!(::ParameterRef))
for more information. Errors if one of the infinite set types is unrecognized.

**Example**
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel(); @infinite_parameter(model, 0 <= x <= 1);)
julia> fill_in_supports!(model, num_supports = 4, sig_fig = 3)

julia> supports(x)
4-element Array{Number,1}:
 0.0
 0.333
 0.667
 1.0
```
"""
function fill_in_supports!(model::InfiniteModel; num_supports::Int = 10,
                           sig_fig::Int = 5)
    for key in keys(model.params)
        pref = ParameterRef(model, key)
        fill_in_supports!(pref, num_supports = num_supports, sig_fig = sig_fig)
    end
    return
end

"""
    fill_in_supports!(pref::IndependentParameterRef; [num_supports::Int = 10,
                                           sig_fig::Int = 5])

Automatically generate support points for a particular infinite parameter `pref`.
Generating `num_supports` for the parameter. The supports are generated uniformly
if the underlying infinite set is an `IntervalSet` or they are generating randomly
accordingly to the distribution if the set is a `DistributionSet`.
User can specify the number of digits kept after decimal point for the
auto-generated supports wtih `sig_fig`. Extensions that use user defined
set types should extend [`generate_and_add_supports!`](@ref) and/or
[`generate_support_values`](@ref) as needed. Errors if the infinite set type is
not recognized.

**Example**
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel(); @infinite_parameter(model, 0 <= x <= 1);)
julia> fill_in_supports!(x, num_supports = 4, sig_fig = 3)

julia> supports(x)
4-element Array{Number,1}:
 0.0
 0.333
 0.667
 1.0

```
"""
function fill_in_supports!(pref::IndependentParameterRef;
                           num_supports::Int = 10, sig_fig::Int = 5)
    p = JuMP.owner_model(pref).params[JuMP.index(pref)]
    if length(p.supports) == 0
        generate_and_add_supports!(pref, p.set, num_supports = num_supports,
                                   sig_fig = sig_fig)
    end
    return
end

"""
    generate_and_add_supports!(pref::IndependentParameterRef,
                               set::AbstractInfiniteSet;
                               [num_supports::Int = 10, sig_fig::Int = 5])

Generate supports for `pref` via [`generate_support_values`](@ref) and add them
to `pref`. This is intended as an extendable internal method for
[`fill_in_supports!`](@ref fill_in_supports!(::ParameterRef)). Note that if
`pref` is part of a `DistributionSet` that features a multivariate distribution,
all the associated parameters with `pref` will also have supports added to them.
Most extensions that empoy user-defined infinite sets can typically enable this
by extending [`generate_support_values`](@ref). However, in some cases it may be
necessary to extend this when more complex operations need to take place then just
adding supports to a single infinite parameter (e.g., how we enable multivariate
distribution sets). Errors if the infinite set type is not recognized.
"""
function generate_and_add_supports!(pref::IndependentParameterRef,
                                    set::AbstractInfiniteSet;
                                    num_supports::Int = 10, sig_fig::Int = 5)
    add_supports(pref, generate_support_values(set, num_supports = num_supports,
                                               sig_fig = sig_fig))
    return
end
#=
"""
    parameter_by_name(model::InfiniteModel, name::String)::Union{ParameterRef,
                                                                 Nothing}

Return the parameter reference assoociated with a parameter name. Errors if
multiple parameters have the same name. Returns nothing if no such name exists.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> parameter_by_name(model, "t")
t
```
"""
function parameter_by_name(model::InfiniteModel,
                           name::String)::Union{ParameterRef, Nothing}
    if model.name_to_param === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_param = Dict{String, Int}()
        for (param, param_name) in model.param_to_name
            if haskey(model.name_to_param, param_name)
                # -1 is a special value that means this string does not map to
                # a unique variable name.
                model.name_to_param[param_name] = -1
            else
                model.name_to_param[param_name] = param
            end
        end
    end
    index = get(model.name_to_param, name, nothing)
    if index isa Nothing
        return
    elseif index == -1
        error("Multiple parameters have the name $name.")
    else
        return ParameterRef(model, index)
    end
end

"""
    all_parameters(model::InfiniteModel)::Vector{ParameterRef}

Return all of the infinite parameter references currently in `model`.

**Example**
```julia-repl
julia> all_parameters(model)
3-element Array{ParameterRef,1}:
 t
 x[1]
 x[2]
```
"""
function all_parameters(model::InfiniteModel)::Vector{ParameterRef}
    pref_list = Vector{ParameterRef}(undef, num_parameters(model))
    indexes = sort([index for index in keys(model.params)])
    counter = 1
    for index in indexes
        pref_list[counter] = ParameterRef(model, index)
        counter += 1
    end
    return pref_list
end
=#
#=
################################################################################
#                               DELETE FUNCTIONS
################################################################################


# Check if parameter is used by measure data and error if it is to prevent bad
# deleting behavior
function _check_param_in_data(pref::ParameterRef, data::AbstractMeasureData)
    prefs = parameter_refs(data)
    if (pref == prefs || pref in prefs)
        error("Unable to delete `$pref` since it is used to evaluate measures.")
    end
    return
end

# Used to update infinite variable when one of its parameters is deleted
function _update_infinite_variable(vref::InfiniteVariableRef)
    JuMP.set_name(vref, _root_name(vref))
    if used_by_measure(vref)
        for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
            JuMP.set_name(MeasureRef(JuMP.owner_model(vref), mindex),
                       _make_meas_name(JuMP.owner_model(vref).measures[mindex]))
        end
    end
    return
end

# Update point variable for which a parameter is deleted
function _update_point_variable(pvref::PointVariableRef)
    # update name if no alias was provided
    if !isa(findfirst(isequal('('), JuMP.name(pvref)), Nothing)
        JuMP.set_name(pvref, "")
    end
    if used_by_measure(pvref)
        for mindex in JuMP.owner_model(pvref).var_to_meas[JuMP.index(pvref)]
            JuMP.set_name(MeasureRef(JuMP.owner_model(pvref), mindex),
                      _make_meas_name(JuMP.owner_model(pvref).measures[mindex]))
        end
    end
    return
end

# Update a reduced variable associated with an infinite variable whose parameter
# was removed
function _update_reduced_variable(vref::ReducedInfiniteVariableRef,
                                  delete_index::Int)
    eval_supps = eval_supports(vref)
    new_supports = Dict{Int, Number}()
    for (index, support) in eval_supps
        if index < delete_index
            new_supports[index] = support
        elseif index > delete_index
            new_supports[index - 1] = support
        end
    end
    new_info = ReducedInfiniteInfo(infinite_variable_ref(vref), new_supports)
    JuMP.owner_model(vref).reduced_info[JuMP.index(vref)] = new_info
    # update measure dependencies
    if used_by_measure(vref)
        for mindex in JuMP.owner_model(vref).reduced_to_meas[JuMP.index(vref)]
            JuMP.set_name(MeasureRef(JuMP.owner_model(vref), mindex),
                       _make_meas_name(JuMP.owner_model(vref).measures[mindex]))
        end
    end
    return
end

"""
    JuMP.delete(model::InfiniteModel, pref::ParameterRef)

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)) to delete
infinite parameters and their dependencies. All variables, constraints, and
measure functions that depend on `pref` are updated to exclude it. Errors if the
parameter is contained in an `AbstractMeasureData` datatype that is employed by
a measure since the measure becomes invalid otherwise. Thus, measures that
contain this dependency must be deleted first. Note that
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) needs to be
extended to allow deletion of parameters when custom `AbstractMeasureData`
datatypes are used.

**Example**
```julia-repl
julia> print(model)
Min measure(g(t, x)*t + x) + z
Subject to
 z >= 0.0
 g(t, x) + z >= 42.0
 g(0.5, x) == 0
 t in [0, 6]
 x in [0, 1]

julia> delete(model, x)

julia> print(model)
Min measure(g(t)*t) + z
Subject to
 g(t) + z >= 42.0
 g(0.5) == 0
 z >= 0.0
 t in [0, 6]
```
"""
function JuMP.delete(model::InfiniteModel, pref::ParameterRef)
    @assert JuMP.is_valid(model, pref) "Parameter reference is invalid."
    # update optimizer model status
    if is_used(pref)
        set_optimizer_model_ready(model, false)
    end
    # update measures
    if used_by_measure(pref)
        # ensure deletion is okay (pref isn't used by measure data)
        for mindex in model.param_to_meas[JuMP.index(pref)]
            _check_param_in_data(pref, model.measures[mindex].data)
        end
        # delete dependence of measures on pref
        for mindex in model.param_to_meas[JuMP.index(pref)]
            if isa(model.measures[mindex].func, ParameterRef)
                model.measures[mindex] = Measure(zero(JuMP.AffExpr),
                                                 model.measures[mindex].data)
            else
                _remove_variable(model.measures[mindex].func, pref)
            end
            JuMP.set_name(MeasureRef(model, mindex),
                          _make_meas_name(model.measures[mindex]))
        end
        # delete mapping
        delete!(model.param_to_meas, JuMP.index(pref))
    end
    # update variables
    if used_by_variable(pref)
        # update infinite variables that depend on pref
        for vindex in model.param_to_vars[JuMP.index(pref)]
            # remove the parameter dependence
            vref = InfiniteVariableRef(model, vindex)
            prefs = raw_parameter_refs(vref)
            delete_index = findfirst(isequal(pref), prefs)
            deleteat!(prefs, delete_index)
            _update_infinite_variable(vref)
            # update any point variables that depend on vref accordingly
            if used_by_point_variable(vref)
                for pindex in model.infinite_to_points[vindex]
                    pvref = PointVariableRef(model, pindex)
                    deleteat!(raw_parameter_values(pvref), delete_index)
                    _update_point_variable(pvref)
                end
            end
            # update any reduced variables that depend on vref accordingly
            if used_by_reduced_variable(vref)
                for rindex in model.infinite_to_reduced[vindex]
                    rvref = ReducedInfiniteVariableRef(model,rindex)
                    _update_reduced_variable(rvref, delete_index)
                end
            end
        end
        # delete mapping
        delete!(model.param_to_vars, JuMP.index(pref))
    end
    # update constraints
    if used_by_constraint(pref)
        # update constraints in mapping to remove the parameter
        for cindex in model.param_to_constrs[JuMP.index(pref)]
            if isa(model.constrs[cindex].func, ParameterRef)
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr),
                                                      model.constrs[cindex].set)
            else
                _remove_variable(model.constrs[cindex].func, pref)
            end
        end
        # delete mapping
        delete!(model.param_to_constrs, JuMP.index(pref))
    end
    # delete parameter information stored in model
    delete!(model.params, JuMP.index(pref))
    delete!(model.param_to_name, JuMP.index(pref))
    delete!(model.param_to_group_id, JuMP.index(pref))
    return
end
=#
