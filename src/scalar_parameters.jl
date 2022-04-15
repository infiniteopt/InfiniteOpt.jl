################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::IndependentParameterIndex
    )
    return IndependentParameterRef(model, index)
end
function dispatch_variable_ref(
    model::InfiniteModel,
    index::FiniteParameterIndex
    )
    return FiniteParameterRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::ScalarParameterData{<:IndependentParameter}
    )
    index =  MOIUC.add_item(model.independent_params, object)
    push!(model.param_object_indices, index)
    return index
end
function _add_data_object(
    model::InfiniteModel,
    object::ScalarParameterData{<:FiniteParameter}
    )
    return MOIUC.add_item(model.finite_params, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel,
    ::Type{IndependentParameter}
    )
    return model.independent_params
end
function _data_dictionary(
    model::InfiniteModel,
    ::Type{FiniteParameter}
    )
    return model.finite_params
end

# Extend _data_dictionary (ref based)
function _data_dictionary(pref::IndependentParameterRef)
    return JuMP.owner_model(pref).independent_params
end
function _data_dictionary(pref::FiniteParameterRef)
    return JuMP.owner_model(pref).finite_params
end

# Extend _data_object
function _data_object(pref::ScalarParameterRef)
    object = Base.get(_data_dictionary(pref), JuMP.index(pref), nothing)
    if isnothing(object)
        error("Invalid scalar parameter reference, cannot find ",
              "corresponding parameter in the model. This is likely ",
              "caused by using the reference of a deleted parameter.")
    end
    return object
end

################################################################################
#                             CORE OBJECT METHODS
################################################################################
# Extend _core_variable_object for IndependentParameterRefs
function _core_variable_object(pref::IndependentParameterRef)
    return _data_object(pref).parameter
end

# Extend _core_variable_object for FiniteParameterRefs
function _core_variable_object(pref::FiniteParameterRef)::FiniteParameter
    return _data_object(pref).parameter
end

# Extend _parameter_number
function _parameter_number(pref::IndependentParameterRef)
    return _data_object(pref).parameter_num
end

# Extend _parameter_numbers
function _parameter_numbers(pref::IndependentParameterRef)
    return [_parameter_number(pref)]
end

# Extend _object_number
function _object_number(pref::IndependentParameterRef)
    return _data_object(pref).object_num
end

# Extend _object_numbers
function _object_numbers(pref::IndependentParameterRef)
    return [_object_number(pref)]
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(
    pref::ScalarParameterRef, 
    param::P, 
    data::ScalarParameterData{P}
    ) where {P <: ScalarParameter}
    data.parameter = param
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(
    pref::ScalarParameterRef, 
    param::P1, 
    data::ScalarParameterData{P2}
    )  where {P1, P2}
    new_data = ScalarParameterData(param, data.object_num, data.parameter_num, 
                                   data.name, data.parameter_func_indices,
                                   data.infinite_var_indices, 
                                   data.derivative_indices, data.measure_indices,
                                   data.constraint_indices, data.in_objective)
    _data_dictionary(pref)[JuMP.index(pref)] = new_data
    return
end

# Extend _set_core_variable_object for ScalarParameterRefs
function _set_core_variable_object(
    pref::ScalarParameterRef,
    param::ScalarParameter
    )
    _adaptive_data_update(pref, param, _data_object(pref))
    return
end

################################################################################
#                       TRANSFORM ATTRIBUTE KEYWORD SUPPORT
################################################################################
# Store the keywords to be accepted by @infinite_parameter
const _InfiniteParameterKeywords = Dict{Symbol, Any}()

# Add registeration method for @infinite_parameter
function register_transform_keyword(
    kw::Symbol, 
    attr::InfiniteParameterAttr, 
    type = nothing
    )
    return _add_transform_keyword(_InfiniteParameterKeywords, kw, attr, type)
end

# Store the keywords to be accepted by @finite_parameter
const _FiniteParameterKeywords = Dict{Symbol, Any}()

# Add registeration method for @finite_parameter
function register_transform_keyword(
    kw::Symbol, 
    attr::FiniteParameterAttr, 
    type = nothing
    )
    return _add_transform_keyword(_FiniteParameterKeywords, kw, attr, type)
end

################################################################################
#                            PARAMETER DEFINITION
################################################################################
"""
    build_parameter(_error::Function, domain::InfiniteScalarDomain; [kwargs...])::IndependentParameter

Returns a [`IndependentParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. This is meant to primarily serve as a
helper method for [`@infinite_parameter`](@ref).

**Example**
```julia-repl
julia> param = build_parameter(error, IntervalDomain(0, 3));
```
"""
function build_parameter(_error::Function, domain::InfiniteScalarDomain; kwargs...)
    return _process_transform_kwargs(_error, _InfiniteParameterKeywords, kwargs, 
                                     IndependentParameter(domain))
end

# Fallback for bad domain types 
function build_parameter(_error::Function, domain::AbstractInfiniteDomain)
    _error("Expected scalar infinite domain for each independent parameter, ",
           "but got a domain of type `$(domain)`. If you are trying to use an ",
           "`InfiniteArrayDomain`, try setting `independent = false`.")
end

"""
    build_parameter(_error::Function, value::Real; [kwargs...])::FiniteParameter

Returns a [`FiniteParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. This is meant to primarily serve as
a helper method for [`@finite_parameter`](@ref).

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> build_finite_parameter(error, 1)
FiniteParameter(1.0)
```
"""
function build_parameter(_error::Function, value::Real; kwargs...)
    return _process_transform_kwargs(_error, _FiniteParameterKeywords, kwargs, 
                                     FiniteParameter(value))
end

# Generic fallback
function build_parameter(_error::Function, arg, kwargs...)
    _error("Unexpected input given for the value of the parameter.")
end

"""
    add_parameter(model::InfiniteModel, p::IndependentParameter,
                  [name::String = ""])::GeneralVariableRef

Returns a [`GeneralVariableRef`](@ref) associated with the parameter `p` that is added 
to `model`. This adds a parameter to the model in a manner similar to 
`JuMP.add_variable`. This is used to add parameters with the use of 
[`@infinite_parameter`](@ref). 
[`build_parameter`](@ref build_parameter(::Function, ::InfiniteScalarDomain)) 
should be used to construct `p`.

**Example**
```julia-repl
julia> p = build_parameter(error, IntervalDomain(0, 3), supports = Vector(0:3));

julia> param_ref = add_parameter(model, p, "name")
name
```
"""
function add_parameter(
    model::InfiniteModel, 
    p::IndependentParameter,
    name::String = ""
    )
    obj_num = length(_param_object_indices(model)) + 1
    param_num = model.last_param_num += 1
    data_object = ScalarParameterData(p, obj_num, param_num, name)
    obj_index = _add_data_object(model, data_object)
    model.name_to_param = nothing
    _update_transform_attributes(model, p)
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

"""
    add_parameter(model::InfiniteModel, p::FiniteParameter,
                  [name::String = ""])::GeneralVariableRef

Returns a [`GeneralVariableRef`](@ref) associated with the parameter `p` that is 
added to `model`. This adds a parameter to the model in a manner similar to
`JuMP.add_variable`. This is to add parameters with the use of 
[`@finite_parameter`](@ref). 
[`build_parameter`](@ref build_parameter(::Function, ::Real)) should be used to
construct `p`.

**Example**
```julia-repl
julia> p = build_parameter(error, 42);

julia> param_ref = add_parameter(model, p, "name")
name
```
"""
function add_parameter(
    model::InfiniteModel, 
    p::FiniteParameter,
    name::String = ""
    )
    data_object = ScalarParameterData(p, -1, -1, name)
    obj_index = _add_data_object(model, data_object)
    model.name_to_param = nothing
    _update_transform_attributes(model, p)
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

"""
    add_parameter(model::InfiniteModel, obj::ObjectWithAttributes{<:ScalarParameter},
                  name::String = "")::GeneralVariableRef

Add a parameter build with `build_parameter` that contains [`InfiniteParameterAttr`](@ref)s 
that need to be added to `model` as well. 
"""
function add_parameter(
    model::InfiniteModel, 
    obj::ObjectWithAttributes{<:ScalarParameter},
    name::String = ""
    )
    pref = add_parameter(model, obj.object, name)
    idx = JuMP.index(pref)
    for (a, v) in obj.attributes
        InfiniteOpt.set(model, idx, a, v)
    end
    return pref
end

################################################################################
#                           PARAMETER DEPENDENCIES
################################################################################
# Extend _infinite_variable_dependencies
function _infinite_variable_dependencies(pref::ScalarParameterRef
    )::Vector{InfiniteVariableIndex}
    return _data_object(pref).infinite_var_indices
end

# Extend _parameter_function_dependencies
function _parameter_function_dependencies(pref::ScalarParameterRef
    )::Vector{ParameterFunctionIndex}
    return _data_object(pref).parameter_func_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(pref::ScalarParameterRef
    )::Vector{DerivativeIndex}
    return _data_object(pref).derivative_indices
end

# Extend _measure_dependencies
function _measure_dependencies(pref::ScalarParameterRef
    )::Vector{MeasureIndex}
    return _data_object(pref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(pref::ScalarParameterRef
    )::Vector{InfOptConstraintIndex}
    return _data_object(pref).constraint_indices
end

################################################################################
#                             USED_BY FUNCTIONS
################################################################################
"""
    used_by_infinite_variable(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by an infinite variable or false otherwise.

**Example**
```julia-repl
julia> used_by_infinite_variable(t)
true
```
"""
function used_by_infinite_variable(pref::IndependentParameterRef)::Bool
    return !isempty(_infinite_variable_dependencies(pref))
end

# FiniteParameter
used_by_infinite_variable(pref::FiniteParameterRef)::Bool = false

"""
    used_by_parameter_function(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by an infinite parameter function or false otherwise.

**Example**
```julia-repl
julia> used_by_parameter_function(t)
false
```
"""
function used_by_parameter_function(pref::IndependentParameterRef)::Bool
    return !isempty(_parameter_function_dependencies(pref))
end

# FiniteParameter
used_by_parameter_function(pref::FiniteParameterRef)::Bool = false

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
    return !isempty(_measure_dependencies(pref))
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
    return !isempty(_constraint_dependencies(pref))
end

"""
    used_by_objective(pref::FiniteParameterRef)::Bool

Return true if `pref` is used by the objective function.

**Example**
```julia-repl
```
"""
function used_by_objective(pref::FiniteParameterRef)::Bool
    return _data_object(pref).in_objective
end

# IndependentParameter
used_by_objective(::IndependentParameterRef)::Bool = false

"""
    used_by_derivative(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by a derivative or false otherwise.

**Example**
```julia-repl
julia> used_by_derivative(t)
false
```
"""
function used_by_derivative(pref::IndependentParameterRef)::Bool
    return !isempty(_derivative_dependencies(pref))
end

# FiniteParameter
used_by_derivative(::FiniteParameterRef)::Bool = false

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
    return used_by_measure(pref) || used_by_constraint(pref) ||
           used_by_infinite_variable(pref) || used_by_objective(pref) || 
           used_by_derivative(pref) || used_by_parameter_function(pref)
end

################################################################################
#                              NAME METHODS
################################################################################

"""
    JuMP.name(pref::Union{IndependentParameterRef, FiniteParameterRef})::String

Extend the `JuMP.name` function to accomodate infinite parameters. Returns the 
name string associated with `pref`.

**Example**
```julia-repl
julia> name(t)
"t"
```
"""
function JuMP.name(pref::ScalarParameterRef)::String
    object = Base.get(_data_dictionary(pref), JuMP.index(pref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(pref::ScalarParameterRef, name::String)

Extend the `JuMP.set_name` function to accomodate infinite parameters. Set a new 
base name to be associated with `pref`.

**Example**
```julia-repl
julia> set_name(t, "time")

julia> name(t)
"time"
```
"""
function JuMP.set_name(pref::ScalarParameterRef, name::String)
    _data_object(pref).name = name
    JuMP.owner_model(pref).name_to_param = nothing
    return
end

# Make a parameter reference
function _make_parameter_ref(model::InfiniteModel, index::AbstractInfOptIndex)
    return GeneralVariableRef(model, MOIUC.key_to_index(index), typeof(index))
end
function _make_parameter_ref(model::InfiniteModel, index::DependentParameterIndex)
    return GeneralVariableRef(model, MOIUC.key_to_index(index.object_index),
                              typeof(index), index.param_index)
end

# Get the name_to_param Dictionary
function _param_name_dict(model::InfiniteModel)
    return model.name_to_param
end

# Update name_to_param
function _update_param_name_dict(
    model::InfiniteModel,
    param_dict::MOIUC.CleverDict{K, V}
    ) where {K, V <: ScalarParameterData}
    name_dict = _param_name_dict(model)
    for (index, data_object) in param_dict
        param_name = data_object.name
        if haskey(name_dict, param_name)
            # IndependentParameterIndex(-1) is a special value that means
            # this string does not map to a unique variable name.
            name_dict[param_name] = IndependentParameterIndex(-1)
        else
            name_dict[param_name] = index
        end
    end
    model.name_to_param = name_dict
    return
end
function _update_param_name_dict(model::InfiniteModel,
    param_dict::MOIUC.CleverDict{K, V}
    ) where {K, V <: MultiParameterData}
    name_dict = _param_name_dict(model)
    for (index, data_object) in param_dict
        param_nums = data_object.parameter_nums
        for i in eachindex(param_nums)
            name = data_object.names[i]
            if haskey(name_dict, name)
                # IndependentParameterIndex(-1) is a special value that means
                # this string does not map to a unique variable name.
                name_dict[name] = IndependentParameterIndex(-1)
            else
                individual_index = DependentParameterIndex(index, i)
                name_dict[name] = individual_index
            end
         end
    end
    model.name_to_param = name_dict
    return
end

"""
    parameter_by_name(model::InfiniteModel,
                      name::String)::Union{GeneralVariableRef, Nothing}

Return the parameter reference assoociated with a parameter name. Errors if
multiple parameters have the same name. Returns nothing if no such name exists.

**Example**
```julia-repl
julia> parameter_by_name(model, "t")
t
```
"""
function parameter_by_name(model::InfiniteModel, name::String)
    if isnothing(_param_name_dict(model))
        model.name_to_param = Dict{String, AbstractInfOptIndex}()
        _update_param_name_dict(model, model.independent_params)
        _update_param_name_dict(model, model.dependent_params)
        _update_param_name_dict(model, model.finite_params)
    end
    index = Base.get(_param_name_dict(model), name, nothing)
    if isnothing(index)
        return nothing
    elseif index == IndependentParameterIndex(-1)
        error("Multiple parameters have the name $name.")
    else
        return _make_parameter_ref(model, index)
    end
end

################################################################################
#                              DOMAIN FUNCTIONS
################################################################################
# Internal functions
function _parameter_domain(pref::IndependentParameterRef)
    return _core_variable_object(pref).domain
end
function _update_parameter_domain(
    pref::IndependentParameterRef,
    domain::AbstractInfiniteDomain
    )
    new_param = IndependentParameter(domain)
    _set_core_variable_object(pref, new_param)
    # TODO do something about the supports and other attributes
    if is_used(pref)
        set_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    infinite_domain(pref::IndependentParameterRef)::InfiniteScalarDomain

Return the infinite domain associated with `pref`.

**Example**
```julia-repl
julia> infinite_domain(t)
[0, 1]
```
"""
function infinite_domain(pref::IndependentParameterRef)
    return _parameter_domain(pref)
end

"""
    set_infinite_domain(pref::IndependentParameterRef,
                        domain::InfiniteScalarDomain)::Nothing

Reset the infinite domain of `pref` with another `InfiniteScalarDomain`. An error will 
be thrown if `pref` is being used by some measure.

**Example**
```julia-repl
julia> set_infinite_domain(t, IntervalDomain(0, 2))

julia> infinite_domain(t)
[0, 2]
```
"""
function set_infinite_domain(
    pref::IndependentParameterRef,
    domain::InfiniteScalarDomain
    )
    if used_by_measure(pref)
        error("$pref is used by a measure so changing its " *
              "infinite domain is not allowed.")
    end
    _update_parameter_domain(pref, domain)
    return
end

################################################################################
#                        LOWER/UPPER BOUND FUNCTIONS
################################################################################
"""
    JuMP.has_lower_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_lower_bound` function to accomodate infinite parameters.
Return true if the domain associated with `pref` has a defined lower bound or if a
lower bound can be found. Extensions with user-defined infinite domain types
should extend `JuMP.has_lower_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_lower_bound(t)
true
```
"""
function JuMP.has_lower_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    return JuMP.has_lower_bound(domain)
end

"""
    JuMP.lower_bound(pref::IndependentParameterRef)::Real

Extend the `JuMP.lower_bound` function to accomodate infinite parameters.
Returns the lower bound associated with the infinite domain. Errors if such a bound
is not well-defined.

**Example**
```julia-repl
julia> lower_bound(t)
0.0
```
"""
function JuMP.lower_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    return JuMP.lower_bound(domain)
end

"""
    JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_lower_bound` function to accomodate infinite parameters.
Updates the infinite domain lower bound if such an operation is supported. Set
extensions that seek to employ this should extend
`JuMP.set_lower_bound(domain::NewType, lower::Number)`.

**Example**
```julia-repl
julia> set_lower_bound(t, -1)

julia> lower_bound(t)
-1.0
```
"""
function JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Real)
    domain = _parameter_domain(pref)
    new_domain = JuMP.set_lower_bound(domain, lower)
    _update_parameter_domain(pref, new_domain)
    return
end

"""
    JuMP.has_upper_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_upper_bound` function to accomodate infinite parameters.
Return true if the domain associated with `pref` has a defined upper bound or if a
upper bound can be found. Extensions with user-defined domains should extend
`JuMP.has_upper_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_upper_bound(t)
true
```
"""
function JuMP.has_upper_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    return JuMP.has_upper_bound(domain)
end

"""
    JuMP.upper_bound(pref::IndependentParameterRef)::Real

Extend the `JuMP.upper_bound` function to accomodate infinite parameters.
Returns the upper bound associated with the infinite domain. Errors if such a bound
is not well-defined. Extensions with user-defined domain types should extend
`JuMP.has_upper_bound(domain::NewType)` and `JuMP.upper_bound(domain::NewType)` if
appropriate.

**Example**
```julia-repl
julia> upper_bound(t)
1.0
```
"""
function JuMP.upper_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    return JuMP.upper_bound(domain)
end

"""
    JuMP.set_upper_bound(pref::IndependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_upper_bound` function to accomodate infinite parameters.
Updates the infinite domain upper bound if and only if it is an IntervalDomain. Errors
otherwise. Extensions with user-defined infinite domains should extend
`JuMP.set_upper_bound(domain::NewType, upper::Number)` if appropriate.

**Example**
```julia-repl
julia> set_upper_bound(t, 2)

julia> upper_bound(t)
2.0
```
"""
function JuMP.set_upper_bound(pref::IndependentParameterRef, upper::Real)
    domain = _parameter_domain(pref)
    new_domain = JuMP.set_upper_bound(domain, upper)
    _update_parameter_domain(pref, new_domain)
    return
end

"""
    parameter_value(pref::FiniteParameterRef)::Float64

Return the value of a finite parameter reference `pref`. Errors if it is
an infinite parameter.

**Example**
```julia-repl
julia> value(cost)
42.0
```
"""
function parameter_value(pref::FiniteParameterRef)
    return _core_variable_object(pref).value
end

"""
    JuMP.set_value(pref::FiniteParameterRef, value::Real)::Nothing

Set the value of `pref` so long as it is a finite parameter. Errors if it is
an infinite parameter.

**Example**
```julia-repl
julia> set_value(cost, 27)

julia> value(cost)
27.0
```
"""
function JuMP.set_value(pref::FiniteParameterRef, value::Real)
    _data_object(pref).parameter = FiniteParameter(value)
    if is_used(pref)
        set_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

################################################################################
#                               DELETE FUNCTIONS
################################################################################
# Check if parameter is used by measure data and error if it is to prevent bad
# deleting behavior
function _check_param_in_data(pref::GeneralVariableRef, data::AbstractMeasureData)
    prefs = parameter_refs(data)
    if isequal(pref, prefs) || any(isequal(pref), prefs)
        error("Unable to delete `$pref` since it is used to evaluate measures.")
    end
    return
end

# Update the dependent measures
function _update_measures(model::InfiniteModel, pref::GeneralVariableRef)
    for mindex in _measure_dependencies(pref)
        mref = dispatch_variable_ref(model, mindex)
        func = measure_function(mref)
        if func isa GeneralVariableRef
            data = measure_data(mref)
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_meas = Measure(new_func, data, Int[], Int[], true)
            _set_core_variable_object(mref, new_meas)
        else
            _remove_variable(func, pref)
        end
    end
    return
end

# Update the dependent constraints
function _update_constraints(model::InfiniteModel, pref::GeneralVariableRef)
    for cindex in copy(_constraint_dependencies(pref))
        cref = _make_constraint_ref(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_constraint_object(cref, new_constr)
            empty!(_object_numbers(cref))
        elseif func isa AbstractArray && any(isequal(pref), func)
            JuMP.delete(model, cref)
        else
            _remove_variable(func, pref)
        end
    end
    return
end

# Remove given object/parameter number and update the list
function _update_number_list(nums::Vector{Int}, list::Vector{Int})
    filter!(e -> !(e in nums), list)
    max_num = maximum(nums)
    for i in eachindex(list)
        if list[i] > max_num
            list[i] -= length(nums)
        end
    end
    return
end

# Update the model with the removed parameter/object numbers
function _update_model_numbers(
    model::InfiniteModel, 
    obj_num::Int,
    param_nums::Vector{Int}
    )
    # update the independent parameters
    for (_, object) in _data_dictionary(model, IndependentParameter)
        if object.object_num > obj_num
            object.object_num -= 1
            object.parameter_num -= length(param_nums)
        end
    end
    # update the dependent parameters
    for (_, object) in _data_dictionary(model, DependentParameters)
        if object.object_num > obj_num
            object.object_num -= 1
            object.parameter_nums = object.parameter_nums .- length(param_nums)
        end
    end
    # update the infinite parameter functions
    for (_, object) in model.param_functions
        _update_number_list([obj_num], object.func.object_nums)
        _update_number_list(param_nums, object.func.parameter_nums)
    end
    # update the infinite variables
    for vref in JuMP.all_variables(model, InfiniteVariable)
        _update_number_list([obj_num], _object_numbers(vref))
        _update_number_list(param_nums, _parameter_numbers(vref))
    end
    # update the semi-infinite variables
    for vref in JuMP.all_variables(model, SemiInfiniteVariable)
        _update_number_list([obj_num], _object_numbers(vref))
        _update_number_list(param_nums, _parameter_numbers(vref))
    end
    # update the measures
    for mref in all_measures(model)
        _update_number_list([obj_num], _object_numbers(mref))
        _update_number_list(param_nums, _parameter_numbers(mref))
    end
    # update the constraints
    for (_, object) in model.constraints
        _update_number_list([obj_num], object.object_nums)
    end
    # update the central info
    deleteat!(_param_object_indices(model), obj_num)
    model.last_param_num -= length(param_nums)
    return
end

"""
    JuMP.delete(model::InfiniteModel, pref::ScalarParameterRef)::Nothing

Extend `JuMP.delete` to delete
scalar parameters and their dependencies. All variables, constraints, and
measure functions that depend on `pref` are updated to exclude it. Errors if the
parameter is used by an infinite variable or if it is contained in an 
`AbstractMeasureData` DataType that is employed by
a measure since the measure becomes invalid otherwise. Thus, measures that
contain this dependency must be deleted first. Note that
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) needs to be
extended to allow deletion of parameters when custom `AbstractMeasureData`
datatypes are used.

**Example**
```julia-repl
julia> delete(model, x)
```
"""
function JuMP.delete(
    model::InfiniteModel, 
    pref::IndependentParameterRef
    )
    @assert JuMP.is_valid(model, pref) "Parameter reference is invalid."
    gvref = _make_parameter_ref(JuMP.owner_model(pref), JuMP.index(pref))
    # ensure deletion is okay (pref isn't used by measure data)
    for mindex in _measure_dependencies(pref)
        data = measure_data(dispatch_variable_ref(model, mindex))
        _check_param_in_data(gvref, data)
    end
    # ensure pref is not used by an infinite variable
    if used_by_infinite_variable(pref)
        error("Cannot delete `$pref` since it is used by an infinite ",
              "variable(s).")
    end
    # ensure pref is not used by a parameter function 
    if used_by_parameter_function(pref)
        error("Cannot delete `$pref` since it is used by a parameter ",
              "function(s).")
    end
    # update backend status
    if is_used(pref)
        set_backend_ready(model, false)
    end
    # delete dependence of measures on pref
    _update_measures(model, gvref)
    # delete any derivatives that use pref 
    for index in _derivative_dependencies(pref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # update constraints in mapping to remove the parameter
    _update_constraints(model, gvref)
    # delete parameter information stored in model
    obj_num = _object_number(pref)
    param_nums = _parameter_numbers(pref)
    _delete_data_object(pref)
    # update the object numbers and parameter numbers
    _update_model_numbers(model, obj_num, param_nums)
    # delete any backend information
    delete!(model.transform_backend, JuMP.index(pref))
    return
end

# FiniteParameterRef
function JuMP.delete(model::InfiniteModel, pref::FiniteParameterRef)::Nothing
    @assert JuMP.is_valid(model, pref) "Parameter reference is invalid."
    # update backend status
    if is_used(pref)
        set_backend_ready(model, false)
    end
    gvref = _make_parameter_ref(model, JuMP.index(pref))
    # delete dependence of measures on pref
    _update_measures(model, gvref)
    # update constraints in mapping to remove the parameter
    _update_constraints(model, gvref)
    # update the objective if necessary
    if used_by_objective(pref)
        func = JuMP.objective_function(model)
        if func isa GeneralVariableRef
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            JuMP.set_objective_function(model, new_func)
            JuMP.set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        else
            _remove_variable(func, gvref)
        end
    end
    # delete parameter information stored in model
    _delete_data_object(pref)
    # delete any backend information
    delete!(model.transform_backend, JuMP.index(pref))
    return
end
