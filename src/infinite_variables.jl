################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::InfiniteVariableIndex
                               )::InfiniteVariableRef
    return InfiniteVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::VariableData{<:InfiniteVariable}
                          )::InfiniteVariableIndex
    return MOIUC.add_item(model.infinite_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{InfiniteVariable}
    )::MOIUC.CleverDict{InfiniteVariableIndex, VariableData{InfiniteVariable{GeneralVariableRef}}}
    return model.infinite_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(vref::InfiniteVariableRef
    )::MOIUC.CleverDict{InfiniteVariableIndex, VariableData{InfiniteVariable{GeneralVariableRef}}}
    return JuMP.owner_model(vref).infinite_vars
end

# Extend _data_object
function _data_object(vref::InfiniteVariableRef
    )::VariableData{InfiniteVariable{GeneralVariableRef}}
  object = get(_data_dictionary(vref), JuMP.index(vref), nothing)
  object === nothing && error("Invalid infinite variable reference, cannot find " *
                        "corresponding variable in the model. This is likely " *
                        "caused by using the reference of a deleted variable.")
  return object
end

# Extend _core_variable_object
function _core_variable_object(vref::InfiniteVariableRef
    )::InfiniteVariable{GeneralVariableRef}
    return _data_object(vref).variable
end

# Extend _object_numbers
function _object_numbers(vref::InfiniteVariableRef)::Vector{Int}
    return _core_variable_object(vref).object_nums
end

# Extend _parameter_numbers
function _parameter_numbers(vref::InfiniteVariableRef)::Vector{Int}
    return _core_variable_object(vref).parameter_nums
end

# Define getter function for var.is_vector_start
function _is_vector_start(vref::InfiniteVariableRef)::Bool
    return _core_variable_object(vref).is_vector_start
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
## Check that each parameter tuple element is formatted correctly
# IndependentParameterRefs
function _check_tuple_element(_error::Function,
                              prefs::Vector{IndependentParameterRef})::Nothing
    return
end

# DependentParameterRefs
function _check_tuple_element(_error::Function,
                              prefs::Vector{DependentParameterRef})::Nothing
    if length(prefs) != _num_parameters(first(prefs))
        _error("Infinite variables cannot depend on a subset of dependent " *
               "parameters.")
    end
    return
end

# Fallback
function _check_tuple_element(_error::Function, prefs)
    _error("Cannot have mixed parameter types in a tuple element and can only " *
           "specify infinite parameters.")
end

# Check parameter tuple, ensure all elements contain parameter references
function _check_parameter_tuple(_error::Function,
                                raw_prefs::VectorTuple{GeneralVariableRef}
                                )::Nothing
    allunique(raw_prefs) || _error("Cannot double specify infinite parameter " *
                                   "references.")
    for i in 1:size(raw_prefs, 1)
        prefs = dispatch_variable_ref.(raw_prefs[i, :])
        _check_tuple_element(_error, prefs)
    end
    return
end

## Check and format the variable info considering functional start values
# Just a number given for the start value
function _check_and_format_infinite_info(_error::Function,
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Real},
    prefs::VectorTuple
    )::Tuple{JuMP.VariableInfo{Float64, Float64, Float64, Function}, Bool}
    # prepare the start value function and return the info
    start_func = (s::Vector{<:Real}) -> info.start
    return JuMP.VariableInfo{Float64, Float64, Float64, Function}(
        info.has_lb, info.lower_bound, info.has_ub, info.upper_bound,
        info.has_fix, info.fixed_value, !isnan(info.start), start_func,
        info.binary, info.integer), true
end

# A function is given for the start value generation
function _check_and_format_infinite_info(_error::Function,
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Function},
    prefs::VectorTuple
    )::Tuple{JuMP.VariableInfo{Float64, Float64, Float64, Function}, Bool}
    # check the function properties
    if length(methods(info.start)) != 1
        _error("Start value function name is not unique.")
    elseif length(first(methods(info.start)).sig.parameters) != size(prefs, 1) + 1
        _error("Start value function must match the formatting of the infinite " *
               "parameter tuple of the infinite variable.")
    end
    # make the info and return
    return JuMP.VariableInfo{Float64, Float64, Float64, Function}(
        info.has_lb, info.lower_bound, info.has_ub, info.upper_bound,
        info.has_fix, info.fixed_value, true, info.start,
        info.binary, info.integer), false
end

# Fallback
function _check_and_format_infinite_info(_error::Function,
    info::JuMP.VariableInfo,
    prefs::VectorTuple)
    _error("Unrecognized formatting for the variable information.")
end

# Define _make_variable (used by JuMP.build_variable)
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Infinite};
                        parameter_refs::Union{GeneralVariableRef,
                                              AbstractArray{<:GeneralVariableRef},
                                              Tuple, Nothing} = nothing,
                        extra_kw_args...)::InfiniteVariable{GeneralVariableRef}
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with infinite variables.")
    end
    # check that we have been given parameter references
    if parameter_refs === nothing
        _error("Parameter references not specified, use the var(params...) " *
               "syntax or the parameter_refs keyword argument.")
    end
    # format parameter_refs into a VectorTuple
    prefs = VectorTuple(parameter_refs)
    # check the VectorTuple for validity and format
    _check_parameter_tuple(_error, prefs)
    # check and format the info (accounting for start value functions)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # get the parameter object numbers
    object_nums = _object_numbers(parameter_list(prefs))
    # make the variable and return
    return InfiniteVariable(new_info, prefs,
                            [_parameter_number(pref) for pref in prefs],
                            object_nums, is_vect_func)
end

# check the pref tuple contains only valid parameters
function _check_parameters_valid(model::InfiniteModel,
                                 prefs::VectorTuple)::Nothing
    for pref in prefs
        JuMP.check_belongs_to_model(pref, model)
    end
    return
end

# Used to update parameter-infinite variable mappings
function _update_param_var_mapping(vref::InfiniteVariableRef,
                                   prefs::VectorTuple)::Nothing
    for pref in prefs
        dependency_list = _infinite_variable_dependencies(pref)
        if !(JuMP.index(vref) in dependency_list)
            push!(dependency_list, JuMP.index(vref))
        end
    end
    return
end

# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(model::InfiniteModel,
                                      v::InfiniteVariable,
                                      name::String)::InfiniteVariableRef
    _check_parameters_valid(model, v.parameter_refs)
    data_object = VariableData(v, name)
    vindex = _add_data_object(model, data_object)
    vref = InfiniteVariableRef(model, vindex)
    _update_param_var_mapping(vref, v.parameter_refs)
    return vref
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _reduced_variable_dependencies
function _reduced_variable_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
     )::Vector{ReducedVariableIndex}
    return _data_object(vref).reduced_var_indices
end

# Extend _point_variable_dependencies
function _point_variable_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Vector{PointVariableIndex}
    return _data_object(vref).point_var_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
     )::Vector{DerivativeIndex}
    return _data_object(vref).derivative_indices
end

"""
    used_by_reduced_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a reduced infinite variable.

**Example**
```julia-repl
julia> used_by_reduced_variable(vref)
false
```
"""
function used_by_reduced_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    return !isempty(_reduced_variable_dependencies(vref))
end

"""
    used_by_point_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a point variable.

**Example**
```julia-repl
julia> used_by_point_variable(vref)
false
```
"""
function used_by_point_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    return !isempty(_point_variable_dependencies(vref))
end

"""
    used_by_derivative(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a derivative.

**Example**
```julia-repl
julia> used_by_derivative(vref)
true
```
"""
function used_by_derivative(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    return !isempty(_derivative_dependencies(vref))
end

"""
    is_used(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia-repl
julia> is_used(vref)
false
```
"""
function is_used(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    if used_by_measure(vref) || used_by_constraint(vref)
        return true
    end
    for vindex in _point_variable_dependencies(vref)
        if is_used(PointVariableRef(JuMP.owner_model(vref), vindex))
            return true
        end
    end
    for vindex in _reduced_variable_dependencies(vref)
        if is_used(ReducedVariableRef(JuMP.owner_model(vref), vindex))
            return true
        end
    end
    for dindex in _derivative_dependencies(vref)
        if is_used(DerivativeRef(JuMP.owner_model(vref), dindex))
            return true
        end
    end
    return false
end

################################################################################
#                           PARAMETER REFERENCES
################################################################################
"""
    raw_parameter_refs(vref::InfiniteVariableRef)::VectorTuple{GeneralVariableRef}

Return the raw [`VectorTuple`](@ref) of the parameter references that `vref`
depends on. This is primarily an internal method where
[`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef))
is intended as the preferred user function.
"""
function raw_parameter_refs(vref::InfiniteVariableRef)::VectorTuple{GeneralVariableRef}
    return _core_variable_object(vref).parameter_refs
end

"""
    parameter_refs(vref::InfiniteVariableRef)::Tuple

Return the parameter references associated with the infinite variable `vref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `vref`.

**Example**
```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> parameter_refs(T)
(t,)
```
"""
function parameter_refs(vref::InfiniteVariableRef)::Tuple
    return Tuple(raw_parameter_refs(vref))
end

"""
    parameter_list(vref::InfiniteVariableRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `vref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef))
is intended as the preferred user function.
"""
function parameter_list(vref::InfiniteVariableRef)::Vector{GeneralVariableRef}
    return raw_parameter_refs(vref).values
end

# get parameter list from raw VectorTuple
function parameter_list(prefs::VectorTuple{GeneralVariableRef}
                        )::Vector{GeneralVariableRef}
    return prefs.values
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for infinite variables
function _update_variable_info(vref::InfiniteVariableRef,
                               info::JuMP.VariableInfo)::Nothing
    new_info = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                                 info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 info.has_start, info.start, info.binary, info.integer)
    prefs = raw_parameter_refs(vref)
    param_nums = _parameter_numbers(vref)
    obj_nums = _object_numbers(vref)
    is_vect_func = _is_vector_start(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums, is_vect_func)
    _set_core_variable_object(vref, new_var)
    return
end

# Specify start_value fallback for infinite variables
function JuMP.start_value(vref::Union{InfiniteVariableRef, DerivativeRef})
    error("`start_value` not defined for infinite variables, consider calling " *
          "`start_value_function` instead.")
end

# Specify set_start_value fallback for infinite variables
function JuMP.set_start_value(vref::Union{InfiniteVariableRef, DerivativeRef}, value::Real)
    error("`set_start_value` not defined for infinite variables, consider calling " *
          "`set_start_value_function` instead.")
end

"""
    start_value_function(vref::Union{InfiniteVariableRef, DerivativeRef})::Union{Nothing, Function}

Return the function that is used to generate the start values of `vref` for
particular support values. Returns `nothing` if no start behavior has been
specified.

**Example**
```julia-repl
julia> start_value_function(vref)
my_start_func
```
"""
function start_value_function(vref::Union{InfiniteVariableRef, DerivativeRef})::Union{Nothing, Function}
    if _variable_info(vref).has_start
        return _variable_info(vref).start
    else
        return
    end
end

"""
    set_start_value_function(vref::InfiniteVariableRef,
                             start::Union{Real, Function})::Nothing

Set the start value function of `vref`. If `start::Real` then a function is
generated to such that the start value will be `start` for the entire infinite
domain. If `start::Function` then this function should map to a scalar start value
given a support value arguments matching the format of the parameter elements in
`parameter_refs(vref)`.

**Example**
```julia-repl
julia> set_start_value_function(vref, 1) # all start values will be 1

julia> set_start_value_function(vref, my_func) # each value will be made via my_func
```
"""
function set_start_value_function(vref::InfiniteVariableRef,
                                  start::Union{Real, Function})::Nothing
    info = _variable_info(vref)
    set_optimizer_model_ready(JuMP.owner_model(vref), false)
    prefs = raw_parameter_refs(vref)
    temp_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 true, start, info.binary, info.integer)
    new_info, is_vect_func = _check_and_format_infinite_info(error, temp_info, prefs)
    obj_nums = _object_numbers(vref)
    param_nums = _parameter_numbers(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums, is_vect_func)
    _set_core_variable_object(vref, new_var)
    return
end

"""
    reset_start_value_function(vref::InfiniteVariableRef)::Nothing

Remove the existing start value function and return to the default. Generally,
this is triggered by deleting an infinite parameter that `vref` depends on.

**Example**
```julia-repl
julia> reset_start_value_function(vref)
```
"""
function reset_start_value_function(vref::InfiniteVariableRef)::Nothing
    info = _variable_info(vref)
    set_optimizer_model_ready(JuMP.owner_model(vref), false)
    start_func = (s::Vector{<:Real}) -> NaN
    new_info = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                                 info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 false, start_func, info.binary, info.integer)
    prefs = raw_parameter_refs(vref)
    obj_nums = _object_numbers(vref)
    param_nums = _parameter_numbers(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums, true)
    _set_core_variable_object(vref, new_var)
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::InfiniteVariableRef)::Nothing
    # remove variable info constraints associated with vref
    _delete_info_constraints(vref)
    # update parameter mapping
    all_prefs = parameter_list(vref)
    for pref in all_prefs
        filter!(e -> e != JuMP.index(vref), _infinite_variable_dependencies(pref))
    end
    model = JuMP.owner_model(vref)
    # delete associated point variables and mapping
    for index in _point_variable_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated reduced variables and mapping
    for index in _reduced_variable_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    return
end
