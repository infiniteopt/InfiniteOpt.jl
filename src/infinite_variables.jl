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
  object = _get(_data_dictionary(vref), JuMP.index(vref), nothing)
  isnothing(object) && error("Invalid infinite variable reference, cannot find " *
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
    # get the parameter object numbers
    object_nums = _object_numbers(parameter_list(prefs))
    # make the variable and return
    return InfiniteVariable(_make_float_info(info), prefs,
                            [_parameter_number(pref) for pref in prefs],
                            object_nums)
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
                                      v::InfiniteVariable)::InfiniteVariableRef
    _check_parameters_valid(model, v.parameter_refs)
    data_object = VariableData(v)
    vindex = _add_data_object(model, data_object)
    vref = InfiniteVariableRef(model, vindex)
    _update_param_var_mapping(vref, v.parameter_refs)
    return vref
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _reduced_variable_dependencies
function _reduced_variable_dependencies(vref::InfiniteVariableRef
                                        )::Vector{ReducedVariableIndex}
    return _data_object(vref).reduced_var_indices
end

# Extend _point_variable_dependencies
function _point_variable_dependencies(vref::InfiniteVariableRef
                                      )::Vector{PointVariableIndex}
    return _data_object(vref).point_var_indices
end

"""
    used_by_reduced_variable(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a reduced infinite variable.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @infinite_variable(m, vref(@infinite_parameter(m, t in [0, 1]))))
julia> used_by_reduced_variable(vref)
false
```
"""
function used_by_reduced_variable(vref::InfiniteVariableRef)::Bool
    return !isempty(_reduced_variable_dependencies(vref))
end

"""
    used_by_point_variable(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a point variable.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @infinite_variable(m, vref(@infinite_parameter(m, t in [0, 1]))))
julia> used_by_point_variable(vref)
false
```
"""
function used_by_point_variable(vref::InfiniteVariableRef)::Bool
    return !isempty(_point_variable_dependencies(vref))
end

"""
    is_used(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @infinite_variable(m, vref(@infinite_parameter(m, t in [0, 1]))))
julia> is_used(vref)
false
```
"""
function is_used(vref::InfiniteVariableRef)::Bool
    if used_by_measure(vref) || used_by_constraint(vref)
        return true
    end
    if used_by_point_variable(vref)
        for vindex in _point_variable_dependencies(vref)
            if is_used(PointVariableRef(JuMP.owner_model(vref), vindex))
                return true
            end
        end
    end
    if used_by_reduced_variable(vref)
        for vindex in _reduced_variable_dependencies(vref)
            if is_used(ReducedVariableRef(JuMP.owner_model(vref), vindex))
                return true
            end
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

# Internal function used to change the parameter reference tuple of an infinite
# variable
function _update_variable_param_refs(vref::InfiniteVariableRef,
                                     prefs::VectorTuple{GeneralVariableRef}
                                     )::Nothing
    # get basic information
    info = _variable_info(vref)
    param_nums = [_parameter_number(pref) for pref in prefs]
    # get the parameter object numbers
    object_nums = _object_numbers(parameter_list(prefs))
    new_var = InfiniteVariable(info, prefs, param_nums, object_nums)
    _set_core_variable_object(vref, new_var)
    JuMP.set_name(vref, _root_name(vref))
    if is_used(vref)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    end
    return
end

"""
    set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)::Nothing

Specify a new parameter reference tuple `prefs` for the infinite variable `vref`.
Note each element must contain a single parameter reference or an array of
parameter references. Errors if a parameter is double specified or if there are
point/reduced variables that depend on `vref`.

**Example**
```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> @infinite_parameter(model, x[1:2] in [-1, 1])
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]

julia> set_parameter_refs(T, (t, x))

julia> parameter_refs(T)
(t, [x[1], x[2]])
```
"""
function set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)::Nothing
    if used_by_point_variable(vref) || used_by_reduced_variable(vref)
        error("Cannot modify parameter dependencies if infinite variable has " *
              "dependent point variables and/or reduced infinite variables.")
    end
    new_prefs = VectorTuple(prefs)
    _check_parameter_tuple(error, new_prefs)
    _update_variable_param_refs(vref, new_prefs)
    return
end

"""
    add_parameter_ref(vref::InfiniteVariableRef,
        pref::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
        )::Nothing

Add additional parameter reference or group of parameter references to be
associated with the infinite variable `vref`. Errors if the parameter references
are already added to the variable, if the added parameters have different
or if `vref` has point/reduced variable dependencies.

```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> @infinite_parameter(model, x[1:2] in [-1, 1])
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]

julia> add_parameter_ref(T, x)

julia> name(T)
"T(t, x)"
```
"""
function add_parameter_ref(vref::InfiniteVariableRef,
                           pref::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
                           )::Nothing
    if used_by_point_variable(vref) || used_by_reduced_variable(vref)
       error("Cannot modify parameter dependencies if infinite variable has " *
             "dependent point variables and/or reduced infinite variables.")
    end
    # check that new group is unique from old ones
    prefs = raw_parameter_refs(vref)
    if !allunique([parameter_list(prefs); pref])
        error("Cannot double specify infinite parameter references.")
    end
    if pref isa GeneralVariableRef
        _check_tuple_element(error, [dispatch_variable_ref(pref)])
    else
        _check_tuple_element(error, dispatch_variable_ref.(pref))
    end
    # add the new parameter(s)
    prefs = push!(prefs, pref)
    _update_variable_param_refs(vref, prefs)
    return
end

################################################################################
#                            VARIABLE NAMING
################################################################################
"""
    JuMP.set_name(vref::InfiniteVariableRef, root_name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of infinite variables. Adds on to `root_name` the ending `(prefs...)`
where the parameter reference names are listed in the same format as input in
the parameter reference tuple.

**Example**
```julia-repl
julia> name(vref)
"vref(t)"

julia> set_name(vref, "new_name")

julia> name(vref)
"new_name(t)"
```
"""
function JuMP.set_name(vref::InfiniteVariableRef, root_name::String)::Nothing
    if length(root_name) == 0
        root_name = "noname"
    end
    prefs = raw_parameter_refs(vref)
    param_name_tuple = "("
    for i in 1:size(prefs, 1)
        element_prefs = prefs[i, :]
        type = _index_type(first(prefs))
        if type == DependentParameterIndex
            param_name = _remove_name_index(first(element_prefs))
        elseif length(element_prefs) == 1
            param_name = JuMP.name(first(element_prefs))
        else
            names = _remove_name_index.(element_prefs)
            if _allequal(names)
                param_name = first(names)
            else
                param_name = string("[", join(element_prefs, ", "), "]")
            end
        end
        if i != size(prefs, 1)
            param_name_tuple *= string(param_name, ", ")
        else
            param_name_tuple *= string(param_name, ")")
        end
    end
    var_name = string(root_name, param_name_tuple)
    _data_object(vref).name = var_name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for infinite variables
function _update_variable_info(vref::InfiniteVariableRef,
                               info::JuMP.VariableInfo)::Nothing
    new_info = _make_float_info(info)
    prefs = raw_parameter_refs(vref)
    param_nums = _parameter_numbers(vref)
    obj_nums = _object_numbers(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums)
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
