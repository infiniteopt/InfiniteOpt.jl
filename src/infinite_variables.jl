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

# Extend _data_dictionary
function _data_dictionary(vref::InfiniteVariableRef)::MOIUC.CleverDict
    return JuMP.owner_model(vref).infinite_vars
end

# Extend _data_object
function _data_object(vref::InfiniteVariableRef)::VariableData{<:InfiniteVariable}
    return _data_dictionary(vref)[JuMP.index(vref)]
end

# Extend _core_variable_object
function _core_variable_object(vref::InfiniteVariableRef)::InfiniteVariable
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
# Check parameter tuple, ensure all elements contain parameter references
function _check_parameter_tuple(_error::Function, prefs::VectorTuple)::Nothing
    if !all(pref.index_type <: InfiniteParameterIndex for pref in prefs)
        _error("Invalid parameter type(s) given.")
    end
    allunique(prefs) ||  _error("Cannot double specify infinite parameter " *
                                "references.")
    return
end

# Define _make_variable (used by JuMP.build_variable)
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Infinite};
                        parameter_refs::Union{GeneralVariableRef,
                                              AbstractArray{<:GeneralVariableRef},
                                              Tuple, Nothing} = nothing,
                        extra_kw_args...)::InfiniteVariable
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
    # make the variable and return
    return InfiniteVariable(info, prefs,
                            [_parameter_number(pref) for pref in prefs],
                            unique!([_object_number(pref) for pref in prefs]))
end

# check the pref tuple contains only valid parameters
function _check_parameters_valid(model::InfiniteModel,
                                 prefs::VectorTuple)::Nothing
    for pref in prefs
        JuMP.is_valid(model, pref) || error("Invalid Parameter reference " *
                                            "provided.")
    end
    return
end

# Used to update parameter-infinite variable mappings
function _update_param_var_mapping(vref::InfiniteVariableRef,
                                   prefs::VectorTuple)::Nothing
    for pref in prefs
        push!(_infinite_variable_dependencies(pref), JuMP.index(vref))
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
                                        )::Vector{ReducedInfiniteVariableIndex}
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
            rvref = ReducedInfiniteVariableRef(JuMP.owner_model(vref), vindex)
            if is_used(rvref)
                return true
            end
        end
    end
    return false
end

################################################################################
#                            VARIABLE NAMING
################################################################################
#=
"""
    JuMP.set_name(vref::InfiniteVariableRef, root_name::String)

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of infinite variables. Adds on to `root_name` the ending `(prefs...)`
where the parameter reference names are listed in the same format as input in
the parameter reference tuple.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @infinite_variable(m, vref(@infinite_parameter(m, t in [0, 1]))))
julia> name(vref)
"vref(t)"

julia> set_name(vref, "new_name")

julia> name(vref)
"new_name(t)"
```
"""
function JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
    if length(root_name) == 0
        root_name = "noname"
    end
    prefs = raw_parameter_refs(vref)
    # TODO need to change this paradigm completely
    param_names = [_root_name(pref) for pref in prefs[:, 1]]
    param_name_tuple = "("
    for i in eachindex(param_names)
        if i != length(param_names)
            param_name_tuple *= string(param_names[i], ", ")
        else
            param_name_tuple *= string(param_names[i])
        end
    end
    param_name_tuple *= ")"
    var_name = string(root_name, param_name_tuple)
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = var_name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    raw_parameter_refs(vref::InfiniteVariableRef)::VectorTuple{ParameterRef}

Return the raw [`VectorTuple`](@ref) of the parameter references that `vref`
depends on. This is primarily an internal method where
[`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef))
is intended as the preferred user function.
"""
function raw_parameter_refs(vref::InfiniteVariableRef)::VectorTuple{ParameterRef}
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_refs
end

"""
    parameter_refs(vref::InfiniteVariableRef)::Tuple

Return the `ParameterRef`(s) associated with the infinite variable `vref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
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
    parameter_list(vref::InfiniteVariableRef)::Vector{ParameterRef}

Return a vector of the parameter references that `vref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef))
is intended as the preferred user function.
"""
function parameter_list(vref::InfiniteVariableRef)::Vector{ParameterRef}
    return raw_parameter_refs(vref).values
end

# get parameter list from raw VectorTuple
function parameter_list(prefs::VectorTuple{ParameterRef})::Vector{ParameterRef}
    return prefs.values
end

# Internal function used to change the parameter reference tuple of an infinite
# variable
function _update_variable_param_refs(vref::InfiniteVariableRef, prefs::VectorTuple)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, prefs)
    return
end

"""
    set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)

Specify a new parameter reference tuple `prefs` for the infinite variable `vref`.
Note each element must contain a single parameter reference or an array of
parameter references. Errors if a parameter is double specified, if an element
contains parameters with different group IDs, or if there are point variables
that depend on `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> @infinite_variable(model, T(t))
T(t)

julia> @infinite_parameter(model, x[1:2] in [-1, 1])
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> set_parameter_refs(T, (t, x))

julia> parameter_refs(T)
(t, [x[1], x[2]])
```
"""
function set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)
    if used_by_point_variable(vref) || used_by_reduced_variable(vref)
        error("Cannot modify parameter dependencies if infinite variable has " *
              "dependent point variables and/or reduced infinite variables.")
    end
    prefs = VectorTuple(prefs)
    _check_parameter_tuple(error, prefs)
    _check_tuple_groups(error, prefs)
    _update_variable_param_refs(vref, prefs)
    JuMP.set_name(vref, _root_name(vref))
    if is_used(vref)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    end
    return
end

"""
    add_parameter_ref(vref::InfiniteVariableRef,
                      pref::Union{ParameterRef, AbstractArray{<:ParameterRef}})

Add additional parameter reference or group of parameter references to be
associated with the infinite variable `vref`. Errors if the parameter references
are already added to the variable, if the added parameters have different
group IDs, or if `vref` has point variable dependencies.

```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> @infinite_variable(model, T(t))
T(t)

julia> @infinite_parameter(model, x[1:2] in [-1, 1])
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> add_parameter_ref(T, x)

julia> name(T)
"T(t, x)"
```
"""
function add_parameter_ref(vref::InfiniteVariableRef,
                       pref::Union{ParameterRef, AbstractArray{<:ParameterRef}})
    if used_by_point_variable(vref) || used_by_reduced_variable(vref)
       error("Cannot modify parameter dependencies if infinite variable has " *
             "dependent point variables and/or reduced infinite variables.")
    end
    # check that array only contains one group ID
    if pref isa AbstractArray
        first_group = group_id(first(pref))
        for i in 2:length(pref)
            if group_id(pref[i]) == first_group
                error("Each parameter tuple element must have contain only " *
                      "infinite parameters with the same group ID.")
            end
        end
    end
    # check that new group is unique from old ones
    prefs = raw_parameter_refs(vref)
    if !allunique(group_id.([prefs[:, 1]; first(pref)]))
        error("Cannot double specify infinite parameter references.")
    end
    # add the new parameter(s)
    prefs = push!(prefs, pref)
    _update_variable_param_refs(vref, prefs)
    JuMP.set_name(vref, _root_name(vref))
    if is_used(vref)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    end
    return
end
=#
