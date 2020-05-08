################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::PointVariableIndex
                               )::PointVariableRef
    return PointVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::VariableData{<:PointVariable}
                          )::PointVariableIndex
    return MOIUC.add_item(model.point_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{PointVariable}
    )::MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}
    return model.point_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(vref::PointVariableRef
    )::MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}
    return JuMP.owner_model(vref).point_vars
end

# Extend _data_object
function _data_object(vref::PointVariableRef
    )::VariableData{PointVariable{GeneralVariableRef}}
    return _data_dictionary(vref)[JuMP.index(vref)]
end

# Extend _core_variable_object
function _core_variable_object(vref::PointVariableRef
    )::PointVariable{GeneralVariableRef}
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
# Ensure parameter values match shape of parameter reference tuple stored in the
# infinite variable reference
function _check_tuple_shape(_error::Function,
                            ivref::InfiniteVariableRef,
                            values::VectorTuple)::Nothing
    prefs = raw_parameter_refs(ivref)
    if !same_structure(prefs, values)
        _error("The dimensions and array formatting of the infinite parameter " *
               "values must match those of the parameter references for the " *
               "infinite variable.")
    end
    return
end

## Dispatch methods for checking the supports of parameters in point variable
# IndependentParameterRefs
function _check_element_support(_error::Function, prefs::Vector{IndependentParameterRef},
                                param_values::Vector{Float64},
                                counter::Int)::Int
    for pref in prefs
        if !supports_in_set(param_values[counter], infinite_set(pref))
            _error("Parameter values violate parameter bounds.")
        end
        counter += 1
    end
    return counter
end

# DependentParameterRefs
function _check_element_support(_error::Function, prefs::Vector{DependentParameterRef},
                                param_values::Vector{Float64},
                                counter::Int)::Int
    len = length(prefs)
    supp = reshape(param_values[counter:counter+len-1], len, 1)
    if !supports_in_set(supp, infinite_set(prefs))
        _error("Parameter values violate parameter bounds.")
    end
    return counter += len
end

# Used to ensure values don't violate parameter bounds
function _check_tuple_values(_error::Function, ivref::InfiniteVariableRef,
                             param_values::Vector{Float64})::Nothing
    raw_prefs = raw_parameter_refs(ivref)
    counter = 1
    for i in 1:size(raw_prefs, 1)
        prefs = map(e -> dispatch_variable_ref(e), raw_prefs[i, :])
        counter = _check_element_support(_error, prefs, param_values, counter)
    end
    return
end

# Update point variable info to consider the infinite variable
function _update_point_info(info::JuMP.VariableInfo,
                            ivref::InfiniteVariableRef)::JuMP.VariableInfo
    if JuMP.has_lower_bound(ivref) && !info.has_fix && !info.has_lb
        info = JuMP.VariableInfo(true, JuMP.lower_bound(ivref),
                                 info.has_ub, info.upper_bound,
                                 info.has_fix, info.fixed_value,
                                 info.has_start, info.start,
                                 info.binary, info.integer)
    end
    if JuMP.has_upper_bound(ivref) && !info.has_fix && !info.has_ub
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 true, JuMP.upper_bound(ivref),
                                 info.has_fix, info.fixed_value,
                                 info.has_start, info.start,
                                 info.binary, info.integer)
    end
    if JuMP.is_fixed(ivref) && !info.has_fix  && !info.has_lb  && !info.has_ub
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 info.has_ub, info.upper_bound,
                                 true, JuMP.fix_value(ivref),
                                 info.has_start, info.start,
                                 info.binary, info.integer)
    end
    if !(JuMP.start_value(ivref) === NaN) && !info.has_start
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 info.has_ub, info.upper_bound,
                                 info.has_fix, info.fixed_value,
                                 true, JuMP.start_value(ivref),
                                 info.binary, info.integer)
    end
    if JuMP.is_binary(ivref) && !info.integer
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 info.has_ub, info.upper_bound,
                                 info.has_fix, info.fixed_value,
                                 info.has_start, info.start,
                                 true, info.integer)
    end
    if JuMP.is_integer(ivref) && !info.binary
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 info.has_ub, info.upper_bound,
                                 info.has_fix, info.fixed_value,
                                 info.has_start, info.start,
                                 info.binary, true)
    end
    return info
end

# Define _make_variable (used by JuMP.build_variable)
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Point};
                        infinite_variable_ref::Union{GeneralVariableRef,
                                                     Nothing} = nothing,
                        parameter_values::Union{Real, AbstractArray{<:Real},
                                                Tuple, Nothing} = nothing,
                        extra_kw_args...)::PointVariable{GeneralVariableRef}
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with point variables.")
    end
    # ensure the needed arguments are given
    if parameter_values === nothing || infinite_variable_ref === nothing
        _error("Must specify the infinite variable and the values of its " *
               "infinite parameters")
    end
    # format tuple as VectorTuple
    pvalues = VectorTuple{Float64}(parameter_values)
    # check information and prepare format
    dispatch_ivref = dispatch_variable_ref(infinite_variable_ref)
    if !(dispatch_ivref isa InfiniteVariableRef)
     _error("Expected an infinite variable reference dependency, but got a " *
            "variable reference of type $(typeof(dispatch_ivref)).")
    end
    _check_tuple_shape(_error, dispatch_ivref, pvalues)
    _check_tuple_values(_error, dispatch_ivref, pvalues.values)
    info = _update_point_info(info, dispatch_ivref)
    # make variable and return
    return PointVariable(_make_float_info(info), infinite_variable_ref,
                         pvalues.values)
end

## Dispatch methods for updating the supports of parameters in point variable
# IndependentParameterRefs
function _add_point_support(prefs::Vector{IndependentParameterRef},
                            param_values::Vector{Float64},
                            counter::Int)::Int
    for pref in prefs
        add_supports(pref, param_values[counter], check = false,
                     label = UserDefined)
        counter += 1
    end
    return counter
end

# DependentParameterRefs
function _add_point_support(prefs::Vector{DependentParameterRef},
                            param_values::Vector{Float64},
                            counter::Int)::Int
    len = length(prefs)
    supp = reshape(param_values[counter:counter+len-1], len, 1)
    add_supports(prefs, supp, check = false, label = UserDefined)
    return counter += len
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(ivref::InfiniteVariableRef,
                                param_values::Vector{Float64})::Nothing
    raw_prefs = raw_parameter_refs(ivref)
    model = JuMP.owner_model(ivref)
    counter = 1
    for i in 1:size(raw_prefs, 1)
        prefs = dispatch_variable_ref.(raw_prefs[i, :])
        counter = _add_point_support(prefs, param_values, counter)
    end
    return
end

# Used to update mapping infinite_to_points
function _update_infinite_point_mapping(pvref::PointVariableRef,
                                        ivref::InfiniteVariableRef)::Nothing
    push!(_point_variable_dependencies(ivref), JuMP.index(pvref))
    return
end

# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(model::InfiniteModel,
                                      v::PointVariable)::PointVariableRef
    ivref = dispatch_variable_ref(v.infinite_variable_ref)
    JuMP.check_belongs_to_model(ivref, model)
    data_object = VariableData(v)
    vindex = _add_data_object(model, data_object)
    vref = PointVariableRef(model, vindex)
    _update_param_supports(ivref, v.parameter_values)
    _update_infinite_point_mapping(vref, ivref)
    return vref
end

################################################################################
#                         PARAMETER VALUE METHODS
################################################################################

"""
    infinite_variable_ref(vref::PointVariableRef)::GeneralVariableRef

Return the `InfiniteVariableRef` associated with the point variable `vref`.

**Example**
```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> infinite_variable_ref(vref)
T(t)
```
"""
function infinite_variable_ref(vref::PointVariableRef)::GeneralVariableRef
    return _core_variable_object(vref).infinite_variable_ref
end


"""
    raw_parameter_values(vref::PointVariableRef)::Vector{Float64}

Return the raw support point values associated with the point variable `vref`.
"""
function raw_parameter_values(vref::PointVariableRef)::Vector{Float64}
    return _core_variable_object(vref).parameter_values
end

"""
    parameter_values(vref::PointVariableRef)::Tuple

Return the support point associated with the point variable `vref`.

**Example**
```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> parameter_values(vref)
(0,)
```
"""
function parameter_values(vref::PointVariableRef)::Tuple
    prefs = raw_parameter_refs(infinite_variable_ref(vref))
    return Tuple(VectorTuple(raw_parameter_values(vref), prefs.ranges,
                             prefs.indices))
end

# Internal function used to change the parameter value tuple of a point variable
function _update_variable_param_values(vref::PointVariableRef,
                                       pref_vals::Vector{<:Real})::Nothing
    info = _variable_info(vref)
    ivref = infinite_variable_ref(vref)
    new_var = PointVariable(info, ivref, pref_vals)
    _set_core_variable_object(vref, new_var)
    return
end

################################################################################
#                         VARIABLE NAMING METHODS
################################################################################
# Get root name of infinite variable
function _root_name(vref::InfiniteVariableRef)
    name = JuMP.name(vref)
    return name[1:findfirst(isequal('('), name)-1]
end

## Return the parameter value as an appropriate string
# Number
function _make_str_value(value)::String
    return JuMP._string_round(value)
end

# Array{<:Number}
function _make_str_value(values::Array)::String
    if length(values) == 1
        return _make_str_value(first(values))
    end
    if length(values) <= 4
        str_value = "["
        for i in eachindex(values)
            if i != length(values)
                str_value *= JuMP._string_round(values[i]) * ", "
            else
                str_value *= JuMP._string_round(values[i]) * "]"
            end
        end
        return str_value
    else
        return string("[", JuMP._string_round(first(values)), ", ..., ",
                      JuMP._string_round(last(values)), "]")
    end
end

"""
    JuMP.set_name(vref::PointVariableRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
the names of point variables.

**Example**
```julia-repl
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> set_name(vref, "new_name")

julia> name(vref)
"new_name"
```
"""
function JuMP.set_name(vref::PointVariableRef, name::String)::Nothing
    if length(name) == 0
        ivref = dispatch_variable_ref(infinite_variable_ref(vref))
        prefs = raw_parameter_refs(ivref)
        name = _root_name(ivref)
        values = raw_parameter_values(vref)
        name = string(name, "(")
        for i in 1:size(prefs, 1)
            if i != size(prefs, 1)
                name *= string(_make_str_value(values[prefs.ranges[i]]), ", ")
            else
                name *= string(_make_str_value(values[prefs.ranges[i]]), ")")
            end
        end
    end
    _data_object(vref).name = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for point variables
function _update_variable_info(vref::PointVariableRef,
                               info::JuMP.VariableInfo)::Nothing
    ivref = infinite_variable_ref(vref)
    param_values = raw_parameter_values(vref)
    new_var = PointVariable(info, ivref, param_values)
    _set_core_variable_object(vref, new_var)
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::PointVariableRef)::Nothing
    # remove variable info constraints associated with vref
    _delete_info_constraints(vref)
    # remove the infinite variable dependency
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref), _point_variable_dependencies(ivref))
    return
end
