################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::PointModel,
                               index::PointVariableIndex
                               )::PointVariableRef
    return PointVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::PointModel,
                          object::VariableData{<:PointVariable}
                          )::PointVariableIndex
    return MOIUC.add_item(model.infinite_vars, object)
end

# Extend _data_dictionary
function _data_dictionary(vref::PointVariableRef)::MOIUC.CleverDict
    return model.infinite_vars
end

# Extend _data_object
function _data_object(vref::PointVariableRef)::VariableData{<:PointVariable}
    return _data_dictionary(vref)[JuMP.index(vref)]
end

# Extend _core_variable_object
function _core_variable_object(vref::PointVariableRef)::PointVariable
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

# Used to ensure values don't violate parameter bounds
function _check_tuple_values(_error::Function, ivref::InfiniteVariableRef,
                             param_values::Vector{Float64})::Nothing
    prefs = raw_parameter_refs(ivref)
    for i in eachindex(prefs)
        if !supports_in_set(param_values[i], infinite_set(prefs[i]))
            _error("Parameter values violate parameter bounds.")
        end
    end
    return
end

# Update point variable info to consider the infinite variable
function _update_point_info(info::JuMP.VariableInfo,
                            ivref::InfiniteVariableRef)::Nothing
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
                        parameter_values::Union{Real,
                                                AbstractArray{<:Real},
                                                Tuple, Nothing} = nothing,
                        extra_kw_args...)::PointVariable
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
    return PointVariable(info, infinite_variable_ref, pvalues.values)
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(ivref::InfiniteVariableRef,
                                param_values::Vector)::Nothing
    prefs = raw_parameter_refs(ivref)
    for i in eachindex(prefs)
        add_supports(prefs[i], param_values[i])
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
    JuMP.is_valid(model, ivref) || error("Invalid infinite variable reference.")
    data_object = VariableData(v)
    vindex = _add_data_object(model, data_object)
    vref = PointVariableRef(model, vindex)
    _update_param_supports(ivref, v.parameter_values)
    _update_infinite_point_mapping(vref, ivref)
    return vref
end

################################################################################
#                            VARIABLE NAMING
################################################################################

"""
    infinite_variable_ref(vref::PointVariableRef)::InfiniteVariableRef

Return the `InfiniteVariableRef` associated with the point variable `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> infinite_variable_ref(vref)
T(t)
```
"""
function infinite_variable_ref(vref::PointVariableRef)::InfiniteVariableRef
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
end

"""
    raw_parameter_values(vref::PointVariableRef)::VectorTuple{<:Number}

Return the raw [`VectorTuple`](@ref) support point associated with the
point variable `vref`.
```
"""
function raw_parameter_values(vref::PointVariableRef)::VectorTuple{<:Number}
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
end

"""
    parameter_values(vref::PointVariableRef)::Tuple

Return the support point associated with the point variable `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> parameter_values(vref)
(0,)
```
"""
function parameter_values(vref::PointVariableRef)::Tuple
    return Tuple(raw_parameter_values(vref))
end

# Internal function used to change the parameter value tuple of a point variable
function _update_variable_param_values(vref::PointVariableRef,
                                       pref_vals::VectorTuple)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    ivref = JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = PointVariable(info, ivref,
                                                                  pref_vals)
    return
end

# Get root name of infinite variable
function _root_name(vref::InfiniteVariableRef)
    name = JuMP.name(vref)
    return name[1:findfirst(isequal('('), name)-1]
end

## Return the parameter value as an appropriate string
# Number
function _make_str_value(value)::String
    return string(JuMP._string_round(value))
end

# Array{<:Number}
function _make_str_value(values::Array)::String
    if length(values) == 1
        return _make_str_value(first(values))
    end
    if length(values) <= 4
        str_value = "["
        counter = 1
        for value in values
            if counter != length(values)
                str_value *= JuMP._string_round(value) * ", "
            else
                str_value *= JuMP._string_round(value) * "]"
            end
            counter += 1
        end
        return string(str_value)
    else
        return string("[", JuMP._string_round(first(values)), ", ..., ",
                      JuMP._string_round(last(values)), "]")
    end
end

"""
    JuMP.set_name(vref::PointVariableRef, name::String)

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
the names of point variables.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> @infinite_variable(model, T(t))
T(t)

julia> vref = @point_variable(model, T(0))
T(0)

julia> set_name(vref, "new_name")

julia> name(vref)
"new_name"
```
"""
function JuMP.set_name(vref::PointVariableRef, name::String)
    if length(name) == 0
        ivref = infinite_variable_ref(vref::PointVariableRef)
        name = _root_name(ivref)
        values = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
        name = string(name, "(")
        for i in 1:size(values, 1)
            if i != size(values, 1)
                name *= _make_str_value(values[i, :]) * ", "
            else
                name *= _make_str_value(values[i, :]) * ")"
            end
        end
    end
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end
