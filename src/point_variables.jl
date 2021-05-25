################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::PointVariableIndex
    )::PointVariableRef
    return PointVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:PointVariable}
    )::PointVariableIndex
    return MOIUC.add_item(model.point_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel, 
    ::Type{PointVariable}
    )::MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}
    return model.point_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(
    vref::PointVariableRef
    )::MOIUC.CleverDict{PointVariableIndex, VariableData{PointVariable{GeneralVariableRef}}}
    return JuMP.owner_model(vref).point_vars
end

# Extend _data_object
function _data_object(
    vref::PointVariableRef
    )::VariableData{PointVariable{GeneralVariableRef}}
    object = get(_data_dictionary(vref), JuMP.index(vref), nothing)
    object === nothing && error("Invalid point variable reference, cannot find " *
                        "corresponding variable in the model. This is likely " *
                        "caused by using the reference of a deleted variable.")
    return object
end

# Extend _core_variable_object
function _core_variable_object(vref::PointVariableRef
    )::PointVariable{GeneralVariableRef}
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
"""
    Point{V, T} <: InfOptVariableType 

A `DataType` to assist in making point variables. This can be passed as an 
extra argument to `@variable` to make such a variable: 
```julia 
@variable(model, var_expr, Point(inf_var, parameter_values...), args..., 
          kwargs...)
```
Here `parameter_values` must match the format of the infinite parameter 
references associated with the infinite variable `inf_var` and can be comprised 
of both real valued supports.

**Fields**
- `infinite_variable_ref::V`
- `parameter_values::VectorTuple{T}`: The infinite parameter support values the 
   variable will depend on.
"""
struct Point{V, T} <: InfOptVariableType 
    infinite_variable_ref::V 
    parameter_values::VectorTuple{T}
    function Point(ivref::V, vals...) where {V}
        vt = VectorTuple(vals)
        T = param_type(vt)
        return new{V, T}(ivref, vt)
    end
end

# Ensure parameter values match shape of parameter reference tuple stored in the
# infinite variable reference
function _check_tuple_shape(
    _error::Function,
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    values::VectorTuple
    )::Nothing
    prefs = raw_parameter_refs(ivref)
    if !same_structure(prefs, values)
        _error("The dimensions and array formatting of the infinite parameter ",
               "values must match those of the parameter references for the ",
               "infinite variable/derivative.")
    end
    return
end

## Dispatch methods for checking the supports of parameters in point variable
# IndependentParameterRefs
function _check_element_support(
    _error::Function, 
    prefs::Vector{IndependentParameterRef},
    param_values::Vector{Float64},
    counter::Int
    )::Int
    for pref in prefs
        if !supports_in_domain(param_values[counter], infinite_domain(pref))
            _error("Parameter values violate parameter bounds.")
        end
        counter += 1
    end
    return counter
end

# DependentParameterRefs
function _check_element_support(
    _error::Function, 
    prefs::Vector{DependentParameterRef},
    param_values::Vector{Float64},
    counter::Int
    )::Int
    len = length(prefs)
    supp = reshape(param_values[counter:counter+len-1], len, 1)
    if !supports_in_domain(supp, infinite_domain(prefs))
        _error("Parameter values violate parameter bounds.")
    end
    return counter += len
end

# Used to ensure values don't violate parameter bounds
function _check_tuple_values(
    _error::Function, 
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    param_values::Vector{Float64}
    )::Nothing
    raw_prefs = raw_parameter_refs(ivref)
    counter = 1
    for i in 1:size(raw_prefs, 1)
        prefs = map(e -> dispatch_variable_ref(e), raw_prefs[i, :])
        counter = _check_element_support(_error, prefs, param_values, counter)
    end
    return
end

# Update point variable info to consider the infinite variable
function _update_point_info(
    info::JuMP.VariableInfo,
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    point::Vector{Float64}
    )::JuMP.VariableInfo
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
    if start_value_function(ivref) !== nothing && !info.has_start
        if _is_vector_start(ivref)
            start = start_value_function(ivref)(point)
        else
            prefs = raw_parameter_refs(ivref)
            start = start_value_function(ivref)(Tuple(point, prefs)...)
        end
        info = JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                 info.has_ub, info.upper_bound,
                                 info.has_fix, info.fixed_value,
                                 true, Float64(start),
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

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, 
                        var_type::Point)::InfiniteVariable{GeneralVariableRef}

Build and return a point variable based on `info` and `var_type`. Errors 
if the information stored in `var_type` is invalid. See [`Point`](@ref) 
for more information.

**Example**
```julia-repl
julia> y
y(t)

julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> pt_var = build_variable(error, info, SemiInfinite(y, 0));
```
"""
function JuMP.build_variable(
    _error::Function, 
    info::JuMP.VariableInfo, 
    var_type::Point;
    extra_kwargs...
    )::PointVariable{GeneralVariableRef}
    # check for unneeded keywords
    for (kwarg, _) in extra_kwargs
        _error("Keyword argument $kwarg is not for use with point variables.")
    end
    # check the infinite variable reference
    ivref = var_type.infinite_variable_ref
    if !(ivref isa GeneralVariableRef)
        _error("Expected an infinite variable/derivative reference dependency ",
               "of type `GeneralVariableRef`, but got an argument of type ",
               "$(typeof(ivref)).")
    end
    dispatch_ivref = dispatch_variable_ref(ivref)
    if !(dispatch_ivref isa Union{InfiniteVariableRef, DerivativeRef})
        _error("Expected an infinite variable/derivative reference dependency,", 
               "but got a variable reference of type $(typeof(dispatch_ivref)).")
    end
    # check and format the values 
    raw_vals = var_type.parameter_values
    if !(param_type(raw_vals) <: Real)
        _error("Expected parameter values consisting of real numbers.")
    end
    pvalues = Vector{Float64}(raw_vals.values)
    _check_tuple_shape(_error, dispatch_ivref, raw_vals)
    _check_tuple_values(_error, dispatch_ivref, pvalues)
    info = _update_point_info(info, dispatch_ivref, pvalues)
    # enforce parameter significant digits on the values
    prefs = parameter_list(dispatch_ivref)
    for i in eachindex(pvalues)
        pvalues[i] = round(pvalues[i], sigdigits = significant_digits(prefs[i]))
    end
    # make variable and return
    return PointVariable(_make_float_info(info), ivref, pvalues)
end

## Dispatch methods for updating the supports of parameters in point variable
# IndependentParameterRefs
function _add_point_support(
    prefs::Vector{IndependentParameterRef},
    param_values::Vector{Float64},
    counter::Int
    )::Int
    for pref in prefs
        add_supports(pref, param_values[counter], check = false,
                     label = UserDefined)
        counter += 1
    end
    return counter
end

# DependentParameterRefs
function _add_point_support(
    prefs::Vector{DependentParameterRef},
    param_values::Vector{Float64},
    counter::Int
    )::Int
    len = length(prefs)
    supp = reshape(param_values[counter:counter+len-1], len, 1)
    add_supports(prefs, supp, check = false, label = UserDefined)
    return counter += len
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    param_values::Vector{Float64}
    )::Nothing
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
function _update_infinite_point_mapping(
    pvref::PointVariableRef,
    ivref::Union{InfiniteVariableRef, DerivativeRef}
    )::Nothing
    push!(_point_variable_dependencies(ivref), JuMP.index(pvref))
    return
end

# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(
    model::InfiniteModel,
    v::PointVariable,
    name::String;
    add_support = true
    )::PointVariableRef
    ivref = dispatch_variable_ref(v.infinite_variable_ref)
    JuMP.check_belongs_to_model(ivref, model)
    data_object = VariableData(v, name)
    vindex = _add_data_object(model, data_object)
    vref = PointVariableRef(model, vindex)
    if add_support
        _update_param_supports(ivref, v.parameter_values)
    end
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
julia> @variable(model, T, Infinite(t))
T(t)

julia> @variable(model, T0, Point(T, 0))
T0

julia> infinite_variable_ref(T0)
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
julia> @variable(model, T, Infinite(t))
T(t)

julia> @variable(model, T0, Point(T, 0))
T0

julia> parameter_values(T0)
(0,)
```
"""
function parameter_values(vref::PointVariableRef)::Tuple
    prefs = raw_parameter_refs(infinite_variable_ref(vref))
    return Tuple(VectorTuple(raw_parameter_values(vref), prefs.ranges,
                             prefs.indices))
end

# Internal function used to change the parameter value tuple of a point variable
function _update_variable_param_values(
    vref::PointVariableRef,
    pref_vals::Vector{<:Real}
    )::Nothing
    info = _variable_info(vref)
    ivref = infinite_variable_ref(vref)
    new_var = PointVariable(info, ivref, pref_vals)
    _set_core_variable_object(vref, new_var)
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for point variables
function _update_variable_info(
    vref::PointVariableRef,
    info::JuMP.VariableInfo
    )::Nothing
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
