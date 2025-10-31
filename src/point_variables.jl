################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel, index::PointVariableIndex)
    return PointVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:PointVariable}
    )
    return MOIUC.add_item(model.point_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{PointVariable})
    return model.point_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(vref::PointVariableRef)
    return JuMP.owner_model(vref).point_vars
end

# Extend _data_object
function _data_object(vref::PointVariableRef)
    object = get(_data_dictionary(vref), JuMP.index(vref), nothing)
    if isnothing(object) 
        error("Invalid point variable reference, cannot find ",
        "corresponding variable in the model. This is likely ",
        "caused by using the reference of a deleted variable.")
    end
    return object
end

"""
    core_object(vref::PointVariableRef)::PointVariable

Retrieve the underlying core [`PointVariable`] object for `vref`. 
This is intended as an advanced method for developers.
"""
function core_object(vref::PointVariableRef)
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
# build point variable without checks
function _build_point_variable(
    ivref::GeneralVariableRef, 
    support::Vector{<:Real},
    pref_list::Vector{GeneralVariableRef},
    info::RestrictedDomainInfo
    )
    # enforce parameter significant digits on the values
    pvalues = Vector{Float64}(support)
    for i in eachindex(pvalues)
        pvalues[i] = round(pvalues[i], sigdigits = significant_digits(pref_list[i]))
    end
    # make variable and return
    return PointVariable(info, ivref, pvalues)
end
function _build_point_variable(
    ivref::GeneralVariableRef, 
    support::Vector{<:Real},
    pref_list::Vector{GeneralVariableRef}
    )
    return _build_point_variable(ivref, support, pref_list, RestrictedDomainInfo())
end

# Ensure parameter values match shape of parameter reference tuple stored in the
# infinite variable reference
function _check_tuple_shape(
    _error::Function,
    prefs::Collections.VectorTuple{GeneralVariableRef},
    values::Collections.VectorTuple
    )
    if !Collections.same_structure(prefs, values)
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
    param_values::Vector{<:Real}
    )
    if !supports_in_domain(only(param_values), infinite_domain(only(prefs)))
        _error("Parameter values violate parameter bounds.")
    end
    return
end

# DependentParameterRefs
function _check_element_support(
    _error::Function, 
    prefs::Vector{DependentParameterRef},
    param_values::Vector{<:Real}
    )
    supp = reshape(param_values, length(prefs), 1)
    if !supports_in_domain(supp, infinite_domain(prefs))
        _error("Parameter values violate parameter bounds.")
    end
    return
end

# Used to ensure values don't violate parameter bounds
function _check_tuple_values(
    _error::Function, 
    prefs::Collections.VectorTuple{GeneralVariableRef},
    param_values::Vector{<:Real}
    )
    for i in 1:size(prefs, 1)
        dprefs = map(e -> dispatch_variable_ref(e), prefs[i, :])
        supp = param_values[prefs.ranges[i]]
        _check_element_support(_error, dprefs, supp)
    end
    return
end

"""
    JuMP.build_variable(
        [_error::Function],
        ivref::GeneralVariableRef, 
        support::Collections.VectorTuple{<:Real},
        [info::RestrictedDomainInfo]
    )::PointVariable

Construct a `PointVariable` that corresponds to the infinite variable `ivref`
at the point defined by `support`. Users are encourage to construct point 
variable via the functional restriction syntax, e.g., `ivref(support...)`
or using `@variable` if you want to override the variable domain information.
"""
function JuMP.build_variable(
    _error::Function,
    ivref::GeneralVariableRef, 
    support::Collections.VectorTuple{<:Real},
    info::RestrictedDomainInfo
    )
    # check the infinite variable reference
    if !(_index_type(ivref) in (InfiniteVariableIndex, DerivativeIndex, ParameterFunctionIndex))
        _error("Expected an infinite variable/derivative reference dependency,", 
               "but got `$(ivref)`.")
    end
    prefs = raw_parameter_refs(ivref)
    _check_tuple_shape(_error, prefs, support)
    _check_tuple_values(_error, prefs, support.values)
    # make the point variable
    return _build_point_variable(ivref, support.values, prefs.values, info)
end

"""
    Point{T <: Real} <: InfOptVariableType 

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
- `infinite_variable_ref::GeneralVariableRef`
- `parameter_values::VectorTuple{T}`: The infinite parameter support values the 
   variable will depend on.
"""
struct Point{T} <: InfOptVariableType 
    infinite_variable_ref::GeneralVariableRef
    parameter_values::Collections.VectorTuple{T}
    function Point(
        vref::GeneralVariableRef,
        vt::Collections.VectorTuple{T}
        ) where {T}
        return new{T}(vref, vt)
    end
end
function Point(ivref::GeneralVariableRef, vals...)
    if isempty(vals)
        error("No point values given for point variable.")
    end
    return Point(ivref, Collections.VectorTuple(vals))
end

# Convert JuMP.VariableInfo to RestrictedDomainInfo
function _process_restricted_info(
    _error::Function,
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Real}
    )
    if info.binary || info.integer
        _error("Cannot set the integrality of a point or semi-infinite " *
               "variable. They only inherit the integrality of the infinite " *
               "variable they originate from.")
    end
    return RestrictedDomainInfo(
        true, info.has_lb ? info.lower_bound : NaN,
        true, info.has_ub ? info.upper_bound : NaN,
        true, info.has_fix ? info.fixed_value : NaN,
        true, info.has_start ? info.start : NaN
        )
end
function _process_restricted_info(_error::Function, info::JuMP.VariableInfo)
    _error("Bounds and start value for point and semi-infinite variables must " *
           "be of type `Real`.")
end

"""
    JuMP.build_variable(
        _error::Function,
        info::JuMP.VariableInfo, 
        var_type::Point
    )::PointVariable

Build and return a point variable based on `info` and `var_type`. Errors 
if the information stored in `var_type` is invalid. See [`Point`](@ref) 
for more information. This is intended to enable the use of `@variable`.
"""
function JuMP.build_variable(
    _error::Function, 
    info::JuMP.VariableInfo, 
    var_type::Point;
    extra_kwargs...
    )
    # check for unneeded keywords
    if !isempty(extra_kwargs)
        _error("Keyword argument `$(first(keys(extra_kwargs)))`` is not " *
               "for use with point variables.")
    end
    # check and format the values 
    raw_vals = var_type.parameter_values
    if !(raw_vals.values isa Vector{<:Real})
        _error("Expected parameter values consisting of real numbers.")
    end
    restricted_info = _process_restricted_info(_error, info)
    # build the parameter
    return JuMP.build_variable(
        _error,
        var_type.infinite_variable_ref,
        raw_vals,
        restricted_info
    )
end

## Dispatch methods for updating the supports of parameters in point variable
# IndependentParameterRefs
function _add_point_support(
    prefs::Vector{IndependentParameterRef},
    param_values::Vector{Float64}
    )
    return add_supports(
        only(prefs),
        only(param_values), 
        check = false, 
        label = UserDefined
    )
end

# DependentParameterRefs
function _add_point_support(
    prefs::Vector{DependentParameterRef},
    param_values::Vector{Float64}
    )
    supp = reshape(param_values, length(prefs), 1)
    add_supports(prefs, supp, check = false, label = UserDefined)
    return
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(
    ivref::Union{InfiniteVariableRef, DerivativeRef, ParameterFunctionRef},
    param_values::Vector{Float64}
    )
    raw_prefs = raw_parameter_refs(ivref)
    for i in 1:size(raw_prefs, 1)
        prefs = dispatch_variable_ref.(raw_prefs[i, :])
        supp = param_values[raw_prefs.ranges[i]]
        _add_point_support(prefs, supp)
    end
    return
end

# Used to update mapping infinite_to_points
function _update_infinite_point_mapping(
    pvref::PointVariableRef,
    ivref::Union{InfiniteVariableRef, DerivativeRef, ParameterFunctionRef}
    )
    push!(_point_variable_dependencies(ivref), JuMP.index(pvref))
    return
end

"""
    JuMP.add_variable(
        model::InfiniteModel,
        var::PointVariable,
        [name::String = ""]
    )::GeneralVariableRef

Extend the `JuMP.add_variable` function to accomodate `PointVariable` variable 
types. Adds a variable to an infinite model `model` and returns a 
[`GeneralVariableRef`](@ref). Primarily intended to be an internal function of 
the constructor macro `@variable`. However, it can be used in combination with
`JuMP.build_variable` to add variables to an infinite model object.
Errors if an invalid infinite variable reference is included in `var`.

**Example**
```julia-repl
julia> @infinite_parameter(m, t in [0, 10]);

julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite(t));

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, ivref, Collections.VectorTuple(0.5));

julia> pvref = add_variable(m, pt_var, "var_alias")
var_alias
```
"""
function JuMP.add_variable(
    model::InfiniteModel,
    v::PointVariable,
    name::String = "";
    add_support = true
    )
    ivref = v.infinite_variable_ref
    info = v.info
    JuMP.check_belongs_to_model(ivref, model)
    existing_index = get(model.point_lookup, (ivref, v.parameter_values), nothing)
    is_active_info = info.active_lower_bound_info || info.active_upper_bound_info || 
                    info.active_fix_info || info.active_start_info
    if isnothing(existing_index)
        divref = dispatch_variable_ref(ivref)
        data_object = VariableData(v, name)
        vindex = _add_data_object(model, data_object)
        vref = PointVariableRef(model, vindex)
        if add_support
            _update_param_supports(divref, v.parameter_values)
        end
        _update_infinite_point_mapping(vref, divref)
        model.point_lookup[(ivref, v.parameter_values)] = vindex
    else
        vref = PointVariableRef(model, existing_index)
        if !isempty(name)
            JuMP.set_name(vref, name)
        end
        if is_active_info
            _delete_info_constraints(vref)
            _set_core_object(vref, v)
        end
    end
    gvref = GeneralVariableRef(vref)
    if is_active_info
        _set_info_constraints(info, gvref, vref)
    end
    model.name_to_var = nothing
    return gvref
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
function infinite_variable_ref(vref::PointVariableRef)
    return core_object(vref).infinite_variable_ref
end


"""
    raw_parameter_values(vref::PointVariableRef)::Vector{Float64}

Return the raw support point values associated with the point variable `vref`.
"""
function raw_parameter_values(vref::PointVariableRef)
    return core_object(vref).parameter_values
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
function parameter_values(vref::PointVariableRef)
    prefs = raw_parameter_refs(infinite_variable_ref(vref))
    return Tuple(raw_parameter_values(vref), prefs)
end

# Internal function used to change the parameter value tuple of a point variable
function _update_variable_param_values(
    vref::PointVariableRef,
    pref_vals::Vector{<:Real}
    )
    info = _variable_info(vref)
    ivref = infinite_variable_ref(vref)
    new_var = PointVariable(info, ivref, pref_vals)
    _set_core_object(vref, new_var)
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for point variables
function _update_variable_info(
    vref::PointVariableRef,
    info::RestrictedDomainInfo
    )
    ivref = infinite_variable_ref(vref)
    param_values = raw_parameter_values(vref)
    new_var = PointVariable(info, ivref, param_values)
    _set_core_object(vref, new_var)
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::PointVariableRef)
    # remove the infinite variable dependency
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref), _point_variable_dependencies(ivref))
    # remove the lookup entry
    delete!(JuMP.owner_model(vref).point_lookup, (ivref, raw_parameter_values(vref)))
    return
end
