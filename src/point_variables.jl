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
"""
    Point{V, VT <: VectorTuple} <: InfOptVariableType 

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
- `parameter_values::VT`: The infinite parameter support values the 
   variable will depend on.
"""
struct Point{V, VT <: Collections.VectorTuple} <: InfOptVariableType 
    infinite_variable_ref::V 
    parameter_values::VT
    function Point(vref::V, vt::VT) where {V, VT <: Collections.VectorTuple}
        return new{V, VT}(vref, vt)
    end
end
function Point(ivref, vals...)
    return Point(ivref, Collections.VectorTuple(vals))
end

# Ensure parameter values match shape of parameter reference tuple stored in the
# infinite variable reference
function _check_tuple_shape(
    _error::Function,
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    values::Collections.VectorTuple
    )
    prefs = raw_parameter_refs(ivref)
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
    param_values::Vector{Float64},
    counter::Int
    )
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
    )
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
    )
    raw_prefs = raw_parameter_refs(ivref)
    counter = 1
    for i in 1:size(raw_prefs, 1)
        prefs = map(e -> dispatch_variable_ref(e), raw_prefs[i, :])
        counter = _check_element_support(_error, prefs, param_values, counter)
    end
    return
end

# Convert JuMP.VariableInfo to RestrictedVariableDomainInfo
function _process_restricted_info(_error::Function, info::JuMP.VariableInfo)
    if info.is_binary || info.is_integer
        _error("Cannot set the integrality of a point or semi-infinite " *
               "variable. They only inherit the integrality of the infinite " *
               "variable they originate from.")
    end
    return RestrictedVariableDomainInfo(
        info.has_lb, info.lower_bound,
        info.has_ub, info.upper_bound,
        info.has_fix, info.fixed_value,
        info.has_start, info.start
        )
end

"""
    JuMP.build_variable(
        _error::Function,
        info::JuMP.VariableInfo, 
        var_type::Point
    )::InfiniteVariable{GeneralVariableRef}

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
    )
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
    if !(raw_vals.values isa Vector{<:Real})
        _error("Expected parameter values consisting of real numbers.")
    end
    pvalues = Vector{Float64}(raw_vals.values)
    _check_tuple_shape(_error, dispatch_ivref, raw_vals)
    _check_tuple_values(_error, dispatch_ivref, pvalues)
    restricted_info = _process_restricted_info(_error, info)
    # enforce parameter significant digits on the values
    prefs = parameter_list(dispatch_ivref)
    for i in eachindex(pvalues)
        pvalues[i] = round(pvalues[i], sigdigits = significant_digits(prefs[i]))
    end
    # make variable and return
    return PointVariable(restricted_info, ivref, pvalues)
end

## Dispatch methods for updating the supports of parameters in point variable
# IndependentParameterRefs
function _add_point_support(
    prefs::Vector{IndependentParameterRef},
    param_values::Vector{Float64},
    counter::Int
    )
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
    )
    len = length(prefs)
    supp = reshape(param_values[counter:counter+len-1], len, 1)
    add_supports(prefs, supp, check = false, label = UserDefined)
    return counter += len
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(
    ivref::Union{InfiniteVariableRef, DerivativeRef},
    param_values::Vector{Float64}
    )
    raw_prefs = raw_parameter_refs(ivref)
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
    )
    push!(_point_variable_dependencies(ivref), JuMP.index(pvref))
    return
end

# TODO CONTINUE FROM HERE

# Update the information constraints for a variable that was already created
function _update_info_constraints(info::JuMP.VariableInfo, gvref, vref)::Nothing
    # extract the preliminaries
    old_info = _variable_info(vref)
    old_info == info && return
    model = JuMP.owner_model(vref)

    # update the lower bound constraint as needed
    if info.has_lb && (!old_info.has_lb || info.lower_bound != old_info.lower_bound)
        constr = JuMP.ScalarConstraint(gvref, MOI.GreaterThan(info.lower_bound))
        if old_info.has_lb
            cindex = _lower_bound_index(vref)
            cref = InfOptConstraintRef(model, cindex)
            _set_core_object(cref, constr)
        else
            cref = JuMP.add_constraint(model, constr, is_info_constr = true)
            _set_lower_bound_index(vref, JuMP.index(cref))
        end
    elseif old_info.has_lb && !info.has_lb
        JuMP.delete(model, JuMP.LowerBoundRef(vref))
        _set_lower_bound_index(vref, nothing)
    end

    # update the upper bound constraint as needed
    if info.has_ub && (!old_info.has_ub || info.upper_bound != old_info.upper_bound)
        constr = JuMP.ScalarConstraint(gvref, MOI.LessThan(info.upper_bound))
        if old_info.has_ub
            cindex = _upper_bound_index(vref)
            cref = InfOptConstraintRef(model, cindex)
            _set_core_object(cref, constr)
        else
            cref = JuMP.add_constraint(model, constr, is_info_constr = true)
            _set_upper_bound_index(vref, JuMP.index(cref))
        end
    elseif old_info.has_ub && !info.has_ub
        JuMP.delete(model, JuMP.UpperBoundRef(vref))
        _set_upper_bound_index(vref, nothing)
    end
    
    # update the fix constraint as needed
    if info.has_fix && (!old_info.has_fix || info.fixed_value != old_info.fixed_value)
        constr = JuMP.ScalarConstraint(gvref, MOI.EqualTo(info.fixed_value))
        if old_info.has_fix
            cindex = _fix_index(vref)
            cref = InfOptConstraintRef(model, cindex)
            _set_core_object(cref, constr)
        else
            cref = JuMP.add_constraint(model, constr, is_info_constr = true)
            _set_fix_index(vref, JuMP.index(cref))
        end
    elseif old_info.has_fix && !info.has_fix
        JuMP.delete(model, JuMP.FixRef(vref))
        _set_fix_index(vref, nothing)
    end
    
    # update the binary constraint as needed
    if info.binary && !old_info.binary
        constr = JuMP.ScalarConstraint(gvref, MOI.ZeroOne())
        cref = JuMP.add_constraint(model, constr, is_info_constr = true)
        _set_binary_index(vref, JuMP.index(cref))
    elseif old_info.binary && !info.binary
        JuMP.delete(model, JuMP.BinaryRef(vref))
        _set_binary_index(vref, nothing)
    end

    # update the integer constraint as needed
    if info.integer && !old_info.integer
        constr = JuMP.ScalarConstraint(gvref, MOI.Integer())
        cref = JuMP.add_constraint(model, constr, is_info_constr = true)
        _set_integer_index(vref, JuMP.index(cref))
    elseif old_info.integer && !info.integer
        JuMP.delete(model, JuMP.IntegerRef(vref))
        _set_integer_index(vref, nothing)
    end

    # finalize the update
    _update_variable_info(vref, info)
    set_transformation_backend_ready(model, false)
    return
end

"""
    JuMP.add_variable(model::InfiniteModel, var::PointVariable,
                      [name::String = ""])::GeneralVariableRef

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

julia> pt_var = build_variable(error, info, Point(ivref, 0.5));

julia> pvref = add_variable(m, pt_var, "var_alias")
var_alias
```
"""
function JuMP.add_variable(
    model::InfiniteModel,
    v::PointVariable,
    name::String = "";
    add_support = true,
    update_info = true
    )
    ivref = v.infinite_variable_ref
    divref = dispatch_variable_ref(ivref)
    JuMP.check_belongs_to_model(divref, model)
    existing_index = get(model.point_lookup, (ivref, v.parameter_values), nothing)
    if isnothing(existing_index)
        data_object = VariableData(v, name)
        vindex = _add_data_object(model, data_object)
        vref = PointVariableRef(model, vindex)
        if add_support
            _update_param_supports(divref, v.parameter_values)
        end
        _update_infinite_point_mapping(vref, divref)
        model.point_lookup[(ivref, v.parameter_values)] = vindex
        gvref = GeneralVariableRef(model, vindex)
        _set_info_constraints(v.info, gvref, vref)
    else
        gvref = GeneralVariableRef(model, existing_index)
        if update_info
            vref = PointVariableRef(model, existing_index)
            _update_info_constraints(v.info, gvref, vref)
        end
        if !isempty(name)
            JuMP.set_name(vref, name)
        end
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
    info::JuMP.VariableInfo
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
function _delete_variable_dependencies(vref::PointVariableRef)::Nothing
    # remove the infinite variable dependency
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref), _point_variable_dependencies(ivref))
    # remove the lookup entry
    delete!(JuMP.owner_model(vref).point_lookup, (ivref, raw_parameter_values(vref)))
    return
end
