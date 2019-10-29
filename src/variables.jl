# Define symbol inputs for different variable types
const Infinite = :Infinite
const Point = :Point
const Hold = :Hold

## Extend Base.copy for new variable types
# GeneralVariableRef
Base.copy(v::GeneralVariableRef) = v

# InfiniteVariableRef
function Base.copy(v::InfiniteVariableRef,
                   new_model::InfiniteModel)::InfiniteVariableRef
    return InfiniteVariableRef(new_model, v.index)
end

# HoldVariableRef
function Base.copy(v::HoldVariableRef,
                   new_model::InfiniteModel)::HoldVariableRef
    return HoldVariableRef(new_model, v.index)
end

# PointVariableRef
function Base.copy(v::PointVariableRef,
                   new_model::InfiniteModel)::PointVariableRef
     return PointVariableRef(new_model, v.index)
 end

# Extend other Base functions
function Base.:(==)(v::T, w::U)::Bool where {T <: GeneralVariableRef,
                                             U <: GeneralVariableRef}
    return v.model === w.model && v.index == w.index && T == U
end
Base.broadcastable(v::GeneralVariableRef) = Ref(v)

# Extend JuMP functions
JuMP.isequal_canonical(v::GeneralVariableRef, w::GeneralVariableRef) = v == w
JuMP.variable_type(model::InfiniteModel) = GeneralVariableRef
function JuMP.variable_type(model::InfiniteModel, type::Symbol)
    if type == Infinite
        return InfiniteVariableRef
    elseif type == Point
        return PointVariableRef
    elseif type == Hold
        return HoldVariableRef
    elseif type == Parameter
        return ParameterRef
    else
        error("Invalid variable type.")
    end
end

# Check parameter tuple, ensure all elements contain parameter references
function _check_parameter_tuple(_error::Function, prefs::Tuple)
    types = [typeof(pref) for pref in prefs]
    num_params = length(types)
    valid_types = zeros(Bool, num_params)
    for i in eachindex(types)
        if types[i] == ParameterRef || types[i] <: AbstractArray{<:ParameterRef}
            valid_types[i] = true
        end
    end
    if !all(valid_types)
        _error("Invalid parameter type(s) given.")
    end
    return
end

# Convert parameter tuple s.t. array elements are SparseAxisArrays
function _make_formatted_tuple(prefs::Tuple)::Tuple
    converted_prefs = ()
    for pref in prefs
        if isa(pref, ParameterRef) || isa(pref, Number)
            converted_prefs = (converted_prefs..., pref)
        else
            converted_prefs = (converted_prefs...,
                               convert(JuMPC.SparseAxisArray, pref))
        end
    end
    return converted_prefs
end

# Ensure each element onyl contains parameters with same group ID
function _check_tuple_groups(_error::Function, prefs::Tuple)
    valid_elements = _only_one_group.(prefs)
    if sum(valid_elements) != length(prefs)
        _error("Each parameter tuple element must have contain only infinite " *
               "parameters with the same group ID.")
    end
    groups = _group.(prefs)
    if length(unique(groups)) != length(groups)
        _error("Cannot double specify infinite parameter references.")
    end
    return
end

# Ensure parameter values match shape of parameter reference tuple stored in the
# infinite variable reference
function _check_tuple_shape(_error::Function,
                            infinite_variable_ref::InfiniteVariableRef,
                            values::Tuple)
    prefs = parameter_refs(infinite_variable_ref)
    container = JuMPC.SparseAxisArray
    if length(prefs) != length(values)
        _error("The dimensions of the infinite parameter values must match " *
               "those defined for the infinite variable.")
    end
    for i in eachindex(values)
        if isa(prefs[i], ParameterRef) && !(isa(values[i], Number))
            _error("The dimensions and array type of the infinite parameter " *
                   "values must match those defined for the infinite variable.")
        elseif isa(prefs[i], container) && !isa(values[i], container)
            _error("The dimensions and array type of the infinite parameter " *
                   "values must match those defined for the infinite variable.")
        elseif isa(prefs[i], container)
            if keys(prefs[i].data) != keys(values[i].data)
                _error("Index keys of infinite parameter values don't match " *
                       "those defined for the infinite variable.")
            end
        end
    end
    return
end

# Used to ensure values don't violate parameter bounds
function _check_tuple_values(_error::Function, inf_vref::InfiniteVariableRef,
                             param_values::Tuple)
    prefs = parameter_refs(inf_vref)
    for i in eachindex(prefs)
        if isa(prefs[i], ParameterRef) && JuMP.has_lower_bound(prefs[i])
            check1 = param_values[i] < JuMP.lower_bound(prefs[i])
            check2 = param_values[i] > JuMP.upper_bound(prefs[i])
            if check1 || check2
                _error("Parameter values violate parameter bounds.")
            end
        else
            for (k, v) in prefs[i].data
                if JuMP.has_lower_bound(v)
                    check1 = param_values[i].data[k] < JuMP.lower_bound(v)
                    check2 = param_values[i].data[k] > JuMP.upper_bound(v)
                    if check1 || check2
                        _error("Parameter values violate parameter bounds.")
                    end
                end
            end
        end
    end
    return
end

# Update point variable info to consider the infinite variable
function _update_point_info(info::JuMP.VariableInfo, ivref::InfiniteVariableRef)
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

# Check that parameter_bounds argument is valid
function _check_bounds(bounds::ParameterBounds; _error = error)
    for (pref, set) in bounds.intervals
        # check that respects lower bound
        if JuMP.has_lower_bound(pref) && (set.lower_bound < JuMP.lower_bound(pref))
                _error("Specified parameter lower bound exceeds that defined " *
                       "for $pref.")
        end
        # check that respects upper bound
        if JuMP.has_upper_bound(pref) && (set.upper_bound > JuMP.upper_bound(pref))
                _error("Specified parameter upper bound exceeds that defined " *
                       "for $pref.")
        end
    end
    return
end

## Check to ensure correct inputs and build variables and return
# InfiniteVariable
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Infinite};
                        parameter_refs::Union{ParameterRef,
                                              AbstractArray{<:ParameterRef},
                                              Tuple, Nothing} = nothing,
                        extra_kw_args...)::InfiniteVariable
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with infinite variables.")
    end
    # check that we have been given parameter references
    if parameter_refs == nothing
        _error("Parameter references not specified, use the var(params...) " *
               "syntax or the parameter_refs keyword argument.")
    end
    # make sure the parameters are in a tuple
    if !isa(parameter_refs, Tuple)
        parameter_refs = (parameter_refs, )
    end
    # check tuple for validity and format
    _check_parameter_tuple(_error, parameter_refs)
    parameter_refs = _make_formatted_tuple(parameter_refs)
    _check_tuple_groups(_error, parameter_refs)
    # make the variable and return
    return InfiniteVariable(info, parameter_refs)
end

# PointVariable
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Point};
                        infinite_variable_ref::Union{InfiniteVariableRef,
                                                     Nothing} = nothing,
                        parameter_values::Union{Number,
                                                AbstractArray{<:Number},
                                                Tuple, Nothing} = nothing,
                        extra_kw_args...)::PointVariable
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with point variables.")
    end
    # ensure the needed arguments are given
    if parameter_values == nothing || infinite_variable_ref == nothing
        _error("Must specify the infinite variable and the values of its " *
               "infinite parameters")
    end
    # format as tuple if needed
    if !isa(parameter_values, Tuple)
        parameter_values = (parameter_values, )
    end
    # check information and prepare format
    parameter_values = _make_formatted_tuple(parameter_values)
    _check_tuple_shape(_error, infinite_variable_ref, parameter_values)
    _check_tuple_values(_error, infinite_variable_ref, parameter_values)
    info = _update_point_info(info, infinite_variable_ref)
    # make variable and return
    return PointVariable(info, infinite_variable_ref, parameter_values)
end

# HoldVariable
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Hold};
                        parameter_bounds::ParameterBounds = ParameterBounds(),
                        extra_kw_args...)::HoldVariable
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with hold variables.")
    end
    # check that the bounds don't violate parameter domains
    _check_bounds(parameter_bounds)
    # make variable and return
    return HoldVariable(info, parameter_bounds)
end

# Fallback method
function _make_variable(_error::Function, info::JuMP.VariableInfo, type;
                        extra_kw_args...)
    _error("Unrecognized variable type $type, should be Infinite, " *
           "Point, or Hold.")
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                        var_type::Symbol;
                        parameter_refs::Union{ParameterRef,
                                              AbstractArray{<:ParameterRef},
                                              Tuple, Nothing} = nothing,
                        infinite_variable_ref::Union{InfiniteVariableRef,
                                                     Nothing} = nothing,
                        parameter_values::Union{Number, AbstractArray{<:Number},
                                                Tuple, Nothing} = nothing,
                        parameter_bounds::Union{Dict{ParameterRef, IntervalSet},
                                                Nothing} = nothing,
                        extra_kw_args...)

Extend the [`JuMP.build_variable`](@ref) function to accomodate `InfiniteOpt`
variable types. Returns the appropriate variable Datatype (i.e.,
[`InfiniteVariable`](@ref), [`PointVariable`](@ref), and
[`HoldVariable`](@ref)). Primarily this method is to be used internally by the
appropriate constructor macros [`@infinite_variable`](@ref),
[`@point_variable`](@ref), and [`@hold_variable`](@ref). However, it can be
called manually to build `InfiniteOpt` variables. Errors if an unneeded keyword
argument is given or if the keywoard arguments are formatted incorrectly (e.g.,
`parameter_refs` contains repeated parameter references when an infinite variable
is defined). Also errors if needed kewword arguments are negated.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel())
julia> @infinite_parameter(m, 0 <= t <= 1)
t

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = t)
InfiniteVariable{Int64,Int64,Int64,Int64}(VariableInfo{Int64,Int64,Int64,Int64}(false, 0, false, 0, false, 0, false, 0, false, false), (t,))

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, info, Point, infinite_variable_ref = ivref,
                               parameter_values = 0.5)
PointVariable{Int64,Int64,Int64,Float64}(VariableInfo{Int64,Int64,Int64,Float64}(false, 0, false, 0, false, 0, true, 0.0, false, false), var_name(t), (0.5,))

julia> hd_var = build_variable(error, info, Hold)
HoldVariable{Int64,Int64,Int64,Int64}(VariableInfo{Int64,Int64,Int64,Int64}(false, 0, false, 0, false, 0, false, 0, false, false), Subdomain bounds (0): )
```
"""
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                             var_type::Symbol;
                             macro_error::Union{Function, Nothing} = nothing,
                             kw_args...)
    if macro_error != nothing
        _error = macro_error # replace with macro error function
    end
    # make the variable and conduct necessary checks
    return _make_variable(_error, info, Val(var_type); kw_args...)
end

# Used to update the model.param_to_vars field
function _update_param_var_mapping(vref::InfiniteVariableRef, prefs::Tuple)
    model = JuMP.owner_model(vref)
    pref_list = _list_parameter_refs(prefs)
    for pref in pref_list
        if haskey(model.param_to_vars, JuMP.index(pref))
            push!(model.param_to_vars[JuMP.index(pref)], JuMP.index(vref))
        else
            model.param_to_vars[JuMP.index(pref)] = [JuMP.index(vref)]
        end
    end
    return
end

# check the pref tuple contains only valid parameters
function _check_parameters_valid(model::InfiniteModel, prefs::Tuple)
    pref_list = _list_parameter_refs(prefs)
    for pref in pref_list
        !JuMP.is_valid(model, pref) && error("Invalid Parameter reference " *
                                             "provided.")
    end
    return
end

# Used to add point variable support to parameter supports if necessary
function _update_param_supports(inf_vref::InfiniteVariableRef,
                                param_values::Tuple)
    prefs = parameter_refs(inf_vref)
    for i in eachindex(prefs)
        if isa(prefs[i], ParameterRef)
            add_supports(prefs[i], param_values[i])
        else
            for (k, v) in prefs[i].data
                add_supports(v, param_values[i].data[k])
            end
        end
    end
    return
end

# Used to update mapping infinite_to_points
function _update_infinite_point_mapping(pvref::PointVariableRef,
                                        ivref::InfiniteVariableRef)
    model = JuMP.owner_model(pvref)
    if haskey(model.infinite_to_points, JuMP.index(ivref))
        push!(model.infinite_to_points[JuMP.index(ivref)], JuMP.index(pvref))
    else
        model.infinite_to_points[JuMP.index(ivref)] = [JuMP.index(pvref)]
    end
    return
end

# Validate parameter bounds and add support(s) if needed
function _validate_bounds(model::InfiniteModel, bounds::ParameterBounds; _error = error)
    for (pref, set) in bounds.intervals
        # check validity
        !JuMP.is_valid(model, pref) && _error("Parameter bound reference " *
                                              "is invalid.")
        # ensure has a support if a point constraint was given
        if set.lower_bound == set.upper_bound
            add_supports(pref, set.lower_bound)
        end
    end
    return
end

## Make the variable reference and do checks/mapping updates
# InfiniteVariable
function _check_make_variable_ref(model::InfiniteModel,
                                  v::InfiniteVariable)::InfiniteVariableRef
    _check_parameters_valid(model, v.parameter_refs)
    vref = InfiniteVariableRef(model, model.next_var_index)
    _update_param_var_mapping(vref, v.parameter_refs)
    return vref
end

# PointVariable
function _check_make_variable_ref(model::InfiniteModel,
                                  v::PointVariable)::PointVariableRef
    ivref = v.infinite_variable_ref
    !JuMP.is_valid(model, ivref) && error("Invalid infinite variable " *
                                          "reference.")
    vref = PointVariableRef(model, model.next_var_index)
    _update_param_supports(ivref, v.parameter_values)
    _update_infinite_point_mapping(vref, ivref)
    return vref
end

# HoldVariable
function _check_make_variable_ref(model::InfiniteModel,
                                  v::HoldVariable)::HoldVariableRef
    _validate_bounds(model, v.parameter_bounds)
    vref = HoldVariableRef(model, model.next_var_index)
    if length(v.parameter_bounds.intervals) != 0
        model.has_hold_bounds = true
    end
    return vref
end

# Fallback
function _check_make_variable_ref(model::InfiniteModel, v)
    error("Invalid variable object type.")
end

"""
    JuMP.add_variable(model::InfiniteModel, var::InfOptVariable, name::String = "")

Extend the [`JuMP.add_variable`](@ref) function to accomodate `InfiniteOpt`
variable types. Adds a variable to an infinite model `model` and returns an
appropriate variable reference (i.e., [`InfiniteVariableRef`](@ref),
[`PointVariableRef`](@ref), or [`HoldVariableRef`](@ref)). Primarily intended
to be an internal function of the constructor macros [`@infinite_variable`](@ref),
[`@point_variable`](@ref), and [`@hold_variable`](@ref). However, it can be used
in combination with [`JuMP.build_variable`](@ref) to add variables to an infinite
model object. Errors if invalid parameters reference(s) or an invalid infinite
variable reference is included in `var`.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel())
julia> @infinite_parameter(m, t in [0, 10]);

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = t);

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, info, Point, infinite_variable_ref = ivref,
                               parameter_values = 0.5);

julia> pvref = add_variable(m, pt_var, "var_alias")
var_alias

julia> hd_var = build_variable(error, info, Hold);

julia> hvref = add_variable(m, hd_var, "var_name")
var_name
```
"""
function JuMP.add_variable(model::InfiniteModel, var::InfOptVariable,
                           name::String = "")
    model.next_var_index += 1
    vref = _check_make_variable_ref(model, var)
    model.vars[JuMP.index(vref)] = var
    JuMP.set_name(vref, name)
    if var.info.has_lb
        newset = MOI.GreaterThan(convert(Float64, var.info.lower_bound))
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, newset))
        _set_lower_bound_index(vref, JuMP.index(cref))
        model.constr_in_var_info[JuMP.index(cref)] = true
    end
    if var.info.has_ub
        newset = MOI.LessThan(convert(Float64, var.info.upper_bound))
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, newset))
        _set_upper_bound_index(vref, JuMP.index(cref))
        model.constr_in_var_info[JuMP.index(cref)] = true
    end
    if var.info.has_fix
        newset = MOI.EqualTo(convert(Float64, var.info.fixed_value))
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(vref, newset))
        _set_fix_index(vref, JuMP.index(cref))
        model.constr_in_var_info[JuMP.index(cref)] = true
    end
    if var.info.binary
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, MOI.ZeroOne()))
        _set_binary_index(vref, JuMP.index(cref))
        model.constr_in_var_info[JuMP.index(cref)] = true
    elseif var.info.integer
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, MOI.Integer()))
        _set_integer_index(vref, JuMP.index(cref))
        model.constr_in_var_info[JuMP.index(cref)] = true
    end
    model.var_in_objective[JuMP.index(vref)] = false
    return vref
end

"""
    JuMP.owner_model(vref::GeneralVariableRef)::InfiniteModel

Extend [`JuMP.owner_model`](@ref) function for `InfiniteOpt` variables. Returns
the infinite model associated with `vref`.

**Example**
```julia
julia> owner_model(vref)
An InfiniteOpt Model
Feasibility problem with:
Variable: 1
`HoldVariableRef`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
`HoldVariableRef`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`HoldVariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
Names registered in the model: vref
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
JuMP.owner_model(vref::GeneralVariableRef)::InfiniteModel = vref.model

"""
    JuMP.index(v::GeneralVariableRef)::Int

Extent [`JuMP.index`](@ref) to return the index of a `InfiniteOpt` variable.

**Example**
```julia
julia> index(vref)
1
```
"""
JuMP.index(v::GeneralVariableRef)::Int = v.index

"""
    used_by_constraint(vref::InfOptVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a constraint.

**Example**
```julia
julia> used_by_constraint(vref)
false
```
"""
function used_by_constraint(vref::InfOptVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).var_to_constrs, JuMP.index(vref))
end

"""
    used_by_measure(vref::InfOptVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a measure.

**Example**
```julia
julia> used_by_measure(vref)
true
```
"""
function used_by_measure(vref::InfOptVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).var_to_meas, JuMP.index(vref))
end

"""
    used_by_objective(vref::InfOptVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by the objective.

**Example**
```julia
julia> used_by_objective(vref)
true
```
"""
function used_by_objective(vref::InfOptVariableRef)::Bool
    return JuMP.owner_model(vref).var_in_objective[JuMP.index(vref)]
end

"""
    is_used(vref::InfOptVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia
julia> is_used(vref)
true
```
"""
function is_used(vref::InfOptVariableRef)::Bool
    return used_by_measure(vref) || used_by_constraint(vref) || used_by_objective(vref)
end

"""
    used_by_point_variable(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a point variable.

**Example**
```julia
julia> used_by_point_variable(vref)
false
```
"""
function used_by_point_variable(vref::InfiniteVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).infinite_to_points, JuMP.index(vref))
end

"""
    used_by_reduced_variable(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a reduced infinite variable.

**Example**
```julia
julia> used_by_reduced_variable(vref)
true
```
"""
function used_by_reduced_variable(vref::InfiniteVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).infinite_to_reduced, JuMP.index(vref))
end

"""
    is_used(vref::InfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia
julia> is_used(vref)
false
```
"""
function is_used(vref::InfiniteVariableRef)::Bool
    if used_by_measure(vref) || used_by_constraint(vref)
        return true
    end
    if used_by_point_variable(vref)
        for vindex in JuMP.owner_model(vref).infinite_to_points[JuMP.index(vref)]
            if is_used(PointVariableRef(JuMP.owner_model(vref), vindex))
                return true
            end
        end
    end
    if used_by_reduced_variable(vref)
        for rindex in JuMP.owner_model(vref).infinite_to_reduced[JuMP.index(vref)]
            rvref = ReducedInfiniteVariableRef(JuMP.owner_model(vref), rindex)
            if used_by_constraint(rvref) || used_by_measure(rvref)
                return true
            end
        end
    end
    return false
end

"""
    JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)

Extend [`JuMP.delete`](@ref) to delete `InfiniteOpt` variables and their
dependencies. Errors if variable is invalid, meaning it has already been
deleted or it belongs to another model.

**Example**
```julia
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z >= 0.0
 g(t) + z >= 42.0
 g(0.5) == 0
 t in [0, 6]

julia> delete(model, g)

julia> print(model)
Min measure(t) + z
Subject to
 z >= 0.0
 z >= 42.0
 t in [0, 6]
```
"""
function JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
    @assert JuMP.is_valid(model, vref) "Variable is invalid."
    # update the optimizer model status
    if is_used(vref)
        set_optimizer_model_ready(model, false)
    end
    # remove variable info constraints associated with vref
    if JuMP.has_lower_bound(vref)
        JuMP.delete_lower_bound(vref)
    end
    if JuMP.has_upper_bound(vref)
        JuMP.delete_upper_bound(vref)
    end
    if JuMP.is_fixed(vref)
        JuMP.unfix(vref)
    end
    if JuMP.is_binary(vref)
        JuMP.unset_binary(vref)
    elseif JuMP.is_integer(vref)
        JuMP.unset_integer(vref)
    end
    # remove dependencies from measures and update them
    if used_by_measure(vref)
        for mindex in model.var_to_meas[JuMP.index(vref)]
            if isa(model.measures[mindex].func, InfOptVariableRef)
                model.measures[mindex] = Measure(zero(JuMP.AffExpr),
                                                 model.measures[mindex].data)
            else
                _remove_variable(model.measures[mindex].func, vref)
            end
            JuMP.set_name(MeasureRef(model, mindex),
                           _make_meas_name(model.measures[mindex]))
        end
        # delete mapping
        delete!(model.var_to_meas, JuMP.index(vref))
    end
    # remove dependencies from measures and update them
    if used_by_constraint(vref)
        for cindex in model.var_to_constrs[JuMP.index(vref)]
            if isa(model.constrs[cindex].func, InfOptVariableRef)
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr),
                                                      model.constrs[cindex].set)
            else
                _remove_variable(model.constrs[cindex].func, vref)
            end
        end
        # delete mapping
        delete!(model.var_to_constrs, JuMP.index(vref))
    end
    # remove from objective if vref is in it
    if used_by_objective(vref)
        if isa(model.objective_function, InfOptVariableRef)
            model.objective_function = zero(JuMP.AffExpr)
        else
            _remove_variable(model.objective_function, vref)
        end
    end
    # do specific updates if vref is infinite
    if isa(vref, InfiniteVariableRef)
        # update parameter mapping
        all_prefs = _list_parameter_refs(parameter_refs(vref))
        for pref in all_prefs
            filter!(e -> e != JuMP.index(vref),
                    model.param_to_vars[JuMP.index(pref)])
            if length(model.param_to_vars[JuMP.index(pref)]) == 0
                delete!(model.param_to_vars, JuMP.index(pref))
            end
        end
        # delete associated point variables and mapping
        if used_by_point_variable(vref)
            for index in model.infinite_to_points[JuMP.index(vref)]
                JuMP.delete(model, PointVariableRef(model, index))
            end
            delete!(model.infinite_to_points, JuMP.index(vref))
        end
        # delete associated reduced variables and mapping
        if used_by_reduced_variable(vref)
            for index in model.infinite_to_reduced[JuMP.index(vref)]
                JuMP.delete(model, ReducedInfiniteVariableRef(model, index))
            end
            delete!(model.infinite_to_reduced, JuMP.index(vref))
        end
    end
    # update mappings if is point variable
    if isa(vref, PointVariableRef)
        ivref = infinite_variable_ref(vref)
        filter!(e -> e != JuMP.index(vref),
                model.infinite_to_points[JuMP.index(ivref)])
        if length(model.infinite_to_points[JuMP.index(ivref)]) == 0
            delete!(model.infinite_to_points, JuMP.index(ivref))
        end
    end
    # delete the variable information
    delete!(model.var_in_objective, JuMP.index(vref))
    delete!(model.vars, JuMP.index(vref))
    delete!(model.var_to_name, JuMP.index(vref))
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, vref::InfOptVariableRef)::Bool

Extend [`JuMP.is_valid`](@ref) to accomodate `InfiniteOpt` variables.

**Example**
```julia
julia> is_valid(model, ivref)
true
```
"""
function JuMP.is_valid(model::InfiniteModel, vref::InfOptVariableRef)::Bool
    return (model === JuMP.owner_model(vref) && JuMP.index(vref) in keys(model.vars))
end

"""
    JuMP.num_variables(model::InfiniteModel)::Int

Extend [`JuMP.num_variables`](@ref) to return the number of `InfiniteOpt`
variables assigned to `model`.

**Example**
```julia
julia> num_variables(model)
3
```
"""
JuMP.num_variables(model::InfiniteModel)::Int = length(model.vars)

# Include all the extension functions for manipulating the properties associated
# with VariableInfo
include("variable_info.jl")

"""
    JuMP.name(vref::InfOptVariableRef)::String

Extend [`JuMP.name`](@ref) to return the names of `InfiniteOpt` variables.

**Example**
```julia
julia> name(vref)
"var_name"
```
"""
function JuMP.name(vref::InfOptVariableRef)::String
    return JuMP.owner_model(vref).var_to_name[JuMP.index(vref)]
end

"""
    JuMP.set_name(vref::HoldVariableRef, name::String)

Extend [`JuMP.set_name`](@ref) to set names of hold variables.

**Example**
```julia
julia> set_name(hvref, "var_name")

julia> name(t)
"var_name"
```
"""
function JuMP.set_name(vref::HoldVariableRef, name::String)
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    infinite_variable_ref(vref::PointVariableRef)::InfiniteVariableRef

Return the `InfiniteVariableRef` associated with the point variable `vref`.

**Example**
```julia
julia> infinite_variable_ref(vref)
T(t, x)
```
"""
function infinite_variable_ref(vref::PointVariableRef)::InfiniteVariableRef
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
end

"""
    parameter_values(vref::PointVariableRef)::Tuple

Return the support point associated with the point variable `vref`.

**Example**
```julia
julia> parameter_values(vref)
(0, )
```
"""
function parameter_values(vref::PointVariableRef)::Tuple
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
end

# Internal function used to change the parameter value tuple of a point variable
function _update_variable_param_values(vref::PointVariableRef, pref_vals::Tuple)
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
function _make_str_value(value::Number)::String
    return string(JuMP._string_round(value))
end

# AbstractArray
function _make_str_value(value::AbstractArray)::String
    if length(keys(value)) <= 4
        str_value = "["
        counter = 1
        for key in sort(collect(keys(value)))
            if counter != length(keys(value))
                str_value *= JuMP._string_round(value[key]) * ", "
            else
                str_value *= JuMP._string_round(value[key]) * "]"
            end
            counter += 1
        end
        return string(str_value)
    else
        value = [value[key] for key in sort(collect(keys(value)))]
        return string("[", JuMP._string_round(first(value)), ", ..., ",
                      JuMP._string_round(last(value)), "]")
    end
end

"""
    JuMP.set_name(vref::PointVariableRef, name::String)

Extend [`JuMP.set_name`](@ref) to set the names of point variables.

**Example**
```julia
julia> name(vref)
old_name

julia> set_name(vref, "new_name")

julia> name(vref)
new_name
```
"""
function JuMP.set_name(vref::PointVariableRef, name::String)
    if length(name) == 0
        inf_var_ref = infinite_variable_ref(vref::PointVariableRef)
        name = _root_name(inf_var_ref)
        values = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
        name = string(name, "(")
        for i in eachindex(values)
            if i != length(values)
                name *= _make_str_value(values[i]) * ", "
            else
                name *= _make_str_value(values[i]) * ")"
            end
        end
    end
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    parameter_refs(vref::InfiniteVariableRef)::Tuple

Return the `ParameterRef`(s) associated with the infinite variable `vref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `vref`.

**Example**
```julia
julia> parameter_refs(vref)
(t,   [2]  =  x[2]
  [1]  =  x[1])
```
"""
function parameter_refs(vref::InfiniteVariableRef)
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_refs
end

# Internal function used to change the parameter reference tuple of an infinite
# variable
function _update_variable_param_refs(vref::InfiniteVariableRef, prefs::Tuple)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, prefs)
    return
end

"""
    set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)

Specify a new parameter reference tuple `prefs` for the infinite variable `vref`.
Note each element must contain a single parameter reference or an array of
parameter references. Errors if a parameter is double specified or if an element
contains parameters with different group IDs.

**Example**
```julia
julia> set_parameter_refs(vref, (t, x))

julia> parameter_refs(vref)
(t,   [2]  =  x[2]
  [1]  =  x[1])
```
"""
function set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)
    _check_parameter_tuple(error, prefs)
    prefs = _make_formatted_tuple(prefs)
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
are already added to the variable or if the added parameters have different
group IDs.

```julia
julia> name(vref)
T(t)

julia> add_parameter_ref(vref, x)

julia> name(vref)
T(t, x)
```
"""
function add_parameter_ref(vref::InfiniteVariableRef,
                       pref::Union{ParameterRef, AbstractArray{<:ParameterRef}})
    set_parameter_refs(vref, (parameter_refs(vref)..., pref))
    return
end

# TODO implement test for all of these new methods and include in docs
"""
    parameter_bounds(vref::HoldVariableRef)::ParameterBounds

Return the [`ParameterBounds`](@ref) object associated with the hold variable
`vref`. It contains a dictionary where each key is a `ParameterRef` which points
to an `IntervalSet` that that defines a sub-domain for `vref` relative to that
parameter reference.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @hold_variable(model, vref, parameter_bounds = (t in [0, 2]))
vref

julia> parameter_bounds(vref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function parameter_bounds(vref::HoldVariableRef)::ParameterBounds
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_bounds
end

"""
    has_parameter_bounds(vref::HoldVariableRef)::Bool

Return a `Bool` indicating if `vref` is limited to a sub-domain as defined
by parameter bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @hold_variable(model, vref, parameter_bounds = (t in [0, 2]))
vref

julia> has_parameter_bounds(vref)
true
```
"""
function has_parameter_bounds(vref::HoldVariableRef)::Bool
    return length(parameter_bounds(vref)) != 0
end

# Other variable types
function has_parameter_bounds(vref::GeneralVariableRef)::Bool
    return false
end

# Internal function used to change the parameter bounds of a hold variable
function _update_variable_param_bounds(vref::HoldVariableRef,
                                       bounds::ParameterBounds)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = HoldVariable(info, bounds)
    return
end

## Check that the bounds dictionary is compadable with existing dependent measures
# DiscreteMeasureData
function _check_meas_bounds(bounds::ParameterBounds, data::DiscreteMeasureData;
                            _error = error)
    pref = data.parameter_ref
    supports = data.supports
    if haskey(bounds.intervals, pref)
        if bounds.intervals[pref].lower_bound > minimum(supports) ||
            bounds.intervals[pref].upper_bound < maximum(supports)
            _error("New bounds don't span existing dependent measure bounds.")
        end
    end
    return
end

# MultiDiscreteMeasureData
function _check_meas_bounds(bounds::ParameterBounds,
                            data::MultiDiscreteMeasureData; _error = error)
    prefs = data.parameter_ref
    supports = data.supports
    mins = minimum(supports)
    maxs = maximum(supports)
    for key in keys(prefs)
        if haskey(bounds.intervals, prefs[key])
            if bounds.intervals[prefs[key]].lower_bound > mins[key] ||
                bounds.intervals[prefs[key]].upper_bound < maxs[key]
                _error("New bounds don't span existing dependent measure bounds.")
            end
        end
    end
    return
end

# Fallback
function _check_meas_bounds(bounds::ParameterBounds, data::AbstractMeasureData;
                            _error = error)
    @warn "Unable to check if hold variables bounds are valid in measure with" *
          " custom measure data type."
    return
end

# Update the current bounds to overlap with the new bounds if possible
function _update_bounds(bounds1::Dict, bounds2::Dict; _error = error)
    # check each new bound
    for (pref, set) in bounds2
        # we have a new bound
        if !haskey(bounds1, pref)
            bounds1[pref] = set
        # the previous set and the new one do not overlap
        elseif set.lower_bound > bounds1[pref].upper_bound || set.upper_bound < bounds1[pref].lower_bound
            _error("Sub-domains of constraint and/or hold variable(s) do not" *
                   " overlap. Consider changing the parameter bounds of the" *
                   " constraint and/or hold variable(s).")
        # we have an existing bound
        else
            # we have a new stricter lower bound to update with
            if set.lower_bound > bounds1[pref].lower_bound
                bounds1[pref] = IntervalSet(set.lower_bound, bounds1[pref].upper_bound)
            end
            # we have a new stricter upper bound to update with
            if set.upper_bound < bounds1[pref].upper_bound
                bounds1[pref] = IntervalSet(bounds1[pref].lower_bound, set.upper_bound)
            end
        end
    end
    return
end

## Check and update the constraint bounds (don't change in case of error)
# BoundedScalarConstraint
function _update_constr_bounds(bounds::ParameterBounds, c::BoundedScalarConstraint;
                               _error = error)
    new_bounds_dict = copy(c.bounds.intervals)
    _update_bounds(new_bounds_dict, bounds.intervals, _error = _error)
    return BoundedScalarConstraint(c.func, c.set, ParameterBounds(new_bounds_dict),
                                   c.orig_bounds)
end

# ScalarConstraint
function _update_constr_bounds(bounds::ParameterBounds, c::JuMP.ScalarConstraint;
                               _error = error)
    return BoundedScalarConstraint(c.func, c.set, bounds, ParameterBounds())
end

"""
    set_parameter_bounds(vref::HoldVariableRef, bounds::ParameterBounds;
                         [force = false])

Specify a new dictionary of parameter bounds `bounds` for the hold variable `vref`.
These are stored in a [`ParameterBounds`](@ref) object which contains a dictionary.
Note the dictionary keys must be `ParameterRef`s and the values must be
`IntervalSet`s that indicate a particular sub-domain for which `vref` is defined.
This is meant to be primarily used by [`@set_parameter_bounds`](@ref) which
provides a more intuitive syntax.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @hold_variable(model, vref)
vref

julia> set_parameter_bounds(vref, ParameterBounds(Dict(t => IntervalSet(0, 2))))

julia> parameter_bounds(vref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function set_parameter_bounds(vref::HoldVariableRef, bounds::ParameterBounds;
                              force = false, _error = error)
    if has_parameter_bounds(vref) && !force
        _error("$vref already has parameter bounds. Consider adding more using " *
               "`add_parameter_bounds` or overwriting them by setting " *
               "the keyword argument `force = true`")
    else
        # check that bounds are valid and add support(s) if necessary
        _check_bounds(bounds, _error = _error)
        # check dependent measures
        meas_cindices = []
        if used_by_measure(vref)
            for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
                meas = JuMP.owner_model(vref).measures[mindex]
                _check_meas_bounds(bounds, meas.data, _error = _error)
                if used_by_constraint(MeasureRef(JuMP.owner_model(vref), mindex))
                    indices = JuMP.owner_model(vref).meas_to_constrs[mindex]
                    meas_cindices = [meas_cindices; indices]
                end
            end
        end
        # check and update dependent constraints
        if used_by_constraint(vref) || length(meas_cindices) != 0
            for cindex in unique([meas_cindices; JuMP.owner_model(vref).var_to_constrs[JuMP.index(vref)]])
                constr = JuMP.owner_model(vref).constrs[cindex]
                new_constr = _update_constr_bounds(bounds, constr, _error = _error)
                JuMP.owner_model(vref).constrs[cindex] = new_constr
            end
        end
        _validate_bounds(JuMP.owner_model(vref), bounds, _error = _error)
        # set the new bounds
        _update_variable_param_bounds(vref, bounds)
        # update status
        JuMP.owner_model(vref).has_hold_bounds = true
        if is_used(vref)
            set_optimizer_model_ready(JuMP.owner_model(vref), false)
        end
    end
    return
end

"""
    add_parameter_bound(vref::HoldVariableRef, pref::ParameterRef,
                        lower::Number, upper::Number)

Add an additional parameter bound to `vref` such that it is defined over the
sub-domain based on `pref` from `lower` to `upper`. This is primarily meant to be
used by [`@add_parameter_bounds`](@ref).

```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @hold_variable(model, vref)
vref

julia> add_parameter_bound(vref, t, 0, 2)

julia> parameter_bounds(vref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function add_parameter_bound(vref::HoldVariableRef, pref::ParameterRef,
                             lower::Number, upper::Number; _error = error)
    # check the new bounds
    new_bounds = ParameterBounds(Dict(pref => IntervalSet(lower, upper)))
    _check_bounds(new_bounds, _error = _error)
    # check dependent measures
    meas_cindices = []
    if used_by_measure(vref)
        for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
            meas = JuMP.owner_model(vref).measures[mindex]
            _check_meas_bounds(new_bounds, meas.data, _error = _error)
            if used_by_constraint(MeasureRef(JuMP.owner_model(vref), mindex))
                indices = JuMP.owner_model(vref).meas_to_constrs[mindex]
                meas_cindices = [meas_cindices; indices]
            end
        end
    end
    # check and update dependent constraints
    if used_by_constraint(vref) || length(meas_cindices) != 0
        for cindex in unique([meas_cindices; JuMP.owner_model(vref).var_to_constrs[JuMP.index(vref)]])
            constr = JuMP.owner_model(vref).constrs[cindex]
            new_constr = _update_constr_bounds(new_bounds, constr, _error = _error)
            JuMP.owner_model(vref).constrs[cindex] = new_constr
        end
    end
    _validate_bounds(JuMP.owner_model(vref), new_bounds, _error = _error)
    # add the bounds
    parameter_bounds(vref).intervals[pref] = IntervalSet(lower, upper)
    # update status
    JuMP.owner_model(vref).has_hold_bounds = true
    if is_used(vref)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    end
    return
end

# TODO add parameter bound deletion


"""
    JuMP.set_name(vref::InfiniteVariableRef, root_name::String)

Extend [`JuMP.set_name`](@ref) to set names of infinite variables. Adds on to
`root_name` the ending `(prefs...)` where the parameter reference names are
listed in the same format as input in the parameter reference tuple.

**Example**
```julia
julia> name(vref)
old_name(t, x)

julia> set_name(vref, "new_name")

julia> name(vref)
new_name(t, x)
```
"""
function JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
    if length(root_name) == 0
        root_name = "noname"
    end
    # TODO do something about SparseAxisArrays (report array of values in order)
    # TODO list as vector and use ... like REPL if there are a lot
    prefs = parameter_refs(vref)
    param_names = _root_names(prefs)
    param_name_tuple = "("
    for i = 1:length(param_names)
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

# Make a variable reference
function _make_variable_ref(model::InfiniteModel, index::Int)::GeneralVariableRef
    if isa(model.vars[index], InfiniteVariable)
        return InfiniteVariableRef(model, index)
    elseif isa(model.vars[index], PointVariable)
        return PointVariableRef(model, index)
    else
        return HoldVariableRef(model, index)
    end
end

"""
    JuMP.variable_by_name(model::InfiniteModel,
                          name::String)::Union{GeneralVariableRef, Nothing}

Extend [`JuMP.variable_by_name`](@ref) for `InfiniteModel` objects. Return the
varaible reference assoociated with a variable name. Errors if multiple
variables have the same name. Returns nothing if no such name exists.

**Examples**
```julia
julia> variable_by_name(m, "var_name")
var_name

julia> variable_by_name(m, "fake_name")

```
"""
function JuMP.variable_by_name(model::InfiniteModel,
                               name::String)::Union{GeneralVariableRef, Nothing}
    if model.name_to_var === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_var = Dict{String, Int}()
        for (var, var_name) in model.var_to_name
            if haskey(model.name_to_var, var_name)
                # -1 is a special value that means this string does not map to
                # a unique variable name.
                model.name_to_var[var_name] = -1
            else
                model.name_to_var[var_name] = var
            end
        end
    end
    index = get(model.name_to_var, name, nothing)
    if index isa Nothing
        return nothing
    elseif index == -1
        error("Multiple variables have the name $name.")
    else
        return _make_variable_ref(model, index)
    end
    return
end

"""
    JuMP.all_variables(model::InfiniteModel)::Vector{GeneralVariableRef}

Extend [`JuMP.all_variables`](@ref) to return a list of all the variable
references associated with `model`.

**Examples**
```julia
julia> all_variables(m)
4-element Array{GeneralVariableRef,1}:
 ivar(test, θ)
 ivar2(test, x)
 name
 z
```
"""
function JuMP.all_variables(model::InfiniteModel)::Vector{GeneralVariableRef}
    vrefs_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model))
    indexes = sort([index for index in keys(model.vars)])
    counter = 1
    for index in indexes
        vrefs_list[counter] = _make_variable_ref(model, index)
        counter += 1
    end
    return vrefs_list
end
