################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::HoldVariableIndex
                               )::HoldVariableRef
    return HoldVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::VariableData{<:HoldVariable}
                          )::HoldVariableIndex
    return MOIUC.add_item(model.hold_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{HoldVariable}
    )::MOIUC.CleverDict{HoldVariableIndex, VariableData{HoldVariable{GeneralVariableRef}}}
    return model.hold_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(vref::HoldVariableRef
    )::MOIUC.CleverDict{HoldVariableIndex, VariableData{HoldVariable{GeneralVariableRef}}}
    return JuMP.owner_model(vref).hold_vars
end

# Extend _data_object
function _data_object(vref::HoldVariableRef)::VariableData{<:HoldVariable}
    return _data_dictionary(vref)[JuMP.index(vref)]
end

# Extend _core_variable_object
function _core_variable_object(vref::HoldVariableRef)::HoldVariable
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
# Check that parameter_bounds argument is valid
function _check_bounds(bounds::ParameterBounds; _error = error)::Nothing
    for (pref, set) in bounds
        # check that respects lower bound
        if JuMP.has_lower_bound(pref) && (JuMP.lower_bound(set) < JuMP.lower_bound(pref))
                _error("Specified parameter lower bound exceeds that defined " *
                       "for $pref.")
        end
        # check that respects upper bound
        if JuMP.has_upper_bound(pref) && (JuMP.upper_bound(set) > JuMP.upper_bound(pref))
                _error("Specified parameter upper bound exceeds that defined " *
                       "for $pref.")
        end
    end
    return
end

# Define _make_variable (used by JuMP.build_variable)
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Hold};
    parameter_bounds::ParameterBounds{GeneralVariableRef} = ParameterBounds(),
    extra_kw_args...)::HoldVariable
    # check for unneeded keywords
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with hold variables.")
    end
    # check that the bounds don't violate parameter domains
    _check_bounds(parameter_bounds)
    # make variable and return
    return HoldVariable(_make_float_info(info), parameter_bounds)
end

# Validate parameter bounds and add support(s) if needed
function _validate_bounds(model::InfiniteModel, bounds::ParameterBounds;
                          _error = error)::Nothing
    for (pref, set) in bounds
        # check validity
        JuMP.check_belongs_to_model(pref, model)
        # ensure has a support if a point constraint was given
        if JuMP.lower_bound(set) == JuMP.upper_bound(set) # TODO resolve for dependent parameters
            add_supports(pref, JuMP.lower_bound(set))
        end
    end
    return
end

# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(model::InfiniteModel,
                                      v::HoldVariable)::HoldVariableRef
    _validate_bounds(model, v.parameter_bounds)
    data_object = VariableData(v)
    vindex = _add_data_object(model, data_object)
    vref = HoldVariableRef(model, vindex)
    if length(v.parameter_bounds) != 0
        model.has_hold_bounds = true
    end
    return vref
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for hold variables
function _update_variable_info(vref::HoldVariableRef,
                               info::JuMP.VariableInfo)::Nothing
    bounds = _core_variable_object(vref).parameter_bounds
    _set_core_variable_object(vref, HoldVariable(info, bounds))
    return
end

################################################################################
#                           PARAMETER BOUND METHODS
################################################################################
#=
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
function _check_meas_bounds(bounds::ParameterBounds, data::AbstractMeasureData;
                            _error = error)
    if !measure_data_in_hold_bounds(data, bounds)
        _error("New bounds don't span existing dependent measure bounds.")
    end
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

## Update the variable bounds if it has any
# GeneralVariableRef
function _update_var_bounds(vref::GeneralVariableRef,
                            constr_bounds::ParameterBounds)
    return
end

# HoldVariableRef
function _update_var_bounds(vref::HoldVariableRef,
                            constr_bounds::ParameterBounds)
    if has_parameter_bounds(vref)
        _update_bounds(constr_bounds.intervals, parameter_bounds(vref).intervals)
    end
    return
end

# MeasureRef
function _update_var_bounds(mref::MeasureRef,
                            constr_bounds::ParameterBounds)
    vrefs = _all_function_variables(measure_function(mref))
    for vref in vrefs
        _update_var_bounds(vref, constr_bounds)
    end
    return
end

## Rebuild the constraint bounds (don't change in case of error)
# BoundedScalarConstraint
function _rebuild_constr_bounds(c::BoundedScalarConstraint,
                                var_bounds::ParameterBounds; _error = error)
    # prepare new constraint
    vrefs = _all_function_variables(c.func)
    c_new = BoundedScalarConstraint(c.func, c.set, copy(c.orig_bounds), c.orig_bounds)
    # look for bounded hold variables and update bounds
    for vref in vrefs
        _update_var_bounds(vref, c_new.bounds)
    end
    # check if the constraint bounds have and update if doesn't
    if length(c_new.bounds) == 0
        c_new = JuMP.ScalarConstraint(c.func, c.set)
    end
    return c_new
end

# ScalarConstraint
function _rebuild_constr_bounds(c::JuMP.ScalarConstraint,
                                var_bounds::ParameterBounds; _error = error)
    return BoundedScalarConstraint(c.func, c.set, var_bounds, ParameterBounds())
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
        cindices = Int[]
        if used_by_measure(vref)
            for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
                meas = JuMP.owner_model(vref).measures[mindex]
                _check_meas_bounds(bounds, meas.data, _error = _error)
                if used_by_constraint(MeasureRef(JuMP.owner_model(vref), mindex))
                    indices = JuMP.owner_model(vref).meas_to_constrs[mindex]
                    push!(cindices, indices...)
                end
            end
        end
        # set the new bounds
        _validate_bounds(JuMP.owner_model(vref), bounds, _error = _error)
        _update_variable_param_bounds(vref, bounds)
        # check and update dependent constraints
        if used_by_constraint(vref)
            union!(cindices, JuMP.owner_model(vref).var_to_constrs[JuMP.index(vref)])
        end
        for cindex in cindices
            constr = JuMP.owner_model(vref).constrs[cindex]
            new_constr = _rebuild_constr_bounds(constr, bounds, _error = _error)
            JuMP.owner_model(vref).constrs[cindex] = new_constr
        end
        # update status
        JuMP.owner_model(vref).has_hold_bounds = true
        if is_used(vref)
            set_optimizer_model_ready(JuMP.owner_model(vref), false)
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

"""
    delete_parameter_bound(vref::HoldVariableRef, pref::ParameterRef)

Delete the parameter bound of the hold variable `vref` associated with the
infinite parameter `pref` if `vref` has such a bound. Note that any other
parameter bounds will be unaffected. Any constraints that employ `vref` will
be updated accordingly.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, x[1:2] in [0, 10])
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> @hold_variable(model, z, parameter_bounds = (x in [0, 1]))
z

julia> delete_parameter_bound(z, x[2])

julia> parameter_bounds(z)
Subdomain bounds (1): x[1] ∈ [0, 1]
```
"""
function delete_parameter_bound(vref::HoldVariableRef, pref::ParameterRef)
    # get the current bounds
    bounds = parameter_bounds(vref)
    # check if there are bounds for pref and act accordingly
    if haskey(bounds.intervals, pref)
        delete!(bounds.intervals, pref)
        # check for dependent measures that are used by constraints
        meas_cindices = []
        if used_by_measure(vref)
            for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
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
                new_constr = _rebuild_constr_bounds(constr, bounds)
                JuMP.owner_model(vref).constrs[cindex] = new_constr
            end
        end
    end
    return
end

"""
    delete_parameter_bounds(vref::HoldVariableRef)

Delete all the parameter bounds of the hold variable `vref`. Any constraints
that employ `vref` will be updated accordingly.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, x[1:2] in [0, 10])
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> @hold_variable(model, z, parameter_bounds = (x in [0, 1]))
z

julia> delete_parameter_bounds(z)

julia> parameter_bounds(z)
Subdomain bounds (0):
```
"""
function delete_parameter_bounds(vref::HoldVariableRef)
    # get the current bounds
    bounds = parameter_bounds(vref)
    # check if there are bounds and act accordingly
    if length(bounds) > 0
        _update_variable_param_bounds(vref, ParameterBounds())
        # check for dependent measures that are used by constraints
        meas_cindices = []
        if used_by_measure(vref)
            for mindex in JuMP.owner_model(vref).var_to_meas[JuMP.index(vref)]
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
                new_constr = _rebuild_constr_bounds(constr, bounds)
                JuMP.owner_model(vref).constrs[cindex] = new_constr
            end
        end
    end
    return
end
=#

################################################################################
#                                 DELETION
################################################################################
