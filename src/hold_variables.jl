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
function _data_object(vref::HoldVariableRef
    )::VariableData{HoldVariable{GeneralVariableRef}}
  object = _get(_data_dictionary(vref), JuMP.index(vref), nothing)
  object === nothing && error("Invalid hold variable reference, cannot find " *
                        "corresponding variable in the model. This is likely " *
                        "caused by using the reference of a deleted variable.")
  return object
end

# Extend _core_variable_object
function _core_variable_object(vref::HoldVariableRef
    )::HoldVariable{GeneralVariableRef}
    return _data_object(vref).variable
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
# Check that parameter_bounds argument is valid
function _check_bounds(bounds::ParameterBounds{GeneralVariableRef};
                       _error = error)::Nothing
    depend_counter = Dict{DependentParameterRef, Int}()
    for (pref, set) in bounds
        # check that pref is an infinite parameter
        if !(_index_type(pref) <: InfiniteParameterIndex)
            _error("Can only specify infinite parameters for parameter bounds.")
        end
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
        # keep track of dependent parameters with equality conditions to ensure completeness
        if (_index_type(pref) == DependentParameterIndex) &&
           (JuMP.lower_bound(set) == JuMP.upper_bound(set))
            index = DependentParameterIndex(DependentParametersIndex(_raw_index(pref)), 1)
            dumby_pref = dispatch_variable_ref(pref.model, index)
            if haskey(depend_counter, dumby_pref)
                depend_counter[dumby_pref] += 1
            else
                depend_counter[dumby_pref] = 1
            end
        end
    end
    # check dimensions of dependent parameter equalities (if any are provided)
    for (pref, count) in depend_counter
        if _num_parameters(pref) != count
            _error("Cannot specify equality parameter bounds for a subset of " *
                   "dependent infinite parameters.")
        end
    end
    return
end

# Define _make_variable (used by JuMP.build_variable)
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Val{Hold};
    parameter_bounds::ParameterBounds{GeneralVariableRef} = ParameterBounds(),
    extra_kw_args...)::HoldVariable{GeneralVariableRef}
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
function _validate_bounds(model::InfiniteModel,
                          bounds::ParameterBounds{GeneralVariableRef};
                          _error = error)::Nothing
    depend_supps = Dict{DependentParametersIndex, Matrix{Float64}}()
    for (pref, set) in bounds
        # check validity
        JuMP.check_belongs_to_model(pref, model)
        # ensure has a support if a point constraint was given
        if (JuMP.lower_bound(set) == JuMP.upper_bound(set))
            if _index_type(pref) == IndependentParameterIndex
                # label will be UserDefined
                add_supports(pref, JuMP.lower_bound(set), check = false)
            else
                index = DependentParametersIndex(_raw_index(pref))
                if !haskey(depend_supps, index)
                    dumby_pref = dispatch_variable_ref(model, DependentParameterIndex(index, 1))
                    depend_supps[index] = Matrix{Float64}(undef, _num_parameters(dumby_pref), 1)
                end
                depend_supps[index][_param_index(pref)] = JuMP.lower_bound(set)
            end
        end
    end
    # add dependent supports if any are given
    for (index, supp) in depend_supps
        prefs = [dispatch_variable_ref(model, DependentParameterIndex(index, i))
                 for i in 1:length(supp)]
        add_supports(prefs, supp, check = false) # label will be UserDefined
    end
    return
end

# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(model::InfiniteModel,
                                      v::HoldVariable,
                                      name::String)::HoldVariableRef
    _validate_bounds(model, v.parameter_bounds)
    data_object = VariableData(v, name)
    vindex = _add_data_object(model, data_object)
    vref = HoldVariableRef(model, vindex)
    if !isempty(v.parameter_bounds)
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
"""
    parameter_bounds(vref::HoldVariableRef)::ParameterBounds

Return the [`ParameterBounds`](@ref) object associated with the hold variable
`vref`. It contains a dictionary where each key is a parameter reference
which points to an `IntervalSet` that that defines a sub-domain for `vref`
relative to that parameter reference.

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
function parameter_bounds(vref::HoldVariableRef
    )::ParameterBounds{GeneralVariableRef}
    return _core_variable_object(vref).parameter_bounds
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
    return !isempty(parameter_bounds(vref))
end

# Other variable types
function has_parameter_bounds(vref::DispatchVariableRef)::Bool
    return false
end

# Internal function used to change the parameter bounds of a hold variable
function _update_variable_param_bounds(vref::HoldVariableRef,
                                       bounds::ParameterBounds{GeneralVariableRef}
                                       )::Nothing
    info = _variable_info(vref)
    _set_core_variable_object(vref, HoldVariable(info, bounds))
    JuMP.owner_model(vref).has_hold_bounds = true
    if is_used(vref)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    end
    return
end

# Check that the bounds dictionary is compadable with existing dependent measures
function _check_meas_bounds(bounds::ParameterBounds{GeneralVariableRef},
                            data::AbstractMeasureData;
                            _error::Function = error
                            )::Nothing
    if !measure_data_in_hold_bounds(data, bounds)
        _error("New bounds don't span existing dependent measure bounds.")
    end
    return
end

# Update the current bounds to overlap with the new bounds if possible
function _update_bounds(bounds1::ParameterBounds{GeneralVariableRef},
                        bounds2::ParameterBounds{GeneralVariableRef};
                        _error::Function = error
                        )::Nothing
    # check each new bound
    for (pref, set) in bounds2
        # we have a new bound
        if !haskey(bounds1, pref)
            bounds1[pref] = set
        # the previous set and the new one do not overlap
        elseif (JuMP.lower_bound(set) > JuMP.upper_bound(bounds1[pref])) ||
               (JuMP.upper_bound(set) < JuMP.lower_bound(bounds1[pref]))
            _error("Sub-domains of constraint and/or hold variable(s) do not" *
                   " overlap. Consider changing the parameter bounds of the" *
                   " constraint and/or hold variable(s).")
        # we have an existing bound
        else
            # we have a new stricter lower bound to update with
            if JuMP.lower_bound(set) > JuMP.lower_bound(bounds1[pref])
                bounds1[pref] = IntervalSet(JuMP.lower_bound(set),
                                            JuMP.upper_bound(bounds1[pref]))
            end
            # we have a new stricter upper bound to update with
            if JuMP.upper_bound(set) < JuMP.upper_bound(bounds1[pref])
                bounds1[pref] = IntervalSet(JuMP.lower_bound(bounds1[pref]),
                                            JuMP.upper_bound(set))
            end
        end
    end
    return
end

## Update the variable bounds if it has any
# DispatchVariableRef
function _update_var_bounds(vref::DispatchVariableRef,
                            constr_bounds::ParameterBounds{GeneralVariableRef}
                            )::Nothing
    return
end

# HoldVariableRef
function _update_var_bounds(vref::HoldVariableRef,
                            constr_bounds::ParameterBounds{GeneralVariableRef}
                            )::Nothing
    if has_parameter_bounds(vref)
        _update_bounds(constr_bounds, parameter_bounds(vref))
    end
    return
end

# MeasureRef
function _update_var_bounds(mref::MeasureRef,
                            constr_bounds::ParameterBounds{GeneralVariableRef}
                            )::Nothing
    vrefs = _all_function_variables(measure_function(mref))
    for vref in vrefs
        _update_var_bounds(vref, constr_bounds)
    end
    return
end

# GeneralVariableRef
function _update_var_bounds(vref::GeneralVariableRef,
                            constr_bounds::ParameterBounds{GeneralVariableRef}
                            )::Nothing
    return _update_var_bounds(dispatch_variable_ref(vref), constr_bounds)
end

## Rebuild the constraint bounds (don't change in case of error)
# BoundedScalarConstraint
function _rebuild_constr_bounds(c::BoundedScalarConstraint,
                                var_bounds::ParameterBounds{GeneralVariableRef};
                                _error::Function = error
                                )::JuMP.AbstractConstraint
    # prepare new constraint
    func = JuMP.jump_function(c)
    set = JuMP.moi_set(c)
    vrefs = _all_function_variables(func)
    orig_bounds = original_parameter_bounds(c)
    c_new = BoundedScalarConstraint(func, set, copy(orig_bounds), orig_bounds)
    # look for bounded hold variables and update bounds
    for vref in vrefs
        _update_var_bounds(vref, parameter_bounds(c_new))
    end
    # check if the constraint has bounds and update if doesn't
    if isempty(parameter_bounds(c_new))
        c_new = JuMP.ScalarConstraint(func, set)
    end
    return c_new
end

# ScalarConstraint
function _rebuild_constr_bounds(c::JuMP.ScalarConstraint,
                                var_bounds::ParameterBounds{GeneralVariableRef};
                                _error::Function = error
                                )::BoundedScalarConstraint
    func = JuMP.jump_function(c)
    set = JuMP.moi_set(c)
    return BoundedScalarConstraint(func, set, var_bounds, ParameterBounds())
end

"""
    set_parameter_bounds(vref::HoldVariableRef,
                         bounds::ParameterBounds{GeneralVariableRef};
                         [force::Bool = false])::Nothing

Specify a new dictionary of parameter bounds `bounds` for the hold variable `vref`.
These are stored in a [`ParameterBounds`](@ref) object which contains a dictionary.
Note the dictionary keys must be infinite parameter references and the values
must be `IntervalSet`s that indicate a particular sub-domain for which `vref`
is defined. This is meant to be primarily used by
[`@set_parameter_bounds`](@ref) which provides a more intuitive syntax.

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
function set_parameter_bounds(vref::HoldVariableRef,
                              bounds::ParameterBounds{GeneralVariableRef};
                              force::Bool = false, _error::Function = error
                              )::Nothing
    if has_parameter_bounds(vref) && !force
        _error("$vref already has parameter bounds. Consider adding more using " *
               "`add_parameter_bounds` or overwriting them by setting " *
               "the keyword argument `force = true`")
    else
        # check that bounds are valid
        _check_bounds(bounds, _error = _error)
        # check dependent measures
        model = JuMP.owner_model(vref)
        cindices = ConstraintIndex[]
        for mindex in _measure_dependencies(vref)
            mref = dispatch_variable_ref(model, mindex)
            _check_meas_bounds(bounds, _core_variable_object(mref).data,
                               _error = _error)
            append!(cindices, _constraint_dependencies(mref))
        end
        # set the new bounds
        _validate_bounds(model, bounds, _error = _error)
        _update_variable_param_bounds(vref, bounds)
        # check and update dependent constraints
        union!(cindices, _constraint_dependencies(vref))
        for cindex in cindices
            cref = _temp_constraint_ref(model, cindex)
            constr = _core_constraint_object(cref)
            new_constr = _rebuild_constr_bounds(constr, bounds, _error = _error)
            _set_core_constraint_object(cref, new_constr)
        end
    end
    return
end

## Check and update the constraint bounds (don't change in case of error)
# BoundedScalarConstraint
function _update_constr_bounds(bounds::ParameterBounds{GeneralVariableRef},
                               c::BoundedScalarConstraint;
                               _error::Function = error
                               )::BoundedScalarConstraint
    new_bounds = copy(parameter_bounds(c))
    _update_bounds(new_bounds, bounds, _error = _error)
    func = JuMP.jump_function(c)
    set = JuMP.moi_set(c)
    orig_bounds = original_parameter_bounds(c)
    return BoundedScalarConstraint(func, set, new_bounds, orig_bounds)
end

# ScalarConstraint
function _update_constr_bounds(bounds::ParameterBounds{GeneralVariableRef},
                               c::JuMP.ScalarConstraint;
                               _error::Function = error
                               )::BoundedScalarConstraint
    func = JuMP.jump_function(c)
    set = JuMP.moi_set(c)
    return BoundedScalarConstraint(func, set, bounds, ParameterBounds())
end

"""
    add_parameter_bounds(vref::HoldVariableRef,
                         new_bounds::ParameterBounds{GeneralVariableRef}
                         )::Nothing

Add additional parameter bounds to `vref` such that it is defined over the
sub-domain based on the intersection of the existing bounds and `new_bounds`.
This is primarily meant to be used by [`@add_parameter_bounds`](@ref).

```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @hold_variable(model, vref)
vref

julia> add_parameter_bounds(vref, ParameterBounds(t => IntervalSet(0, 2)))

julia> parameter_bounds(vref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function add_parameter_bounds(vref::HoldVariableRef,
                              new_bounds::ParameterBounds{GeneralVariableRef};
                              _error::Function = error
                              )::Nothing
    # if overriding any existing parameter bounds, dispatch to set new ones
    if any(haskey(parameter_bounds(vref), pref) for pref in keys(new_bounds))
        merged_bounds = merge(parameter_bounds(vref), new_bounds)
        return set_parameter_bounds(vref, merged_bounds, force = true,
                                    _error = _error)
    end
    # check the new bounds
    _check_bounds(new_bounds, _error = _error)
    # check dependent measures
    model = JuMP.owner_model(vref)
    cindices = ConstraintIndex[]
    for mindex in _measure_dependencies(vref)
        mref = dispatch_variable_ref(model, mindex)
        _check_meas_bounds(new_bounds, _core_variable_object(mref).data,
                           _error = _error)
        append!(cindices, _constraint_dependencies(mref))
    end
    # check and update dependent constraints
    union!(cindices, _constraint_dependencies(vref))
    for cindex in cindices
        cref = _temp_constraint_ref(model, cindex)
        constr = _core_constraint_object(cref)
        new_constr = _update_constr_bounds(new_bounds, constr, _error = _error)
        _set_core_constraint_object(cref, new_constr)
    end
    _validate_bounds(model, new_bounds, _error = _error)
    # add the bounds
    merge!(parameter_bounds(vref), new_bounds)
    # update status
    JuMP.owner_model(vref).has_hold_bounds = true
    if is_used(vref)
        set_optimizer_model_ready(model, false)
    end
    return
end

"""
    delete_parameter_bounds(vref::HoldVariableRef)::Nothing

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
function delete_parameter_bounds(vref::HoldVariableRef)::Nothing
    # get the current bounds
    bounds = parameter_bounds(vref)
    # check if there are bounds and act accordingly
    if !isempty(bounds)
        _update_variable_param_bounds(vref, ParameterBounds())
        # check for dependent measures that are used by constraints
        model = JuMP.owner_model(vref)
        cindices = ConstraintIndex[]
        for mindex in _measure_dependencies(vref)
            mref = dispatch_variable_ref(model, mindex)
            _check_meas_bounds(bounds, _core_variable_object(mref).data,
                               _error = error)
            append!(cindices, _constraint_dependencies(mref))
        end
        # check and update dependent constraints
        union!(cindices, _constraint_dependencies(vref))
        for cindex in cindices
            cref = _temp_constraint_ref(model, cindex)
            constr = _core_constraint_object(cref)
            new_constr = _rebuild_constr_bounds(constr, bounds, _error = error)
            _set_core_constraint_object(cref, new_constr)
        end
        # update status
        if is_used(vref)
            set_optimizer_model_ready(model, false)
        end
    end
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::HoldVariableRef)::Nothing
    # remove variable info constraints associated with vref
    _delete_info_constraints(vref)
    # update dependent constraint bounds by deleting the variable bounds
    delete_parameter_bounds(vref)
    return
end
