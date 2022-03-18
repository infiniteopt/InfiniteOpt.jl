################################################################################
#                               VARIABLE DEFINITION
################################################################################
"""
    InfOptVariableType

An abstract `DataType` for variable type objects used to create `InfiniteOpt` 
variables via `JuMP.@variable`.
"""
abstract type InfOptVariableType end 

# helper function for setting the info constraints (TODO add test)
function _set_info_constraints(
    info::JuMP.VariableInfo, 
    gvref::GeneralVariableRef, 
    dvref::DispatchVariableRef
    )::Nothing
    model = JuMP.owner_model(gvref)
    if info.has_lb
        newset = MOI.GreaterThan(info.lower_bound)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_lower_bound_index(dvref, JuMP.index(cref))
    end
    if info.has_ub
        newset = MOI.LessThan(info.upper_bound)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_upper_bound_index(dvref, JuMP.index(cref))
    end
    if info.has_fix
        newset = MOI.EqualTo(info.fixed_value)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_fix_index(dvref, JuMP.index(cref))
    end
    if info.binary
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.ZeroOne()),
                                   is_info_constr = true)
        _set_binary_index(dvref, JuMP.index(cref))
    elseif info.integer
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.Integer()),
                                   is_info_constr = true)
        _set_integer_index(dvref, JuMP.index(cref))
    end
    return
end

# Add VariableConstrainedOnCreation (TODO change paradigm to leverage MOI.add_constrained_variable)
function JuMP.add_variable(
    model::InfiniteModel, 
    variable::JuMP.VariableConstrainedOnCreation,
    name::String
    )::GeneralVariableRef
    vref = JuMP.add_variable(model, variable.scalar_variable, name)
    con = JuMP.ScalarConstraint(vref, variable.set)
    JuMP.add_constraint(model, con)
    return vref
end

# Add VariablesConstrainedOnCreation (TODO change paradigm to leverage MOI.add_constrained_variables)
function JuMP.add_variable(
    model::InfiniteModel, 
    variables::JuMP.VariablesConstrainedOnCreation,
    names::AbstractArray{<:String}
    )
    shape = variables.shape
    vector_names = JuMP.vectorize(names, shape)
    vrefs = JuMP.add_variable.(model, variables.scalar_variables, vector_names)
    shaped_vrefs = JuMP.reshape_vector(vrefs, shape)
    shaped_set = JuMP.reshape_set(variables.set, shape)
    con = JuMP.build_constraint(error, shaped_vrefs, shaped_set)
    JuMP.add_constraint(model, con)
    return shaped_vrefs
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _measure_dependencies
function _measure_dependencies(vref::DecisionVariableRef)::Vector{MeasureIndex}
    return _data_object(vref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(vref::DecisionVariableRef)::Vector{InfOptConstraintIndex}
    return _data_object(vref).constraint_indices
end

"""
    used_by_measure(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a measure.

**Example**
```julia-repl
julia> used_by_measure(vref)
true
```
"""
function used_by_measure(vref::DecisionVariableRef)::Bool
    return !isempty(_measure_dependencies(vref))
end

"""
    used_by_constraint(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a constraint.

**Example**
```julia-repl
julia> used_by_constraint(vref)
false
```
"""
function used_by_constraint(vref::DecisionVariableRef)::Bool
    return !isempty(_constraint_dependencies(vref))
end

"""
    used_by_objective(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by the objective.

**Example**
```julia-repl
julia> used_by_objective(vref)
true
```
"""
function used_by_objective(vref::DecisionVariableRef)::Bool
    return _data_object(vref).in_objective
end

"""
    is_used(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia-repl
julia> is_used(vref)
true
```
"""
function is_used(vref::DecisionVariableRef)::Bool
    return used_by_measure(vref) || used_by_constraint(vref) || used_by_objective(vref)
end

################################################################################
#                                VARIABLE NAMING
################################################################################
"""
    JuMP.name(vref::DecisionVariableRef)::String

Extend `JuMP.name` to return the names of `InfiniteOpt` variables.

**Example**
```julia-repl
julia> name(vref)
"var_name"
```
"""
function JuMP.name(vref::DecisionVariableRef)::String
    object = Base.get(_data_dictionary(vref), JuMP.index(vref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(vref::DecisionVariableRef, name::String)::Nothing

Extend `JuMP.set_name` to set names of decision variables.

**Example**
```julia-repl
julia> set_name(vref, "var_name")

julia> name(vref)
"var_name"
```
"""
function JuMP.set_name(vref::DecisionVariableRef, name::String)::Nothing
    _data_object(vref).name = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

# Make a variable reference
function _make_variable_ref(
    model::InfiniteModel, 
    index::ObjectIndex
    )::GeneralVariableRef
    return GeneralVariableRef(model, index.value, typeof(index))
end

# Get the name_to_var Dictionary
function _var_name_dict(
    model::InfiniteModel
    )::Union{Nothing, Dict{String, ObjectIndex}}
    return model.name_to_var
end

# Update name_to_var
function _update_var_name_dict(
    model::InfiniteModel, 
    var_dict::MOIUC.CleverDict
    )::Nothing
    name_dict = _var_name_dict(model)
    for (index, data_object) in var_dict
        var_name = data_object.name
        if haskey(name_dict, var_name)
            name_dict[var_name] = FiniteVariableIndex(-1) # dumby value
        else
            name_dict[var_name] = index
        end
    end
    model.name_to_var = name_dict
    return
end

"""
    JuMP.variable_by_name(model::InfiniteModel,
                          name::String)::Union{GeneralVariableRef, Nothing}

Extend `JuMP.variable_by_name` for `InfiniteModel` objects. Return the variable 
reference assoociated with a variable name. Errors if multiple variables have the 
same name. Returns nothing if no such name exists.

**Examples**
```julia-repl
julia> variable_by_name(m, "var_name")
var_name

julia> variable_by_name(m, "fake_name")

```
"""
function JuMP.variable_by_name(
    model::InfiniteModel,
    name::String
    )::Union{GeneralVariableRef, Nothing}
    if isnothing(_var_name_dict(model))
        model.name_to_var = Dict{String, ObjectIndex}()
        _update_var_name_dict(model, model.infinite_vars)
        _update_var_name_dict(model, model.semi_infinite_vars)
        _update_var_name_dict(model, model.point_vars)
        _update_var_name_dict(model, model.finite_vars)
    end
    index = Base.get(_var_name_dict(model), name, nothing)
    if isnothing(index)
        return nothing
    elseif index == FiniteVariableIndex(-1)
        error("Multiple variables have the name $name.")
    else
        return _make_variable_ref(model, index)
    end
end

################################################################################
#                              PARAMETER REFERENCES
################################################################################
# Extend parameter_refs for variables (this serves as a fallback for finite types)
function parameter_refs(vref::FiniteRef)::Tuple
    return ()
end

################################################################################
#                         VARIABLE OBJECT MODIFICATION
################################################################################
# Extend _set_core_variable_object
function _set_core_variable_object(
    vref::DecisionVariableRef,
    object::JuMP.AbstractVariable
    )::Nothing
    _data_object(vref).variable = object
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Get info
function _variable_info(vref::UserDecisionVariableRef)::JuMP.VariableInfo
    return _core_variable_object(vref).info
end

"""
    JuMP.has_lower_bound(vref::UserDecisionVariableRef)::Bool

Extend `JuMP.has_lower_bound` to return a `Bool` whether an `InfiniteOpt` 
variable has a lower bound.

**Example**
```julia-repl
julia> has_lower_bound(vref)
true
```
"""
function JuMP.has_lower_bound(vref::UserDecisionVariableRef)::Bool 
    return _variable_info(vref).has_lb
end

"""
    JuMP.lower_bound(vref::UserDecisionVariableRef)::Float64

Extend `JuMP.lower_bound` to return the lower bound of an `InfiniteOpt` variable. 
Errors if `vref` doesn't have a lower bound.

**Example**
```julia-repl
julia> lower_bound(vref)
0.0
```
"""
function JuMP.lower_bound(vref::UserDecisionVariableRef)::Float64
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return _variable_info(vref).lower_bound
end

# Extend to return the index of the lower bound constraint associated with `vref`.
function _lower_bound_index(vref::UserDecisionVariableRef)::InfOptConstraintIndex
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return _data_object(vref).lower_bound_index
end

# Extend to specify the index `cindex` of a lower bound constraint for `vref`.
function _set_lower_bound_index(
    vref::UserDecisionVariableRef,
    cindex::Union{InfOptConstraintIndex, Nothing}
    )::Nothing
    _data_object(vref).lower_bound_index = cindex
    return
end

"""
    JuMP.set_lower_bound(vref::UserDecisionVariableRef, lower::Real)::Nothing

Extend `JuMP.set_lower_bound` to specify the lower bound of an `InfiniteOpt` 
variable `vref`. Errors if `vref` is fixed.

**Example**
```julia-repl
julia> set_lower_bound(vref, -1)

julia> lower_bound(vref)
-1.0
```
"""
function JuMP.set_lower_bound(
    vref::UserDecisionVariableRef,
    lower::Real
    )::Nothing
    newset = MOI.GreaterThan(convert(Float64, lower))
    model = JuMP.owner_model(vref)
    gvref = _make_variable_ref(model, JuMP.index(vref))
    new_constr = JuMP.ScalarConstraint(gvref, newset)
    if JuMP.has_lower_bound(vref)
        cindex = _lower_bound_index(vref)
        cref = _make_constraint_ref(model, cindex)
        _set_core_constraint_object(cref, new_constr)
        set_optimizer_model_ready(model, false)
    else
        @assert !JuMP.is_fixed(vref) "$vref is fixed, cannot set lower bound."
        cref = JuMP.add_constraint(model, new_constr, is_info_constr = true)
        _set_lower_bound_index(vref, JuMP.index(cref))
    end
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(true, convert(Float64, lower),
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.LowerBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef

Extend `JuMP.LowerBoundRef` to extract a constraint reference for the lower 
bound of `vref`.

**Example**
```julia-repl
var ≥ 0.0
```
"""
function JuMP.LowerBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = _lower_bound_index(vref)
    model = JuMP.owner_model(vref)
    return _make_constraint_ref(model, cindex)
end

"""
    JuMP.delete_lower_bound(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.delete_lower_bound` to delete lower bound of `vref`. Errors if it 
doesn't have a lower bound.

**Example**
```julia-repl
julia> delete_lower_bound(vref)

julia> has_lower_bound(vref)
false
```
"""
function JuMP.delete_lower_bound(vref::UserDecisionVariableRef)::Nothing
    JuMP.delete(JuMP.owner_model(vref), JuMP.LowerBoundRef(vref))
    _set_lower_bound_index(vref, nothing)
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(false, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.has_upper_bound(vref::UserDecisionVariableRef)::Bool

Extend `JuMP.has_upper_bound` to return a `Bool` whether an `InfiniteOpt` 
variable has an upper bound.

**Example**
```julia-repl
julia> has_upper_bound(vref)
true
```
"""
function JuMP.has_upper_bound(vref::UserDecisionVariableRef)::Bool 
    return _variable_info(vref).has_ub
end

"""
    JuMP.upper_bound(vref::UserDecisionVariableRef)::Float64

Extend `JuMP.upper_bound` to return the upper bound of an `InfiniteOpt` variable. 
Errors if `vref` doesn't have a upper bound.

**Example**
```julia-repl
julia> upper_bound(vref)
0.0
```
"""
function JuMP.upper_bound(vref::UserDecisionVariableRef)::Float64
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return _variable_info(vref).upper_bound
end

# Extend to return the index of the upper bound constraint associated with `vref`.
function _upper_bound_index(vref::UserDecisionVariableRef)::InfOptConstraintIndex
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return _data_object(vref).upper_bound_index
end

# Extend to specify the index `cindex` of a upper bound constraint for `vref`
function _set_upper_bound_index(
    vref::UserDecisionVariableRef,
    cindex::Union{InfOptConstraintIndex, Nothing}
    )::Nothing
    _data_object(vref).upper_bound_index = cindex
    return
end

"""
    JuMP.set_upper_bound(vref::UserDecisionVariableRef, upper::Real)::Nothing

Extend `JuMP.set_upper_bound` to specify the upper bound of an `InfiniteOpt` 
variable `vref`. Errors if `vref` is fixed.

**Example**
```julia-repl
julia> set_upper_bound(vref, 1)

julia> upper_bound(vref)
1.0
```
"""
function JuMP.set_upper_bound(
    vref::UserDecisionVariableRef,
    upper::Real
    )::Nothing
    newset = MOI.LessThan(convert(Float64, upper))
    model = JuMP.owner_model(vref)
    gvref = _make_variable_ref(model, JuMP.index(vref))
    new_constr = JuMP.ScalarConstraint(gvref, newset)
    if JuMP.has_upper_bound(vref)
        cindex = _upper_bound_index(vref)
        cref = _make_constraint_ref(model, cindex)
        _set_core_constraint_object(cref, new_constr)
        set_optimizer_model_ready(model, false)
    else
        @assert !JuMP.is_fixed(vref) "$vref is fixed, cannot set upper bound."
        cref = JuMP.add_constraint(model, new_constr, is_info_constr = true)
        _set_upper_bound_index(vref, JuMP.index(cref))
    end
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  true, convert(Float64, upper),
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.UpperBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef

Extend `JuMP.UpperBoundRef` to extract a constraint reference for the upper 
bound of `vref`.

**Example**
```julia-repl
julia> cref = UpperBoundRef(vref)
var ≤ 1.0
```
"""
function JuMP.UpperBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = _upper_bound_index(vref)
    model = JuMP.owner_model(vref)
    return _make_constraint_ref(model, cindex)
end

"""
    JuMP.delete_upper_bound(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.delete_upper_bound` to delete the upper bound of `vref`. Errors if 
it doesn't have an upper bound.

**Example**
```julia-repl
julia> delete_upper_bound(vref)

julia> has_upper_bound(vref)
false
```
"""
function JuMP.delete_upper_bound(vref::UserDecisionVariableRef)::Nothing
    JuMP.delete(JuMP.owner_model(vref), JuMP.UpperBoundRef(vref))
    _set_upper_bound_index(vref, nothing)
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  false, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.is_fixed(vref::UserDecisionVariableRef)::Bool

Extend `JuMP.is_fixed` to return `Bool` whether an `InfiniteOpt` variable is 
fixed.

**Example**
```julia-repl
julia> is_fixed(vref)
true
```
"""
JuMP.is_fixed(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).has_fix

"""
    JuMP.fix_value(vref::UserDecisionVariableRef)::Float64

Extend `JuMP.fix_value` to return the fix value of an `InfiniteOpt` variable. 
Errors if variable is not fixed.

**Example**
```julia-repl
julia> fix_value(vref)
0.0
```
"""
function JuMP.fix_value(vref::UserDecisionVariableRef)::Float64
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return _variable_info(vref).fixed_value
end

# Extend to return the index of the fix constraint associated with `vref`.
function _fix_index(vref::UserDecisionVariableRef)::InfOptConstraintIndex
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return _data_object(vref).fix_index
end

# Extend to set the index of the fix constraintassociated with `vref`.
function _set_fix_index(
    vref::UserDecisionVariableRef,
    cindex::Union{InfOptConstraintIndex, Nothing}
    )::Nothing
    _data_object(vref).fix_index = cindex
    return
end

"""
    JuMP.fix(vref::UserDecisionVariableRef, value::Real;
             force::Bool = false)::Nothing

Extend `JuMP.fix` to fix the value of an `InfiniteOpt` variable. Errors if 
variable has a lower and/or an upper bound(s) unless `force = true`.

**Examples**
```julia-repl
julia> fix(vref, 3)

julia> fix_value(vref)
3.0

julia> fix(vref2, 2, force = true)

julia> fix_value(vref2)
2.0
```
"""
function JuMP.fix(
    vref::UserDecisionVariableRef, 
    value::Real;
    force::Bool = false
    )::Nothing
    new_set = MOI.EqualTo(convert(Float64, value))
    model = JuMP.owner_model(vref)
    gvref = _make_variable_ref(model, JuMP.index(vref))
    new_constr = JuMP.ScalarConstraint(gvref, new_set)
    if JuMP.is_fixed(vref)  # Update existing fixing constraint.
        cindex = _fix_index(vref)
        cref = _make_constraint_ref(model, cindex)
        _set_core_constraint_object(cref, new_constr)
        set_optimizer_model_ready(model, false)
    else  # Add a new fixing constraint.
        if  JuMP.has_upper_bound(vref) ||  JuMP.has_lower_bound(vref)
            if !force
                error("Unable to fix $(vref) to $(value) because it has " *
                      "existing variable bounds. Consider calling " *
                      "`JuMP.fix(variable, value; force=true)` which will " *
                      "delete existing bounds before fixing the variable.")
            end
            if  JuMP.has_upper_bound(vref)
                 JuMP.delete_upper_bound(vref)
            end
            if  JuMP.has_lower_bound(vref)
                 JuMP.delete_lower_bound(vref)
            end
        end
        cref = JuMP.add_constraint(model, new_constr, is_info_constr = true)
        _set_fix_index(vref, JuMP.index(cref))
    end
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(false, info.lower_bound,
                                                  false, info.upper_bound,
                                                  true, convert(Float64, value),
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.FixRef(vref::UserDecisionVariableRef)::InfOptConstraintRef

Extend `JuMP.FixRef` to return the constraint reference of the fix constraint 
associated with `vref`. Errors `vref` is not fixed.

**Examples**
```julia-repl
julia> cref = FixRef(vref)
var = 1.0
```
"""
function JuMP.FixRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = _fix_index(vref)
    model = JuMP.owner_model(vref)
    return _make_constraint_ref(model, cindex)
end

"""
    JuMP.unfix(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.unfix` to unfix `vref`. Errors if it is not fixed.

**Example**
```julia-repl
julia> unfix(vref)

julia> is_fixed(vref)
false
```
"""
function JuMP.unfix(vref::UserDecisionVariableRef)::Nothing
    JuMP.delete(JuMP.owner_model(vref), JuMP.FixRef(vref))
    _set_fix_index(vref, nothing)
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  false, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.start_value(vref::UserDecisionVariableRef)::Union{Nothing, Float64}

Extend `JuMP.start_value` to return starting value of `InfiniteOpt` variable if 
it has one. Returns `nothing` otherwise.

**Example**
```julia-repl
julia> start_value(vref)
0.0
```
"""
function JuMP.start_value(vref::UserDecisionVariableRef)::Union{Nothing, Float64}
    if _variable_info(vref).has_start
        return _variable_info(vref).start
    else
        return
    end
end

"""
    JuMP.set_start_value(vref::UserDecisionVariableRef, value::Real)::Nothing

Extend `JuMP.set_start_value` to specify the start value of `InfiniteOpt` 
variables.

**Example**
```julia-repl
julia> set_start_value(vref, 1)

julia> start_value(vref)
1.0
```
"""
function JuMP.set_start_value(
    vref::UserDecisionVariableRef,
    value::Real
    )::Nothing
    info = _variable_info(vref)
    set_optimizer_model_ready(JuMP.owner_model(vref), false)
    _update_variable_info(vref,
                          JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                            info.has_ub, info.upper_bound,
                                            info.has_fix, info.fixed_value,
                                            true, Float64(value),
                                            info.binary, info.integer))
    return
end

"""
    JuMP.is_binary(vref::UserDecisionVariableRef)::Bool

Extend `JuMP.is_binary` to return `Bool` whether an `InfiniteOpt` variable is 
binary.

**Example**
```julia-repl
true
```
"""
JuMP.is_binary(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).binary

# Extend to return the index of the binary constraint associated with `vref`.
function _binary_index(vref::UserDecisionVariableRef)::InfOptConstraintIndex
    if !JuMP.is_binary(vref)
        error("Variable $(vref) is not binary.")
    end
    return _data_object(vref).zero_one_index
end

# Extend to specify the index of the binary constraint associated with `vref`.
function _set_binary_index(
    vref::UserDecisionVariableRef,
    cindex::Union{InfOptConstraintIndex, Nothing}
    )::Nothing
    _data_object(vref).zero_one_index = cindex
    return
end

"""
    JuMP.set_binary(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.set_binary` to specify an `InfiniteOpt` variable as a binary 
variable. Errors if `vref` is an integer variable.

**Example**
```julia-repl
julia> set_binary(vref)

julia> is_binary(vref)
true
```
"""
function JuMP.set_binary(vref::UserDecisionVariableRef)::Nothing
    if JuMP.is_binary(vref)
        return
    elseif JuMP.is_integer(vref)
        error("Cannot set the variable_ref $(vref) to binary as it " *
              "is already integer.")
    end
    gvref = _make_variable_ref(JuMP.owner_model(vref), JuMP.index(vref))
    cref = JuMP.add_constraint(JuMP.owner_model(vref),
                               JuMP.ScalarConstraint(gvref, MOI.ZeroOne()),
                               is_info_constr = true)
    _set_binary_index(vref, JuMP.index(cref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  true, info.integer))
    return
end

"""
    JuMP.BinaryRef(vref::UserDecisionVariableRef)::InfOptConstraintRef

Extend `JuMP.BinaryRef` to return a constraint reference to the constraint 
constrainting `vref` to be binary. Errors if one does not exist.

**Example**
```julia-repl
julia> cref = BinaryRef(vref)
var binary
```
"""
function JuMP.BinaryRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = _binary_index(vref)
    model = JuMP.owner_model(vref)
    return _make_constraint_ref(model, cindex)
end

"""
    JuMP.unset_binary(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.unset_binary` to unset `vref` as a binary variable. Errors if it is 
not binary.

```julia-repl
julia> unset_binary(vref)

julia> is_binary(vref)
false
```
"""
function JuMP.unset_binary(vref::UserDecisionVariableRef)::Nothing
    JuMP.delete(JuMP.owner_model(vref), JuMP.BinaryRef(vref))
    _set_binary_index(vref, nothing)
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  false, info.integer))
    return
end

"""
    JuMP.is_integer(vref::UserDecisionVariableRef)::Bool

Extend `JuMP.is_integer` to return `Bool` whether an `InfiniteOpt` variable is 
integer.

**Example**
```julia-repl
julia> is_integer(vref)
true
```
"""
function JuMP.is_integer(vref::UserDecisionVariableRef)::Bool 
    return _variable_info(vref).integer
end

# Extend to return the index of the integer constraintassociated with `vref`.
function _integer_index(vref::UserDecisionVariableRef)::InfOptConstraintIndex
    if !JuMP.is_integer(vref)
        error("Variable $(vref) is not an integer.")
    end
    return _data_object(vref).integrality_index
end

# Extend to specify the index of the integer constraint associated with `vref`.
function _set_integer_index(
    vref::UserDecisionVariableRef,
    cindex::Union{InfOptConstraintIndex, Nothing}
    )::Nothing
    _data_object(vref).integrality_index = cindex
    return
end

"""
    JuMP.set_integer(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.set_integer` to specify an `InfiniteOpt` variable as a integer 
variable. Errors if `vref` is an binary variable.

**Example**
```julia-repl
julia> set_integer(vref)

julia> is_integer(vref)
true
```
"""
function JuMP.set_integer(vref::UserDecisionVariableRef)::Nothing
    if JuMP.is_integer(vref)
        return
    elseif JuMP.is_binary(vref)
        error("Cannot set the variable_ref $(vref) to integer as it " *
              "is already binary.")
    end
    gvref = _make_variable_ref(JuMP.owner_model(vref), JuMP.index(vref))
    cref = JuMP.add_constraint(JuMP.owner_model(vref),
                               JuMP.ScalarConstraint(gvref, MOI.Integer()),
                               is_info_constr = true)
    _set_integer_index(vref, JuMP.index(cref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, true))
    return
end

"""
    JuMP.IntegerRef(vref::UserDecisionVariableRef)::InfOptConstraintRef

Extend `JuMP.IntegerRef` to return a constraint reference to the constraint 
constrainting `vref` to be integer. Errors if one does not exist.

**Example**
```julia-repl
julia> cref = IntegerRef(vref)
var integer
```
"""
function JuMP.IntegerRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = _integer_index(vref)
    model = JuMP.owner_model(vref)
    return _make_constraint_ref(model, cindex)
end

"""
    JuMP.unset_integer(vref::UserDecisionVariableRef)::Nothing

Extend `JuMP.unset_integer` to unset `vref` as an integer variable. Errors if it 
is not an integer variable.

```julia-repl
julia> unset_integer(vref)

julia> is_integer(vref)
false
```
"""
function JuMP.unset_integer(vref::UserDecisionVariableRef)::Nothing
    JuMP.delete(JuMP.owner_model(vref), JuMP.IntegerRef(vref))
    _set_integer_index(vref, nothing)
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, false))
    return
end

################################################################################
#                          MODEL VARIABLE QUERIES
################################################################################
# Define finite variable struct for queries
struct FiniteVariable end

"""
    JuMP.num_variables(model::InfiniteModel, [type])::Int

Extend `JuMP.num_variables` to return the number of `InfiniteOpt` variables 
assigned to `model`. By default, the total number of infinite, semi-infinite, 
point, and finite variables is returned. The amount of a particular type is 
obtained by specifying the concrete variable type via `type`. Type options 
include:
 - `InfiniteVariable`: all infinite variables
 - `SemiInfiniteVariable`: all semi-infinite variables
 - `PointVariable`: all point variables
 - `FiniteVariable`: all finite variables

**Example**
```julia-repl
julia> num_variables(model)
3

julia> num_variables(model, InfiniteVariable)
2
```
"""
function JuMP.num_variables(model::InfiniteModel)::Int
    num_vars = JuMP.num_variables(model, InfiniteVariable)
    num_vars += JuMP.num_variables(model, SemiInfiniteVariable)
    num_vars += JuMP.num_variables(model, PointVariable)
    num_vars += JuMP.num_variables(model, FiniteVariable)
    return num_vars
end

# Particular variable types
function JuMP.num_variables(model::InfiniteModel,
                            type::Type{C})::Int where {C}
    return length(_data_dictionary(model, type))
end

"""
    JuMP.all_variables(model::InfiniteModel, [type])::Vector{GeneralVariableRef}

Extend `JuMP.all_variables`] to return a list of all the variable references 
associated with `model`. By default, all of the infinite, semi-infinite, point, 
and finite variables is returned. Those of a particular type is obtained by 
specifying the concrete variable type via `type`. Type options include:
 - `InfiniteVariable`: all infinite variables
 - `SemiInfiniteVariable`: all semi-infinite variables
 - `PointVariable`: all point variables
 - `FiniteVariable`: all finite variables

**Examples**
```julia-repl
julia> all_variables(model)
4-element Array{GeneralVariableRef,1}:
 y(t)
 w(t, x)
 y(0)
 z

julia> all_variables(model, PointVariable)
1-element Array{GeneralVariableRef,1}:
 y(0)
```
"""
function JuMP.all_variables(model::InfiniteModel)::Vector{GeneralVariableRef}
    vrefs_list = JuMP.all_variables(model, InfiniteVariable)
    append!(vrefs_list, JuMP.all_variables(model, SemiInfiniteVariable))
    append!(vrefs_list, JuMP.all_variables(model, PointVariable))
    append!(vrefs_list, JuMP.all_variables(model, FiniteVariable))
    return vrefs_list
end

# Particular variable types
function JuMP.all_variables(model::InfiniteModel,
                            type::Type{C}
                            )::Vector{GeneralVariableRef} where {C}
    vrefs_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model, type))
    for (i, (index, _)) in enumerate(_data_dictionary(model, type))
        vrefs_list[i] = _make_variable_ref(model, index)
    end
    return vrefs_list
end

################################################################################
#                                 DELETION
################################################################################
# Helper function to delete the info constraints
function _delete_info_constraints(vref::UserDecisionVariableRef)::Nothing
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
    return
end

"""
    JuMP.delete(model::InfiniteModel, vref::DecisionVariableRef)::Nothing

Extend `JuMP.delete` to delete `InfiniteOpt` variables and their dependencies. 
Errors if variable is invalid, meaning it has already been deleted or it belongs 
to another model.

**Example**
```julia-repl
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z ≥ 0.0
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]
 g(0.5) = 0

julia> delete(model, g)

julia> print(model)
Min measure(t) + z
Subject to
 z ≥ 0.0
 z ≥ 42.0
```
"""
function JuMP.delete(model::InfiniteModel, vref::DecisionVariableRef)::Nothing
    @assert JuMP.is_valid(model, vref) "Variable is invalid."
    # update the optimizer model status
    if is_used(vref)
        set_optimizer_model_ready(model, false)
    end
    # delete attributes specific to the variable type
    _delete_variable_dependencies(vref)
    gvref = _make_variable_ref(model, JuMP.index(vref))
    # remove from measures if used
    for mindex in _measure_dependencies(vref)
        mref = dispatch_variable_ref(model, mindex)
        func = measure_function(mref)
        data = measure_data(mref)
        if func isa GeneralVariableRef
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_meas = Measure(new_func, data, Int[], Int[], true)
        else
            _remove_variable(func, gvref)
            new_meas = build_measure(func, data)
        end
        _set_core_variable_object(mref, new_meas)
    end
    # remove from constraints if used
    for cindex in copy(_constraint_dependencies(vref)) # copy in case a constraint is deleted
        cref = _make_constraint_ref(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_constraint_object(cref, new_constr)
            empty!(_object_numbers(cref))
        elseif func isa AbstractArray && any(isequal(gvref), func)
            JuMP.delete(model, cref)
        else
            _remove_variable(func, gvref)
            # update the object numbers if vref is infinite
            if vref isa Union{InfiniteVariableRef, SemiInfiniteVariableRef}
                _data_object(cref).object_nums = sort(_object_numbers(func))
            end
        end
    end
    # remove from objective if vref is in it
    if used_by_objective(vref)
        if JuMP.objective_function(model) isa GeneralVariableRef
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            JuMP.set_objective_function(model, new_func)
            JuMP.set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        else
            _remove_variable(JuMP.objective_function(model), gvref)
        end
    end
    # delete the variable information
    _delete_data_object(vref)
    return
end

################################################################################
#                           INTEGRALITY RELAXATION
################################################################################
"""
    JuMP.relax_integrality(model::InfiniteModel)::Function

Modifies `model` to "relax" all binary and integrality constraints on
variables. Specifically,

- Binary constraints are deleted, and variable bounds are tightened if
  necessary to ensure the variable is constrained to the interval ``[0, 1]``.
- Integrality constraints are deleted without modifying variable bounds.
- All other constraints are ignored (left in place). This includes discrete
  constraints like SOS and indicator constraints.

Returns a function that can be called without any arguments to restore the
original model. The behavior of this function is undefined if additional
changes are made to the affected variables in the meantime.

**Example**
```julia-repl
julia> undo_relax = relax_integrality(model);

julia> print(model)
Min x + ∫{t ∈ [0, 10]}(y(t))
Subject to
 x ≥ 0.0
 y(t) ≥ 1.0
 x ≤ 1.0
 y(t) ≤ 10.0

julia> undo_relax()

julia> print(model)
Min x + ∫{t ∈ [0, 10]}(y(t))
Subject to
 y(t) ≥ 1.0
 y(t) ≤ 10.0
 y(t) integer
 x binary
```
"""
function JuMP.relax_integrality(model::InfiniteModel)
    # TODO ensure variables are not semi-continous/integer
    # get the variable info 
    info_pre_relaxation = map(v -> (dispatch_variable_ref(v), 
                              _variable_info(dispatch_variable_ref(v))),
                              JuMP.all_variables(model))
    # relax the appropriate variables
    for (vref, info) in info_pre_relaxation
        if info.integer
            JuMP.unset_integer(vref)
        elseif info.binary
            JuMP.unset_binary(vref)
            if isnan(info.lower_bound)
                JuMP.set_lower_bound(vref, 0.0)
            else
                JuMP.set_lower_bound(vref, max(0.0, info.lower_bound))
            end 
            if isnan(info.upper_bound)
                JuMP.set_upper_bound(vref, 1.0)
            else 
                JuMP.set_upper_bound(vref, min(1.0, info.upper_bound))
            end
        end
    end
    function unrelax()::Nothing
        for (vref, info) in info_pre_relaxation
            if info.integer
                JuMP.set_integer(vref)
            elseif info.binary
                JuMP.set_binary(vref)
                if info.has_lb
                    JuMP.set_lower_bound(vref, info.lower_bound)
                else
                    JuMP.delete_lower_bound(vref)
                end
                if info.has_ub
                    JuMP.set_upper_bound(vref, info.upper_bound)
                else
                    JuMP.delete_upper_bound(vref)
                end
            end
        end
        return
    end
    return unrelax
end
