################################################################################
#                               VARIABLE DEFINITION
################################################################################
# Define symbol inputs for different variable types
const Infinite = :Infinite
const Point = :Point
const Hold = :Hold

# Fallback _make_variable (the methods are defined in the variable files)
function _make_variable(_error::Function, info::JuMP.VariableInfo, type;
                        extra_kw_args...)
    _error("Unrecognized variable type `$(typeof(type).parameters[1])`, " *
           "should be Infinite, Point, or Hold.")
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                        var_type::Symbol;
                        [parameter_refs::Union{GeneralVariableRef,
                                              AbstractArray{<:GeneralVariableRef},
                                              Tuple, Nothing} = nothing,
                        infinite_variable_ref::Union{GeneralVariableRef,
                                                     Nothing} = nothing,
                        parameter_values::Union{Number, AbstractArray{<:Real},
                                                Tuple, Nothing} = nothing,
                        parameter_bounds::Union{ParameterBounds{GeneralVariableRef},
                                                Nothing} = nothing]
                        )::InfOptVariable

Extend the `JuMP.build_variable` function to accomodate `InfiniteOpt`
variable types. Returns the appropriate variable Datatype (i.e.,
[`InfiniteVariable`](@ref), [`PointVariable`](@ref), and
[`HoldVariable`](@ref)). Primarily, this method is to be used internally by the
appropriate constructor macros [`@infinite_variable`](@ref),
[`@point_variable`](@ref), and [`@hold_variable`](@ref). However, it can be
called manually to build `InfiniteOpt` variables. Errors if an unneeded keyword
argument is given or if the keywoard arguments are formatted incorrectly (e.g.,
`parameter_refs` contains repeated parameter references when an infinite variable
is defined). Also errors if needed keyword arguments are negated.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel())
julia> @independent_parameter(m, t in [0, 1])
t

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = t)
InfiniteVariable{GeneralVariableRef}(VariableInfo{Float64,Float64,Float64,Float64}(false, 0.0, false, 0.0, false, 0.0, false, 0.0, false, false), (t,), Int64[], Int64[])

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, info, Point, infinite_variable_ref = ivref,
                               parameter_values = 0.5)
PointVariable{GeneralVariableRef}(VariableInfo{Float64,Float64,Float64,Float64}(false, 0.0, false, 0.0, false, 0.0, true, 0.0, false, false), var_name(t), [0.5])

julia> hd_var = build_variable(error, info, Hold)
HoldVariable{GeneralVariableRef}(VariableInfo{Float64,Float64,Float64,Float64}(false, 0.0, false, 0.0, false, 0.0, false, 0.0, false, false), Subdomain bounds (0): )
```
"""
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                             var_type::Symbol;
                             macro_error::Union{Function, Nothing} = nothing,
                             kw_args...)::InfOptVariable
    if macro_error != nothing
        _error = macro_error # replace with macro error function
    end
    # make the variable and conduct necessary checks
    return _make_variable(_error, info, Val(var_type); kw_args...)
end

# Fallback
function _check_and_make_variable_ref(model::InfiniteModel, v::T) where {T}
    throw(ArgumentError("Invalid variable object type `$T`."))
end

"""
    JuMP.add_variable(model::InfiniteModel, var::InfOptVariable,
                      [name::String = ""])::GeneralVariableRef

Extend the [`JuMP.add_variable`](@ref JuMP.add_variable(::JuMP.Model, ::JuMP.ScalarVariable, ::String))
function to accomodate `InfiniteOpt` variable types. Adds a variable to an
infinite model `model` and returns a [`GeneralVariableRef`](@ref).
Primarily intended to be an internal function of the
constructor macros [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@hold_variable`](@ref). However, it can be used in combination with
[`JuMP.build_variable`](@ref) to add variables to an infinite model object.
Errors if invalid parameters reference(s) or an invalid infinite variable
reference is included in `var`.

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
                           name::String = "")::GeneralVariableRef
    dvref = _check_and_make_variable_ref(model, var)
    JuMP.set_name(dvref, name)
    vindex = JuMP.index(dvref)
    gvref = _make_variable_ref(model, vindex)
    if var.info.has_lb
        newset = MOI.GreaterThan(var.info.lower_bound)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_lower_bound_index(dvref, JuMP.index(cref))
    end
    if var.info.has_ub
        newset = MOI.LessThan(var.info.upper_bound)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_upper_bound_index(dvref, JuMP.index(cref))
    end
    if var.info.has_fix
        newset = MOI.EqualTo(var.info.fixed_value)
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_fix_index(dvref, JuMP.index(cref))
    end
    if var.info.binary
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.ZeroOne()),
                                   is_info_constr = true)
        _set_binary_index(dvref, JuMP.index(cref))
    elseif var.info.integer
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.Integer()),
                                   is_info_constr = true)
        _set_integer_index(dvref, JuMP.index(cref))
    end
    return gvref
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _measure_dependencies
function _measure_dependencies(vref::DecisionVariableRef)::Vector{MeasureIndex}
    return _data_object(vref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(vref::DecisionVariableRef)::Vector{ConstraintIndex}
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

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to return the names of
`InfiniteOpt` variables.

**Example**
```julia-repl
julia> name(vref)
"var_name"
```
"""
function JuMP.name(vref::DecisionVariableRef)::String
    return _data_object(vref).name
end

"""
    JuMP.set_name(vref::DecisionVariableRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of decision variables.

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
function _make_variable_ref(model::InfiniteModel, index::ObjectIndex)::GeneralVariableRef
    return GeneralVariableRef(model, index.value, typeof(index))
end

# Get the name_to_var Dictionary
function _var_name_dict(model::InfiniteModel)::Union{Nothing, Dict{String, ObjectIndex}}
    return model.name_to_var
end

# Update name_to_var
function _update_var_name_dict(model::InfiniteModel, var_dict::MOIUC.CleverDict)::Nothing
    name_dict = _var_name_dict(model)
    for (index, data_object) in var_dict
        var_name = data_object.name
        if haskey(name_dict, var_name)
            name_dict[var_name] = HoldVariableIndex(-1) # dumby value
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

Extend [`JuMP.variable_by_name`](@ref JuMP.variable_by_name(::JuMP.Model, ::String))
for `InfiniteModel` objects. Return the variable reference assoociated with a
variable name. Errors if multiple variables have the same name. Returns nothing
if no such name exists.

**Examples**
```julia-repl
julia> variable_by_name(m, "var_name")
var_name

julia> variable_by_name(m, "fake_name")

```
"""
function JuMP.variable_by_name(model::InfiniteModel,
                               name::String)::Union{GeneralVariableRef, Nothing}
    if _var_name_dict(model) === nothing
        model.name_to_var = Dict{String, ObjectIndex}()
        _update_var_name_dict(model, model.infinite_vars)
        _update_var_name_dict(model, model.reduced_vars)
        _update_var_name_dict(model, model.point_vars)
        _update_var_name_dict(model, model.hold_vars)
    end
    index = get(_var_name_dict(model), name, nothing)
    if index isa Nothing
        return nothing
    elseif index == HoldVariableIndex(-1)
        error("Multiple variables have the name $name.")
    else
        return _make_variable_ref(model, index)
    end
end

################################################################################
#                         VARIABLE OBJECT MODIFICATION
################################################################################
# Extend _set_core_variable_object
function _set_core_variable_object(vref::DecisionVariableRef,
                                   object::InfOptVariable)::Nothing
    _data_object(vref).variable = object
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Get info
function _variable_info(vref::UserDecisionVariableRef
    )::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    return _core_variable_object(vref).info
end

# Make a temporary constraint reference to use the data accessors
function _temp_constraint_ref(model::InfiniteModel,
                              cindex::ConstraintIndex)::FiniteConstraintRef
    return FiniteConstraintRef(model, cindex, JuMP.ScalarShape())
end

"""
    JuMP.has_lower_bound(vref::UserDecisionVariableRef)::Bool

Extend [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::JuMP.VariableRef)) to
return a `Bool` whether an `InfiniteOpt` variable has a lower bound.

**Example**
```julia-repl
julia> has_lower_bound(vref)
true
```
"""
JuMP.has_lower_bound(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).has_lb

"""
    JuMP.lower_bound(vref::UserDecisionVariableRef)::Float64

Extend [`JuMP.lower_bound`](@ref JuMP.lower_bound(::JuMP.VariableRef)) to return the
lower bound of an `InfiniteOpt` variable. Errors if `vref` doesn't have a lower
bound.

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
function JuMP._lower_bound_index(vref::UserDecisionVariableRef)::ConstraintIndex
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return _data_object(vref).lower_bound_index
end

# Extend to specify the index `cindex` of a lower bound constraint for `vref`.
function _set_lower_bound_index(vref::UserDecisionVariableRef,
                                cindex::Union{ConstraintIndex, Nothing})::Nothing
    _data_object(vref).lower_bound_index = cindex
    return
end

"""
    JuMP.set_lower_bound(vref::UserDecisionVariableRef, lower::Real)::Nothing

Extend [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::JuMP.VariableRef, ::Number))
to specify the lower bound of an `InfiniteOpt` variable `vref`. Errors if `vref`
is fixed.

**Example**
```julia-repl
julia> set_lower_bound(vref, -1)

julia> lower_bound(vref)
-1.0
```
"""
function JuMP.set_lower_bound(vref::UserDecisionVariableRef,
                              lower::Real)::Nothing
    newset = MOI.GreaterThan(convert(Float64, lower))
    new_constr = JuMP.ScalarConstraint(vref, newset)
    model = JuMP.owner_model(vref)
    if JuMP.has_lower_bound(vref)
        cindex = JuMP._lower_bound_index(vref)
        cref = _temp_constraint_ref(model, cindex)
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

Extend [`JuMP.LowerBoundRef`](@ref JuMP.LowerBoundRef(::JuMP.VariableRef)) to extract
a constraint reference for the lower bound of `vref`.

**Example**
```julia-repl
var ≥ 0.0
```
"""
function JuMP.LowerBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = JuMP._lower_bound_index(vref)
    model = JuMP.owner_model(vref)
    cref = _temp_constraint_ref(model, cindex)
    if isempty(_object_numbers(cref))
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(_core_constraint_object(cref)))
    else
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(_core_constraint_object(cref)))
    end
end

"""
    JuMP.delete_lower_bound(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.delete_lower_bound`](@ref JuMP.delete_lower_bound(::JuMP.VariableRef))
to delete lower bound of `vref`. Errors if it doesn't have a lower bound.

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

Extend [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::JuMP.VariableRef)) to
return a `Bool` whether an `InfiniteOpt` variable has an upper bound.

**Example**
```julia-repl
julia> has_upper_bound(vref)
true
```
"""
JuMP.has_upper_bound(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).has_ub

"""
    JuMP.upper_bound(vref::UserDecisionVariableRef)::Float64

Extend [`JuMP.upper_bound`](@ref JuMP.upper_bound(::JuMP.VariableRef)) to return the
upper bound of an `InfiniteOpt` variable. Errors if `vref` doesn't have a upper bound.

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
function JuMP._upper_bound_index(vref::UserDecisionVariableRef)::ConstraintIndex
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return _data_object(vref).upper_bound_index
end

# Extend to specify the index `cindex` of a upper bound constraint for `vref`
function _set_upper_bound_index(vref::UserDecisionVariableRef,
                                cindex::Union{ConstraintIndex, Nothing})::Nothing
    _data_object(vref).upper_bound_index = cindex
    return
end

"""
    JuMP.set_upper_bound(vref::UserDecisionVariableRef, upper::Real)::Nothing

Extend [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::JuMP.VariableRef, ::Number))
to specify the upper bound of an `InfiniteOpt` variable `vref`. Errors if `vref`
is fixed.

**Example**
```julia-repl
julia> set_upper_bound(vref, 1)

julia> upper_bound(vref)
1.0
```
"""
function JuMP.set_upper_bound(vref::UserDecisionVariableRef,
                              upper::Real)::Nothing
    newset = MOI.LessThan(convert(Float64, upper))
    new_constr = JuMP.ScalarConstraint(vref, newset)
    model = JuMP.owner_model(vref)
    if JuMP.has_upper_bound(vref)
        cindex = JuMP._upper_bound_index(vref)
        cref = _temp_constraint_ref(model, cindex)
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

Extend [`JuMP.UpperBoundRef`](@ref JuMP.UpperBoundRef(::JuMP.VariableRef)) to extract
a constraint reference for the upper bound of `vref`.

**Example**
```julia-repl
julia> cref = UpperBoundRef(vref)
var ≤ 1.0
```
"""
function JuMP.UpperBoundRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = JuMP._upper_bound_index(vref)
    model = JuMP.owner_model(vref)
    cref = _temp_constraint_ref(model, cindex)
    if isempty(_object_numbers(cref))
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(_core_constraint_object(cref)))
    else
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(_core_constraint_object(cref)))
    end
end

"""
    JuMP.delete_upper_bound(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.delete_upper_bound`](@ref JuMP.delete_upper_bound(::JuMP.VariableRef))
to delete the upper bound of `vref`. Errors if it doesn't have an upper bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
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

Extend [`JuMP.is_fixed`](@ref JuMP.is_fixed(::JuMP.VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is fixed.

**Example**
```julia-repl
julia> is_fixed(vref)
true
```
"""
JuMP.is_fixed(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).has_fix

"""
    JuMP.fix_value(vref::UserDecisionVariableRef)::Float64

Extend [`JuMP.fix_value`](@ref JuMP.fix_value(::JuMP.VariableRef)) to return the fix
value of an `InfiniteOpt` variable. Errors if variable is not fixed.

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
function JuMP._fix_index(vref::UserDecisionVariableRef)::ConstraintIndex
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return _data_object(vref).fix_index
end

# Extend to set the index of the fix constraintassociated with `vref`.
function _set_fix_index(vref::UserDecisionVariableRef,
                        cindex::Union{ConstraintIndex, Nothing})::Nothing
    _data_object(vref).fix_index = cindex
    return
end

"""
    JuMP.fix(vref::UserDecisionVariableRef, value::Real;
             force::Bool = false)::Nothing

Extend [`JuMP.fix`](@ref JuMP.fix(::JuMP.VariableRef, ::Number)) to fix the value of
an `InfiniteOpt` variable. Errors if variable has a lower and/or an upper
bound(s) unless `force = true`.

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
function JuMP.fix(vref::UserDecisionVariableRef, value::Real;
                  force::Bool = false)::Nothing
    new_set = MOI.EqualTo(convert(Float64, value))
    model = JuMP.owner_model(vref)
    new_constr = JuMP.ScalarConstraint(vref, new_set)
    if JuMP.is_fixed(vref)  # Update existing fixing constraint.
        cindex = JuMP._fix_index(vref)
        cref = _temp_constraint_ref(model, cindex)
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

Extend [`JuMP.FixRef`](@ref JuMP.FixRef(::JuMP.VariableRef)) to return the constraint
reference of the fix constraint associated with `vref`. Errors `vref` is not
fixed.

**Examples**
```julia-repl
julia> cref = FixRef(vref)
var = 1.0
```
"""
function JuMP.FixRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = JuMP._fix_index(vref)
    model = JuMP.owner_model(vref)
    cref = _temp_constraint_ref(model, cindex)
    if isempty(_object_numbers(cref))
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(_core_constraint_object(cref)))
    else
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(_core_constraint_object(cref)))
    end
end

"""
    JuMP.unfix(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.unfix`](@ref JuMP.unfix(::JuMP.VariableRef)) to unfix `vref`. Errors
if it is not fixed.

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

Extend [`JuMP.start_value`](@ref JuMP.start_value(::JuMP.VariableRef)) to return
starting value of `InfiniteOpt` variable if it has one. Returns `nothing` otherwise.

**Example**
```julia-repl
julia> start_value(vref)
0.0
```
"""
function JuMP.start_value(vref::UserDecisionVariableRef)::Union{Nothing, Float64}
    return _variable_info(vref).start
end

"""
    JuMP.set_start_value(vref::UserDecisionVariableRef, value::Real)::Nothing

Extend [`JuMP.set_start_value`](@ref JuMP.set_start_value(::JuMP.VariableRef, ::Number))
to specify the start value of `InfiniteOpt` variables.

**Example**
```julia-repl
julia> set_start_value(vref, 1)

julia> start_value(vref)
1.0
```
"""
function JuMP.set_start_value(vref::UserDecisionVariableRef,
                              value::Real)::Nothing
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

Extend [`JuMP.is_binary`](@ref JuMP.is_binary(::JuMP.VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is binary.

**Example**
```julia-repl
true
```
"""
JuMP.is_binary(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).binary

# Extend to return the index of the binary constraint associated with `vref`.
function JuMP._binary_index(vref::UserDecisionVariableRef)::ConstraintIndex
    if !JuMP.is_binary(vref)
        error("Variable $(vref) is not binary.")
    end
    return _data_object(vref).zero_one_index
end

# Extend to specify the index of the binary constraint associated with `vref`.
function _set_binary_index(vref::UserDecisionVariableRef,
                           cindex::Union{ConstraintIndex, Nothing})::Nothing
    _data_object(vref).zero_one_index = cindex
    return
end

"""
    JuMP.set_binary(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.set_binary`](@ref JuMP.set_binary(::JuMP.VariableRef)) to specify an
`InfiniteOpt` variable as a binary variable. Errors if `vref` is an integer
variable.

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
    cref = JuMP.add_constraint(JuMP.owner_model(vref),
                               JuMP.ScalarConstraint(vref, MOI.ZeroOne()),
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

Extend [`JuMP.BinaryRef`](@ref JuMP.BinaryRef(::JuMP.VariableRef)) to return a
constraint reference to the constraint constrainting `vref` to be binary. Errors
if one does not exist.

**Example**
```julia-repl
julia> cref = BinaryRef(vref)
var binary
```
"""
function JuMP.BinaryRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = JuMP._binary_index(vref)
    model = JuMP.owner_model(vref)
    cref = _temp_constraint_ref(model, cindex)
    if isempty(_object_numbers(cref))
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(_core_constraint_object(cref)))
    else
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(_core_constraint_object(cref)))
    end
end

"""
    JuMP.unset_binary(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.unset_binary`](@ref JuMP.unset_binary(::JuMP.VariableRef)) to unset
`vref` as a binary variable. Errors if it is not binary.

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

Extend [`JuMP.is_integer`](@ref JuMP.is_integer(::JuMP.VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is integer.

**Example**
```julia-repl
julia> is_integer(vref)
true
```
"""
JuMP.is_integer(vref::UserDecisionVariableRef)::Bool = _variable_info(vref).integer

# Extend to return the index of the integer constraintassociated with `vref`.
function JuMP._integer_index(vref::UserDecisionVariableRef)::ConstraintIndex
    if !JuMP.is_integer(vref)
        error("Variable $(vref) is not an integer.")
    end
    return _data_object(vref).integrality_index
end

# Extend to specify the index of the integer constraint associated with `vref`.
function _set_integer_index(vref::UserDecisionVariableRef,
                            cindex::Union{ConstraintIndex, Nothing})::Nothing
    _data_object(vref).integrality_index = cindex
    return
end

"""
    JuMP.set_integer(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.set_integer`](@ref JuMP.set_integer(::JuMP.VariableRef)) to specify an
`InfiniteOpt` variable as a integer variable. Errors if `vref` is an binary
variable.

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
    cref = JuMP.add_constraint(JuMP.owner_model(vref),
                               JuMP.ScalarConstraint(vref, MOI.Integer()),
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

Extend [`JuMP.IntegerRef`](@ref JuMP.IntegerRef(::JuMP.VariableRef)) to return a
constraint reference to the constraint constrainting `vref` to be integer.
Errors if one does not exist.

**Example**
```julia-repl
julia> cref = IntegerRef(vref)
var integer
```
"""
function JuMP.IntegerRef(vref::UserDecisionVariableRef)::InfOptConstraintRef
    cindex = JuMP._integer_index(vref)
    model = JuMP.owner_model(vref)
    cref = _temp_constraint_ref(model, cindex)
    if isempty(_object_numbers(cref))
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(_core_constraint_object(cref)))
    else
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(_core_constraint_object(cref)))
    end
end

"""
    JuMP.unset_integer(vref::UserDecisionVariableRef)::Nothing

Extend [`JuMP.unset_integer`](@ref JuMP.unset_integer(::JuMP.VariableRef)) to unset
`vref` as an integer variable. Errors if it is not an integer variable.

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
"""
    JuMP.num_variables(model::InfiniteModel,
                       [type::Type{InfOptVariable} = InfOptVariable])::Int

Extend [`JuMP.num_variables`](@ref JuMP.num_variables(::JuMP.Model)) to return the
number of `InfiniteOpt` variables assigned to `model`. By default, the total
number of infinite, reduced, point, and hold variables is returned. The amount
of a particular type is obtained by specifying the concrete variable type
of [`InfOptVariable`](@ref) via `type`. Type options include:
 - `InfOptVariable`: all variables
 - `InfiniteVariable`: all infinite variables
 - `ReducedInfiniteVariable`: all reduced infinite variables
 - `PointVariable`: all point variables
 - `HoldVariable`: all hold variables

**Example**
```julia-repl
julia> num_variables(model)
3

julia> num_variables(model, InfiniteVariable)
2
```
"""
function JuMP.num_variables(model::InfiniteModel,
                            type::Type{InfOptVariable} = InfOptVariable
                            )::Int
    num_vars = JuMP.num_variables(model, InfiniteVariable)
    num_vars += JuMP.num_variables(model, ReducedInfiniteVariable)
    num_vars += JuMP.num_variables(model, PointVariable)
    num_vars += JuMP.num_variables(model, HoldVariable)
    return num_vars
end

# Particular variable types
function JuMP.num_variables(model::InfiniteModel,
                            type::Type{C})::Int where {C <: InfOptVariable}
    return length(_data_dictionary(model, type))
end

"""
    JuMP.all_variables(model::InfiniteModel,
                       type::Type{InfOptVariable} = InfOptVariable
                       )::Vector{GeneralVariableRef}

Extend [`JuMP.all_variables`](@ref JuMP.all_variables(::JuMP.Model)) to return a
list of all the variable references associated with `model`. By default, all
of the infinite, reduced, point, and hold variables is returned. Those
of a particular type is obtained by specifying the concrete variable type
of [`InfOptVariable`](@ref) via `type`. Type options include:
 - `InfOptVariable`: all variables
 - `InfiniteVariable`: all infinite variables
 - `ReducedInfiniteVariable`: all reduced infinite variables
 - `PointVariable`: all point variables
 - `HoldVariable`: all hold variables

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
function JuMP.all_variables(model::InfiniteModel,
                            type::Type{InfOptVariable} = InfOptVariable
                            )::Vector{GeneralVariableRef}
    vrefs_list = JuMP.all_variables(model, InfiniteVariable)
    append!(vrefs_list, JuMP.all_variables(model, ReducedInfiniteVariable))
    append!(vrefs_list, JuMP.all_variables(model, PointVariable))
    append!(vrefs_list, JuMP.all_variables(model, HoldVariable))
    return vrefs_list
end

# Particular variable types
function JuMP.all_variables(model::InfiniteModel,
                            type::Type{C}
                            )::Vector{GeneralVariableRef} where {C <: InfOptVariable}
    vrefs_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model, type))
    counter = 1
    for (index, object) in _data_dictionary(model, type)
        vrefs_list[counter] = _make_variable_ref(model, index)
        counter += 1
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

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)) to delete
`InfiniteOpt` variables and their dependencies. Errors if variable is invalid,
meaning it has already been deleted or it belongs to another model.

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
    # remove from measures if used
    for mindex in _measure_dependencies(vref)
        mref = dispatch_variable_ref(model, mindex)
        func = measure_function(mref)
        if func isa GeneralVariableRef
            data = measure_data(mref)
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_meas = Measure(new_func, data, [], [])
            _set_core_variable_object(mref, new_meas)
        else
            _remove_variable(func, vref)
            # TODO update object/param numbers via measure specifc functions
        end
        meas = _core_variable_object(mref)
        JuMP.set_name(mref, _make_meas_name(meas))
    end
    # remove from constraints if used
    for cindex in _constraint_dependencies(vref)
        cref = _temp_constraint_ref(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.AffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_constraint_object(cref, new_constr)
            empty!(_object_numbers(cref))
        else
            _remove_variable(func, vref)
            # update the object numbers if vref is infinite
            if vref isa Union{InfiniteVariableRef, ReducedInfiniteVariableRef}
                new_obj_nums = _object_numbers(func)
                filter!(e -> e in new_obj_nums, _object_numbers(cref))
            end
        end
    end
    # remove from objective if vref is in it
    if used_by_objective(vref)
        if JuMP.objective_function(model) isa GeneralVariableRef
            new_func = zero(JuMP.AffExpr{Float64, GeneralVariableRef})
            JuMP.set_objective_function(model, new_func)
        else
            _remove_variable(JuMP.objective_function(model), vref)
        end
    end
    # delete the variable information
    _delete_data_object(vref)
    return
end
