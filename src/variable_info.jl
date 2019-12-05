## Internal functions for getting and setting variable info
# Get info
function _variable_info(vref::InfOptVariableRef)::JuMP.VariableInfo
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].info
end

# Set info for infinite variables
function _update_variable_info(vref::InfiniteVariableRef,
                               info::JuMP.VariableInfo)
    parameter_refs = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_refs
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info,
                                                                 parameter_refs)
    return
end

# Set info for point variables
function _update_variable_info(vref::PointVariableRef, info::JuMP.VariableInfo)
    infinite_variable_ref = JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
    parameter_values = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = PointVariable(info,
                                                          infinite_variable_ref,
                                                          parameter_values)
    return
end

# Set info for hold variables
function _update_variable_info(vref::HoldVariableRef, info::JuMP.VariableInfo)
    parameter_bounds = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_bounds
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = HoldVariable(info,
                                                               parameter_bounds)
    return
end

"""
    JuMP.has_lower_bound(vref::InfOptVariableRef)::Bool

Extend [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::VariableRef)) to
return a `Bool` whether an `InfiniteOpt` variable has a lower bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> has_lower_bound(vref)
true
```
"""
JuMP.has_lower_bound(vref::InfOptVariableRef)::Bool = _variable_info(vref).has_lb

"""
    JuMP.lower_bound(vref::InfOptVariableRef)::Float64

Extend [`JuMP.lower_bound`](@ref JuMP.lower_bound(::VariableRef)) to return the
lower bound of an `InfiniteOpt` variable. Errors if `vref` doesn't have a lower
bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> lower_bound(vref)
0.0
```
"""
function JuMP.lower_bound(vref::InfOptVariableRef)::Float64
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return _variable_info(vref).lower_bound
end

# TODO modify to eliminate dictionaries if possible like SingleVariable in MOI
# Extend to return the index of the lower bound constraint associated with `vref`.
function JuMP._lower_bound_index(vref::InfOptVariableRef)::Int
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return JuMP.owner_model(vref).var_to_lower_bound[JuMP.index(vref)]
end

# Extend to specify the index `cindex` of a lower bound constraint for `vref`.
function _set_lower_bound_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_lower_bound[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_lower_bound(vref::InfOptVariableRef, lower::Number)

Extend [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::VariableRef, ::Number))
to specify the lower bound of an `InfiniteOpt` variable `vref`. Errors if `vref`
is fixed.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> set_lower_bound(vref, -1)

julia> lower_bound(vref)
-1.0
```
"""
function JuMP.set_lower_bound(vref::InfOptVariableRef, lower::Number)
    newset = MOI.GreaterThan(convert(Float64, lower))
    if JuMP.has_lower_bound(vref)
        cindex = JuMP._lower_bound_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref,
                                                                       newset)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    else
        @assert !JuMP.is_fixed(vref) "$vref is fixed, cannot set lower bound."
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, newset))
        _set_lower_bound_index(vref, JuMP.index(cref))
        JuMP.owner_model(vref).constr_in_var_info[JuMP.index(cref)] = true
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
    JuMP.LowerBoundRef(vref::InfOptVariableRef)::GeneralConstraintRef

Extend [`JuMP.LowerBoundRef`](@ref JuMP.LowerBoundRef(::VariableRef)) to extract
a constraint reference for the lower bound of `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> cref = LowerBoundRef(vref)
var ≥ 0.0
```
"""
function JuMP.LowerBoundRef(vref::InfOptVariableRef)::GeneralConstraintRef
    index = JuMP._lower_bound_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.delete_lower_bound(vref::InfOptVariableRef)

Extend [`JuMP.delete_lower_bound`](@ref JuMP.delete_lower_bound(::VariableRef))
to delete lower bound of `vref`. Errors if it doesn't have a lower bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> delete_lower_bound(vref)

julia> has_lower_bound(vref)
false
```
"""
function JuMP.delete_lower_bound(vref::InfOptVariableRef)
    JuMP.delete(JuMP.owner_model(vref), JuMP.LowerBoundRef(vref))
    delete!(JuMP.owner_model(vref).var_to_lower_bound, JuMP.index(vref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(false, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.has_upper_bound(vref::InfOptVariableRef)::Bool

Extend [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::VariableRef)) to
return a `Bool` whether an `InfiniteOpt` variable has an upper bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
julia> has_upper_bound(vref)
true
```
"""
JuMP.has_upper_bound(vref::InfOptVariableRef)::Bool = _variable_info(vref).has_ub

"""
    JuMP.upper_bound(vref::InfOptVariableRef)::Float64

Extend [`JuMP.upper_bound`](@ref JuMP.upper_bound(::VariableRef)) to return the
upper bound of an `InfiniteOpt` variable. Errors if `vref` doesn't have a upper bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
julia> upper_bound(vref)
0.0
```
"""
function JuMP.upper_bound(vref::InfOptVariableRef)::Float64
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return _variable_info(vref).upper_bound
end

# Extend to return the index of the upper bound constraint associated with `vref`.
function JuMP._upper_bound_index(vref::InfOptVariableRef)::Int
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return JuMP.owner_model(vref).var_to_upper_bound[JuMP.index(vref)]
end

# Extend to specify the index `cindex` of a upper bound constraint for `vref`
function _set_upper_bound_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_upper_bound[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_upper_bound(vref::InfOptVariableRef, upper::Number)

Extend [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::VariableRef, ::Number))
to specify the upper bound of an `InfiniteOpt` variable `vref`. Errors if `vref`
is fixed.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
julia> set_upper_bound(vref, 1)

julia> upper_bound(vref)
1.0
```
"""
function JuMP.set_upper_bound(vref::InfOptVariableRef, upper::Number)
    newset = MOI.LessThan(convert(Float64, upper))
    if JuMP.has_upper_bound(vref)
        cindex = JuMP._upper_bound_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref,
                                                                       newset)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
    else
        @assert !JuMP.is_fixed(vref)
        cref = JuMP.add_constraint(JuMP.owner_model(vref),
                                   JuMP.ScalarConstraint(vref, newset))
        _set_upper_bound_index(vref, JuMP.index(cref))
        JuMP.owner_model(vref).constr_in_var_info[JuMP.index(cref)] = true
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
    JuMP.UpperBoundRef(vref::InfOptVariableRef)::GeneralConstraintRef

Extend [`JuMP.UpperBoundRef`](@ref JuMP.UpperBoundRef(::VariableRef)) to extract
a constraint reference for the upper bound of `vref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
julia> cref = UpperBoundRef(vref)
var ≤ 1.0
```
"""
function JuMP.UpperBoundRef(vref::InfOptVariableRef)::GeneralConstraintRef
    index = JuMP._upper_bound_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.delete_upper_bound(vref::InfOptVariableRef)

Extend [`JuMP.delete_upper_bound`](@ref JuMP.delete_upper_bound(::VariableRef))
to delete the upper bound of `vref`. Errors if it doesn't have an upper bound.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 >= vref))
julia> delete_upper_bound(vref)

julia> has_upper_bound(vref)
false
```
"""
function JuMP.delete_upper_bound(vref::InfOptVariableRef)
    JuMP.delete(JuMP.owner_model(vref), JuMP.UpperBoundRef(vref))
    delete!(JuMP.owner_model(vref).var_to_upper_bound, JuMP.index(vref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  false, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.is_fixed(vref::InfOptVariableRef)::Bool

Extend [`JuMP.is_fixed`](@ref JuMP.is_fixed(::VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is fixed.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref == 0))
julia> is_fixed(vref)
true
```
"""
JuMP.is_fixed(vref::InfOptVariableRef)::Bool = _variable_info(vref).has_fix

"""
    JuMP.fix_value(vref::InfOptVariableRef)::Float64

Extend [`JuMP.fix_value`](@ref JuMP.fix_value(::VariableRef)) to return the fix
value of an `InfiniteOpt` variable. Errors if variable is not fixed.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref == 0))
julia> fix_value(vref)
0.0
```
"""
function JuMP.fix_value(vref::InfOptVariableRef)::Float64
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return _variable_info(vref).fixed_value
end

# Extend to return the index of the fix constraint associated with `vref`.
function JuMP._fix_index(vref::InfOptVariableRef)::Int
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return JuMP.owner_model(vref).var_to_fix[JuMP.index(vref)]
end

# Extend to set the index of the fix constraintassociated with `vref`.
function _set_fix_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_fix[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.fix(vref::InfOptVariableRef, value::Number; force::Bool = false)

Extend [`JuMP.fix`](@ref JuMP.fix(::VariableRef, ::Number)) to fix the value of
an `InfiniteOpt` variable. Errors if variable has a lower and/or an upper
bound(s) unless `force = true`.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref == 0); @hold_variable(m, vref2 <= 0))
julia> fix(vref, 3)

julia> fix_value(vref)
3.0

julia> fix(vref2, 2, force = true)

julia> fix_value(vref2)
2.0
```
"""
function JuMP.fix(vref::InfOptVariableRef, value::Number; force::Bool = false)
    new_set = MOI.EqualTo(convert(Float64, value))
    model = JuMP.owner_model(vref)
    if JuMP.is_fixed(vref)  # Update existing fixing constraint.
        cindex = JuMP._fix_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref,
                                                                       new_set)
        set_optimizer_model_ready(JuMP.owner_model(vref), false)
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
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(vref, new_set))
        _set_fix_index(vref, JuMP.index(cref))
        JuMP.owner_model(vref).constr_in_var_info[JuMP.index(cref)] = true
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
    JuMP.FixRef(vref::InfOptVariableRef)::GeneralConstraintRef

Extend [`JuMP.FixRef`](@ref JuMP.FixRef(::VariableRef)) to return the constraint
reference of the fix constraint associated with `vref`. Errors `vref` is not
fixed.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref == 0))
julia> cref = FixRef(vref)
var = 1.0
```
"""
function JuMP.FixRef(vref::InfOptVariableRef)::GeneralConstraintRef
    index = JuMP._fix_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unfix(vref::InfOptVariableRef)

Extend [`JuMP.unfix`](@ref JuMP.unfix(::VariableRef)) to unfix `vref`. Errors
if it is not fixed.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref == 0))
julia> unfix(vref)

julia> is_fixed(vref)
false
```
"""
function JuMP.unfix(vref::InfOptVariableRef)
    JuMP.delete(JuMP.owner_model(vref), JuMP.FixRef(vref))
    delete!(JuMP.owner_model(vref).var_to_fix, JuMP.index(vref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  false, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, info.integer))
    return
end

"""
    JuMP.start_value(vref::InfOptVariableRef)::Union{Nothing, Float64}

Extend [`JuMP.start_value`](@ref JuMP.start_value(::VariableRef)) to return
starting value of `InfiniteOpt` variable if it has one. Returns `nothing` otherwise.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, start = 0))
julia> start_value(vref)
0.0
```
"""
function JuMP.start_value(vref::InfOptVariableRef)::Union{Nothing, Float64}
    return _variable_info(vref).start
end

"""
    JuMP.set_start_value(vref::InfOptVariableRef, value::Number)

Extend [`JuMP.set_start_value`](@ref JuMP.set_start_value(::VariableRef, ::Number))
to specify the start value of `InfiniteOpt` variables.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, start = 0))
julia> set_start_value(vref, 1)

julia> start_value(vref)
1.0
```
"""
function JuMP.set_start_value(vref::InfOptVariableRef, value::Number)
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
    JuMP.is_binary(vref::InfOptVariableRef)::Bool

Extend [`JuMP.is_binary`](@ref JuMP.is_binary(::VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is binary.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Bin))
julia> is_binary(vref)
true
```
"""
JuMP.is_binary(vref::InfOptVariableRef)::Bool = _variable_info(vref).binary

# Extend to return the index of the binary constraint associated with `vref`.
function JuMP._binary_index(vref::InfOptVariableRef)::Int
    if !JuMP.is_binary(vref)
        error("Variable $(vref) is not binary.")
    end
    return JuMP.owner_model(vref).var_to_zero_one[JuMP.index(vref)]
end

# Extend to specify the index of the binary constraint associated with `vref`.
function _set_binary_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_zero_one[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_binary(vref::InfOptVariableRef)

Extend [`JuMP.set_binary`](@ref JuMP.set_binary(::VariableRef)) to specify an
`InfiniteOpt` variable as a binary variable. Errors if `vref` is an integer
variable.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Bin))
julia> set_binary(vref)

julia> is_binary(vref)
true
```
"""
function JuMP.set_binary(vref::InfOptVariableRef)
    if JuMP.is_binary(vref)
        return
    elseif JuMP.is_integer(vref)
        error("Cannot set the variable_ref $(vref) to binary as it " *
              "is already integer.")
    end
    cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, MOI.ZeroOne()))
    _set_binary_index(vref, JuMP.index(cref))
    JuMP.owner_model(vref).constr_in_var_info[JuMP.index(cref)] = true
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  true, info.integer))
    return
end

"""
    JuMP.BinaryRef(vref::InfOptVariableRef)::GeneralConstraintRef

Extend [`JuMP.BinaryRef`](@ref JuMP.BinaryRef(::VariableRef)) to return a
constraint reference to the constraint constrainting `vref` to be binary. Errors
if one does not exist.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Bin))
julia> cref = BinaryRef(vref)
var binary
```
"""
function JuMP.BinaryRef(vref::InfOptVariableRef)::GeneralConstraintRef
    index = JuMP._binary_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unset_binary(vref::InfOptVariableRef)

Extend [`JuMP.unset_binary`](@ref JuMP.unset_binary(::VariableRef)) to unset
`vref` as a binary variable. Errors if it is not binary.

```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Bin))
julia> unset_binary(vref)

julia> is_binary(vref)
false
```
"""
function JuMP.unset_binary(vref::InfOptVariableRef)
    JuMP.delete(JuMP.owner_model(vref), JuMP.BinaryRef(vref))
    delete!(JuMP.owner_model(vref).var_to_zero_one, JuMP.index(vref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  false, info.integer))
    return
end

"""
    JuMP.is_integer(vref::InfOptVariableRef)::Bool

Extend [`JuMP.is_integer`](@ref JuMP.is_integer(::VariableRef)) to return `Bool`
whether an `InfiniteOpt` variable is integer.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Int))
julia> is_integer(vref)
true
```
"""
JuMP.is_integer(vref::InfOptVariableRef)::Bool = _variable_info(vref).integer

# Extend to return the index of the integer constraintassociated with `vref`.
function JuMP._integer_index(vref::InfOptVariableRef)::Int
    if !JuMP.is_integer(vref)
        error("Variable $(vref) is not an integer.")
    end
    return JuMP.owner_model(vref).var_to_integrality[JuMP.index(vref)]
end

# Extend to specify the index of the integer constraint associated with `vref`.
function _set_integer_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_integrality[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_integer(vref::InfOptVariableRef)

Extend [`JuMP.set_integer`](@ref JuMP.set_integer(::VariableRef)) to specify an
`InfiniteOpt` variable as a integer variable. Errors if `vref` is an binary
variable.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Int))
julia> set_integery(vref)

julia> is_integer(vref)
true
```
"""
function JuMP.set_integer(vref::InfOptVariableRef)
    if JuMP.is_integer(vref)
        return
    elseif JuMP.is_binary(vref)
        error("Cannot set the variable_ref $(vref) to integer as it " *
              "is already binary.")
    end
    cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref,
                               MOI.Integer()))
    _set_integer_index(vref, JuMP.index(cref))
    JuMP.owner_model(vref).constr_in_var_info[JuMP.index(cref)] = true
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, true))
    return
end

"""
    JuMP.IntegerRef(vref::InfOptVariableRef)::GeneralConstraintRef

Extend [`JuMP.IntegerRef`](@ref JuMP.IntegerRef(::VariableRef)) to return a
constraint reference to the constraint constrainting `vref` to be integer.
Errors if one does not exist.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Int))
julia> cref = IntegerRef(vref)
var integer
```
"""
function JuMP.IntegerRef(vref::InfOptVariableRef)::GeneralConstraintRef
    index = JuMP._integer_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unset_integer(vref::InfOptVariableRef)

Extend [`JuMP.unset_integer`](@ref JuMP.unset_integer(::VariableRef)) to unset
`vref` as an integer variable. Errors if it is not an integer variable.

```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, Int))
julia> unset_integer(vref)

julia> is_integer(vref)
false
```
"""
function JuMP.unset_integer(vref::InfOptVariableRef)
    JuMP.delete(JuMP.owner_model(vref), JuMP.IntegerRef(vref))
    delete!(JuMP.owner_model(vref).var_to_integrality, JuMP.index(vref))
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                                  info.has_ub, info.upper_bound,
                                                  info.has_fix, info.fixed_value,
                                                  info.has_start, info.start,
                                                  info.binary, false))
    return
end
