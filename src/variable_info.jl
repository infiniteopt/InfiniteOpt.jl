# Internal functions
_variable_info(vref::InfOptVariableRef) = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
function _update_variable_info(vref::InfiniteVariableRef, info::JuMP.VariableInfo)
    parameter_refs = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_refs
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, parameter_refs)
    return
end
function _update_variable_info(vref::PointVariableRef, info::JuMP.VariableInfo)
    infinite_variable_ref = JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
    parameter_values = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = PointVariable(info, infinite_variable_ref, parameter_values)
    return
end
function _update_variable_info(vref::GlobalVariableRef, info::JuMP.VariableInfo)
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = GlobalVariable(info)
    return
end

"""
    JuMP.has_lower_bound(vref::InfOptVariableRef)
Extend the `JuMP.has_lower_bound` function to accomodate our new variable types.
"""
JuMP.has_lower_bound(vref::InfOptVariableRef) = _variable_info(vref).has_lb

"""
    JuMP.lower_bound(vref::InfOptVariableRef)::Float64
Extend the `JuMP.lower_bound` function to accomodate our new variable types.
"""
function JuMP.lower_bound(vref::InfOptVariableRef)::Float64
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return _variable_info(vref).lower_bound
end

"""
    JuMP.lower_bound_index(vref::InfOptVariableRef)
Extend the `JuMP.lower_bound_index` function to accomodate our new variable types.
"""
function JuMP.lower_bound_index(vref::InfOptVariableRef)
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return JuMP.owner_model(vref).var_to_lower_bound[JuMP.index(vref)]
end

"""
    JuMP.set_lower_bound_index(vref::InfOptVariableRef, cindex::Int)
Extend the `JuMP.set_lower_bound_index` function to accomodate our new variable types.
"""
function JuMP.set_lower_bound_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_lower_bound[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_lower_bound(vref::InfOptVariableRef, lower::Number)
Extend the `JuMP.set_lower_bound` function to accomodate our new variable types.
"""
function JuMP.set_lower_bound(vref::InfOptVariableRef, lower::Number)
    newset = MOI.GreaterThan(convert(Float64, lower))
    if JuMP.has_lower_bound(vref)
        cindex = JuMP.lower_bound_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref, newset)
        set_optimizer_model_status(JuMP.owner_model(vref), false)
    else
        @assert !JuMP.is_fixed(vref)
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_lower_bound_index(vref, JuMP.index(cref))
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
    JuMP.LowerBoundRef(vref::InfOptVariableRef)
Extract a constraint reference for the lower bound.
"""
function JuMP.LowerBoundRef(vref::InfOptVariableRef)
    index = JuMP.lower_bound_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.delete_lower_bound(vref::InfOptVariableRef)
Extend the `JuMP.delete_lower_bound` function to accomodate our new variable types.
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
    JuMP.has_upper_bound(vref::InfOptVariableRef
Extend the `JuMP.has_upper_bound` function to accomodate our new variable types.
"""
JuMP.has_upper_bound(vref::InfOptVariableRef) = _variable_info(vref).has_ub

"""
    JuMP.upper_bound(vref::InfOptVariableRef)
Extend the `JuMP.upper_bound` function to accomodate our new variable types.
"""
function JuMP.upper_bound(vref::InfOptVariableRef)::Float64
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return _variable_info(vref).upper_bound
end

"""
    JuMP.upper_bound_index(vref::InfOptVariableRef)
Extend the `JuMP.upper_bound_index` function to accomodate our new variable types.
"""
function JuMP.upper_bound_index(vref::InfOptVariableRef)
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return JuMP.owner_model(vref).var_to_upper_bound[JuMP.index(vref)]
end

"""
    JuMP.set_upper_bound_index(vref::InfOptVariableRef, cindex::Int)
Extend the `JuMP.set_upper_bound_index` function to accomodate our new variable types.
"""
function JuMP.set_upper_bound_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_upper_bound[JuMP.index(vref)] = cindex
    return
end

"""
    JuMP.set_upper_bound(vref::InfOptVariableRef, upper::Number)
Extend the `JuMP.set_upper_bound` function to accomodate our new variable types.
"""
function JuMP.set_upper_bound(vref::InfOptVariableRef, upper::Number)
    newset = MOI.LessThan(convert(Float64, upper))
    if JuMP.has_upper_bound(vref)
        cindex = JuMP.upper_bound_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref, newset)
        set_optimizer_model_status(JuMP.owner_model(vref), false)
    else
        @assert !JuMP.is_fixed(vref)
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_upper_bound_index(vref, JuMP.index(cref))
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
    JuMP.UpperBoundRef(vref::InfOptVariableRef)
Extract a constraint reference for the upper bound.
"""
function JuMP.UpperBoundRef(vref::InfOptVariableRef)
    index = JuMP.upper_bound_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.delete_upper_bound(vref::InfOptVariableRef)
Extend the `JuMP.delete_upper_bound` function to accomodate our new variable types.
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
    JuMP.is_fixed(vref::InfOptVariableRef)
Extend the `is_fixed` function to accomodate our new variable types.
"""
JuMP.is_fixed(vref::InfOptVariableRef) = _variable_info(vref).has_fix

"""
    JuMP.fix_value(vref::InfOptVariableRef)
Extend the `JuMP.fix_value` function to accomodate our new variable types.
"""
function JuMP.fix_value(vref::InfOptVariableRef)::Float64
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return _variable_info(vref).fixed_value
end

"""
    JuMP.fix_index(vref::InfOptVariableRef)
Extend the `JuMP.fix_index` function to accomodate our new variable types.
"""
function JuMP.fix_index(vref::InfOptVariableRef)
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return JuMP.owner_model(vref).var_to_fix[JuMP.index(vref)]
end

"""
    JuMP.set_fix_index(vref::InfOptVariableRef, cindex::Int)
Extend the `JuMP.set_fix_index` function to accomodate our new variable types.
"""
function JuMP.set_fix_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_fix[JuMP.index(vref)] = cindex
end

"""
    JuMP.fix(vref::InfOptVariableRef, value::Number; force::Bool = false)
Extend the `JuMP.fix` function to accomodate our new variable types.
"""
function JuMP.fix(vref::InfOptVariableRef, value::Number; force::Bool = false)
    new_set = MOI.EqualTo(convert(Float64, value))
    model = JuMP.owner_model(vref)
    if JuMP.is_fixed(vref)  # Update existing fixing constraint.
        cindex = JuMP.fix_index(vref)
        JuMP.owner_model(vref).constrs[cindex] = JuMP.ScalarConstraint(vref, new_set)
        set_optimizer_model_status(JuMP.owner_model(vref), false)
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
        JuMP.set_fix_index(vref, JuMP.index(cref))
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
    JuMP.FixRef(vref::InfOptVariableRef)
Extend the `JuMP.FixRef` function for our new variable types.
"""
function JuMP.FixRef(vref::InfOptVariableRef)
    index = JuMP.fix_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unfix(vref::InfOptVariableRef)
Extend the `JuMP.unfix` function to accomodate our new variable types.
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
    JuMP.start_value(vref::InfOptVariableRef)
Extend the `JuMP.start_value` function to accomodate our new variable types.
"""
function JuMP.start_value(vref::InfOptVariableRef)::Union{Nothing, Float64}
    return _variable_info(vref).start
end

"""
    JuMP.set_start_value(vref::InfOptVariableRef, value::Number)
Extend the `JuMP.set_start_value` function to accomodate our new variable types.
"""
function JuMP.set_start_value(vref::InfOptVariableRef, value::Number)
    info = _variable_info(vref)
    set_optimizer_model_status(JuMP.owner_model(vref), false)
    _update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           true, Float64(value),
                                           info.binary, info.integer))
    return
end

"""
    JuMP.is_binary(vref::InfOptVariableRef)
Extend the `JuMP.is_binary` function to accomodate our new variable types.
"""
JuMP.is_binary(vref::InfOptVariableRef) = _variable_info(vref).binary

"""
    JuMP.binary_index(vref::InfOptVariableRef)
Extend the `JuMP.binary_index` function to accomodate our new variable types.
"""
function JuMP.binary_index(vref::InfOptVariableRef)
    if !JuMP.is_binary(vref)
        error("Variable $(vref) is not binary.")
    end
    return JuMP.owner_model(vref).var_to_zero_one[JuMP.index(vref)]
end

"""
    JuMP.set_binary_index(vref::InfOptVariableRef, cindex::Int)
Extend the `JuMP.set_binary_index` function to accomodate our new variable types.
"""
function JuMP.set_binary_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_zero_one[JuMP.index(vref)] = cindex
end

"""
    JuMP.set_binary(vref::InfOptVariableRef)
Extend the `JuMP.set_binary` function to accomodate our new variable types.
"""
function JuMP.set_binary(vref::InfOptVariableRef)
    if JuMP.is_binary(vref)
        return
    elseif JuMP.is_integer(vref)
        error("Cannot set the variable_ref $(vref) to binary as it " *
              "is already integer.")
    end
    cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, MOI.ZeroOne()))
    JuMP.set_binary_index(vref, JuMP.index(cref))
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
    JuMP.BinaryRef(vref::InfOptVariableRef)
Return a constraint reference to the constraint constrainting `vref` to be binary.
Errors if one does not exist.
"""
function JuMP.BinaryRef(vref::InfOptVariableRef)
    index = JuMP.binary_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unset_binary(vref::InfOptVariableRef)
Extend the `JuMP.unset_binary` function to accomodate our new variable types.
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
    JuMP.is_integer(vref::InfOptVariableRef)
Extend the `JuMP.is_integer` function to accomodate our new variable types.
"""
JuMP.is_integer(vref::InfOptVariableRef) = _variable_info(vref).integer

"""
    JuMP.integer_index(vref::InfOptVariableRef)
Extend the `JuMP.integer_index` function to accomodate our new variable types.
"""
function JuMP.integer_index(vref::InfOptVariableRef)
    if !JuMP.is_integer(vref)
        error("Variable $(vref) is not an integer.")
    end
    return JuMP.owner_model(vref).var_to_integrality[JuMP.index(vref)]
end

"""
    JuMP.set_integer_index(vref::InfOptVariableRef, cindex::Int)
Extend the `JuMP.set_integer_index` function to accomodate our new variable types.
"""
function JuMP.set_integer_index(vref::InfOptVariableRef, cindex::Int)
    JuMP.owner_model(vref).var_to_integrality[JuMP.index(vref)] = cindex
end

"""
    JuMP.set_integer(vref::InfOptVariableRef)
Extend the `JuMP.set_integer` function to accomodate our new variable types.
"""
function JuMP.set_integer(vref::InfOptVariableRef)
    if JuMP.is_integer(vref)
        return
    elseif JuMP.is_binary(vref)
        error("Cannot set the variable_ref $(vref) to integer as it " *
              "is already binary.")
    end
    cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, MOI.Integer()))
    JuMP.set_integer_index(vref, JuMP.index(cref))
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
    JuMP.IntegerRef(vref::InfOptVariableRef)
Return a constraint reference to the constraint constrainting `vref` to be integer.
Errors if one does not exist.
"""
function JuMP.IntegerRef(vref::InfOptVariableRef)
    index = JuMP.integer_index(vref)
    model = JuMP.owner_model(vref)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.unset_integer(vref::InfOptVariableRef)
Extend the `JuMP.unset_integer` function to accomodate our new variable types.
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
