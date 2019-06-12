# Define symbol inputs for different variable types
const Infinite = :Infinite
const Point = :Point
const Global = :Global

# Extend Base.copy for new variable types
Base.copy(v::GeneralVariableRef) = v
Base.copy(v::InfiniteVariableRef, new_model::InfiniteModel) = InfiniteVariableRef(new_model, v.index)
Base.copy(v::GlobalVariableRef, new_model::InfiniteModel) = GlobalVariableRef(new_model, v.index)
Base.copy(v::PointVariableRef, new_model::InfiniteModel) = PointVariableRef(new_model, v.index)

# Extend other Base functions
Base.:(==)(v::GeneralVariableRef, w::GeneralVariableRef) = v.model === w.model && v.index == w.index
Base.broadcastable(v::GeneralVariableRef) = Ref(v)

# Extend JuMP functions
JuMP.isequal_canonical(v::GeneralVariableRef, w::GeneralVariableRef) = v == w
JuMP.variable_type(model::InfiniteModel) = GeneralVariableRef
function JuMP.variable_type(model::InfiniteModel, type::Symbol)
    if type == Infinite
        return InfiniteVariableRef
    elseif type == Point
        return PointVariableRef
    elseif type == Global
        return GlobalVariableRef
    elseif type == Parameter
        return ParameterRef
    else
        error("Invalid variable type.")
    end
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, test_arg::Int;
                        param_set::Union{AbstractInfiniteSet, Vector{AbstractInfiniteSet}} = EmptySet(),
                        extra_kw_args...)
Extend the `JuMP.build_variable` function to accomodate our new variable types.
"""
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, var_type::Symbol;
                             param_refs::Union{ParameterRef, Tuple, Nothing} = nothing,
                             param_set::Union{AbstractInfiniteSet, Nothing} = nothing,
                             inf_var_ref::Union{InfiniteVariableRef, Nothing} = nothing,
                             param_values::Union{Tuple, Nothing} = nothing,
                             extra_kw_args...)
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    if !(var_type in [Infinite, Point, Global, Parameter])
        _error("Unrecognized variable type $var_type, should be Inf, Point, Global, or Parameter.")
    end
    if var_type != Parameter && param_set != nothing
        _error("Can only use the keyword argument 'param_set' with infinite parameters.")
    elseif var_type != Infinite && param_refs != nothing
        _error("Can only use the keyword argument 'parameters' with infinite variables.")
    elseif var_type != Point && (inf_var_ref != nothing || param_values != nothing)
        _error("Can only use the keyword arguments 'infinite_var' and 'param_values' with point variables.")
    elseif var_type == Parameter
        if info.has_lb || info.has_lb || info.has_fix
            _error("Cannot specify bounds or fix a parameter.")
        end
        return InfOptParameter(param_set)
    elseif var_type == Infinite
        if param_refs isa ParameterRef
            param_refs = (param_refs, )
        end
        return InfiniteVariable(info, param_refs)
    elseif var_type == Point
        return PointVariable(info, inf_var_ref, param_values)
    else
        return GlobalVariable(info)
    end
end

"""
    JuMP.add_variable(model::InfiniteModel, v::InfOptVariable, name::String="")
Extend the `JuMP.add_variable` function to accomodate our new variable types.
"""
function JuMP.add_variable(model::InfiniteModel, v::InfOptVariable, name::String="")
    model.next_var_index += 1
    if isa(v, InfiniteVariable)
        vref = InfiniteVariableRef(model, model.next_var_index)
    elseif isa(v, PointVariable)
        vref = PointVariableRef(model, model.next_var_index)
    else
        vref = GlobalVariableRef(model, model.next_var_index)
    end
    model.vars[vref.index] = v
    if isa(vref, PointVariableRef)
        first_bracket = findfirst(isequal('('), JuMP.name(v.inf_var_ref))
        name = JuMP.name(v.inf_var_ref)[1:first_bracket-1]
    end
    JuMP.set_name(vref, name)
    if v.info.has_lb
        newset = MOI.GreaterThan(convert(Float64, v.info.lower_bound))
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_lower_bound_index(vref, JuMP.index(cref))
    end
    if v.info.has_ub
        newset = MOI.LessThan(convert(Float64, v.info.upper_bound))
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_upper_bound_index(vref, JuMP.index(cref))
    end
    if v.info.has_fix
        newset = MOI.EqualTo(convert(Float64, v.info.fixed_value))
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(vref, newset))
        JuMP.set_fix_index(vref, JuMP.index(cref))
    end
    if v.info.binary
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, MOI.ZeroOne()))
        JuMP.set_binary_index(vref, JuMP.index(cref))
    elseif v.info.integer
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, MOI.Integer()))
        JuMP.set_integer_index(vref, JuMP.index(cref))
    end
    return vref
end

"""
    JuMP.owner_model(vref::GeneralVariableRef)
Extend the `JuMP.owner_model` function for our new variable types.
"""
JuMP.owner_model(vref::GeneralVariableRef) = vref.model

"""
    JuMP.index(v::GeneralVariableRef::Int
Extent `JuMP.index` to return the index of a InfOpt variable.
"""
JuMP.index(v::GeneralVariableRef) = v.index

"""
    JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
Extend the `JuMP.delete` function to accomodate our new variable types.
"""
function JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
    @assert JuMP.is_valid(model, vref)
    delete!(model.vars, JuMP.index(vref))
    delete!(model.var_to_name, JuMP.index(vref))
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, vref::InfOptVariableRef)
Extend the `JuMP.is_valid` function to accomodate our new variable types.
"""
function JuMP.is_valid(model::InfiniteModel, vref::InfOptVariableRef)
        return (model === JuMP.owner_model(vref) && JuMP.index(vref) in keys(model.vars))
end

"""
    JuMP.num_variables(model::InfiniteModel)
Extend the `JuMP.num_variables` function to accomodate our new variable types.
"""
JuMP.num_variables(model::InfiniteModel) = length(model.vars)

# Internal functions
_variable_info(vref::InfOptVariableRef) = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
function _update_variable_info(vref::InfiniteVariableRef, info::JuMP.VariableInfo)
    param_refs = JuMP.owner_model(vref).vars[JuMP.index(vref)].param_refs
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, param_refs)
    return
end
function _update_variable_info(vref::PointVariableRef, info::JuMP.VariableInfo)
    inf_var_ref = JuMP.owner_model(vref).vars[JuMP.index(vref)].inf_var_ref
    param_values = JuMP.owner_model(vref).vars[JuMP.index(vref)].param_values
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = PointVariable(info, inf_var_ref, param_values)
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
    else
        @assert !JuMP.is_fixed(vref)
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_lower_bound_index(vref, JuMP.index(cref))
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
    else
        @assert !JuMP.is_fixed(vref)
        cref = JuMP.add_constraint(JuMP.owner_model(vref), JuMP.ScalarConstraint(vref, newset))
        JuMP.set_upper_bound_index(vref, JuMP.index(cref))
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
    end
    info = _variable_info(vref)
    _update_variable_info(vref, JuMP.VariableInfo(
        false, info.lower_bound, false, info.upper_bound, true, convert(Float64, value),
        info.has_start, info.start, info.binary, info.integer)
    )
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

"""
    JuMP.name(vref::InfOptVariableRef)
Extend the `JuMP.name` function to accomodate our new variable types.
"""
JuMP.name(vref::InfOptVariableRef) = JuMP.owner_model(vref).var_to_name[JuMP.index(vref)]

"""
    JuMP.set_name(vref::GlobalVariableRef, name::String)
Extend the `JuMP.set_name` function to accomodate point and global variables.
"""
function JuMP.set_name(vref::GlobalVariableRef, name::String)
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    JuMP.set_name(vref::PointVariableRef, name::String)
Extend the `JuMP.set_name` function to accomodate point and global variables.
"""
function JuMP.set_name(vref::PointVariableRef, root_name::String)
    param_values = JuMP.owner_model(vref).vars[JuMP.index(vref)].param_values
    name = string(root_name, param_values)
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    get_parameter_refs(vref::InfiniteVariableRef)
Get the `ParameterRef`(s) associated with the infinite variable `vref`.
"""
function get_parameter_refs(vref::InfiniteVariableRef)
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].param_refs
end

# Internal function
function _update_variable_param_refs(vref::InfiniteVariableRef, param_refs::Tuple)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, param_refs)
    return
end

"""
    set_parameter_refs(vref::InfiniteVariableRef, params::Tuple; force = false)
Set the `ParameterRef`(s) associated with the infinite variable `vref`.
"""
function set_parameter_refs(vref::InfiniteVariableRef, params::Tuple; force = false)
    _check_parameter_tuple(params)
    if get_parameter_refs(vref) != nothing && !force
        error("A parameter is already assigned to $vref, use the keyword argument
               force = true to overwrite it.")
    end
    _update_variable_param_refs(vref, params)
    return
end

"""
    add_parameter_ref(vref::InfiniteVariableRef, param::Union{ParameterRef, AbstractArray{<:ParameterRef}})
Add additional `ParameterRef` to be associated with the infinite
variable `vref`.
"""
function add_parameter_ref(vref::InfiniteVariableRef, param::Union{ParameterRef, AbstractArray{<:ParameterRef}})
    return set_parameter_refs(vref, (get_parameter_refs(vref)..., param), force = true)
end

"""
    JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
Extend the `JuMP.set_name` function to accomodate infinite variables.
"""
function JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
    param_refs = get_parameter_refs(vref)
    if isa(param_refs, Nothing) # no parameters given
        JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    else
        param_names = _get_root_names(param_refs)
        param_name_tuple = "("
        for i = 1:length(param_names)
            if i != length(param_names)
                param_name_tuple *= string(param_names[i], ", ")
            else
                param_name_tuple *= string(param_names[i], ")")
            end
        end
        var_name = string(root_name, param_name_tuple)
        JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = var_name
    end
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    JuMP.variable_by_name(model::InfiniteModel, name::String)
Extend the `JuMP.variable_by_name` function to accomodate `InfiniteModel` objects.
"""
function JuMP.variable_by_name(model::InfiniteModel, name::String)
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
        if isa(model.vars[index], InfiniteVariable)
            return InfiniteVariableRef(model, index)
        elseif isa(model.vars[index], PointVariable)
            return PointVariableRef(model, index)
        else
            return GlobalVariableRef(model, index)
        end
    end
    return
end

"""
    JuMP.all_variables(model::InfiniteModel)
Extend the `JuMP.all_variables` function to accomodate `InfiniteModel` objects.
"""
function JuMP.all_variables(model::InfiniteModel)
    vars_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model))
    indexes = sort([index for index in keys(model.vars)])
    counter = 1
    for index in indexes
        if isa(model.vars[index], InfiniteVariable)
            vars_list[counter] = InfiniteVariableRef(model, index)
        elseif isa(model.vars[index], PointVariable)
            vars_list[counter] = PointVariableRef(model, index)
        else
            vars_list[counter] = GlobalVariableRef(model, index)
        end
        counter += 1
    end
    return vars_list
end
