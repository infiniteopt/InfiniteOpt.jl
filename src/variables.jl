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

# Make check tuple check functions
function _check_parameter_tuple(_error::Function, prefs::Tuple)
    types = [typeof(pref) for pref in prefs]
    num_params = length(types)
    valid_types = zeros(Bool, num_params)
    for i = 1:num_params
        if types[i] == ParameterRef || types[i] <: AbstractArray{<:ParameterRef}
            valid_types[i] = true
        end
    end
    if sum(valid_types) != num_params
        _error("Invalid parameter type(s) given.")
    end
    return
end

function _make_formatted_tuple(prefs::Tuple)
    converted_prefs = ()
    for pref in prefs
        if isa(pref, ParameterRef) || isa(pref, Number)
            converted_prefs = (converted_prefs..., pref)
        else
            converted_prefs = (converted_prefs...,
                               convert(JuMP.Containers.SparseAxisArray, pref))
        end
    end
    return converted_prefs
end

function _check_tuple_names(_error::Function, prefs::Tuple)
    valid_elements = [_only_one_name(prefs[i]) for i = 1:length(prefs)]
    if sum(valid_elements) != length(prefs)
        _error("Each parameter tuple element must have contain only one infinite parameter name.")
    end
    root_names = _get_root_names(prefs)
    if length(unique(root_names)) != length(root_names)
        _error("Cannot double specify infinite parameter references.")
    end
    return
end

function _check_tuple_shape(_error::Function, infinite_variable_ref::InfiniteVariableRef, values::Tuple)
    prefs = get_parameter_refs(infinite_variable_ref)
    if length(prefs) != length(values)
        _error("The dimensions of the infinite parameter values must match those defined for the infinite variable.")
    end
    for i = 1:length(values)
        if isa(prefs[i], ParameterRef) && !(typeof(values[i]) <: Number)
            _error("The dimensions and array type of the infinite parameter values must match those defined for the infinite variable.")
        elseif isa(prefs[i], JuMP.Containers.SparseAxisArray) && !isa(values[i], JuMP.Containers.SparseAxisArray)
            _error("The dimensions and array type of the infinite parameter values must match those defined for the infinite variable.")
        elseif isa(prefs[i], JuMP.Containers.SparseAxisArray)
            if keys(prefs[i].data) != keys(values[i].data)
                _error("Index keys of infinite parameter values don't match those defined for the infinite variable.")
            end
        end
    end
    return
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, test_arg::Int;
                        param_set::Union{AbstractInfiniteSet, Vector{AbstractInfiniteSet}} = EmptySet(),
                        extra_kw_args...)
Extend the `JuMP.build_variable` function to accomodate our new variable types.
"""
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, var_type::Symbol;
                             parameter_refs::Union{ParameterRef, AbstractArray{<:ParameterRef}, Tuple, Nothing} = nothing,
                             infinite_variable_ref::Union{InfiniteVariableRef, Nothing} = nothing,
                             parameter_values::Union{Number, AbstractArray{<:Number}, Tuple, Nothing} = nothing,
                             error::Union{Function, Nothing} = nothing,
                             extra_kw_args...)
    if error == nothing
        _error = error # replace with macro error function
    end
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    if !(var_type in [Infinite, Point, Global])
        _error("Unrecognized variable type $var_type, should be Inf, Point, or Global.")
    end
    if var_type != Infinite && parameter_refs != nothing
        _error("Can only use the keyword argument 'parameter_refs' with infinite variables.")
    elseif var_type != Point && (infinite_variable_ref != nothing || parameter_values != nothing)
        _error("Can only use the keyword arguments 'infinite_var' and 'parameter_values' with point variables.")
    elseif var_type == Infinite
        if parameter_refs == nothing
            _error("Parameter references not specified, use the var(params...) syntax or the parameter_refs keyword argument.")
        end
        if !isa(parameter_refs, Tuple)
            parameter_refs = (parameter_refs, )
        end
        _check_parameter_tuple(_error, parameter_refs)
        parameter_refs = _make_formatted_tuple(parameter_refs)
        _check_tuple_names(_error, parameter_refs)
        return InfiniteVariable(info, parameter_refs)
    elseif var_type == Point
        if parameter_values == nothing || infinite_variable_ref == nothing
            _error("Must specify the infinite variable and the values of its infinite parameters")
        end
        if !isa(parameter_values, Tuple)
            parameter_values = (parameter_values, )
        end
        parameter_values = _make_formatted_tuple(parameter_values)
        _check_tuple_shape(_error, infinite_variable_ref, parameter_values)
        return PointVariable(info, infinite_variable_ref, parameter_values)
    else
        return GlobalVariable(info)
    end
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

"""
    JuMP.add_variable(model::InfiniteModel, v::InfOptVariable, name::String="")
Extend the `JuMP.add_variable` function to accomodate our new variable types.
"""
function JuMP.add_variable(model::InfiniteModel, v::InfOptVariable, name::String="")
    model.next_var_index += 1
    if isa(v, InfiniteVariable)
        vref = InfiniteVariableRef(model, model.next_var_index)
        _update_param_var_mapping(vref, v.parameter_refs)
    elseif isa(v, PointVariable)
        vref = PointVariableRef(model, model.next_var_index)
    else
        vref = GlobalVariableRef(model, model.next_var_index)
    end
    model.vars[JuMP.index(vref)] = v
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
    model.var_in_objective[JuMP.index(vref)] = false
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
    used_by_constraint(vref::InfOptVariableRef)::Bool
Return a Boolean indicating if `vref` is used by a constraint.
"""
function used_by_constraint(vref::InfOptVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).var_to_constrs, JuMP.index(vref))
end

"""
    used_by_measure(vref::InfOptVariableRef)::Bool
Return a Boolean indicating if `vref` is used by a measure.
"""
function used_by_measure(vref::InfOptVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).var_to_meas, JuMP.index(vref))
end

"""
    used_by_objective(vref::InfOptVariableRef)::Bool
Return a Boolean indicating if `vref` is used by the objective.
"""
function used_by_objective(vref::InfOptVariableRef)::Bool
    return JuMP.owner_model(vref).var_in_objective[JuMP.index(vref)]
end

"""
    is_used(vref::InfOptVariableRef)::Bool
Return a Boolean indicating if `vref` is used in the model.
"""
function is_used(vref::InfOptVariableRef)::Bool
    return used_by_measure(vref) || used_by_constraint(vref) || used_by_objective(vref)
end

"""
    JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
Extend the `JuMP.delete` function to accomodate our new variable types.
"""
function JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
    @assert JuMP.is_valid(model, vref)
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
    if used_by_measure(vref)
        for mindex in model.var_to_meas[JuMP.index(vref)]
            if isa(model.measures[mindex].func, InfOptVariableRef)
                data = model.measures[mindex].data
                model.measures[mindex] = Measure(zero(JuMP.AffExpr), data)
            else
                _remove_variable(model.measures[mindex].func, vref)
            end
            measure = model.measures[mindex]
            mref = MeasureRef(model, mindex)
            JuMP.set_name(mref, _make_meas_name(measure))
        end
        delete!(model.var_to_meas, JuMP.index(vref))
    end
    if used_by_constraint(vref)
        for cindex in model.var_to_constrs[JuMP.index(vref)]
            if isa(model.constrs[cindex].func, InfOptVariableRef)
                set = model.constrs[cindex].set
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr), set)
            else
                _remove_variable(model.constrs[cindex].func, vref)
            end
        end
        delete!(model.var_to_constrs, JuMP.index(vref))
    end
    if used_by_objective(vref)
        if isa(model.objective_function, InfOptVariableRef)
            model.objective_function = zero(JuMP.AffExpr)
        else
            _remove_variable(model.objective_function, vref)
        end
    end
    if isa(vref, InfiniteVariableRef)
        all_prefs = _list_parameter_refs(get_parameter_refs(vref))
        for pref in all_prefs
            filter!(e -> e != JuMP.index(vref), model.param_to_vars[JuMP.index(pref)])
            if length(model.param_to_vars[JuMP.index(pref)]) == 0
                delete!(model.param_to_vars, JuMP.index(pref))
            end
        end
    end
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
    get_infinite_variable_ref(vref::PointVariableRef)
Get the `InfiniteVariableRef`associated with the point variable `vref`.
"""
function get_infinite_variable_ref(vref::PointVariableRef)
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].infinite_variable_ref
end

"""
    JuMP.set_name(vref::PointVariableRef, name::String)
Extend the `JuMP.set_name` function to accomodate point and global variables.
"""
function JuMP.set_name(vref::PointVariableRef, name::String)
    if length(name) == 0
        inf_var_ref = get_infinite_variable_ref(vref::PointVariableRef)
        first_bracket = findfirst(isequal('('), JuMP.name(inf_var_ref))
        name = JuMP.name(inf_var_ref)[1:first_bracket-1]
        # TODO do something about SparseAxisArrays
        parameter_values = JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_values
        name = string(name, parameter_values)
    end
    JuMP.owner_model(vref).var_to_name[JuMP.index(vref)] = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

"""
    get_parameter_refs(vref::InfiniteVariableRef)
Get the `ParameterRef`(s) associated with the infinite variable `vref`.
"""
function get_parameter_refs(vref::InfiniteVariableRef)
    return JuMP.owner_model(vref).vars[JuMP.index(vref)].parameter_refs
end

# Internal function
function _update_variable_param_refs(vref::InfiniteVariableRef, prefs::Tuple)
    info = JuMP.owner_model(vref).vars[JuMP.index(vref)].info
    JuMP.owner_model(vref).vars[JuMP.index(vref)] = InfiniteVariable(info, prefs)
    return
end

"""
    set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)
Set the `ParameterRef`(s) associated with the infinite variable `vref`.
"""
function set_parameter_refs(vref::InfiniteVariableRef, prefs::Tuple)
    _check_parameter_tuple(error, prefs)
    prefs = _make_formatted_tuple(prefs)
    _check_tuple_names(error, prefs)
    _update_variable_param_refs(vref, prefs)
    return
end

"""
    add_parameter_ref(vref::InfiniteVariableRef, pref::Union{ParameterRef, AbstractArray{<:ParameterRef}})
Add additional `ParameterRef` to be associated with the infinite
variable `vref`.
"""
function add_parameter_ref(vref::InfiniteVariableRef, pref::Union{ParameterRef, AbstractArray{<:ParameterRef}})
    return set_parameter_refs(vref, (get_parameter_refs(vref)..., pref))
end

"""
    JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
Extend the `JuMP.set_name` function to accomodate infinite variables.
"""
function JuMP.set_name(vref::InfiniteVariableRef, root_name::String)
    if length(root_name) == 0
        root_name = "noname"
    end
    prefs = get_parameter_refs(vref)
    param_names = _get_root_names(prefs)
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
    vrefs_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model))
    indexes = sort([index for index in keys(model.vars)])
    counter = 1
    for index in indexes
        if isa(model.vars[index], InfiniteVariable)
            vrefs_list[counter] = InfiniteVariableRef(model, index)
        elseif isa(model.vars[index], PointVariable)
            vrefs_list[counter] = PointVariableRef(model, index)
        else
            vrefs_list[counter] = GlobalVariableRef(model, index)
        end
        counter += 1
    end
    return vrefs_list
end
