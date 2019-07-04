# Extend Base.copy for new variable types
Base.copy(v::MeasureRef, new_model::InfiniteModel) = MeasureRef(new_model, v.index)

# Parse the string for displaying a measure
function _make_meas_name(meas::Measure)
    return string(meas.data.name, "(", JuMP.function_string(JuMP.REPLMode, meas.func), ")")
end

# Used to update the model.var_to_meas field
function _update_var_meas_mapping(vrefs::Vector{<:GeneralVariableRef}, mindex::Int)
    for vref in vrefs
        model = JuMP.owner_model(vref)
        if isa(vref, InfOptVariableRef)
            if haskey(model.var_to_meas, JuMP.index(vref))
                push!(model.var_to_meas[JuMP.index(vref)], mindex)
            else
                model.var_to_meas[JuMP.index(vref)] = [mindex]
            end
        elseif isa(vref, ParameterRef)
            if haskey(model.param_to_meas, JuMP.index(vref))
                push!(model.param_to_meas[JuMP.index(vref)], mindex)
            else
                model.param_to_meas[JuMP.index(vref)] = [mindex]
            end
        end
    end
    return
end

# Set a default weight function
_w(t) = 1

"""
    DiscreteMeasureData(parameter_ref::ParameterRef, coefficients::Vector{<:Number},
                        supports::Vector{<:Number}; name::String = "measure",
                        weight_function::Function = w(t) = 1)
Constructor for DiscreteMeasureData.
"""
function DiscreteMeasureData(parameter_ref::ParameterRef,
                             coefficients::Vector{<:Number},
                             supports::Vector{<:Number};
                             name::String = "measure",
                             weight_function::Function = _w)
    return DiscreteMeasureData(parameter_ref, coefficients, supports, name,
                               weight_function)
end

"""
    DiscreteMeasureData(parameter_ref::AbstractArray{<:ParameterRef},
                        coefficients::Vector{<:Number},
                        supports::Vector{<:AbstractArray{<:Number}};
                        name::String = "measure",
                        weight_function::Function = w(t) = 1)
Constructor for MultiDiscreteMeasureData.
"""
function DiscreteMeasureData(parameter_ref::AbstractArray{<:ParameterRef},
                             coefficients::Vector{<:Number},
                             supports::Vector{<:AbstractArray};
                             name::String = "measure",
                             weight_function::Function = _w)
    supports = [convert(JuMP.Containers.SparseAxisArray, s) for s in supports]
    parameter_ref = convert(JuMP.Containers.SparseAxisArray, parameter_ref)
    return MultiDiscreteMeasureData(parameter_ref, coefficients, supports, name,
                                    weight_function)
end

"""
    add_measure(model::InfiniteModel, v::Measure)
Add a measure to the `InfiniteModel` object in an analagous way to `JuMP.add_variable`.
"""
function add_measure(model::InfiniteModel, meas::Measure)
    model.next_meas_index += 1
    index = model.next_meas_index
    vrefs = _all_function_variables(meas.func)
    _update_var_meas_mapping(vrefs, index)
    mref = MeasureRef(model, model.next_meas_index)
    model.measures[mref.index] = meas
    JuMP.set_name(mref, _make_meas_name(meas))
    return mref
end

# Parse the model pertaining to an expression
function _model_from_expr(expr::JuMP.AbstractJuMPScalar)
    all_vrefs = _all_function_variables(expr)
    if length(all_vrefs) > 0
        return JuMP.owner_model(all_vrefs[1])
    else
        return
    end
end

# check a measure function for a particular parameter
function _has_parameter(vrefs::Vector{<:GeneralVariableRef}, pref::ParameterRef)
    if _has_variable(vrefs, pref)
        return true
    end
    model = JuMP.owner_model(pref)
    relavent_vindices = model.param_to_vars[JuMP.index(pref)]
    relavent_vrefs = [InfiniteVariableRef(model, vindex) for vindex in relavent_vindices]
    for vref in relavent_vrefs
        if _has_variable(vrefs, vref)
            return true
        end
    end
    return false
end

function _check_has_parameter(expr::Union{InfiniteExpr, MeasureExpr},
                              pref::ParameterRef)
    vrefs = _all_function_variables(expr)
    if !_has_parameter(vrefs, pref)
        error("Measure expression is not parameterized by the parameter " *
              "specified in the measure data.")
    end
    return
end

function _check_has_parameter(expr::Union{InfiniteExpr, MeasureExpr},
                              pref::JuMP.Containers.SparseAxisArray{<:ParameterRef})
    vrefs = _all_function_variables(expr)
    for key in keys(pref.data)
        if !_has_parameter(vrefs, pref.data[key])
            error("Measure expression is not parameterized by the parameter " *
                  "specified in the measure data.")
        end
    end
    return
end

# Internal function for adding measure data supports directly to the parameter supports
function _add_supports_to_parameters(pref::ParameterRef, supports::Vector{<:Number})
    add_supports(pref, supports)
    return
end

function _add_supports_to_parameters(pref::JuMP.Containers.SparseAxisArray{<:ParameterRef},
                                     supports::Array{<:JuMP.Containers.SparseAxisArray{<:Number}})
    for i = 1:length(supports)
        for key in keys(pref.data)
            add_supports(pref.data[key], supports[i].data[key])
        end
    end
    return
end

"""
    measure(expr::Union{InfiniteExpr, MeasureRef}, data::AbstractMeasureData)
Implement a measure in an expression in a similar fashion to the `sum` method in JuMP.
"""
function measure(expr::Union{InfiniteExpr, MeasureExpr}, data::AbstractMeasureData)
    pref = data.parameter_ref
    _check_has_parameter(expr, pref)
    meas = Measure(expr, data)
    model = _model_from_expr(expr)
    if model == nothing # --> Might not be necessary with above checks
        error("Expression contains no variables.")
    end
    _add_supports_to_parameters(pref, data.supports)
    return add_measure(model, meas)
end

"""
    JuMP.name(mref::MeasureRef)
Extend the `JuMP.name` function to accomodate measure references.
"""
JuMP.name(mref::MeasureRef) = mref.model.meas_to_name[mref.index]

"""
    JuMP.set_name(mref::MeasureRef, name::String)
Extend the `JuMP.set_name` function to accomodate measure references.
"""
function JuMP.set_name(mref::MeasureRef, name::String)
    JuMP.owner_model(mref).meas_to_name[JuMP.index(mref)] = name
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, mref::MeasureRef)
Extend the `JuMP.is_valid` function to accomodate measures.
"""
function JuMP.is_valid(model::InfiniteModel, mref::MeasureRef)
    return (model === JuMP.owner_model(mref) && JuMP.index(mref) in keys(model.measures))
end

"""
    JuMP.delete(model::InfiniteModel, mref::MeasureRef)
Extend the `JuMP.delete` function to accomodate measures
"""
function JuMP.delete(model::InfiniteModel, mref::MeasureRef)
    @assert JuMP.is_valid(model, mref)
    set_optimizer_model_status(model, false)
    delete!(model.measures, JuMP.index(mref))
    delete!(model.meas_to_name, JuMP.index(mref))
    return
end

"""
    measure_function(mref::MeasureRef)
Return the function associated with a `mref`
"""
function measure_function(mref::MeasureRef)
    return JuMP.owner_model(mref).measures[JuMP.index(mref)].func
end

"""
    measure_data(mref::MeasureRef)
Return the data associated with a `mref`
"""
function measure_data(mref::MeasureRef)
    return JuMP.owner_model(mref).measures[JuMP.index(mref)].data
end

# Helper function for making place holder point variables
function _make_point_variable(ivref::InfiniteVariableRef)
    inf_model = JuMP.owner_model(ivref)
    index = inf_model.next_var_index += 1
    return PointVariableRef(inf_model, index)
end

# Helper function for making place holder infinite variables
function _make_reduced_variable(ivref::InfiniteVariableRef, removed_index::Int,
                                support::Union{Number, JuMP.Containers.SparseAxisArray{<:Number}})
    inf_model = JuMP.owner_model(ivref)
    index = inf_model.next_var_index += 1
    return _ReducedInfiniteRef(inf_model, index, ivref, Dict(removed_index => support))
end

function _make_reduced_variable(ivref::InfiniteVariableRef, supports::Dict)
    inf_model = JuMP.owner_model(ivref)
    index = inf_model.next_var_index += 1
    return _ReducedInfiniteRef(inf_model, index, ivref, copy(supports))
end

## Implement functions for expanding measures into regular expressions
# InfiniteVariableRef
function _expand_measure(ivref::InfiniteVariableRef,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model)
    # figure out the parameter groups
    prefs = parameter_refs(ivref)
    group = _groups((data.parameter_ref, ))[1]
    groups = _groups(prefs)
    # treat variable as constant if doesn't have measure parameter
    if !(group in groups)
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, ivref)
        end
    # convert variable into point variables if its only parameter is the the measure parameter
    elseif length(prefs) == 1
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i = 1:length(data.supports)
            pvref = _make_point_variable(ivref)
            support = (data.supports[i],)
            TranscriptionOpt._update_point_mapping(trans_model, pvref, ivref, support)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pvref)
        end
    # make reduced variables if the variable contains other parameters
    else
        tuple_loc = findfirst(isequal(group), groups)
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i = 1:length(data.supports)
            support = data.supports[i]
            rvref = _make_reduced_variable(ivref, tuple_loc, support)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, rvref)
        end
    end
    return aff
end

# _ReducedInfiniteRef
function _expand_measure(rvref::_ReducedInfiniteRef,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model)
    # figure out the parameters used by the reduced infinite variable
    orig_prefs = parameter_refs(rvref.original)
    prefs = parameter_refs(rvref)
    # figure out the parameter groups
    group = _groups((data.parameter_ref, ))[1]
    groups = _groups(prefs)
    # treat variable as constant if doesn't have measure parameter
    if !(group in groups)
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, rvref)
        end
    # convert variable into point variables if its only parameter is the the measure parameter
    elseif length(prefs) == 1
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        tuple_loc = findfirst(isequal(group), _groups(orig_prefs))
        for i = 1:length(data.supports)
            pvref = _make_point_variable(rvref.original)
            rvref.supports[tuple_loc] = data.supports[i]
            support = ()
            for j = 1:length(rvref.supports)
                support = (support..., rvref.supports[j])
            end
            TranscriptionOpt._update_point_mapping(trans_model, pvref, rvref.original, support)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pvref)
        end
    # make reduced variables if the variable contains other parameters
    else
        tuple_loc = findfirst(isequal(group), _groups(orig_prefs))
        aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i = 1:length(data.supports)
            new_rvref = _make_reduced_variable(rvref.original, rvref.supports)
            new_rvref.supports[tuple_loc] = data.supports[i]
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, new_rvref)
        end
    end
    return aff
end

# FiniteVariableRef
function _expand_measure(vref::FiniteVariableRef,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    # treat the variable as a contant
    for i = 1:length(data.supports)
        coef = data.coefficients[i]
        weight = data.weight_function(data.supports[i])
        JuMP.add_to_expression!(aff, coef * weight, vref)
    end
    return aff
end

# ParameterRef with scalar data
function _expand_measure(pref::ParameterRef,
                         data::DiscreteMeasureData,
                         trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    # replace the parameter with its value if it is the measure parameter
    if data.parameter_ref == pref
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight * data.supports[i])
        end
    # treat the parameter as a constant otherwise
    else
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pref)
        end
    end
    return aff
end

# ParameterRef with vector data
function _expand_measure(pref::ParameterRef,
                         data::MultiDiscreteMeasureData,
                         trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    # determine if pref is part of the measure parameters
    prefs = collect(values(data.parameter_ref.data))
    pref_index = findfirst(isequal(pref), prefs)
    # replace the parameter with its value if it is the measure parameter
    if pref_index != nothing
        key = collect(keys(data.parameter_ref.data))[pref_index]
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight * data.supports[i][key])
        end
    # treat the parameter as a constant otherwise
    else
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pref)
        end
    end
    return aff
end

# GenericAffExpr
function _expand_measure(expr::JuMP.GenericAffExpr,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model) where {V}
    # need to use a quadratic expression in case contains measures with quadratic expressions
    quad = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
    # expand each variable independently and add all together
    for (var, coef) in expr.terms
        var_expr = _expand_measure(var, data, trans_model)
        JuMP.add_to_expression!(quad, coef, var_expr)
    end
    # expand over the cantant
    if expr.constant != 0
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(quad, coef * weight * expr.constant)
        end
    end
    if length(quad.terms) == 0
        return quad.aff
    else
        return quad
    end
end

# GenericQuadExpr
function _expand_measure(expr::JuMP.GenericQuadExpr,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model) where {V}
    quad = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
    # convert the GenericAffExpr
    quad.aff = _expand_measure(expr.aff, data, trans_model)
    for (pair, coef) in expr.terms
        var_a = pair.a
        var_b = pair.b
        # expand on both variables
        expr_a = _expand_measure(var_a, data, trans_model)
        expr_b = _expand_measure(var_b, data, trans_model)
        vars_a = collect(keys(expr_a.terms))
        vars_b = collect(keys(expr_b.terms))
        # combine both variable expressions using the coefficients from one of them
        if length(vars_a) == length(vars_b)
            # are same length therefore have same coefficients
            for i = 1:length(vars_a)
                JuMP.add_to_expression!(quad, coef * expr_a.terms[vars_a[i]], vars_a[i], vars_b[i])
            end
        elseif length(vars_a) == 1
            # var_a was effectively a constant and var_b was't
            for i = 1:length(vars_b)
                JuMP.add_to_expression!(quad, coef * expr_b.terms[vars_b[i]], vars_a[1], vars_b[i])
            end
        else
            # var_b was effectively a constant and var_a was't
            for i = 1:length(vars_a)
                JuMP.add_to_expression!(quad, coef * expr_a.terms[vars_a[i]], vars_a[i], vars_b[1])
            end
        end
    end
    return quad
end

# MeasureRef
function _expand_measure(mref::MeasureRef,
                         data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                         trans_model::JuMP.Model) where {V}
    # determine function and data of the inner measure
    deeper_func = measure_function(mref)
    deeper_data = measure_data(mref)
    # expand the inner measure (note this is recursive for nested measures)
    new_func = _expand_measure(deeper_func, deeper_data, trans_model)
    # expand current level with the inner measure now expanded
    return _expand_measure(new_func, data, trans_model)
end
