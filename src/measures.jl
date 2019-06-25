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
                             supports::Vector{<:AbstractArray{<:Number}};
                             name::String = "measure",
                             weight_function::Function = _w)
    supports = [convert(JuMP.Containers.SparseAxisArray, s) for s in supports]
    parameter_ref = convert(JuMP.Containers.SparseAxisArray, parameter_ref)
    return MultiDiscreteMeasureData(parameter_ref, coefficients, supports, name,
                                    weight_function)
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
function measure(expr::Union{InfiniteExpr, MeasureExpr}, data::DiscreteMeasureData)
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
    return _ReducedInfiniteRef(inf_model, index, ivref, supports)
end


# Implement functions for expanding measures into regular expressions
function _expand_measure_function(ivref::InfiniteVariableRef,
                                  data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                                  trans_model::JuMP.Model)
    prefs = parameter_refs(ivref)
    group = _groups((data.parameter_ref, ))[1]
    groups = _groups(prefs)
    if !(group in groups)
        aff = zero(JuMP.GenericAffExpr{Float64, InfiniteVariableRef})
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, ivref)
        end
    elseif length(prefs) == 1
        aff = zero(JuMP.GenericAffExpr{Float64, PointVariableRef})
        for i = 1:length(data.supports)
            pvref = _make_point_variable(ivref)
            support = (data.supports[i],)
            _update_point_mapping(trans_model, pvref, ivref, support)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pvref)
        end
    else
        tuple_loc = findfirst(isequal(group), groups)
        aff = zero(JuMP.GenericAffExpr{Float64, _ReducedInfiniteRef})
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

function _expand_measure_function(rvref::_ReducedInfiniteRef,
                                  data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                                  trans_model::JuMP.Model)
    orig_prefs = parameter_refs(rvref.original)
    prefs = ()
    for i = 1:length(orig_prefs)
        if !haskey(rvref.supports, i)
            prefs = (prefs..., orig_prefs[i])
        end
    end
    group = _groups((data.parameter_ref, ))[1]
    groups = _groups(prefs)
    if !(group in groups)
        aff = zero(JuMP.GenericAffExpr{Float64, _ReducedInfiniteRef})
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, rvref)
        end
    elseif length(prefs) == 1
        aff = zero(JuMP.GenericAffExpr{Float64, PointVariableRef})
        tuple_loc = findfirst(isequal(group), _groups(orig_prefs))
        for i = 1:length(data.supports)
            pvref = _make_point_variable(rvref.original)
            rvref.supports[tuple_loc] = data.supports[i]
            support = ()
            for j = 1:length(rvref.supports)
                support = (support..., rvref.supports[j])
            end
            _update_point_mapping(trans_model, pvref, rvref.original, support)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pvref)
        end
    else
        tuple_loc = findfirst(isequal(group), _groups(orig_prefs))
        aff = zero(JuMP.GenericAffExpr{Float64, _ReducedInfiniteRef})
        for i = 1:length(data.supports)
            supports = Dict(k => rvref.supports[k] for k in keys(rvref.supports))
            new_rvref = _make_reduced_variable(rvref.original, supports)
            support = data.supports[i]
            new_rvref.supports[tuple_loc] = data.supports[i]
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, new_rvref)
        end
    end
    return aff
end

function _expand_measure_function(vref::FiniteVariableRef,
                                  data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                                  trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, FiniteVariableRef})
    for i = 1:length(data.supports)
        coef = data.coefficients[i]
        weight = data.weight_function(data.supports[i])
        JuMP.add_to_expression!(aff, coef * weight, vref)
    end
    return aff
end

function _expand_measure_function(pref::ParameterRef,
                                  data::DiscreteMeasureData,
                                  trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, PointVariableRef})
    if data.parameter_ref == pref
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight * data.supports[i])
        end
    else
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pref)
        end
    end
    return aff
end

function _expand_measure_function(pref::ParameterRef,
                                  data::MultiDiscreteMeasureData,
                                  trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, FiniteVariableRef})
    pref_vals = collect(values(data.parameter_ref.data))
    pref_loc = findfirst(isequal(pref), pref_vals)
    if pref_loc != nothing
        key = collect(keys(data.parameter_ref.data))[pref_loc]
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight * data.supports[i][key])
        end
    else
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight, pref)
        end
    end
    return aff
end

function _expand_measure_function(expr::JuMP.GenericAffExpr,
                                  data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                                  trans_model::JuMP.Model) where {V}
    aff = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    for k in keys(expr.terms)
        var_expr = _expand_measure_function(k, data, trans_model)
        JuMP.add_to_expression!(aff, expr.terms[k], var_expr)
    end
    if expr.constant != 0
        for i = 1:length(data.supports)
            coef = data.coefficients[i]
            weight = data.weight_function(data.supports[i])
            JuMP.add_to_expression!(aff, coef * weight * expr.constant)
        end
    end
    return aff
end

function _expand_measure_function(mref::MeasureRef,
                                  data::Union{DiscreteMeasureData, MultiDiscreteMeasureData},
                                  trans_model::JuMP.Model) where {V}
    deeper_func = measure_function(mref)
    deeper_data = measure_data(mref)
    new_func = _expand_measure_function(deeper_func, deeper_data, trans_model)
    return _expand_measure_function(new_func, data, trans_model)
end

# Expand for more functions
