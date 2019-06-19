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
function _get_model_from_expr(expr::JuMP.AbstractJuMPScalar)
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
    DiscreteMeasureData(parameter_ref::Union{ParameterRef, AbstractArray{<:ParameterRef}},
                        coefficients::Vector{<:Number},
                        supports::Union{Vector{<:Number}, Array{<:AbstractArray{<:Number}}};
                        name::String = "measure",
                        weight_function::Function = w(t) = 1)
Constructor for DiscreteMeasureData.
"""
function DiscreteMeasureData(parameter_ref::Union{ParameterRef, AbstractArray{<:ParameterRef}},
                             coefficients::Vector{<:Number},
                             supports::Union{Vector{<:Number}, Array{<:AbstractArray{<:Number}}};
                             name::String = "measure",
                             weight_function::Function = _w)
    return DiscreteMeasureData(name, weight_function, coefficients, supports, parameter_ref)
end

"""
    measure(expr::Union{InfiniteExpr, MeasureRef}, data::AbstractMeasureData)
Implement a measure in an expression in a similar fashion to the `sum` method in JuMP.
"""
function measure(expr::Union{InfiniteExpr, MeasureRef}, data::AbstractMeasureData)
    #TODO Check that variables have the correct parameters
    meas = Measure(expr, data)
    model = _get_model_from_expr(expr)
    if model == nothing
        error("Expression contains no variables.")
    end
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
