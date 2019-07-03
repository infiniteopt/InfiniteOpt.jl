"""
    TranscriptionData
A DataType for storing the data mapping an `InfiniteModel` to a regular JuMP
`Model`.
**Fields**

"""
mutable struct TranscriptionData
    # Variable mapping
    infinite_to_vars::Dict{InfiniteOpt.InfiniteVariableRef, Vector{JuMP.VariableRef}}
    global_to_var::Dict{InfiniteOpt.GlobalVariableRef, JuMP.VariableRef}
    point_to_var::Dict{InfiniteOpt.PointVariableRef, JuMP.VariableRef}

    # Variable support data
    infvar_to_supports::Dict{InfiniteOpt.InfiniteVariableRef, Dict}

    # Constraint mapping
    infinite_to_constr::Dict{InfiniteOpt.InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}
    measure_to_constr::Dict{InfiniteOpt.MeasureConstraintRef, Vector{JuMP.ConstraintRef}}
    finite_to_constr::Dict{InfiniteOpt.FiniteConstraintRef, JuMP.ConstraintRef}

    # Constraint support data
    infconstr_to_supports::Dict{InfiniteOpt.InfiniteConstraintRef, Dict}
    measconstr_to_supports::Dict{InfiniteOpt.MeasureConstraintRef, Dict}
    infconstr_to_params::Dict{InfiniteOpt.InfiniteConstraintRef, Tuple}
    measconstr_to_params::Dict{InfiniteOpt.MeasureConstraintRef, Tuple}

    # Default constructor
    function TranscriptionData()
        return new(Dict{InfiniteOpt.InfiniteVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{InfiniteOpt.GlobalVariableRef, JuMP.VariableRef}(),
                   Dict{InfiniteOpt.PointVariableRef, JuMP.VariableRef}(),
                   Dict{InfiniteOpt.InfiniteVariableRef, Dict}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.MeasureConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.FiniteConstraintRef, JuMP.ConstraintRef}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef, Dict}(),
                   Dict{InfiniteOpt.MeasureConstraintRef, Dict}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef, Tuple}(),
                   Dict{InfiniteOpt.MeasureConstraintRef, Tuple}())
    end
end

"""
    TranscriptionModel(args...)
Return a JuMP `Model` with `TranscriptionData` included in the extension
data field.
"""
function TranscriptionModel(; kwargs...)
    model = JuMP.Model(; kwargs...)
    model.ext[:TransData] = TranscriptionData()
    return model
end

"""
    TranscriptionModel(args...)
Return a JuMP `Model` with `TranscriptionData` included in the extension
data field.
"""
function TranscriptionModel(optimizer_factory::JuMP.OptimizerFactory; kwargs...)
    model = JuMP.Model(optimizer_factory; kwargs...)
    model.ext[:TransData] = TranscriptionData()
    return model
end

"""
    is_transcription_model(model::JuMP.Model)
Return true if JuMP `Model` is a `TranscriptionModel` or false otherwise.
"""
function is_transcription_model(model::JuMP.Model)
    return haskey(model.ext, :TransData)
end

"""
    transcription_data(model::JuMP.Model)
Return the `TranscriptionData` from a `TranscriptionModel`. Error if it is not a
a `TranscriptionModel`.
"""
function transcription_data(model::JuMP.Model)
    !is_transcription_model(model) && error("Model is not a transcription model.")
    return model.ext[:TransData]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfiniteOpt.GlobalVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfiniteOpt.GlobalVariableRef)
    !haskey(transcription_data(model).global_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).global_to_var[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfiniteOpt.InfiniteVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfiniteOpt.InfiniteVariableRef)
    !haskey(transcription_data(model).infinite_to_vars, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infinite_to_vars[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfiniteOpt.PointVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfiniteOpt.PointVariableRef)
    !haskey(transcription_data(model).point_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).point_to_var[vref]
end

"""
    InfiniteOpt.supports(model::JuMP.Model, vref::InfiniteOpt.InfiniteVariableRef)
Return the supports associated `vref` in the transcribed model.
"""
function InfiniteOpt.supports(model::JuMP.Model, vref::InfiniteOpt.InfiniteVariableRef)
    !haskey(transcription_data(model).infvar_to_supports, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infvar_to_supports[vref]
end

"""
    InfiniteOpt.supports(vref::InfiniteOpt.InfiniteVariableRef)
Return the supports associated `vref` in the transcribed model.
"""
function InfiniteOpt.supports(vref::InfiniteOpt.InfiniteVariableRef)
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(vref))
    !haskey(transcription_data(model).infvar_to_supports, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infvar_to_supports[vref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infinite_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infinite_to_constr[cref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measure_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).measure_to_constr[cref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.FiniteConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfiniteOpt.FiniteConstraintRef)
    !haskey(transcription_data(model).finite_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).finite_to_constr[cref]
end

"""
    InfiniteOpt.supports(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.supports(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_supports[cref]
end

"""
    InfiniteOpt.supports(cref::InfiniteOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.supports(cref::InfiniteOpt.InfiniteConstraintRef)
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).infconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_supports[cref]
end

"""
    InfiniteOpt.supports(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.supports(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have supports.")
    return transcription_data(model).measconstr_to_supports[cref]
end

"""
    InfiniteOpt.supports(cref::InfiniteOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.supports(cref::InfiniteOpt.MeasureConstraintRef)
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).measconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have supports.")
    return transcription_data(model).measconstr_to_supports[cref]
end

"""
    InfiniteOpt.parameter_refs(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.parameter_refs(model::JuMP.Model, cref::InfiniteOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_params[cref]
end

"""
    InfiniteOpt.parameter_refs(cref::InfiniteOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.parameter_refs(cref::InfiniteOpt.InfiniteConstraintRef)
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).infconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_params[cref]
end

"""
    InfiniteOpt.parameter_refs(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.parameter_refs(model::JuMP.Model, cref::InfiniteOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have parameter references.")
    return transcription_data(model).measconstr_to_params[cref]
end

"""
    InfiniteOpt.parameter_refs(cref::InfiniteOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfiniteOpt.parameter_refs(cref::InfiniteOpt.MeasureConstraintRef)
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).measconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have parameter references.")
    return transcription_data(model).measconstr_to_params[cref]
end
