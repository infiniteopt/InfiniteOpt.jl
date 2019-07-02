"""
    TranscriptionData
A DataType for storing the data mapping an `InfiniteModel` to a regular JuMP
`Model`.
**Fields**

"""
mutable struct TranscriptionData
    # Variable mapping
    infinite_to_vars::Dict{InfOpt.InfiniteVariableRef, Vector{JuMP.VariableRef}}
    global_to_var::Dict{InfOpt.GlobalVariableRef, JuMP.VariableRef}
    point_to_var::Dict{InfOpt.PointVariableRef, JuMP.VariableRef}

    # Variable support data
    infvar_to_supports::Dict{InfOpt.InfiniteVariableRef, Dict}

    # Constraint mapping
    infinite_to_constr::Dict{InfOpt.InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}
    measure_to_constr::Dict{InfOpt.MeasureConstraintRef, Vector{JuMP.ConstraintRef}}
    finite_to_constr::Dict{InfOpt.FiniteConstraintRef, JuMP.ConstraintRef}

    # Constraint support data
    infconstr_to_supports::Dict{InfOpt.InfiniteConstraintRef, Dict}
    measconstr_to_supports::Dict{InfOpt.MeasureConstraintRef, Dict}
    infconstr_to_params::Dict{InfOpt.InfiniteConstraintRef, Tuple}
    measconstr_to_params::Dict{InfOpt.MeasureConstraintRef, Tuple}

    # Default constructor
    function TranscriptionData()
        return new(Dict{InfOpt.InfiniteVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{InfOpt.GlobalVariableRef, JuMP.VariableRef}(),
                   Dict{InfOpt.PointVariableRef, JuMP.VariableRef}(),
                   Dict{InfOpt.InfiniteVariableRef, Dict}(),
                   Dict{InfOpt.InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfOpt.MeasureConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfOpt.FiniteConstraintRef, JuMP.ConstraintRef}(),
                   Dict{InfOpt.InfiniteConstraintRef, Dict}(),
                   Dict{InfOpt.MeasureConstraintRef, Dict}(),
                   Dict{InfOpt.InfiniteConstraintRef, Tuple}(),
                   Dict{InfOpt.MeasureConstraintRef, Tuple}())
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
    transcription_variable(model::JuMP.Model, vref::InfOpt.GlobalVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfOpt.GlobalVariableRef)
    !haskey(transcription_data(model).global_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).global_to_var[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfOpt.InfiniteVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfOpt.InfiniteVariableRef)
    !haskey(transcription_data(model).infinite_to_vars, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infinite_to_vars[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfOpt.PointVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfOpt.PointVariableRef)
    !haskey(transcription_data(model).point_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).point_to_var[vref]
end

"""
    InfOpt.supports(model::JuMP.Model, vref::InfOpt.InfiniteVariableRef)
Return the supports associated `vref` in the transcribed model.
"""
function InfOpt.supports(model::JuMP.Model, vref::InfOpt.InfiniteVariableRef)
    !haskey(transcription_data(model).infvar_to_supports, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infvar_to_supports[vref]
end

"""
    InfOpt.supports(vref::InfOpt.InfiniteVariableRef)
Return the supports associated `vref` in the transcribed model.
"""
function InfOpt.supports(vref::InfOpt.InfiniteVariableRef)
    model = InfOpt.optimizer_model(JuMP.owner_model(vref))
    !haskey(transcription_data(model).infvar_to_supports, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infvar_to_supports[vref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infinite_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infinite_to_constr[cref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measure_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).measure_to_constr[cref]
end

"""
    transcription_constraint(model::JuMP.Model, cref::InfOpt.FiniteConstraintRef)
Return the transcribed constraint reference(s) corresponding to `cref`.
"""
function transcription_constraint(model::JuMP.Model, cref::InfOpt.FiniteConstraintRef)
    !haskey(transcription_data(model).finite_to_constr, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).finite_to_constr[cref]
end

"""
    InfOpt.supports(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.supports(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_supports[cref]
end

"""
    InfOpt.supports(cref::InfOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.supports(cref::InfOpt.InfiniteConstraintRef)
    model = InfOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).infconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_supports[cref]
end

"""
    InfOpt.supports(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.supports(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have supports.")
    return transcription_data(model).measconstr_to_supports[cref]
end

"""
    InfOpt.supports(cref::InfOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.supports(cref::InfOpt.MeasureConstraintRef)
    model = InfOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).measconstr_to_supports, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have supports.")
    return transcription_data(model).measconstr_to_supports[cref]
end

"""
    InfOpt.parameter_refs(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.parameter_refs(model::JuMP.Model, cref::InfOpt.InfiniteConstraintRef)
    !haskey(transcription_data(model).infconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_params[cref]
end

"""
    InfOpt.parameter_refs(cref::InfOpt.InfiniteConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.parameter_refs(cref::InfOpt.InfiniteConstraintRef)
    model = InfOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).infconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model.")
    return transcription_data(model).infconstr_to_params[cref]
end

"""
    InfOpt.parameter_refs(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.parameter_refs(model::JuMP.Model, cref::InfOpt.MeasureConstraintRef)
    !haskey(transcription_data(model).measconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have parameter references.")
    return transcription_data(model).measconstr_to_params[cref]
end

"""
    InfOpt.parameter_refs(cref::InfOpt.MeasureConstraintRef)
Return the supports associated with `cref`.
"""
function InfOpt.parameter_refs(cref::InfOpt.MeasureConstraintRef)
    model = InfOpt.optimizer_model(JuMP.owner_model(cref))
    !haskey(transcription_data(model).measconstr_to_params, cref) && error("Constraint reference $cref not used in transcription model and/or is finite and doesn't have parameter references.")
    return transcription_data(model).measconstr_to_params[cref]
end
