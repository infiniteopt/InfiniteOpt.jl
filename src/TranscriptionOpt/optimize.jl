"""
    transcription_model(model::InfOpt.InfiniteModel)
Return the transcription model stored in `model` if that is what is stored in
`model.optimizer_model`.
"""
function transcription_model(model::InfOpt.InfiniteModel)
    !is_transcription_model(InfOpt.optimizer_model(model)) && error("The model does not contain a transcription model.")
    return InfOpt.optimizer_model(model)
end

"""
    InfOpt.build_optimizer_model!(model::InfOpt.InfiniteModel, key::Val{:TransData})
Transcribe `model` and store it as a `TranscriptionModel` in the
`model.optimizer_model` field which can be accessed with `transcription_model`.
"""
function InfOpt.build_optimizer_model!(model::InfOpt.InfiniteModel, key::Val{:TransData})
    mode = JuMP.backend(InfOpt.optimizer_model(model)).mode
    trans_model = TranscriptionModel(model, caching_mode = mode)
    if !isa(model.optimizer_factory, Nothing)
        bridge_constrs = JuMP.bridge_constraints(model)
        JuMP.set_optimizer(trans_model, model.optimizer_factory,
                           bridge_constraints = bridge_constrs)
    end
    InfOpt.set_optimizer_model(model, trans_model)
    InfOpt.set_optimizer_model_status(model, true)
    return
end
