"""
    transcription_model(model::InfiniteOpt.InfiniteModel)::JuMP.Model

Return the transcription model stored in `model` if that is what is stored in
`model.optimizer_model`.
"""
function transcription_model(model::InfiniteOpt.InfiniteModel)::JuMP.Model
    if !is_transcription_model(InfiniteOpt.optimizer_model(model))
        error("The model does not contain a transcription model.")
    end
    return InfiniteOpt.optimizer_model(model)
end

"""
    InfiniteOpt.build_optimizer_model!(model::InfiniteOpt.InfiniteModel,
                                       key::Val{:TransData})

Transcribe `model` and store it as a `TranscriptionModel` in the
`model.optimizer_model` field which can be accessed with `transcription_model`.
"""
function InfiniteOpt.build_optimizer_model!(model::InfiniteOpt.InfiniteModel,
                                            key::Val{:TransData})
    # clear the optimzier model contents
    trans_model = InfiniteOpt.clear_optimizer_model_build!(model)
    # build the transcription model based on model
    _build_transcription_model!(trans_model, model)
    # update the optimizer model status
    InfiniteOpt.set_optimizer_model_ready(model, true)
    return
end
