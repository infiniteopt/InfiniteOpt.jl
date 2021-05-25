"""
    transcription_model(model::InfiniteOpt.InfiniteModel)::JuMP.Model

Return the transcription model stored in `model` if that is what is stored in
`model.optimizer_model`.
"""
function transcription_model(model::InfiniteOpt.InfiniteModel)::JuMP.Model
    trans_model = InfiniteOpt.optimizer_model(model)
    if !is_transcription_model(trans_model)
        error("The model does not contain a transcription model.")
    end
    return trans_model
end

"""
    InfiniteOpt.build_optimizer_model!(model::InfiniteOpt.InfiniteModel,
                                       key::Val{:TransData};
                                       check_support_dims::Bool = true)::Nothing

Transcribe `model` and store it as a `TranscriptionModel` in the
`model.optimizer_model` field which can be accessed with `transcription_model`.
Ths clears the existing `TranscriptionModel` via
[`InfiniteOpt.clear_optimizer_model_build!`](@ref) and then builds a new one
using [`build_transcription_model!`](@ref).
"""
function InfiniteOpt.build_optimizer_model!(
    model::InfiniteOpt.InfiniteModel,
    key::Val{:TransData};
    check_support_dims::Bool = true
    )::Nothing
    # clear the optimzier model contents
    trans_model = InfiniteOpt.clear_optimizer_model_build!(model)
    # build the transcription model based on model
    build_transcription_model!(trans_model, model, 
                               check_support_dims = check_support_dims)
    # update the optimizer model status
    InfiniteOpt.set_optimizer_model_ready(model, true)
    return
end
