"""
    InfiniteOpt.map_value(fvref::InfiniteOpt.FiniteVariableRef,
                          key::Val{:TransData})

Map the value of the appropriate transcription variable in the transcription
model to `fvref`.
"""
function InfiniteOpt.map_value(fvref::InfiniteOpt.FiniteVariableRef,
                               key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fvref))
    return JuMP.value(transcription_variable(trans_model, fvref))
end

"""
    InfiniteOpt.map_value(ivref::InfiniteOpt.InfiniteVariableRef,
                          key::Val{:TransData})

Map the value of the appropriate transcription variable in the transcription
model to `ivref`.
"""
function InfiniteOpt.map_value(ivref::InfiniteOpt.InfiniteVariableRef,
                               key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(ivref))
    return JuMP.value.(transcription_variable(trans_model, ivref))
end
