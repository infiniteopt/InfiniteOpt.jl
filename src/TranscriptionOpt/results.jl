"""
    InfOpt.map_value(fvref::InfOpt.FiniteVariableRef, key::Val{:TransData})
Map the value of the appropriate transcription variable in the transcription
model to `fvref`.
"""
function InfOpt.map_value(fvref::InfOpt.FiniteVariableRef, key::Val{:TransData})
    trans_model = InfOpt.optimizer_model(JuMP.owner_model(fvref))
    return JuMP.value(transcription_variable(trans_model, fvref))
end

"""
    InfOpt.map_value(ivref::InfOpt.InfiniteVariableRef, key::Val{:TransData})
Map the value of the appropriate transcription variable in the transcription
model to `ivref`.
"""
function InfOpt.map_value(ivref::InfOpt.InfiniteVariableRef, key::Val{:TransData})
    trans_model = InfOpt.optimizer_model(JuMP.owner_model(ivref))
    return JuMP.value.(transcription_variable(trans_model, ivref))
end
