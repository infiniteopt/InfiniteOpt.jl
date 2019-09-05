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

Map the value of the appropriate transcription variables in the transcription
model to `ivref`.
"""
function InfiniteOpt.map_value(ivref::InfiniteOpt.InfiniteVariableRef,
                               key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(ivref))
    return JuMP.value.(transcription_variable(trans_model, ivref))
end

"""
    InfiniteOpt.map_value(fcref::InfiniteOpt.FiniteConstraintRef,
                          key::Val{:TransData})

Map the value of the appropriate transcription constraint function in the
transcription model to `fcref`.
"""
function InfiniteOpt.map_value(fcref::InfiniteOpt.FiniteConstraintRef,
                               key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fcref))
    return JuMP.value(transcription_constraint(trans_model, fcref))
end

"""
    InfiniteOpt.map_value(icref::InfiniteOpt.GeneralConstraintRef,
                          key::Val{:TransData})

Map the value of the appropriate transcription constraint functions in the
transcription model to `icref`.
"""
function InfiniteOpt.map_value(icref::InfiniteOpt.GeneralConstraintRef,
                               key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(icref))
    return JuMP.value.(transcription_constraint(trans_model, icref))
end

"""
    InfiniteOpt.map_optimizer_index(fvref::InfiniteOpt.FiniteVariableRef,
                                    key::Val{:TransData})

Map the optimizer model index of the appropriate transcription variable in the
transcription model to `fvref`.
"""
function InfiniteOpt.map_optimizer_index(fvref::InfiniteOpt.FiniteVariableRef,
                                         key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fvref))
    return JuMP.optimizer_index(transcription_variable(trans_model, fvref))
end

"""
    InfiniteOpt.map_optimizer_index(ivref::InfiniteOpt.InfiniteVariableRef,
                                    key::Val{:TransData})

Map the optimizer model index of the appropriate transcription variables in the
transcription model to `ivref`.
"""
function InfiniteOpt.map_optimizer_index(ivref::InfiniteOpt.InfiniteVariableRef,
                                         key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(ivref))
    return JuMP.optimizer_index.(transcription_variable(trans_model, ivref))
end

"""
    InfiniteOpt.map_optimizer_index(fcref::InfiniteOpt.FiniteConstraintRef,
                                    key::Val{:TransData})

Map the optimizer model index of the appropriate transcription constraint in the
transcription model to `fcref`.
"""
function InfiniteOpt.map_optimizer_index(fcref::InfiniteOpt.FiniteConstraintRef,
                                         key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fcref))
    return JuMP.optimizer_index(transcription_constraint(trans_model, fcref))
end

"""
    InfiniteOpt.map_optimizer_index(cref::InfiniteOpt.GeneralConstraintRef,
                                    key::Val{:TransData})

Map the optimizer model index of the appropriate transcription constraints in the
transcription model to `cref`.
"""
function InfiniteOpt.map_optimizer_index(cref::InfiniteOpt.GeneralConstraintRef,
                                         key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return JuMP.optimizer_index.(transcription_constraint(trans_model, cref))
end

"""
    InfiniteOpt.map_dual(fcref::InfiniteOpt.FiniteConstraintRef,
                         key::Val{:TransData})

Map the dual of the appropriate transcription constraint in the
transcription model to `fcref`.
"""
function InfiniteOpt.map_dual(fcref::InfiniteOpt.FiniteConstraintRef,
                              key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fcref))
    return JuMP.dual(transcription_constraint(trans_model, fcref))
end

"""
    InfiniteOpt.map_dual(icref::InfiniteOpt.InfiniteConstraintRef,
                         key::Val{:TransData})

Map the duals of the appropriate transcription constraints in the
transcription model to `icref`.
"""
function InfiniteOpt.map_dual(icref::InfiniteOpt.InfiniteConstraintRef,
                              key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(icref))
    return JuMP.dual.(transcription_constraint(trans_model, icref))
end

"""
    InfiniteOpt.map_shadow_price(fcref::InfiniteOpt.FiniteConstraintRef,
                                 key::Val{:TransData})

Map the shadow price of the appropriate transcription constraint in the
transcription model to `fcref`.
"""
function InfiniteOpt.map_shadow_price(fcref::InfiniteOpt.FiniteConstraintRef,
                                      key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(fcref))
    return JuMP.shadow_price(transcription_constraint(trans_model, fcref))
end

"""
    InfiniteOpt.map_shadow_price(icref::InfiniteOpt.InfiniteConstraintRef,
                                 key::Val{:TransData})

Map the shadow prices of the appropriate transcription constraints in the
transcription model to `icref`.
"""
function InfiniteOpt.map_shadow_price(icref::InfiniteOpt.InfiniteConstraintRef,
                                      key::Val{:TransData})
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(icref))
    return JuMP.shadow_price.(transcription_constraint(trans_model, icref))
end
