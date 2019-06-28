"""
    optimize_model(model::InfiniteModel)
Return the JuMP model stored in `model` that is used to solve it.
"""
optimize_model(model::InfiniteModel) = model.optimize_model

"""
    transcription_model(model::InfiniteModel)
Return the transcription model stored in `model` if that is what is stored in
`model.optimize_model`.
"""
function transcription_model(model::InfiniteModel)
    !is_transcription_model(optimize_model(model)) && error("The model does not contain a transcription model.")
    return optimize_model(model)
end

"""
    set_optimize_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
Set the JuMP model that is used to solve `inf_model`.
"""
function set_optimize_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
    inf_model.optimize_model = opt_model
    return
end

"""
    transcribe!(model::InfiniteModel)
Transcribe `model` and store it as a `TranscriptionModel` in the
`model.optimize_model` field which can be accessed with `transcription_model`.
"""
function transcribe!(model::InfiniteModel)
    trans_model = TranscriptionModel(model, model.kwargs...)
    if !isa(model.optimizer_factory, Nothing)
        JuMP.set_optimizer(trans_model, model.optimizer_factory)
    end
    return
end

"""
    JuMP.set_optimizer(model::InfiniteModel, optimizer_factory::JuMP.OptimizerFactory;
                       bridge_constraints::Bool=true)
Extend `JuMP.set_optimizer` to accomodate infinite models.
"""
function JuMP.set_optimizer(model::InfiniteModel, optimizer_factory::JuMP.OptimizerFactory;
                            bridge_constraints::Bool=true)
    JuMP.set_optimizer(optimize_model(model), optimizer_factory,
                       bridge_constraints = bridge_constraints)
    return
end

# Extend the solve error function
JuMP.solve(model::InfiniteModel) = JuMP.solve(optimize_model(model))

# TODO Make more modular to accomodate extensions
"""
    JuMP.optimize!(model::InfiniteModel,
                   optimizer_factory::Union{Nothing, OptimizerFactory}=nothing;
                   bridge_constraints::Bool=true,
                   ignore_optimize_hook=(model.optimize_hook === nothing),
                   kwargs...)
Extend `JuMP.optimize!` to accomodate infinite models.
"""
function JuMP.optimize!(model::InfiniteModel,
                        optimizer_factory::Union{Nothing, OptimizerFactory}=nothing;
                        bridge_constraints::Bool=true,
                        ignore_optimize_hook=(model.optimize_hook === nothing),
                        kwargs...)
    transcribe!(model)
    JuMP.optimize!(optimize_model(model), optimizer_factory,
                   bridge_constraints = bridge_constraints,
                   ignore_optimize_hook = ignore_optimize_hook,
                   kwargs...)
    return
end
