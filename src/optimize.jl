"""
    optimizer_model(model::InfiniteModel)
Return the JuMP model stored in `model` that is used to solve it.
"""
optimizer_model(model::InfiniteModel) = model.optimizer_model

"""
    JuMP.bridge_constraints(model::InfiniteModel)
Extend `JuMP.bridge_constraints` to accomodate infinite models.
"""
function JuMP.bridge_constraints(model::InfiniteModel)
    return JuMP.bridge_constraints(optimizer_model(model))
end

"""
    JuMP.add_bridge(model::Model, BridgeType::Type{<:MOI.Bridges.AbstractBridge})
Extend `JuMP.add_bridge` to accomodate infinite models.
"""
function JuMP.add_bridge(model::InfiniteModel,
                    BridgeType::Type{<:MOI.Bridges.AbstractBridge})
    return JuMP.add_bridge(optimizer_model(model), BridgeType)
end

"""
    optimizer_model_status(model::InfiniteModel)
Return true if the optimizer model is up to date with `model` or false otherwise.
"""
optimizer_model_status(model::InfiniteModel) = model.ready_to_optimize

"""
    set_optimizer_model_status(model::InfiniteModel, status::Bool)
Set the status of the optimizer model to whether it is up to date or not.
"""
function set_optimizer_model_status(model::InfiniteModel, status::Bool)
     model.ready_to_optimize = status
     return
end

"""
    set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
Set the JuMP model that is used to solve `inf_model`.
"""
function set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
    inf_model.optimizer_model = opt_model
    set_optimizer_model_status(inf_model, false)
    return
end

"""
    JuMP.set_optimizer(model::InfiniteModel, optimizer_factory::JuMP.OptimizerFactory;
                       bridge_constraints::Bool=true)
Extend `JuMP.set_optimizer` to accomodate infinite models.
"""
function JuMP.set_optimizer(model::InfiniteModel,
                            optimizer_factory::JuMP.OptimizerFactory;
                            bridge_constraints::Bool=true)
    JuMP.set_optimizer(optimizer_model(model), optimizer_factory,
                       bridge_constraints = bridge_constraints)
    return
end

"""
    JuMP.set_silent(model::InfiniteModel)
Takes precedence over any other attribute controlling verbosity
and requires the solver to produce no output.
"""
# TODO include with next JuMP version
# function JuMP.set_silent(model::InfiniteModel)
#     return MOI.set(optimizer_model(model), MOI.Silent(), true)
# end

"""
    JuMP.unset_silent(model::InfiniteModel)
Neutralize the effect of the `set_silent` function and let the solver
attributes control the verbosity.
"""
# TODO include with next JuMP version
# function JuMP.unset_silent(model::InfiniteModel)
#     return MOI.set(optimizer_model(model), MOI.Silent(), false)
# end

# Extend the solve error function
JuMP.solve(model::InfiniteModel) = JuMP.solve(optimizer_model(model))

"""
    set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
Set the JuMP model that is used to solve `inf_model`.
"""
function optimizer_model_key(model::InfiniteModel)
    key = collect(keys(optimizer_model(model).ext))
    if length(key) != 1
        error("Optimizer models should only have 1 extension key.")
    end
    return key[1]
end

"""
    build_optimizer_model!(model::InfiniteModel, key; kwargs...)
Build the optimizer model stored in `model` such that it can be
treated as a normal JuMP model, where the `Model.ext` field contains a key
that points to a datastructure that appropriately mapps the data between the
two models. The key argument should be be typed to `Val{ext_key_name}`.
"""
 function build_optimizer_model! end

 """
     build_optimizer_model!(model::InfiniteModel)
 Build the optimizer model stored in `model` such that it can be
 treated as a normal JuMP model.
 """
  function build_optimizer_model!(model::InfiniteModel)
      key = optimizer_model_key(model)
      return build_optimizer_model!(model, Val(key))
  end

"""
    JuMP.optimize!(model::InfiniteModel,
                   optimizer_factory::Union{Nothing, OptimizerFactory}=nothing;
                   bridge_constraints::Bool=true,
                   ignore_optimize_hook=(model.optimize_hook === nothing),
                   kwargs...)
Extend `JuMP.optimize!` to accomodate infinite models.
"""
function JuMP.optimize!(model::InfiniteModel,
                        optimizer_factory::Union{Nothing, JuMP.OptimizerFactory} = nothing;
                        bridge_constraints::Bool = true,
                        kwargs...)
    key = optimizer_model_key(model)
    if !optimizer_model_status(model)
        build_optimizer_model!(model, Val(key), kwargs...)
    end
    JuMP.optimize!(optimizer_model(model), optimizer_factory,
                   bridge_constraints = bridge_constraints)
    return
end
