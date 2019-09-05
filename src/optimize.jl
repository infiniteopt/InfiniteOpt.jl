"""
    optimizer_model(model::InfiniteModel)::JuMP.Model

Return the JuMP model stored in `model` that is used to solve it.

**Example**
```julia
julia> opt_model = optimizer_model(model)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
optimizer_model(model::InfiniteModel)::JuMP.Model = model.optimizer_model

"""
    JuMP.bridge_constraints(model::InfiniteModel)::Bool

Extend [`JuMP.bridge_constraints`](@ref) to return if an infinite model `model`
has an optimizer model where the optimizer is set and unsupported constraints
are automatically bridged to equivalent supported constraints when an
appropriate transformation is available.

**Example**
```julia
julia> bridge_constraints(model)
false
```
"""
function JuMP.bridge_constraints(model::InfiniteModel)::Bool
    return JuMP.bridge_constraints(optimizer_model(model))
end

"""
    JuMP.add_bridge(model::Model, BridgeType::Type{<:MOI.Bridges.AbstractBridge})

Extend [`JuMP.add_bridge`] to add `BridgeType` to the list of bridges that can
be used by the optimizer model to transform unsupported constraints into an
equivalent formulation using only constraints supported by the optimizer.
"""
function JuMP.add_bridge(model::InfiniteModel,
                    BridgeType::Type{<:MOI.Bridges.AbstractBridge})
    JuMP.add_bridge(optimizer_model(model), BridgeType)
    return
end

"""
    optimizer_model_ready(model::InfiniteModel)::Bool

Return `Bool` if the optimizer model is up to date with `model`.

**Example**
```julia
julia> optimizer_model_ready(model)
false
```
"""
optimizer_model_ready(model::InfiniteModel)::Bool = model.ready_to_optimize

"""
    set_optimizer_model_ready(model::InfiniteModel, status::Bool)

Set the status of the optimizer model to whether it is up to date or not. Note
is more intended as an internal function, but is useful for extensions.

**Example**
```julia
julia> set_optimizer_model_ready(model, true)

julia> optimizer_model_ready(model)
true
```
"""
function set_optimizer_model_ready(model::InfiniteModel, status::Bool)
     model.ready_to_optimize = status
     return
end

"""
    set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model)

Specify the JuMP model that is used to solve `inf_model`. This is intended for
internal use and extensions. Note that `opt_model` should contain extension
data to allow it to map to `inf_model` in a manner similar to
[`TranscriptionModel`](@ref).

**Example**
```julia
julia> set_optimizer_model(model, TranscriptionModel())

julia> optimizer_model(model)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
function set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model)
    inf_model.optimizer_model = opt_model
    set_optimizer_model_ready(inf_model, false)
    return
end

"""
    JuMP.set_optimizer(model::InfiniteModel,
                       optimizer_factory::JuMP.OptimizerFactory;
                       bridge_constraints::Bool=true)

Extend [`JuMP.set_optimizer`](@ref) to set optimizer of infinite models.
Specifically, the optimizer of the optimizer model is modified.

**Example**
```julia
julia> set_optimizer(model, with_optimizer(Clp.Optimizer))

julia> optimizer_model(model)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: SolverName() attribute not implemented by the optimizer.
```
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

Extend `JuMP.set_silent` for infinite models to take precedence over any other
attribute controlling verbosity and requires the solver to produce no output.
"""
function JuMP.set_silent(model::InfiniteModel)
    return JuMP.set_silent(optimizer_model(model))
end

"""
    JuMP.unset_silent(model::InfiniteModel)

Extend `JuMP.unset_silent` for infinite models to Neutralize the effect of the
`set_silent` function and let the solver attributes control the verbosity.
"""
function JuMP.unset_silent(model::InfiniteModel)
    return JuMP.unset_silent(optimizer_model(model))
end

"""
    JuMP.set_parameter(model::Model, name, value)

Sets solver-specific parameter identified by `name` to `value`.
"""
function JuMP.set_parameter(model::InfiniteModel, name, value)
    return JuMP.set_parameter(optimizer_model(model), name, value)
end

# Extend the solve error function
JuMP.solve(model::InfiniteModel) = JuMP.solve(optimizer_model(model))

"""
    optimizer_model_key(model::InfiniteModel)::Any

Return the extension key used in the optimizer model of `model`. Errors if
`optimizer_model.ext` contains more than one key. This is intended for internal
use and extensions. For extensions this is used to dispatch to the appropriate
optmizer model functions such as extensions to [`build_optimizer_model!`](@ref).

**Example**
```julia
julia> optimizer_model_key(model)
:TransData
```
"""
function optimizer_model_key(model::InfiniteModel)::Any
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
that points to a datastructure that appropriately maps the data between the
two models. The key argument should be be typed to `Val{ext_key_name}`.
"""
 function build_optimizer_model! end

 """
     build_optimizer_model!(model::InfiniteModel)

 Build the optimizer model stored in `model` such that it can be
 treated as a normal JuMP model. Specifically, translate the variables and
 constraints stored in `model` into ones that are stored in the optimizer model
 and can be solved. This is build generally to accomodate extensions that use
 custom optimizer model types in accordance with [`optimizer_model_key`](@ref).
 Extensions will need to implement their own version of the function
 `build_optimizer_model!(model::InfiniteModel, key::Val{ext_key_name})`.

 **Example**
 ```julia
julia> build_optimizer_model!(model)

julia> optimizer_model_ready(model)
true
 ```
 """
  function build_optimizer_model!(model::InfiniteModel)
      key = optimizer_model_key(model)
      build_optimizer_model!(model, Val(key))
      return
  end

"""
    JuMP.optimize!(model::InfiniteModel,
                   optimizer_factory::Union{Nothing, OptimizerFactory} = nothing;
                   bridge_constraints::Bool=true, kwargs...)

Extend [`JuMP.optimize!`](@ref) to optimize infinite models using the internal
optimizer model. Will call [`build_optimizer_model!`](@ref) if the optimizer
model isn't up to date. The `kwargs` correspond to keyword arguments passed to
[`build_optimizer_model!`](@ref) if any are defined.

**Example**
```julia
julia> optimize!(model, with_optimizer(Clp.Optimizer))

julia> has_values(model)
true
```
"""
function JuMP.optimize!(model::InfiniteModel,
                        optimizer_factory::Union{Nothing,
                                               JuMP.OptimizerFactory} = nothing;
                        bridge_constraints::Bool = true,
                        kwargs...)
    key = optimizer_model_key(model)
    if !optimizer_model_ready(model)
        build_optimizer_model!(model, Val(key), kwargs...)
    end
    JuMP.optimize!(optimizer_model(model), optimizer_factory,
                   bridge_constraints = bridge_constraints)
    return
end
