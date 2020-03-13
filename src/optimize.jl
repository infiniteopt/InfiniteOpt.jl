"""
    optimizer_model(model::InfiniteModel)::JuMP.Model

Return the JuMP model stored in `model` that is used to solve it.

**Example**
```julia-repl
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

Extend [`JuMP.bridge_constraints`](@ref JuMP.bridge_constraints(::JuMP.Model))
to return if an infinite model `model`
has an optimizer model where the optimizer is set and unsupported constraints
are automatically bridged to equivalent supported constraints when an
appropriate transformation is available.

**Example**
```julia-repl
julia> bridge_constraints(model)
false
```
"""
function JuMP.bridge_constraints(model::InfiniteModel)::Bool
    return JuMP.bridge_constraints(optimizer_model(model))
end

"""
    JuMP.add_bridge(model::InfiniteModel,
                    BridgeType::Type{<:MOI.Bridges.AbstractBridge})

Extend [`JuMP.add_bridge`](@ref JuMP.add_bridge(::JuMP.Model, ::Type{<:MOI.Bridges.AbstractBridge}))
to add `BridgeType` to the list of bridges that can
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
```julia-repl
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
```julia-repl
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
```julia-repl
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
                       [optimizer_constructor;
                       bridge_constraints::Bool = true])

Extend `JuMP.set_optimizer` to set optimizer of infinite models.
Specifically, the optimizer of the optimizer model is modified.

**Example**
```julia-repl
julia> set_optimizer(model, Clp.Optimizer)

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
                            optimizer_constructor;
                            bridge_constraints::Bool = true)
    JuMP.set_optimizer(optimizer_model(model), optimizer_constructor,
                       bridge_constraints = bridge_constraints)
    _set_optimizer_constructor(model, optimizer_constructor)
    return
end

"""
    JuMP.set_silent(model::InfiniteModel)

Extend [`JuMP.set_silent`](@ref JuMP.set_silent(::JuMP.Model)) for infinite
models to take precedence over any other
attribute controlling verbosity and requires the solver to produce no output.

**Example**
```julia-repl
julia> set_silent(model)
true
```
"""
function JuMP.set_silent(model::InfiniteModel)
    return JuMP.set_silent(optimizer_model(model))
end

"""
    JuMP.unset_silent(model::InfiniteModel)

Extend [`JuMP.unset_silent`](@ref JuMP.unset_silent(::JuMP.Model)) for infinite
models to neutralize the effect of the
`set_silent` function and let the solver attributes control the verbosity.

**Example**
```julia-repl
julia> unset_silent(model)
false
```
"""
function JuMP.unset_silent(model::InfiniteModel)
    return JuMP.unset_silent(optimizer_model(model))
end

"""
    JuMP.set_time_limit_sec(model::InfiniteModel, limit)

Extend [`set_time_limit_sec`](@ref JuMP.set_time_limit_sec(::JuMP.Model, ::Any))
to set the time limit (in seconds) of the solver.
Can be unset using `unset_time_limit_sec` or with `limit` set to `nothing`.

**Example**
```julia-repl
julia> set_time_limit_sec(model, 100)
100
```
"""
function JuMP.set_time_limit_sec(model::InfiniteModel, limit)
    return JuMP.set_time_limit_sec(optimizer_model(model), limit)
end

"""
    JuMP.unset_time_limit_sec(model::InfiniteModel)

Extend [`unset_time_limit_sec`](@ref JuMP.unset_time_limit_sec(::JuMP.Model)) to
unset the time limit of the solver. Can be set using `set_time_limit_sec`.

**Example**
```julia-repl
julia> unset_time_limit_sec(model)
```
"""
function JuMP.unset_time_limit_sec(model::InfiniteModel)
    return JuMP.unset_time_limit_sec(optimizer_model(model))
end

"""
    JuMP.time_limit_sec(model::InfiniteModel)

Extend [`time_limit_sec`](@ref JuMP.time_limit_sec(::JuMP.Model) to get the time
limit (in seconds) of the solve used by the optimizer model (`nothing` if unset).
Can be set using `set_time_limit_sec`.

**Example**
```julia-repl
julia> time_limit_sec(model)
100
```
"""
function JuMP.time_limit_sec(model::InfiniteModel)
    return JuMP.time_limit_sec(optimizer_model(model))
end

"""
    JuMP.set_optimizer_attribute(model::InfiniteModel, name::String, value)

Extend [`set_optimizer_attribute`](@ref JuMP.set_optimizer_attribute(::JuMP.Model, ::String, ::Any))
to specify a solver-specific attribute identified by `name` to `value`.

**Example**
```julia-repl
julia> set_optimizer_attribute(model, "SolverSpecificAttributeName", true)
true
```
"""
function JuMP.set_optimizer_attribute(model::InfiniteModel, name::String, value)
    return JuMP.set_optimizer_attribute(optimizer_model(model), name, value)
end

"""
    JuMP.set_optimizer_attribute(model::InfiniteModel,
                                 attr::MOI.AbstractOptimizerAttribute,
                                 value)

Extend [`set_optimizer_attribute`](@ref JuMP.set_optimizer_attribute(::JuMP.Model, ::MOI.AbstractOptimizerAttribute, ::Any))
to set the solver-specific attribute `attr` in `model` to `value`.

**Example**
```julia-repl
julia> set_optimizer_attribute(model, MOI.Silent(), true)
true
```
"""
function JuMP.set_optimizer_attribute(model::InfiniteModel,
                                      attr::MOI.AbstractOptimizerAttribute,
                                      value)
    return MOI.set(optimizer_model(model), attr, value)
end

"""
    JuMP.set_optimizer_attributes(model::InfiniteModel, pairs::Pair...)

Extend [`set_optimizer_attributes`](@ref JuMP.set_optimizer_attributes(::JuMP.Model, ::Pair))
to set multiple solver attributes given a list of `attribute => value` pairs.
Calls `set_optimizer_attribute(model, attribute, value)` for each pair.

**Example**
```julia-repl
julia> model = Model(Ipopt.Optimizer);

julia> set_optimizer_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
is equivalent to:
```julia-repl
julia> set_optimizer_attribute(model, "tol", 1e-4);

julia> set_optimizer_attribute(model, "max_iter", 100);
```
"""
function JuMP.set_optimizer_attributes(model::InfiniteModel, pairs::Pair...)
    for (name, value) in pairs
        JuMP.set_optimizer_attribute(model, name, value)
    end
    return
end

"""
    JuMP.get_optimizer_attribute(model::InfiniteModel, name::String)

Extend [`get_optimizer_attribute`](@ref JuMP.get_optimizer_attribute(::JuMP.Model, ::String))
to return the value associated with the solver-specific attribute named `name`.

**Example**
```julia-repl
julia> get_optimizer_attribute(model, "tol")
0.0001
````
"""
function JuMP.get_optimizer_attribute(model::InfiniteModel, name::String)
    return JuMP.get_optimizer_attribute(optimizer_model(model), name)
end

"""
    JuMP.get_optimizer_attribute(model::InfiniteModel,
                                 attr::MOI.AbstractOptimizerAttribute)

Extend [`get_optimizer_attribute`](@ref JuMP.get_optimizer_attribute(::JuMP.Model, ::MOI.AbstractOptimizerAttribute))
to return the value of the solver-specific attribute `attr` in `model`.

**Example**
```julia-repl
julia> get_optimizer_attribute(model, MOI.Silent())
true
````
"""
function JuMP.get_optimizer_attribute(model::InfiniteModel,
                                      attr::MOI.AbstractOptimizerAttribute)
    return MOI.get(optimizer_model(model), attr)
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
```julia-repl
julia> optimizer_model_key(model)
:TransData
```
"""
function optimizer_model_key(model::InfiniteModel)::Any
    key = collect(keys(optimizer_model(model).ext))
    if length(key) != 1
        error("Optimizer models should have 1 and only 1 extension key.")
    end
    return key[1]
end

"""
    build_optimizer_model!(model::InfiniteModel, key::Val{ext_key_name};
                           [kwargs...])

Build the optimizer model stored in `model` such that it can be
treated as a normal JuMP model, where the `Model.ext` field contains a key
that points to a datastructure that appropriately maps the data between the
two models. The key argument should be be typed to `Val{ext_key_name}`.
"""
 function build_optimizer_model! end

 """
     build_optimizer_model!(model::InfiniteModel; [kwargs...])

 Build the optimizer model stored in `model` such that it can be
 treated as a normal JuMP model. Specifically, translate the variables and
 constraints stored in `model` into ones that are stored in the optimizer model
 and can be solved. This is provided generally to accomodate extensions that use
 custom optimizer model types in accordance with [`optimizer_model_key`](@ref).
 However, it may be useful in certain applications when the user desires to
 force a build without calling `optimize!`.
 Extensions will need to implement their own version of the function
 `build_optimizer_model!(model::InfiniteModel, key::Val{ext_key_name}; kwargs...)`.

**Example**
```julia-repl
julia> build_optimizer_model!(model)

julia> optimizer_model_ready(model)
true
```
"""
function build_optimizer_model!(model::InfiniteModel; kwargs...)
  key = optimizer_model_key(model)
  build_optimizer_model!(model, Val(key); kwargs...)
  return
end

"""
    optimizer_model_variable(vref::InfOptVariableRef, key::Val{ext_key_name};
                             [kwargs...])

Return the reformulation variable(s) stored in the optimizer model that correspond
to `vref`. This needs to be defined for extensions that implement a custom
optimizer model type. Principally, this is accomplished by typed the `key`
argument to `Val{ext_key_name}`. Keyword arguments can be added as needed.
"""
function optimizer_model_variable end

# Fallback for unextended keys
function optimizer_model_variable(vref::InfOptVariableRef, key; kwargs...)
    error("`optimizer_model_variable` not implemented for optimizer model
          key `$key`.")
end

"""
    optimizer_model_variable(vref::InfOptVariableRef; [kwargs...])

Return the reformulation variable(s) stored in the optimizer model that correspond
to `vref`. By default, no keyword arguments `kwargs` are employed by
`TranscriptionOpt`, but extensions may employ `kwargs` in accordance with
their implementation of [`optimizer_model_variable`](@ref). Errors if such an
extension has not been written. Also errors if no such variable can be found in
the optimizer model.

**Example**
```julia-repl
julia> optimizer_model_variable(x) # infinite variable
2-element Array{VariableRef,1}:
 x(support: 1)
 x(support: 2)

julia> optimizer_model_variable(z) # hold variable
z
```
"""
function optimizer_model_variable(vref::InfOptVariableRef; kwargs...)
    key = optimizer_model_key(JuMP.owner_model(vref))
    return optimizer_model_variable(vref, Val(key); kwargs...)
end

"""
    variable_supports(optimizer_model::JuMP.Model, vref::InfiniteVariableRef,
                      key::Val{ext_key_name}; [kwargs...])::Vector

Return the supports associated with the mappings of `vref` in `optimizer_model`.
This dispatches off of `key` which permits optimizer model extensions. This
should throw an error if `vref` is not associated with the variable mappings
stored in `optimizer_model`. Keyword arguments can be added as needed.
"""
function variable_supports end

# fallback for unextended keys
function variable_supports(optimizer_model::JuMP.Model, vref::InfiniteVariableRef,
                           key; kwargs...)
  error("`variable_supports` not implemented for optimizer model key `$key`.")
end

"""
    supports(vref::InfiniteVariableRef; [kwargs...])::Vector

Return the supports associated with `vref` in the optimizer
model. Errors if [`variable_supports`](@ref) has not been extended for the
optimizer model type or if `vref` is not be reformulated in the optimizer model.
By default, keyword arugments are not used, but may employed by extensions.

**Example**
```julia-repl
julia> supports(vref)
Dict{Int64,Tuple{Float64}} with 2 entries:
  2 => (1.0,)
  1 => (0.0,)
```
"""
function supports(vref::InfiniteVariableRef; kwargs...)::Vector
    model = optimizer_model(JuMP.owner_model(vref))
    key = optimizer_model_key(JuMP.owner_model(vref))
    return variable_supports(model, vref, Val(key); kwargs...)
end

"""
    optimizer_model_constraint(cref::GeneralConstraintRef,
                               key::Val{ext_key_name}; [kwargs...])

Return the reformulation constraint(s) stored in the optimizer model that correspond
to `cref`. This needs to be defined for extensions that implement a custom
optimizer model type. Principally, this is accomplished by typed the `key`
argument to `Val{ext_key_name}`. Keyword arguments can be added as needed.
"""
function optimizer_model_constraint end

# Fallback for unextended keys
function optimizer_model_constraint(cref::GeneralConstraintRef, key; kwargs...)
    error("`optimizer_model_constraint` not implemented for optimizer model
          key `$key`.")
end

"""
    optimizer_model_constraint(cref::GeneralConstraintRef; [kwargs...])

Return the reformulation constraint(s) stored in the optimizer model that correspond
to `cref`. By default, no keyword arguments `kwargs` are employed by
`TranscriptionOpt`, but extensions may employ `kwargs` in accordance with
their implementation of [`optimizer_model_constraint`](@ref). Errors if such an
extension has not been written. Also errors if no such constraint can be found in
the optimizer model.

**Example**
```julia-repl
julia> optimizer_model_constraint(c1) # finite constraint
c1 : x(support: 1) - y <= 3.0
```
"""
function optimizer_model_constraint(cref::GeneralConstraintRef; kwargs...)
    key = optimizer_model_key(JuMP.owner_model(cref))
    return optimizer_model_constraint(cref, Val(key); kwargs...)
end

"""
    constraint_supports(optimizer_model::JuMP.Model, cref::GeneralConstraintRef,
                        key::Val{ext_key_name}; [kwargs...])::Vector

Return the supports associated with the mappings of `cref` in `optimizer_model`.
This dispatches off of `key` which permits optimizer model extensions. This
should throw an error if `cref` is not associated with the variable mappings
stored in `optimizer_model`. Keyword arguments can be added as needed.
"""
function constraint_supports end

# fallback for unextended keys
function constraint_supports(optimizer_model::JuMP.Model,
                             cref::GeneralConstraintRef,
                             key; kwargs...)
  error("`constraint_supports` not implemented for optimizer model key `$key` " *
        "and/or constraint type `$(typeof(cref))`.")
end

"""
    supports(cref::GeneralConstraintRef; [kwargs...])::Vector

Return the support associated with `cref`. Errors if `cref` is
not associated with the constraint mappings stored in `optimizer_model` or if
[`constraint_supports`](@ref) has not been extended. By default, no keyword
arguments are accepted, but extensions may employ some.

**Example**
```julia-repl
julia> supports(cref)
Dict{Int64,Tuple{Float64}} with 2 entries:
  2 => (1.0,)
  1 => (0.0,)
```
"""
function supports(cref::GeneralConstraintRef; kwargs...)::Vector
    model = optimizer_model(JuMP.owner_model(cref))
    key = optimizer_model_key(JuMP.owner_model(cref))
    return constraint_supports(model, cref, Val(key); kwargs...)
end

"""
    constraint_parameter_refs(optimizer_model::JuMP.Model,
                              cref::GeneralConstraintRef,
                              key::Val{ext_key_name}; [kwargs...])::Tuple

Return the infinite parameter references associated with the mappings of `cref`
in `optimizer_model`. This dispatches off of `key` which permits optimizer model
extensions. This should throw an error if `cref` is not associated with the
variable mappings stored in `optimizer_model`. Keyword arguments can be added
as needed.
"""
function constraint_parameter_refs end

# fallback for unextended keys
function constraint_parameter_refs(optimizer_model::JuMP.Model,
                                   cref::GeneralConstraintRef,
                                   key; kwargs...)
  error("`constraint_parameter_refs` not implemented for optimizer model key `$key` " *
        "and/or constraint type `$(typeof(cref))`.")
end

"""
    parameter_refs(cref::GeneralConstraintRef; [kwargs...])::Tuple

Return the infinite parameters associated with `cref`. Errors if `cref` is
not associated with the constraint mappings stored in `optimizer_model` or if
[`constraint_parameter_refs`](@ref) has not been extended. By default, no keyword
arguments are accepted, but extensions may employ some.

**Example**
```julia-repl
julia> parameter_refs(cref)
(t, x)
```
"""
function parameter_refs(cref::GeneralConstraintRef; kwargs...)::Tuple
    model = optimizer_model(JuMP.owner_model(cref))
    key = optimizer_model_key(JuMP.owner_model(cref))
    return constraint_parameter_refs(model, cref, Val(key); kwargs...)
end

"""
    JuMP.solver_name(model::InfiniteModel)

Extend [`solver_name`](@ref JuMP.solver_name(::JuMP.Model)) to return the name
of the solver being used if there is an optimizer selected and it has a name
attribute. Otherwise, an error is thrown.

**Example**
```julia-repl
julia> solver_name(model)
"Gurobi"
```
"""
function JuMP.solver_name(model::InfiniteModel)
    return JuMP.solver_name(optimizer_model(model))
end

"""
    JuMP.backend(model::InfiniteModel)

Extend [`backend`](@ref JuMP.backend(::JuMP.Model)) to return the
`MathOptInterface` backend associated with the optimizer model. Note this will
be empty if the optimizer model has not been build yet.

**Example**
```julia-repl
julia> moi_model = backend(model);
```
"""
function JuMP.backend(model::InfiniteModel)
    return JuMP.backend(optimizer_model(model))
end

"""
    JuMP.mode(model::InfiniteModel)

Extend [`mode`](@ref JuMP.mode(::JuMP.Model)) to return the `MathOptInterface`
mode the optimizer model is in.

**Example**
```julia-repl
julia> mode(model)
AUTOMATIC::ModelMode = 0
```
"""
function JuMP.mode(model::InfiniteModel)
    return JuMP.mode(optimizer_model(model))
end

"""
    JuMP.optimize!(model::InfiniteModel;
                   bridge_constraints::Bool=true, kwargs...])

Extend [`JuMP.optimize!`](@ref JuMP.optimize!(::JuMP.Model, ::Any))
to optimize infinite models using the internal
optimizer model. Will call [`build_optimizer_model!`](@ref) if the optimizer
model isn't up to date. The `kwargs` correspond to keyword arguments passed to
[`build_optimizer_model!`](@ref) if any are defined.

**Example**
```julia-repl
julia> optimize!(model)

julia> has_values(model)
true
```
"""
function JuMP.optimize!(model::InfiniteModel;
                        bridge_constraints::Bool = true,
                        kwargs...)
    if !optimizer_model_ready(model)
        build_optimizer_model!(model; kwargs...)
    end
    JuMP.optimize!(optimizer_model(model))
    return
end

"""
    JuMP.result_count(model::InfiniteModel)

Extend [`result_count`](@ref JuMP.result_count(::JuMP.Model)) to return the
number of results available to query after a call to `optimize!`.

**Example**
```julia-repla
julia> result_count(model)
1
```
"""
function JuMP.result_count(model::InfiniteModel)::Int
    return MOI.get(optimizer_model(model), MOI.ResultCount())
end
