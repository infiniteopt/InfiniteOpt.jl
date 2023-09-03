################################################################################
#                              OPTIMIZER MODEL BASICS
################################################################################
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

Extend `JuMP.bridge_constraints` to return if an infinite model `model`
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

Extend `JuMP.add_bridge` to add `BridgeType` to the list of bridges that can
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
backend_ready(model::InfiniteModel) = model.ready_to_optimize

"""
    set_backend_ready(model::InfiniteModel, status::Bool)

Set the status of the optimizer model to whether it is up to date or not. Note
is more intended as an internal function, but is useful for extensions.

**Example**
```julia-repl
julia> set_backend_ready(model, true)

julia> backend_ready(model)
true
```
"""
function set_backend_ready(model::InfiniteModel, status::Bool)
    model.backend_ready = status
    return
end

"""
    add_infinite_model_optimizer(opt_model::JuMP.Model, inf_model::InfiniteModel)

Parse the current optimizer and its attributes associated with `inf_model` and load
them into `opt_model`. This is intended to be used as an internal method
for [`set_optimizer_model`](@ref).
"""
function add_infinite_model_optimizer(opt_model::JuMP.Model,
                                      inf_model::InfiniteModel)
    if !isa(inf_model.optimizer_constructor, Nothing)
        bridge_constrs = JuMP.bridge_constraints(inf_model)
        JuMP.set_optimizer(opt_model, inf_model.optimizer_constructor,
                           add_bridges =  bridge_constrs)
    end
    # parse the attributes (this is a hacky workaround)
    for (attr, val) in JuMP.backend(inf_model).model_cache.optattr
        MOI.set(opt_model, attr, val)
    end
    return
end

"""
    set_optimizer_model(inf_model::InfiniteModel, opt_model::JuMP.Model;
                        inherit_optimizer::Bool = true)

Specify the JuMP model that is used to solve `inf_model`. This is intended for
internal use and extensions. Note that `opt_model` should contain extension
data to allow it to map to `inf_model` in a manner similar to
[`TranscriptionModel`](@ref). `inherit_optimizer` indicates whether
[`add_infinite_model_optimizer`](@ref) should be invoked on the new optimizer
mode to inherit the optimizer constuctor and attributes currently stored in
`inf_model`.

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
function set_optimizer_model(
    inf_model::InfiniteModel, 
    opt_model::JuMP.Model;
    inherit_optimizer::Bool = true
    )
    if inherit_optimizer
        add_infinite_model_optimizer(opt_model, inf_model)
    end
    inf_model.optimizer_model = opt_model
    set_optimizer_model_ready(inf_model, false)
    return
end

"""
    optimizer_model_key(model::JuMP.Model)::Any

Return the extension key used in the optimizer model `model`. Errors if
`model.ext` contains more than one key. This is intended for internal
use and extensions. For extensions this is used to dispatch to the appropriate
optmizer model functions such as extensions to [`build_optimizer_model!`](@ref).
This is intended as an internal method. See 
[`optimizer_model_key`](@ref optimizer_model_key(::InfiniteModel)) 
for the public method
"""
function optimizer_model_key(model::JuMP.Model)
    if length(model.ext) != 1
        error("Optimizer models should have 1 and only 1 extension key of the " *
              "form `Model.ext[:my_ext_key] = MyExtData`.")
    end
    return first(keys(model.ext))
end

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
    return optimizer_model_key(optimizer_model(model))
end

################################################################################
#                         OPTIMIZER METHOD EXTENSIONS
################################################################################
"""
    JuMP.set_optimizer(model::InfiniteModel,
                       [optimizer_constructor;
                       add_bridges::Bool = true])

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
function JuMP.set_optimizer(
    model::InfiniteModel,
    optimizer_constructor;
    add_bridges::Bool = true
    )
    JuMP.set_optimizer(optimizer_model(model), optimizer_constructor,
                       add_bridges = add_bridges)
    _set_optimizer_constructor(model, optimizer_constructor)
    return
end

"""
    JuMP.set_silent(model::InfiniteModel)

Extend `JuMP.set_silent` for infinite models to take precedence over any other
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

Extend `JuMP.unset_silent` for infinite models to neutralize the effect of the
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

Extend `set_time_limit_sec` to set the time limit (in seconds) of the solver.
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

Extend `unset_time_limit_sec` to unset the time limit of the solver. Can be set 
using `set_time_limit_sec`.

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

Extend `time_limit_sec` to get the time limit (in seconds) of the solve used by 
the optimizer model (`nothing` if unset). Can be set using `set_time_limit_sec`.

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

Extend `set_optimizer_attribute` to specify a solver-specific attribute 
identified by `name` to `value`.

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

Extend `set_optimizer_attribute` to set the solver-specific attribute `attr` in 
`model` to `value`.

**Example**
```julia-repl
julia> set_optimizer_attribute(model, MOI.Silent(), true)
true
```
"""
function JuMP.set_optimizer_attribute(
    model::InfiniteModel,
    attr::MOI.AbstractOptimizerAttribute,
    value
    )
    return MOI.set(optimizer_model(model), attr, value)
end

"""
    JuMP.set_optimizer_attributes(model::InfiniteModel, pairs::Pair...)

Extend `set_optimizer_attributes` to set multiple solver attributes given a 
list of `attribute => value` pairs. Calls 
`set_optimizer_attribute(model, attribute, value)` for each pair.

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

Extend `get_optimizer_attribute` to return the value associated with the 
solver-specific attribute named `name`.

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

Extend `get_optimizer_attribute` to return the value of the solver-specific 
attribute `attr` in `model`.

**Example**
```julia-repl
julia> get_optimizer_attribute(model, MOI.Silent())
true
````
"""
function JuMP.get_optimizer_attribute(
    model::InfiniteModel,
    attr::MOI.AbstractOptimizerAttribute
    )
    return MOI.Base.get(optimizer_model(model), attr)
end

"""
    JuMP.solver_name(model::InfiniteModel)

Extend `solver_name` to return the name of the solver being used if there is an 
optimizer selected and it has a name attribute. Otherwise, an error is thrown.

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

Extend `backend` to return the `MathOptInterface` backend associated with the 
optimizer model. Note this will be empty if the optimizer model has not been 
build yet.

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

Extend `mode` to return the `MathOptInterface` mode the optimizer model is in.

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
    JuMP.result_count(model::InfiniteModel)

Extend `result_count` to return the number of results available to query after a 
call to `optimize!`.

**Example**
```julia-repla
julia> result_count(model)
1
```
"""
function JuMP.result_count(model::InfiniteModel)::Int
    return MOI.Base.get(optimizer_model(model), MOI.ResultCount())
end

################################################################################
#                       OPTIMIZER MODEL BUILD METHODS
################################################################################
"""
    build_optimizer_model!(model::InfiniteModel, key::Val{ext_key_name};
                           [kwargs...])

Build the optimizer model stored in `model` such that it can be
treated as a normal JuMP model, where the `Model.ext` field contains a key
that points to a datastructure that appropriately maps the data between the
two models. The key argument should be be typed to `Val{ext_key_name}`. This
should also use [`clear_optimizer_model_build!`](@ref) to empty the out the current
optimizer model. Ultimately, [`set_optimizer_model`](@ref) should be called
to insert the build optimizer model into `model` and [`set_optimizer_model_ready`](@ref)
should be used to update the optimizer model's status.
"""
function build_optimizer_model! end

"""
    clear_optimizer_model_build!(model::JuMP.Model)::JuMP.Model

Empty the optimizer model using appropriate calls of `Base.empty!`. This
effectively resets `model` except the optimizer, its attributes, and an an emptied
optimizer model data struct are maintained. This is intended as an internal
method for use by [`build_optimizer_model!`](@ref).
"""
function clear_optimizer_model_build!(model::JuMP.Model)
    key = optimizer_model_key(model)
    data_type = typeof(model.ext[key])
    empty!(model)
    model.ext[key] = data_type()
    model.operator_counter = 0
    return model
end

"""
    clear_optimizer_model_build!(model::InfiniteModel)::JuMP.Model

Empty the optimizer model using appropriate calls of `Base.empty!`. This
effectively resets `model.optimizer_model` except the optimizer, its attributes,
and an an emptied optimizer model data struct are maintained. This is intended
as an internal method for use by [`build_optimizer_model!`](@ref).
"""
function clear_optimizer_model_build!(model::InfiniteModel)::JuMP.Model
    return clear_optimizer_model_build!(optimizer_model(model))
end

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
    if num_parameters(model, InfiniteParameter) == 0
        @warn("Finite models (i.e., `InfiniteModel`s with no infinite " * 
              "parameters) should be modeled directly via a `Model` in JuMP.jl.")
    end
    key = optimizer_model_key(model)
    build_optimizer_model!(model, Val(key); kwargs...)
    return
end

################################################################################
#                  OPTIMIZER MODEL MAPPING METHODS (VARIABLES)
################################################################################
"""
    optimizer_model_variable(vref::GeneralVariableRef, key::Val{ext_key_name};
                             [kwargs...])

Return the reformulation variable(s) stored in the optimizer model that correspond
to `vref`. This needs to be defined for extensions that implement a custom
optimizer model type. Principally, this is accomplished by typed the `key`
argument to `Val{ext_key_name}`. Keyword arguments can be added as needed.
"""
function optimizer_model_variable end

# Fallback for unextended keys
function optimizer_model_variable(vref::GeneralVariableRef, key; kwargs...)
    error("`optimizer_model_variable` not implemented for optimizer model " *
          "key `$(typeof(key).parameters[1])`.")
end

"""
    optimizer_model_variable(vref::GeneralVariableRef; 
                             [label::Type{<:AbstractSupportLabel} = PublicLabel, 
                             ndarray::Bool = false,
                             kwargs...])

Return the reformulation variable(s) stored in the optimizer model that correspond
to `vref`. Also errors if no such variable can be found in
the optimizer model.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of [`optimizer_model_variable`](@ref). Errors if such an
extension has not been written. 

By default only the variables associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, infinite variables are 
returned as a list corresponding to their supports. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the variable has multiple 
infinite parameter dependencies. The corresponding supports are obtained via 
`supports` using the same keyword arguments.

**Example**
```julia-repl
julia> optimizer_model_variable(x) # infinite variable
2-element Array{VariableRef,1}:
 x(support: 1)
 x(support: 2)

julia> optimizer_model_variable(z) # finite variable
z
```
"""
function optimizer_model_variable(vref::GeneralVariableRef; kwargs...)
    key = optimizer_model_key(JuMP.owner_model(vref))
    return optimizer_model_variable(vref, Val(key); kwargs...)
end

"""
    variable_supports(optimizer_model::JuMP.Model, vref,
                      key::Val{ext_key_name}; 
                      [kwargs...])::Vector

Return the supports associated with the mappings of `vref` in `optimizer_model`.
This dispatches off of `key` which permits optimizer model extensions. This
should throw an error if `vref` is not associated with the variable mappings
stored in `optimizer_model`. Keyword arguments can be added as needed. Note that
no extension is necessary for point or finite variables. 
"""
function variable_supports end

# fallback for unextended keys
function variable_supports(optimizer_model::JuMP.Model, vref, key; kwargs...)
    error("`variable_supports` not implemented for optimizer model key " *
          "`$(typeof(key).parameters[1])` and/or variable type $(typeof(vref)).")
end

# FiniteRef
function variable_supports(optimizer_model::JuMP.Model, vref::FiniteRef,
                           key; kwargs...)
    return ()
end

"""
    supports(vref::DecisionVariableRef; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel, 
             ndarray::Bool = false,
             kwargs...])

Return the supports associated with `vref` in the optimizer
model. Errors if [`InfiniteOpt.variable_supports`](@ref) has not been extended for the
optimizer model type or if `vref` is not be reformulated in the optimizer model.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `variable_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the supports of infinite 
variables are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the variable has multiple 
infinite parameter dependencies.

**Example**
```julia-repl
julia> supports(vref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(
    vref::Union{DecisionVariableRef, MeasureRef, ParameterFunctionRef}; 
    kwargs...
    )
    model = optimizer_model(JuMP.owner_model(vref))
    key = optimizer_model_key(JuMP.owner_model(vref))
    return variable_supports(model, vref, Val(key); kwargs...)
end

################################################################################
#                 OPTIMIZER MODEL MAPPING METHODS (EXPRESSIONS)
################################################################################
"""
    optimizer_model_expression(expr, key::Val{ext_key_name}; [kwargs...])

Return the reformulation expression(s) stored in the optimizer model that correspond
to `expr`. This needs to be defined for extensions that implement a custom
optimizer model type. Principally, this is accomplished by typed the `key`
argument to `Val{ext_key_name}`. Keyword arguments can be added as needed.
Note that if `expr` is a `GeneralVariableRef` this just dispatches to
`optimizer_model_variable`.
"""
function optimizer_model_expression end

# Fallback for unextended keys
function optimizer_model_expression(expr, key; kwargs...)
    error("`optimizer_model_expression` not defined for optimizer model " *
          "key `$(typeof(key).parameters[1])` and expression type " *
          "`$(typeof(expr))`.")
end

# Define for variable reference expressions
function optimizer_model_expression(expr::GeneralVariableRef, key; kwargs...)
    return optimizer_model_variable(expr, key; kwargs...)
end

"""
    optimizer_model_expression(expr::JuMP.AbstractJuMPScalar; 
                               [label::Type{<:AbstractSupportLabel} = PublicLabel,
                               ndarray::Bool = false, 
                               kwargs...])

Return the reformulation expression(s) stored in the optimizer model that correspond
to `expr`. Also errors if no such expression can be found in
the optimizer model (meaning one or more of the underlying variables have not
been transcribed).

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of [`optimizer_model_expression`](@ref). Errors if such an
extension has not been written. 

By default only the expressions associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, infinite expressions are 
returned as a list corresponding to their supports. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the expression has multiple 
infinite parameter dependencies. The corresponding supports are obtained via 
`supports` using the same keyword arguments.

**Example**
```julia-repl
julia> optimizer_model_expression(my_expr) # finite expression
x(support: 1) - y
```
"""
function optimizer_model_expression(expr::JuMP.AbstractJuMPScalar; kwargs...)
    model = _model_from_expr(expr)
    if isnothing(model)
        return zero(JuMP.AffExpr) + JuMP.constant(expr)
    else
        key = optimizer_model_key(model)
        return optimizer_model_expression(expr, Val(key); kwargs...)
    end
end

"""
    expression_supports(optimizer_model::JuMP.Model, expr,
                        key::Val{ext_key_name}; [kwargs...])

Return the supports associated with the mappings of `expr` in `optimizer_model`.
This dispatches off of `key` which permits optimizer model extensions. This
should throw an error if `expr` is not associated with the variable mappings
stored in `optimizer_model`. Keyword arguments can be added as needed. Note that
if `expr` is a `GeneralVariableRef` this just dispatches to `variable_supports`.
"""
function expression_supports end

# fallback for unextended keys
function expression_supports(optimizer_model::JuMP.Model, expr, key; kwargs...)
  error("`constraint_supports` not implemented for optimizer model key " *
        "`$(typeof(key).parameters[1])` and/or expressions of type " *
        "`$(typeof(expr))`.")
end

# Variable reference expressions
function expression_supports(model::JuMP.Model, vref::GeneralVariableRef, key;
                             kwargs...)
    return variable_supports(model, dispatch_variable_ref(vref), key; kwargs...)
end

"""
    supports(expr::JuMP.AbstractJuMPScalar; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel,
             ndarray::Bool = false,
             kwargs...])

Return the support associated with `expr`. Errors if `expr` is
not associated with the constraint mappings stored in `optimizer_model`.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `expression_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the supports of infinite 
expressions are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the expression has multiple 
infinite parameter dependencies.

**Example**
```julia-repl
julia> supports(cref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(expr::JuMP.AbstractJuMPScalar; kwargs...)
    model = _model_from_expr(expr)
    if isnothing(model)
        return ()
    else
        key = optimizer_model_key(model)
        opt_model = optimizer_model(model)
        return expression_supports(opt_model, expr, Val(key); kwargs...)
    end
end

################################################################################
#                OPTIMIZER MODEL MAPPING METHODS (CONSTRAINTS)
################################################################################
"""
    optimizer_model_constraint(cref::InfOptConstraintRef,
                               key::Val{ext_key_name}; [kwargs...])

Return the reformulation constraint(s) stored in the optimizer model that correspond
to `cref`. This needs to be defined for extensions that implement a custom
optimizer model type. Principally, this is accomplished by typed the `key`
argument to `Val{ext_key_name}`. Keyword arguments can be added as needed.
"""
function optimizer_model_constraint end

# Fallback for unextended keys
function optimizer_model_constraint(
    cref::InfOptConstraintRef,
    key; 
    kwargs...
    )
    error("`optimizer_model_constraint` not implemented for optimizer model " *
          "key `$(typeof(key).parameters[1])`.")
end

"""
    optimizer_model_constraint(cref::InfOptConstraintRef; 
                               [label::Type{<:AbstractSupportLabel} = PublicLabel, 
                               ndarray::Bool = false,
                               kwargs...])

Return the reformulation constraint(s) stored in the optimizer model that correspond
to `cref`. Errors if no such constraint can be found in
the optimizer model.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of [`optimizer_model_constraint`](@ref). Errors if such an
extension has not been written. 

By default only the constraints associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, infinite constraints are 
returned as a list corresponding to their supports. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies. The corresponding supports are obtained via 
`supports` using the same keyword arguments.

**Example**
```julia-repl
julia> optimizer_model_constraint(c1) # finite constraint
c1 : x(support: 1) - y <= 3.0
```
"""
function optimizer_model_constraint(
    cref::InfOptConstraintRef; 
    kwargs...
    )
    key = optimizer_model_key(JuMP.owner_model(cref))
    return optimizer_model_constraint(cref, Val(key); kwargs...)
end

"""
    constraint_supports(optimizer_model::JuMP.Model, 
                        cref::InfOptConstraintRef,
                        key::Val{ext_key_name}; [kwargs...])

Return the supports associated with the mappings of `cref` in `optimizer_model`.
This dispatches off of `key` which permits optimizer model extensions. This
should throw an error if `cref` is not associated with the variable mappings
stored in `optimizer_model`. Keyword arguments can be added as needed.
"""
function constraint_supports end

# fallback for unextended keys
function constraint_supports(optimizer_model::JuMP.Model,
                             cref::InfOptConstraintRef,
                             key; kwargs...)
  error("`constraint_supports` not implemented for optimizer model key " *
        "`$(typeof(key).parameters[1])`.")
end

"""
    supports(cref::InfOptConstraintRef; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel,
             ndarray::Bool = false,
             kwargs...])

Return the support associated with `cref`. Errors if `cref` is
not associated with the constraint mappings stored in `optimizer_model`.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `constraint_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the supports of infinite 
constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

**Example**
```julia-repl
julia> supports(cref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(cref::InfOptConstraintRef; kwargs...)
    model = optimizer_model(JuMP.owner_model(cref))
    key = optimizer_model_key(JuMP.owner_model(cref))
    return constraint_supports(model, cref, Val(key); kwargs...)
end

################################################################################
#                             OPTIMIZATION METHODS
################################################################################
"""
    JuMP.optimize!(model::InfiniteModel; [kwargs...])

Extend `JuMP.optimize!` to optimize infinite models using the internal
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
function JuMP.optimize!(model::InfiniteModel; kwargs...)
    if !optimizer_model_ready(model)
        build_optimizer_model!(model; kwargs...)
    end
    JuMP.optimize!(optimizer_model(model))
    return
end
