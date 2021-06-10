# [Optimization](@id opt_manual)
A technical manual for optimizing (solving) `InfiniteOpt` models. See the 
respective [guide](@ref opt_docs) for more information.

## Optimize
```@docs
JuMP.optimize!(::InfiniteModel)
```

## Optimizer Settings
```@docs
JuMP.set_optimizer(::InfiniteModel, ::Any)
JuMP.set_silent(::InfiniteModel)
JuMP.unset_silent(::InfiniteModel)
JuMP.set_time_limit_sec(::InfiniteModel, ::Any)
JuMP.unset_time_limit_sec(::InfiniteModel)
JuMP.time_limit_sec(::InfiniteModel)
JuMP.set_optimizer_attribute(::InfiniteModel, ::String, ::Any)
JuMP.set_optimizer_attribute(::InfiniteModel,::MOI.AbstractOptimizerAttribute,::Any)
JuMP.set_optimizer_attributes(::InfiniteModel, ::Pair)
JuMP.get_optimizer_attribute(::InfiniteModel, ::String)
JuMP.get_optimizer_attribute(::InfiniteModel,::MOI.AbstractOptimizerAttribute)
JuMP.add_bridge(::InfiniteModel, ::Type{<:MOI.Bridges.AbstractBridge})
```

## Optimizer Queries
```@docs
JuMP.solver_name(model::InfiniteModel)
JuMP.backend(model::InfiniteModel)
JuMP.mode(model::InfiniteModel)
JuMP.bridge_constraints(::InfiniteModel)
```

## Optimizer Model API 
```@docs
optimizer_model
set_optimizer_model
optimizer_model_key(::InfiniteModel)
optimizer_model_key(::JuMP.Model)
build_optimizer_model!(::InfiniteModel)
build_optimizer_model!
clear_optimizer_model_build!(::InfiniteModel)
clear_optimizer_model_build!(::JuMP.Model)
InfiniteOpt.add_infinite_model_optimizer
optimizer_model_variable(::GeneralVariableRef)
optimizer_model_variable
supports(::Union{DecisionVariableRef, MeasureRef})
InfiniteOpt.variable_supports
optimizer_model_expression(::JuMP.AbstractJuMPScalar)
optimizer_model_expression
supports(::JuMP.AbstractJuMPScalar)
InfiniteOpt.expression_supports
InfiniteOpt.optimizer_model_constraint(::InfOptConstraintRef)
optimizer_model_constraint
supports(::InfOptConstraintRef)
InfiniteOpt.constraint_supports
optimizer_model_ready
set_optimizer_model_ready
```
