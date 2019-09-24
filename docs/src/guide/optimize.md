# Optimization
A guide and manual for optimizing (solving) `InfiniteOpt` models.

## Overview


## Basic Usage


## Optimizer Models


## Optimizer Settings


## Methods
```@index
Pages   = ["optimize.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:function]
```
```@docs
optimizer_model
set_optimizer_model
optimizer_model_key
build_optimizer_model!(::InfiniteModel)
build_optimizer_model!
optimizer_model_ready
set_optimizer_model_ready
JuMP.bridge_constraints(::InfiniteModel)
JuMP.add_bridge(::InfiniteModel, ::Type{<:MOI.Bridges.AbstractBridge})
JuMP.set_optimizer(::InfiniteModel, ::JuMP.OptimizerFactory)
JuMP.set_silent(::InfiniteModel)
JuMP.unset_silent(::InfiniteModel)
JuMP.set_parameter(::InfiniteModel, ::Any, ::Any)
JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory})
```
