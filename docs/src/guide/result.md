# Results
A guide and manual to querying optimized `InfiniteOpt` models.

## Overview


## Basic Usage


## Termination Queries


## Variable Queries


## Constraint Queries


## Methods
```@index
Pages   = ["optimize.md"]
Modules = [JuMP, InfiniteOpt, InfiniteOpt.TranscriptionOpt]
Order   = [:function]
```
```@docs
JuMP.termination_status(::InfiniteModel)
JuMP.raw_status(::InfiniteModel)
JuMP.primal_status(::InfiniteModel)
JuMP.dual_status(::InfiniteModel)
JuMP.solve_time(::InfiniteModel)
JuMP.has_values(::InfiniteModel)
JuMP.objective_bound(::InfiniteModel)
JuMP.objective_value(::InfiniteModel)
JuMP.value(::GeneralVariableRef)
JuMP.value(::GeneralConstraintRef)
map_value
InfiniteOpt.map_value(::InfiniteOpt.FiniteVariableRef, ::Val{:TransData})
InfiniteOpt.map_value(::InfiniteOpt.InfiniteVariableRef, ::Val{:TransData})
InfiniteOpt.map_value(::InfiniteOpt.FiniteConstraintRef, ::Val{:TransData})
InfiniteOpt.map_value(::InfiniteOpt.GeneralConstraintRef, ::Val{:TransData})
JuMP.optimizer_index(::GeneralVariableRef)
JuMP.optimizer_index(::GeneralConstraintRef)
map_optimizer_index
InfiniteOpt.map_optimizer_index(::InfiniteOpt.FiniteVariableRef, ::Val{:TransData})
InfiniteOpt.map_optimizer_index(::InfiniteOpt.InfiniteVariableRef, ::Val{:TransData})
InfiniteOpt.map_optimizer_index(::InfiniteOpt.FiniteConstraintRef, ::Val{:TransData})
InfiniteOpt.map_optimizer_index(::InfiniteOpt.GeneralConstraintRef, ::Val{:TransData})
JuMP.dual(::GeneralConstraintRef)
map_dual
InfiniteOpt.map_dual(::InfiniteOpt.FiniteConstraintRef, ::Val{:TransData})
InfiniteOpt.map_dual(::InfiniteOpt.InfiniteConstraintRef, ::Val{:TransData})
JuMP.shadow_price(::GeneralConstraintRef)
map_shadow_price
InfiniteOpt.map_shadow_price(::InfiniteOpt.FiniteConstraintRef, ::Val{:TransData})
InfiniteOpt.map_shadow_price(::InfiniteOpt.InfiniteConstraintRef, ::Val{:TransData})
```
