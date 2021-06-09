# [Results](@id result_manual)
A technical manual for querying optimized `InfiniteOpt` models. See the 
respective [guide](@ref result_manual) for more information.

## Statuses 
```@docs
JuMP.termination_status(::InfiniteModel)
JuMP.raw_status(::InfiniteModel)
JuMP.primal_status(::InfiniteModel)
JuMP.dual_status(::InfiniteModel)
```

## General
```@docs
JuMP.solve_time(::InfiniteModel)
JuMP.simplex_iterations(::InfiniteModel)
JuMP.barrier_iterations(::InfiniteModel)
JuMP.node_count(::InfiniteModel)
JuMP.result_count(::InfiniteModel)
```

## Objective
```@docs
JuMP.objective_bound(::InfiniteModel)
JuMP.objective_value(::InfiniteModel)
JuMP.dual_objective_value(::InfiniteModel)
```

## Variables
```@docs
JuMP.has_values(::InfiniteModel)
JuMP.value(::GeneralVariableRef)
JuMP.reduced_cost(::GeneralVariableRef)
JuMP.optimizer_index(::GeneralVariableRef)
InfiniteOpt.map_value
InfiniteOpt.map_reduced_cost
InfiniteOpt.map_optimizer_index
```

## Constraints
```@docs
JuMP.has_duals(::InfiniteModel)
JuMP.value(::InfOptConstraintRef)
JuMP.optimizer_index(::InfOptConstraintRef)
JuMP.dual(::InfOptConstraintRef)
JuMP.shadow_price(::InfOptConstraintRef)
InfiniteOpt.map_dual
InfiniteOpt.map_shadow_price
```

## Expressions
```@docs
JuMP.value(::Union{JuMP.GenericAffExpr{<:Any, <:GeneralVariableRef}, JuMP.GenericQuadExpr{<:Any, <:GeneralVariableRef}})
```

## LP Sensitivity
```@docs
JuMP.lp_sensitivity_report(::InfiniteModel)
InfOptSensitivityReport 
```
