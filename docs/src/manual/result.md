# [Results](@id result_manual)
A technical manual for querying optimized `InfiniteOpt` models. See the 
respective [guide](@ref result_manual) for more information.

## Statuses 
```@docs
JuMP.termination_status(::InfiniteModel)
JuMP.raw_status(::InfiniteModel)
JuMP.primal_status(::InfiniteModel)
JuMP.dual_status(::InfiniteModel)
JuMP.is_solved_and_feasible(::InfiniteModel)
```

## General
```@docs
JuMP.solve_time(::InfiniteModel)
JuMP.relative_gap(::InfiniteModel)
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
```

## [InfiniteInterpolate](@id infiniteInterpolate)
```@docs
JuMP.value(::GeneralVariableRef, ::Union{typeof(Interpolations.constant_interpolation), typeof(Interpolations.cubic_spline_interpolation), typeof(Interpolations.linear_interpolation)})
JuMP.value(::GeneralVariableRef, ::Interpolations.Degree)
```

## Constraints
```@docs
JuMP.has_duals(::InfiniteModel)
JuMP.value(::InfOptConstraintRef)
JuMP.optimizer_index(::InfOptConstraintRef)
JuMP.dual(::InfOptConstraintRef)
JuMP.shadow_price(::InfOptConstraintRef)
```

## Expressions
```@docs
JuMP.value(::Union{JuMP.GenericAffExpr{Float64, GeneralVariableRef}, JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, JuMP.GenericNonlinearExpr{GeneralVariableRef}})
```

## LP Sensitivity
```@docs
JuMP.lp_sensitivity_report(::InfiniteModel)
InfOptSensitivityReport 
```

## Transformation Backend Extension API
```@docs
map_value(::Any, ::AbstractTransformationBackend)
map_infinite_parameter_value
map_reduced_cost(::GeneralVariableRef, ::AbstractTransformationBackend)
map_optimizer_index(::GeneralVariableRef, ::AbstractTransformationBackend)
map_dual(::InfOptConstraintRef, ::AbstractTransformationBackend)
map_shadow_price(::InfOptConstraintRef, ::AbstractTransformationBackend)
map_optimizer_index(::InfOptConstraintRef, ::AbstractTransformationBackend)
JuMP.lp_sensitivity_report(::AbstractTransformationBackend)
```