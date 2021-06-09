# [Objectives](@id obj_manual)
A technical manual for objective functions in `InfiniteOpt`. See the respective 
[guide](@ref obj_docs) for more information.

## Queries
```@docs
JuMP.objective_sense(::InfiniteModel)
JuMP.objective_function_type(::InfiniteModel)
JuMP.objective_function(::InfiniteModel)
objective_has_measures
```

## Modification
```@docs
JuMP.set_objective_function(::InfiniteModel, ::JuMP.AbstractJuMPScalar)
JuMP.set_objective_function(::InfiniteModel, ::Real)
JuMP.set_objective_sense(::InfiniteModel, ::MOI.OptimizationSense)
JuMP.set_objective(::InfiniteModel, ::MOI.OptimizationSense, ::Union{JuMP.AbstractJuMPScalar, Real})
JuMP.set_objective_coefficient(::InfiniteModel, ::GeneralVariableRef, ::Real)
```
