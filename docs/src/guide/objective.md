# Objectives
A guide and manual for specifying and modifying objective functions in
`InfiniteOpt`.

## Overview


## Basic Usage


## Queries


## Modification


## Methods
```@index
Pages   = ["objective.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
JuMP.set_objective_function(::InfiniteModel, ::JuMP.AbstractJuMPScalar)
JuMP.set_objective_function(::InfiniteModel, ::Real)
JuMP.set_objective_sense(::InfiniteModel, ::MOI.OptimizationSense)
JuMP.set_objective(::InfiniteModel, ::MOI.OptimizationSense, ::Union{JuMP.AbstractJuMPScalar, Real})
JuMP.objective_sense(::InfiniteModel)
JuMP.objective_function_type(::InfiniteModel)
JuMP.objective_function(::InfiniteModel)
JuMP.set_objective_coefficient(::InfiniteModel, ::GeneralVariableRef, ::Real)
```
