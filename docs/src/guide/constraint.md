# Constraints
A guide and manual for defining and manipulating constraints in `InfiniteOpt`.

## Overview


## Basic Usage


## Data Structure


## Modification


## Datatypes
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
BoundedScalarConstraint
GeneralConstraintRef
InfiniteConstraintRef
FiniteConstraintRef
MeasureConstraintRef
```

## Methods/Macros
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
JuMP.build_constraint(::Function, ::InfiniteExpr, ::MOI.AbstractScalarSet)
JuMP.add_constraint(::InfiniteModel, ::JuMP.AbstractConstraint)
JuMP.owner_model(::GeneralConstraintRef)
JuMP.index(::GeneralConstraintRef)
JuMP.constraint_object(::GeneralConstraintRef)
JuMP.name(::GeneralConstraintRef)
JuMP.set_name(::GeneralConstraintRef, ::String)
JuMP.is_valid(::InfiniteModel, ::GeneralConstraintRef)
JuMP.delete(::InfiniteModel, ::GeneralConstraintRef)
JuMP.set_normalized_rhs(::GeneralConstraintRef, ::Real)
JuMP.normalized_rhs(::GeneralConstraintRef)
JuMP.add_to_function_constant(::GeneralConstraintRef, ::Real)
JuMP.set_normalized_coefficient(::GeneralConstraintRef, ::GeneralVariableRef, ::Real)
JuMP.normalized_coefficient(::GeneralConstraintRef, ::GeneralVariableRef)
JuMP.constraint_by_name(::InfiniteModel, ::String)
JuMP.list_of_constraint_types(::InfiniteModel)
JuMP.num_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet})
JuMP.num_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar})
JuMP.num_constraints(::InfiniteModel, ::Type{<:MOI.AbstractSet})
JuMP.num_constraints(::InfiniteModel)
JuMP.all_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet})
JuMP.all_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar})
JuMP.all_constraints(::InfiniteModel, ::Type{<:MOI.AbstractSet})
JuMP.all_constraints(::InfiniteModel)
```
