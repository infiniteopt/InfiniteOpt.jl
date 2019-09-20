# Variables
A guide and manual for the definition and use of variables in `InfiniteOpt`.

## Overview


## Basic Usage


## Infinite Variable Definition


## Point Variable Definition


## Global Variable Definition


## Manipulation


## Datatypes
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
InfOptVariable
InfiniteVariable
PointVariable
GlobalVariable
GeneralVariableRef
MeasureFiniteVariableRef
FiniteVariableRef
InfiniteVariableRef
PointVariableRef
GlobalVariableRef
```

## Methods
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Symbol)
JuMP.add_variable(::InfiniteModel, ::InfOptVariable, ::String)
@infinite_variable
@point_variable
@global_variable
JuMP.owner_model(::GeneralVariableRef)
JuMP.index(::GeneralVariableRef)
used_by_constraint(::InfOptVariableRef)
used_by_measure(::InfOptVariableRef)
used_by_objective(::InfOptVariableRef)
is_used(::InfOptVariableRef)
used_by_point_variable(::InfiniteVariableRef)
used_by_reduced_variable(::InfiniteVariableRef)
is_used(::InfiniteVariableRef)
JuMP.delete(::InfiniteModel, ::InfOptVariableRef)
JuMP.is_valid(::InfiniteModel, ::InfOptVariableRef)
JuMP.num_variables(::InfiniteModel)
JuMP.has_lower_bound(::InfOptVariableRef)
JuMP.lower_bound(::InfOptVariableRef)
JuMP.set_lower_bound(::InfOptVariableRef, ::Number)
JuMP.LowerBoundRef(::InfOptVariableRef)
JuMP.delete_lower_bound(::InfOptVariableRef)
JuMP.has_upper_bound(::InfOptVariableRef)
JuMP.upper_bound(::InfOptVariableRef)
JuMP.set_upper_bound(::InfOptVariableRef, ::Number)
JuMP.UpperBoundRef(::InfOptVariableRef)
JuMP.delete_upper_bound(::InfOptVariableRef)
JuMP.is_fixed(::InfOptVariableRef)
JuMP.fix_value(::InfOptVariableRef)
JuMP.fix(::InfOptVariableRef, ::Number; ::Bool)
JuMP.FixRef(::InfOptVariableRef)
JuMP.unfix(::InfOptVariableRef)
JuMP.start_value(::InfOptVariableRef)
JuMP.set_start_value(::InfOptVariableRef, ::Number)
JuMP.is_binary(::InfOptVariableRef)
JuMP.set_binary(::InfOptVariableRef)
JuMP.BinaryRef(::InfOptVariableRef)
JuMP.unset_binary(::InfOptVariableRef)
JuMP.is_integer(::InfOptVariableRef)
JuMP.set_integer(::InfOptVariableRef)
JuMP.IntegerRef(::InfOptVariableRef)
JuMP.unset_integer(::InfOptVariableRef)
JuMP.name(::InfOptVariableRef)
JuMP.set_name(::InfiniteVariableRef, ::String)
JuMP.set_name(::PointVariableRef, ::String)
JuMP.set_name(::GlobalVariableRef, ::String)
parameter_refs(::InfiniteVariableRef)
set_parameter_refs(::InfiniteVariableRef, ::Tuple)
add_parameter_ref(::InfiniteVariableRef,::Union{ParameterRef, AbstractArray{<:ParameterRef}})
infinite_variable_ref(::PointVariableRef)
parameter_values(::PointVariableRef)
JuMP.variable_by_name(::InfiniteModel, ::String)
JuMP.all_variables(::InfiniteModel)
```
