# Infinite Parameters
A guide and manual to the definition and use of infinite parameters in
`InfiniteOpt`.

## Overview


## Basic Usage


## Infinite Sets


## Parameter Definition


## Parameter Use  


## Datatypes
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractInfiniteSet
IntervalSet
DistributionSet
InfOptParameter
ParameterRef
```

## Methods/Macros
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
IntervalSet(::Number, ::Number)
add_parameter
build_parameter
@infinite_parameter
used_by_constraint(::ParameterRef)
used_by_measure(::ParameterRef)
used_by_variable(::ParameterRef)
is_used(::ParameterRef)
JuMP.delete(::InfiniteModel, ::ParameterRef)
JuMP.is_valid(::InfiniteModel, ::ParameterRef)
JuMP.name(::ParameterRef)
JuMP.set_name(::ParameterRef, ::String)
num_parameters(::InfiniteModel)
infinite_set(::ParameterRef)
set_infinite_set(::ParameterRef, ::AbstractInfiniteSet)
JuMP.has_lower_bound(::ParameterRef)
JuMP.lower_bound(::ParameterRef)
JuMP.set_lower_bound(::ParameterRef, ::Number)
JuMP.has_upper_bound(::ParameterRef)
JuMP.upper_bound(::ParameterRef)
JuMP.set_upper_bound(::ParameterRef, ::Number)
num_supports(::ParameterRef)
has_supports(::ParameterRef)
supports(::ParameterRef)
supports(::AbstractArray{<:ParameterRef})
set_supports(::ParameterRef, ::Vector{<:Number})
add_supports(::ParameterRef, ::Union{Number, Vector{<:Number}})
delete_supports(::ParameterRef)
group_id(::ParameterRef)
group_id(::AbstractArray{<:ParameterRef})
is_independent(::ParameterRef)
set_independent(::ParameterRef)
unset_independent(::ParameterRef)
parameter_by_name(::InfiniteModel, ::String)
all_parameters(::InfiniteModel)
```
