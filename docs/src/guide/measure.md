# Measures
A guide and manual for defining and using measures in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview


## Basic Usage


## Theoretical Abstraction


## Preset Evaluations


## Custom Evaluations


## Expansion


## Reduced Infinite Variables


## Datatypes
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractMeasureData
DiscreteMeasureData
MultiDiscreteMeasureData
Measure
MeasureRef
AbstractReducedInfo
ReducedInfiniteInfo
ReducedInfiniteVariableRef
```

## Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:function]
```
```@docs
expect
support_sum
measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})
DiscreteMeasureData(::ParameterRef, ::Vector{<:Number}, ::Vector{<:Number})
DiscreteMeasureData(::AbstractArray{<:ParameterRef}, ::Vector{<:Number}, ::Vector{<:AbstractArray})
measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)
add_measure
measure_function
measure_data
used_by_constraint(::MeasureRef)
used_by_measure(::MeasureRef)
used_by_objective(::MeasureRef)
is_used(::MeasureRef)
JuMP.is_valid(::InfiniteModel, ::MeasureRef)
JuMP.delete(::InfiniteModel, ::MeasureRef)
JuMP.name(::MeasureRef)
JuMP.set_name(::MeasureRef, ::String)
expand
expand_all_measures!
infinite_variable_ref(::ReducedInfiniteVariableRef)
eval_supports(::ReducedInfiniteVariableRef)
parameter_refs(::ReducedInfiniteVariableRef)
JuMP.name(::ReducedInfiniteVariableRef)
JuMP.has_lower_bound(::ReducedInfiniteVariableRef)
JuMP.lower_bound(::ReducedInfiniteVariableRef)
JuMP.LowerBoundRef(::ReducedInfiniteVariableRef)
JuMP.has_upper_bound(::ReducedInfiniteVariableRef)
JuMP.upper_bound(::ReducedInfiniteVariableRef)
JuMP.UpperBoundRef(::ReducedInfiniteVariableRef)
JuMP.is_fixed(::ReducedInfiniteVariableRef)
JuMP.fix_value(::ReducedInfiniteVariableRef)
JuMP.FixRef(::ReducedInfiniteVariableRef)
JuMP.start_value(::ReducedInfiniteVariableRef)
JuMP.is_binary(::ReducedInfiniteVariableRef)
JuMP.BinaryRef(::ReducedInfiniteVariableRef)
JuMP.is_integer(::ReducedInfiniteVariableRef)
JuMP.IntegerRef(::ReducedInfiniteVariableRef)
used_by_constraint(::ReducedInfiniteVariableRef)
used_by_measure(::ReducedInfiniteVariableRef)
JuMP.is_valid(::InfiniteModel, ::ReducedInfiniteVariableRef)
JuMP.delete(::InfiniteModel, ::ReducedInfiniteVariableRef)
```

## MeasureEvalMethods Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt.MeasureEvalMethods]
Order   = [:function]
```

```@docs
InfiniteOpt.MeasureEvalMethods.generate_measure_data
InfiniteOpt.MeasureEvalMethods.mc_sampling(::Number, ::Number, ::Int)
InfiniteOpt.MeasureEvalMethods.mc_sampling(::Distributions.UnivariateDistribution, ::InfiniteOpt.ParameterRef, ::Int)
InfiniteOpt.MeasureEvalMethods.gauss_legendre
InfiniteOpt.MeasureEvalMethods.gauss_hermite
InfiniteOpt.MeasureEvalMethods.gauss_laguerre
InfiniteOpt.MeasureEvalMethods.infinite_transform
```
