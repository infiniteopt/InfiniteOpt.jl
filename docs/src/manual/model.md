# [Infinite Models](@id infinite_model_manual)
A technical manual for infinite dimensional models. See the respective 
[guide](@ref infinite_model_docs) for more information.

## Models
```@docs
InfiniteModel
InfiniteModel()
JuMP.object_dictionary(::InfiniteModel)
has_internal_supports
Base.empty!(::InfiniteModel)
JuMP.set_optimize_hook(::InfiniteModel, ::Union{Function, Nothing})
parameter_group_indices
```

## Abstract Dependencies
```@docs
AbstractDataObject
AbstractInfOptIndex
ObjectIndex
```
