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
```

## Abstract Dependencies
```@docs
AbstractDataObject
AbstractInfOptIndex
ObjectIndex
```
