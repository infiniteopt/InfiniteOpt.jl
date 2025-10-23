# [Finite Parameters](@id finite_param_manual)
A technical manual for finite parameters in `InfiniteOpt`. See the respective 
[guide](@ref finite_param_docs) for more information.

## Definition
```@docs
@finite_parameter
FiniteParameter
build_parameter(::Function, ::Real)
add_parameter(::InfiniteModel, ::FiniteParameter, ::String)
FiniteParameterIndex
FiniteParameterRef
```

## Methods
Many methods are shared with independent infinite parameters since both 
finite and independent infinite parameters are scalar. See the infinite 
parameter [technical manual](@ref inf_par_manual) for the remainder of the 
methods available for finite parameters (i.e., any method typed for 
`ScalarParameterRef`s)
```@docs
JuMP.parameter_value(::FiniteParameterRef)
JuMP.set_parameter_value(::FiniteParameterRef, ::Real)
used_by_objective(::FiniteParameterRef)
core_object(::FiniteParameterRef)
```
