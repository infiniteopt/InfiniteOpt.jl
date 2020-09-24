# Derivatives
A guide and manual for the definition and use of derivatives in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview


## Basic Usage


## Advanced Information


## Numerical Evaluation Methods


## Datatypes
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
DerivativeIndex
DerivativeRef
AbstractDerivativeMethod
Integral
Derivative
```

## [Methods/Macros] (@id deriv_methods)
```@index
Pages   = ["derivative.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
deriv
build_derivative
add_derivative
default_eval_method
set_default_eval_method
derivative_argument(::DerivativeRef)
operator_parameter(::DerivativeRef)
eval_method(::DerivativeRef)
set_eval_method(::DerivativeRef, ::AbstractDerivativeMethod)
set_all_eval_methods
raw_parameter_refs(::DerivativeRef)
parameter_refs(::DerivativeRef)
parameter_list(::DerivativeRef)
set_start_value_function(::DerivativeRef,::Union{Real, Function})
reset_start_value_function(::DerivativeRef)
num_derivatives
all_derivatives
evaluate 
evaluate_all_derivatives
evaluate_derivative
```