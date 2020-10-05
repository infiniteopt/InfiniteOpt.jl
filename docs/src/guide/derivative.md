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
Derivative
AbstractDerivativeMethod
GenerativeDerivativeMethod
OrthogonalCollocation
OCQuadrature
Lobatto
NonGenerativeDerivateMethod
FiniteDifference
FDTechniques
FDForward
FDCentral
FDBackward
```

## [Methods/Macros] (@id deriv_methods)
```@index
Pages   = ["derivative.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@deriv
deriv
@derivative_variable
build_derivative
add_derivative
derivative_argument(::DerivativeRef)
operator_parameter(::DerivativeRef)
derivative_method(::DerivativeRef)
raw_parameter_refs(::DerivativeRef)
parameter_refs(::DerivativeRef)
parameter_list(::DerivativeRef)
set_start_value_function(::DerivativeRef,::Union{Real, Function})
reset_start_value_function(::DerivativeRef)
num_derivatives
all_derivatives
set_derivative_method(::IndependentParameterRef, ::AbstractDerivativeMethod)
set_derivative_method(::DependentParameterRef, ::AbstractDerivativeMethod)
set_all_derivative_methods
evaluate(::DerivativeRef)
evaluate_all_derivatives!
evaluate_derivative
InfiniteOpt.support_label(::AbstractDerivativeMethod)
InfiniteOPt.generate_derivative_supports
InfiniteOpt.add_derivative_supports(::IndependentParameterRef)
InfiniteOpt.make_reduced_expr
```