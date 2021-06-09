# [Derivative Operators](@id deriv_manual)
A technical manual for derivatives in `InfiniteOpt`. See the respective 
[guide](@ref deriv_docs) for more information.

## Definition
```@docs
deriv
∂
@deriv
@∂
Deriv
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Deriv)
build_derivative
Derivative
add_derivative
DerivativeIndex
DerivativeRef
```

## Queries
```@docs
derivative_argument(::DerivativeRef)
operator_parameter(::DerivativeRef)
num_derivatives
all_derivatives
parameter_refs(::DerivativeRef)
parameter_list(::DerivativeRef)
raw_parameter_refs(::DerivativeRef)
```

## Modification
```@docs
set_start_value_function(::DerivativeRef,::Union{Real, Function})
reset_start_value_function(::DerivativeRef)
```

## Evaluation 
```@docs
AbstractDerivativeMethod
GenerativeDerivativeMethod
OrthogonalCollocation
NonGenerativeDerivativeMethod
FiniteDifference
FDTechnique
Forward
Central
Backward
derivative_method(::DerivativeRef)
set_derivative_method(::IndependentParameterRef, ::NonGenerativeDerivativeMethod)
set_derivative_method(::DependentParameterRef, ::AbstractDerivativeMethod)
set_all_derivative_methods
evaluate(::DerivativeRef)
evaluate_all_derivatives!
has_derivative_constraints(::DerivativeRef)
derivative_constraints(::DerivativeRef)
delete_derivative_constraints(::DerivativeRef)
evaluate_derivative
generative_support_info(::AbstractDerivativeMethod)
support_label(::AbstractDerivativeMethod)
InfiniteOpt.make_reduced_expr
```
