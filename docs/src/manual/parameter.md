# [Infinite Parameters](@id inf_par_manual)
A technical manual for infinite parameters in `InfiniteOpt`. See the respective 
[guide](@ref inf_par_docs) for more information.

## Definition
### Macro
```@docs 
@infinite_parameter
```

### Independent Parameters
```@docs 
InfOptParameter
ScalarParameter
IndependentParameter
build_parameter(::Function, ::InfiniteScalarDomain)
add_parameter(::InfiniteModel, ::IndependentParameter, ::String)
ScalarParameterData
IndependentParameterIndex
IndependentParameterRef
```

### Dependent Parameters
```@docs
DependentParameters
add_parameters
MultiParameterData
DependentParametersIndex
DependentParameterIndex
DependentParameterRef
```

## Queries 
### General 
```@docs
parameter_by_name(::InfiniteModel,::String)
num_parameters
all_parameters
```

### Independent Parameters
```@docs
JuMP.name(::ScalarParameterRef)
infinite_domain(::IndependentParameterRef)
JuMP.has_lower_bound(::IndependentParameterRef)
JuMP.lower_bound(::IndependentParameterRef)
JuMP.has_upper_bound(::IndependentParameterRef)
JuMP.upper_bound(::IndependentParameterRef)
has_supports(::IndependentParameterRef)
num_supports(::IndependentParameterRef)
supports(::IndependentParameterRef)
has_internal_supports(::Union{IndependentParameterRef, DependentParameterRef})
significant_digits(::IndependentParameterRef)
derivative_method(::IndependentParameterRef)
is_used(::ScalarParameterRef)
used_by_infinite_variable(::IndependentParameterRef)
used_by_parameter_function(::IndependentParameterRef)
used_by_derivative(::IndependentParameterRef)
used_by_measure(::ScalarParameterRef)
used_by_constraint(::ScalarParameterRef)
parameter_group_int_index(::IndependentParameterRef)
core_object(::IndependentParameterRef)
```

### Dependent Parameters
```@docs
JuMP.name(::DependentParameterRef)
infinite_domain(::DependentParameterRef)
infinite_domain(::AbstractArray{<:DependentParameterRef})
JuMP.has_lower_bound(::DependentParameterRef)
JuMP.lower_bound(::DependentParameterRef)
JuMP.has_upper_bound(::DependentParameterRef)
JuMP.upper_bound(::DependentParameterRef)
has_supports(::DependentParameterRef)
has_supports(::AbstractArray{<:DependentParameterRef})
num_supports(::DependentParameterRef)
num_supports(::AbstractArray{<:DependentParameterRef})
supports(::DependentParameterRef)
supports(::AbstractArray{<:DependentParameterRef})
significant_digits(::DependentParameterRef)
derivative_method(::DependentParameterRef)
is_used(::DependentParameterRef)
used_by_infinite_variable(::DependentParameterRef)
used_by_parameter_function(::DependentParameterRef)
used_by_derivative(::DependentParameterRef)
used_by_measure(::DependentParameterRef)
used_by_constraint(::DependentParameterRef)
parameter_group_int_index(::DependentParameterRef)
core_object(::DependentParameterRef)
```

## Modification
### General
```@docs
fill_in_supports!(::InfiniteModel)
```

### Independent Parameters 
```@docs
JuMP.set_name(::ScalarParameterRef, ::String)
set_infinite_domain(::IndependentParameterRef,::InfiniteScalarDomain)
JuMP.set_lower_bound(::IndependentParameterRef, ::Real)
JuMP.set_upper_bound(::IndependentParameterRef,::Real)
add_supports(::IndependentParameterRef,::Union{Real, Vector{<:Real}})
set_supports(::IndependentParameterRef, ::Vector{<:Real})
delete_supports(::IndependentParameterRef)
generate_and_add_supports!(::IndependentParameterRef,::AbstractInfiniteDomain)
fill_in_supports!(::IndependentParameterRef)
JuMP.delete(::InfiniteModel, ::IndependentParameterRef)
```

### Dependent Parameters
```@docs
JuMP.set_name(::DependentParameterRef, ::String)
set_infinite_domain(::DependentParameterRef,::InfiniteScalarDomain)
set_infinite_domain(::AbstractArray{<:DependentParameterRef},::InfiniteArrayDomain)
JuMP.set_lower_bound(::DependentParameterRef,::Real)
JuMP.set_upper_bound(::DependentParameterRef,::Real)
add_supports(::AbstractArray{<:DependentParameterRef},::Vector{<:AbstractArray{<:Real}})
set_supports(::AbstractArray{<:DependentParameterRef},::Vector{<:AbstractArray{<:Real}})
delete_supports(::AbstractArray{<:DependentParameterRef})
generate_and_add_supports!(::AbstractArray{<:DependentParameterRef},::InfiniteArrayDomain)
fill_in_supports!(::AbstractArray{<:DependentParameterRef})
JuMP.delete(::InfiniteModel,::AbstractArray{<:DependentParameterRef})
```

## Generative Supports
```@docs
AbstractGenerativeInfo
NoGenerativeSupports
UniformGenerativeInfo
has_generative_supports(::IndependentParameterRef)
support_label(::AbstractGenerativeInfo)
generative_support_info(::IndependentParameterRef)
make_generative_supports
add_generative_supports
```
