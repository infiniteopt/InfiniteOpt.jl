# [Variables](@id var_manual)
A technical manual for variables in `InfiniteOpt`. See the respective 
[guide](@ref var_docs) for more information.


## Definition
Note that the principle way for defining variables is by using 
[`JuMP.@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
which originates from `JuMP.jl`.

### Infinite Variables
```@docs
InfOptVariableType
Infinite
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Infinite)
JuMP.add_variable(::InfiniteModel, ::InfiniteVariable, ::String)
InfiniteVariable
restrict
VariableData
InfiniteVariableIndex
InfiniteVariableRef
InfiniteOpt.Collections.VectorTuple
```

### Semi-Infinite Variables
```@docs
SemiInfinite
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::SemiInfinite)
JuMP.build_variable(::Function, ::GeneralVariableRef, ::Vector{Float64}, ::RestrictedDomainInfo)
JuMP.add_variable(::InfiniteModel, ::SemiInfiniteVariable, ::String)
SemiInfiniteVariable
RestrictedDomainInfo
SemiInfiniteVariableIndex
SemiInfiniteVariableRef
```

### Point Variables
```@docs
Point
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Point)
JuMP.add_variable(::InfiniteModel, ::PointVariable, ::String)
PointVariable
PointVariableIndex
PointVariableRef
```

### Finite Variables
Note that finite variables simply correspond to using 
[`JuMP.ScalarVariable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.ScalarVariable) 
which originates from `JuMP.jl` as well. In other words, these are defined via 
`JuMP.@variable` without specifying any [`InfOptVariableType`](@ref).
```@docs
JuMP.add_variable(::InfiniteModel, ::JuMP.ScalarVariable, ::String)
FiniteVariableIndex
FiniteVariableRef
```

## Queries
### General
```@docs
JuMP.name(::DecisionVariableRef)
JuMP.variable_by_name(::InfiniteModel, ::String)
JuMP.num_variables(::InfiniteModel)
JuMP.all_variables(::InfiniteModel)
JuMP.has_lower_bound(::UserDecisionVariableRef)
JuMP.lower_bound(::UserDecisionVariableRef)
JuMP.LowerBoundRef(::UserDecisionVariableRef)
JuMP.has_upper_bound(::UserDecisionVariableRef)
JuMP.upper_bound(::UserDecisionVariableRef)
JuMP.UpperBoundRef(::UserDecisionVariableRef)
JuMP.is_fixed(::UserDecisionVariableRef)
JuMP.fix_value(::UserDecisionVariableRef)
JuMP.FixRef(::UserDecisionVariableRef)
JuMP.start_value(::UserDecisionVariableRef)
JuMP.is_binary(::UserDecisionVariableRef)
JuMP.BinaryRef(::UserDecisionVariableRef)
JuMP.is_integer(::UserDecisionVariableRef)
JuMP.IntegerRef(::UserDecisionVariableRef)
is_used(::DecisionVariableRef)
used_by_constraint(::DecisionVariableRef)
used_by_measure(::DecisionVariableRef)
used_by_objective(::DecisionVariableRef)
```

### Infinite Variables
```@docs
parameter_refs(::InfiniteVariableRef)
parameter_list(::InfiniteVariableRef)
raw_parameter_refs(::InfiniteVariableRef)
is_used(::Union{InfiniteVariableRef, DerivativeRef})
used_by_derivative(::Union{DerivativeRef, InfiniteVariableRef})
used_by_point_variable(::Union{InfiniteVariableRef, DerivativeRef})
used_by_semi_infinite_variable(::Union{InfiniteVariableRef, DerivativeRef})
parameter_group_int_indices(::InfiniteVariableRef)
core_object(::InfiniteVariableRef)
```

### Semi-Infinite Variables
```@docs
infinite_variable_ref(::SemiInfiniteVariableRef)
parameter_refs(::SemiInfiniteVariableRef)
parameter_list(::SemiInfiniteVariableRef)
raw_parameter_refs(::SemiInfiniteVariableRef)
eval_support(::SemiInfiniteVariableRef)
is_used(::SemiInfiniteVariableRef)
used_by_derivative(::SemiInfiniteVariableRef)
parameter_group_int_indices(::SemiInfiniteVariableRef)
core_object(::SemiInfiniteVariableRef)
```

### Point Variables
```@docs
infinite_variable_ref(::PointVariableRef)
parameter_values(::PointVariableRef)
raw_parameter_values(::PointVariableRef)
core_object(::PointVariableRef)
```

### Finite Variables
```@docs
core_object(::FiniteVariableRef)
```

## Modification
### General
```@docs
JuMP.set_name(::DecisionVariableRef, ::String)
JuMP.set_lower_bound(::UserDecisionVariableRef, ::Union{<:Real, <:Function})
JuMP.delete_lower_bound(::UserDecisionVariableRef)
JuMP.set_upper_bound(::UserDecisionVariableRef, ::Union{<:Real, <:Function})
JuMP.delete_upper_bound(::UserDecisionVariableRef)
JuMP.fix(::UserDecisionVariableRef, ::Union{<:Real, <:Function}; ::Bool)
JuMP.unfix(::UserDecisionVariableRef)
JuMP.set_start_value(::UserDecisionVariableRef, ::Union{<:Real, <:Function})
JuMP.set_binary(::UserDecisionVariableRef)
JuMP.unset_binary(::UserDecisionVariableRef)
JuMP.set_integer(::UserDecisionVariableRef)
JuMP.unset_integer(::UserDecisionVariableRef)
JuMP.relax_integrality(::InfiniteModel)
JuMP.set_start_values(::InfiniteModel)
JuMP.delete(::InfiniteModel, ::DecisionVariableRef)
```

### Infinite Variables
```@docs
constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)
```

### Semi-Infinite Variables
```@docs
JuMP.set_name(::SemiInfiniteVariableRef,::String)
```
