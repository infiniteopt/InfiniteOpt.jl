# [Variables](@id var_manual)
A technical manual for variables in `InfiniteOpt`. See the respective 
[guide](@ref var_docs) for more information.


## Definition
Note that the principle way for defining variables is by using 
[`JuMP.@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
which originates from `JuMP.jl`.

### Infinite
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

### Semi-Infinite
```@docs
SemiInfinite
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::SemiInfinite)
JuMP.build_variable(::Function, ::GeneralVariableRef, ::Dict{Int, Float64})
JuMP.add_variable(::InfiniteModel, ::SemiInfiniteVariable, ::String)
SemiInfiniteVariable
SemiInfiniteVariableIndex
SemiInfiniteVariableRef
```

### Point
```@docs
Point
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Point)
JuMP.add_variable(::InfiniteModel, ::PointVariable, ::String)
PointVariable
PointVariableIndex
PointVariableRef
```

### Finite 
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

### Infinite
```@docs
start_value_function(::Union{InfiniteVariableRef, DerivativeRef})
parameter_refs(::InfiniteVariableRef)
parameter_list(::InfiniteVariableRef)
raw_parameter_refs(::InfiniteVariableRef)
is_used(::Union{InfiniteVariableRef, DerivativeRef})
used_by_derivative(::Union{DerivativeRef, InfiniteVariableRef})
used_by_point_variable(::Union{InfiniteVariableRef, DerivativeRef})
used_by_semi_infinite_variable(::Union{InfiniteVariableRef, DerivativeRef})
```

### Semi-Infinite
```@docs
JuMP.has_lower_bound(::SemiInfiniteVariableRef)
JuMP.lower_bound(::SemiInfiniteVariableRef)
JuMP.LowerBoundRef(::SemiInfiniteVariableRef)
JuMP.has_upper_bound(::SemiInfiniteVariableRef)
JuMP.upper_bound(::SemiInfiniteVariableRef)
JuMP.UpperBoundRef(::SemiInfiniteVariableRef)
JuMP.is_fixed(::SemiInfiniteVariableRef)
JuMP.fix_value(::SemiInfiniteVariableRef)
JuMP.FixRef(::SemiInfiniteVariableRef)
start_value_function(::SemiInfiniteVariableRef)
JuMP.is_binary(::SemiInfiniteVariableRef)
JuMP.BinaryRef(::SemiInfiniteVariableRef)
JuMP.is_integer(::SemiInfiniteVariableRef)
JuMP.IntegerRef(::SemiInfiniteVariableRef)
infinite_variable_ref(::SemiInfiniteVariableRef)
parameter_refs(::SemiInfiniteVariableRef)
parameter_list(::SemiInfiniteVariableRef)
raw_parameter_refs(::SemiInfiniteVariableRef)
eval_supports(::SemiInfiniteVariableRef)
is_used(::SemiInfiniteVariableRef)
used_by_derivative(::SemiInfiniteVariableRef)
```

### Point
```@docs
infinite_variable_ref(::PointVariableRef)
parameter_values(::PointVariableRef)
raw_parameter_values(::PointVariableRef)
```

## Modification
### General
```@docs
JuMP.set_name(::DecisionVariableRef, ::String)
JuMP.set_lower_bound(::UserDecisionVariableRef, ::Real)
JuMP.delete_lower_bound(::UserDecisionVariableRef)
JuMP.set_upper_bound(::UserDecisionVariableRef, ::Real)
JuMP.delete_upper_bound(::UserDecisionVariableRef)
JuMP.fix(::UserDecisionVariableRef, ::Real; ::Bool)
JuMP.unfix(::UserDecisionVariableRef)
JuMP.set_start_value(::UserDecisionVariableRef, ::Real)
JuMP.set_binary(::UserDecisionVariableRef)
JuMP.unset_binary(::UserDecisionVariableRef)
JuMP.set_integer(::UserDecisionVariableRef)
JuMP.unset_integer(::UserDecisionVariableRef)
JuMP.relax_integrality(::InfiniteModel)
JuMP.delete(::InfiniteModel, ::DecisionVariableRef)
```

### Infinite
```@docs
set_start_value_function(::InfiniteVariableRef, ::Union{Real, Function})
reset_start_value_function(::InfiniteVariableRef)
constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)
```

### Semi-Infinite
```@docs
JuMP.set_name(::SemiInfiniteVariableRef,::String)
```
