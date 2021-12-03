# [Constraints](@id constr_manual)
A technical manual for constraints in `InfiniteOpt`. See the respective 
[guide](@ref constr_docs) for more information.

## Definition
Note that constraints are defined principally with 
[`JuMP.@constraint`](https://jump.dev/JuMP.jl/v0.22/reference/constraints/#JuMP.@constraint) 
which originates from `JuMP.jl`. Below we show build methods for 
`DomainRestrictedConstraint`s, but any `JuMP.AbstractConstraint` can be used.
```@docs
DomainRestrictions
JuMP.build_constraint(::Function, ::Any, ::Any, ::DomainRestrictions)
DomainRestrictedConstraint
JuMP.add_constraint(::InfiniteModel, ::JuMP.AbstractConstraint, ::String)
ConstraintData
InfOptConstraintIndex
InfOptConstraintRef
```

## Queries
```@docs
JuMP.owner_model(::InfOptConstraintRef)
JuMP.index(::InfOptConstraintRef)
JuMP.constraint_object(::InfOptConstraintRef)
JuMP.name(::InfOptConstraintRef)
JuMP.constraint_by_name(::InfiniteModel, ::String)
JuMP.list_of_constraint_types(::InfiniteModel)
JuMP.num_constraints(::InfiniteModel, ::Any, ::Any)
JuMP.all_constraints(::InfiniteModel, ::Any, ::Any)
JuMP.is_valid(::InfiniteModel, ::InfOptConstraintRef)
parameter_refs(::InfOptConstraintRef)
has_domain_restrictions
domain_restrictions
JuMP.normalized_rhs(::InfOptConstraintRef)
JuMP.normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef)
```

## Modification
```@docs
JuMP.set_name(::InfOptConstraintRef, ::String)
set_domain_restrictions
add_domain_restrictions
delete_domain_restrictions
JuMP.set_normalized_rhs(::InfOptConstraintRef, ::Real)
JuMP.add_to_function_constant(::InfOptConstraintRef, ::Real)
JuMP.set_normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef, ::Real)
JuMP.delete(::InfiniteModel, ::InfOptConstraintRef)
```
