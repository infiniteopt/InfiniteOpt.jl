```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

# JuMP Documentation
The docstrings below are sourced from JuMP.jl and provided for completeness.


## Model
```@docs
JuMP.Model
JuMP.Model()
JuMP.Model(::Any)
JuMP.NoOptimizer
JuMP.optimizer_with_attributes(::Any, ::Pair)
JuMP.direct_model
JuMP.unregister
```

## Variables
```@docs
JuMP.AbstractVariableRef
JuMP.add_variable
JuMP.@variable
JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)
JuMP.is_valid(::JuMP.Model, ::JuMP.VariableRef)
JuMP.set_name(::JuMP.VariableRef, ::String)
JuMP.owner_model(::AbstractVariableRef)
JuMP.index(::JuMP.VariableRef)
JuMP.num_variables(::JuMP.Model)
JuMP.name(::JuMP.VariableRef)
JuMP.variable_by_name(::JuMP.Model, ::String)
JuMP.all_variables(::JuMP.Model)
JuMP.has_lower_bound(::JuMP.VariableRef)
JuMP.lower_bound(::JuMP.VariableRef)
JuMP.set_lower_bound(::JuMP.VariableRef, ::Number)
JuMP.LowerBoundRef(::JuMP.VariableRef)
JuMP.delete_lower_bound(::JuMP.VariableRef)
JuMP.has_upper_bound(::JuMP.VariableRef)
JuMP.upper_bound(::JuMP.VariableRef)
JuMP.set_upper_bound(::JuMP.VariableRef, ::Number)
JuMP.UpperBoundRef(::JuMP.VariableRef)
JuMP.delete_upper_bound(::JuMP.VariableRef)
JuMP.is_fixed(::JuMP.VariableRef)
JuMP.fix_value(::JuMP.VariableRef)
JuMP.fix(::JuMP.VariableRef, ::Number)
JuMP.FixRef(::JuMP.VariableRef)
JuMP.unfix(::JuMP.VariableRef)
JuMP.start_value(::JuMP.VariableRef)
JuMP.set_start_value(::JuMP.VariableRef, ::Number)
JuMP.is_binary(::JuMP.VariableRef)
JuMP.set_binary(::JuMP.VariableRef)
JuMP.BinaryRef(::JuMP.VariableRef)
JuMP.unset_binary(::JuMP.VariableRef)
JuMP.is_integer(::JuMP.VariableRef)
JuMP.set_integer(::JuMP.VariableRef)
JuMP.IntegerRef(::JuMP.VariableRef)
JuMP.unset_integer(::JuMP.VariableRef)
```

## Containers 
```@docs
JuMP.Containers.SparseAxisArray
JuMP.Containers.DenseAxisArray
```

## Expressions
```@docs
JuMP.@expression
JuMP.add_to_expression!
```

## Objectives
```@docs
JuMP.@objective
JuMP.set_objective_function
JuMP.set_objective_sense(::JuMP.Model, ::MOI.OptimizationSense)
JuMP.objective_function(::JuMP.Model)
JuMP.objective_function_type(::JuMP.Model)
JuMP.objective_sense(::JuMP.Model)
JuMP.set_objective_coefficient(::JuMP.Model, ::JuMP.VariableRef, ::Real)
```

## [Constraints] (@id jump_constrs)
```@docs
JuMP.ScalarConstraint
JuMP.add_constraint(::JuMP.Model, ::JuMP.AbstractConstraint, ::String)
JuMP.@constraint
JuMP.owner_model(::JuMP.ConstraintRef)
JuMP.index(::JuMP.ConstraintRef)
JuMP.delete(::JuMP.Model, ::JuMP.ConstraintRef{JuMP.Model})
JuMP.is_valid(::JuMP.Model, ::JuMP.ConstraintRef{JuMP.Model})
JuMP.constraint_object
JuMP.name(::JuMP.ConstraintRef{JuMP.Model,<:JuMP._MOICON})
JuMP.set_name(::JuMP.ConstraintRef{JuMP.Model,<:JuMP._MOICON}, ::String)
JuMP.constraint_by_name
JuMP.num_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet})
JuMP.all_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet})
JuMP.list_of_constraint_types(::JuMP.Model)
JuMP.SecondOrderCone
JuMP.RotatedSecondOrderCone
JuMP.PSDCone
```

## Optimization
```@docs
JuMP.optimize!(::JuMP.Model, ::Any)
JuMP.set_silent(::JuMP.Model)
JuMP.unset_silent(::JuMP.Model)
JuMP.set_time_limit_sec(::JuMP.Model, ::Any)
JuMP.unset_time_limit_sec(::JuMP.Model)
JuMP.time_limit_sec(::JuMP.Model)
JuMP.bridge_constraints(::JuMP.Model)
JuMP.add_bridge(::JuMP.Model, ::Type{<:MOI.Bridges.AbstractBridge})
JuMP.backend(::JuMP.Model)
JuMP.mode(::JuMP.Model)
JuMP.solver_name(::JuMP.Model)
JuMP.set_optimizer_attribute(::JuMP.Model, ::String, ::Any)
JuMP.set_optimizer_attribute(::JuMP.Model, ::MOI.AbstractOptimizerAttribute, ::Any)
JuMP.set_optimizer_attributes(::JuMP.Model, ::Pair)
JuMP.get_optimizer_attribute(::JuMP.Model, ::String)
JuMP.get_optimizer_attribute(::JuMP.Model, ::MOI.AbstractOptimizerAttribute)
JuMP.set_optimizer(::JuMP.Model, ::Any)
JuMP.result_count(::JuMP.Model)
```

## Queries
```@docs
JuMP.termination_status(::JuMP.Model)
JuMP.raw_status(::JuMP.Model)
JuMP.primal_status(::JuMP.Model)
JuMP.dual_status(::JuMP.Model)
JuMP.solve_time(::JuMP.Model)
JuMP.has_values(::JuMP.Model)
JuMP.has_duals(::JuMP.Model)
JuMP.objective_bound(::JuMP.Model)
JuMP.objective_value(::JuMP.Model)
JuMP.dual_objective_value(::JuMP.Model)
JuMP.value(::JuMP.VariableRef)
JuMP.value(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON})
JuMP.dual(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON})
JuMP.shadow_price(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON})
JuMP.optimizer_index(::JuMP.VariableRef)
JuMP.optimizer_index(::JuMP.ConstraintRef{JuMP.Model})
JuMP.lp_sensitivity_report(::JuMP.Model)
JuMP.SensitivityReport
JuMP.reduced_cost(::JuMP.VariableRef)

```

# JuMP Section References
## Containers in macros
The `container` function encodes the logic for how containers are
constructed in JuMP's macros. The `@container` macro is available to create
containers independently of any JuMP model.
```@docs
JuMP.Containers.container
JuMP.Containers.default_container
JuMP.Containers.VectorizedProductIterator
JuMP.Containers.NestedIterator
JuMP.Containers.@container
```

In the [`@variable`](@ref) (resp. [`@constraint`](@ref)) macro, containers of
variables (resp. constraints) can be created with the following syntax:

* `name[index_set_1, index_set_2, ..., index_set_n]` creating an `n`-dimensional
  container of name `name`; or
* `[index_set_1, index_set_2, ..., index_set_n]` creating an *anonymous*
  `n`-dimensional container.

Each expression `index_set_i` can either be

* of the form `index_set` specifying that the `i`th index set of the container
  is `index_set`; or
* of the form `index_name=index_set` specifying that the `i`th index set of the
  container is `index_set` and allowing values used in the macro expression and
  keyword arguments to be expressions depending on the `index_name`.

The macro then creates the container using the
[`JuMP.Containers.container`](@ref) function with the following
arguments:

1. A function taking as argument the value of the indices and returning the
   value to be stored in the container, e.g. a variable for the
   [`@variable`](@ref) macro and a constraint for the [`@constraint`](@ref)
   macro.
2. An iterator over the indices of the container.
4. The value of the `container` keyword argument if given.

## Variables constrained on creation
!!! info
    When using JuMP in Direct mode, it may be required to constrain
    variables on creation instead of constraining free variables as the solver
    may only support variables constrained on creation. In Automatic and Manual
    modes, both ways of adding constraints on variables are equivalent.
    Indeed, during the copy of the cache to the optimizer, the choice of the
    constraints on variables that are copied as variables constrained on creation
    does not depend on how it was added to the cache.

All uses of the `@variable` macro documented so far translate to a separate
call for variable creation and adding of constraints.

For example, `@variable(model, x >= 0, Int)`, is equivalent to:
```julia
@variable(model, x)
set_lower_bound(x, 0.0)
@constraint(model, x in MOI.Integer())
```
Importantly, the bound and integrality constraints are added _after_ the
variable has been created.

However, some solvers require a constraining set _at creation time_.
We say that these variables are _constrained on creation_.

Use `in` within `@variable` to access the special syntax for constraining
variables on creation. For example, the following creates a vector of variables
constrained on creation to belong to the `SecondOrderCone`:
```julia-repl
julia> @variable(model, y[1:3] in SecondOrderCone())
3-element Array{VariableRef,1}:
 y[1]
 y[2]
 y[3]
```

For contrast, the more standard approach is as follows:
```julia-repl
julia> @variable(model, x[1:3])
3-element Array{VariableRef,1}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SecondOrderCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.SecondOrderCone(3)
```

The technical difference between the former and the latter is that the former
calls `MOI.add_constrained_variables` while the latter calls `MOI.add_variables`
and then `MOI.add_constraint`. This distinction is important only in
Direct mode, depending on the solver being used. It's often not
possible to delete the `SecondOrderCone` constraint if it was specified
at variable creation time.
