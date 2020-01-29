```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

# JuMP Documentation

## Model
```@docs
JuMP.Model
JuMP.Model()
JuMP.Model(::OptimizerFactory)
JuMP.with_optimizer
JuMP.NoOptimizer
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

## Constraints
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
JuMP.optimize!(::JuMP.Model, ::Union{Nothing, JuMP.OptimizerFactory})
JuMP.set_parameter(::JuMP.Model, ::Any, ::Any)
JuMP.set_silent(::JuMP.Model)
JuMP.unset_silent(::JuMP.Model)
JuMP.bridge_constraints(::JuMP.Model)
JuMP.add_bridge(::JuMP.Model, ::Type{<:MOI.Bridges.AbstractBridge})
```
