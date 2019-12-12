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
JuMP.ScalarConstraint
JuMP.@constraint
JuMP.SecondOrderCone
JuMP.RotatedSecondOrderCone
JuMP.PSDCone
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
