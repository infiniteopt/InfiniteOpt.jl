# JuMP Documentation

```@docs
JuMP.Model
JuMP.Model()
JuMP.Model(::OptimizerFactory)
JuMP.with_optimizer
JuMP.@expression
JuMP.add_to_expression!
JuMP.@objective
JuMP.ScalarConstraint
JuMP.@constraint
JuMP.SecondOrderCone
JuMP.RotatedSecondOrderCone
JuMP.PSDCone
JuMP.NoOptimizer
JuMP.delete(::Model, ::VariableRef)
JuMP.is_valid(::JuMP.Model, ::JuMP.VariableRef)
JuMP.set_name(::VariableRef, ::String)
JuMP.add_variable
JuMP.@variable
JuMP.owner_model(::AbstractVariableRef)
JuMP.index(::VariableRef)
JuMP.num_variables(::Model)
JuMP.name(::VariableRef)
JuMP.variable_by_name(::Model, ::String)
all_variables(::Model)
JuMP.has_lower_bound(::VariableRef)
JuMP.lower_bound(::VariableRef)
JuMP.set_lower_bound(::VariableRef, ::Number)
JuMP.LowerBoundRef(::VariableRef)
JuMP.delete_lower_bound(::VariableRef)
JuMP.has_upper_bound(::VariableRef)
JuMP.upper_bound(::VariableRef)
JuMP.set_upper_bound(::VariableRef, ::Number)
JuMP.UpperBoundRef(::VariableRef)
JuMP.delete_upper_bound(::VariableRef)
JuMP.is_fixed(::VariableRef)
JuMP.fix_value(::VariableRef)
JuMP.fix(::VariableRef, ::Number)
JuMP.FixRef(::VariableRef)
JuMP.unfix(::VariableRef)
JuMP.start_value(::VariableRef)
JuMP.set_start_value(::VariableRef, ::Number)
JuMP.is_binary(::VariableRef)
JuMP.set_binary(::VariableRef)
JuMP.BinaryRef(::VariableRef)
JuMP.unset_binary(::VariableRef)
JuMP.is_integer(::VariableRef)
JuMP.set_integer(::VariableRef)
JuMP.IntegerRef(::VariableRef)
JuMP.unset_integer(::VariableRef)
```
