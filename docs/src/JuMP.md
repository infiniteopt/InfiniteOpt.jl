# JuMP Documentation

```@docs
JuMP.Model
JuMP.Model()
JuMP.Model(::OptimizerFactory)
JuMP.with_optimizer
JuMP.fix(::JuMP.VariableRef, ::Number)
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
JuMP.@variable
```
