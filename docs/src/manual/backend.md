# [Transformation Backends](@id opt_manual)
A technical manual for optimizing (solving) `InfiniteOpt` models via 
transformation backends. See the respective [guide](@ref opt_docs) 
for more information.

## Optimize
```@docs
JuMP.optimize!(::InfiniteModel)
```

## Backend Settings/Queries
```@docs
transformation_model(::InfiniteModel)
transformation_data(::InfiniteModel)
transformation_backend
set_transformation_backend
build_transformation_backend!(::InfiniteModel)
transformation_variable(::GeneralVariableRef)
transformation_expression(::JuMP.AbstractJuMPScalar)
transformation_constraint(::InfOptConstraintRef)
supports(::Union{DecisionVariableRef, MeasureRef, ParameterFunctionRef})
supports(::JuMP.AbstractJuMPScalar)
supports(::InfOptConstraintRef)
JuMP.set_optimizer(::InfiniteModel, ::Any)
JuMP.set_silent(::InfiniteModel)
JuMP.unset_silent(::InfiniteModel)
JuMP.set_time_limit_sec(::InfiniteModel, ::Real)
JuMP.unset_time_limit_sec(::InfiniteModel)
JuMP.time_limit_sec(::InfiniteModel)
JuMP.solver_name(model::InfiniteModel)
JuMP.mode(::InfiniteModel)
JuMP.compute_conflict!(::InfiniteModel)
JuMP.copy_conflict(::InfiniteModel)
JuMP.bridge_constraints(::InfiniteModel)
JuMP.add_bridge(::InfiniteModel, ::Any)
JuMP.print_active_bridges(::IO, ::InfiniteModel, ::Vararg{Any})
JuMP.print_bridge_graph(::IO, ::InfiniteModel)
JuMP.set_attribute(::InfiniteModel, ::Any, ::Any)
JuMP.set_attributes
JuMP.get_attribute(::InfiniteModel, ::Any)
JuMP.backend(::InfiniteModel)
JuMP.unsafe_backend(::InfiniteModel)
```

## Transformation Backend API 
```@docs
AbstractTransformationBackend
JuMPBackend
AbstractJuMPTag
transformation_model(::AbstractTransformationBackend)
transformation_data(::AbstractTransformationBackend)
JuMP.get_attribute(::AbstractTransformationBackend, ::Any)
JuMP.set_attribute(::AbstractTransformationBackend, ::Any, ::Any)
Base.empty!(::AbstractTransformationBackend)
build_transformation_backend!(::InfiniteModel, ::AbstractTransformationBackend)
JuMP.optimize!(::AbstractTransformationBackend)
JuMP.set_optimizer(::AbstractTransformationBackend, ::Any)
JuMP.mode(::AbstractTransformationBackend)
JuMP.compute_conflict!(::AbstractTransformationBackend)
JuMP.copy_conflict(::AbstractTransformationBackend)
JuMP.bridge_constraints(::AbstractTransformationBackend)
JuMP.add_bridge(::AbstractTransformationBackend, ::Any)
JuMP.print_active_bridges(::IO, ::AbstractTransformationBackend, ::Vararg{Any})
JuMP.print_bridge_graph(::IO, ::AbstractTransformationBackend)
JuMP.variable_ref_type(::AbstractTransformationBackend)
JuMP.backend(::AbstractTransformationBackend)
JuMP.unsafe_backend(::AbstractTransformationBackend)
transformation_variable(::GeneralVariableRef, ::AbstractTransformationBackend)
transformation_expression(::Any, ::AbstractTransformationBackend)
transformation_constraint(::InfOptConstraintRef, ::AbstractTransformationBackend)
variable_supports(::Any, ::AbstractTransformationBackend)
expression_supports(::Any, ::AbstractTransformationBackend)
constraint_supports(::InfOptConstraintRef, ::AbstractTransformationBackend)
transformation_backend_ready
set_transformation_backend_ready
update_parameter_value(::AbstractTransformationBackend, ::Any, ::Any)
warmstart_backend_start_values(::AbstractTransformationBackend)
map_value_to_start
```
