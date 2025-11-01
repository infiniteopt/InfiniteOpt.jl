# [Expressions](@id expr_manual)
A technical manual for variable expressions in `InfiniteOpt`. See the respective 
[guide](@ref expr_docs) for more information.

## [Parameter Functions](@id par_func_manual)
### Definition
```@docs
@parameter_function
parameter_function
build_parameter_function
ParameterFunction
add_parameter_function
ParameterFunctionData
ParameterFunctionIndex
ParameterFunctionRef
```

### Queries
```@docs
JuMP.name(::ParameterFunctionRef)
raw_function(::ParameterFunctionRef)
JuMP.parameter_value(::ParameterFunctionRef)
call_function
parameter_refs(::ParameterFunctionRef)
parameter_list(::ParameterFunctionRef)
raw_parameter_refs(::ParameterFunctionRef)
is_used(::ParameterFunctionRef)
used_by_semi_infinite_variable(::ParameterFunctionRef)
used_by_point_variable(::ParameterFunctionRef)
used_by_measure(::ParameterFunctionRef)
used_by_constraint(::ParameterFunctionRef)
parameter_group_int_indices(::ParameterFunctionRef)
num_parameter_functions
all_parameter_functions
```

### Modification
```@docs
JuMP.set_name(::ParameterFunctionRef, ::String)
JuMP.set_parameter_value(::ParameterFunctionRef, ::Function)
JuMP.delete(::InfiniteModel, ::ParameterFunctionRef)
```

## [Nonlinear Expressions](@id nlp_manual)
### DataTypes
```@docs
NLPOperator
```

### Methods
```@docs
JuMP.add_nonlinear_operator
all_nonlinear_operators
name_to_operator
added_nonlinear_operators
add_operators_to_jump
```

## Expression Methods
```@docs
parameter_refs(::Union{JuMP.GenericAffExpr, JuMP.GenericNonlinearExpr, JuMP.GenericQuadExpr})
restrict(::JuMP.AbstractJuMPScalar)
map_expression
map_expression_to_ast
all_expression_variables
parameter_group_int_indices(::Any)
```

## GeneralVariableRef User Methods
```@docs
GeneralVariableRef
DispatchVariableRef
FiniteRef
JuMP.owner_model(::GeneralVariableRef)
JuMP.owner_model(::DispatchVariableRef)
JuMP.index(::GeneralVariableRef)
JuMP.index(::DispatchVariableRef)
dispatch_variable_ref(::GeneralVariableRef)
dispatch_variable_ref
JuMP.name(::GeneralVariableRef)
JuMP.set_name(::GeneralVariableRef, ::String)
JuMP.is_valid(::InfiniteModel,::GeneralVariableRef)
JuMP.is_valid(::InfiniteModel, ::DispatchVariableRef)
used_by_infinite_variable(::GeneralVariableRef)
used_by_point_variable(::GeneralVariableRef)
used_by_semi_infinite_variable(::GeneralVariableRef)
used_by_parameter_function(::GeneralVariableRef)
used_by_derivative(::GeneralVariableRef)
used_by_measure(::GeneralVariableRef)
used_by_objective(::GeneralVariableRef)
used_by_constraint(::GeneralVariableRef)
is_used(::GeneralVariableRef)
has_derivative_constraints(::GeneralVariableRef)
JuMP.parameter_value(::GeneralVariableRef)
JuMP.set_parameter_value(::GeneralVariableRef, ::Any)
infinite_domain(::GeneralVariableRef)
infinite_domain(::Array{<:GeneralVariableRef})
set_infinite_domain(::GeneralVariableRef, ::InfiniteScalarDomain)
set_infinite_domain(::Array{<:GeneralVariableRef}, ::InfiniteArrayDomain)
num_supports(::GeneralVariableRef)
num_supports(::Array{<:GeneralVariableRef})
has_supports(::GeneralVariableRef)
has_supports(::Array{<:GeneralVariableRef})
supports(::GeneralVariableRef)
supports(::Array{<:GeneralVariableRef})
set_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
set_supports(::Array{<:GeneralVariableRef},::Union{Array{<:Real, 2}, Vector{<:Array{<:Real}}})
add_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
add_supports(::Array{<:GeneralVariableRef},::Union{Array{<:Real, 2}, Vector{<:Array{<:Real}}})
delete_supports(::GeneralVariableRef)
delete_supports(::Array{<:GeneralVariableRef})
fill_in_supports!(::GeneralVariableRef)
fill_in_supports!(::Array{<:GeneralVariableRef})
raw_parameter_refs(::GeneralVariableRef)
parameter_refs(::GeneralVariableRef)
parameter_list(::GeneralVariableRef)
raw_function(::GeneralVariableRef)
infinite_variable_ref(::GeneralVariableRef)
eval_support(::GeneralVariableRef)
raw_parameter_values(::GeneralVariableRef)
parameter_values(::GeneralVariableRef)
significant_digits(::GeneralVariableRef)
measure_function(::GeneralVariableRef)
measure_data(::GeneralVariableRef)
is_analytic(::GeneralVariableRef)
derivative_argument(::GeneralVariableRef)
operator_parameter(::GeneralVariableRef)
derivative_order(::GeneralVariableRef)
derivative_method(::GeneralVariableRef)
evaluate(::GeneralVariableRef)
derivative_constraints(::GeneralVariableRef)
delete_derivative_constraints(::GeneralVariableRef)
InfiniteOpt.add_generative_supports(::GeneralVariableRef)
set_derivative_method(::GeneralVariableRef, ::AbstractDerivativeMethod)
has_generative_supports(::GeneralVariableRef)
has_internal_supports(::GeneralVariableRef)
JuMP.delete(::InfiniteModel, ::GeneralVariableRef)
JuMP.delete(::InfiniteModel,::Array{<:GeneralVariableRef})
JuMP.has_lower_bound(::GeneralVariableRef)
JuMP.lower_bound(::GeneralVariableRef)
JuMP.set_lower_bound(::GeneralVariableRef,::Real)
JuMP.LowerBoundRef(::GeneralVariableRef)
JuMP.delete_lower_bound(::GeneralVariableRef)
JuMP.has_upper_bound(::GeneralVariableRef)
JuMP.upper_bound(::GeneralVariableRef)
JuMP.set_upper_bound(::GeneralVariableRef,::Real)
JuMP.UpperBoundRef(::GeneralVariableRef)
JuMP.delete_upper_bound(::GeneralVariableRef)
JuMP.is_fixed(::GeneralVariableRef)
JuMP.fix_value(::GeneralVariableRef)
JuMP.fix(::GeneralVariableRef, ::Real)
JuMP.FixRef(::GeneralVariableRef)
JuMP.unfix(::GeneralVariableRef)
JuMP.start_value(::GeneralVariableRef)
JuMP.set_start_value(::GeneralVariableRef, ::Real)
JuMP.is_binary(::GeneralVariableRef)
JuMP.set_binary(::GeneralVariableRef)
JuMP.BinaryRef(::GeneralVariableRef)
JuMP.unset_binary(::GeneralVariableRef)
JuMP.is_integer(::GeneralVariableRef)
JuMP.set_integer(::GeneralVariableRef)
JuMP.IntegerRef(::GeneralVariableRef)
JuMP.unset_integer(::GeneralVariableRef)
constant_over_collocation(::GeneralVariableRef, ::GeneralVariableRef)
core_object
core_object(::GeneralVariableRef)
parameter_group_int_indices(::GeneralVariableRef)
InfiniteOpt.parameter_group_int_index
InfiniteOpt.parameter_group_int_index(::GeneralVariableRef)
```

## Developer Internal Methods
```@docs
InfiniteOpt._add_data_object
InfiniteOpt._data_dictionary
InfiniteOpt._data_object
InfiniteOpt._delete_data_object
InfiniteOpt._set_core_object
InfiniteOpt._infinite_variable_dependencies
InfiniteOpt._infinite_variable_dependencies(::GeneralVariableRef)
InfiniteOpt._semi_infinite_variable_dependencies
InfiniteOpt._semi_infinite_variable_dependencies(::GeneralVariableRef)
InfiniteOpt._point_variable_dependencies
InfiniteOpt._point_variable_dependencies(::GeneralVariableRef)
InfiniteOpt._derivative_dependencies
InfiniteOpt._derivative_dependencies(::GeneralVariableRef)
InfiniteOpt._measure_dependencies
InfiniteOpt._measure_dependencies(::GeneralVariableRef)
InfiniteOpt._generative_measures
InfiniteOpt._generative_measures(::GeneralVariableRef)
InfiniteOpt._constraint_dependencies
InfiniteOpt._constraint_dependencies(::GeneralVariableRef)
InfiniteOpt._derivative_constraint_dependencies
InfiniteOpt._derivative_constraint_dependencies(::GeneralVariableRef)
InfiniteOpt._parameter_number
InfiniteOpt._parameter_number(::GeneralVariableRef)
```
