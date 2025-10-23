# [Finite Parameters](@id finite_param_docs)
A guide for the definition and use of finite parameters in `InfiniteOpt`. See 
the respective [technical manual](@ref finite_param_manual) for more details.

## Overview
Often a mathematical model needs to be optimized several times in accordance
with a set of fixed parameter values. In such cases, it is typically preferable
to modify these values in place without having to redefine the entire model. This
ability is provided in `InfiniteOpt` via [`@finite_parameter`](@ref) which
permits users to define finite parameters whose values can later be modified
as needed. Furthermore, at the optimization step these parameters are replaced
with their numeric values. Thus, not adding unnecessary decision variables as is
typically done in `JuMP` models using `JuMP.fix` on placeholder variables.  

!!! warning 
    In some cases, using [`@finite_parameter`](@ref) can unexpectedly make 
    the underlying transformation backend contain nonlinear constraints/objectives. This 
    occurs when a quadratic expression is mutliplied by a finite parameter 
    (making a `GenericNonlinearExpr`):
    ```julia-repl
    julia> model = InfiniteModel(); @variable(model, z); @finite_parameter(model, p == 2);

    julia> @objective(model, Min,  p * z^2) # becomes a nonlinear objective 
    p * (z²)
    ```
    In these cases, a nonlinear solver like `Ipopt` should be used or the 
    finite parameter syntax should be avoided if a quadratic solver like 
    `Gurobi` is needed.

## Basic Usage
Once an `InfiniteModel` `model` has been defined we can add a finite parameter
via [`@finite_parameter`](@ref). For example, let's define a maximum cost
parameter called `max_cost` with an initial value of `42`:
```jldoctest fpar; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @finite_parameter(model, max_cost == 42)
max_cost
```
Notice that a `Julia` variable called `max_cost` is returned that contains a
[`GeneralVariableRef`](@ref) that points to the finite parameter we have just created.
An array of parameters can also be defined following standard `JuMP` macro syntax:
```jldoctest fpar
julia> values = [2, 3.2, 1];

julia> @finite_parameter(model, params[i = 1:3] == values[i])
3-element Vector{GeneralVariableRef}:
 params[1]
 params[2]
 params[3]
```
The `@finite_parameter` macro emulates all typical `JuMP` functionality and can
define anonymous parameters, use `JuMP` containers and more. We refer to its
documentation below to learn more. Once a finite parameter is defined the
corresponding `GeneralVariableRef` can be used in expressions, objectives, measures,
and constraints just like infinite parameters.

The value of a finite parameter can be checked using
[`parameter_value`](@ref JuMP.parameter_value(::FiniteParameterRef)) and can modified using
[`set_parameter_value`](@ref JuMP.set_parameter_value(::FiniteParameterRef, ::Real)). For example,
let's update the value of `max_cost` to be now be `10.2`:
```jldoctest fpar
julia> parameter_value(max_cost)
42.0

julia> set_parameter_value(max_cost, 10.2)

julia> parameter_value(max_cost)
10.2
```

## Advanced Details
The ability to implement finite parameters stems from its ability to support 
mixed variable types using by using the `GeneralVariableRef` buffer. As such, 
finite parameters will be treated as variables until the model is transcribed. 
For example, this means that the expression `max_cost * x` will be treated as a 
quadratic expression when it is expressed in its `InfiniteOpt` form, however it is 
converted into the appropriate affine expression when transcribed. 
