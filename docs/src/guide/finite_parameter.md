# Finite Parameters
A guide and manual to the definition and use of finite parameters in
`InfiniteOpt`. The Datatypes and Methods sections at the end comprise the manual,
and the above sections comprise the guide.  

## Overview
Often a mathematical model needs to be optimized several times in accordance
with a set of fixed parameter values. In such cases, it is typically preferable
to modify these values in place without having to redefine the entire model. This
ability is provided in `InfiniteOpt` via [`@finite_parameter`](@ref) which
permits users to define finite parameters whose values can later be modified
as needed. Furthermore, at the optimization step these parameters are replaced
with their numeric values. Thus, not adding unnecessary decision variables as is
typically done in `JuMP` models using
[`JuMP.fix`](@ref JuMP.fix(::JuMP.VariableRef, ::Number)) on placeholder
variables.  

## Basic Usage
Once an `InfiniteModel` `model` has been defined we can add a finite parameter
via [`@finite_parameter`](@ref). For example, let's define a maximum cost
parameter called `max_cost` with an initial value of `42`:
```jldoctest fpar; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @finite_parameter(model, max_cost, 42)
max_cost
```
Notice that a `Julia` variable called `max_cost` is returned that contains a
[`ParameterRef`](@ref) that points to the finite parameter we have just created.
An array of parameters can also be defined following standard `JuMP` macro syntax:
```jldoctest fpar
julia> values = [2, 3.2, 1];

julia> @finite_parameter(model, params[i = 1:3], values[i])
3-element Array{ParameterRef,1}:
 params[1]
 params[2]
 params[3]
```
The `@finite_parameter` macro emulates all typical `JuMP` functionality and can
define anonymous parameters, use `JuMP` containers and more. We refer to its
documentation below to learn more. Once a finite parameter is defined the
corresponding `ParameterRef` can be used in expressions, objectives, measures,
and constraints just like infinite parameters.

The value of a finite parameter can be checked using
[`JuMP.value`](@ref JuMP.value(::ParameterRef)) and can modified using
[`JuMP.set_value`](@ref JuMP.set_value(::ParameterRef, ::Number)). For example,
let's update the value of `max_cost` to be now be `10.2`:
```jldoctest fpar
julia> value(max_cost)
42

julia> set_value(max_cost, 10.2)

julia> value(max_cost)
10.2
```

## Advanced Details
The ability to implement finite parameters simply stems from `InfiniteOpt`'s
framework for infinite parameters. In reality finite parameters are simply
infinite parameters defined with an [`IntervalSet`](@ref) whose lower and upper
bounds equal to the desired value. Thus, [`@finite_parameter`](@ref) simply
processes the names, values, and other arguments and then makes the appropriate
[`@infinite_parameter`](@ref) call. For example, the following calls are
equivalent:
```julia
julia> @finite_parameter(model, param[1:2], 42)
2-element Array{ParameterRef,1}:
 param[1]
 param[2]

julia> @infinite_parameter(model, param[1:2] in [42, 42], supports = [42])
2-element Array{ParameterRef,1}:
 param[1]
 param[2]
```
Furthermore, the [`JuMP.value`](@ref JuMP.value(::ParameterRef)) checks the
bound value (also ensuring the bounds match) and
[`JuMP.set_value`](@ref JuMP.set_value(::ParameterRef, ::Number)) updates the
the `IntervalSet` and updates the support value.

## Methods/Macros
```@index
Pages   = ["finite_parameter.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@finite_parameter
is_finite_parameter
JuMP.value(::ParameterRef)
JuMP.set_value(::ParameterRef, ::Number)
```
