# [Finite Parameters] (@id finite_param_docs)
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
```jldoctest fpar; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @finite_parameter(model, max_cost, 42)
max_cost
```
Notice that a `Julia` variable called `max_cost` is returned that contains a
[`GeneralVariableRef`](@ref) that points to the finite parameter we have just created.
An array of parameters can also be defined following standard `JuMP` macro syntax:
```jldoctest fpar
julia> values = [2, 3.2, 1];

julia> @finite_parameter(model, params[i = 1:3], values[i])
3-element Array{GeneralVariableRef,1}:
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
[`parameter_value`](@ref parameter_value(::FiniteParameterRef)) and can modified using
[`JuMP.set_value`](@ref JuMP.set_value(::FiniteParameterRef, ::Real)). For example,
let's update the value of `max_cost` to be now be `10.2`:
```jldoctest fpar
julia> parameter_value(max_cost)
42.0

julia> set_value(max_cost, 10.2)

julia> parameter_value(max_cost)
10.2
```

## Advanced Details
The ability to implement finite parameters stems from its ability to support 
mixed variable types using by using the `GeneralVariableRef` buffer. As such, 
finite parameters will be treated as variables until the model is transcribed. 
For example, this means that the expression `max_cost * x` will be treated as a 
quadratic expression when it is expressed in its `InfiniteOpt` form, however it is 
converted into the appropriate affine expression when transcripted. 

!!! note 
    In previous versions finite parameters were just special cases of infinite 
    parameters. However, they now have their own distinct underlying data structure. 

!!! warning 
    `InfiniteOpt`'s implementation of finite parameters should not be a reason to 
    use `InfiniteOpt` to model non-infinite-dimensional problems, since the added 
    overhead will make it slower than just iteratively building `JuMP` models. For 
    this behavior, we recommend looking into using `ParameterJuMP`.

## Datatypes
```@index
Pages   = ["finite_parameter.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
FiniteParameter
FiniteParameterIndex
FiniteParameterRef
```

## Methods/Macros
```@index
Pages   = ["finite_parameter.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@finite_parameter
parameter_value(::FiniteParameterRef)
JuMP.set_value(::FiniteParameterRef, ::Real)
```
