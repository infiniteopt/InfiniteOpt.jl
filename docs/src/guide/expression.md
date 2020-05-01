# [Expressions] (@id expr_page)
A guide for the defining and understanding the variable expressions
used in `InfiniteOpt`.

## Overview
Expressions in `InfiniteOpt` (also called functions) refer to mathematical
statements involving variables and numbers. Thus, these comprise the
mathematical expressions used that are used in measures, objectives, and
constraints. Programmatically, `InfiniteOpt` simply extends `JuMP` expression
types and methods principally pertaining to affine and quadratic mathematical
expressions. An natively supported abstraction for general nonlinear expressions
is currently under development since that of `JuMP` is not readily extendable.

## Variable Hierarchy
Expressions employ variable reference types inherited from
[`JuMP.AbstractVariableRef`](@ref) to form expression objects. `InfiniteOpt`
uses a hierarchy of such types to organize the complexities associated with
modeling infinite dimensional programs. The figure below summarizes this
hierarchy of variable reference types where the abstract types are depicted in
green and the concrete types are shown blue.

![tree](../assets/variable_tree.png)

Following `JuMP`, expression objects are parameterized by the variable reference
type that is present in the expression. In `InfiniteOpt` expressions
automatically, select the most specific variable reference type possible in
accordance with the above figure. For instance, an expression that only contains
hold variables will be classified as a [`HoldVariableRef`](@ref) expression object,
whereas an expression containing hold variables and a measure would be classified
as a [`MeasureFiniteVariableRef`](@ref) expression object. This hierarchical
classification becomes convenient to guide infinite program reformulation
schemes in how to treat different expressions. The default transcription
methodology employed by `InfiniteOpt.TranscriptionOpt` uses these classifications
to efficiently differentiate between finite and infinite expressions.

## Affine Expressions
An affine expression pertains to a mathematical function of the form:
```math
f_a(x) = a_1x_1 + ... + a_nx_n + b
```
where ``x \in \mathbb{R}^n`` denote variables, ``a \in \mathbb{R}^n`` denote
coefficients, and ``b \in \mathbb{R}`` denotes a constant value. Such
expressions, are prevalent in any problem than involves linear constraints
and/or objectives.

In `InfiniteOpt`, affine expressions can be defined directly
using `Julia`'s arithmetic operators (i.e., `+`, `-`, `*`, etc.) or using
[`@expression`](@ref).  For example, let's define the expression
``2y(t) + z - 3t`` noting that the following methods are equivalent:
```jldoctest affine; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_variable(model, y(t))
y(t)

julia> @hold_variable(model, z)
z

julia> expr = 2y + z - 3t
2 y(t) + z - 3 t

julia> expr = 2 * y + z - 3 * t
2 y(t) + z - 3 t

julia> expr = @expression(model, 2y + z - 3t)
2 y(t) + z - 3 t

julia> typeof(expr)
GenericAffExpr{Float64,GeneralVariableRef}
```
Notice that coefficients to variables can simply be put alongside variables
without having to use the `*` operator. Also, note that all of these expressions
are stored in a container referred to as a `GenericAffExpr` which is a `JuMP`
object for storing affine expressions. Furthermore, this object is parameterized
by [`GeneralVariableRef`](@ref) since it is the lowest common variable reference
type in common between hold variables and infinite variables.

!!! note
    Where possible, it is preferable to use [`@expression`](@ref) for defining
    expressions as it is much more efficient than explicitly using the
    standard operators.

`GenericAffExpr` objects contain 2 fields which are:
- `constant::CoefType` The constant value of the affine expression.
- `terms::OrderDict{VarType, CoefType}` A dictionary mapping variables to coefficients.
For example, let's see what these fields look like in the above example:
```jldoctest affine
julia> expr.terms
OrderedCollections.OrderedDict{GeneralVariableRef,Float64} with 3 entries:
  y(t) => 2.0
  z    => 1.0
  t    => -3.0

julia> expr.constant
0.0
```
Notice that the ordered dictionary preserves the order in which the variables
appear in the expression.

More information can be found in the documentation for affine expressions in
[`JuMP`](http://www.juliaopt.org/JuMP.jl/stable/expressions/#Affine-expressions-1).

## Quadratic Expressions
A quadratic function pertains to a mathematical function of the form:
```math
f_q(x) = a_1x_1^2 + a_2 x_1 x_2 + ... + a_m x_n^2 + f_a(x)
```
where ``x \in \mathbb{R}^n`` are the variables,
``f_a(x): \mathbb{R}^n \mapsto \mathbb{R}`` is an affine function, and
``m = n(n+1)/2`` is the number of unique combinations of variables ``x``.
Like affine expressions, quadratic expressions can be defined via `Julia`'s
arithmetic operators or via [`@expression`](@ref). For example, let's define
``2y^2(t) - zy(t) + 42t - 3`` using the following equivalent methods:
```jldoctest affine
julia> expr = 2y^2 - z * y + 42t - 3
2 y(t)² - z*y(t) + 42 t - 3

julia> expr = @expression(model, 2y^2 - z * y + 42t - 3)
2 y(t)² - z*y(t) + 42 t - 3

julia> typeof(expr)
GenericQuadExpr{Float64,GeneralVariableRef}
```
Again, notice that coefficients need not employ `*`. Also, the object used to
store the expression is a `GenericQuadExpr` which is a `JuMP` object used for
storing quadratic expressions. Again, this expression container is parameterized
by [`GeneralVariableRef`](@ref) since that is the common variable reference type
between the hold variable `z` and the infinite variable `y(t)`.

`GenericQuadExpr` object contains 2 data fields which are:
- `aff::GenericAffExpr{CoefType,VarType}` An affine expression
- `terms::OrderedDict{UnorderedPair{VarType}, CoefType}` A dictionary mapping quadratic variable pairs to coefficients.
Here the `UnorderedPair` type is unique to `JuMP` and contains the fields:
- `a::AbstractVariableRef` One variable in a quadratic pair
- `b::AbstractVariableRef` The other variable in a quadratic pair.
Thus, this form can be used to store arbitrary quadratic expressions. For
example, let's look at what these fields look like in the above example:
```jldoctest affine
julia> expr.aff
42 t - 3

julia> typeof(expr.aff)
GenericAffExpr{Float64,GeneralVariableRef}

julia> expr.terms
OrderedCollections.OrderedDict{UnorderedPair{GeneralVariableRef},Float64} with 2 entries:
  UnorderedPair{GeneralVariableRef}(y(t), y(t)) => 2.0
  UnorderedPair{GeneralVariableRef}(z, y(t))    => -1.0
```
Notice again that the ordered dictionary preserves the order.

!!! tip
    Polynomial expressions can be represented by introducing dumby variables
    and nested quadratic/affine expressions. For instance, ``z^3 + 2`` can be
    expressed by introducing a dumby variable ``x = z^2``:
    ```jldoctest affine
    julia> @hold_variable(model, x)
    x

    julia> @constraint(model, x == z^2)
    -z² + x = 0.0

    julia> expr = @expression(model, z * x + 2)
    z*x + 2
    ```

More information can be found in the documentation for quadratic expressions in
[`JuMP`](http://www.juliaopt.org/JuMP.jl/stable/expressions/#Quadratic-expressions-1).

## Nonlinear Expressions
General nonlinear expressions as generated via `@NLexpression` and similar
methods in `JuMP` are not yet extended for `InfiniteOpt`. This is because
`JuMP` no longer readily supports nonlinear extensions, but a native nonlinear
implementation is currently under development and should be released in the near
future.

## DataTypes
```@index
Pages   = ["expression.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
GeneralVariableRef
DispatchVariableRef
MeasureFiniteVariableRef
FiniteVariableRef
```

## User Methods
```@docs
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
used_by_measure(::GeneralVariableRef)
used_by_objective(::GeneralVariableRef)
used_by_constraint(::GeneralVariableRef)
is_used(::GeneralVariableRef)
JuMP.value(::GeneralVariableRef)
JuMP.set_value(::GeneralVariableRef, ::Real)
infinite_set(::GeneralVariableRef)
infinite_set(::AbstractArray{<:GeneralVariableRef})
set_infinite_set(::GeneralVariableRef, ::InfiniteScalarSet)
set_infinite_set(::AbstractArray{<:GeneralVariableRef}, ::InfiniteArraySet)
num_supports(::GeneralVariableRef)
num_supports(::AbstractArray{<:GeneralVariableRef})
has_supports(::GeneralVariableRef)
has_supports(::AbstractArray{<:GeneralVariableRef})
supports(::GeneralVariableRef)
supports(::AbstractArray{<:GeneralVariableRef})
set_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
set_supports(::AbstractArray{<:GeneralVariableRef},::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}})
add_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
add_supports(::AbstractArray{<:GeneralVariableRef},::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}})
delete_supports(::GeneralVariableRef)
delete_supports(::AbstractArray{<:GeneralVariableRef})
fill_in_supports!(::GeneralVariableRef)
fill_in_supports!(::AbstractArray{<:GeneralVariableRef})
raw_parameter_refs(::GeneralVariableRef)
parameter_refs(::GeneralVariableRef)
parameter_list(::GeneralVariableRef)
set_parameter_refs(::GeneralVariableRef, ::Tuple)
add_parameter_ref(::GeneralVariableRef, ::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}})
infinite_variable_ref(::GeneralVariableRef)
eval_supports(::GeneralVariableRef)
raw_parameter_values(::GeneralVariableRef)
parameter_values(::GeneralVariableRef)
parameter_bounds(::GeneralVariableRef)
has_parameter_bounds(::GeneralVariableRef)
set_parameter_bounds(::GeneralVariableRef,::ParameterBounds{GeneralVariableRef})
add_parameter_bounds(::GeneralVariableRef,::ParameterBounds{GeneralVariableRef})
delete_parameter_bounds(::GeneralVariableRef)
measure_function(::GeneralVariableRef)
measure_data(::GeneralVariableRef)
JuMP.delete(::InfiniteModel, ::GeneralVariableRef)
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
```

## Developer Internal Methods
```@index
Pages   = ["expression.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
_add_data_object
_data_dictionary
_data_object
_delete_data_object
_core_variable_object
_core_variable_object(::GeneralVariableRef)
_set_core_variable_object
_infinite_variable_dependencies
_infinite_variable_dependencies(::GeneralVariableRef)
_reduced_variable_dependencies
_reduced_variable_dependencies(::GeneralVariableRef)
_point_variable_dependencies
_point_variable_dependencies(::GeneralVariableRef)
_measure_dependencies
_measure_dependencies(::GeneralVariableRef)
_constraint_dependencies
_constraint_dependencies(::GeneralVariableRef)
_parameter_number
_parameter_number(::GeneralVariableRef)
_object_number
_object_number(::GeneralVariableRef)
```
