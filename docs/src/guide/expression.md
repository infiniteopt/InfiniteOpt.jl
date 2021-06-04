```@meta
DocTestFilters = [r"≤|<=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI",
                  r" for all | ∀ ", r"d|∂", r"integral|∫"]
```

# [Expressions] (@id expr_page)
A guide for the defining and understanding the variable expressions
used in `InfiniteOpt`.

!!! note 
    Nonlinear objects as defined by `JuMP.@NL[macro_name]` are not currently 
    supported by `InfiniteOpt`. See [Nonlinear Expressions](@ref) for more 
    information and possible workarounds. 

## Overview
Expressions in `InfiniteOpt` (also called functions) refer to mathematical
statements involving variables and numbers. Thus, these comprise the
mathematical expressions used that are used in measures, objectives, and
constraints. Programmatically, `InfiniteOpt` simply extends `JuMP` expression
types and methods principally pertaining to affine and quadratic mathematical
expressions. A natively supported abstraction for general nonlinear expressions
is planned for development since that of `JuMP` is not readily extendable.

## Parameter Functions
As described further below, InfiniteOpt.jl only supports affine and quadratic 
expressions in its current rendition. However, there several use cases where we 
might want to provide a more complex known function of infinite parameter(s) (e.g., 
nonlinear setpoint tracking). Thus, we provide parameter function objects 
that given a particular realization of infinite parameters will output a scalar 
value. Note that this can be interpreted as an infinite variable that is 
constrained to a particular known function. This is accomplished via 
[`@parameter_function`](@ref) or [`parameter_function`](@ref) and is exemplified 
by defining a parameter function `f(t)` that uses `sin(t)`:
```jldoctest param_func
julia> using InfiniteOpt;

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @parameter_function(model, f == sin(t))
f(t)
```
Here we created an parameter function object, added it to `model`, and 
then created a Julia variable `f` that serves as a `GeneralVariableRef` that points 
to it. From here we can treat `f` as a normal infinite variable and use it with 
measures, derivatives, and constraints. For example, we can do the following:
```jldoctest param_func
julia> @variable(model, y, Infinite(t));

julia> df = deriv(f, t)
∂/∂t[f(t)]

julia> meas = integral(y - f, t)
∫{t ∈ [0, 10]}[y(t) - f(t)]

julia> @constraint(model, y - f <= 0)
y(t) - f(t) ≤ 0.0, ∀ t ∈ [0, 10]
```
We can also define parameter functions that depend on multiple infinite 
parameters even use an anonymous function if prefer:
```jldoctest param_func
julia> @infinite_parameter(model, x[1:2] in [-1, 1]);

julia> @parameter_function(model, myname == (t, x) -> t + sum(x))
myname(t, x)
```
In many applications, we may also desire to define an array of parameter functions 
that each use a different realization of some parent function by varying some 
additional positional/keyword arguments. We readily support this behavior since 
parameter functions can be defined with additional known arguments:
```jldoctest param_func
julia> @parameter_function(model, pfunc_alt[i = 1:3] == t -> mysin(t, as[i], b = 0))
3-element Array{GeneralVariableRef,1}:
 pfunc_alt[1](t)
 pfunc_alt[2](t)
 pfunc_alt[3](t)
```
The main recommended use case for [`parameter_function`](@ref) is that it is 
amendable to defining complex anonymous functions via a do-block which is useful 
for applications like defining a time-varied setpoint:
```jldoctest param_func
julia> setpoint = parameter_function(t, name = "setpoint") do t_supp
                    if t_supp <= 5
                        return 2.0
                    else 
                        return 10.2
                    end
                 end
setpoint(t)
```
Please consult the following links for more information about defining parameter 
functions: [`@parameter_function`](@ref) and [`parameter_function`](@ref).

Beyond this, there are number of query and modification methods that can be 
employed for parameter functions and these are detailed in the 
[Parameter Function Methods](@ref) Section below.

## Variable Hierarchy
Expressions employ variable reference types inherited from
`JuMP.AbstractVariableRef` to form expression objects. `InfiniteOpt`
uses a hierarchy of such types to organize the complexities associated with
modeling infinite dimensional programs. The figure below summarizes this
hierarchy of variable reference types where the abstract types are depicted in
green and the concrete types are shown blue.

![tree](../assets/variable_tree.png)

In consistently with `JuMP` expression support, [`GeneralVariableRef`](@ref)
exists as a variable reference type that is able to represent any of the above
concrete subtypes of [`DispatchVariableRef`](@ref). This allows the expression
containers to be homogeneous in variable type. This is a paradigm shift from
previous versions of `InfiniteOpt` that used the hierarchy of types directly
to construct expressions. This behavior led to stability and performance
limitations and thus a has been discontinued.

However, the variable hierarchy is still used to create for variable methods.
To accomplish this appropriate `GeneralVariableRef` dispatch methods are implemented
(which are detailed in User Methods section at the bottom of this page) that
utilize [`dispatch_variable_ref`](@ref) to create the appropriate concrete
subtype of `DispatchVariableRef` and call the appropriate underlying method.
These dispatch methods have been implemented for all public methods and the
underlying methods are what are documented in the method manuals throughout the
User Guide pages.

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
`@expression`.  For example, let's define the expression
``2y(t) + z - 3t`` noting that the following methods are equivalent:
```jldoctest affine; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @variable(model, y, Infinite(t))
y(t)

julia> @variable(model, z)
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
object for storing affine expressions.

!!! note
    Where possible, it is preferable to use 
    [`@expression`](https://jump.dev/JuMP.jl/v0.21.8/reference/expressions/#JuMP.@expression) 
    for defining expressions as it is much more efficient than explicitly using 
    the standard operators.

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
[`JuMP`](https://jump.dev/JuMP.jl/v0.21.8/reference/expressions/#Affine-expressions).

## Quadratic Expressions
A quadratic function pertains to a mathematical function of the form:
```math
f_q(x) = a_1x_1^2 + a_2 x_1 x_2 + ... + a_m x_n^2 + f_a(x)
```
where ``x \in \mathbb{R}^n`` are the variables,
``f_a(x): \mathbb{R}^n \mapsto \mathbb{R}`` is an affine function, and
``m = n(n+1)/2`` is the number of unique combinations of variables ``x``.
Like affine expressions, quadratic expressions can be defined via `Julia`'s
arithmetic operators or via `@expression`. For example, let's define
``2y^2(t) - zy(t) + 42t - 3`` using the following equivalent methods:
```jldoctest affine
julia> expr = 2y^2 - z * y + 42t - 3
2 y(t)² - z*y(t) + 42 t - 3

julia> expr = @expression(model, 2y^2 - z * y + 42t - 3)
2 y(t)² - y(t)*z + 42 t - 3

julia> typeof(expr)
GenericQuadExpr{Float64,GeneralVariableRef}
```
Again, notice that coefficients need not employ `*`. Also, the object used to
store the expression is a `GenericQuadExpr` which is a `JuMP` object used for
storing quadratic expressions.

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
  UnorderedPair{GeneralVariableRef}(y(t), z)    => -1.0
```
Notice again that the ordered dictionary preserves the order.

!!! tip
    Polynomial expressions can be represented by introducing dumby variables
    and nested quadratic/affine expressions. For instance, ``z^3 + 2`` can be
    expressed by introducing a dumby variable ``x = z^2``:
    ```jldoctest affine
    julia> @variable(model, x)
    x

    julia> @constraint(model, x == z^2)
    -z² + x = 0.0

    julia> expr = @expression(model, z * x + 2)
    z*x + 2
    ```

More information can be found in the documentation for quadratic expressions in
[`JuMP`](https://jump.dev/JuMP.jl/v0.21.8/reference/expressions/#Quadratic-expressions).

## Nonlinear Expressions
General nonlinear expressions as generated via `JuMP.@NLexpression`, 
`JuMP.@NLobjective`, and/or `JuMP.@NLconstraint` macros in `JuMP` are not 
currently for `InfiniteOpt`. This is because `JuMP` does not readily support 
nonlinear extensions. However, a fundamental overhaul is planned to resolve this 
problem (check the status on 
[GitHub](https://github.com/jump-dev/MathOptInterface.jl/issues/846)).

### Workarounds
In the meantime, this limitation can often be overcome by reformulating the 
problem formulation. 

One common case involves expressions that entail integer powers that are greater 
than 2. This can readily be remedied by adding placeholder variables. For example, 
consider the expression ``z^2x - 3z``. We can reformulate by introducing 
``z' = z^2``:
```jldoctest affine
julia> @variable(model, z_squared)
z_squared

julia> @constraint(model, z_squared == z^2)
-z² + z_squared = 0.0

julia> @expression(model, z_squared * x - 3z)
z_squared*x - 3 z
```

We can also reformulate for a variety of nonlinear function types:

| Function     | Example                        | Reformulation Method                                                                                                               |
|:------------:|:------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------:|
| Square Root  | ``\sqrt{z}``                   | Make a squared reformulation variable                                                                                              |
| Indicator    | ``\mathbb{1}_{z \geq \alpha}`` | Use [`JuMP`'s indicator constraint syntax](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/#Indicator-constraints)             | 
| Indicator    | ``\mathbb{1}_{z \geq \alpha}`` | Replace with big-M constraints ([reference](https://www.gurobi.com/documentation/9.1/refman/dealing_with_big_m_constra.html))      | 
| Max/Min      | ``\max(z, a)``                 | Linear programming cuts or big-m constraints ([reference](https://or.stackexchange.com/questions/711/how-to-formulate-linearize-a-maximum-function-in-a-constraint/712#712)) |

In other cases, it may be possible to use a formulation that uses vector 
constraints. For example, it might be possible to model your problem using 
semi-definite and/or conic constraints. 

Also note that any nonlinearites that only involve infinite parameters (i.e., 
no decision variables) are enabled via parameter functions. See 
[Parameter Functions](@ref) for more information.

For problems that cannot be readily reformulated, `JuMP` can be used directly. In 
this case the user will need to first transform the formulation into a finite 
representation (e.g., discretize it). 

## DataTypes
```@index
Pages   = ["expression.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
GeneralVariableRef
DispatchVariableRef
FiniteRef
ParameterFunctionRef
ParameterFunctionIndex
ParameterFunction
ParameterFunctionData
```

## Parameter Function Methods 
```@docs
@parameter_function
parameter_function
build_parameter_function
add_parameter_function
JuMP.name(::ParameterFunctionRef)
JuMP.set_name(::ParameterFunctionRef, ::String)
parameter_refs(::ParameterFunctionRef)
raw_parameter_refs(::ParameterFunctionRef)
parameter_list(::ParameterFunctionRef)
raw_function(::ParameterFunctionRef)
call_function
used_by_semi_infinite_variable(::ParameterFunctionRef)
used_by_derivative(::ParameterFunctionRef)
used_by_measure(::ParameterFunctionRef)
used_by_constraint(::ParameterFunctionRef)
is_used(::ParameterFunctionRef)
JuMP.delete(::InfiniteModel, ::ParameterFunctionRef)
```

## Expression Methods
```@docs
parameter_refs(::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr})
```

## GeneralVariableRef User Methods
```@index
Pages   = ["expression.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
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
used_by_semi_infinite_variable(::GeneralVariableRef)
used_by_derivative(::GeneralVariableRef)
used_by_measure(::GeneralVariableRef)
used_by_objective(::GeneralVariableRef)
used_by_constraint(::GeneralVariableRef)
is_used(::GeneralVariableRef)
has_derivative_constraints(::GeneralVariableRef)
parameter_value(::GeneralVariableRef)
JuMP.set_value(::GeneralVariableRef, ::Real)
infinite_domain(::GeneralVariableRef)
infinite_domain(::AbstractArray{<:GeneralVariableRef})
set_infinite_domain(::GeneralVariableRef, ::InfiniteScalarDomain)
set_infinite_domain(::AbstractArray{<:GeneralVariableRef}, ::InfiniteArrayDomain)
num_supports(::GeneralVariableRef)
num_supports(::AbstractArray{<:GeneralVariableRef})
has_supports(::GeneralVariableRef)
has_supports(::AbstractArray{<:GeneralVariableRef})
supports(::GeneralVariableRef)
supports(::AbstractArray{<:GeneralVariableRef})
set_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
set_supports(::AbstractArray{<:GeneralVariableRef},::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}})
add_supports(::GeneralVariableRef,::Union{Real, Vector{<:Real}})
add_supports(::AbstractArray{<:GeneralVariableRef},::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}})
delete_supports(::GeneralVariableRef)
delete_supports(::AbstractArray{<:GeneralVariableRef})
fill_in_supports!(::GeneralVariableRef)
fill_in_supports!(::AbstractArray{<:GeneralVariableRef})
raw_parameter_refs(::GeneralVariableRef)
parameter_refs(::GeneralVariableRef)
parameter_list(::GeneralVariableRef)
raw_function(::GeneralVariableRef)
infinite_variable_ref(::GeneralVariableRef)
eval_supports(::GeneralVariableRef)
raw_parameter_values(::GeneralVariableRef)
parameter_values(::GeneralVariableRef)
significant_digits(::GeneralVariableRef)
measure_function(::GeneralVariableRef)
measure_data(::GeneralVariableRef)
is_analytic(::GeneralVariableRef)
derivative_argument(::GeneralVariableRef)
operator_parameter(::GeneralVariableRef)
derivative_method(::GeneralVariableRef)
evaluate(::GeneralVariableRef)
derivative_constraints(::GeneralVariableRef)
delete_derivative_constraints(::GeneralVariableRef)
InfiniteOpt.add_generative_supports(::GeneralVariableRef)
set_derivative_method(::GeneralVariableRef, ::AbstractDerivativeMethod)
has_generative_supports(::GeneralVariableRef)
has_internal_supports(::GeneralVariableRef)
JuMP.delete(::InfiniteModel, ::GeneralVariableRef)
JuMP.delete(::InfiniteModel,::AbstractArray{<:GeneralVariableRef})
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
start_value_function(::GeneralVariableRef)
set_start_value_function(::GeneralVariableRef, ::Any)
reset_start_value_function(::GeneralVariableRef)
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
```@docs
InfiniteOpt._add_data_object
InfiniteOpt._data_dictionary
InfiniteOpt._data_object
InfiniteOpt._delete_data_object
InfiniteOpt._core_variable_object
InfiniteOpt._core_variable_object(::GeneralVariableRef)
InfiniteOpt._set_core_variable_object
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
InfiniteOpt._object_number
InfiniteOpt._object_number(::GeneralVariableRef)
```
