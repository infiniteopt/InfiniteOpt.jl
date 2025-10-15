```@meta
DocTestFilters = [r"≤|<=", r" == | = ", r" ∈ | in ", r" for all | ∀ ", r"sin \(.*"]
```

# [Variables](@id var_docs)
A guide for variables in `InfiniteOpt`. See the respective 
[technical manual](@ref var_manual) for more details.

## Overview
Decision variables are at the very core of `InfiniteOpt` as its name alludes
to mathematical programs that entail infinite decision spaces (i.e., contain
infinite decision variables). Principally, 4 variable types are employed:
infinite, semi-infinite, point, and finite. Infinite variables encompass any 
decision variable that is parameterized by an infinite parameter(s) (e.g., 
space-time variables and stochastic recourse variables). Semi-infinite variables 
denote infinite variables where certain infinite parameters are restricted to 
point values. Point variables are infinite variables at a particular point. 
Finally, finite variables are decisions that are made irrespective of the 
infinite domain (e.g., first stage variables and design variables).

## Basic Usage
Infinite, semi-infinite, point, and finite variables are summarized in the 
following table:

| Name          | Variable Type Object   | Description                            | Example                |
|:-------------:|:----------------------:|:--------------------------------------:|:----------------------:|
| Infinite      | [`Infinite`](@ref)     | decision functions                     | ``y(t, x, \xi)``       |
| Semi-Infinite | [`SemiInfinite`](@ref) | partially evaluated decision functions | ``y(t_0, x, \xi)``     |
| Point         | [`Point`](@ref)        | fully evaluated decision functions     | ``y(t_0, x_0, \xi_k)`` |
| Finite        | NA                     | classical decision variables           | ``z``                  |

Infinite, semi-infinite, point, and finite variables are defined via 
[`@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
(inherited from `JuMP`) with their respective variable type 
object arguments: [`Infinite`](@ref), [`SemiInfinite`](@ref), and [`Point`](@ref) 
(finite variables don't use a variable type object).

Let's first set up a simple space-time model with infinite parameters time `t` and
spatial position `x`:
```jldoctest var_basic
julia> using InfiniteOpt

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_parameter(model, x[1:2] in [-1, 1], independent = true)
2-element Vector{GeneralVariableRef}:
 x[1]
 x[2]
```

### Infinite Variables
Now let's define a time dependent infinite variable `y(t)` with a lower bound of
0:
```jldoctest var_basic
julia> @variable(model, y >= 0, Infinite(t))
y(t)
```
This creates a Julia variable `y` that points to the decision variable `y(t)`
that is stored in `model` which is added to include a lower bound of 0. Another
useful case is that of defining an array of variables `w` that depend on both
position and time:
```jldoctest var_basic
julia> @variable(model, w[i = 1:3], Infinite(t, x...), start = [0, 2, 1][i])
3-element Vector{GeneralVariableRef}:
 w[1](t, x[1], x[2])
 w[2](t, x[1], x[2])
 w[3](t, x[1], x[2])
```
Thus, we create a Julia array variable `w` whose elements `w[i]` point to their
respective infinite variables `w[i](t, x)` stored in `model`. Note that the `i`
used in the array definition can be used to index attributes assigned to each
variable in the array. In this case, we used `i` to assign different initial
guess values for each variable via the `start` keyword argument. Also note that
we use `Infinite(t, x...)` instead of `Infinite(t, x)` since each scalar independent
infinite parameter must be in its own argument.

!!! note
    Previous versions of InfiniteOpt supported arrays of independent infinite
    parameters for infinite variables (e.g., `Infinite(x)` where `x` is an array
    of independent infinite parameters). Now each independent infinite parameter
    must be an individual argument. This can be readily accomplished by splatting
    the infinite parameters (e.g., `Infinite(x...)`).

Moreover, for infinite variables a function can be given to determine the start
values over a range of support points (e.g., a guess trajectory). This is
discussed further below in the Macro Definition section.

### Semi-Infinite Variables
Now let's restrict the above infinite variables `w[i](t, x)` to a particular 
time via semi-infinite variables:
```jldoctest var_basic
julia> @variable(model, w0[i = 1:3], SemiInfinite(w[i], 0, x...))
3-element Vector{GeneralVariableRef}:
 w0[1]
 w0[2]
 w0[3]
```
Thus, we create a Julia array variable `w0` whose elements `w0[i]` point to their
respective semi-infinite variables `w[i](0, x)` stored in `model`. Alternatively, 
we can make a semi-infinite variable via our restriction syntax:
```jldoctest var_basic
julia> [w[i](0, x...) for i in 1:3]
3-element Vector{GeneralVariableRef}:
 w0[1]
 w0[2]
 w0[3]
```
These are often useful to define semi-infinite variables directly in constraint 
expressions. See [Restricted Variables](@ref) to learn about symbolic inline 
definition of semi-infinite variables.

### Point Variables
Now let's add some point variables. These allow us to consider an infinite
variable evaluated at a certain infinite parameter point. For example, let's
define a point variable for `y(0)` with the alias `y0` that is fixed at a value
of 0:
```jldoctest var_basic
julia> @variable(model, y0 == 0, Point(y, 0))
y0
```
Here we create a Julia variable `y0` which points to the point variable `y(0)`.
Notice that in the second argument we specify the infinite variable indexed at
the appropriate parameter value(s). Point variables automatically inherit
attributes of the infinite variable (e.g., bounds, start values, etc.), but
these are overwritten with properties specified for the point variable. In this
case the lower bound inherited from `y(t)` is overwritten by instead fixing
`y(0)` to a value of 0.  

Alternatively, we can use the convenient restriction syntax:
```jldoctest var_basic
julia> y(0)
y0
```
Again this is very useful when embedded directly in constraint expressions 
(e.g., when defining boundary conditions). See [Restricted Variables](@ref) to 
learn about symbolic inline definition of point variables.

### Finite Variables
Finally, we can add finite variables to our model. These denote variables that
hold a single value over the infinite domain or some portion of it (e.g.,
design variables, first stage variables, etc.). Let's add a finite variable
``0 \leq d \leq 42`` that is an integer variable and defined over all infinite
domains (i.e., time and space):
```jldoctest var_basic
julia> @variable(model, 0 <= d <= 42, Int)
d
```
This creates a Julia variable `d` that points to the finite variable `d` which has
a lower bound of 0, an upper bound of 42, and is an integer variable. Thus, finite 
variables are equivalent to those employed in `JuMP`.

Now we have defined variables that we can use in the objective, measures, and
constraints. Please note that the above tutorial only shows a small portion of
the capabilities and options available in defining variables. A full description
is provided in the documentation below.

## Variable Definition Methodology
Defining/initializing a variable (what happens behind the scenes of the variable 
macros) principally involves the following steps:
1. Define the variable information pertaining to `JuMP.VariableInfo` (e.g., 
   bounds, indicate if it is integer, etc.)
2. Construct a concrete subtype of [`InfOptVariableType`](@ref) to specify the 
   desired type and its required additional information if appropriate
3. Build the variable object via `JuMP.build_variable`
4. Add the variable object to an `InfiniteModel` and assign a name via 
   `JuMP.add_variable`
5. Create a [`GeneralVariableRef`](@ref) that points to the variable object 
   stored in the model

!!! note 
    This methodology is presented for those wanting to learn more about the ins 
    and outs of variable definition. We recommend that all variables be created 
    via `@variable`. See [Macro Variable Definition](@ref).

The `JuMP.VariableInfo` data structure stores the following variable information:
- `has_lb::Bool`: Specifies a `Bool` it has a lower bound
- `lower_bound::Real`: Specifies lower bound value
- `has_ub::Bool`: Specifies a `Bool` it has a upper bound
- `upper_bound::Real`: Specifies upper bound value
- `has_fix::Bool`: Specifies a `Bool` it is fixed
- `fixed_value::Real`: Specifies the fixed value
- `has_start::Bool`: Specifies a `Bool` it has a start value
- `start::Union{Real, Function}`: Specifies the start guess value, this can be a
                                  function for infinite variables that intakes a
                                  support and maps it to a guess value (allowing
                                  to specify guess trajectories)
- `binary`: Specifies `Bool` if it is binary
- `integer`: Specifies `Bool` if it is integer.
Thus, the user specifies this information to prepare such an object:
```jldoctest genvar_define; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> info = VariableInfo(true, 0., true, 42., false, 0., false, 0., false, true);
```
Here we specified a lower bound of 0, an upper bound of 42, and that it is
integer valued.

The variable type objects (`InfOptVariableType` subtypes) are used with 
`build_variable` to specify the desired variable type along with any additional 
information needed for that type. For example, let's build an infinite variable 
`y(t)` that has an lower bound of 0, an upper bound of 42, and is integer valued:
```jldoctest genvar_define
julia> @infinite_parameter(model, t in [0, 10])
t

julia> info = VariableInfo(true, 0, true, 42, false, 0, false, 0, false, true);

julia> inf_var = build_variable(error, info, Infinite(t));
```
Thus, we create an [`InfiniteVariable`](@ref) object with the desired properties.

Once a variable has been built, it needs to be added to our `model` and a Julia
variable should be defined to reference it. Variables are added via
[`add_variable`](@ref JuMP.add_variable(::InfiniteModel, ::JuMP.AbstractVariable, ::String))
which adds a variable object to the model, assigns a name to the variable,
adds any constraints associated with the `JuMP.VariableInfo`, and returns
an appropriate variable reference variable (a [`GeneralVariableRef`](@ref)).
For example, let's add `inf_var` to `model`:
```jldoctest genvar_define
julia> var_ref = add_variable(model, inf_var, "y")
y(t)
```
Thus, we have added an infinite variable `y` that is parameterized by `t` with the
variable information mentioned above and now have a `GeneralVariableRef`
called `var_ref` that can be used in defining our infinite model.

Note that the use of `GeneralVariableRef`s and the corresponding concrete subtypes
of [`DispatchVariableRef`](@ref)s is discussed on the [Expressions](@ref expr_docs)
page.

## Macro Variable Definition
The [`@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
macro automates the variable definition process discussed above in the 
[Variable Definition Methodology](@ref) section via a straightforward symbolic 
syntax. The only key difference is that non-anonymous macro calls will register 
variable names to ensure they are not repeated. Anonymous macro calls forgo this 
step and exactly follow the process described above. This section will highlight 
the details of using this macro which is the recommended way to define variables.

!!! tip
    `JuMP`'s [documentation on variables](https://jump.dev/JuMP.jl/v1/manual/variables/) 
    is a good place to start since `InfiniteOpt` simply extends `JuMP` to 
    accommodate our additional variable types.

We directly build upon 
[`JuMP.@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
to create all of our decision variable types. To illustrate this via example, 
let's setup a model with a variety of infinite parameters ``t \in [0,10]``, 
``x \in [-1, 1]^3``, and ``\xi \in \mathcal{N}(0, 1)``:
```jldoctest var_macro
julia> using InfiniteOpt, Distributions

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_parameter(model, x[1:3] in [-1, 1], independent = true);

julia> @infinite_parameter(model, ξ ~ Normal());
```

### Variable Types
We specify the variable type by providing a subtype of [`InfOptVariableType`](@ref) 
as an extra positional argument:
```jldoctest var_macro
julia> @variable(model, y, Infinite(t, x..., ξ)) # explicit infinite variable
y(t, x[1], x[2], x[3], ξ)

julia> @variable(model, ys, SemiInfinite(y, 0, x..., ξ)) # explicit semi-infinite variable
ys

julia> @variable(model, yp, Point(y, 0, [1, 1, 1]..., 0)) # explicit point variable
yp

julia> @variable(model, z) # explicit finite variable
z
```
For anonymous definition, we use the `variable_type` keyword argument instead:
```jldoctest var_macro
julia>  y = @variable(model, variable_type = Infinite(t, x..., ξ)) # anon infinite variable
noname(t, x, ξ)

julia> ys = @variable(model, variable_type = SemiInfinite(y, 0, x..., ξ)) # anon semi-infinite variable
noname(0, x[1], x[2], x[3], ξ)

julia> yp = @variable(model, variable_type = Point(y, 0, [1, 1, 1]..., 0)) # anon point variable
noname(0, 1, 1, 1, 0)

julia> z = @variable(model) # anon finite variable
noname
```
Please refer to [`Infinite`](@ref), [`SemiInfinite`](@ref), and [`Point`](@ref) 
for more information.

### Variable Names
Variable inherit their names from the symbolic literal given with explicit 
definitions:
```jldoctest var_macro
julia> @variable(model, myname, Infinite(t))
myname(t)
```
This creates an infinite variable with name `"myname"` that is added to `model` 
and creates a Julia variable `myname` that stores a [`GeneralVariableRef`](@ref) 
which points to the infinite variable in `model`.

We can overwrite the inherited name using the `base_name` keyword argument:
```jldoctest var_macro
julia> @variable(model, myjlvar, Infinite(t), base_name = "myname")
myname(t)
```
This creates an infinite variable with name `"myname"` that is added to `model` 
and creates a Julia variable `myjlvar` that stores a [`GeneralVariableRef`](@ref) 
which points to the infinite variable in `model`.

This syntax is particularly useful for anonymous variables to have meaningful 
names:
```jldoctest var_macro
julia> myjlvar = @variable(model, variable_type = Infinite(t), base_name = "myname")
myname(t)
```

See the Queries and Modification sections further below for more information on 
how to query/modify variable names.

### Variable Bounds
We can specify variable bounds in like manner to `JuMP` variables. Let's 
demonstrate this with infinite variables: 
```jldoctest var_macro
julia> @variable(model, y_lb >= 0, Infinite(t, x...)) # add w/ lower bound
y_lb(t, x[1], x[2], x[3])

julia> @variable(model, y_ub <= 10, Infinite(t, x...)) # add w/ upper bound
y_ub(t, x[1], x[2], x[3])

julia> @variable(model, 0 <= y_bd <= 10, Infinite(t, x...)) # add w/ bounds
y_bd(t, x[1], x[2], x[3])

julia> @variable(model, y_fix == 42, Infinite(t, x...)) # add w/ fixed value 
y_fix(t, x[1], x[2], x[3])
```

!!! warning
    When creating a variable with only a single bound and the value of the bound 
    is not an explicit numeric literal, the name of the variable must appear on 
    the left-hand side. Otherwise, the macro will error.
    ```julia
    @variable(model, 0 <= y, Infinite(t)) # okay

    a = 0
    @variable(model, a <= y, Infinite(t)) # bad 
    @variable(model, y >= a, Infinite(t)) # okay
    ```

For anonymous definition, we use the `lower_bound` and `upper_bound`. Let's use 
finite variables for example:
```jldoctest var_macro
julia> z_lb = @variable(model, lower_bound = 0, base_name = "z_lb") # add w/ lower bound
z_lb

julia> z_ub = @variable(model, upper_bound = 10, base_name = "z_ub") # add w/ upper bound
z_ub

julia> z_bd = @variable(model, lower_bound = 0, upper_bound = 10, 
                        base_name = "z_bd") # add w/ bounds
z_bd

julia> z_fix = @variable(model, lower_bound = 10, upper_bound = 10, 
                         base_name = "z_fix") # ~add w/ fixed value 
z_fix
```
Note that there isn't a keyword for fixing variables. Instead, 
[`fix`](@ref JuMP.fix(::UserDecisionVariableRef, ::Real)) should be used. 

See the Queries and Modification sections further below for more information on 
how to query/modify variable bounds.

!!! note 
    Point variables inherit all the bounds of their respective infinite variables 
    by default. This can be overwritten by specifying different ones at creation.
    ```julia
    @variable(model, y >= 0, Infinite(t, x...)) # has lower bound
    @variable(model, yp == 0, Point(w, 0, 0, 0, 0)) # forces the point to be fixed
    ```

!!! note 
    Bounds cannot be specified on creation for semi-infinite variables. Note that 
    they will inherit these from the infinite variable they depend on. Additional 
    bound be created by directly adding constraints. For example:
    ```julia
    @variable(model, y >= 0, Infinite(t, x...)) # has lower bound
    @variable(model, ys, SemiInfinite(w, 0, x...)) # inherits the lower bound
    @constraint(model, ys <= 10) # add upper bound to ys
    ```

### Variable Integrality 
We can constrain the integrality of decision variables in like manner to `JuMP` 
using the `Bin` and `Int` positional arguments for explicit macro definition:
```jldoctest var_macro
julia> @variable(model, y_bin, Infinite(t, x...), Bin) # add as binary variable
y_bin(t, x[1], x[2], x[3])

julia> @variable(model, y_int, Infinite(t, x...), Int) # add as integer variable
y_int(t, x[1], x[2], x[3])
```

For anonymous definition, we use the `binary` and `integer` keyword arguments:
```jldoctest var_macro
julia> y_bin = @variable(model, variable_type = Infinite(t, x...), binary = true)
noname(t, x[1], x[2], x[3])

julia> y_int = @variable(model, variable_type = Infinite(t, x...), integer = true)
noname(t, x[1], x[2], x[3])
```

Moreover, we can add bounds as needed to constrain the domain of integer variables:
```jldoctest var_macro
julia> @variable(model, 0 <= y_int2 <= 10, Infinite(t, x...), Int)
y_int2(t, x[1], x[2], x[3])
```

See the Queries and Modification sections further below for more information on 
how to query/modify variable integralities.

!!! note 
    Point variables inherit the integrality of their respective infinite variables 
    by default. This can be overwritten by specifying different ones at creation.
    ```julia
    @variable(model, y, Infinite(t, x...), Bin) # is binary
    @variable(model, yp, Point(w, 0, 0, 0, 0), Int) # is integer
    ```

!!! note 
    Integrality cannot be specified for semi-infinite variables. Note that 
    they will inherit these from the infinite variable they depend on.
    ```

### Variable Start Values
Optimization solvers often benefit from giving initial guesses for the optimal 
decision variable values. Following `JuMP` vernacular, these are called start 
values. We use the keyword `start` to specify these at variable creation:
```jldoctest var_macro
julia> @variable(model, z_start, start = 42)
z_start
```

Moreover, infinite variables can accept a function that specifies the start 
value of over the range of its infinite parameters (e.g., a function that provides 
an initial guess trajectory). For example, consider the difference between these 
two infinite variables:
```jldoctest var_macro
julia> @variable(model, y_uniform, Infinite(t), start = 0) # start with y(t) = 0
y_uniform(t)

julia> @variable(model, y_sin, Infinite(t), start = sin) # start with y(t) = sin(t)
y_sin(t)
```
Note that such start functions must be able to accept parameter values as 
arguments that exactly match the format of the infinite parameters given in 
`Infinite(params...)`.

!!! note 
    Start values be specified for semi-infinite variables. Note that 
    they will inherit these from the infinite variable they depend on.

See the Queries and Modification sections further below for more information on 
how to query/modify variable names.

### Variable Containers
Optimization problems often involve multi-dimensional decision variables. Luckily, 
`JuMP` provides a versatile syntax for specifying collections (i.e., containers) 
of variables. See 
[JuMP's container documentation](https://jump.dev/JuMP.jl/v1/manual/containers/) 
for a thorough tutorial on the syntax. It uses `Array`s, `DenseAxisArray`s, and 
`SparseAxisArray`s to contain the variable references created. Here 
`DenseAxisArray`s and `SparseAxisArray`s allow the use of nontraditional indices 
(i.e., can use indices that are not sequential integers).

To illustrate what this means, consider the two equivalent ways to define 
a 3-dimensional vector of variables with indices `[1, 2, 3]`:
```jldoctest var_macro
julia> s = [0, 2, 1];

julia> var_refs = @variable(model, [i = 1:3], start = s[i], base_name = "z")
3-element Vector{GeneralVariableRef}:
 z[1]
 z[2]
 z[3]

julia> var_refs = Vector{GeneralVariableRef}(undef, 3);

julia> for i in eachindex(var_refs)
          var_refs[i] = @variable(model, start = s[i], base_name = "z")
       end
```

Moreover, here are a few illustrative examples:
```jldoctest var_macro
julia> @variable(model, z_dense[2:4])
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 2:4
And data, a 3-element Vector{GeneralVariableRef}:
 z_dense[2]
 z_dense[3]
 z_dense[4]

julia> @variable(model, z_named[[:A, :C, :Z]])
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, [:A, :C, :Z]
And data, a 3-element Vector{GeneralVariableRef}:
 z_named[A]
 z_named[C]
 z_named[Z]

julia> @variable(model, z_sparse[i = 1:2, j = 1:2; i + j <= 3])
JuMP.Containers.SparseAxisArray{GeneralVariableRef, 2, Tuple{Int64, Int64}} with 3 entries:
  [1, 1]  =  z_sparse[1,1]
  [1, 2]  =  z_sparse[1,2]
  [2, 1]  =  z_sparse[2,1]
```

The variable macro will by default automatically detect which container type 
should be used. However, the user can specify a particular container type using 
the `container` keyword. For example, if we want to use indices `a:b` where 
`a = 1` and `b = 3`, a `DenseAxisArray` will be used by default, but we can 
force it to be a regular `Array`:
```jldoctest var_macro
julia> a = 1; b = 3;

julia> var_refs1 = @variable(model, [a:b], base_name = "z")
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:3
And data, a 3-element Vector{GeneralVariableRef}:
 z[1]
 z[2]
 z[3]

julia> var_refs2 = @variable(model, [a:b], base_name = "z", container = Array)
3-element Vector{GeneralVariableRef}:
 z[1]
 z[2]
 z[3]
```

### Variable Sets
Like `JuMP` variables, we can constrain variables on creation to lie in 
particular sets. This allows us to make semi-definite variables, cone constrained 
variables, and more. 

For example:
```jldoctest var_macro
julia> @variable(model, z_psd[1:2, 1:2], PSD) # positive semi-definite variable matrix
2×2 LinearAlgebra.Symmetric{GeneralVariableRef, Matrix{GeneralVariableRef}}:
 z_psd[1,1]  z_psd[1,2]
 z_psd[1,2]  z_psd[2,2]

julia> @variable(model, z_cone[1:3] in SecondOrderCone()) # 2nd order cone variables
3-element Vector{GeneralVariableRef}:
 z_cone[1]
 z_cone[2]
 z_cone[3]
```
Typically, variable sets can be defined symbolically using the syntax 
`var in set`. For anonymous variables, the `set` keyword argument must be used:
```jldoctest var_macro
julia> z_cone = @variable(model, [1:3], set = SecondOrderCone())
3-element Vector{GeneralVariableRef}:
 noname
 noname
 noname
```

For a more thorough tutorial please see 
[JuMP's semi-definite documentation](https://jump.dev/JuMP.jl/v1/manual/variables/#Semidefinite-variables) 
and/or [JuMP's variables constrained on creation documentation](https://jump.dev/JuMP.jl/v1/manual/variables/#Variables-constrained-on-creation).

### Anonymous Variables
Above we talked showed the syntax for both explicit and anonymous variable 
creation. Anonymous creation is typically helpful in the following situations:
- defining multiple variables with the same name
- creating variables in user defined extensions
- using nontraditional naming

For anonymous variables, the only accepted positional arguments are the `model` 
and the container expression `[indices...]`. Everything else must be specified 
via keyword arguments `kwargs...` as shown in the subsections above.
```julia
@variable(model, [indices...], kwargs...)
```

For more information, see 
[JuMP's anonymous variable documentation](https://jump.dev/JuMP.jl/v1/manual/variables/#Anonymous-JuMP-variables).

### The `@variables` Macro
When using many `@variable` calls, we can instead use 
[`@variables`](https://jump.dev/JuMP.jl/v1/manual/variables/#variables) to 
enhance the readability:
```jldoctest var_macro
julia> @variables(model, begin
           y1, Infinite(t, x....)
           y2[i=1:2] >= i, Infinite(t), (start = i, base_name = "Y_$i")
           z2, Bin
       end)
(y1(t, x[1], x[2], x[3]), GeneralVariableRef[Y_1[1](t), Y_2[2](t)], z2)
```

## Restricted Variables
To define point and semi-infinite variables, we can also use [`restrict`](@ref) 
for convenient inline definitions.

For example, let's consider restricting the infinite variable `y(t, x)`:
```jldoctest restrict_vars; setup = :(using InfiniteOpt)
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 1]);

julia> @infinite_parameter(model, x[1:2] in [-1, 1]);

julia> @variable(model, y, Infinite(t, x))
y(t, x)

julia> pt = restrict(y, 0, [-1, 1]) # make point variable y(0, [-1, 1])
y(0, [-1, 1])

julia> semi = restrict(y, 0, x) # make semi-infinite variable y(0, x)
y(0, [x[1], x[2]])
```

We can also, even more conveniently, treat the infinite variable as a function 
to accomplish this in a more intuitive syntax:
```jldoctest restrict_vars
julia> pt = y(0, [-1, 1]) # make point variable y(0, [-1, 1])
y(0, [-1, 1])

julia> semi = y(0, x) # make semi-infinite variable y(0, x)
y(0, [x[1], x[2]])
```

## Queries
`InfiniteOpt` contains a large suite of methods to query information about
variables. This suite comprises extensions to all current `JuMP` query
methods and many more that are specific to `InfiniteOpt`. A number of the more
commonly used ones are explained in this section, but all the available methods
are explained in the [technical manual](@ref var_manual).

### General Information
Here we describe some methods used to query general variable information such as
the name. Variable names can be extracted via [`name`](@ref JuMP.name(::DecisionVariableRef))
which returns the name of a variable. The index of a variable (where it is stored 
in the infinite model) is accessed via 
[`index`](@ref JuMP.index(::GeneralVariableRef)) and the infinite model it 
belongs to is given by 
[`owner_model`](@ref JuMP.owner_model(::GeneralVariableRef)). These methods
are demonstrated below:
```jldoctest var_macro; setup = :(set_name(y, "y"))
julia> name(y)
"y"

julia> index(y)
InfiniteVariableIndex(2)

julia> model_where_stored = owner_model(y);
```

Also, [`num_variables`](@ref JuMP.all_variables(::InfiniteModel))
is useful in returning the total number of decision variables currently stored
in an infinite model:
```jldoctest var_macro
julia> num_variables(model)
61

julia> num_variables(model, PointVariable)
2
```
Similarly, [`all_variables`](@ref JuMP.all_variables(::InfiniteModel))
returns a list of all the variables currently added to the model.

Finally, [`variable_by_name`](@ref JuMP.variable_by_name(::InfiniteModel,::String))
can be employed to return the appropriate [`GeneralVariableRef`](@ref) based off 
of the variable name if it is unique. Returns `nothing` if such a name cannot be 
found and errors if it is not unique. For example, we can request the reference 
associated with `"y_ub"`:
```jldoctest var_macro
julia> variable_by_name(model, "y_ub")
y_ub(t, x[1], x[2], x[3])
```

### Variable Constraint Info
As described above, variables in `InfiniteOpt` can have constraints associated
with them like `JuMP` variables. These constraints include:
- lower bounds
- upper bounds
- fixed values
- binary valued
- integer valued.
Thus, a number of methods exist to query information about these constraints.

First, the ```[has/is]_[variable constraint type]``` methods indicate whether 
a variable has that particular constraint type. For example, to query if a 
variable `y_lb` has a lower bound we can use
[`has_lower_bound`](@ref JuMP.has_lower_bound(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> has_lower_bound(y_bd)
true
```
Thus, `y_bd` does have a lower bound. The other methods are
[`has_upper_bound`](@ref JuMP.has_upper_bound(::UserDecisionVariableRef)),
[`is_fixed`](@ref JuMP.is_fixed(::UserDecisionVariableRef)),
[`is_binary`](@ref JuMP.is_binary(::UserDecisionVariableRef)), and
[`is_integer`](@ref JuMP.is_integer(::UserDecisionVariableRef)).

Next, the ```[ConstraintType]Ref``` methods return an appropriate explicit type
[`InfOptConstraintRef`](@ref) that points to the constraint (errors if no
such constraint exists). For example, the upper bound constraint of `y_bd` can be
obtained via [`UpperBoundRef`](@ref JuMP.UpperBoundRef(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> UpperBoundRef(y_bd)
y_bd(t, x[1], x[2], x[3]) ≤ 10, ∀ t ∈ [0, 10], x[1] ∈ [-1, 1], x[2] ∈ [-1, 1], x[3] ∈ [-1, 1]
```
The other methods are [`LowerBoundRef`](@ref JuMP.LowerBoundRef(::UserDecisionVariableRef)),
[`FixRef`](@ref JuMP.FixRef(::UserDecisionVariableRef)),
[`BinaryRef`](@ref JuMP.BinaryRef(::UserDecisionVariableRef)), and
[`IntegerRef`](@ref JuMP.IntegerRef(::UserDecisionVariableRef)).

Finally, variable constraints that entail values (i.e., lower bounds, upper
bounds, and fixed values) have their values queried via the appropriate method.
For example, the lower bound value of `y_bd` is obtained via
[`lower_bound`](@ref JuMP.lower_bound(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> lower_bound(y_bd)
0.0
```
Note these methods error when no such constraint is associated with the variable.
The other methods are [`upper_bound`](@ref JuMP.upper_bound(::UserDecisionVariableRef))
and [`fix_value`](@ref JuMP.fix_value(::UserDecisionVariableRef)).

The start value can also be queried via 
[`start_value`](@ref JuMP.start_value(::UserDecisionVariableRef)) where nothing 
is returned if not start value is specified:
```jldoctest var_macro
julia> start_value(var_refs[1])
0.0

julia> start_value(yp)
```

For infinite and semi-infinite variables, the [`start_value_function`](@ref) 
should be used instead:
```jldoctest var_macro
julia> start_value_function(y_sin)
sin (generic function with 18 methods)
```

### Variable Use
`InfiniteOpt` defines a number of methods to track if and how variables are used
in an infinite model. For example,
[`used_by_constraint`](@ref used_by_constraint(::DecisionVariableRef)) is used to
determine if a variable is used by a constraint. For example, let's see if `y_bd`
is used by a constraint:
```jldoctest var_macro
julia> used_by_constraint(y_bd)
true
```
Other methods include [`used_by_measure`](@ref used_by_measure(::DecisionVariableRef))
and [`used_by_objective`](@ref used_by_objective(::DecisionVariableRef)).
For infinite variables, [`used_by_point_variable`](@ref) can also be used similarly.

Finally, in general [`is_used`](@ref is_used(::DecisionVariableRef)) can be used
to determine if a variable is used at all in the infinite model or not. For 
example, if we check `yp` using `is_used` we find that it isn't:
```jldoctest var_macro
julia> is_used(yp)
false
```

### Type Specific
`InfiniteOpt` also employs a few methods for specific variable types that return
information pertaining to that particular variable type. For infinite variables 
and semi-infinite variables,
[`parameter_refs`](@ref parameter_refs(::InfiniteVariableRef)) returns the tuple
of infinite parameters that the variable depends on. For example, consider
`y(t, x)`:
```jldoctest var_macro
julia> parameter_refs(y)
(t, x[1], x[2], x[3], ξ)
```

For point variables,
[`infinite_variable_ref`](@ref infinite_variable_ref(::PointVariableRef)) and
[`parameter_values`](@ref parameter_values(::PointVariableRef))
return the infinite variable it depends on and the infinite parameter point
values, respectively. For example, consider the point variable `yp`:
```jldoctest var_macro
julia> infinite_variable_ref(yp)
y(t, x[1], x[2], x[3], ξ)

julia> parameter_values(yp)
(0.0, 1.0, 1.0, 1.0, 0.0)
```

## Modification
`InfiniteOpt` employs a wide variety of methods to modify/delete variables.
These are comprised of `JuMP` extensions and methods native only to `InfiniteOpt`.
This section will highlight some of the more commonly used ones. All the
methods/macros are detailed in the [technical manual](@ref var_manual).

### Deletion
Like `JuMP v0.19+`, `InfiniteOpt` fully supports deletion throughout its data
types. Any variable and its dependencies can be deleted via
[`delete`](@ref JuMP.delete(::InfiniteModel, ::DecisionVariableRef)). Thus, when 
`delete` is invoked any bound/type constraints associated with the variable will 
be removed and it will be removed from any other constraints, measures, and/or 
objectives. For example, if we delete `y(t, x, ξ)` it will be removed along with 
its  bounds and the point variable `yp` will also be removed since it is a 
dependent:
```jldoctest var_macro
julia> delete(model, y)

julia> is_valid(model, yp)
false
```

Another class of deletion methods correspond to variable constraints. For example,
[`delete_lower_bound`](@ref JuMP.delete_lower_bound(::UserDecisionVariableRef)) is
used to delete a lower bound associated with a
variable if it has one. Let's illustrate this by deleting the lower bound of
`y_bd`:
```jldoctest var_macro
julia> delete_lower_bound(y_bd)

julia> has_lower_bound(y_bd)
false
```
Other similar methods are
[`delete_upper_bound`](@ref JuMP.delete_upper_bound(::UserDecisionVariableRef)),
[`unfix`](@ref JuMP.unfix(::UserDecisionVariableRef)),
[`unset_binary`](@ref JuMP.unset_binary(::UserDecisionVariableRef)), and
[`unset_integer`](@ref JuMP.unset_integer(::UserDecisionVariableRef)).

### Variable Constraints
Another class of methods seek to add/modify variable constraints such as
bounds. For example, [`set_lower_bound`](@ref JuMP.set_lower_bound(::UserDecisionVariableRef, ::Real))
specifies the lower bound of a
variable. We can add a lower bound of 0 to `z` by:
```jldoctest var_macro
julia> set_lower_bound(z, 0)

julia> lower_bound(z)
0.0
```
Thus, adding a lower bound to `z`. Furthermore, we can later modify the lower
bound using the same method:
```jldoctest var_macro
julia> set_lower_bound(z, -2)

julia> lower_bound(z)
-2.0
```
Other similar methods are
[`set_upper_bound`](@ref JuMP.set_upper_bound(::UserDecisionVariableRef, ::Real)),
[`fix`](@ref JuMP.fix(::UserDecisionVariableRef, ::Real)),
[`set_binary`](@ref JuMP.set_binary(::UserDecisionVariableRef)), and
[`set_integer`](@ref JuMP.set_integer(::UserDecisionVariableRef)).

### Start Values
We can update the start value of a variable using 
[`set_start_value`](@ref JuMP.set_start_value(::UserDecisionVariableRef, ::Real)).
For example:
```jldoctest var_macro
julia> set_start_value(z, 0)

julia> start_value(z)
0.0
```

For infinite variables, this should be done using 
[`set_start_value_function`](@ref). FOr example:
```jldoctest var_macro
julia> set_start_value_function(myname, sin)

julia> start_value_function(myname)
sin (generic function with 18 methods)
```
Again note that such start functions must be able to accept parameter values as 
arguments that exactly match the format of the infinite parameters given in 
`Infinite(params...)`.

A number of other techniques exist for the various variable types can be found in 
the manual below.
