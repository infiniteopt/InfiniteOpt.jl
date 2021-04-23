# Variables
A guide and manual for the definition and use of variables in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview
Decision variables are at the very core of `InfiniteOpt` as its name alludes
to mathematical programs that entail infinite decision spaces (i.e., contain
infinite decision variables). Principally, three variable types are employed:
infinite, point, and finite. Infinite variables encompass any decision variable
that is parameterized by an infinite parameter (e.g., space-time variables and
recourse variables). Point variables are infinite variables at a particular
infinite parameter value (point). Finally, finite variables are decisions that
are made irrespective of the infinite domain (e.g., first stage variables and
design variables). Or in other words, they hold a particular value over the
infinite domain or some sub-domain of it. We also employ semi-infinite variables
which correspond to partially transcribed infinite variables which are artifacts
of measure evaluation and thus are discussed on the [Measures](@ref measure_page) page.

## Basic Usage
Infinite, point, and finite variables are summarized in the following table:

| Variable Type | Description                                    | Examples                               |
|:-------------:|:----------------------------------------------:|:--------------------------------------:|
| Infinite      | Functional variables with infinite co-domain   | ``y(t)``, ``y(\xi)``, ``y(t, x)``      |
| Point         | Infinite variable evaluated at parameter point | ``y(0)``, ``y(t_0, x_0)``              |
| Finite          | Held constant over infinite domain             | ``z`` (design and 1st stage variables) |

Infinite, point, and finite variables are typically defined via their respective
macros: [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@finite_variable`](@ref). These macros generally emulate [`JuMP.@variable`](@ref)
except that they each employ additional syntax capabilities to employ their
respective variable type.

!!! warning
    The "hold variable" nomenclature used by previous versions of `InfiniteOpt` has 
    been removed in favor of "finite variable". Hence, `@hold_variable` is now 
    deprecated in favor of `@finite_variable`.

Let's first setup a simple space-time model with infinite parameters time `t` and
spatial position `x`:
```jldoctest var_basic
julia> using InfiniteOpt

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_parameter(model, x[1:2] in [-1, 1], independent = true)
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]
```

### Infinite Variables
Now let's define a time dependent infinite variable `y(t)` with a lower bound of
0:
```jldoctest var_basic
julia> @infinite_variable(model, y(t) >= 0)
y(t)
```
This creates a Julia variable `y` that points to the decision variable `y(t)`
that is stored in `model` which is added to include a lower bound of 0. Another
useful case is that of defining an array of variables `w` that depend on both
position and time:
```jldoctest var_basic
julia> @infinite_variable(model, w[i = 1:3](t, x), start = [0, 2, 1][i])
3-element Array{GeneralVariableRef,1}:
 w[1](t, x)
 w[2](t, x)
 w[3](t, x)
```
Thus we create a Julia array variable `w` whose elements `w[i]` point to their
respective infinite variables `w[i](t, x)` stored in `model`. Note that the `i`
used in the array definition can be used to index attributes assigned to each
variable in the array. In this case, we used `i` to assign different initial
guess values for each variable via the `start` keyword argument.

Moreover, for infinite variables a function can be given to determine the start
values over a range of support points (e.g., a guess trajectory). This is
discussed further below in the Macro Definition section.

### Point Variables
Now let's add some point variables. These allow us to consider an infinite
variable evaluated at a certain infinite parameter point. For example, let's
define a point variable for `y(0)` with the alias `y0` that is fixed at a value
of 0:
```jldoctest var_basic
julia> @point_variable(model, y(0), y0 == 0)
y0
```
Here we create a Julia variable `y0` which points to the point variable `y(0)`.
Notice that in the second argument we specify the infinite variable indexed at
the appropriate parameter value(s). Point variables automatically inherit
attributes of the infinite variable (e.g., bounds, start values, etc.), but
these are overwritten with properties specified for the point variable. In this
case the lower bound inherited from `y(t)` is overwritten by instead fixing
`y(0)` to a value of 0.  

!!! note
    Point variables are provided for enhancing the generality of
    `InfiniteOpt`, but typically can be avoided by using infinite variables in
    combination with [`@BDconstraint`](@ref) to define bounded constraints (
    e.g., initial conditions).

### Finite Variables
Finally, we can add finite variables to our model. These denote variables that
hold a single value over the infinite domain or some portion of it (e.g.,
design variables, first stage variables, etc.). Let's add a finite variable
``0 \leq d \leq 42`` that is an integer variable and defined over all infinite
domains (i.e., time and space):
```jldoctest var_basic
julia> @finite_variable(model, 0 <= d <= 42, Int)
d
```
This creates a Julia variable `d` that points to the finite variable `d` which has
a lower bound of 0, an upper bound of 42, and is an integer variable. Thus,
[`@finite_variable`](@ref) follows the same exact syntax as [`JuMP.@variable`](@ref)
except that it also allows the user to specify a subdomain over which the finite
variable is valid. For example, let's add a finite variable `z` that is only valid
over the subdomain ``t \in [0, 5]`` via the `parameter_bounds` keyword argument:
```jldoctest var_basic
julia> anon = @finite_variable(model, parameter_bounds = (t in [0, 5]),
                             base_name = "z")
z
```
Here we make an anonymous variable for the sake of example whose reference is
stored to the Julia variable `anon` and points to a finite variable `z` which is
only valid for ``t \in [0, 5]``. Thus, this will be enforced in any constraints that
involve `anon`, meaning they will automatically be bounded to such a subdomain.
Any number of parameters bounds (bounds on the parameters of the infinite
domain) can be added in a tuple like argument as explained in the documentation
for [`@finite_variable`](@ref).

Now we have defined variables that we can use in the objective, measures, and
constraints. Please note that the above tutorial only shows a small portion of
the capabilities and options available in defining variables. A full description
is provided in the documentation below.

## Variable Definition Methodology
The [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@finite_variable`](@ref) macros all follow a similar methodology behind the
scenes and these commonalities are discussed in this section for conciseness.
Defining/initializing a variable principally involves the following steps:
1. Define the variable information pertaining to `JuMP.VariableInfo` (e.g., bounds, indicate if it is integer, etc.)
2. Construct a concrete subtype of [`InfOptVariable`](@ref) to store the variable information
3. Add the `InfOptVariable` object to an `InfiniteModel` and assign a name
4. Create a [`GeneralVariableRef`](@ref) that points to the variable object stored in the model

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
julia> info = VariableInfo(true, 0., true, 42., false, 0., false, 0., false, true)
VariableInfo{Float64,Float64,Float64,Float64}(true, 0.0, true, 42.0, false, 0.0, false, 0.0, false, true)
```
Here we specified a lower bound of 0, an upper bound of 42, and that it is
integer.

The variable objects (`InfOptVariable` subtypes) are defined via
[`build_variable`](@ref JuMP.build_variable(::Function,::JuMP.VariableInfo,::Symbol))
which requires that the user provide a
`JuMP.VariableInfo` object, the variable type to be defined (`Infinite`,
`Point`, or `Finite`), and any necessary keyword arguments required for that
variable type (i.e., `parameter_refs`, `infinite_variable_ref`, and/or
`parameter_values`). For example, let's build an infinite variable `y(t)` that
has an lower bound of 0, an upper bound of 42, and is integer:
```jldoctest genvar_define
julia> @infinite_parameter(model, t in [0, 10])
t

julia> info = VariableInfo(true, 0, true, 42, false, 0, false, 0, false, true);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = (t));
```
Thus, we create an [`InfiniteVariable`](@ref) object with the desired properties.
Note that in this case the `parameter_refs` keyword argument is required to
indicate which infinite parameter(s) this infinite variable will depend on.

Once a variable has been built, it needs to be added to our `model` and a Julia
variable should be defined to reference it. Variables are added via
[`add_variable`](@ref JuMP.add_variable(::InfiniteModel, ::InfOptVariable, ::String))
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
of [`DispatchVariableRef`](@ref)s is discussed on the [Expressions](@ref expr_page)
page.

## Macro Variable Definition
The [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@finite_variable`](@ref) macros automate the variable definition process
discussed above in the [Variable Definition Methodology](@ref) section
via a straightforward symbolic syntax. The only key difference
is that non-anonymous macro calls will register variable names to ensure they
are not repeated. Anonymous macro calls forgo this step and exactly follow the
process described above. This section will highlight the details of using these
macros.

### [General Usage] (@id var_macro_gen_usage)
Here we discuss the features that the variable macros have in common (generally
these pertain to `JuMP`-like features). To illustrate this via example, let's
setup a model with a variety of infinite parameters ``t \in [0,10]``,
``x \in [-1, 1]^3``, and ``\xi \in \mathcal{N}(0, 1)``:
```jldoctest var_macro
julia> using InfiniteOpt, Distributions

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_parameter(model, x[1:3] in [-1, 1], independent = true);

julia> @infinite_parameter(model, ξ in Normal());
```

We will first consider anonymous variable macro calls which generally are less
convenient than non-anonymous macro calls which offer a much more intuitive
mathematical syntax. However, anonymous variables can be useful and provide a
good foundation to understanding non-anonymous variables.
Furthermore, we'll use finite variables as the motivating examples since they
best exemplify commonalities between the macros. First, let's consider single
anonymous definition a finite variable:
```jldoctest var_macro
julia> var_ref = @finite_variable(model)
noname
```
Here we just added a nameless finite variable to `model` and defined `var_ref`
as a `GeneralVariableRef` that points to it. We can add a name via the
`base_name` keyword argument:
```jldoctest var_macro
julia> var_ref1 = @finite_variable(model, base_name = "d")
d

julia> var_ref2 = @finite_variable(model, base_name = "d")
d
```
Now we've made 2 more finite variables both called `d`. Thus, the anonymous syntax
allows us to define variables with the same name. Moreover, any variable
information can be specified via the appropriate keywords which include:
- `lower_bound::Real`: specifies lower bound
- `upper_bound::Real`: specifies upper bound
- `start::Union{Real, Function}`: specifies the initial guess value the solver will use,
                                  for infinite variables this can be a function that
                                  takes as input a realization of the infinite parameters
                                  (using the same format as `parameter_refs`) and return the
                                  guess value for such as realization (e.g., functionally
                                  define a trajectory).
- `binary::Bool`: specifies if is binary variable
- `integer::Bool`: specifies if is integer variable.
Anonymous variables must use these keyword arguments since symbolic definition
is only permitted for non-anonymous macro calls. For example, let's define a
finite variable ``0 \leq d \leq 5`` that is integer:
```jldoctest var_macro
julia> var_ref = @finite_variable(model, base_name = "d", lower_bound = 0,
                                upper_bound = 5, integer = true)
d
```

We can also define arrays of variables using any indices of our choice. For
example, let's define a 3-dimensional vector with indices `[1, 2, 3]`:
```jldoctest var_macro
julia> var_refs = @finite_variable(model, [i = 1:3], start = [0, 2, 1][i],
                                 base_name = "d")
3-element Array{GeneralVariableRef,1}:
 d[1]
 d[2]
 d[3]
```
Thus, we define 3 variables named `d[i]` and each with different start values
and define `var_refs` which is a vector of `GeneralVariableRef`s that uses the
indices we specified. Note the syntax `i = indices` is used to define an
iteration variable to use with the keyword arguments to assign different values
for each variable being defined. Note the above example is equivalent to:
```jldoctest var_macro
julia> starts = [0, 2, 1];

julia> var_refs = Vector{GeneralVariableRef}(undef, 3);

julia> for i = eachindex(var_refs)
          var_refs[i] = @finite_variable(model, base_name = "d", start = starts[i])
       end
```
Other non-standard indices can also be used such as the following examples:
```jldoctest var_macro
julia> var_refs2 = @finite_variable(model, [2:4], base_name = "d")
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 2:4
And data, a 3-element Array{GeneralVariableRef,1}:
 d[2]
 d[3]
 d[4]

julia> var_refs3 = @finite_variable(model, [[:A, :C, :Z]], base_name = "d")
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, Symbol[:A, :C, :Z]
And data, a 3-element Array{GeneralVariableRef,1}:
 d[A]
 d[C]
 d[Z]

julia> var_refs3 = @finite_variable(model, [i=1:2, j=i:2], base_name = "d")
JuMP.Containers.SparseAxisArray{GeneralVariableRef,2,Tuple{Int64,Int64}} with 3 entries:
  [1, 2]  =  d[1,2]
  [2, 2]  =  d[2,2]
  [1, 1]  =  d[1,1]
```
Here we see that a variety of indices can be used and this is explained more
fully in the documentation of `JuMP`.

`JuMP` employs 2 special array container types: `DenseAxisArray`s and
`SparseAxisArray`s which help facilitate this special indexing. The variable
macros will by default automatically detect which container type should be used.
However, the user can specify a particular container type using the
`container` keyword. For example, if we want to use indices `a:b` where `a = 1`
and `b = 3`, a `DenseAxisArray` will be used by default, but we can force it to
be a regular `Array`:
```jldoctest var_macro
julia> a = 1; b = 3;

julia> var_refs1 = @finite_variable(model, [a:b], base_name = "d")
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:3
And data, a 3-element Array{GeneralVariableRef,1}:
 d[1]
 d[2]
 d[3]

julia> var_refs2 = @finite_variable(model, [a:b], base_name = "d", container = Array)
3-element Array{GeneralVariableRef,1}:
 d[1]
 d[2]
 d[3]
```
For more information on `JuMP` containers please visit their page
[here](http://www.juliaopt.org/JuMP.jl/stable/containers/).

Now that we have a foundation with anonymous variable macro calls, let's focus
on non-anonymous calls which offer a much more straightforward syntax. These calls
can still implement all of the same keyword arguments. Moreover, they automatically
create a Julia variable with the variable name provided and register this name
to ensure subsequent automatic Julia variables do not overwrite it.

The supported symbolic syntax principally implements the following keyword
arguments:
- `base_name`
- `lower_bound`
- `upper_bound`
- `integer`
- `binary`
- `parameter_refs` (for infinite variables)
- `infinite_variable_ref` (for point variables)
- `parameter_values` (for point variables).
These are implemented via the syntax
`@[type]_variable(model, expr, integrality_arg, keyword_args...)`. Here `expr`
specifies the name, bounds, and/or variable specific keyword arguments. It can
use the following forms (note that in the following the symbol `<=` can be used
instead of `≤` and the symbol `>=`can be used instead of `≥`):
- `varexpr` creating variables described by `varexpr`
- `varexpr ≤ ub` (resp. `varexpr ≥ lb`) creating variables described by
  `varexpr` with upper bounds given by `ub` (resp. lower bounds given by `lb`)
- `varexpr == value` creating variables described by `varexpr` with fixed values
   given by `value`
- `lb ≤ varexpr ≤ ub` or `ub ≥ varexpr ≥ lb` creating variables described by
  `varexpr` with lower bounds given by `lb` and upper bounds given by `ub`
Thus, providing an intuitive means to specify bounds. The expressions `varexpr`
specifies the name, dimensions, and/or type specific keywords and can be of the
form:
- `varname` creating a scalar real variable of name `varname`
- `varname[...]` creating a container of variables with indices `...`
- `varname(params)` creating an infinite variable dependent on `params`
- `varname[...](params)` creating infinite variables dependent on `params`.
The `integrality_arg` optionally is used to indicate if the variable(s) is/are
integer or binary using `Int` or `Bin`, respectively.  

For example, let's define a finite variable ``0 \leq d \leq 3`` that is integer:
```jldoctest var_macro
julia> @finite_variable(model, 0 <= d <= 3, Int)
d
```
Note this is equivalent to
```jldoctest var_macro
julia> d = @finite_variable(model, base_name = "d", lower_bound = 0, upper_bound = 3,
                          integer = true)
d
```
with the exception that the non-anonymous definition registers `d` as variable
name that cannot be duplicated. For one more example let's define a vector of
variables ``a \in \mathbb{R}_+^3`` with starting values of 0:
```jldoctest var_macro
julia> @finite_variable(model, a[1:3] >= 0, start = 0)
3-element Array{GeneralVariableRef,1}:
 a[1]
 a[2]
 a[3]
```
Thus, highlighting how keyword arguments can still be used.

### Infinite Variables
Infinite variables entail decision variables that depend on infinite parameter(s).
Thus, [`@infinite_variable`](@ref) follows the general definition methodology
with this additional consideration.

Let's first consider a basic anonymous definition of an infinite variable
``y(t, x)``:
```jldoctest var_macro
julia> y = @infinite_variable(model, parameter_refs = (t, x), base_name = "y")
y(t, x)
```
Here we created an infinite variable with the base name ``y`` that depends on
``t`` and ``x``. Notice that each group of parameters are specified in a
particular element of the `parameter_refs` tuple, this is the required format.
Moreover, the keyword argument `parameter_refs` is required to
specify what parameterizes the infinite variable. Because we used an anonymous
call, we can still make another variable with the same name. For example let's
define another infinite variable also called ``y`` that only depends on ``t``:
```jldoctest var_macro
julia> y2 = @infinite_variable(model, parameter_refs = (t), base_name = "y")
y(t)
```

More conveniently, we can equivalently define ``y(t, x)`` symbolically:
```jldoctest var_macro
julia> @infinite_variable(model, y(t, x))
y(t, x)
```
We can also use this symbolic syntax to add constraint information as described
in the previous section. For example, let's define a vector of infinite variables
``z(t) \in \{0, 1, 2\}^3``:
```jldoctest var_macro
julia> @infinite_variable(model, 0 <= z[1:3](t) <= 2, Int)
3-element Array{GeneralVariableRef,1}:
 z[1](t)
 z[2](t)
 z[3](t)
```
In many cases, it is important to specify a guess trajectory for an infinite
variable. Currently, this can be done by defining a start value function and
specifying via the `start` keyword. In particular, the arguments of this function
must match the format of the infinite parameters and return a start value for a
given realization (point) of the infinite parameters:
```jldoctest var_macro
julia> my_traj_func(t_val, x_vals) = t_val + sum(x_vals) # replace with actual functionality
my_traj_func (generic function with 1 method)

julia> @infinite_variable(model, y2(t, x), start = (t_val, x_vals) -> my_traj_func(t_val, x_vals))
y2(t, x)
```
Notice that above we use an anonymous function in the macro call, this is done
to avoid ambiguity in case we have multiple `my_traj_func` definitions. In
general it is recommended that anonymous functions be used.

### Point Variables
Point variables denote infinite variables (or derivatives) evaluated at a particular point in the
infinite decision space (i.e., particular infinite parameter values). These
are commonly employed when using initial and terminal conditions, and to build
discrete characterizations of complex operators such as derivatives. The
[`@point_variable`](@ref) macro is employed to define such variables. Principally,
it follows the general variable definition paradigm, but allows us to specify
the infinite variable it refers to and the parameter values it is evaluated at.
Also, note that by default it inherits the characteristics of the infinite variable
(e.g., bounds), but these are overwritten as specified in the point variable
macro.

To begin let's consider defining a boundary point ``y(0, -1)`` based on the
infinite variable ``y(t, x)`` and enforce that it be nonnegative. Note that
-1 need be a 3 element vector in this case to match the dimensions of ``x``.
The anonymous syntax would be:
```jldoctest var_macro
julia> y0 = @point_variable(model, infinite_variable_ref = y,
                            parameter_values = (0, [-1, -1, -1]), lower_bound = 0)
y(0, [-1, -1, -1])
```
This creates a point variable `y(0, [-1, -1, -1]) ≥ 0` that is added to `model`
and assigns to the associated `GeneralVariableRef` to the Julia variable
`y0`. Equivalently, this can accomplished much more conveniently via:
```jldoctest var_macro
julia> @point_variable(model, y(0, [-1, -1, -1]), y0 >= 0)
y0
```
Here the 2nd argument specifies the infinite variable and parameter values, and
the next argument is used to provide a convenient alias that can be used in
combination with the typical symbolic variable syntax described in the previous
sections. Let's also demonstrate how this works for multi-dimensional infinite
variables. For example consider defining ``z(0) = 0`` for
``z(t) \in \{0, 1, 2\}^3``:
```jldoctest var_macro
julia> @point_variable(model, z[i](0), z0[i = 1:3] == 0)
3-element Array{GeneralVariableRef,1}:
 z0[1]
 z0[2]
 z0[3]
```

### Finite Variables
Finite variables denote decision variables that are constant (agnostic) over the
infinite domain or some sub-domain of it. This is accomplished via the
[`@finite_variable`](@ref) macro as demonstrated in the
[General Usage](@ref var_macro_gen_usage) section. By default and as shown in
the above examples, finite variables are valid over the entire infinite domain.
However, this scope can be limited via the `parameter_bounds` keyword argument.

For example, let's define a finite variable ``b \in \{0, 1\}`` that is valid over
the entire infinite domain (any infinite parameter value):
```jldoctest var_macro
julia> @finite_variable(model, b, Bin)
b
```
Again, this follows the methodology outlined above. Now let's suppose we want to
define a finite variable ``0 \leq c \leq 42`` that is only valid over the time
interval ``t \in [0, 5]`` which is a subset of the entire range being considered.
This can be accomplished via `parameter_bounds`:
```jldoctest var_macro
julia> @finite_variable(model, 0 <= c <= 42, parameter_bounds = (t in [0, 5]))
c
```
Thus, we defined `c` and it can only be used in constraints and measures in
accordance with this limited sub-domain. When such a limited finite variable is
used in a constraint, the constraint parameter bounds be overlapped with those of
`c` if possible. Otherwise, an error will be thrown. This is further explained
on the [Constraints](@ref constr_page) page.

Any number of parameters can be specified in a finite variable's sub-domain. For
example, let's define `e` such over the domain ``t \in [0, 1]``, ``x = -1``:
```jldoctest var_macro
julia> @finite_variable(model, e, parameter_bounds = (t in [0, 1], x == -1))
e
```

## Queries
`InfiniteOpt` contains a large suite of methods to query information about
variables. This suite is comprised of extensions to all current `JuMP` query
methods and many more that are specific to `InfiniteOpt`. A number of the more
commonly used ones are explained in this section, but all of the available methods
are explained in the [Methods/Macros](@ref var_methods) section (i.e., the
manual) below.

### General Information
Here we describe some methods used to query general variable information such as
the name. Variable names can be extracted via [`name`](@ref JuMP.name(::DecisionVariableRef))
which returns the
name of a variable. The index of a variable (where it is stored in the infinite
model) is accessed via [`index`](@ref JuMP.index(::GeneralVariableRef)) and the
infinite model it belongs to is
given by [`owner_model`](@ref JuMP.owner_model(::GeneralVariableRef)). These methods
are demonstrated below:
```jldoctest var_macro
julia> name(y)
"y"

julia> index(y)
InfiniteVariableIndex(3)

julia> model_where_stored = owner_model(y);
```

Also, [`num_variables`](@ref JuMP.all_variables(::InfiniteModel,::Type{InfOptVariable}))
is useful in returning the total number of decision variables currently stored
in an infinite model:
```jldoctest var_macro
julia> num_variables(model)
45

julia> num_variables(model, PointVariable)
5
```
Similarly, [`all_variables`](@ref JuMP.all_variables(::InfiniteModel,::Type{InfOptVariable}))
returns a list of all the variables currently added to the model.

Finally, [`variable_by_name`](@ref JuMP.variable_by_name(::InfiniteModel,::String))
can be employed to return the appropriate
[`GeneralVariableRef`](@ref) based off of the variable name if
it is unique. Returns `nothing` if such a name cannot be found and errors if it
is not unique. For example, we can request the reference associated with `"c"`:
```jldoctest var_macro
julia> variable_by_name(model, "c")
c
```

### Variable Constraint Info
As described above, variables in `InfiniteOpt` can have constraints associated
with them like `JuMP` variables. These constraints include:
- lower bounds
- upper bounds
- fixed values
- binary specifications
- integer specifications.
Thus, a number of methods exist to query information about these constraints.

First, the ```[has/is]_[variable constraint type]``` methods indicate whether or not a
variable has that particular constraint type. For example, to query if a variable
`d` has a lower bound we can use
[`has_lower_bound`](@ref JuMP.has_lower_bound(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> has_lower_bound(d)
true
```
Thus, `d` does have a lower bound. The other methods are
[`has_upper_bound`](@ref JuMP.has_upper_bound(::UserDecisionVariableRef)),
[`is_fixed`](@ref JuMP.is_fixed(::UserDecisionVariableRef)),
[`is_binary`](@ref JuMP.is_binary(::UserDecisionVariableRef)), and
[`is_integer`](@ref JuMP.is_integer(::UserDecisionVariableRef)).

Next, the ```[ConstraintType]Ref``` methods return an appropriate explicit type
[`InfOptConstraintRef`](@ref) that points to the constraint (errors if no
such constraint exists). For example, the upper bound constraint of `d` can be
obtained via [`UpperBoundRef`](@ref JuMP.UpperBoundRef(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> UpperBoundRef(d)
d ≤ 3.0
```
The other methods are [`LowerBoundRef`](@ref JuMP.LowerBoundRef(::UserDecisionVariableRef)),
[`FixRef`](@ref JuMP.FixRef(::UserDecisionVariableRef)),
[`BinaryRef`](@ref JuMP.BinaryRef(::UserDecisionVariableRef)), and
[`IntegerRef`](@ref JuMP.IntegerRef(::UserDecisionVariableRef)).

Finally, variable constraints that entail values (i.e., lower bounds, upper
bounds, and fixed values) have their values queried via the appropriate method.
For example, the lower bound value of `d` is obtained via
[`lower_bound`](@ref JuMP.lower_bound(::UserDecisionVariableRef)):
```jldoctest var_macro
julia> lower_bound(d)
0.0
```
Note these methods error when no such constraint is associated with the variable.
The other methods are [`upper_bound`](@ref JuMP.upper_bound(::UserDecisionVariableRef))
and [`fix_value`](@ref JuMP.fix_value(::UserDecisionVariableRef)).

### Variable Use
`InfiniteOpt` defines a number of methods to track if and how variables are used
in an infinite model. For example,
[`used_by_constraint`](@ref used_by_constraint(::DecisionVariableRef)) is used to
determine if a variable is used by a constraint. For example, let's see if `c`
is used by a constraint:
```jldoctest var_macro
julia> used_by_constraint(c)
true
```
Other methods include [`used_by_measure`](@ref used_by_measure(::DecisionVariableRef))
and [`used_by_objective`](@ref used_by_objective(::DecisionVariableRef)).
For infinite variables, [`used_by_point_variable`](@ref) can also be used in a
similar manner.

Finally, in general [`is_used`](@ref is_used(::DecisionVariableRef)) can be used
to determine if a variable is
used at all in the infinite model or not. For example, if we check `e` using
`is_used` we find that it isn't:
```jldoctest var_macro
julia> is_used(e)
false
```

### Type Specific
`InfiniteOpt` also employs a few methods for specific variable types that return
information pertaining to that particular variable type. For infinite variables,
[`parameter_refs`](@ref parameter_refs(::InfiniteVariableRef)) returns the tuple
of infinite parameters that the variable depends on. For example, consider
`y(t, x)`:
```jldoctest var_macro
julia> parameter_refs(y)
(t, GeneralVariableRef[x[1], x[2], x[3]])
```

For point variables,
[`infinite_variable_ref`](@ref infinite_variable_ref(::PointVariableRef)) and
[`parameter_values`](@ref parameter_values(::PointVariableRef))
return the infinite variable it depends on and the infinite parameter point
values, respectively. For example, consider the point variable `y0`:
```jldoctest var_macro
julia> infinite_variable_ref(y0)
y(t, x)

julia> parameter_values(y0)
(0.0, [-1.0, -1.0, -1.0])
```

For finite variables,
[`has_parameter_bounds`](@ref has_parameter_bounds(::FiniteVariableRef)) returns
if a finite variable has parameter bounds (i.e., a specified sub-domain) and
[`parameter_bounds`](@ref parameter_bounds(::FiniteVariableRef))
returns those bounds if there are any. For example, consider `c`:
```jldoctest var_macro
julia> has_parameter_bounds(c)
true

julia> parameter_bounds(c)
Subdomain bounds (1): t ∈ [0, 5]
```

## Modification
`InfiniteOpt` employs a wide variety of methods to modify/delete variables.
These are comprised of `JuMP` extensions and methods native only to `InfiniteOpt`.
This section will highlight some of the more commonly used ones. All of the
methods/macros are detailed in the [Methods/Macros](@ref var_methods) section
(i.e., the manual) below.

### Deletion
Like `JuMP v0.19+`, `InfiniteOpt` fully supports deletion throughout its data
types. Any variable and its dependencies can be deleted via
[`delete`](@ref JuMP.delete(::InfiniteModel, ::DecisionVariableRef)). Thus, when `delete`
is invoked any bound/type constraints associated with the variable will be removed
and it will be removed from any other constraints, measures, and/or objectives.
For example, if we delete `y(t, x)` it will be removed along with its bounds and
the point variable `y0` will also be removed since it is a dependent:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel(); @finite_variable(model, y))
julia> delete(model, y)
```

Another class of deletion methods correspond to variable constraints. For example,
[`delete_lower_bound`](@ref JuMP.delete_lower_bound(::UserDecisionVariableRef)) is
used to delete a lower bound associated with a
variable if it has one. Let's illustrate this by deleting the lower bound of
`d`:
```jldoctest var_macro
julia> delete_lower_bound(d)

julia> has_lower_bound(d)
false
```
Other similar methods are
[`delete_upper_bound`](@ref JuMP.delete_upper_bound(::UserDecisionVariableRef)),
[`unfix`](@ref JuMP.unfix(::UserDecisionVariableRef)),
[`unset_binary`](@ref JuMP.unset_binary(::UserDecisionVariableRef)), and
[`unset_integer`](@ref JuMP.unset_integer(::UserDecisionVariableRef)).

Finally, [`delete_parameter_bounds`](@ref delete_parameter_bounds(::FiniteVariableRef))
can be used on finite variables to delete all of their parameter bounds. For
example, let's delete all of the parameter bounds associated with `c`:
```jldoctest var_macro
julia> parameter_bounds(c)
Subdomain bounds (1): t ∈ [0, 5]

julia> delete_parameter_bounds(c)

julia> has_parameter_bounds(c)
false
```

### Variable Constraints
Another class of methods seek to add/modify variable constraints such as
bounds. For example, [`set_lower_bound`](@ref JuMP.set_lower_bound(::UserDecisionVariableRef, ::Number))
specifies the lower bound of a
variable. We can add a lower bound of 0 to `c` by:
```jldoctest var_macro
julia> set_lower_bound(c, 0)

julia> lower_bound(c)
0.0
```
Thus, adding a lower bound to `c`. Furthermore, we can later modify the lower
bound using the same method:
```jldoctest var_macro
julia> set_lower_bound(c, -2)

julia> lower_bound(c)
-2.0
```
Other similar methods are
[`set_upper_bound`](@ref JuMP.set_upper_bound(::UserDecisionVariableRef, ::Real)),
[`fix`](@ref JuMP.fix(::UserDecisionVariableRef, ::Real)),
[`set_binary`](@ref JuMP.set_binary(::UserDecisionVariableRef)), and
[`set_integer`](@ref JuMP.set_integer(::UserDecisionVariableRef)).

### Type Specific
Finally, we consider methods unique to `InfiniteOpt` that exist to modify
specific variable types.

For finite variables, the parameter bounds can be modified via
[`@add_parameter_bounds`](@ref) and [`@set_parameter_bounds`](@ref) which
facilitate an intuitive symbolic syntax to add and/or overwrite existing
parameter bounds for a finite variable. For example, let's add the bounds
``t \in [0, 5]`` and ``x = 0`` to `c`:
```jldoctest var_macro
julia> @add_parameter_bounds(c, (t in [0, 5], x == 0))
```

A number of other techniques exist for the various variable types can be found in 
the manual below.

## Datatypes
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt, InfiniteOpt.Collections]
Order   = [:type]
```
```@docs
InfOptVariable
InfiniteVariable
PointVariable
ParameterBounds
FiniteVariable
VariableData
InfiniteVariableIndex
PointVariableIndex
FiniteVariableIndex
InfiniteVariableRef
PointVariableRef
FiniteVariableRef
InfiniteOpt.Collections.VectorTuple
```

## [Methods/Macros] (@id var_methods)
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@infinite_variable
@point_variable
@finite_variable
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Type{<:InfOptVariableType})
JuMP.add_variable(::InfiniteModel, ::InfOptVariable, ::String)
used_by_constraint(::DecisionVariableRef)
used_by_measure(::DecisionVariableRef)
used_by_objective(::DecisionVariableRef)
is_used(::DecisionVariableRef)
used_by_point_variable(::Union{InfiniteVariableRef, DerivativeRef})
used_by_semi_infinite_variable(::Union{InfiniteVariableRef, DerivativeRef})
is_used(::Union{InfiniteVariableRef, DerivativeRef})
JuMP.delete(::InfiniteModel, ::DecisionVariableRef)
JuMP.num_variables(::InfiniteModel, ::Type{InfOptVariable})
JuMP.all_variables(::InfiniteModel, ::Type{InfOptVariable})
JuMP.name(::DecisionVariableRef)
JuMP.set_name(::DecisionVariableRef, ::String)
JuMP.variable_by_name(::InfiniteModel, ::String)
JuMP.has_lower_bound(::UserDecisionVariableRef)
JuMP.lower_bound(::UserDecisionVariableRef)
JuMP.set_lower_bound(::UserDecisionVariableRef, ::Real)
JuMP.LowerBoundRef(::UserDecisionVariableRef)
JuMP.delete_lower_bound(::UserDecisionVariableRef)
JuMP.has_upper_bound(::UserDecisionVariableRef)
JuMP.upper_bound(::UserDecisionVariableRef)
JuMP.set_upper_bound(::UserDecisionVariableRef, ::Real)
JuMP.UpperBoundRef(::UserDecisionVariableRef)
JuMP.delete_upper_bound(::UserDecisionVariableRef)
JuMP.is_fixed(::UserDecisionVariableRef)
JuMP.fix_value(::UserDecisionVariableRef)
JuMP.fix(::UserDecisionVariableRef, ::Real; ::Bool)
JuMP.FixRef(::UserDecisionVariableRef)
JuMP.unfix(::UserDecisionVariableRef)
JuMP.start_value(::UserDecisionVariableRef)
JuMP.set_start_value(::UserDecisionVariableRef, ::Real)
start_value_function(::Union{InfiniteVariableRef, DerivativeRef})
set_start_value_function(::InfiniteVariableRef, ::Union{Real, Function})
reset_start_value_function(::InfiniteVariableRef)
JuMP.is_binary(::UserDecisionVariableRef)
JuMP.set_binary(::UserDecisionVariableRef)
JuMP.BinaryRef(::UserDecisionVariableRef)
JuMP.unset_binary(::UserDecisionVariableRef)
JuMP.is_integer(::UserDecisionVariableRef)
JuMP.set_integer(::UserDecisionVariableRef)
JuMP.IntegerRef(::UserDecisionVariableRef)
JuMP.unset_integer(::UserDecisionVariableRef)
parameter_refs(::InfiniteVariableRef)
parameter_list(::InfiniteVariableRef)
raw_parameter_refs(::InfiniteVariableRef)
@set_parameter_bounds
@add_parameter_bounds
has_parameter_bounds(::FiniteVariableRef)
parameter_bounds(::FiniteVariableRef)
set_parameter_bounds(::FiniteVariableRef, ::ParameterBounds{GeneralVariableRef})
add_parameter_bounds(::FiniteVariableRef, ::ParameterBounds{GeneralVariableRef})
delete_parameter_bounds(::FiniteVariableRef)
infinite_variable_ref(::PointVariableRef)
parameter_values(::PointVariableRef)
raw_parameter_values(::PointVariableRef)
JuMP.relax_integrality(::InfiniteModel)
```
