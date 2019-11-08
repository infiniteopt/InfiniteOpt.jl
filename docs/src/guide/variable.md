# Variables
A guide and manual for the definition and use of variables in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview
Decision variables are at the very core of `InfiniteOpt` as its name alludes
to mathematical programs that entail infinite decision spaces (i.e., contain
infinite decision variables). Principally, three variable types are employed:
infinite, point, and hold. Infinite variables encompass any decision variable
that is parameterized by an infinite parameter (e.g., space-time variables and
recourse variables). Point variables are infinite variables at a particular
infinite parameter value (point). Finally, hold variables are decisions that
are made irrespective of the infinite domain (e.g., first stage variables and
design variables). Or in other words they hold a particular value over the
infinite domain or some sub-domain of it.

## Basic Usage
Infinite, point, and hold variables are typically defined via their respective
macros: [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@hold_variable`](@ref). These macros generally emulate [`JuMP.@variable`](@ref)
except that they each employ additional syntax capabilities to employ their
respective variable type.

Let's first setup a simple space-time model with infinite parameters time `t` and
spacial position `x`:
```jldoctest var_basic
julia> using InfiniteOpt, JuMP

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_parameter(model, x[1:2] in [-1, 1], independent = true)
2-element Array{ParameterRef,1}:
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
3-element Array{InfiniteVariableRef,1}:
 w[1](t, x)
 w[2](t, x)
 w[3](t, x)
```
Thus we create a Julia array variable `w` whose elements `w[i]` point to their
respective infinite variables `w[i](t, x)` stored in `model`. Note that the `i`
used in the array definition can be used to index attributes assigned to each
variable in the array. In this case, we used `i` to assign different initial
guess values for each variable via the `start` keyword argument.

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

### Hold Variables
Finally, we can add hold variables to our model. These denote variables that
hold a single value over the infinite domain or some portion of it (e.g.,
design variables, first stage variables, etc.). Let's add a hold variable
``0 \leq d \leq 42`` that is an integer variable and defined over all infinite
domains (i.e., time and space):
```jldoctest var_basic
julia> @hold_variable(model, 0 <= d <= 42, Int)
d
```
This creates a Julia variable `d` that points to the hold variable `d` which has
a lower bound of 0, an upper bound of 42, and is an integer variable. Thus,
[`@hold_variable`](@ref) follows the same exact syntax as [`JuMP.@variable`](@ref)
except that it also allows the user to specify a subdomain over which the hold
variable is valid. For example, let's add a hold variable `z` that is only valid
over the subdomain ``t \in [0, 5]`` via the `parameter_bounds` keyword argument:
```jldoctest var_basic
julia> anon = @hold_variable(model, parameter_bounds = (t in [0, 5]),
                             base_name = "z")
z
```
Here we make an anonymous variable for the sake of example whose reference is
stored to the Julia variable `anon` and points to a hold variable `z` which is
only valid for ``t \in [0, 5]``. Thus, this will be enforced in any constraints that
involve `anon`, meaning they will automatically be bounded to such a subdomain.
Any number of parameters bounds (bounds on the parameters of the infinite
domain) can be added in a tuple like argument as explained in the documentation
for [`@hold_variable`](@ref).

Now we have defined variables that we can use in the objective, measures, and
constraints. Please note that the above tutorial only shows a small portion of
the capabilities and options available in defining variables. A full description
is provided in the documentation below.

## Variable Definition Methodology
The [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@hold_variable`](@ref) macros all follow a similar methodology behind the
scenes and these commonalities are discussed in this section for conciseness.
Defining/initializing a variable principally involves the following steps:
- Define the variable information pertaining to `JuMP.VariableInfo` (e.g., bounds, indicate if it is integer, etc.)
- Construct a concrete subtype of [`InfOptVariable`](@ref) to store the variable information
- Add the `InfOptVariable` object to an `InfiniteModel` and assign a name
- Create a concrete subtype of [`InfOptVariableRef`](@ref) that points to the variable object stored in the model

The `JuMP.VariableInfo` data structure stores the following variable information:
- `has_lb::Bool`: Specifies a `Bool` it has a lower bound
- `lower_bound::Number`: Specifies lower bound value
- `has_ub::Bool`: Specifies a `Bool` it has a upper bound
- `upper_bound::Number`: Specifies upper bound value
- `has_fix::Bool`: Specifies a `Bool` it is fixed
- `fixed_value::Number`: Specifies the fixed value
- `has_start::Bool`: Specifies a `Bool` it has a start value
- `start::Number`: Specifies the start value
- `binary`: Specifies `Bool` if it is binary
- `integer`: Specifies `Bool` if it is integer.
Thus, the user specifies this information to prepare such an object:
```jldoctest genvar_define; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> info = VariableInfo(true, 0, true, 42, false, 0, false, 0, false, true)
VariableInfo{Int64,Int64,Int64,Int64}(true, 0, true, 42, false, 0, false, 0, false, true)
```
Here we specified a lower bound of 0, an upper bound of 42, and that it is
integer.

The variable objects (`InfOptVariable` subtypes) are defined via
[`build_variable`](@ref) which requires that the user provide a
`JuMP.VariableInfo` object, the variable type to be defined (`Infinite`,
`Point`, or `Hold`), and any necessary keyword arguments required for that
variable type (i.e., `parameter_refs`, `infinite_variable_ref`, and/or
`parameter_values`). For example, let's build an infinite variable `y(t)` that
has an lower bound of 0, an upper bound of 42, and is integer:
```jldoctest genvar_define
julia> @infinite_parameter(model, t in [0, 10])
t

julia> info = VariableInfo(true, 0, true, 42, false, 0, false, 0, false, true);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = (t))
InfiniteVariable{Int64,Int64,Int64,Int64}(VariableInfo{Int64,Int64,Int64,Int64}(true, 0, true, 42, false, 0, false, 0, false, true), (t,))
```
Thus, we create an [`InfiniteVariable`](@ref) object with the desired properties.
Note that in this case the `parameter_refs` keyword argument is required to
indicate which infinite parameter(s) this infinite variable will depend on.

Once a variable has been built, it needs to be added to our `model` and a Julia
variable should be defined to reference it. Variables are added via
[`add_variable`](@ref JuMP.add_variable(::InfiniteModel, ::InfOptVariable))
which adds a variable object to the model, assigns a name to the variable,
adds any constraints associated with the `JuMP.VariableInfo`, and returns
an appropriate variable reference variable (subtyped from [`InfOptVariableRef`](@ref)).
For example, let's add `inf_var` to `model`:
```jldoctest genvar_define
julia> var_ref = add_variable(model, inf_var, "y")
y(t)
```
Thus, we have added an infinite variable `y` that is parameterized by `t` with the
variable information mentioned above and now have a [`InfiniteVariableRef`](@ref)
called `var_ref` that can be used in defining our infinite model.

## Macro Variable Definition
The [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@hold_variable`](@ref) macros automate the variable definition process
discussed above in the [Variable Definition Methodology](@ref) section
via a straightforward symbolic syntax. The only key difference
is that non-anonymous macro calls will register variable names to ensure they
are not repeated. Anonymous macro calls forgo this step and exactly follow the
process described above. This section will highlight the details of using these
macros.

### General Usage
Here we discuss the features that the variable macros have in common (generally
these pertain to `JuMP`-like features). To illustrate this via example, let's
setup a model with a variety of infinite parameters ``t \in [0,10]``,
``x \in [-1, 1]^3``, and ``\xi \in \mathcal{N}(0, 1)``:
```jldoctest var_macro
julia> using InfiniteOpt, JuMP, Distributions

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_parameter(model, x[1:3] in [-1, 1], independent = true);

julia> @infinite_parameter(model, ξ in Normal());
```

We will first consider anonymous variable macro calls which generally are less
convenient than non-anonymous macro calls which offer a much more intuitive
mathematical syntax. However, anonymous variables can be useful and provide a
good foundation to understanding non-anonymous variables.
Furthermore, we'll use hold variables as the motivating examples since they
best exemplify commonalities between the macros. First, let's consider single
anonymous definition a hold variable:
```jldoctest var_macro
julia> var_ref = @hold_variable(model)
noname
```
Here we just added a nameless hold variable to `model` and defined `var_ref`
as a [`HoldVariableRef`](@ref) that points to it. We can add a name via the
`base_name` keyword argument:
```jldoctest var_macro
julia> var_ref1 = @hold_variable(model, base_name = "d")
d

julia> var_ref2 = @hold_variable(model, base_name = "d")
d
```
Now we've made 2 more hold variables both called `d`. Thus, the anonymous syntax
allows us to define variables with the same name. Moreover, any variable
information can be specified via the appropriate keywords which include:
- `lower_bound::Number`: specifies lower bound
- `upper_bound::Number`: specifies upper bound
- `start::Number`: specifies the initial guess value the solver will use
- `binary::Bool`: specified if is binary variable
- `integer::Bool`: specifies if is integer variable.
Anonymous variables must use these keyword arguments since symbolic definition
is only permitted for non-anonymous macro calls. For example, let's define a
hold variable ``0 \leq d \leq 5`` that is integer:
```jldoctest var_macro
julia> var_ref = @hold_variable(model, base_name = "d", lower_bound = 0,
                                upper_bound = 5, integer = true)
d
```

We can also define arrays of variables using any indices of our choice. For
example, let's define a 3 dimensional vector with indices `[1, 2, 3]`:
```jldoctest var_macro
julia> var_refs = @hold_variable(model, [i = 1:3], start = [0, 2, 1][i],
                                 base_name = "d")
3-element Array{HoldVariableRef,1}:
 d[1]
 d[2]
 d[3]
```
Thus, we define 3 variables named `d[i]` and each with different start values
and define `var_refs` which is a vector of `HoldVariableRef`s that uses the
indices we specified. Note the syntax `i = indices` is used to define an
iteration variable to use with the keyword arguments to assign different values
for each variable being defined. Note the above example is equivalent to:
```jldoctest var_macro
julia> starts = [0, 2, 1];

julia> var_refs = Vector{HoldVariableRef}(undef, 3);

julia> for i = eachindex(var_refs)
          var_refs[i] = @hold_variable(model, base_name = "d", start = starts[i])
       end
```
Other non-standard indices can also be used such as the following examples:
```jldoctest var_macro
julia> var_refs2 = @hold_variable(model, [2:4], base_name = "d")
1-dimensional DenseAxisArray{HoldVariableRef,1,...} with index sets:
    Dimension 1, 2:4
And data, a 3-element Array{HoldVariableRef,1}:
 d[2]
 d[3]
 d[4]

julia> var_refs3 = @hold_variable(model, [[:A, :C, :Z]], base_name = "d")
1-dimensional DenseAxisArray{HoldVariableRef,1,...} with index sets:
    Dimension 1, Symbol[:A, :C, :Z]
And data, a 3-element Array{HoldVariableRef,1}:
 d[A]
 d[C]
 d[Z]

julia> var_refs3 = @hold_variable(model, [i=1:2, j=i:2], base_name = "d")
JuMP.Containers.SparseAxisArray{HoldVariableRef,2,Tuple{Any,Any}} with 3 entries:
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

julia> var_refs1 = @hold_variable(model, [a:b], base_name = "d")
1-dimensional DenseAxisArray{HoldVariableRef,1,...} with index sets:
    Dimension 1, 1:3
And data, a 3-element Array{HoldVariableRef,1}:
 d[1]
 d[2]
 d[3]

julia> var_refs2 = @hold_variable(model, [a:b], base_name = "d", container = Array)
3-element Array{HoldVariableRef,1}:
 d[1]
 d[2]
 d[3]
```

Now that we have a foundation with anonymous variable macro calls, let's focus
on non-anonymous calls which offer a much more straightforward syntax.

TO BE CONTINUED...

### Infinite Variables
Infinite variables entail decision variables that depend on infinite parameter(s).
Thus, [`@infinite_variable`](@ref) follows the general definition methodology
with this additional consideration.

Let's first consider a basic anonymous definition of an infinite variable
``y(t, x)``:
```jldoctest var_macro
julia> yref = @infinite_variable(model, parameter_refs = (t, x), base_name = "y")
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
julia> yref2 = @infinite_variable(model, parameter_refs = (t), base_name = "y")
y(t)
```

### Point Variables


### Hold Variables


## Manipulation


## Datatypes
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
InfOptVariable
InfiniteVariable
PointVariable
ParameterBounds
HoldVariable
InfOptVariableRef
GeneralVariableRef
MeasureFiniteVariableRef
FiniteVariableRef
InfiniteVariableRef
PointVariableRef
HoldVariableRef
```

## Methods/Macros
```@index
Pages   = ["variable.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
JuMP.build_variable(::Function, ::JuMP.VariableInfo, ::Symbol)
JuMP.add_variable(::InfiniteModel, ::InfOptVariable, ::String)
@infinite_variable
@point_variable
@hold_variable
JuMP.owner_model(::GeneralVariableRef)
JuMP.index(::GeneralVariableRef)
used_by_constraint(::InfOptVariableRef)
used_by_measure(::InfOptVariableRef)
used_by_objective(::InfOptVariableRef)
is_used(::InfOptVariableRef)
used_by_point_variable(::InfiniteVariableRef)
used_by_reduced_variable(::InfiniteVariableRef)
is_used(::InfiniteVariableRef)
JuMP.delete(::InfiniteModel, ::InfOptVariableRef)
JuMP.is_valid(::InfiniteModel, ::InfOptVariableRef)
JuMP.num_variables(::InfiniteModel)
JuMP.has_lower_bound(::InfOptVariableRef)
JuMP.lower_bound(::InfOptVariableRef)
JuMP.set_lower_bound(::InfOptVariableRef, ::Number)
JuMP.LowerBoundRef(::InfOptVariableRef)
JuMP.delete_lower_bound(::InfOptVariableRef)
JuMP.has_upper_bound(::InfOptVariableRef)
JuMP.upper_bound(::InfOptVariableRef)
JuMP.set_upper_bound(::InfOptVariableRef, ::Number)
JuMP.UpperBoundRef(::InfOptVariableRef)
JuMP.delete_upper_bound(::InfOptVariableRef)
JuMP.is_fixed(::InfOptVariableRef)
JuMP.fix_value(::InfOptVariableRef)
JuMP.fix(::InfOptVariableRef, ::Number; ::Bool)
JuMP.FixRef(::InfOptVariableRef)
JuMP.unfix(::InfOptVariableRef)
JuMP.start_value(::InfOptVariableRef)
JuMP.set_start_value(::InfOptVariableRef, ::Number)
JuMP.is_binary(::InfOptVariableRef)
JuMP.set_binary(::InfOptVariableRef)
JuMP.BinaryRef(::InfOptVariableRef)
JuMP.unset_binary(::InfOptVariableRef)
JuMP.is_integer(::InfOptVariableRef)
JuMP.set_integer(::InfOptVariableRef)
JuMP.IntegerRef(::InfOptVariableRef)
JuMP.unset_integer(::InfOptVariableRef)
JuMP.name(::InfOptVariableRef)
JuMP.set_name(::InfiniteVariableRef, ::String)
JuMP.set_name(::PointVariableRef, ::String)
JuMP.set_name(::HoldVariableRef, ::String)
parameter_refs(::InfiniteVariableRef)
set_parameter_refs(::InfiniteVariableRef, ::Tuple)
add_parameter_ref(::InfiniteVariableRef,::Union{ParameterRef, AbstractArray{<:ParameterRef}})
@set_parameter_bounds
@add_parameter_bounds
has_parameter_bounds(::HoldVariableRef)
parameter_bounds(::HoldVariableRef)
set_parameter_bounds(::HoldVariableRef, ::ParameterBounds)
add_parameter_bound(::HoldVariableRef, ::ParameterRef, ::Number, ::Number)
delete_parameter_bound(::HoldVariableRef, ::ParameterRef)
delete_parameter_bounds(::HoldVariableRef)
infinite_variable_ref(::PointVariableRef)
parameter_values(::PointVariableRef)
JuMP.variable_by_name(::InfiniteModel, ::String)
JuMP.all_variables(::InfiniteModel)
```
