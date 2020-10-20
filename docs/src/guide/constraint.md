# [Constraints] (@id constr_page)
A guide and manual for defining and manipulating constraints in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.

## Overview
Constraints are an integral part of infinite dimensional problems and serve
as a fundamental utility of `InfiniteOpt`. In particular, `InfiniteOpt` supports
finite constraints that entail finite variables and/or measures that fully remove any
infinite parameter dependencies (e.g., first stage constraints), infinite
constraints that are enforced over the entire domain of its infinite parameter
dependencies (e.g., path constraints), and bounded infinite constraints which
are enforced over some specified sub-domain of its infinite parameter
dependencies (e.g., boundary conditions). This page will highlight how to
implement these types of constraints in `InfiniteOpt`.

## Basic Usage
Principally, the [`@constraint`](@ref) and [`@BDconstraint`](@ref) macros should
be employed to specify constraints. Here `@constraint` is used to specify
typical finite and infinite constraints and `@BDconstraint` is used to specify
bounded constraints (i.e., infinite constraints with a constrained sub-domain
of its full infinite domain). First, let's setup an infinite model with
variables that we can add constraints to:
```jldoctest constrs; setup = :(using InfiniteOpt, JuMP, Distributions)
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_parameter(model, x[1:2] in Normal());

julia> @infinite_variable(model, T(t, x));

julia> @infinite_variable(model, g(t));

julia> @hold_variable(model, z[1:2]);

julia> @hold_variable(model, w, parameter_bounds = (t in [0, 5]));
```

### Unbounded Constraints
Unbounded constraints don't entail user-defined curtailing of dependent infinite
parameter domains and are defined using [`@constraint`](@ref). This is done
following typical `JuMP` syntax. Anonymous constraints (constraints without
associated register names) are defined via 2 arguments where the first is the
infinite model and the second is the constraint expression. For example, let's
define the constraint ``2T(t, x) + ||z||^2 \leq 0``:
```jldoctest constrs
julia> cref = @constraint(model, 2T + sum(z[i]^2 for i = 1:2) <= 0)
z[1]² + z[2]² + 2 T(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ~ Normal, x[2] ~ Normal
```
Thus, we added an infinite constraint (which infinite with respect to `t` and `x`)
to `model` and stored the corresponding constraint reference to `cref`. The
allowed constraint operators are `==`, `<=`, `≤`, `>=`, and `≥`.

!!! note
    Currently, `InfiniteOpt` only supports scalar constraint expressions. Other
    types such as vector constraints are not currently extended.

Similarly, we can define an array of constraints with varied indexes by including
an additional argument before the constraint expression. For example,
let's define ``3z[i] - 14 == 0, \forall i \in \{1,2\}``:
```jldoctest constrs
julia> crefs = @constraint(model, [i = 1:2], 3z[i] - 14 == 0)
2-element Array{InfOptConstraintRef,1}:
 3 z[1] = 14.0
 3 z[2] = 14.0
```
Thus, we added two constraints to `model` and stored a vector of the corresponding
constraint references to the `Julia` variable `crefs`.

Named constraints are defined by including a name as part of the second argument.
For example, let's add the constraint ``\int_0^10 g(t) dt == 4``:
```jldoctest constrs
julia> @constraint(model, measure_constr, integral(g, t) == 4)
measure_constr : integral{t ∈ [0, 10]}[g(t)] = 4.0
```
Thus, we added another constraint named `measure_constr` and created a `Julia`
variable `measure_constr` that contains a reference to that constraint.

!!! note
    Linear algebra constraints can also be used when defining constraints
    when `.` is added in front of the constraint operators (e.g., `.<=`). This
    behavior is further explained in `JuMP`'s documentation 
    [here](https://jump.dev/JuMP.jl/stable/constraints/#Vectorized-constraints-1). 
    However, note that that vector array sets such as `MOI.Nonnegatives` are not 
    currently supported.

### Bounded Constraints
Bounded constraints denote constraints with a particular sub-domain of an infinite
domain. Such constraints might pertain to point constraints evaluated at a
particular infinite parameter values, constraints limited to a certain range of
infinite parameter values, and/or entail bounded hold variables whose
sub-domains are also taken into consideration. These types of constraints are
defined using [`@BDconstraint`](@ref). For example, let's add the initial
condition ``g(0) == 0``:
```jldoctest constrs
julia> @BDconstraint(model, initial(t == 0), g == 0)
initial : g(t) = 0.0, ∀ t = 0
```
Thus, we have added a constraint to `model` defined over the sub-domain ``t = 0``
in accordance with the initial condition. This is referred to as a bounded
constraint named `initial` where a `Julia` variable `initial` has been defined to
store a reference to that constraint.

More complex sub-domains can be specified by simply adding more conditions in
the second argument. To illustrate this, let's define an anonymous constraint
for ``2T^2(t, x) + w \geq 3, \ \forall t = 0, \ x \in [-1, 1]^2``:
```jldoctest constrs
julia> cref = @BDconstraint(model, (t == 0, x in [-1, 1]), 2T^2 + w >= 3)
2 T(t, x)² + w ≥ 3.0, ∀ t = 0, x[1] ~ Normal ∩ [-1, 1], x[2] ~ Normal ∩ [-1, 1]
```
where `cref` contains the corresponding constraint reference.

Now we have added constraints to our model and it is ready to be transcripted!

## Data Structure
Here we detail the data structures used to store constraints in `InfiniteOpt`.
In general, constraints in `JuMP` are of the form: `[function] in [set]` where
`function` corresponds to a `JuMP` expression and `set` corresponds to a `MOI`
set. Since `InfiniteOpt` only supports scalar constraints currently, expressions
must be inherited from `JuMP.AbstractJuMPScalar` and the supported `MOI` sets are
`MOI.EqualTo`, `MOI.LessThan`, `MOI.GreaterThan`, `MOI.Integer`, and `MOI.ZeroOne`.
Where the last 2 are intended for single variable constraints.  

Furthermore, constraints in the form mentioned above are stored in appropriate
constraint object inherited from `JuMP.AbstractConstraint`. Typical scalar
constraints use [`JuMP.ScalarConstraint`](@ref) which has the fields:
- `func::JuMP.AbstractJuMPScalar` The constraint expression
- `set::MOI.AbstractScalarSet` The MOI set.

Similarly, bounded constraints are stored in [`BoundedScalarConstraint`](@ref)
which has `bounds` and `orig_bounds` fields in addition to `func` and `set`.
These additional fields store the current constraint bounds and the bounds
originally given, respectively. This distinction is needed to facilitate
deletion methods such as deleting a bounded hold variable.

These constraint objects are what store constraints in `InfiniteModel`s. And
these are referred to by an [`InfOptConstraintRef`](@ref).

## Definition
In this section, we describe the ins and outs of defining constraints. Note that
this process is analogous to the manner in which variables are defined and added
to the model.

### Manual Definition
The [`@constraint`](@ref) and [`@BDconstraint`](@ref) both follow a similar
methodology behind the scenes and these commonalities are discussed in this
section for conciseness. Defining/initializing a constraint principally involves
the following steps:
- Define the constraint information (i.e., function, set, and parameter bounds)
- Construct a concrete subtype of `JuMP.AbstractConstraint` to store the constraint information
- Add the `AbstractConstraint` object to an `InfiniteModel` and assign a name
- Create an [`InfOptConstraintRef`](@ref) that points to the constraint object stored in the model.

The constraint objects are specified via
[`JuMP.build_constraint`](@ref) which requires that the user provides
a `JuMP.AbstractJuMPScalar`, a `MOI.AbstractScalarSet`, and any keyword
arguments such as `ParameterBound`. For example, let's build a scalar constraint
for ``3T(t, x) - g^2(t) \leq 0``:
```jldoctest constrs; setup = :(using MathOptInterface; const MOI = MathOptInterface)
julia> constr = build_constraint(error, 3T - g^2, MOI.LessThan(0.0))
ScalarConstraint{GenericQuadExpr{Float64,GeneralVariableRef},MathOptInterface.LessThan{Float64}}(-g(t)² + 3 T(t, x), MathOptInterface.LessThan{Float64}(0.0))
```

Now the built constraint object can be added to the infinite model via
[`add_constraint`](@ref JuMP.add_constraint(::InfiniteModel, ::JuMP.AbstractConstraint, ::String)).
Let's do so with our example and assign it the name of `c1` (note that adding a
name is optional):
```jldoctest constrs
julia> cref = add_constraint(model, constr, "c1")
c1 : -g(t)² + 3 T(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ~ Normal, x[2] ~ Normal
```

Thus, we have made our constraint and added it `model` and now have a constraint
reference `cref` that we can use to access it.

The [`@constraint`](@ref) and [`@BDconstraint`](@ref) automate the above steps
where `@BDconstraint` is specifically used to implement the `parameter_bounds`
keyword argument in `build_constraint`.

### Macro Definition
As mentioned above in the Basic Usage section, the [`@constraint`](@ref) macro
should be used for defining all constraints that do not entail specified
parameter bounds. In general, the syntax follows the form
`@constraint([InfiniteModel], [name][indexing expr], [scalar constr expr])`.
The second argument is optional and is used to assign a name and/or define
indexing variables to be used in the constraint expr. When a name is provided it
is registered and cannot be used again for another constraint or variable name.
The indexing expression can be used to produce an array of constraints as shown
below (notice this is equivalent to looping over individual `@constraint` calls):
```jldoctest constrs
julia> crefs = @constraint(model, [i = 1:2], 2z[i] - g == 0)
2-element Array{InfOptConstraintRef,1}:
 2 z[1] - g(t) = 0.0, ∀ t ∈ [0, 10]
 2 z[2] - g(t) = 0.0, ∀ t ∈ [0, 10]

julia> crefs = Vector{InfOptConstraintRef}(undef, 2);

julia> for i = 1:2
           crefs[i] = @constraint(model, 2z[i] - g == 0)
       end

julia> crefs
2-element Array{InfOptConstraintRef,1}:
 2 z[1] - g(t) = 0.0, ∀ t ∈ [0, 10]
 2 z[2] - g(t) = 0.0, ∀ t ∈ [0, 10]
```
Again, please note that only scalar constraints are currently supported and thus
the `[scalar constr expr]` must be scalar.

An interesting corollary is that `@constraint` in certain cases will produce
bounded constraints even if the user doesn't specify any parameter bounds. This,
occurs when a bounded hold variable is included. For example, if we want to
make the constraint ``g(t) + 2.5w \leg 2`` it would only be valid over the sub-domain
``t \in [0, 5] \subsetneq [0, 10]`` due to the restrictions on `w`. Thus,
`@constraint` automatically takes this into account and produces a bounded
constraint:
```jldoctest constrs
julia> @constraint(model, bounded_example, g + 2.5w <= 2)
bounded_example : g(t) + 2.5 w ≤ 2.0, ∀ t ∈ [0, 5]
```
If this behavior is not desired, the parameter bounds associated with the
variable should be deleted and then should be managed individually for each
constraint.

The [`@BDconstraint`](@ref) is very similar except that it adds the capability
of symbolically specifying parameter bounds. Thus the syntax is of the form
`@BDconstraint([InfiniteModel], [name][indexing expr](bound expr), [scalar constr expr])`.
Note that here the `name` and `[indexing expr]` are optional, but the
`(param bounds)` tuple is required. Thus, three arguments must always be given.
The `(bound expr)` can be of the form:
- `(param in [lb, ub], ...)` enforcing `param` to be in a sub-domain from `lb`
                             to `ub` (note `∈` can be used in place of `in`)
- `(params in [lb, ub], ...)` enforcing that all parameter references in `params`
                              each be a in sub-domain from `lb` to `ub`
- `(lb <= param <= ub, ...)` enforcing `param` to be in a sub-domain from `lb`
                             to `ub`
- `(lb <= params <= ub, ...)` enforcing that all parameter references in `params`
                              each be a in sub-domain from `lb` to `ub`
- `(param == value, ...)` enforcing `param` to be equal to `value`
- `(params == value, ...)` enforcing that all parameter references in `params`
                            each to be equal to `value`
- Any combination of the above forms. Must be inside parentheses and comma
  separated.
Please refer to the Basic Usage section and the manual for examples.

Finally, the [`@constraint`](@ref) and [`@BDconstraint`](@ref) macros allow
the user to specify the `container` keyword argument when defining an array
of constraints. For example, we can force a group of bounded constraint
references to be stored in a `JuMP.SparseAxisArray`:
```jldoctest constrs
julia> @BDconstraint(model, [i = 1:2](x[i] in [0, 1]), T^2 + z[i] <= 1,
                     container = SparseAxisArray)
JuMP.Containers.SparseAxisArray{InfOptConstraintRef,1,Tuple{Int64}} with 2 entries:
  [2]  =  T(t, x)² + z[2] ≤ 1.0, ∀ t ∈ [0, 10], x[1] ~ Normal, x[2] ~ Normal ∩ [0, 1] 
  [1]  =  T(t, x)² + z[1] ≤ 1.0, ∀ t ∈ [0, 10], x[1] ~ Normal ∩ [0, 1], x[2] ~ Normal
```
For more information on `JuMP` containers please visit their page
[here](http://www.juliaopt.org/JuMP.jl/stable/containers/).

## Queries
In this section, we describe a variety of methods to extract constraint
information.

### Basic
A number of constraint properties can be extracted via constraint references.
Principally, the validity, name, model, index, and constraint object can be queried
via [`is_valid`](@ref JuMP.is_valid(::InfiniteModel, ::InfOptConstraintRef)),
[`name`](@ref JuMP.name(::InfOptConstraintRef)),
[`owner_model`](@ref JuMP.owner_model(::InfOptConstraintRef)),
[`index`](@ref JuMP.index(::InfOptConstraintRef)),
and [`constraint_object`](@ref JuMP.constraint_object(::InfOptConstraintRef)),
respectively. These methods all constitute extensions of `JuMP` methods and
follow exactly the same behavior. Let's try them out with the following example:
```jldoctest constrs
julia> is_valid(model, measure_constr) # check if contained in model
true

julia> name(measure_constr) # get the name
"measure_constr"

julia> m = owner_model(measure_constr); # get the model it is added to

julia> index(measure_constr) # get the constraint's index
ConstraintIndex(4)

julia> constraint_object(measure_constr) # get the raw constraint datatype
ScalarConstraint{GenericAffExpr{Float64,GeneralVariableRef},MathOptInterface.EqualTo{Float64}}(integral{t ∈ [0, 10]}[g(t)], MathOptInterface.EqualTo{Float64}(4.0))
```

Also, [`constraint_by_name`](@ref JuMP.constraint_by_name(::InfiniteModel, ::String))
can be used to retrieve a constraint reference if only the name is known and its
name is unique. For example, let's extract the reference for `"c1"`:
```jldoctest constrs
julia> cref = constraint_by_name(model, "c1")
c1 : -g(t)² + 3 T(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ~ Normal, x[2] ~ Normal
```

### Parameter Bounds
As explained above, bounded constraints serve as an integral capability of
`InfiniteOpt`. Information about parameter bounds can be obtained via
[`has_parameter_bounds`](@ref has_parameter_bounds(::InfOptConstraintRef)) and
[`parameter_bounds`](@ref parameter_bounds(::InfOptConstraintRef)) which indicate
if a constraint is bounded and what its [`ParameterBounds`](@ref) are,
respectively. These are exemplified below:
```jldoctest constrs
julia> has_parameter_bounds(measure_constr) # check if constraint is bounded
false

julia> has_parameter_bounds(initial)
true

julia> parameter_bounds(initial)
Subdomain bounds (1): t = 0
```
Note that `parameter_bounds` will error if the constraint is not bounded.

### General
Constraints can be defined in a number of ways symbolically that differ from
how it is actually stored in the model. This principally occurs since like terms
and constants are combined together where possible with the variable terms on the
left hand side and the constant on the right hand side. For instance, the
constraint ``2g(t) + 3g(t) - 2 \leq 1 + z_1`` would be normalized ``5g(t) - z_1 \leq 3``. In
accordance with this behavior [`normalized_rhs`](@ref JuMP.normalized_rhs(::InfOptConstraintRef))
and [`normalized_coefficient`](@ref JuMP.normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef))
can be used to query the normalized right hand side and the coefficient of a
particular variable reference, respectively. Let's employ the above example to
illustrate this:
```jldoctest constrs
julia> @constraint(model, constr, 2g + 3g - 2 <= 1 + z[1])
constr : 5 g(t) - z[1] ≤ 3.0, ∀ t ∈ [0, 10]

julia> normalized_rhs(constr)
3.0

julia> normalized_coefficient(constr, g)
5.0
```

There also exist a number of methods for querying an infinite model about what
constraints it contains.
[`list_of_constraint_types`](@ref JuMP.list_of_constraint_types(::InfiniteModel))
can be used query what types of constraints have been added to a model. This
is provided as a list of tuples where the first element is the expression type
and the second element is the set type (recall that constraints are stored in
the form `func-in-set`). Thus, for our current model we obtain:
```julia-repl
julia> list_of_constraint_types(model)
4-element Array{Tuple{DataType,DataType},1}:
 (GenericQuadExpr{Float64,GeneralVariableRef}, MathOptInterface.LessThan{Float64})
 (GenericQuadExpr{Float64,GeneralVariableRef}, MathOptInterface.GreaterThan{Float64})
 (GenericAffExpr{Float64,GeneralVariableRef}, MathOptInterface.LessThan{Float64})
 (GenericAffExpr{Float64,GeneralVariableRef}, MathOptInterface.EqualTo{Float64})
```
This information is useful when in combination with the
[`num_constraints`](@ref JuMP.num_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet}))
and [`all_constraints`](@ref JuMP.all_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet}))
methods which can take the expression type and/or the set type as inputs. Here
`num_constraints` provides the number of constraints that match a certain type 
and `all_constraints` returns a list of constraint references matching the criteria
provided. These  have been extended beyond `JuMP` functionality such additional
methods have been provided for the cases in which one wants to query solely off
of set or off expression type. Let's illustrate this with `num_constraints`:
```jldoctest constrs
julia> num_constraints(model) # total number of constraints
15

julia> num_constraints(model, GenericQuadExpr{Float64,GeneralVariableRef})
5

julia> num_constraints(model, MathOptInterface.LessThan{Float64})
6

julia> num_constraints(model, GenericQuadExpr{Float64,GeneralVariableRef},
                       MathOptInterface.LessThan{Float64})
4                   
```

## Modification
In this section, we highlight a number of methods that can be used to modify
existing constraints.

### Deletion
All constraints in `InfiniteOpt` can be removed in like manner to typical `JuMP`
constraints with the appropriate extension of
[`delete`](@ref JuMP.delete(::InfiniteModel, ::InfOptConstraintRef)). This will
remove the corresponding constraint object from the model. However, please note
any registered names will remain registered in the infinite model. This means
that a constraint with a registered name cannot be repeatedly added and removed
using the same name. To exemplify this, let's delete the constraint `c1`:
```jldoctest constrs
julia> cref = constraint_by_name(model, "c1")
c1 : -g(t)² + 3 T(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ~ Normal, x[2] ~ Normal

julia> delete(model, cref)
```

### General
There also are a number of ways to modify information and characteristics of
constraints. First, [`set_name`](@ref JuMP.set_name(::InfOptConstraintRef, ::String))
can be used to specify a new name for a particular constraint. For instance,
let's update the name of `initial` to `"init_cond"`:
```jldoctest constrs
julia> set_name(initial, "init_cond")

julia> initial
init_cond : g(t) = 0.0, ∀ t = 0
```

We can also update the normalized right hand side constant value or normalized
left hand side variable coefficient value using
[`set_normalized_rhs`](@ref JuMP.set_normalized_rhs(::InfOptConstraintRef, ::Real))
and [`set_normalized_coefficient`](@ref JuMP.set_normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef, ::Real))
, respectively. Let's again consider the constraint ``5g(t) - z_1 \leq 3`` as an
example. Let's change the constant term to -1 and the `g(t)` coefficient to 2.5:
```jldoctest constrs
julia> set_normalized_rhs(constr, -1)

julia> set_normalized_coefficient(constr, g, 2.5)

julia> constr
constr : 2.5 g(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 10]
```

!!! note
    In some cases, it may be more convenient to dynamically modify coefficients
    and other values via the use of finite parameters. This provides an avenue
    to update parameters without having to be concerned about the normalized form.
    For more information, see the [Finite Parameters](@ref finite_param_docs) page.

### Parameter Bounds
Parameter bounds can be added to, modified, or removed from any constraint in
`InfiniteOpt`. Principally, this is accomplished via [`@add_parameter_bounds`](@ref),
[`@set_parameter_bounds`](@ref),
and [`delete_parameter_bounds`](@ref delete_parameter_bounds(::InfOptConstraintRef))
in like manner to hold variables.

First, parameter bounds can be added to a constraint in an intuitive symbolic
syntax via [`@add_parameter_bounds`](@ref) which follows form
`@add_parameter_bounds(ref, bound_expr)` where `(bound_expr)` can be of the form:

- `(param in [lb, ub], ...)` enforcing `param` to be in a sub-domain from `lb`
                             to `ub` (note `∈` can be used in place of `in`)
- `(params in [lb, ub], ...)` enforcing that all parameter references in `params`
                              each be a in sub-domain from `lb` to `ub`
- `(lb <= param <= ub, ...)` enforcing `param` to be in a sub-domain from `lb`
                             to `ub`
- `(lb <= params <= ub, ...)` enforcing that all parameter references in `params`
                              each be a in sub-domain from `lb` to `ub`
- `(param == value, ...)` enforcing `param` to be equal to `value`
- `(params == value, ...)` enforcing that all parameter references in `params`
                            each to be equal to `value`
- Any combination of the above forms. Must be inside parentheses and comma
  separated.

For example, let's add the bound ``t \in [0, 1]`` to `constr`:
```jldoctest constrs
julia> @add_parameter_bounds(constr, (t in [0, 1]))

julia> constr
constr : 2.5 g(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 1]
```

In similar manner, [`@set_parameter_bounds`](@ref) can be employed to specify
what bounds a constraint has (overwriting any existing ones if forced). It follows
the same syntax, so let's use it to change the bounds on `t` to ``t = 0``:
```jldoctest constrs
julia> @set_parameter_bounds(constr, (t == 0), force = true)

julia> constr
constr : 2.5 g(t) - z[1] ≤ -1.0, ∀ t = 0
```

!!! note
    Constraint parameters bounds are intersected with those of any hold variables
    that are part of the constraint. If no such intersection exists then an error is
    thrown.

Finally, constraint bounds can be deleted via
[`delete_parameter_bounds`](@ref delete_parameter_bounds(::InfOptConstraintRef))
which deletes all parameter bounds. Again, note the parameter bounds associated
with hold variables will be unaffected and can only be removed by deleting them
from the variables directly. Now let's delete the parameter bounds associated
with our example:
```jldoctest constrs
julia> delete_parameter_bounds(constr)

julia> constr
constr : 2.5 g(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 10]
```

## Datatypes
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
BoundedScalarConstraint
ConstraintData
ConstraintIndex
InfOptConstraintRef
```

## Methods/Macros
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@BDconstraint
JuMP.build_constraint(::Function, ::Union{JuMP.GenericAffExpr{Any, <:GeneralVariableRef}, JuMP.GenericQuadExpr{Any, <:GeneralVariableRef}}, ::MOI.AbstractScalarSet)
JuMP.add_constraint(::InfiniteModel, ::JuMP.AbstractConstraint, ::String)
JuMP.owner_model(::InfOptConstraintRef)
JuMP.index(::InfOptConstraintRef)
JuMP.constraint_object(::InfOptConstraintRef)
JuMP.name(::InfOptConstraintRef)
JuMP.set_name(::InfOptConstraintRef, ::String)
JuMP.is_valid(::InfiniteModel, ::InfOptConstraintRef)
JuMP.delete(::InfiniteModel, ::InfOptConstraintRef)
parameter_refs(::InfOptConstraintRef)
has_parameter_bounds(::InfOptConstraintRef)
parameter_bounds(::InfOptConstraintRef)
set_parameter_bounds(::InfOptConstraintRef, ::ParameterBounds{GeneralVariableRef})
add_parameter_bounds(::InfOptConstraintRef, ::ParameterBounds{GeneralVariableRef})
delete_parameter_bounds(::InfOptConstraintRef)
JuMP.set_normalized_rhs(::InfOptConstraintRef, ::Real)
JuMP.normalized_rhs(::InfOptConstraintRef)
JuMP.add_to_function_constant(::InfOptConstraintRef, ::Real)
JuMP.set_normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef, ::Real)
JuMP.normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef)
JuMP.constraint_by_name(::InfiniteModel, ::String)
JuMP.list_of_constraint_types(::InfiniteModel)
JuMP.num_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet})
JuMP.num_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar})
JuMP.num_constraints(::InfiniteModel, ::Type{<:MOI.AbstractSet})
JuMP.num_constraints(::InfiniteModel)
JuMP.all_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet})
JuMP.all_constraints(::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar})
JuMP.all_constraints(::InfiniteModel, ::Type{<:MOI.AbstractSet})
JuMP.all_constraints(::InfiniteModel)
```
