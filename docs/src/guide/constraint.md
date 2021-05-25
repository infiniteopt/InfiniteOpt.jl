```@meta
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", 
                  r"MathOptInterface|MOI", r" for all | ∀ "]
```

# [Constraints] (@id constr_page)
A guide and manual for defining and manipulating constraints in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.

## Overview
Constraints are a key part of infinite dimensional problems and serve
as a fundamental utility of `InfiniteOpt`. In particular, `InfiniteOpt` supports
finite constraints that entail finite variables and/or measures that fully remove any
infinite parameter dependencies (e.g., first stage constraints), infinite
constraints that are enforced over the entire domain of its infinite parameter
dependencies (e.g., path constraints), and restricted constraints which
are enforced over some specified sub-domain of its infinite parameter
dependencies (e.g., boundary conditions). This page will highlight how to
implement these types of constraints in `InfiniteOpt`.

!!! note 
    Nonlinear constraints as defined by `JuMP.@NLconstraint` are not currently 
    supported by `InfiniteOpt`. See [Nonlinear Expressions](@ref) for more 
    information and possible workarounds. 

## Basic Usage
Principally, the 
[`@constraint`](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.@constraint) 
macro is used to define constraints. First, let's setup an infinite model with
variables that we can add constraints to:
```jldoctest constrs; setup = :(using InfiniteOpt)
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_parameter(model, x[1:2] in [-2, 2]);

julia> @variable(model, ya, Infinite(t, x));

julia> @variable(model, yb, Infinite(t));

julia> @variable(model, z[1:2]);
```

!!! note
    Unlike previous versions, `InfiniteOpt` now supports all of the constraints 
    offered by `JuMP`, including vector and semi-definite constraints! Please 
    see [JuMP's constraint documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/#Constraints) 
    for a thorough explanation of the supported types and syntax.

### Scalar Constraints
Scalar constraints use scalar functions of variables. For example, let's define 
the constraint 
``||z||^2 + 2y_a(t, x) \leq 0, \ \forall t \in [0, 10], x \in [-2, 2]^2`` 
using `@constraint`:
```jldoctest constrs
julia> @constraint(model, c1, sum(z[i]^2 for i = 1:2) + 2ya <= 0)
c1 : z[1]² + z[2]² + 2 ya(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ∈ [-2, 2], x[2] ∈ [-2, 2]
```
Thus, we added an infinite constraint (which infinite with respect to `t` and `x`)
to `model` and stored the corresponding constraint reference to `c1`. Note that 
this is enforced over the full infinite domains of the infinite parameters `t` 
and `x` which are implicitly used by `c1`. For scalar constraints like this one, 
the allowed constraint operators are `==`, `<=`, `≤`, `>=`, and `≥`.

!!! note
    Linear algebra constraints can also be used when defining constraints
    when `.` is added in front of the constraint operators (e.g., `.<=`). This
    behavior is further explained in 
    [`JuMP`'s constraint documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/#Vectorized-constraints). 

Similarly, we can define an array of constraints with varied indexes by including
an additional argument before the constraint expression. For example,
let's define ``3z_i - 14 = 0, \ \forall i \in \{1,2\}``:
```jldoctest constrs
julia> @constraint(model, c2[i = 1:2], 3z[i] - 14 == 0)
2-element Array{InfOptConstraintRef,1}:
 c2[1] : 3 z[1] = 14.0
 c2[2] : 3 z[2] = 14.0
```
Thus, we added two constraints to `model` and stored a vector of the corresponding
constraint references to the `Julia` variable `c2`. To learn more about building 
containers of constraints please see 
[`JuMP`'s constraint container documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/#Constraint-containers).

### Multi-Dimensional Constraints
Building upon `JuMP` we support a variety of multi-dimensional constraint types. 
For example, we can define the vector constraint:
```jldoctest constrs
julia> A = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4

julia> b = [5, 6]
2-element Array{Int64,1}:
 5
 6

julia> @constraint(model, A * z - b in MOI.Nonnegatives(2))
[z[1] + 2 z[2] - 5, 3 z[1] + 4 z[2] - 6] ∈ MathOptInterface.Nonnegatives(2)
```
See [`JuMP`'s constraint documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/) 
for a thorough tutorial on the accepted syntax and constraint types.

### Restricted Constraints
Restricted constraints entail an infinite domain (determined by the infinite 
parameters they explicitly/implicitly depend on) that is restricted to a certain 
sub-domain. Such constraints are common for enforcing initial/boundary conditions 
and for enforcing path constraints over a certain sub-domain.

!!! warning 
    Previous versions of `InfiniteOpt` referred to restricted constraints as 
    "bounded constraints" and used `@BDconstraint` to define them. This has been 
    deprecated in favor of the more intuitive domain restricted nomenclature.
    ```julia
    # Old syntax
    @BDconstraint(model, name_expr(restricts...), constr_expr)

    # New syntax
    @constraint(model, name_expr, constr_expr, DomainRestrictions(restricts...))
    ```

These types of constraints are defined adding [`DomainRestrictions`](@ref). For 
example, let's add the initial condition ``y_b(0) = 0``:
```jldoctest constrs
julia> @constraint(model, initial, yb == 0, DomainRestrictions(t => 0))
initial : yb(t) = 0.0, ∀ t = 0
```
Thus, we have added a constraint to `model` defined over the sub-domain ``t = 0``
in accordance with the initial condition.

More complex sub-domains can be specified by simply adding more restrictions. To 
illustrate this, let's define the constraint 
``2y_b^2(t, x) + z_1 \geq 3, \ \forall t = 0, \ x \in [-1, 1]^2``:
```jldoctest constrs
julia> @constraint(model, 2ya^2 + z[1] >= 3, DomainRestrictions(t => 0, x => [-1, 1]))
2 ya(t, x)² + z[1] ≥ 3.0, ∀ t = 0, x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
```

Now we have added constraints to our model and it is ready to be solved!

## Data Structure
Here we detail the data structures used to store constraints in `InfiniteOpt`.
In general, constraints in `JuMP` are of the form: `function in set` where
`function` corresponds to a `JuMP` expression and `set` corresponds to a `MOI`
set. This leads to the following data structures:

| Constraint Type | Function Type                       | Set Type                                                                                                                            |
|:---------------:|:-----------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------:|
| Scalar          | `JuMP.AbstractJuMPScalar`           | [`MOI.AbstractScalarSet`](https://jump.dev/MathOptInterface.jl/v0.9.22/reference/standard_form/#MathOptInterface.AbstractScalarSet) |
| Vector          | `Vector{<:JuMP.AbstractJuMPScalar}` | [`MOI.AbstractVectorSet`](https://jump.dev/MathOptInterface.jl/v0.9.22/reference/standard_form/#MathOptInterface.AbstractVectorSet) |
| Matrix          | `Matrix{<:JuMP.AbstractJuMPScalar}` | `MOI.AbstractVectorSet` [via vectorization](https://jump.dev/MathOptInterface.jl/v0.9.22/reference/standard_form/#Matrix-sets)      |

The above combos are then stored in 
[`JuMP.ScalarConstraint`](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.ScalarConstraint)s 
and [`JuMP.VectorConstraint](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.VectorConstraint)s. 

Restricted constraints are built upon this data structure where the underlying 
constraint is created in the same manner. Then the specified 
[`DomainRestrictions`](@ref) are added by creating a 
[`DomainRestrictedConstraint`](@ref) which stores the `JuMP.AbstractConstraint` 
and the restrictions.

These constraint objects are what store constraints in `InfiniteModel`s. And
these are pointed to by [`InfOptConstraintRef`](@ref)s.

## Definition
In this section, we describe the ins and outs of defining constraints. Note that
this process is analogous to the manner in which variables are defined and added
to the model.

### Manual Definition
Defining a constraint principally involves the following steps:
- Define the constraint information (i.e., function, set, and domain restrictions)
- Construct a concrete subtype of `JuMP.AbstractConstraint` to store the 
  constraint information
- Add the `AbstractConstraint` object to an `InfiniteModel` and assign a name
- Create an [`InfOptConstraintRef`](@ref) that points to the constraint object 
  stored in the model.

The constraint objects are specified via `JuMP.build_constraint` which requires 
that the user provides a function, set, and optionally include domain 
restrictions. For example, let's build a scalar constraint 
``3y_a(t, x) - y_b^2(t) \leq 0, \ \forall t \in [0, 10], x \in [-2, 2]^2`` over 
its full infinite domain (i.e., have no `DomainRestrictions`):
```jldoctest constrs
julia> constr = build_constraint(error, 3ya - yb^2, MOI.LessThan(0.0));
```

Now the built constraint object can be added to the infinite model via
[`add_constraint`](@ref). Let's do so with our example and assign it the name of 
`c3` (note that adding a name is optional):
```jldoctest constrs
julia> cref = add_constraint(model, constr, "c3")
c3 : -yb(t)² + 3 ya(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ∈ [-2, 2], x[2] ∈ [-2, 2]
```

Thus, we have made our constraint and added it `model` and now have a constraint
reference `cref` that we can use to access it.

The [`@constraint`](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.@constraint) 
and [`@SDconstraint`](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.@SDconstraint) 
macros automate the above steps.

### Macro Definition
As mentioned above in the Basic Usage section, the 
[`@constraint`](https://jump.dev/JuMP.jl/v0.21.8/reference/constraints/#JuMP.@constraint) 
macro should be used to define constraints with the syntax: 
`@constraint(model::InfiniteModel, [container/name_expr], constr_expr, [rs::DomainRestrictions])`.

The second argument is optional and is used to assign a name and/or define
indexing variables to be used in the constraint expr. When a name is provided it
is registered and cannot be used again for another constraint or variable name.
The indexing expression can be used to produce an array of constraints as shown
below (notice this is equivalent to looping over individual `@constraint` calls):
```jldoctest constrs
julia> crefs = @constraint(model, [i = 1:2], 2z[i] - yb == 0)
2-element Array{InfOptConstraintRef,1}:
 2 z[1] - yb(t) = 0.0, ∀ t ∈ [0, 10]
 2 z[2] - yb(t) = 0.0, ∀ t ∈ [0, 10]

julia> crefs = Vector{InfOptConstraintRef}(undef, 2);

julia> for i = 1:2
           crefs[i] = @constraint(model, 2z[i] - yb == 0)
       end

julia> crefs
2-element Array{InfOptConstraintRef,1}:
 2 z[1] - yb(t) = 0.0, ∀ t ∈ [0, 10]
 2 z[2] - yb(t) = 0.0, ∀ t ∈ [0, 10]
```
Please refer to 
[`JuMP`'s constraint container documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/#Constraint-containers) 
for a thorough tutorial on creating containers of constraints.

Any constraint type supported by `JuMP` can be specified in the `constr_expr` 
argument. This includes a wealth of constraint types including:
- Variable constraints
- Scalar constraints
- Vector constraints
- Conic constraints 
- Indicator constraints
- Semi-definite constraints
For example, we could define the following semi-definite constraint using 
`@SDconstraint`:
```jldoctest constrs
julia> @SDconstraint(model, [yb 2yb; 3yb 4yb] >= ones(2, 2))
[yb(t) - 1    2 yb(t) - 1;
 3 yb(t) - 1  4 yb(t) - 1] ∈ PSDCone(), ∀ t ∈ [0, 10]
```
See [`JuMP`'s constraint documentation](https://jump.dev/JuMP.jl/v0.21.8/manual/constraints/) 
for a thorough tutorial on the accepted syntax and constraint types.

Finally, restrictions on the inherent infinite domain of a constraint can be 
specified via [`DomainRestrictions`](@ref) with the `rs` argument. The accepted 
syntax is `DomainRestrictions(restricts...)` where each argument of `restricts` 
can be any of the following forms:
- `pref => value`
- `pref => [lb, ub]`
- `pref => IntervalDomain(lb, ub)`
- `prefs => value`
- `prefs => [lb, ub]`
- `prefs => IntervalDomain(lb, ub)`.
Note that `pref` and `prefs` must correspond to infinite parameters.

For example, we can define the constraint ``y_a^2(t, x) + z_i \leq 1`` and 
restrict the infinite domain of ``x_i`` to be ``[0, 1]``:
 ```jldoctest constrs
julia> @constraint(model, [i = 1:2], ya^2 + z[i] <= 1, DomainRestrictions(x[i] => [0, 1]))
2-element Array{InfOptConstraintRef,1}:
 ya(t, x)² + z[1] ≤ 1.0, ∀ t ∈ [0, 10], x[1] ∈ [0, 1], x[2] ∈ [-2, 2]
 ya(t, x)² + z[2] ≤ 1.0, ∀ t ∈ [0, 10], x[1] ∈ [-2, 2], x[2] ∈ [0, 1]
```

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
and [`constraint_object`](@ref), respectively. These methods all constitute 
extensions of `JuMP` methods and follow exactly the same behavior. Let's try 
them out with the following example:
```jldoctest constrs
julia> is_valid(model, c1) # check if contained in model
true

julia> name(c1) # get the name
"c1"

julia> m = owner_model(c1); # get the model it is added to

julia> index(c1) # get the constraint's index
InfOptConstraintIndex(1)

julia> constr = constraint_object(c1); # get the raw constraint datatype
```

Also, [`constraint_by_name`](@ref) can be used to retrieve a constraint reference 
if only the name is known and its name is unique. For example, let's extract the 
reference for `"c1"`:
```jldoctest constrs
julia> cref = constraint_by_name(model, "c1")
c1 : z[1]² + z[2]² + 2 ya(t, x) ≤ 0.0, ∀ t ∈ [0, 10], x[1] ∈ [-2, 2], x[2] ∈ [-2, 2]
```

### Domain Restrictions
As explained above, restricted constraints serve as a key capability of
`InfiniteOpt`. Information about domain restrictions can be obtained via
[`has_domain_restrictions`](@ref) and [`domain_restrictions`](@ref) which indicate
if a constraint is restricted and what its [`DomainRestrictions`](@ref) are,
respectively. These are exemplified below:
```jldoctest constrs
julia> has_domain_restrictions(c1) # check if constraint is bounded
false

julia> has_domain_restrictions(initial)
true

julia> domain_restrictions(initial)
Subdomain restrictions (1): t = 0
```

### General
Constraints can be defined in a number of ways symbolically that differ from
how it is actually stored in the model. This principally occurs since like terms
and constants are combined together where possible with the variable terms on the
left hand side and the constant on the right hand side. For instance, the
constraint ``2y_b(t) + 3y_b(t) - 2 \leq 1 + z_1`` would be normalized 
``5y_b(t) - z_1 \leq 3``. In accordance with this behavior, 
[`normalized_rhs`](@ref) and [`normalized_coefficient`](@ref)
can be used to query the normalized right hand side and the coefficient of a
particular variable reference, respectively. Let's employ the above example to
illustrate this:
```jldoctest constrs
julia> @constraint(model, constr, 2yb + 3yb - 2 <= 1 + z[1])
constr : 5 yb(t) - z[1] ≤ 3.0, ∀ t ∈ [0, 10]

julia> normalized_rhs(constr)
3.0

julia> normalized_coefficient(constr, yb)
5.0
```

There also exist a number of methods for querying an infinite model about what
constraints it contains.
[`list_of_constraint_types`](@ref) can be used query what types of constraints 
have been added to a model. This is provided as a list of tuples where the first 
element is the expression type and the second element is the set type (recall 
that constraints are stored in the form `func-in-set`). Thus, for our current 
model we obtain:
```julia-repl
julia> list_of_constraint_types(model)
4-element Array{Tuple{DataType,DataType},1}:
 (GenericQuadExpr{Float64,GeneralVariableRef}, MathOptInterface.LessThan{Float64})
 (GenericQuadExpr{Float64,GeneralVariableRef}, MathOptInterface.GreaterThan{Float64})
 (GenericAffExpr{Float64,GeneralVariableRef}, MathOptInterface.LessThan{Float64})
 (GenericAffExpr{Float64,GeneralVariableRef}, MathOptInterface.EqualTo{Float64})
```
This information is useful when used in combination with the
[`num_constraints`](@ref) and [`all_constraints`](@ref) methods which can take 
the expression type and/or the set type as inputs. Here `num_constraints` 
provides the number of constraints that match a certain type  and `all_constraints` 
returns a list of constraint references matching the criteria provided. These have 
been extended beyond `JuMP` functionality such additional methods have been 
provided for the cases in which one wants to query solely off of set or off 
expression type. Let's illustrate this with `num_constraints`:
```jldoctest constrs
julia> num_constraints(model) # total number of constraints
15

julia> num_constraints(model, GenericQuadExpr{Float64,GeneralVariableRef})
5

julia> num_constraints(model, MOI.LessThan{Float64})
5

julia> num_constraints(model, GenericQuadExpr{Float64,GeneralVariableRef},
                       MOI.LessThan{Float64})
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
julia> delete(model, c1)
```

### General
There also are a number of ways to modify information and characteristics of
constraints. First, [`set_name`](@ref JuMP.set_name(::InfOptConstraintRef, ::String))
can be used to specify a new name for a particular constraint. For instance,
let's update the name of `initial` to `"init_cond"`:
```jldoctest constrs
julia> set_name(initial, "init_cond")

julia> initial
init_cond : yb(t) = 0.0, ∀ t = 0
```

We can also update the normalized right hand side constant value or normalized
left hand side variable coefficient value using
[`set_normalized_rhs`](@ref) and [`set_normalized_coefficient`](@ref), 
respectively. Let's again consider the constraint ``5y_b(t) - z_1 \leq 3`` as an
example. Let's change the constant term to -1 and the `y_b(t)` coefficient to 2.5:
```jldoctest constrs
julia> set_normalized_rhs(constr, -1)

julia> set_normalized_coefficient(constr, yb, 2.5)

julia> constr
constr : 2.5 yb(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 10]
```

!!! note
    In some cases, it may be more convenient to dynamically modify coefficients
    and other values via the use of finite parameters. This provides an avenue
    to update parameters without having to be concerned about the normalized form.
    For more information, see the [Finite Parameters](@ref finite_param_docs) page.

### Domain Restrictions
Domain Restrictions can be added to, modified, or removed from any constraint in
`InfiniteOpt`. Principally, this is accomplished via 
[`add_domain_restrictions`](@ref), [`set_domain_restrictions`](@ref),
and [`delete_domain_restrictions`](@ref).

!!! note 
    Previous versions of `InfiniteOpt` used `@[set/add]_parameter_bounds` which 
    have been deprecated in favor of using [`DomainRestrictions`](@ref) with the 
    the methods described used in this section.

First, domain restrictions can be added to a constraint via 
[`add_domain_restrictions`](@ref). For example, let's add the bound 
``t \in [0, 1]`` to `constr`:
```jldoctest constrs
julia> add_domain_restrictions(constr, DomainRestrictions(t => [0, 1]))

julia> constr
constr : 2.5 yb(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 1]
```

In similar manner, [`set_domain_restrictions`](@ref) can be employed to specify
what restrictions a constraint has (overwriting any existing ones if forced). It 
follows the same syntax, so let's use it to change the bounds on `t` to ``t = 0``:
```jldoctest constrs
julia> set_domain_restrictions(constr, DomainRestrictions(t => 0), force = true)

julia> constr
constr : 2.5 yb(t) - z[1] ≤ -1.0, ∀ t = 0
```

Finally, constraint restrictions can be deleted via
[`delete_domain_restrictions`](@ref). Now let's delete the domain restrictions 
associated with our example:
```jldoctest constrs
julia> delete_domain_restrictions(constr)

julia> constr
constr : 2.5 yb(t) - z[1] ≤ -1.0, ∀ t ∈ [0, 10]
```

## Datatypes
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
DomainRestrictions
DomainRestrictedConstraint
ConstraintData
InfOptConstraintIndex
InfOptConstraintRef
```

## Methods/Macros
```@index
Pages   = ["constraint.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:function]
```
```@docs
JuMP.build_constraint(::Function, ::Any, ::Any, ::DomainRestrictions)
JuMP.add_constraint(::InfiniteModel, ::JuMP.AbstractConstraint, ::String)
JuMP.owner_model(::InfOptConstraintRef)
JuMP.index(::InfOptConstraintRef)
JuMP.constraint_object(::InfOptConstraintRef)
JuMP.name(::InfOptConstraintRef)
JuMP.set_name(::InfOptConstraintRef, ::String)
JuMP.is_valid(::InfiniteModel, ::InfOptConstraintRef)
JuMP.delete(::InfiniteModel, ::InfOptConstraintRef)
parameter_refs(::InfOptConstraintRef)
has_domain_restrictions
domain_restrictions
set_domain_restrictions
add_domain_restrictions
delete_domain_restrictions
JuMP.set_normalized_rhs(::InfOptConstraintRef, ::Real)
JuMP.normalized_rhs(::InfOptConstraintRef)
JuMP.add_to_function_constant(::InfOptConstraintRef, ::Real)
JuMP.set_normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef, ::Real)
JuMP.normalized_coefficient(::InfOptConstraintRef, ::GeneralVariableRef)
JuMP.constraint_by_name(::InfiniteModel, ::String)
JuMP.list_of_constraint_types(::InfiniteModel)
JuMP.num_constraints(::InfiniteModel, ::Any, ::Any)
JuMP.all_constraints(::InfiniteModel, ::Any, ::Any)
```
