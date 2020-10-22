# [Derivative Operators](@id deriv_page)
A guide and manual for the definition and use of derivatives in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview
Derivative operators commonly arise in many infinite-dimensional problems, 
particularly in space-time optimization. `InfiniteOpt.jl` provides a simple yet 
powerful interface to model these objects for derivatives of any order, including 
partial derivatives. Derivatives can be used in defining measures and constraints. 

## Basic Usage
Derivative operators can defined a few different ways in `InfiniteOpt`. To motivate 
these, let's first define an `InfiniteModel` along with some parameters and variables:
```jldoctest deriv_basic
julia> using InfiniteOpt, JuMP, Distributions;

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10], 
                           derivative_method = OrthogonalCollocation(3));

julia> @infinite_parameter(model, ξ in Uniform(-1, 1));

julia> @infinite_variable(model, y(t, ξ));

julia> @infinite_variable(model, q(t));
```
Notice that we used the `derivative_method` keyword argument to specify which 
numerical method will be used to evaluate any derivatives that depend on that 
infinite parameter `t`. In this case we, specified to use orthogonal collocation 
over finite elements using 3 nodes. We'll come back to this just a little further 
below to more fully describe the various methods we can use. 

First, let's discuss how to define derivatives in `InfiniteOpt.jl`. Principally, 
this is accomplished via [`@deriv`](@ref) which will operate on a particular 
`InfiniteOpt` expression (containing parameters, variables, and/or measures) with 
respect to infinite parameters specified with their associated orders. Behind the 
scenes all the appropriate calculus will be applied, creating derivative variables 
as needed. For example, we can define the following:
```jldoctest deriv_basic 
julia> d1 = @deriv(y, t)
∂/∂t[y(t, ξ)]

julia> d2 = @deriv(y, t, ξ)
∂/∂ξ[∂/∂t[y(t, ξ)]]

julia> d3 = @deriv(q, t^2)
∂/∂t[∂/∂t[q(t)]]

julia> d_expr = @deriv(y * q - 2t, t)
∂/∂t[y(t, ξ)]*q(t) + ∂/∂t[q(t)]*y(t, ξ) - 2
```
Thus, we can define derivatives in a variety of forms according to the problem at 
hand. The last example even shows hwo the product rule is correctly applied. 

Also, notice that the appropriate analytic calculus is applied to infinite 
parameters. For example, we could also compute:
```jldoctest deriv_basic 
julia> @deriv(3t^2 - 2t, t)
6 t - 2
```

Conveniently, `@deriv` can be called within any measure and constraint. However, 
in certain cases we may need to define an initial guess (initial guess trajectory). 
This can be accomplished in 2 ways:
 - Call [`set_start_value_function`](@ref set_start_value_function(::DerivativeRef,::Union{Real, Function})) 
   using the individual derivative (e.g., `d1` above)
 - Define the derivative using [`@derivative_variable`](@ref) and use the `start` keyword argument.
In either case, a single value can be given or a start value function that will
generate a value in accordance with the support values (i.e., following the same 
syntax as infinite variables). For example, we can specify the starting value of 
`d1` to `0` via the following:
```jldoctest deriv_basic 
julia> set_start_value_function(d1, 0)
```

Now let's return to our discussion on derivative evaluation methods. These are the 
methods that can/will be invoked to transcript the derivatives when solving the 
model. The methods native to `InfiniteOpt` are described in the table below:

| Method                         | Type               | Needed Boundary Conditions| Creates Supports |
|:------------------------------:|:------------------:|:-------------------------:|:----------------:|
|[`FiniteDifference`](@ref)      | [`Forward`](@ref)  | Final & optional Initial  | No               |
|[`FiniteDifference`](@ref)      | [`Central`](@ref)  | Initial & Final           | No               |
|[`FiniteDifference`](@ref)      | [`Backward`](@ref) | Initial & optional Final  | No               |
|[`OrthogonalCollocation`](@ref) | [`Lobatto`](@ref)  | Initial                   | Yes              |

Here the default method is backward finite difference. These are enforced on an 
infinite parameter basis (i.e., the parameter the differential operator is taken 
with respect to). Thus, in the above examples any derivatives taken with respect to 
`t` will use orthogonal collocation on finite elements since that is what we 
specified as our derivative method. More information is provided in the 
[Derivative Methods](@ref) Section below. However, we note here that 
[`set_derivative_method`](@ref) can be invoked anytime after parameter definition 
to specify/modify the derivative method used. More conveniently, we can call 
[`set_all_derivative_methods`](@ref):
```jldoctest deriv_basic 
julia> set_all_derivative_methods(model, FiniteDifference(Forward))

```

!!! warn 
    `InfiniteOpt` does not ensure proper boundary conditions are provided by the 
    user. Thus, it is imperative that the user ensure these are provided appropriately 
    with the derivative evaluation method that is used. We recommend specifying 
    such conditions via [`@BDconstraint`](@ref).

## Advanced Definition 
This section will detail the inner-workings and more advanced details behind 
defining derivatives in `InfiniteOpt`.

### Manual Definition 
The workflow for derivative definition mirrors that of variable definition as 
summarized in the following steps:
 1. Define the variable information via a `JuMP.VariableInfo`.
 2. Build the derivative using [`build_derivative`](@ref).
 3. Add the derivative to the model via [`add_derivative`](@ref).

To exemplify this process, let's first define appropriate variable information:
```jldoctest deriv_basic 
julia> info = VariableInfo(true, 0., true, 42., false, 0., false, 0., false, false);
```
More detailed information on `JuMP.VariableInfo` is provided in the 
[Variable Definition Methodology](@ref) section. 

!!! warn 
    Instances of `JuMP.VariableInfo` used to define derivatives should have 
    `info.binary = false` and `info.integer = false`, since most derivative 
    evaluation methods require that derivatives be continuous.

Now that we have our variable information we can make a derivative using 
[`build_derivative`](@ref):
```jldoctest deriv_basic 
julia> d = build_derivative(error, info, y, ξ);

julia> d isa Derivative
true
```
Here the argument variable can be an infinite variable, reduced variable, 
derivative, or measure that depends on the infinite parameter provided. This will 
error to the contrary or if such a derivative has already been to the model 
associated with the infinite parameter. 

Now we can add the derivative to the model via [`add_derivative`](@ref) which 
will add the [`Derivative`](@ref) object and return `GeneralVariableRef` pointing 
to it that we can use in `InfiniteOpt` expressions:
```jldoctest deriv_basic 
julia> dref = add_derivative(model, d)
∂/∂ξ[y(t, ξ)]
```
This will also create any appropriate information based constraints (e.g., lower 
bounds).

Finally, we note that higher order derivatives are made by simply nesting this 
process.

### Macro Definition 
There are two macros we provide for defining derivatives: 
[`@derivative_variable`](@ref) and [`@deriv`](@ref). 

First, [`@derivative_variable`](@ref)
simply automates the process described above in a manner inspired the by the syntax 
of the variable macros. As such it will support all the same keywords and 
constraint syntax used with the variable macros. For example, we can define the 
derivative ``\frac{\partial^2 y(t, \xi)}{\partial t^2}`` using `d1` (defined in the 
a Basic Usage section) enforcing a lower bound of 1 with an initial guess of 0 and 
assign it to an alias `GeneralVariableRef` called `dydt2`:
```jldoctest deriv_basic 
julia> @derivative_variable(model, d(d1)/d(t), dydt2 >= 1, start = 0)
dydt2(t, ξ)
```
This will also support anonymous definition and multi-dimensional definition,  
please refer to [`@derivative_variable`](@ref) in the manual for the full details.

Second, for more convenient definition we use [`@deriv`](@ref) as shown in the 
Basic Usage section above. Unlike `@derivative_variable` this can handle any 
`InfiniteOpt` expression as the argument input and will automatically take care of 
any redundant derivative creation by using the existing derivatives as appropriate. 
It also can build derivatives that depend on multiple infinite parameters and/or 
are taken to higher orders. This is accomplished via recursive derivative 
definition, handling the nesting as appropriate. For example, we can "define" 
``\frac{\partial^2 y(t, \xi)}{\partial t^2}`` again:
```jldoctest deriv_basic 
julia> @deriv(d1, t)
dydt2(t, ξ)

julia> @deriv(y, t^2)
dydt2(t, ξ)
```
Notice that no error is thrown (which would have occurred if we called 
`@derivative_variable` again) and that the derivative references all point to the 
same derivative object we defined up above with its alias name `dydt2`. This macro 
can also tackle complex expressions using the appropriate calculus such as:
```jldoctest deriv_basic 
julia> @deriv(integral(y, ξ) * q, t)
∂/∂t[integral{ξ ∈ [-1, 1]}[y(t, ξ)]]*q(t) + ∂/∂t[q(t)]*integral{ξ ∈ [-1, 1]}[y(t, ξ)]
```
Thus, demonstrating the convenience of using `@deriv`.

With all this in mind, we recommend using `@deriv` as the defacto method, but then 
using `@derivative_variable` as a convenient way to specify information constraints 
and an initial guess value/trajectory. 

## Derivative Evaluation
In this section, we detail how derivatives are evaluated in `InfiniteOpt` to then 
be used in reformulating the model for solution. 

### Theory 
To motivate the principles behind numerical derivative evaluation/transcription, 
let's first consider the initial value problem:
```math 
\frac{d y(t)}{dt} = f(t, y(t)), \ \ \ y(t_0) = y_0
```
With a finite support set ``\{t_0, t_1, \dots, t_k\}`` we can numerically 
approximate the value of ``\frac{d y(t_n)}{dt}`` at each time point ``t_n`` via 
the Euler method (i.e., forward finite difference). We thus obtain a system of 
equations:
```math
\begin{aligned}
&&& y(t_{n+1}) = y(t_n) + (t_{n+1} - t_n) \frac{d y(t_n)}{dt}, && \forall n = 0, \dots, k-1\\
&&& \frac{d y(t_n)}{dt} = f(t_n, y(t_n)), && \forall n = 0, \dots, k \\ 
&&& y(t_0) = y_0
\end{aligned}
```
Thus, we obtain 3 sets of equations: 
 1. constraint transcriptions
 2. auxiliary derivative equations 
 3. boundary conditions. 

In the case above, we could reduce the number of equations by substituting out the 
point derivatives in the constraint transcriptions since we have explicit 
relationships in the auxiliary equations. However, this is not possible in general, 
such as when we encounter more complex partial differential equations. 

Thus, in `InfiniteOpt` derivatives are treated as variables which can be contained 
implicitly in constraints and/or measures. This allows us to support implicit 
dependencies and higher order derivatives. This means that when the model is 
reformulated, its constraints and measures can be reformulated as normal 
(treating any derivative dependencies as variables). We then can apply the 
appropriate derivative evaluation technique to derive the necessary set of 
auxiliary derivative equations to properly characterize the derivative variables. 
This can be formalized as:
```math 
\begin{aligned}
&&& f_j(y(\lambda), Dy(\lambda)) \leq 0, && \forall j \in J, \lambda \in \Lambda \\
&&& h_i(y(\lambda), Dy(\lambda)) == 0, && \forall i \in I, \lambda \in \Lambda \\
&&& g_k(y(\hat{\lambda}), Dy(\hat{\lambda})) == 0, && \forall k \in K, \hat{\lambda} \in \hat{\Lambda}
\end{aligned}
```
where ``y(\lambda)`` and ``Dy(\lambda)`` denote all the variables and derivatives 
in the problem and ``\lambda`` the denote all the problem's infinite parameters. 
With this let the constraints ``f_j`` denote the problem constraints which can 
contain any variables, parameters, derivatives, and/or measures associated with 
the problem. The constraints ``h_i`` denote the auxiliary derivative equations 
formed by the appropriate numerical method to implicitly define the behavior of 
the derivative variables present in ``f_j``. Finally, the necessary boundary 
conditions are provided in the constraints ``g_k``.

Note that this general paradigm captures a wide breadth of problems and 
derivative evaluation techniques. Higher order derivatives are dealt with naturally 
since such techniques can be applied to nested derivative operators recursively. 
For example, consider the second-order partial derivative:
```math
\frac{\partial^2 y(t, \xi)}{\partial t^2} = \frac{\partial}{\partial t}\left(\frac{\partial y(t, \xi)}{\partial t}\right)
```
The 2 forms are equivalent thus when we apply the Euler method we obtain the 
following auxiliary equations:
```math
\begin{aligned}
&&& y(t_{n+1}, \xi) = y(t_n, \xi) + (t_{n+1} - t_n) \frac{\partial y(t_n, \xi)}{\partial t}, && \forall \xi \in \mathcal{D}_\xi, n = 0, \dots, k-1\\
&&& \frac{\partial y(t_{n+1}, \xi)}{\partial t} = \frac{\partial y(t_n, \xi)}{\partial t} + (t_{n+1} - t_n) \frac{\partial^2 y(t_n, \xi)}{\partial t^2}, && \forall \xi \in \mathcal{D}_\xi, n = 0, \dots, k-1\\
\end{aligned}
```

In the section below we detail the derivative evaluation methods that `InfiniteOpt` 
natively implements.

### Derivative Methods
As discussed briefly above in the Basic Usage section, we natively employ 4 
derivative methods in `InfiniteOpt` which are summarized:

| Method                         | Type               | Needed Boundary Conditions| Creates Supports |
|:------------------------------:|:------------------:|:-------------------------:|:----------------:|
|[`FiniteDifference`](@ref)      | [`Forward`](@ref)  | Final & optional Initial  | No               |
|[`FiniteDifference`](@ref)      | [`Central`](@ref)  | Initial & Final           | No               |
|[`FiniteDifference`](@ref)      | [`Backward`](@ref) | Initial & optional Final  | No               |
|[`OrthogonalCollocation`](@ref) | [`Lobatto`](@ref)  | Initial                   | Yes              |

These methods are defined in association with individual infinite parameters and 
will be applied to any derivatives that are taken with respect to that parameter. 
These methods are specified via the `derivative_method` keyword argument in the 
[`@infinite_parameter`](@ref) macro and can also be defined by invoking 
[`set_derivative_method`](@ref) or [`set_all_derivative_methods`](@ref):
```jldoctest deriv_basic
julia> set_derivative_method(t, FiniteDifference(Forward))

```
In this example, we set `t`'s derivative evaluation method to use forward finite 
difference. This will also reset any changes that were made with the old method 
(e.g., removing old collocation points). Now let's describe the ins and outs of 
these methods.

The first class of methods pertain to finite difference techniques. The syntax 
for specifying these techniques is described in [`FiniteDifference`](@ref) and 
exemplified here:
```jldoctest deriv_basic
julia> FiniteDifference(Forward, true)
FiniteDifference(Forward, true)
```
where the first argument indicates the type of finite difference we wish to employ 
and the second argument indicates if this method should be enforced on boundary 
points. By default, we have `FiniteDifference(Backward, true)` which is the default 
for all infinite parameters. 

Forward finite difference (i.e., explicit Euler) is exemplified by approximating first 
order derivative ``\frac{d y(t)}{dt}`` via 
```math
y(t_{n+1}) = y(t_n) + (t_{n+1} - t_{n})\frac{d y(t_n)}{dt}, \ \forall n = 0, 1, \dots, k-1
```
Note that in this case, the boundary relation corresponds to ``n = 0`` and would 
be included if we set `FiniteDifference(Forward, true)` or would excluded if we 
let the second argument be `false`. We recommend, selecting `false` when an initial 
condition is provided. Also, note that a terminal condition should be provided 
when using this method since an auxiliary equation for the derivative at the 
terminal point cannot be made. Thus, if a terminal condition is not given terminal 
point derivative will be a free variable.

Central finite difference is exemplified by approximating the first order derivative 
``\frac{d y(t)}{dt}`` via
```math 
y(t_{n+1}) = y(t_{n-1}) + (t_{n+1} - t_{n-1})\frac{d y(t_n)}{dt}, \ \forall n = 1, 2, \dots, k-1
```
Note that this form cannot be invoked at ``n = 0`` or ``n = k`` and cannot 
an equation at either boundary. With this in mind the syntax is `FiniteDifference(Central)` 
where the second argument is omitted since it doesn't apply to this scheme. As a 
result both initial and terminal conditions should be specified otherwise the 
derivatives at those points will be free variables.

Backward finite difference (i.e., implicit euler) is our last (and default) 
finite difference method and is exemplified by approximating the first order 
derivative ``\frac{d y(t)}{dt}`` via
```math
y(t_{n}) = y(t_{n-1}) + (t_{n} - t_{n-1})\frac{d y(t_{n})}{dt}, \ \forall n = 1, 2, \dots, k
```
Here the boundary case corresponds to ``n = k`` and would be included if we set 
`FiniteDifference(Backward, true)` (the default) or excluded if we set the second 
argument to `false`. We recommend, selecting `false` when a terminal condition is 
provided. Also, note that an initial condition should always be given otherwise 
the derivative at the first point will be free.

Finally, we employ orthogonal collocation on finite elements via the 
[`OrthogonalCollocation`](@ref) object (please refer to it in the manual for 
complete syntax details). In general terms, this technique fits an ``m`` degree 
polynomial to each finite element (i.e., sequential support pair) and this fit is 
done via ``m+1`` collocation nodes (supports) which include the finite element 
supports along with ``m-1`` additional internal collocation nodes chosen at 
orthogonal points to the polynomial. The typical syntax for specifying this method 
is `OrthogonalCollocation(num_nodes)` where `num_nodes` indicates the number 
collocation nodes to be used for each finite element. For example, we can specify 
to use 3 collocation nodes (i.e., 1 internal node per finite element) corresponding 
to a 2nd degree polynomial via
```jldoctest deriv_basic
julia> OrthogonalCollocation(3)
OrthogonalCollocation(1, Lobatto)
```
Notice that the 2nd attribute is `Lobatto` which indicates that we are using 
collocation nodes selected via Lobatto quadrature. This is currently the only 
supported technique employed by `OrthogonalCollocation` although more may be added 
in future versions. Please note that an initial condition must be provided otherwise 
the corresponding derivative will be free variable. For more information on 
orthogonal collocation over finite elements, this 
[page](http://apmonitor.com/do/index.php/Main/OrthogonalCollocation) provides a 
good reference.

Other methods can be employed via user-defined extensions. Please visit our 
Extensions page for more information.

### User-Invoked Evaluation
Typically, derivative evaluation is handled when the model is reformulated in such 
a way that the `InfiniteModel` is unmodified such that modifications and repeated 
solutions can be done efficiently and seamlessly. This is also the recommended 
workflow. However, we do provide user accessible derivative evaluation methods 
that generate the auxiliary derivative equations and add them to the `InfiniteModel`.
This can be useful for visualizing how these techniques work and can be helpful for 
user-defined reformulation extensions (i.e., optimizer model extensions).

We can build these relations for a particular derivative via [`evaluate`](@ref). 
For example, let's build evaluation equations for `d1`:
```jldoctest deriv_basic
julia> d1 
∂/∂t[y(t, ξ)]

julia> fill_in_supports!(t, num_supports = 3) # add supports first

julia> evaluate(d1)

julia> derivative_constraints(d1)
2-element Array{InfOptConstraintRef,1}:
 5 ∂/∂t[y(t, ξ)](5, ξ) - y(10, ξ) + y(5, ξ) = 0.0, ∀ ξ ~ Uniform
 5 ∂/∂t[y(t, ξ)](0, ξ) - y(5, ξ) + y(0, ξ) = 0.0, ∀ ξ ~ Uniform
```
Note that we made sure `t` had supports first over which we could carry out the 
evaluation, otherwise an error would have been thrown. Moreover, once the 
evaluation was completed we were able to access the auxiliary equations via 
[`derivative_constraints`](@ref). 

We can also, add the necessary auxiliary equations for all the derivatives in the 
model if we call [`evaluate_all_derivatives!`](@ref):
```jldoctest deriv_basic
julia> fill_in_supports!(ξ, num_supports = 4) # add supports first

julia> evaluate_all_derivatives!(model)

julia> derivative_constraints(dydt2)
2-element Array{InfOptConstraintRef,1}:
 5 dydt2(5, ξ) - ∂/∂t[y(t, ξ)](10, ξ) + ∂/∂t[y(t, ξ)](5, ξ) = 0.0, ∀ ξ ~ Uniform
 5 dydt2(0, ξ) - ∂/∂t[y(t, ξ)](5, ξ) + ∂/∂t[y(t, ξ)](0, ξ) = 0.0, ∀ ξ ~ Uniform
```

Finally, we note that once derivative constraints have been added to the 
`InfiniteModel` any changes to the respective infinite parameter sets, supports, 
or derivative method will necessitate the deletion of these auxiliary constraints 
and a warning will be thrown to indicate such:
```jldoctest deriv_basic
julia> derivative_constraints(d1)
2-element Array{InfOptConstraintRef,1}:
 5 ∂/∂t[y(t, ξ)](5, ξ) - y(10, ξ) + y(5, ξ) = 0.0, ∀ ξ ~ Uniform
 5 ∂/∂t[y(t, ξ)](0, ξ) - y(5, ξ) + y(0, ξ) = 0.0, ∀ ξ ~ Uniform

julia> add_supports(t, 0.2)
┌ Warning: Support/method changes will invalidate existing derivative evaluation constraints that have been added to the InfiniteModel. Thus, these are being deleted.
└ @ InfiniteOpt ~/build/pulsipher/InfiniteOpt.jl/src/scalar_parameters.jl:656

julia> has_derivative_constraints(d1)
false
```

## Query Methods
Here we describe the various query techniques that we can employ on derivatives 
in `InfiniteOpt`.

### Basic Queries 
First, let's overview the basic object inquiries: [`derivative_argument`](@ref), 
[`operator_parameter`](@ref), [`derivative_method`](@ref), and
[`name`](@ref JuMP.name(::DecisionVariableRef)):
```jldoctest deriv_basic
julia> derivative_argument(dydt2) # get the variable the derivative operates on
∂/∂t[y(t, ξ)]

julia> operator_parameter(dydt2) # get the parameter the operator is taken with respect to
t

julia> derivative_method(dydt2) # get the numerical derivative evaluation method
FiniteDifference(Forward, true)

julia> name(dydt2) # get the name if there is one
"dydt2"
```
These all work as exemplified above. We note that `derivative_method` simply 
queries the derivative method associated with the operator parameter.

Derivatives also inherit all the usage methods employed by infinite variables. 
For example:
```jldoctest deriv_basic
julia> is_used(d1)
true

julia> used_by_measure(dydt2)
false

julia> used_by_reduced_variable(d2)
true
```

Also, since derivatives are analogous to infinite variables, they inherit many 
of the same queries including [`parameter_refs`](@ref):
```jldoctest deriv_basic
julia> parameter_refs(d1)
(t, ξ)

julia> parameter_refs(derivative_argument(d1))
(t, ξ)
```
Since derivatives simply inherit their infinite parameter dependencies from the 
argument variable, the above lines are equivalent.

### Variable Information 
Again, since derivatives are essentially a special case of infinite variables, they 
inherit all the same methods for querying variable information. For example, 
consider the following queries:
```jldoctest deriv_basic
julia> has_lower_bound(dydt2)
true

julia> lower_bound(dydt2)
1.0

julia> LowerBoundRef(dydt2)
dydt2(t, ξ) ≥ 1.0, ∀ t ∈ [0, 10], ξ ~ Uniform

julia> has_upper_bound(dydt2)
false 

julia> start_value_function(dydt2)
#139 (generic function with 1 method)
```

### Model Queries 
We can also determine the number of derivatives a model contains and obtain a list 
of them via [`num_derivatives`](@ref) and [`all_derivatives`](@ref), respectively:
```jldoctest deriv_basic
julia> num_derivatives(model)
7

julia> all_derivatives(model)
7-element Array{GeneralVariableRef,1}:
 ∂/∂t[y(t, ξ)]
 ∂/∂ξ[∂/∂t[y(t, ξ)]]
 ∂/∂t[q(t)]
 ∂/∂t[∂/∂t[q(t)]]
 ∂/∂ξ[y(t, ξ)]
 dydt2(t, ξ)
 ∂/∂t[integral{ξ ∈ [-1, 1]}[y(t, ξ)]]
```

## Modification Methods 
In this section, we'll highlight some of the modification methods that can be 
used on derivatives in `InfiniteOpt`.

### Variable Information 
As discussed above, derivatives inherit the same variable methods as infinite 
variables. Thus we can modify/delete bounds and starting values for derivatives 
using the same methods. For example:
```jldoctest deriv_basic
julia> set_lower_bound(dydt2, 0)

julia> lower_bound(dydt2)
0.0

julia> set_upper_bound(dydt2, 2)

julia> upper_bound(dydt2)
2.0

julia> fix(dydt2, 42, force = true)

julia> fix_value(dydt2) 
42.0

julia> set_start_value_function(dydt2, (t, xi) -> t + xi)

julia> unfix(dydt2)

```

### Deletion
Finally, the are 2 deletion methods we can employ apart from deleting variable 
information. First, we can employ [`delete_derivative_constraints`](@ref) to 
delete any derivative evaluation constraints associated with a particular 
derivative:
```jldoctest deriv_basic
julia> delete_derivative_constraints(d2)

julia> has_derivative_constraints(d2)
false
```

Lastly, we can employ `delete` to delete a particular derivative and all its 
dependencies:
```jldoctest deriv_basic
julia> delete(model, d2)

julia> is_valid(model, d2)
false
```

## Datatypes
```@index
Pages   = ["derivative.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
DerivativeIndex
DerivativeRef
Derivative
AbstractDerivativeMethod
GenerativeDerivativeMethod
OrthogonalCollocation
OCTechnique
Lobatto
NonGenerativeDerivativeMethod
FiniteDifference
FDTechnique
Forward
Central
Backward
```

## [Methods/Macros] (@id deriv_methods)
```@index
Pages   = ["derivative.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@deriv
deriv
@derivative_variable
build_derivative
add_derivative
derivative_argument(::DerivativeRef)
operator_parameter(::DerivativeRef)
derivative_method(::DerivativeRef)
raw_parameter_refs(::DerivativeRef)
parameter_refs(::DerivativeRef)
parameter_list(::DerivativeRef)
set_start_value_function(::DerivativeRef,::Union{Real, Function})
reset_start_value_function(::DerivativeRef)
num_derivatives
all_derivatives
set_derivative_method(::IndependentParameterRef, ::AbstractDerivativeMethod)
set_derivative_method(::DependentParameterRef, ::AbstractDerivativeMethod)
set_all_derivative_methods
evaluate(::DerivativeRef)
evaluate_all_derivatives!
has_derivative_constraints(::DerivativeRef)
derivative_constraints(::DerivativeRef)
delete_derivative_constraints(::DerivativeRef)
evaluate_derivative
InfiniteOpt.support_label(::AbstractDerivativeMethod)
InfiniteOpt.generate_derivative_supports
InfiniteOpt.add_derivative_supports(::IndependentParameterRef)
InfiniteOpt.make_reduced_expr
```
