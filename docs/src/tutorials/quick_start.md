```@meta
DocTestFilters = [r"≤|<=", r" == | = ", r" ∈ | in ", r" for all | ∀ ", r"d|∂", 
                  r"E|𝔼", r"integral|∫"]
```

# Quick Start Guide
Below we exemplify and briefly explain the very basics behind defining and solving 
an infinite-dimensional optimization problem in `InfiniteOpt`. Please refer to the 
Guide on our subsequent pages for more complete information. The Basic Usage sections 
on each guide page are good places to start from. Also, the syntax of `InfiniteOpt` 
is inspired by `JuMP` thus we recommend new users that haven't used `JuMP`, first 
consult their tutorials starting 
[here](https://jump.dev/JuMP.jl/v1/tutorials/getting_started/getting_started_with_JuMP/).

## Preliminaries 
### Software Setup
First, we need to make sure everything is installed. This will include:
- installing Julia 
- installing the `InfiniteOpt.jl`, `JuMP.jl`, and `Distributions.jl` packages
- installing wanted optimizers e.g., `Ipopt.jl` and `HiGHS.jl`
See [Installation](@ref) for more information.

### Problem Formulation
Now we need to formulate the problem we want to solve mathematically. For example, 
let's define a simple optimal control model:
```math
\begin{aligned}
	&&\underset{x_i(t, \xi), v_i(t, \xi), y_w(\xi), u_i(t)}{\text{min}} &&& \int_{t \in \mathcal{D}_t} \sum_{i \in I} u_i^2(t) dt \\
	&&\text{s.t.} &&& x_i(0, \xi) = x0_i, && \forall i \in I, \xi \in \mathcal{D}_\xi\\
    &&&&& v_i(0, \xi) = v0_i, && \forall i \in I, \xi \in \mathcal{D}_\xi \\
	&&&&& \frac{\partial x_i(t, \xi)}{\partial t} = v_i(t, \xi), && \forall i \in I, t \in \mathcal{D}_t, \xi \in \mathcal{D}_\xi\\
    &&&&& \xi\frac{\partial v_i(t, \xi)}{\partial t} = u_i(t), && \forall i \in I, t \in \mathcal{D}_t, \xi \in \mathcal{D}_\xi\\
    &&&&& y_{w}(\xi) = \sum_{i \in I}(x_i(t_w, \xi) - p_{iw})^2, && \forall w \in W, \xi \in \mathcal{D}_\xi \\
    &&&&& y_{w}(\xi) \geq 0, && \forall w \in W, \xi \in \mathcal{D}_\xi \\
    &&&&& \mathbb{E}_{\xi}\left[\sum_{w \in W} y_w(\xi) \right] \leq \epsilon \\
    &&&&& \xi \sim \mathcal{N}(\mu, \sigma^2) \\
    &&&&& t \in \mathcal{D}_t
\end{aligned}
```
Notice this model is both dynamic with time ``t`` and random with respect to ``\xi``.

### Parameter Specification
Before moving on we'll need to define the necessary constants and problem 
parameters. Thus, continuing with our example we define the following in our 
Julia session (these could also be put into a script as is shown at the bottom 
of this page):
```jldoctest quick
julia> μ = 1; σ = 0.2; # set the distribution parameters 

julia> x0 = [0, 0]; v0 = [0, 0]; # set the initial conditions

julia> p = [1 4 6 1; 1 3 0 1]; tw = [0, 25, 50, 60]; # set waypoint specifications

julia> I = 1:2; W = 1:4; # set the finite domains
```

## Model Definition 
### Model Initialization 
The first thing we need to do is initialize our `InfiniteModel` and assign an 
appropriate optimizer that will be used to solve its transcripted variant. For our 
little example let's choose to use Ipopt:
```jldoctest quick
julia> using InfiniteOpt, Distributions, Ipopt;

julia> model = InfiniteModel(Ipopt.Optimizer)
An InfiniteOpt Model
Feasibility problem with:
  Finite parameters: 0
  Infinite parameters: 0
  Variables: 0
  Derivatives: 0
  Measures: 0
Transformation backend information:
  Backend type: TranscriptionBackend
  Solver: Ipopt
  Transformation built and up-to-date: false
```
Learn more about `InfiniteModel`s and optimizers on our 
[Infinite Models](@ref infinite_model_docs) page.

Before moving on, let's go ahead make a finite parameter via [`@finite_parameter`](@ref) 
for ``\epsilon`` since this likely a constant we'll want to update repeatedly 
(e.g., to determine a tradeoff curve by varying it): 
```jldoctest quick
julia> @finite_parameter(model, ϵ == 10)
ϵ
```
Learn more about finite parameters on our [Finite Parameters](@ref finite_param_docs) 
page.

### Infinite Parameters 
The next thing we need to do is identify the infinite domains our problem contains 
and define an infinite parameter(s) for each one via [`@infinite_parameter`]. For 
this problem we have the time domain ``t \in \mathcal{D}_t`` and the random domain 
``\xi \in \mathcal{D}_\xi`` where ``\xi \sim \mathcal{N}(\mu, \sigma^2)``:
```jldoctest quick
julia> @infinite_parameter(model, t in [0, 60], num_supports = 61, 
                           derivative_method = OrthogonalCollocation(3))
t

julia> @infinite_parameter(model, ξ ~ Normal(μ, σ^2), num_supports = 10)
ξ
```
Notice we specify the domain/distribution the parameter depends on via `in`.
Here we also specify the number of finite supports we desire for each parameter 
that will ultimately be used to reformulate and solve the problem (i.e., discretize). 
We also specify the derivative evaluation method associated with `t` that will be 
used evaluate the derivatives numerically. See more information about parameters 
on our [Infinite Parameters] (@ref inf_par_docs) page. Also learn more about 
derivative methods on our [Derivative Operators](@ref deriv_docs) page.

### Variables 
Now that we have an `InfiniteModel` and infinite parameters let's define our 
decision variables. First, infinite variables (ones that depend on infinite 
parameters) are defined via 
[`@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
with the addition of the [`Infinite`](@ref) variable type argument to specify the 
infinite parameters it depends on:
 ```jldoctest quick
julia> @variable(model, x[I], Infinite(t, ξ), start = 0)
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{GeneralVariableRef}:
 x[1](t, ξ)
 x[2](t, ξ)

julia> @variable(model, v[I], Infinite(t, ξ), start = 0)
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{GeneralVariableRef}:
 v[1](t, ξ)
 v[2](t, ξ)

julia> @variable(model, u[I], Infinite(t), start = 0)
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{GeneralVariableRef}:
 u[1](t)
 u[2](t)

julia> @variable(model, y[W] >= 0, Infinite(ξ), start = 0)
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, 1:4
And data, a 4-element Vector{GeneralVariableRef}:
 y[1](ξ)
 y[2](ξ)
 y[3](ξ)
 y[4](ξ)
```
Notice that we specify the initial guess for all of them via `start`. We also 
can symbolically define variable conditions like the lower bound on `y`.

That does it for this example, but other problems might also employ the following:
- Finite variables: variables that do not depend on infinite parameters 
  (defined using `@variable`)
- Semi-infinite variables: infinite variables where 1 or more parameters are 
  set a particular point (defined via [Restricted Variables](@ref)).
- Point variables: infinite variables at a particular point (defined via [Restricted Variables](@ref)).

### Objective & Constraints 
Now that the variables and parameters are ready to go, let's define our problem. 
First, we can define the objective using 
[`@objective`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@objective):
 ```jldoctest quick
julia> @objective(model, Min, integral(sum(u[i]^2 for i in I), t))
∫{t ∈ [0, 60]}[u[1](t)² + u[2](t)²]
```
Notice that we also employ [`integral`](@ref) to define the integral. Note that 
objectives must evaluate over all included infinite domains. 

Now let's define the initial conditions using 
[`@constraint`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@constraint) 
in combination with [Restricted Variables](@ref) which will restrict the domain 
of the variables to only be enforced at the initial time:
 ```jldoctest quick
julia> @constraint(model, [i in I], x[i](0, ξ) == x0[i])
1-dimensional DenseAxisArray{InfOptConstraintRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{InfOptConstraintRef}:
 x[1](0, ξ) = 0, ∀ ξ ~ Normal
 x[2](0, ξ) = 0, ∀ ξ ~ Normal

julia> @constraint(model, [i in I], v[i](0, ξ) == v0[i])
1-dimensional DenseAxisArray{InfOptConstraintRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{InfOptConstraintRef}:
 v[1](0, ξ) = 0, ∀ ξ ~ Normal
 v[2](0, ξ) = 0, ∀ ξ ~ Normal
```
Note it is important that we include appropriate boundary conditions when using 
derivatives in our model. For more information please see 
[Derivative Operators](@ref deriv_docs).

Next, we can add our model constraints that have derivatives using 
[`@constraint`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@constraint) 
and [`deriv`](@ref):
 ```jldoctest quick
julia> @constraint(model, c1[i in I], deriv(x[i], t) == v[i])
1-dimensional DenseAxisArray{InfOptConstraintRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{InfOptConstraintRef}:
 c1[1] : ∂/∂t[x[1](t, ξ)] - v[1](t, ξ) = 0, ∀ t ∈ [0, 60], ξ ~ Normal
 c1[2] : ∂/∂t[x[2](t, ξ)] - v[2](t, ξ) = 0, ∀ t ∈ [0, 60], ξ ~ Normal

julia> @constraint(model, c2[i in I], ξ * deriv(v[i], t) == u[i])
1-dimensional DenseAxisArray{InfOptConstraintRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{InfOptConstraintRef}:
 c2[1] : ξ*∂/∂t[v[1](t, ξ)] - u[1](t) = 0, ∀ t ∈ [0, 60], ξ ~ Normal
 c2[2] : ξ*∂/∂t[v[2](t, ξ)] - u[2](t) = 0, ∀ t ∈ [0, 60], ξ ~ Normal
```

Next, we can define our last 2 constraints:
 ```jldoctest quick
julia> @constraint(model, c3[w in W], y[w] == sum((x[i](tw[w], ξ) - p[i, w])^2 for i in I))
1-dimensional DenseAxisArray{InfOptConstraintRef,1,...} with index sets:
    Dimension 1, 1:4
And data, a 4-element Vector{InfOptConstraintRef}:
 c3[1] : -x[1](0, ξ)² - x[2](0, ξ)² + y[1](ξ) + 2 x[1](0, ξ) + 2 x[2](0, ξ) = 2, ∀ ξ ~ Normal
 c3[2] : -x[1](25, ξ)² - x[2](25, ξ)² + y[2](ξ) + 8 x[1](25, ξ) + 6 x[2](25, ξ) = 25, ∀ ξ ~ Normal
 c3[3] : -x[1](50, ξ)² - x[2](50, ξ)² + y[3](ξ) + 12 x[1](50, ξ) = 36, ∀ ξ ~ Normal
 c3[4] : -x[1](60, ξ)² - x[2](60, ξ)² + y[4](ξ) + 2 x[1](60, ξ) + 2 x[2](60, ξ) = 2, ∀ ξ ~ Normal

julia> @constraint(model, c4, expect(sum(y[w] for w in W), ξ) <= ϵ)
c4 : 𝔼{ξ}[y[1](ξ) + y[2](ξ) + y[3](ξ) + y[4](ξ)] - ϵ ≤ 0
```
Notice we are able to invoke an expectation simply by calling [`expect`](@ref).

Finally, to address any unwanted degrees of freedom introduced by internal collocation 
nodes with [`OrthogonalCollocation`](@ref). We should call [`constant_over_collocation`](@ref constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)) 
on any control variables:
```jldoctest quick
julia> constant_over_collocation.(u, t);

``` 

That's it, now we have our problem defined in `InfiniteOpt`!

## Solution & Queries
### Optimize 
Now that our model is defined, let's optimize it via [`optimize!`](@ref):
```jldoctest quick; setup = :(set_attribute(model, "print_level", 0))
julia> optimize!(model)

```
We can check the solution status via 
[`termination_status`](@ref JuMP.termination_status(::InfiniteModel)):
```jldoctest quick
julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4
```
Thus, our model was solved successfully! For more information please see our 
[Optimization](@ref opt_docs) and [Results](@ref result_docs) pages.

### Query the Solution
Finally, we can query a wide variety of information about our solution. Perhaps 
most commonly we'll want to know the objective value and the optimal primal values 
of decision variables. This is accomplished via 
[`objective_value`](@ref JuMP.objective_value(::InfiniteModel)) and 
[`value`](@ref JuMP.value(::GeneralVariableRef)):
```jldoctest quick
julia> opt_obj = objective_value(model);

julia> u_opt = value.(u);
```
Note that u_opt will be multi-dimensional combination with the support values used 
to transcribe `u(t)` along the domain of `t`. We can query those corresponding 
support values via `supports`:
```jldoctest quick
julia> u_ts = supports.(u)
1-dimensional DenseAxisArray{Vector{Tuple},1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{Vector{Tuple}}:
 [(0.0,), (1.0,), (2.0,), (3.0,), (4.0,), (5.0,), (6.0,), (7.0,), (8.0,), (9.0,)  …  (51.0,), (52.0,), (53.0,), (54.0,), (55.0,), (56.0,), (57.0,), (58.0,), (59.0,), (60.0,)]
 [(0.0,), (1.0,), (2.0,), (3.0,), (4.0,), (5.0,), (6.0,), (7.0,), (8.0,), (9.0,)  …  (51.0,), (52.0,), (53.0,), (54.0,), (55.0,), (56.0,), (57.0,), (58.0,), (59.0,), (60.0,)]
```
Please see the [Results](@ref result_docs) page for more information. 

## Summary Script 
The example used in the sections above is summarized in the script below:
```julia
using InfiniteOpt, Distributions, Ipopt

# DEFINE THE PROBLEM CONSTANTS
μ = 1; σ = 0.2
x0 = [0, 0]; v0 = [0, 0]
p = [1 4 6 1; 1 3 0 1]; tw = [0, 25, 50, 60]
I = 1:2; W = 1:4

# INITIALIZE THE MODEL
model = InfiniteModel(Ipopt.Optimizer)

# INITIALIZE THE PARAMETERS
@finite_parameter(model, ϵ == 10)
@infinite_parameter(model, t in [0, 60], num_supports = 61, 
                    derivative_method = OrthogonalCollocation(3))
@infinite_parameter(model, ξ ~ Normal(μ, σ^2), num_supports = 10)

# INITIALIZE THE VARIABLES
@variable(model, x[I], Infinite(t, ξ), start = 0)
@variable(model, v[I], Infinite(t, ξ), start = 0)
@variable(model, u[I], Infinite(t), start = 0)
@variable(model, y[W] >= 0, Infinite(ξ), start = 0)

# SET THE OBJECTIVE
@objective(model, Min, integral(sum(u[i]^2 for i in I), t))

# SET THE INITIAL CONDITIONS
@constraint(model, [i in I], x[i](0, ξ) == x0[i])
@constraint(model, [i in I], v[i](0, ξ) == v0[i])

# SET THE PROBLEM CONSTRAINTS
@constraint(model, c1[i in I], @deriv(x[i], t) == v[i])
@constraint(model, c2[i in I], ξ * @deriv(v[i], t) == u[i])
@constraint(model, c3[w in W], y[w] == sum((x[i](tw[w], ξ) - p[i, w])^2 for i in I))
@constraint(model, c4, expect(sum(y[w] for w in W), ξ) <= ϵ)

# ADJUST DEGREES OF FREEDOM FOR CONTROL VARIABLES
constant_over_collocation.(u, t)

# SOLVE THE MODEL
optimize!(model)

# GET THE RESULTS
termination_status(model)
opt_obj = objective_value(model)
u_opt = value.(u)
u_ts = supports.(u)
```
