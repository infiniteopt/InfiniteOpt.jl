# Results
A guide and manual to querying optimized `InfiniteOpt` models. The Methods
section at the bottom comprises the manual and the above sections form the
guide.

## Overview
So far we have covered defining, transforming, and optimizing `InfiniteModel`s.
Now comes the point to extract information from our optimized model. This is done
following extended versions of `JuMP`s querying functions in combination with
the mapping information stored in the optimizer model. Thus, this page will
walk through the use of these result query functions.

## Basic Usage
Let's revisit the example from the optimization page to get us started:
```jldoctest results
julia> using InfiniteOpt, Ipopt;

julia> model = InfiniteModel(Ipopt.Optimizer);

julia> set_optimizer_attribute(model, "print_level", 0);

julia> @infinite_parameter(model, t in [0, 10], num_supports = 10);

julia> @infinite_variable(model, x(t) >= 0);

julia> @finite_variable(model, z >= 0);

julia> @objective(model, Min, 2z);

julia> @constraint(model, c1, z >= x);

julia> @BDconstraint(model, c2(t == 0), x == 42);

julia> print(model)
Min 2 z
Subject to
 x(t) ≥ 0.0, ∀ t ∈ [0, 10]
 z ≥ 0.0
 c1 : z - x(t) ≥ 0.0, ∀ t ∈ [0, 10]
 c2 : x(t) = 42.0, ∀ t = 0

julia> optimize!(model)

```
Now that the model has been optimized, let's find out what happened. To determine
why the optimizer stopped, we can use
[`termination_status`](@ref JuMP.termination_status(::InfiniteModel)) to report
the corresponding `MathOptInterface` termination code (possible codes are explained
[here](https://jump.dev/JuMP.jl/stable/solutions/#MathOptInterface.TerminationStatusCode).
```jldoctest results
julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4
```
Here we see that our model was locally solved via Ipopt and that is why it
stopped. Furthermore, we can query the primal and dual problem optimalities via
[`primal_status`](@ref JuMP.primal_status(::InfiniteModel)) and
[`dual_status`](@ref JuMP.dual_status(::InfiniteModel)), respectively.
```jldoctest results
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1

julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
The possible statuses are detailed [here](https://jump.dev/JuMP.jl/stable/solutions/#MathOptInterface.ResultStatusCode).
These results are useful in knowing if information can be drawn from the primal
and/or dual and what it means. We can also verify that we indeed have answers
via [`has_values`](@ref JuMP.has_values(::InfiniteModel)) which indicates if our
model has optimized variable values.
```jldoctest results
julia> has_values(model)
true
```
And indeed we do have values.

Now we can query the objective value via
[`objective_value`](@ref JuMP.objective_value(::InfiniteModel)) which reports
the optimal objective value.
```jldoctest results
julia> objective_value(model)
83.99999998250514
```
Great now we can inquire about variable values via
[`value`](@ref JuMP.value(::GeneralVariableRef)). First, let's retrieve the value
of `z`:
```jldoctest results
julia> value(z)
41.99999999125257
```
We get a single value since `z` is a `FiniteVariable` and therefore finite. Now
let's retrieve the "value" of `x(t)` which is infinite with respect to `t`:
```jldoctest results
julia> value(x)
10-element Array{Float64,1}:
 42.0
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
 20.999999995620495
```
Notice here we obtain an array of values since these correspond to the
transcribed finite (discretized) variables used to solve the problem. We obtain
the corresponding support (discretized `t`) values via `supports`:
```jldoctest results
julia> supports(x)
10-element Array{Tuple,1}:
 (0.0,)
 (1.11111111111,)
 (2.22222222222,)
 (3.33333333333,)
 (4.44444444444,)
 (5.55555555556,)
 (6.66666666667,)
 (7.77777777778,)
 (8.88888888889,)
 (10.0,)
```
There is 1-to-1 correspondence between these supports and the values reported
above. Note that these are stored in tuples to facilitate multiple infinite
parameter dependencies.

!!! note
    The values for an array of variables is obtained via the vectorized call
    of `value` following the syntax:
    ```julia
    value.(::AbstractArray{<:GeneralVariableRef})
    ```
    This also holds true for many other methods in `InfiniteOpt`. For example,
    adding the dot also vectorizes `dual` and `set_binary`.

We can also query the dual of a constraint via
[`dual`](@ref JuMP.dual(::InfOptConstraintRef)) if a model has duals available
as indicated by [`has_duals`](@ref JuMP.has_duals(::InfiniteModel)):
```jldoctest results
julia> has_duals(model)
true

julia> dual(c1)
10-element Array{Float64,1}:
 1.9999999988666093
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
```
`c1` is an infinite constraint and thus we obtain the duals of its transcribed
versions. The underlying infinite parameter(s) and support values are queried
via `parameter_refs` and `supports`:
```jldoctest results
julia> parameter_refs(c1)
(t,)

julia> supports(c1)
10-element Array{Tuple,1}:
 (0.0,)
 (1.11111111111,)
 (2.22222222222,)
 (3.33333333333,)
 (4.44444444444,)
 (5.55555555556,)
 (6.66666666667,)
 (7.77777777778,)
 (8.88888888889,)
 (10.0,)
```
These again all have a 1-to-1 correspondence.

!!! note
    In the case that our variables/constraints depend on multiple infinite 
    parameter it is typically convenient to add the keyword statement 
    `ndarray = true` when calling any variable/constraint queries (e.g., `value` 
    and `dual`). This will reformat the output vector into a n-dimensional array 
    whose dimensions correspond to the supports of the infinite parameters. 

## Termination Queries
Termination queries are those that question about how the infinite model was
solved and what its optimized state entails. Programmatically, such queries on
the `InfiniteModel` are simply routed to its optimizer model.

The commonly used queries include
[`termination_status`](@ref JuMP.termination_status(::InfiniteModel)),
[`primal_status`](@ref JuMP.primal_status(::InfiniteModel)),
[`dual_status`](@ref JuMP.dual_status(::InfiniteModel)),
[`objective_value`](@ref JuMP.objective_value(::InfiniteModel)),
[`result_count`](@ref JuMP.result_count(::InfiniteModel))
[`solve_time`](@ref JuMP.solve_time(::InfiniteModel)). The first four are well
exemplified in the Basic Usage section above and are helpful in quickly
understanding the optimality status of a given model following the many possible
statuses reported by `MathOptInterface` which are documented
[here](http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual/#Solving-and-retrieving-the-results-1).
We use `result_count` to determine how many solutions are recorded in the
optimizer.
```jldoctest results
julia> result_count(model)
1
```
This is useful since it informs what results there are which can be specified
via the `result` keyword argument in many methods such as `primal_status`,
`dual_status`, `objective_value`, `value`, `dual`, and more.

We use `solve_time` to determine the time in seconds used by the optimizer until
it terminated its search.
```julia-repl
julia> solve_time(model)
0.004999876022338867
```
Note that this query might not be supported with all solvers.

The above status queries are designed to report information in a consistent
format irrespective of the chosen optimizer. However,
[`raw_status`](@ref JuMP.raw_status(::InfiniteModel)) will provide the optimality
status verbatim as reported by the optimizer. Thus, following our example with
Ipopt we obtain:
```jldoctest results
julia> raw_status(model)
"Solve_Succeeded"
```

Also, we obtain the best objective bound via
[`objective_bound`](@ref JuMP.objective_bound(::InfiniteModel)) which becomes
particularly useful solutions that are suboptimal. However, this method is not
supported by all optimizers and in this case Ipopt is one such optimizer.

Finally, we get the best dual objective value via
[`dual_objective_value`](@ref JuMP.dual_objective_value(::InfiniteModel)) if the
optimizer supplies this information which again Ipopt does not.

## Variable Queries
Information about the optimized variables is gathered consistently in comparison
to typical `JuMP` models. With `InfiniteModel`s this is done by querying the
optimizer model and using its stored variable mappings to return the correct
information. Thus, here the queries are extended to work with the specifics of
the optimizer model to return the appropriate info.

!!! note 
    1. Like `supports` the all variable based query methods below also employ the 
       `label::Type{AbstractSupportLabel} = PublicLabel` keyword argument that by 
       default will return the desired information associated with public 
       supports. The full set (e.g., ones corresponding to internal collocation nodes) 
       is obtained via `label = All`.
    2. These methods also employ the `ndarray::Bool` keyword argument that will cause the 
       output to be formatted as a n-dimensional array where the dimensions 
       correspond to the infinite parameter dependencies. For example, if we have an 
       infinite variable `y(t, ξ)` and we invoke a query method with `ndarray = true` 
       then we'll get a matrix whose dimensions correspond to the supports of `t` and 
       `ξ`, respectively. Also, if `ndarray = true` then `label` correspond to the 
       intersection of supports labels in contrast to its default of invoking the union 
       of the labels.

First, we should verify that the optimized model in fact has variable values
via [`has_values`](@ref JuMP.has_values(::InfiniteModel)). In our example,
we have:
```jldoctest results
julia> has_values(model)
true
```
So we have values readily available to be extracted.

Now [`value`](@ref JuMP.value(::GeneralVariableRef)) can be used to query the
values as shown above in the Basic Usage section. This works by calling the
appropriate [`map_value`](@ref InfiniteOpt.map_value) defined by the optimizer model. By default this,
employs the `map_value` fallback which uses `optimizer_model_variable` to do the
mapping. Details on how to extend these methods for user-defined optimizer
models is explained on the Extensions page.

We also, support call to `value` that use an expression of variables as input.

Finally, the optimizer index of a variable is queried via
[`optimizer_index`](@ref JuMP.optimizer_index(::GeneralVariableRef)) which
reports back the index of the variable as used in the `MathOptInterface`
backend:
```jldoctest results
julia> optimizer_index(z)
MathOptInterface.VariableIndex(1)

julia> optimizer_index(x)
10-element Array{MathOptInterface.VariableIndex,1}:
 MathOptInterface.VariableIndex(2)
 MathOptInterface.VariableIndex(3)
 MathOptInterface.VariableIndex(4)
 MathOptInterface.VariableIndex(5)
 MathOptInterface.VariableIndex(6)
 MathOptInterface.VariableIndex(7)
 MathOptInterface.VariableIndex(8)
 MathOptInterface.VariableIndex(9)
 MathOptInterface.VariableIndex(10)
 MathOptInterface.VariableIndex(11)
```
As noted previously, an array is returned for `x(t)` in accordance with its
transcription variables. In similar manner to `value`, this is enabled by
appropriate versions of [`map_optimizer_index`](@ref InfiniteOpt.map_optimizer_index).

## Constraint Queries
Like variables, a variety of information can be queried about constraints.

!!! note 
    1. Like `supports` the all constraint query methods below also employ the 
       `label::Type{AbstractSupportLabel} = PublicLabel` keyword argument that by 
       default will return the desired information associated with public 
       supports. The full set (e.g., ones corresponding to internal collocation nodes) 
       is obtained via `label = All`.
    2. These methods also employ the `ndarray::Bool` keyword argument that will cause the 
       output to be formatted as a n-dimensional array where the dimensions 
       correspond to the infinite parameter dependencies. For example, if we have an 
       infinite constraint that depends on `t` and `ξ)` and we invoke a query method 
       with `ndarray = true` then we'll get a matrix whose dimensions correspond to 
       the supports of `t` and `ξ`, respectively. Also, if `ndarray = true` then 
       `label` correspond to the intersection of supports labels in contrast to its 
       default of invoking the union of the labels.

First, recall that constraints are stored in the form `function-in-set` where
generally `function` contains the variables and coefficients and the set contains
the relational operator and the constant value. With this understanding, we
query the value of a constraint's `function` via
[`value`](@ref JuMP.value(::InfOptConstraintRef)):
```jldoctest results
julia> constraint_object(c1).func # show the function expression of c1
z - x(t)

julia> value(c1)
10-element Array{Float64,1}:
 -8.747427671096375e-9
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
 20.999999995632077
```
Again, we obtain an array of values since `c1` is infinite due to its dependence
on `x(t)`. Behind the scenes this is implemented via the appropriate extensions
of [`map_value`](@ref InfiniteOpt.map_value).

Next the optimizer index(es) of the transcribed constraints in the
`MathOptInterface` backend provided via
[`optimizer_index`](@ref JuMP.optimizer_index(::InfOptConstraintRef)).
```jldoctest results
julia> optimizer_index(c1)
10-element Array{MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}},1}:
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(1)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(2)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(3)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(4)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(5)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(6)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(7)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(8)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(9)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(10)
```
Here 10 indices are given in accordance with the transcription constraints.
The mapping between these and the original infinite constraints is managed via
the appropriate extensions of [`map_optimizer_index`](@ref InfiniteOpt.map_optimizer_index).

We can also query dual information from our constraints if it is available.
First, we should verify that dual information is available via
[`has_duals`](@ref JuMP.has_duals(::InfiniteModel)):
```jldoctest results
julia> has_duals(model)
true
```
Now we can query the duals via [`dual`](@ref JuMP.dual(::InfOptConstraintRef)).
```jldoctest results
julia> dual(c1)
10-element Array{Float64,1}:
 1.9999999988666093
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
 1.1930560126841273e-10
```
Here we obtain the optimal dual values for each transcribed version of `c1`. This
is enabled via the proper extensions of [`map_dual`](@ref InfiniteOpt.map_dual).

Finally, we query the shadow price of a constraint via
[`shadow_price`](@ref JuMP.shadow_price(::InfOptConstraintRef)). This denotes
the change in the objective value due to an infinitesimal relaxation of the
constraint. For `c1` we get:
```jldoctest results
julia> shadow_price(c1)
10-element Array{Float64,1}:
 -1.9999999988666093
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
 -1.1930560126841273e-10
```
Similarly, the mapping to the transcription constraints is enabled via the
appropriate version of [`map_shadow_price`](@ref InfiniteOpt.map_shadow_price).

## LP Sensitivity
We also conduct sensitivity analysis for linear problems using
[`lp_sensitivity_report`](@ref JuMP.lp_sensitivity_report(::InfiniteModel)). This 
will generate a [`InfOptSensitivityReport`](@ref) which contains mapping to the 
ranges indicating how much a constraint RHS constant or a objective
coefficient can be changed without violating the feasibility of the solution.
This is further explained in the JuMP documentation
[here](https://jump.dev/JuMP.jl/stable/solutions/#Sensitivity-analysis-for-LP-1).
Furthermore, these analysis can only be employed for a solver that implements
`MOI.ConstraintBasisStatus`. In our running example up above, `Ipopt.jl` does not
support this A solver like `Gurobi.jl` does.
```julia-repl
julia> report = lp_sensitivity_report(model);

julia> report[c1]
10-element Array{Tuple{Float64,Float64},1}:
 (-42.0, Inf)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)

julia> report[z]
(-2.0, Inf)
```
Note that like other query methods, an array of ranges will be provided with
testing the sensitivity of an infinite constraint RHS in accordance with the
discretization scheme. Also, keyword arguments (like `ndarray` and `label`) can 
be invoked when indexing the report:
```julia-repl
julia> report[c1, label = All]
10-element Array{Tuple{Float64,Float64},1}:
 (-42.0, Inf)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
```

## Methods/DataTypes
```@index
Pages   = ["result.md"]
Modules = [JuMP, InfiniteOpt, InfiniteOpt.TranscriptionOpt]
Order   = [:function, :type]
```
```@docs
JuMP.termination_status(::InfiniteModel)
JuMP.raw_status(::InfiniteModel)
JuMP.primal_status(::InfiniteModel)
JuMP.dual_status(::InfiniteModel)
JuMP.solve_time(::InfiniteModel)
JuMP.has_values(::InfiniteModel)
JuMP.has_duals(::InfiniteModel)
JuMP.simplex_iterations(::InfiniteModel)
JuMP.barrier_iterations(::InfiniteModel)
JuMP.node_count(::InfiniteModel)
JuMP.objective_bound(::InfiniteModel)
JuMP.objective_value(::InfiniteModel)
JuMP.dual_objective_value(::InfiniteModel)
JuMP.result_count(::InfiniteModel)
JuMP.value(::GeneralVariableRef)
JuMP.value(::InfOptConstraintRef)
JuMP.value(::Union{JuMP.GenericAffExpr{<:Any, <:GeneralVariableRef}, JuMP.GenericQuadExpr{<:Any, <:GeneralVariableRef}})
JuMP.reduced_cost(::GeneralVariableRef)
JuMP.optimizer_index(::GeneralVariableRef)
JuMP.optimizer_index(::InfOptConstraintRef)
JuMP.dual(::InfOptConstraintRef)
JuMP.shadow_price(::InfOptConstraintRef)
JuMP.lp_sensitivity_report(::InfiniteModel)
InfOptSensitivityReport 
InfiniteOpt.map_value
InfiniteOpt.map_optimizer_index
InfiniteOpt.map_dual
InfiniteOpt.map_shadow_price
InfiniteOpt.map_reduced_cost
```
