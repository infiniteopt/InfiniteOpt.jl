```@meta
DocTestFilters = [r"≥|>=", r" == | = ", r" ∈ | in ", r" for all | ∀ "]
```

# [Optimization](@id opt_docs)
A guide for optimizing (solving) `InfiniteOpt` models. See the respective 
[technical manual](@ref opt_manual) for more details.

## Overview
`InfiniteOpt` offers a general and intuitive platform to optimize infinite 
optimization models. This is accomplished by applying a *transformation* to 
build a transformation backend that can be solved. By default, this is done via a 
[`TranscriptionBackend`](@ref) as described on the previous page. However, user-defined 
reformulation strategies can readily be implemented as described in the 
[Transformation Backends](@ref extend_backends) section on the extensions page. 

## Basic Usage
For most users, [`optimize!`](@ref JuMP.optimize!(::InfiniteModel)) is the only 
method required to optimize an `InfiniteModel`. This is exactly analogous 
to that of any `JuMP.Model` and is designed to provide a similar user experience. 
Let's first define an `InfiniteModel` with an appropriate optimizer:
```jldoctest optimize
julia> using InfiniteOpt, Ipopt;

julia> model = InfiniteModel(Ipopt.Optimizer);

julia> set_attribute(model, "print_level", 0);

julia> @infinite_parameter(model, t in [0, 10], num_supports = 10);

julia> @variable(model, y >= 0, Infinite(t));

julia> @variable(model, z >= 0);

julia> @objective(model, Min, 2z);

julia> @constraint(model, c1, z >= y);

julia> @constraint(model, c2, y(0) == 42);

julia> print(model)
Min 2 z
Subject to
 y(t) ≥ 0, ∀ t ∈ [0, 10]
 z ≥ 0
 c1 : z - y(t) ≥ 0, ∀ t ∈ [0, 10]
 y(0) ≥ 0
 c2 : y(0) = 42
```
Now we optimize the model using `optimize!`:
```jldoctest optimize
julia> optimize!(model);

julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4
```
Now our model has been solved, and we can query the solution. How to query the 
solution is explained on the [Results](@ref result_docs) page.

If no optimizer has been specified for the transformation backend, one can be provided 
via [`set_optimizer`](@ref):
```jldoctest; setup = :(using InfiniteOpt, Ipopt; model = InfiniteModel())
julia> set_optimizer(model, Ipopt.Optimizer)
```

A number of methods also exist to adjust the optimizer settings such as 
suppressing output. This is explained below in the 
[Optimizer Settings](@ref opt_settings) section.

## [Transformation Backends](@id opt_transform_backends)
As discussed previously, `InfiniteModel`s contain a transformation backend to 
store and solve a transformed version of the model. This backend typically
contains a data object (that stores a mapping between the transformed model and 
the infinite model). By default, a `TranscriptionBackend` is used which creates 
a discretized `JuMP.Model`. This transformation backend is used to optimize the 
infinite model, and it provides the information exacted by solution queries 
mapped back to the infinite model using the mapping data structure.

The process for optimizing an `InfiniteModel` is summarized in the following 
steps:
 1. fully define the `InfiniteModel`
 2. build the transformation backend via [`build_transformation_backend!`](@ref)
 3. optimize the `transformation_backend` via [`optimize!`](@ref JuMP.optimize!(::AbstractTransformationBackend)).

Here, `build_transformation_backend!` empties any existing backend information via 
[`empty!`](@ref Base.empty!(::AbstractTransformationBackend)), creates a reformulated 
of the `InfiniteModel`, stores it in in-place, and indicates that the transformation 
backend is ready via [`set_transformation_backend_ready`](@ref). These steps are all 
automated when [`optimize!`](@ref JuMP.optimize!(::InfiniteModel)) is invoked on the 
`InfiniteModel`.

The `transformation_backend` can be queried/extracted at any time from an `InfiniteModel` 
via [`transformation_backend`](@ref) and its underlying model (if there is one) is extracted
via [`transformation_model`](@ref transformation_model(::InfiniteModel)). For example, 
let's extract the transformation backend from the example above in the basic usage section: 
```jldoctest optimize
julia> backend = transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
Minimization problem with:
Variables: 11
Objective function type: AffExpr
`AffExpr`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`AffExpr`-in-`MathOptInterface.GreaterThan{Float64}`: 10 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 11 constraints
Model mode: AUTOMATIC
CachingOptimizer state: ATTACHED_OPTIMIZER
Solver name: Ipopt

julia> tmodel = transformation_model(model)
A JuMP Model
Minimization problem with:
Variables: 11
Objective function type: AffExpr
`AffExpr`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`AffExpr`-in-`MathOptInterface.GreaterThan{Float64}`: 10 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 11 constraints
Model mode: AUTOMATIC
CachingOptimizer state: ATTACHED_OPTIMIZER
Solver name: Ipopt
```

The `JuMP` variable(s) stored in the transformation backend that correspond to a 
particular `InfiniteOpt` variable can be queried via 
[`transformation_variable`](@ref transformation_variable(::GeneralVariableRef)). 
Thus, using the going example we get:
```jldoctest optimize
julia> transformation_variable(y) # infinite variable
10-element Vector{VariableRef}:
 y(support: 1)
 y(support: 2)
 y(support: 3)
 y(support: 4)
 y(support: 5)
 y(support: 6)
 y(support: 7)
 y(support: 8)
 y(support: 9)
 y(support: 10)

julia> transformation_variable(z) # finite variable
z
```
In like manner, we get the `JuMP` constraints corresponding to a particular 
`InfiniteOpt` constraint via 
[`transformation_constraint`](@ref transformation_constraint(::InfOptConstraintRef)). 
Thus, using going example we get: 
```jldoctest optimize
julia> transformation_constraint(c1) # infinite constraint
10-element Vector{ConstraintRef}:
 c1(support: 1) : z - y(support: 1) ≥ 0
 c1(support: 2) : z - y(support: 2) ≥ 0
 c1(support: 3) : z - y(support: 3) ≥ 0
 c1(support: 4) : z - y(support: 4) ≥ 0
 c1(support: 5) : z - y(support: 5) ≥ 0
 c1(support: 6) : z - y(support: 6) ≥ 0
 c1(support: 7) : z - y(support: 7) ≥ 0
 c1(support: 8) : z - y(support: 8) ≥ 0
 c1(support: 9) : z - y(support: 9) ≥ 0
 c1(support: 10) : z - y(support: 10) ≥ 0
```
We can also query the expressions via 
[`transformation_expression`](@ref transformation_expression(::JuMP.AbstractJuMPScalar)):
```jldoctest optimize
julia> transformation_expression(z - y^2 + 3) # infinite expression
10-element Vector{AbstractJuMPScalar}:
 -y(support: 1)² + z + 3
 -y(support: 2)² + z + 3
 -y(support: 3)² + z + 3
 -y(support: 4)² + z + 3
 -y(support: 5)² + z + 3
 -y(support: 6)² + z + 3
 -y(support: 7)² + z + 3
 -y(support: 8)² + z + 3
 -y(support: 9)² + z + 3
 -y(support: 10)² + z + 3
```

!!! note 
    1. Like `supports` the `transformation_[obj]` methods also employ the 
       `label::Type{AbstractSupportLabel} = PublicLabel` keyword argument that by 
       default will return variables/expressions/constraints associated with public 
       supports. The full set (e.g., ones corresponding to internal collocation nodes) 
       is obtained via `label = All`.
    2. These methods also employ the `ndarray::Bool` keyword argument that will cause the 
       output to be formatted as an n-dimensional array where the dimensions 
       correspond to the infinite parameter dependencies. For example, if we have an 
       infinite variable `y(t, ξ)`, and we invoke a query method with `ndarray = true` 
       then we'll get a matrix whose dimensions correspond to the supports of `t` and 
       `ξ`, respectively. Also, if `ndarray = true` then `label` correspond to the 
       intersection of supports labels in contrast to its default of invoking the union 
       of the labels.

The purpose of this `transformation_backend` abstraction is to readily enable user-defined 
reformulation extensions (e.g., using polynomial chaos expansion theory). However, 
this is all handled behind the scenes such that most users can interact with 
`InfiniteModel`s like any `JuMP.Model`.

## [Optimizer Settings](@id opt_settings)
A few optimizer settings can be set consistently agnostic of particular 
solver keywords. One such setting is that of suppressing
optimizer verbose output. This is accomplished via 
[`set_silent`](@ref JuMP.set_silent(::InfiniteModel)) and 
[`unset_silent`](@ref JuMP.unset_silent(::InfiniteModel)). The syntax is 
exemplified below:
```jldoctest optimize
julia> set_silent(model)

julia> unset_silent(model)
```

We can also adjust the time limit in a solver independent fashion via 
[`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref), and 
[`time_limit_sec`](@ref). These methods are illustrated below:
```jldoctest optimize
julia> set_time_limit_sec(model, 100)

julia> time_limit_sec(model)
100.0

julia> unset_time_limit_sec(model)
```

Other optimizer specific settings can be set via [`set_attribute`](@ref). 
For example, let's set the maximum CPU time for Ipopt:
```jldoctest optimize
julia> set_attribute(model, "max_cpu_time", 60.)
```
Multiple settings can be specified via [`set_attributes`](@ref). For 
example, let's specify the tolerance and the maximum number of iterations:
```jldoctest optimize
julia> set_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```

Finally, we can query optimizer settings via [`get_attribute`](@ref). 
For example, let's query the maximum number of iterations:
```jldoctest optimize
julia> get_attribute(model, "max_iter")
100
```
Note this only works if the attribute has been previously specified.
