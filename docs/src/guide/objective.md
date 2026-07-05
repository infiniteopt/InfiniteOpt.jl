```@meta
DocTestFilters = [r"E|𝔼", r"integral|∫"]
```

# [Objectives](@id obj_docs)
A guide specifying and modifying objective functions in `InfiniteOpt`. See the 
respective [technical manual](@ref obj_manual) for more details.

## Overview
Naturally, objective functions serve as a key aspect of optimization problems in 
general and this is certainly the case with infinite dimensional ones. In 
`InfiniteOpt` objectives are defined in much the same way they are in `JuMP`. 
One key idea is that the objective must evaluate to a finite expression which means 
it must only explicitly contain finite variables and point variables. Infinite 
expressions must be summarized by a measure (e.g., taking the expectation of a random variable).


## [Basic Usage](@id obj_basic)
Principally, the objective function is specified via 
[`@objective`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@objective) 
as is done in `JuMP`. For example, let's define the stochastic objective to 
minimize ``0.5 x_1 + 0.5 x_2 + \mathbb{E}_\xi [y^2 - y]``:
```jldoctest obj; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> @infinite_parameter(model, ξ ~ Normal())
ξ

julia> @variable(model, y, Infinite(ξ))
y(ξ)

julia> @variable(model, x[1:2])
2-element Vector{GeneralVariableRef}:
 x[1]
 x[2]

julia> @objective(model, Min, 0.5x[1] + 0.5x[2] + 𝔼(y^2 - y, ξ))
0.5 x[1] + 0.5 x[2] + 𝔼{ξ}[y(ξ)² - y(ξ)]
```
Note that the second argument indicates the objective sense which can be 
expressed `Min` for minimization problems and `Max` for maximization problems. 
The objective function (expression) must be finite containing only finite variables, 
point variables, and/or measures. Also, any included measures must fully summarize 
all infinite parameters contained in the expression they summarize. 
For example, if we define an infinite variable `z(ξ, t)`, then the measure 
`𝔼(z, ξ)` could not be included since the resulting expression would still 
be infinite with respect to `t`. However, adding a measure for `t` would result 
in a valid object to add to an objective: `∫(𝔼(z, ξ), t)`.

!!! note 
    Nonlinear objectives are defined simply by using `@objective` and not 
    using `JuMP.@NLobjective`. See [Nonlinear Expressions](@ref nlp_guide) for 
    more information. 

Now we can add objectives to our infinite models. For more detailed information, 
please review the information below.  

## Queries
This section will highlight the available methods for extracting objective 
information. These are all derived from extensions to `JuMP` functions and thus 
follow syntax.

Principally, these methods correspond to 
[`objective_sense`](@ref JuMP.objective_sense(::InfiniteModel)), 
[`objective_function`](@ref JuMP.objective_function(::InfiniteModel)), and 
[`objective_function_type`](@ref JuMP.objective_function_type(::InfiniteModel)) 
which return the objective sense (a subtype of `MOI.OptimizationSense`), the 
objective function (expression), and the objective function type, respectively. 
These methods are demonstrated in accordance with the example presented above in 
the [Basic Usage](@ref obj_basic) section:
```jldoctest obj
julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0

julia> objective_function(model)
0.5 x[1] + 0.5 x[2] + 𝔼{ξ}[y(ξ)² - y(ξ)]

julia> objective_function_type(model)
GenericAffExpr{Float64, GeneralVariableRef}
```
The objective sense can be one of three possibilities: `MIN_SENSE`, `MAX_SENSE`, 
or `FEASIBILITY_SENSE`. The later sense applies to models that contain no 
objective function.

## Modification
This section will review the methods that can be used to modify the objective. 
First, we'll consider the useful 
[`set_objective_coefficient`](@ref JuMP.set_objective_coefficient(::InfiniteModel, ::GeneralVariableRef, ::Real)) 
method, and then we'll explain the methods that enable `@objective`.

The coefficient of a particular variable in an objective can be readily updated 
via [`set_objective_coefficient`](@ref JuMP.set_objective_coefficient(::InfiniteModel, ::GeneralVariableRef, ::Real)). 
This is useful repeatedly optimizing an infinite model with varied objective 
coefficients (e.g., varied tradeoff parameters). For example, let's consider 
updating the coefficient of `x[1]` in the previous example from 0.5 to 0.25:
```jldoctest obj
julia> set_objective_coefficient(model, x[1], 0.25)

julia> objective_function(model)
0.25 x[1] + 0.5 x[2] + 𝔼{ξ}[y(ξ)² - y(ξ)]
```

Now let's consider the modification methods that enable the `@objective` macro. 
The objective function is specified via 
[`set_objective_function`](@ref JuMP.set_objective_function(::InfiniteModel, ::JuMP.AbstractJuMPScalar)) 
which simply updates the expression stored in the objective. For example, 
let's update out objective to simply be ``0.5x_1 + 0.5x_2``:
```jldoctest obj
julia> set_objective_function(model, 0.5x[1] + 0.5x[2])

julia> objective_function(model)
0.5 x[1] + 0.5 x[2]
```

The objective sense is updated via 
[`set_objective_sense`](@ref JuMP.set_objective_sense(::InfiniteModel, ::MOI.OptimizationSense)) 
which can specify the sense as one of the `MOI.OptimizationSense` subtypes. For 
example, let's change the current objective to be maximization problem:
```jldoctest obj
julia> set_objective_sense(model, MOI.MAX_SENSE)

julia> objective_sense(model)
MAX_SENSE::OptimizationSense = 1
```

The above 2 methods are both called via 
[`set_objective`](@ref JuMP.set_objective(::InfiniteModel, ::MOI.OptimizationSense, ::Union{JuMP.AbstractJuMPScalar, Real})). 
This is the function that enables `@objective` behind the scenes. Thus, the 
previous 2 examples could have been implemented equivalently in the following 
ways:
```jldoctest obj
julia> set_objective(model, MOI.MAX_SENSE, 0.5x[1] + 0.5x[2])

julia> @objective(model, Max, 0.5x[1] + 0.5x[2])
0.5 x[1] + 0.5 x[2]
```
Notice that `@objective` offers a more intuitive syntax and is also 
more efficient at parsing expressions.

!!! note
    When possible, the 
    [`@objective`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@objective) 
    since it is more stable and efficient than the `set_objective_[aspect]` 
    methods due to its enhanced methodology for parsing expressions.

### Multi-Objective Problems
`InfiniteOpt` also supports multi-objective problems by supplying the model with a vector of scalar values.
Instead of a `JuMP.AbstractJuMPScalar`, these become `AbstractVector` within the model.

These models are defined the same way as scalar, single-objective models with `@objective`, with some slight differences in definition. Looking at the aforementioned example and modifying it a bit, it becomes:

```jldoctest obj
julia> set_objective(model, MOI.MAX_SENSE, [0.5x[1], 0.5x[2]])
 0.5 x[1]
 0.5 x[2]
```
If a scalarization method of solving the model is the goal, the individual `@expression` scalars within the vector can be compiled into a single objective. The weight- or priority- of these expressions can be changed by multiplying each component in the vector by fractions, with the sum of the fractions being 1.
In hierarchical model solving, the order of the scalars determines the priority of the model.
For more solvers that support multiple objective optimization, please refer to `JuMP` solvers, such as [`@MultiObjectiveAlgorithms`](https://jump.dev/JuMP.jl/stable/packages/MultiObjectiveAlgorithms/#MultiObjectiveAlgorithms.jl) and [`@Gurobi`](https://jump.dev/JuMP.jl/stable/packages/Gurobi/#Gurobi.jl).