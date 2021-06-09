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
One key idea to keep in mind is that the objective must evaluate to a finite 
expression. Note this means that objectives can only explicitly contain 
finite variables and point variables. Infinite expressions must be evaluated in a 
measure to be included (e.g., evaluate the expectation of a random variable).

!!! note 
    Nonlinear objectives as defined by `JuMP.@NLobjective` are not currently 
    supported by `InfiniteOpt`. See [Nonlinear Expressions](@ref) for more 
    information and possible workarounds. 

## [Basic Usage] (@id obj_basic)
Principally, the objective function is specified via 
[`@objective`](https://jump.dev/JuMP.jl/v0.21.8/reference/objectives/#JuMP.@objective) 
as is done in `JuMP`. For example, let's define the stochastic objective to 
minimize ``0.5 x_1 + 0.5 x_2 + \mathbb{E}_\xi [y^2 - y]``:
```jldoctest obj; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> @infinite_parameter(model, ξ ~ Normal())
ξ

julia> @variable(model, y, Infinite(ξ))
y(ξ)

julia> @variable(model, x[1:2])
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]

julia> @objective(model, Min, 0.5x[1] + 0.5x[2] + 𝔼(y^2 - y, ξ))
0.5 x[1] + 0.5 x[2] + 𝔼{ξ}[y(ξ)² - y(ξ)]
```
Thus, we have defined an objective using `InfiniteOpt`'s straightforward syntax. 
Note that the second argument indicates the objective sense which can be 
expressed `Min` for minimization problems and `Max` for maximization problems. 
The objective function (expression) must be finite containing only finite variables, 
point variables, and/or measures. Also, any included measures must fully 
integrate over all the infinite parameters contained in its input function. 
For example, if we define had an infinite variable `z(ξ, t)` then the measure 
`𝔼(z, ξ)` could not be included since the resulting expression would still 
be infinite with respect to `t`. However, adding a measure for `t` would result 
in a valid object to add to an objective: `∫(𝔼(z, ξ), t)`.

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
GenericAffExpr{Float64,GeneralVariableRef}
```
The objective sense can be one of three possibilities: `MIN_SENSE`, `MAX_SENSE`, 
or `FEASIBILITY_SENSE`. The later sense applies to models that contain no 
objective function.

## Modification
This section will review the methods that can be used to modify the objective. 
First, we'll consider the useful 
[`set_objective_coefficient`](@ref JuMP.set_objective_coefficient(::InfiniteModel, ::GeneralVariableRef, ::Real)) 
method and then we'll explain the methods that enable `@objective`.

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
    [`@objective`](https://jump.dev/JuMP.jl/v0.21.8/reference/objectives/#JuMP.@objective) 
    since it is more stable and efficient than the `set_objective_[aspect]` 
    methods due to its enhanced methodology for parsing expressions.
