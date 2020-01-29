# Optimization
A guide and manual for optimizing (solving) `InfiniteOpt` models. The Methods
section at the end comprise the manual, and the above sections comprise the guide.

## Overview
Fundamentally, we seek to optimize a given infinite optimization model that
we have defined and this is the very reason why `InfiniteOpt` was created. Thus,
`InfiniteOpt` offers a general and intuitive platform to do just this. This
is made up of transforming the `InfiniteModel` into a standard optimization
problem stored as a `JuMP.Model` (referred to as the `optimizer_model`) that is
then optimized via a compatible optimizer. By default this is done via a
`TranscriptionModel` as described on the previous page. However, user-defined
reformulation strategies can readily be implemented as described in the
[Optimizer Models](@ref extend_optimizer_model) section on the extensions page.

## Basic Usage
For most users, [`optimize!`](@ref JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory}))
is the method required to optimize an `InfiniteModel`. This is exactly analogous
to that of any `JuMP.Model` and is designed to provide a similar user experience.
Let's first define an `InfiniteModel` with an appropriate optimizer:
```jldoctest optimize
julia> using InfiniteOpt, JuMP, Ipopt;

julia> model = InfiniteModel(with_optimizer(Ipopt.Optimizer, print_level = 0));

julia> @infinite_parameter(model, t in [0, 10], num_supports = 10);

julia> @infinite_variable(model, x(t) >= 0);

julia> @hold_variable(model, z >= 0);

julia> @objective(model, Min, 2z);

julia> @constraint(model, z >= x);

julia> @BDconstraint(model, t == 0, x == 42);

julia> print(model)
Min 2 z
Subject to
 x(t) ≥ 0.0
 z ≥ 0.0
 z - x(t) ≥ 0.0
 x(t) = 42.0, ∀ t = 0
 t ∈ [0, 10]
```
Now we optimize the model using `optimize!`:
```jldoctest optimize
julia> optimize!(model)

julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4
```
Now our model has been solved and we can query the solution. How to query the
solution is explained on the [Results](@ref) page.

If no optimizer has been specified for the `InfiniteModel`, one can be provided
via [`set_optimizer`](@ref):
```jldoctest; setup = :(using InfiniteOpt, JuMP, Ipopt; model = InfiniteModel())
julia> set_optimizer(model, with_optimizer(Ipopt.Optimizer))
```

A number of methods also exist to adjust the optimizer settings such as
suppressing output. This is explained below in the [Optimizer Settings](@ref)
section.

## Optimizer Models
As discussed previously, `InfiniteModel`s contain an `optimizer_model` field
which stores a transformed finite version of the model in a `JuMP.Model` that
contains a data struct (that stores a mapping between the transformed model and
the infinite model) in the `Model.ext` dictionary with an associated key. By
default a `JuMP.Model` using [`TranscriptionData`](@ref) stored under the key
`:TransData` is used and is referred to as a `TranscriptionModel`. The
optimizer model is then what is used to optimize the infinite model and it provides
the information exacted by solution queries mapped backed back to the infinite
model using the mapping data structure.

The process for optimizing an `InfiniteModel` is summarized in the following
steps:
 - fully define the `InfiniteModel`
 - build the optimizer model via [`build_optimizer_model!`](@ref)
 - optimize the `optimizer_model` via [`optimize!`](@ref JuMP.optimize!(::JuMP.Model, ::Union{Nothing, JuMP.OptimizerFactory})).

Here `build_optimizer_model!` creates a reformulated finite version of the
`InfiniteModel`, stores it in `InfiniteModel.optimizer_model` via
[`set_optimizer_model`](@ref), and indicates that the optimizer model is ready
via [`set_optimizer_model_ready`](@ref). These steps are all automated when
[`optimize!`](@ref JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory}))
is invoked on the `InfiniteModel`.

The `optimizer_model` can be queried/extracted at any time from an `InfiniteModel`
via [`optimizer_model`](@ref). For example, let's extract the optimizer model
from the example above in the basic usage section:
```jldoctest optimize
julia> trans_model = optimizer_model(model)
A JuMP Model
Minimization problem with:
Variables: 11
Objective function type: GenericAffExpr{Float64,VariableRef}
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.GreaterThan{Float64}`: 10 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 11 constraints
Model mode: AUTOMATIC
CachingOptimizer state: ATTACHED_OPTIMIZER
Solver name: Ipopt
```

The purpose of this `optimizer_model` structure is to readily enable user-defined
reformulation extensions (e.g., using polynomial chaos expansion theory). However,
this is all handled behind the scenes such that most users can interact with
`InfiniteModel`s like any `JuMP.Model`.

## Optimizer Settings
A few optimizer settings can be set in a consistent way agnostic of particular
solver keywords. One such setting is that of suppressing and unsuppressing
optimizer verbose output. This is accomplished via
[`set_silent`](@ref JuMP.set_silent(::InfiniteModel)) and
[`unset_silent`](@ref JuMP.unset_silent(::InfiniteModel)). The syntax is
exemplified below:
```jldoctest optimize
julia> set_silent(model)
true

julia> unset_silent(model)
false
```

Other optimizer specific settings can be set via
[`set_parameter`](@ref JuMP.set_parameter(::InfiniteModel, ::Any, ::Any)).
For example, let's set the maximum CPU time for Ipopt:
```jldoctest optimize
julia> set_parameter(model, "max_cpu_time", 60.)
60.0
```

## Methods
```@index
Pages   = ["optimize.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:function]
```
```@docs
optimizer_model
set_optimizer_model
optimizer_model_key
build_optimizer_model!(::InfiniteModel)
build_optimizer_model!
optimizer_model_ready
set_optimizer_model_ready
JuMP.bridge_constraints(::InfiniteModel)
JuMP.add_bridge(::InfiniteModel, ::Type{<:MOI.Bridges.AbstractBridge})
JuMP.set_optimizer(::InfiniteModel, ::JuMP.OptimizerFactory)
JuMP.set_silent(::InfiniteModel)
JuMP.unset_silent(::InfiniteModel)
JuMP.set_parameter(::InfiniteModel, ::Any, ::Any)
JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory})
```
