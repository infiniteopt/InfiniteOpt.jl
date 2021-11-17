# [Infinite Models](@id infinite_model_docs)
A guide for creating infinite models and understanding how they work. See the 
respective [technical manual](@ref infinite_model_manual) for more details.

## Overview
Infinite models are expressed via the [`InfiniteModel`](@ref) datatype which is at the
core of `InfiniteOpt`. These model objects are designed to emulate the behavior
of [`Model`](https://jump.dev/JuMP.jl/v0.22/reference/models/#JuMP.Model) 
objects in `JuMP`. These data objects store the parameters, variables,
measures, objective, constraints, and all other data used in `InfiniteOpt`. This
differs from `JuMP` models which store such information in a `MathOptInterface`
model backend.

!!! note 
    `InfiniteOpt`'s `InfiniteModel`s are intended to be used for 
    infinite-dimensional optimization problems. Finite problems (e.g., 
    directly modeling a discrete time model) should instead be modeled using 
    `Model`'s in [`JuMP`](https://jump.dev/JuMP.jl/stable/).

## Basic Usage
Infinite models can be initialized with no arguments by default:
```jldoctest
julia> using InfiniteOpt

julia> model = InfiniteModel()
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```

The optimizer that will be used to solve the model can also be specified at
model definition:
```jldoctest
julia> using InfiniteOpt, Ipopt

julia> model = InfiniteModel(Ipopt.Optimizer)
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```
Note that any optimizer currently supported by `JuMP v0.19.0` or newer is
supported for use in `InfiniteOpt`. For completeness, the table of currently
supported optimizers is provided below in [Supported Optimizers](@ref).

We can also specify optimizer attributes via
[`optimizer_with_attributes`](https://jump.dev/JuMP.jl/v0.22/reference/models/#JuMP.optimizer_with_attributes)
which allows us to append as many attributes as we like, for example:
```jldoctest
julia> using InfiniteOpt, Ipopt

julia> model = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer,
                                                       "output_level" => 0))
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```

Now you have an initialized `InfiniteModel` that is ready for your mathematical
model to be defined and optimized!

## Advanced Definition Information
As you may have noticed in the above examples, `InfiniteModel`s contain an
optimizer model backend which simply corresponds to a `JuMP.Model` that
will be used to store and optimize the reformulation of the infinite mathematical
model stored in `InfiniteModel`. It also will contain a mapping between its
optimization model and that of the `InfiniteModel` (e.g., a mapping between the
variables and constraints). By default, `InfiniteModel`s use a
[`TranscriptionModel`](@ref) optimizer model backend which will store a
transcribed (discretized) version of the infinite model. More information on
the internal use of `TranscriptionModel`s is provided in
[Model Transcription](@ref transcription_docs).

All the arguments used with the `InfiniteModel` constructor (e.g., the optimizer)
are simply passed on and stored in the optimizer model backend. Thus, any
argument supported by `JuMP.Model` can be passed on to the optimizer
model by including it in the `InfiniteModel` constructor. For example, we can
specify the `add_bridges` keyword argument in the `InfiniteModel` call to use
in the definition of the optimizer model:
```jldoctest
julia> using InfiniteOpt, Ipopt

julia> model = InfiniteModel(Ipopt.Optimizer,
                             add_bridges = false)
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```

Moreover, alternative optimizer model types (i.e., not a `TranscriptionModel`) can be 
specified via the `OptimizerModel` keyword argument when initializing the 
`InfiniteModel`. Thus, to redundantly specify a `TranscriptionModel` we would call:
```jldoctest model_fun
julia> using InfiniteOpt

julia> model = InfiniteModel(OptimizerModel = TranscriptionModel)
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
More information on implementing custom optimizer models is located on the 
Extensions page.

## Supported Optimizers
`InfiniteOpt` can use any optimizer that is supported by `JuMP v0.19.0` or newer 
(i.e., has a `MathOptInterface` implementation). Please refer to `JuMP`'s current
[solver documentation](https://jump.dev/JuMP.jl/v0.22/installation/#Supported-solvers) 
to learn what solvers are supported and how to install them.

## Object Dictionaries
Like `JuMP.Model`s, `InfiniteModel`s register the name symbols of macro defined 
objects. This enables us to access such objects by indexing the `InfiniteModel` 
with the appropriate symbol. This is particularly useful for function defined 
models. For example:
```jldoctest; setup = :(using InfiniteOpt)
julia> function make_model(num_supports)
        model = InfiniteModel()
        @infinite_parameter(model, t ∈ [0, 10], num_supports = num_supports)
        @variable(model, y >= 42, Infinite(t))
        @objective(model, Min, ∫(y, t))
        return model
       end
make_model (generic function with 1 method)

julia> model1 = make_model(2); model2 = make_model(4);

julia> y1 = model1[:y]
y(t)
```
Note that when macro defined objects are deleted from an `InfiniteModel` that the 
corresponding symbols in the object dictionary are not removed by default. This 
can be accomplished by use of 
[`JuMP.unregister`](https://jump.dev/JuMP.jl/v0.22/reference/models/#JuMP.unregister) 
(please click on its link for usage information).
