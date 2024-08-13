# [Infinite Models](@id infinite_model_docs)
A guide for creating infinite models and understanding how they work. See the 
respective [technical manual](@ref infinite_model_manual) for more details.

## Overview
Infinite models are expressed via the [`InfiniteModel`](@ref) datatype which is at the
core of `InfiniteOpt`. These model objects are designed to emulate the behavior
of [`Model`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.Model) 
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
  Finite parameters: 0
  Infinite parameters: 0
  Variables: 0
  Derivatives: 0
  Measures: 0
Transformation backend information:
  Backend type: TranscriptionBackend
  Solver: none
  Transformation built and up-to-date: false
```
Ultimately, `model` will be solved via a transformation backend. By default, 
we see that a [`TranscriptionBackend`](@ref) (a JuMP `Model` that will contain 
a transcribed, i.e., a discretized model) is used. To specify, a backend 
of our choice, we use the syntax:
```jldoctest
julia> using InfiniteOpt

julia> model = InfiniteModel(TranscriptionBackend())
An InfiniteOpt Model
Feasibility problem with:
  Finite parameters: 0
  Infinite parameters: 0
  Variables: 0
  Derivatives: 0
  Measures: 0
Transformation backend information:
  Backend type: TranscriptionBackend
  Solver: none
  Transformation built and up-to-date: false
```

Since `TranscriptionBackend`s are a common choice, we can just pass a JuMP 
compatible optimizer (i.e., solver) to the model and a `TranscriptionBackend` 
that uses that optimizer will be initialized:
```jldoctest
julia> using InfiniteOpt, Ipopt

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
For completeness, the table of currently supported JuMP compatible optimizers 
is provided below in [Supported Optimizers](@ref).

We can also specify optimizer attributes via
[`optimizer_with_attributes`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.optimizer_with_attributes)
which allows us to append as many attributes as we like, for example:
```jldoctest
julia> using InfiniteOpt, Ipopt

julia> jump_opt = optimizer_with_attributes(Ipopt.Optimizer, "output_level" => 0);

julia> model = InfiniteModel(jump_opt)
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

Now you have an initialized `InfiniteModel` that is ready for your mathematical
model to be defined and optimized!

## Advanced Definition Information
As noted above, `InfiniteModel`s contain a transformation backend that will ultimately 
be used to optimize the `InfiniteModel` via a transformed version of it. Such backends 
typically have methods to transform an `InfiniteModel` into a transformed model that 
can be optimized; moreover, they store necessary data to map back to the `InfiniteModel`. 

By default, `InfiniteModel`s use a [`TranscriptionBackend`](@ref) which will store a
transcribed (i.e., discretized) version of the infinite model. More information on
`TranscriptionBackends`s is provided in [Model Transcription](@ref transcription_docs).
Notably, the main argument `TranscriptionBackend` is an appropriate JuMP compatible 
optimizer:
```jldoctest
julia> using InfiniteOpt, Ipopt

julia> backend = TranscriptionBackend(Ipopt.Optimizer)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: Ipopt
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```

We query the underlying transformation backend, transformation model, and transformation
data via [`transformation_backend`](@ref), 
[`transformation_model`](@ref transformation_model(::InfiniteModel)), and
[`transformation_data`](@ref transformation_data(::InfiniteModel)), respectively:
```jldoctest
julia> using InfiniteOpt; model = InfiniteModel();

julia> tbackend = transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> tmodel = transformation_model(model)
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> data = transformation_data(model);
```

A new transformation backend is specified via [`set_transformation_backend`](@ref):
```jldoctest
julia> using InfiniteOpt, Ipopt; model = InfiniteModel();

julia> set_transformation_backend(model, TranscriptionBackend(Ipopt.Optimizer))

julia> tbackend = transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: Ipopt
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```
Again, since `TranscriptionBackend` is the default, the following models are equivalent:
```jldoctest
julia> using InfiniteOpt, Ipopt; 

julia> model1 = InfiniteModel();

julia> set_transformation_backend(model1, TranscriptionBackend(Ipopt.Optimizer, add_bridges = false))

julia> model2 = InfiniteModel(Ipopt.Optimizer, add_bridges = false)
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

More information on implementing custom transformation backends is located on the 
Extensions page.

## Supported Optimizers
Supported optimizers (e.g., solvers) depend on the transformation backend being 
used. For [`JuMPBackend`](@ref)s such as [`TranscriptionBackend`](@ref), any 
JuMP compatible optimizer (i.e., has a `MathOptInterface` implementation) can be 
used. Please refer to `JuMP`'s current
[solver documentation](https://jump.dev/JuMP.jl/v1/installation/#Supported-solvers) 
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
[`JuMP.unregister`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.unregister) 
(please click on its link for usage information).
