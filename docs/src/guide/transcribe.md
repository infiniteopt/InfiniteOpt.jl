# Model Transcription
A guide and manual for transcribing infinite models using `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.

## Overview
All infinite models need to be reformulated in such a way that they can be solved
using traditional optimization methods. Typically, this involves discretization
of the infinite domain via particular parameter support points. By default,
`InfiniteOpt` employs this methodology via the use of transcription models (which
comprise the `optimizer_model` as discussed in the [Infinite Models](@ref) section).
`InfiniteOpt` is built modularly to readily accept other user defined techniques
and this is discussed in further detail on the [Extensions](@ref) page. This
page will detail transcription models based in `InfiniteOpt.TranscriptionOpt`
which provide the default transcription (reformulation) capabilities of
`InfiniteOpt`.

## Basic Usage
Most users will not need to employ the capabilities of `TranscriptionOpt` directly
since they are employed implicitly with the call of
[`optimize!`](@ref JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory}))
on an infinite model. This occurs since `TranscriptionModel`s are the default
optimizer model type that is employed.

However, some users may wish to use `TranscriptionOpt` to extract a fully
discretized/transcribed version of an infinite model that is conveniently output
as a typical `JuMP` model and can then be treated as such. This is principally
accomplished via the [`TranscriptionModel`](@ref TranscriptionModel(::InfiniteModel))
constructor. To illustrate how this is done, let's first define a basic infinite
model with a simple support structure for the sake of example:
```jldoctest transcribe
julia> using InfiniteOpt, JuMP

julia> inf_model = InfiniteModel();

julia> @infinite_parameter(inf_model, t in [0, 10], supports = [0, 5, 10])
t

julia> @infinite_variable(inf_model, g(t) >= 0)
g(t)

julia> @hold_variable(inf_model, z, Bin)
z

julia> @objective(inf_model, Min, 2z + support_sum(g, t))
2 z + sum(g(t))

julia> @BDconstraint(inf_model, initial(t == 0), g == 1)
initial : g(t) = 1.0, ∀ t = 0

julia> @constraint(inf_model, constr, g^2 - z <= 42)
constr : g(t)² - z ≤ 42.0

julia> print(inf_model)
Min 2 z + sum(g(t))
Subject to
 g(t) ≥ 0.0
 z binary
 g(t) = 1.0, ∀ t = 0
 g(t)² - z ≤ 42.0
 t ∈ [0, 10]
```
Now we can make `JuMP` model containing the transcribed version of `inf_model`
via [`TranscriptionModel`](@ref TranscriptionModel(::InfiniteModel)):
```jldoctest transcribe
julia> trans_model = TranscriptionModel(inf_model)
A JuMP Model
Minimization problem with:
Variables: 4
Objective function type: GenericAffExpr{Float64,VariableRef}
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`GenericQuadExpr{Float64,VariableRef}`-in-`MathOptInterface.LessThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.ZeroOne`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> print(trans_model)
Min 2 z + g(support: 1) + g(support: 2) + g(support: 3)
Subject to
 g(support: 1) = 1.0
 g(support: 1)² - z ≤ 42.0
 g(support: 2)² - z ≤ 42.0
 g(support: 3)² - z ≤ 42.0
 g(support: 1) ≥ 0.0
 g(support: 2) ≥ 0.0
 g(support: 3) ≥ 0.0
 z binary
```
Thus, we have a transcribed `JuMP` model. To be precise this actually a
`TranscriptionModel` which is a `JuMP.Model` with some extra data stored in the
`ext` field that retains the mapping between the transcribed variables and
constraints and their infinite counterparts. Notice, that multiple finite variables
have been introduced to discretize `g(t)` at supports 1, 2, and 3 which correspond
to 0, 5, and 10 as can be queried by
[`supports`](@ref InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteVariableRef)):
```jldoctest transcribe
julia> supports(trans_model, g)
3-element Array{Tuple{Float64},1}:
 (0.0,)
 (5.0,)
 (10.0,)
```
Also, notice how the constraints are transcribed in accordance with these supports
except the initial condition which naturally is only invoked for the first support
point. Furthermore, the transcription variable(s) of any variable associated with
the infinite model can be determined via [`transcription_variable`](@ref):
```jldoctest transcribe
julia> transcription_variable(trans_model, g)
3-element Array{VariableRef,1}:
 g(support: 1)
 g(support: 2)
 g(support: 3)

julia> transcription_variable(trans_model, z)
z
```
Similarly, the transcription constraints associated with infinite model constraints
can be queried via [`transcription_constraint`](@ref) and the associated supports
and infinite parameters can be found via
[`supports`](@ref InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef))
and [`parameter_refs`](@ref InfiniteOpt.parameter_refs(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef)):
```jldoctest transcribe
julia> transcription_constraint(trans_model, initial)
1-element Array{ConstraintRef,1}:
 initial(Support: 1) : g(support: 1) = 1.0

julia> transcription_constraint(trans_model, constr)
3-element Array{ConstraintRef,1}:
 constr(Support: 1) : g(support: 1)² - z ≤ 42.0
 constr(Support: 2) : g(support: 2)² - z ≤ 42.0
 constr(Support: 3) : g(support: 3)² - z ≤ 42.0

julia> supports(trans_model, constr)
3-element Array{Tuple{Float64},1}:
 (0.0,)
 (5.0,)
 (10.0,)

julia> parameter_refs(trans_model, constr)
(t,)
```
Note the parameter reference tuple corresponds to the support tuples.

Now we have a transcribed `JuMP` model that can be optimized via traditional
`JuMP` methods whose variables and constraints can be accessed using the methods
mentioned above.

## Transcription Theory
A given infinite dimensional optimization problem is parameterized according to
infinite parameters following our abstraction. In general, most solution strategies
transcribe the problem according to certain finite parameter values (supports) and
thus represent the problem in terms of these supports (e.g., using discrete time
points in dynamic optimization). This methodology can be generalized into the
following steps:
 - define supports for each infinite parameter if not already defined
 - expand any measures according to their underlying numerical representation
 - replace any remaining infinite variables with transcribed variables supported
   over each unique combination of the underlying parameter supports
 - replace any remaining infinite constraints with transcribed ones supported over
   all the unique support combinations stemming from the infinite parameters they
   depend on

TO BE CONTINUED WITH AN EXAMPLE

## TranscriptionOpt
`InfiniteOpt.TranscriptionOpt` is a submodule which principally implements
`TranscriptionModel`s and its related access/modification methods. Thus,
this section will detail what these are and how they work.

### TranscriptionModels
A `TranscriptionModel` is simply a `JuMP.Model` whose `ext` field contains
[`TranscriptionData`](@ref) which acts to map the transcribed model back to the
original infinite model (e.g., map the variables and constraints). Such models
are constructed via the [`TranscriptionModel`](@ref) constructors:
```jldoctest constrs
julia> model1 = TranscriptionModel() # make an empty model
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> model2 = TranscriptionModel(inf_model) # generate from an InfiniteModel
A JuMP Model
Minimization problem with:
Variables: 4
Objective function type: GenericAffExpr{Float64,VariableRef}
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`GenericQuadExpr{Float64,VariableRef}`-in-`MathOptInterface.LessThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.ZeroOne`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
Note that the all the normal `JuMP.Model` arguments can be used with both
constructors such as specifying the optimizer. The first constructor is what
`InfiniteOpt` uses to initialize the default `optimizer_model` attribute of
`InfiniteModel`s. The second constructor is used to build the optimizer model
when [`build_optimizer_model!`](@ref) is called directly or by
[`optimize!`](@ref JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory})).
Thus, second constructor serves as the principle tool for transcribing infinite
models as it encapsulates all of the methods to transcribe measures, variables,
and constraints.

### Queries
In this section we highlight a number of query methods that pertain
`TranscriptionModel`s and their mappings. First, if the `optimizer_model` of an
`InfiniteModel` is a `TranscriptionModel` it can be extracted via
[`transcription_model`](@ref):
```jldoctest constrs
julia> transcription_model(inf_model)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
Here we observe that such a model is currently empty and hasn't been populated
yet. Furthermore, we check that a `Model` is an `TranscriptionModel` via
[`is_transcription_model`](@ref):
```jldoctest constrs
julia> is_transcription_model(model2)
true

julia> is_transcription_model(Model())
false
```
We can also extract the raw [`TranscriptionData`](@ref) object from a
`TranscriptionModel` via [`transcription_data`](@ref).
```jldoctest constrs
julia> transcription_data(trans_model);
```

Next we can retrieve the `JuMP` variable(s) for a particular `InfiniteOpt`
variable via [`transcription_variable`](@ref). For finite variables, this will
be a one to one mapping, and for infinite variables a list of supported variables
will be returned in the order of the supports. Following the initial example in
the basic usage section, this is done:
```jldoctest transcribe
julia> transcription_variable(trans_model, g)
3-element Array{VariableRef,1}:
 g(support: 1)
 g(support: 2)
 g(support: 3)

julia> transcription_variable(trans_model, z)
z
```
Note that if the `TranscriptionModel` is stored as the current `optimizer_model`
then the first argument (specifying the `TranscriptionModel` can be omitted).
However, in this case the argument is required since `trans_model` was built
externally.

Similarly, the parameter supports corresponding to the transcription variables
(in the case of transcribed infinite variables) can be queried via
[`supports`](@ref InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteVariableRef)):
```jldoctest transcribe
julia> supports(trans_model, g)
3-element Array{Tuple{Float64},1}:
 (0.0,)
 (5.0,)
 (10.0,)
```
Again, the first argument can be dropped if this the `TranscriptionModel` of
interest is stored in the `optimizer_model` field of the `InfiniteModel` as is
the case when [`build_optimizer_model!`](@ref) or
[`optimize!`](@ref JuMP.optimize!(::InfiniteModel, ::Union{Nothing, JuMP.OptimizerFactory}))
is invoked.

Likewise, [`transcription_constraint`](@ref) and
[`supports`](@ref InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef))
can be used with constraints to find their transcribed equivalents in the
`JuMP` model and determine their supports. In the case of infinite constraints,
their parameter references can be determined
[`parameter_refs`](@ref InfiniteOpt.parameter_refs(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef))
just like infinite variables.

## Datatypes
```@index
Pages   = ["transcribe.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
TranscriptionData
```

## Methods
```@index
Pages   = ["transcribe.md"]
Modules = [InfiniteOpt, InfiniteOpt.TranscriptionOpt]
Order   = [:function]
```
```@docs
transcription_model
build_optimizer_model!(::InfiniteModel,::Val{:TransData})
TranscriptionModel()
TranscriptionModel(::InfiniteModel)
is_transcription_model
transcription_data
transcription_variable
InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteVariableRef)
InfiniteOpt.supports(::InfiniteOpt.InfiniteVariableRef)
transcription_constraint
InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.supports(::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.parameter_refs(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.parameter_refs(::InfiniteOpt.InfiniteConstraintRef)
```
