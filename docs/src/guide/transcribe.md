```@meta
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI", 
                  r" for all | ∀ ", r"d|∂"]
```

# [Model Transcription](@id transcription_docs)
A guide for transcribing infinite models in `InfiniteOpt`. See the respective 
[technical manual](@ref transcription_manual) for more details.

## Overview
All infinite models need to be transformed into a form that can be solved. A
common approach is direct transcription (e.g., discretize-then-optimize) where 
the domain of an infinite parameter is approximated by a collection of support 
points. This is the idea behind `TranscriptionOpt` (which implements 
[`TranscriptionBackend`](@ref)) which is InfiniteOpt's default transformation
backend as discussed in the [Infinite Models](@ref infinite_model_docs) section.
This page will detail the transcription models based in `InfiniteOpt.TranscriptionOpt`.

## Basic Usage
Most users will not need to employ the capabilities of `TranscriptionOpt` directly 
since they are employed implicitly with the call of 
[`optimize!`](@ref JuMP.optimize!(::InfiniteModel)) on an infinite model. This 
occurs since `TranscriptionBackend`s are the default backend. 

However, some users may wish to use `TranscriptionOpt` to extract a fully 
discretized/transcribed version of an infinite model that is conveniently output 
as a typical `JuMP` model and can then be treated as such. This is principally 
accomplished via [`build_transformation_backend!`](@ref). To illustrate how this is done, 
let's first define a basic infinite model with a simple support structure for the 
sake of example:
```jldoctest transcribe
julia> using InfiniteOpt

julia> inf_model = InfiniteModel();

julia> @infinite_parameter(inf_model, t in [0, 10], supports = [0, 5, 10])
t

julia> @variable(inf_model, y >= 0, Infinite(t))
y(t)

julia> @variable(inf_model, z, Bin)
z

julia> @objective(inf_model, Min, 2z + support_sum(y, t))
2 z + support_sum{t}[y(t)]

julia> @constraint(inf_model, initial, y(0) == 1)
initial : y(0) = 1

julia> @constraint(inf_model, constr, y^2 - z <= 42)
constr : y(t)² - z ≤ 42, ∀ t ∈ [0, 10]

julia> print(inf_model)
Min 2 z + support_sum{t}[y(t)]
Subject to
 y(t) ≥ 0, ∀ t ∈ [0, 10]
 z binary
 y(0) ≥ 0
 initial : y(0) = 1
 constr : y(t)² - z ≤ 42, ∀ t ∈ [0, 10]
```
Now we can make `JuMP` model containing the transcribed version of `inf_model` 
via [`build_transformation_backend!`](@ref) and then extract it via 
[`transformation_model`](@ref transformation_model(::InfiniteModel)):
```jldoctest transcribe
julia> build_transformation_backend!(inf_model)

julia> trans_model = transformation_model(inf_model)
A JuMP Model
Minimization problem with:
Variables: 4
Objective function type: AffExpr
`AffExpr`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`QuadExpr`-in-`MathOptInterface.LessThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.ZeroOne`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> print(trans_model)
Min 2 z + y(0.0) + y(5.0) + y(10.0)
Subject to
 initial : y(0.0) = 1
 constr[1] : y(0.0)² - z ≤ 42
 constr[2] : y(5.0)² - z ≤ 42
 constr[3] : y(10.0)² - z ≤ 42
 y(0.0) ≥ 0
 y(5.0) ≥ 0
 y(10.0) ≥ 0
 z binary
```

Thus, we have a transcribed `JuMP` model. To be precise, data on the mapping between 
the transcribed variables/constraints and their infinite counterparts is also generated 
as part of the `TranscriptionBackend` that `trans_model` is part of. Notice, that 
multiple finite variables have been introduced to discretize `y(t)` at supports 0, 5, and 10 which we can can also query via `supports`:
```jldoctest transcribe
julia> supports(y)
3-element Vector{Tuple}:
 (0.0,)
 (5.0,)
 (10.0,)
```
Also, notice how the constraints are transcribed in accordance with these supports 
except the initial condition which naturally is only invoked for the first support 
point. Furthermore, the transcription variable(s) of any variable associated with 
the infinite model can be determined via [`transformation_variable`](@ref):
```jldoctest transcribe
julia> transformation_variable(y)
3-element Vector{VariableRef}:
 y(0.0)
 y(5.0)
 y(10.0)

julia> transformation_variable(z)
z
```
Similarly, the transcription constraints associated with infinite model constraints 
can be queried via [`transformation_constraint`](@ref) and the associated supports 
and infinite parameters can be found via `supports` and `parameter_refs`:
```jldoctest transcribe
julia> transformation_constraint(initial)
initial : y(0.0) = 1

julia> transformation_constraint(constr)
3-element Vector{ConstraintRef}:
 constr[1] : y(0.0)² - z ≤ 42
 constr[2] : y(5.0)² - z ≤ 42
 constr[3] : y(10.0)² - z ≤ 42

julia> supports(constr)
3-element Vector{Tuple}:
 (0.0,)
 (5.0,)
 (10.0,)

julia> parameter_refs(constr)
(t,)
```
Note the parameter reference tuple corresponds to the support tuples. 

!!! note 
    Methods that query the transcription objects (e.g., `transcformation_variable`) 
    and the respective supports via `supports` also accept the keyword argument 
    `label` to specify which that transcription objects are desired in accordance 
    to the support labels that are inherited from and/or are equal to `label`. By 
    default, this will return any supports that are public (i.e., will hide anything 
    solely associated with internal supports). The full query response can always 
    be obtained via `label = All`.

Now we have a transcribed JuMP model that can be optimized via traditional 
JuMP compatible optimizers whose variables and constraints can be accessed using 
the methods mentioned above.

## Transcription Theory
A given infinite-dimensional optimization problem is parameterized according to 
infinite parameters following our abstraction. In general, most solution strategies 
transcribe the problem according to certain finite parameter values (supports) and 
thus represent the problem in terms of these supports (e.g., using discrete time 
points in dynamic optimization). This methodology can be generalized into the 
following steps:
 1. define supports for each infinite parameter if not already defined,
 2. add any additional support needed for derivative evaluation,
 3. expand any measures according to their underlying numerical representation 
    using transcribed infinite variables as appropriate,
 4. replace any remaining infinite variables/derivatives with transcribed 
    variables supported over each unique combination of the underlying parameter 
    supports,
 5. replace any remaining infinite constraints with transcribed ones supported over 
    all the unique support combinations stemming from the infinite parameters they 
    depend on,
 6. and add on the transcribed versions of the auxiliary derivative evaluation 
    equations. 

For example, let's consider a space-time optimization problem of the form:
```math
\begin{aligned}
	&&\min_{y(t), g(t, x)} &&& \int_0^{10} y^2(t) dt \\
	&&\text{s.t.} &&& y(0) = 1 \\
	&&&&& \int_{x \in [-1, 1]^2} \frac{\partial g(t, x)}{\partial t} dx = 42, && \forall t \in [0, 10] \\
  &&&&& 3g(t, x) + 2y^2(t) \leq 2, && \forall t \in T, \ x \in [-1, 1]^2. \\
\end{aligned}
```
Thus, we have an optimization problem whose decision space is infinite with 
respect to time ``t`` and position ``x``. Now let's transcribe it following the 
above steps. First, we need to specify the infinite parameter supports and for 
simplicity let's choose the following sparse sets:
 - ``t \in \{0, 5, 10\}``
 - ``x \in \{[-1, -1]^T, [-1, 1]^T, [1, -1]^T, [1, 1]^T\}``.
 To handle the derivative ``\frac{\partial g(t, x)}{\partial t}``, we'll use  
 backward finite difference, so no additional supports will need to be added.

Now we expand the two integrals (measures) via a finite approximation using only 
the above supports and term coefficients of 1 (note this is not numerically 
correct but is done for conciseness in example). Doing this, we obtain the 
form:
```math
\begin{aligned}
	&&\min_{y(t), g(t, x)} &&& y^2(0) + y^2(5) + y^2(10) \\
	&&\text{s.t.} &&& y(0) = 1 \\
  &&&&& g(0, x) = 0 \\
	&&&&& \frac{\partial g(t, [-1, -1])}{\partial t} + \frac{\partial g(t, [-1, 1])}{\partial t} + \frac{\partial g(t, [1, -1])}{\partial t} + \frac{\partial g(t, [1, 1])}{\partial t} = 42, && \forall t \in [0, 10] \\
  &&&&& 3g(t, x) + 2y^2(t) \leq 2, && \forall t \in T, \ x \in [-1, 1]^2. \\
\end{aligned}
```
Notice that the infinite variable ``y(t)`` in the objective measure has been 
replaced with finite transcribed variables ``y(0)``, ``y(5)``, ``y(10)``. Also, the 
infinite derivative ``\frac{\partial g(t, x)}{\partial t}`` was replaced with  
partially transcribed variables in the second constraint in accordance with the 
measure over the positional domain ``x``.

Now we need to transcribe the remaining infinite and semi-infinite variables 
with finite variables and duplicate the remaining infinite constraints accordingly. 
This means that the second constraint needs to be transcribed over the time domain 
and the third constraint needs to be transcribed for each unique combination 
of the time and position supports. Applying this transcription yields: 
```math
\begin{aligned}
	&&\min_{y(t), g(t, x)} &&& y^2(0) + y^2(5) + y^2(10) \\
	&&\text{s.t.} &&& y(0) = 1 \\
  &&&&& g(0, [-1, -1]) = 0 \\
  &&&&& g(0, [-1, 1]) = 0 \\
  &&&&& g(0, [1, -1]) = 0 \\
  &&&&& g(0, [1, 1]) = 0 \\
	&&&&& \frac{\partial g(0, [-1, -1])}{\partial t} + \frac{\partial g(0, [-1, 1])}{\partial t} + \frac{\partial g(0, [1, -1])}{\partial t} + \frac{\partial g(0, [1, 1])}{\partial t} = 42\\
  &&&&& \frac{\partial g(5, [-1, -1])}{\partial t} + \frac{\partial g(5, [-1, 1])}{\partial t} + \frac{\partial g(5, [1, -1])}{\partial t} + \frac{\partial g(5, [1, 1])}{\partial t} = 42\\
  &&&&& \frac{\partial g(10, [-1, -1])}{\partial t} + \frac{\partial g(10, [-1, 1])}{\partial t} + \frac{\partial g(10, [1, -1])}{\partial t} + \frac{\partial g(10, [1, 1])}{\partial t} = 42\\
  &&&&& 3g(0, [-1, -1]) + 2y^2(0) \leq 2 \\
  &&&&& 3g(0, [-1, 1]) + 2y^2(0) \leq 2 \\
  &&&&& \vdots \\
  &&&&& 3g(10, [1, 1]) + 2y^2(10) \leq 2.
\end{aligned}
```

Now that the variables and constraints are transcribed, all that remains is 
to add relations to define the behavior of the transcribed partial derivatives. 
We can accomplish this via backward finite difference which will just add one 
infinite equation in this case this we only have 2 supports in the time domain 
is then transcribed over the spatial domain to yield:
```math
\begin{aligned}
&&& g(5, [-1, -1]) = g(0, [-1, -1]) + 5\frac{\partial g(5, [-1, -1])}{\partial t} \\
&&& g(5, [-1, 1]) = g(0, [-1, 1]) + 5\frac{\partial g(5, [-1, 1])}{\partial t} \\
&&& g(5, [1, -1]) = g(0, [1, -1]) + 5\frac{\partial g(5, [1, -1])}{\partial t} \\
&&& g(5, [1, 1]) = g(0, [1, 1]) + 5\frac{\partial g(5, [1, 1])}{\partial t} \\
&&& g(10, [-1, -1]) = g(5, [-1, -1]) + 5\frac{\partial g(10, [-1, -1])}{\partial t} \\
&&& g(10, [-1, 1]) = g(5, [-1, 1]) + 5\frac{\partial g(10, [-1, 1])}{\partial t} \\
&&& g(10, [1, -1]) = g(5, [1, -1]) + 5\frac{\partial g(10, [1, -1])}{\partial t} \\
&&& g(10, [1, 1]) = g(5, [1, 1]) + 5\frac{\partial g(10, [1, 1])}{\partial t}
\end{aligned}
```

Now the problem is fully transcribed (discretized) and can be solved as a 
standard optimization problem. Note that with realistic measure evaluation 
schemes more supports might be added to the support sets and these will need to 
be incorporated when transcribing variables and constraints.

It is easy to imagine how the above procedure can get quite involved to do manually, 
but this is precisely what `InfiniteOpt` automates behind the scenes. Let's 
highlight this by repeating the same example using `InfiniteOpt` (again using 
the incorrect simple representation for the integrals for conciseness).
```jldoctest trans_example
using InfiniteOpt

# Initialize model
inf_model = InfiniteModel()

# Define parameters and supports
@infinite_parameter(inf_model, t in [0, 10], supports = [0, 5, 10])
@infinite_parameter(inf_model, x[1:2] in [-1, 1], supports = [-1, 1], independent = true)

# Define variables
@variable(inf_model, y, Infinite(t))
@variable(inf_model, g, Infinite(t, x))

# Set the objective (using support_sum for the integral given our simple example)
# Note: In real problems integral should be used
@objective(inf_model, Min, support_sum(y^2, t))

# Define the constraints
@constraint(inf_model, y(0) == 1)
@constraint(inf_model, g(0, x) == 0)
@constraint(inf_model, support_sum(deriv(g, t), x) == 42) # support_sum for simplicity
@constraint(inf_model, 3g + y^2 <= 2)

# Print the infinite model
print(inf_model)

# output
Min support_sum{t}[y(t)²]
Subject to
 y(0) = 1
 g(0, [x[1], x[2]]) = 0, ∀ x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
 support_sum{x}[∂/∂t[g(t, x)]] = 42, ∀ t ∈ [0, 10]
 y(t)² + 3 g(t, x) ≤ 2, ∀ t ∈ [0, 10], x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
```
Thus, we obtain the infinite problem in `InfiniteOpt`. As previously noted, 
transcription would be handled automatically behind the scenes when the model is 
optimized. However, we can directly extract the transcribed version by building a 
`TranscriptionBackend`:
```jldoctest trans_example
julia> build_transformation_backend!(inf_model)

julia> trans_model = transformation_model(inf_model);

julia> print(trans_model)
Min y(0.0)² + y(5.0)² + y(10.0)²
Subject to
 y(0.0) = 1
 g(0.0, [-1.0, -1.0]) = 0
 g(0.0, [1.0, -1.0]) = 0
 g(0.0, [-1.0, 1.0]) = 0
 g(0.0, [1.0, 1.0]) = 0
 d/dt[g(t, x)](0.0, [-1.0, -1.0]) + d/dt[g(t, x)](0.0, [1.0, -1.0]) + d/dt[g(t, x)](0.0, [-1.0, 1.0]) + d/dt[g(t, x)](0.0, [1.0, 1.0]) = 42
 d/dt[g(t, x)](5.0, [-1.0, -1.0]) + d/dt[g(t, x)](5.0, [1.0, -1.0]) + d/dt[g(t, x)](5.0, [-1.0, 1.0]) + d/dt[g(t, x)](5.0, [1.0, 1.0]) = 42
 d/dt[g(t, x)](10.0, [-1.0, -1.0]) + d/dt[g(t, x)](10.0, [1.0, -1.0]) + d/dt[g(t, x)](10.0, [-1.0, 1.0]) + d/dt[g(t, x)](10.0, [1.0, 1.0]) = 42
 g(0.0, [-1.0, -1.0]) - g(5.0, [-1.0, -1.0]) + 5 d/dt[g(t, x)](5.0, [-1.0, -1.0]) = 0
 g(5.0, [-1.0, -1.0]) - g(10.0, [-1.0, -1.0]) + 5 d/dt[g(t, x)](10.0, [-1.0, -1.0]) = 0
 g(0.0, [1.0, -1.0]) - g(5.0, [1.0, -1.0]) + 5 d/dt[g(t, x)](5.0, [1.0, -1.0]) = 0
 g(5.0, [1.0, -1.0]) - g(10.0, [1.0, -1.0]) + 5 d/dt[g(t, x)](10.0, [1.0, -1.0]) = 0
 g(0.0, [-1.0, 1.0]) - g(5.0, [-1.0, 1.0]) + 5 d/dt[g(t, x)](5.0, [-1.0, 1.0]) = 0
 g(5.0, [-1.0, 1.0]) - g(10.0, [-1.0, 1.0]) + 5 d/dt[g(t, x)](10.0, [-1.0, 1.0]) = 0
 g(0.0, [1.0, 1.0]) - g(5.0, [1.0, 1.0]) + 5 d/dt[g(t, x)](5.0, [1.0, 1.0]) = 0
 g(5.0, [1.0, 1.0]) - g(10.0, [1.0, 1.0]) + 5 d/dt[g(t, x)](10.0, [1.0, 1.0]) = 0
 y(0.0)² + 3 g(0.0, [-1.0, -1.0]) ≤ 2
 y(5.0)² + 3 g(5.0, [-1.0, -1.0]) ≤ 2
 y(10.0)² + 3 g(10.0, [-1.0, -1.0]) ≤ 2
 y(0.0)² + 3 g(0.0, [1.0, -1.0]) ≤ 2
 y(5.0)² + 3 g(5.0, [1.0, -1.0]) ≤ 2
 y(10.0)² + 3 g(10.0, [1.0, -1.0]) ≤ 2
 y(0.0)² + 3 g(0.0, [-1.0, 1.0]) ≤ 2
 y(5.0)² + 3 g(5.0, [-1.0, 1.0]) ≤ 2
 y(10.0)² + 3 g(10.0, [-1.0, 1.0]) ≤ 2
 y(0.0)² + 3 g(0.0, [1.0, 1.0]) ≤ 2
 y(5.0)² + 3 g(5.0, [1.0, 1.0]) ≤ 2
 y(10.0)² + 3 g(10.0, [1.0, 1.0]) ≤ 2
```
This precisely matches what we found analytically. Note that the unique support 
combinations are determined automatically.

## TranscriptionOpt
`InfiniteOpt.TranscriptionOpt` is a sub-module which principally implements 
[`TranscriptionBackend`](@ref)s and its related access/modification methods. Thus, 
this section will detail what these are and how they work.

### TranscriptionBackends
A `TranscriptionBackend` is simply a [`JuMPBackend`](@ref) that uses the 
[`Transcription`](@ref InfiniteOpt.TranscriptionOpt.Transcription) `AbstractJuMPTag`
and [`TranscriptionData`](@ref InfiniteOpt.TranscriptionOpt.TranscriptionData) which 
acts to map the transcribed model back to the original infinite model (e.g., map the 
variables and constraints). Such models are constructed via a default version of 
[`build_transformation_backend!`](@ref InfiniteOpt.build_transformation_backend!(::InfiniteOpt.InfiniteModel,::TranscriptionBackend)) 
which wraps [`build_transcription_backend!`](@ref InfiniteOpt.TranscriptionOpt.build_transcription_backend!):
```jldoctest transcribe
julia> backend1 = TranscriptionBackend() # make an empty backend
A TranscriptionBackend that uses a
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> build_transformation_backend!(inf_model); 

julia> backend2 = transformation_backend(inf_model) # generate from an InfiniteModel
A TranscriptionBackend that uses a
A JuMP Model
Minimization problem with:
Variables: 4
Objective function type: AffExpr
`AffExpr`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`QuadExpr`-in-`MathOptInterface.LessThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 3 constraints
`VariableRef`-in-`MathOptInterface.ZeroOne`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
The call to `build_transformation_backend!` is the backbone 
behind infinite model transformation and is what encapsulates all the methods to 
transcribe measures, variables, derivatives, and constraints. This is also the 
method that enables the use of [`optimize!`](@ref JuMP.optimize!(::InfiniteModel)).

### Queries
In this section we highlight a number of query methods that pertain to 
`TranscriptionBackend`s and their mappings. First, we can retrieve the JuMP `Model`
via [`transformation_model`](@ref):
```jldoctest transcribe; setup = :(empty!(inf_model.backend))
julia> transformation_model(inf_model)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
Here we observe that such a model is currently empty and hasn't been populated 
yet.

Next we can retrieve the `JuMP` variable(s) for a particular `InfiniteOpt` 
variable via [`transformation_variable`](@ref). For finite variables, this will 
be a one to one mapping, and for infinite variables an array will be returned that corresponds
to the underlying supports. Following the initial example in 
the basic usage section, this is done:
```jldoctest transcribe
julia> build_transformation_backend!(inf_model); backend = transformation_backend(inf_model);

julia> transformation_variable(y, backend)
3-element Vector{VariableRef}:
 y(0.0)
 y(5.0)
 y(10.0)

julia> transformation_variable(z, backend)
z
```
Note that if the `TranscriptionBackend` is the current backend then, 
then the 2nd argument can be omitted.

Similarly, the parameter supports corresponding to the transcription variables 
(in the case of transcribed infinite variables) can be queried via 
[`supports`](@ref):
```jldoctest transcribe
julia> supports(y)
3-element Vector{Tuple}:
 (0.0,)
 (5.0,)
 (10.0,)
```

!!! note 
    Note that like `supports`, the `transformation_[obj]` methods also employ the 
    `label::Type{AbstractSupportLabel} = PublicLabel` keyword argument that by 
    default will return variables/expressions/constraints associated with public 
    supports. The full set (e.g., ones corresponding to internal collocation nodes) 
    is obtained via `label = All`. 

Likewise, [`transformation_constraint`](@ref transformation_constraint(::InfOptConstraintRef)) and 
[`supports`](@ref supports(::InfOptConstraintRef)) can be used with constraints to find their transcribed equivalents in the `JuMP` model and determine their supports.

We can also do this with measures and expressions:
```jldoctest transcribe
julia> meas = support_sum(y^2, t)
support_sum{t}[y(t)²]

julia> build_transformation_backend!(inf_model)

julia> transformation_variable(meas)
y(0.0)² + y(5.0)² + y(10.0)²

julia> supports(meas)
()

julia> transformation_expression(y^2 + z - 42)
3-element Vector{AbstractJuMPScalar}:
 y(0.0)² + z - 42
 y(5.0)² + z - 42
 y(10.0)² + z - 42

julia> supports(y^2 + z - 42)
3-element Vector{Tuple}:
 (0.0,)
 (5.0,)
 (10.0,)

julia> parameter_refs(y^2 + z - 42)
(t,)
```
