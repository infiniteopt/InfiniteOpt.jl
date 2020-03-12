# Measures
A guide and manual for the definition and use of measures in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview

Measures are objects that capture the integration of an expression with respect
to parameters, which is a distinct feature of optimization problems with
infinite decision spaces. In dynamic optimization measures can represent integral
terms such as the total cost over time, and in stochastic optimization measures
can represent integrals over the uncertain parameters, such as expectation. In
`InfiniteOpt`, measures are evaluated by some discretization scheme, which
evaluates the expression at a set of points over the parameter space and
approximates the measures based on the expression values at these points.

## [Basic Usage] (@id measure_basic_usage)

First, we consider a dynamic optimization problem with the time parameter `t`
from 0 to 10. We also consider a state variable `y(t)` and a control variable
`u(t)` that are parameterized by `t`:
```jldoctest meas_basic; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(seed = true))
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_variable(model, y(t))
y(t)

julia> @infinite_variable(model, u(t))
u(t)
```

Now suppose we want to evaluate the integral ``\int_{2}^{8}y(t)^2 + u(t)^2 dt``.
We can construct a measure to represent this integral using the
[`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})) function
```jldoctest meas_basic
julia> mref1 = measure(y^2 + u^2, t, 2, 8)
measure(y(t)² + u(t)²)
```
The four positional arguments of [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})) are the integrand expression, the the parameter of integration, the lower bound, and the upper
bound. If the lower and upper bounds are not specified, then the integration will
be over the entire domain, which is [0, 10] in this case. In addition, if the
parameter of integration is not specified, the measure will search for the
parameter in the expression and choose that as the measure of integration.
Specifying the parameter of integration is required if the expression contains
multiple groups of infinite parameters.

The measure function uses Monte Carlo (MC) sampling method as the default
discretization scheme. However, for integration over univariate parameter,
the user can also use quadrature methods by setting the keyword argument
`eval_method` as `Quad`:
```jldoctest meas_basic
julia> mref2 = measure(y^2 + u^2, eval_method = Quad)
measure(y(t)² + u(t)²)
```

The measure function also for specifying the number of points for the
discretization scheme using the keyword argument `num_supports`. The default
value of `num_supports` is 50.
```jldoctest meas_basic
julia> mref3 = measure(y^2 + u^2, num_supports = 10)
measure(y(t)² + u(t)²)
```

Every time a measure function is called for an infinite parameter, the points
at which the integral is evaluated will be added to the list of supports of the
relevant infinite parameter. The use can choose to evaluate a new integral over
points that already exist in the support list by setting the keyword argument
`use_existing_supports` as `true`:
```jldoctest meas_basic
julia> mref4 = measure(y^2 + u^2, use_existing_supports = true)
measure(y(t)² + u(t)²)
```
Note that by setting `use_existing_supports = true`, the measure function will
ignore the values of `eval_method` and `num_supports`.

Once a measure is created, the evaluation of that measure is stored in a
measure data object. Users can query the measure data object using the
[`measure_data`](@ref) function as follows
```jldoctest meas_basic
julia> measure_data(mref3)
DiscreteMeasureData(t, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.96407, 8.73581, 6.54922, 5.86712, 6.32696, 7.35004, 5.48342, 4.76883, 6.02375, 7.91346], "measure", InfiniteOpt._w)
```
For more details on the measure data object, refer to [Measure Data Generation](@ref).

Similarly, one can also query the integrand function using [`measure_function`](@ref):
```jldoctest meas_basic
julia> measure_function(mref3)
y(t)² + u(t)²
```

Now suppose we want to create multiple measures for the same model that share
the same keyword argument values. We don't have to input the keyword argument
values every time we construct a new measure. Instead, we can modify the
default values of measure keyword arguments for that model, and construct
measures using the new default values. To do that, we use the functions
[`set_measure_defaults`](@ref) and [`measure_defaults`](@ref). For example,
suppose we want to change the default number of supports to 10, and for some
reason we want to add a new keyword argument called `a` to the measure function,
with value of 1. Then we can do the following:
```jldoctest meas_basic
julia> measure_defaults(model)
Dict{Symbol,Any} with 7 entries:
  :num_supports          => 50
  :call_from_expect      => false
  :eval_method           => :Sampling
  :name                  => "measure"
  :check_method          => true
  :weight_func           => _w
  :use_existing_supports => false

julia> set_measure_defaults(model, num_supports = 10, a = 1)

julia> measure_defaults(model)
Dict{Symbol,Any} with 8 entries:
  :a                     => 1
  :num_supports          => 10
  :call_from_expect      => false
  :eval_method           => :Sampling
  :name                  => "measure"
  :check_method          => true
  :weight_func           => _w
  :use_existing_supports => false
```
[`measure_defaults`](@ref) returns the current default keyword argument
values for measures for the model, and [`set_measure_defaults`](@ref) set new
keyword argument values. The first [`measure_defaults`](@ref) call above
shows the default values at model initialization. Note that
[`set_measure_defaults`](@ref) can not only modify values of existing keyword
arguments, but also add new keyword arguments with a specified value. This is
very useful if users want to extend the measure functions with their custom
discretization/evaluation schemes that take additional arguments.
See the extension page for more details.

Now that we change the default number of supports to 10, all measures
created later will have 10 supports by default. The users can still overwrite
the default for specific measures as follows, shown as follows:
```jldoctest meas_basic; setup = :(delete!(measure_defaults(model), :a))
julia> mref1 = measure(y^2)
measure(y(t)²)

julia> length(measure_data(mref1).supports)
10

julia> mref2 = measure(y^2, num_supports = 20)
measure(y(t)²)

julia> length(measure_data(mref2).supports)
20
```

Now we can add integrals to the constraints and objective functions in our
model using these measures. For more detailed information, please review the
information below.

## Theoretical Abstraction
In `InfiniteOpt`, measures represent integrals of the form
```math
\int_{\tau \in \mathcal{T}} f(\tau)w(\tau) d\tau
```
where ``\tau`` is a
(possibly multivariate) infinite parameter, ``f(\tau)`` is an expression
parameterized by ``\tau``, ``w(\tau)`` is a weight function, and ``\mathcal{T}``
is a subset of the domain of ``\tau``. The measures approximate the integrals
by taking a discretization scheme
```math
\int_{\tau \in \mathcal{T}} f(\tau)w(\tau) d\tau \approx \sum_{i=1}^N \alpha_i f(\tau_i) w(\tau_i)
```
where ``\tau_i`` are the grid points where the expression ``f(\tau)`` is
evaluated, and ``N`` is the total number of points taken.

## Measure Data Generation
The most general form of [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)) takes two arguments: the expression to integrate and
a measure data object that contains the details of the discretization scheme.
Measure data objects can be constructed using [`DiscreteMeasureData`](@ref DiscreteMeasureData(::ParameterRef, ::Vector{<:Number}, ::Vector{<:Number})),
where the parameter of integration, the coefficients ``\alpha_i``, and the
support points need to be defined explicitly. For example, if we want to
evaluate a function at each integer time point between 0 and 10, we
can construct the following measure data object to record this discretization
scheme:
```jldoctest meas_basic
julia> md_t = DiscreteMeasureData(t, ones(10), [i for i in 1:10])
DiscreteMeasureData(t, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], "measure", InfiniteOpt._w)
```
The arguments of [`DiscreteMeasureData`](@ref DiscreteMeasureData(::ParameterRef, ::Vector{<:Number}, ::Vector{<:Number})) are parameter, coefficients, and
supports. Users can define a name for the measure data object using the `name`
keyword argument. The default weight function is ``w(\tau) = 1`` for
any ``\tau``, which can be overwritten by the keyword argument `weight_function`.
The `weight_function` should take a function that returns a number for any
value that is well defined for the integrated infinite parameter. The data type
is [`DiscreteMeasureData`](@ref), which is a subtype of the abstract data type
[`AbstractMeasureData`](@ref).

With [`DiscreteMeasureData`](@ref), a measure can be generated in a custom and
quick manner. For example, using the measure data above, we can define a measure
for ``y^2`` as follows:
```jldoctest meas_basic
julia> mref = measure(y^2, md_t)
measure(y(t)²)
```
In the same way, we can define measure data for multivariate infinite parameter.
For example, we can define a discretization scheme for a 2D position parameter
``x \in [0, 1] \times [0, 1]`` as follows:
```jldoctest meas_basic
julia> @infinite_parameter(model, x[1:2] in [0, 1])
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> md_x = DiscreteMeasureData(x, 0.25 * ones(4), [[0.25, 0.25], [0.25, 0.75], [0.75, 0.25], [0.75, 0.75]])
MultiDiscreteMeasureData(  [2]  =  x[2]
  [1]  =  x[1], [0.25, 0.25, 0.25, 0.25], JuMP.Containers.SparseAxisArray{Float64,1,Tuple{Int64}}[  [2]  =  0.25
  [1]  =  0.25,   [2]  =  0.75
  [1]  =  0.25,   [2]  =  0.25
  [1]  =  0.75,   [2]  =  0.75
  [1]  =  0.75], "measure", InfiniteOpt._w)
```
where `md_x` cuts the domain into four 0.5-by-0.5 squares, and evaluates the
integrand on the center of these squares.

Note that for multivariate parameters, each support should be an `AbstractArray`
that records the value at each dimension, and the measure data type is
[`MultiDiscreteMeasureData`](@ref), again a subtype of
[`AbstractMeasureData`](@ref) that's specifically for multivariate parameters.

The [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})) function called
in [Basic Usage](@ref measure_basic_usage) above is a dispatch of the basic [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)) function which does
not require explicit construction of the measure data object. Instead, the
function constructs the appropriate measure data object according to the values
of the positional and keyword arguments.

## Evaluation Methods
If the users call [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})) without giving the measure data object but giving enough
information, the function calls [`generate_measure_data`](@ref) under the hood
to construct the measure data object. [`generate_measure_data`](@ref) takes as
positional arguments the integrated parameter, number of supports, lower bound
and upper bound, and returns a measure data object of type
[`AbstractMeasureData`](@ref).

[`generate_measure_data`](@ref) calls integral evaluation methods that are
encoded in this package, which is specified through the `eval_method` keyword
argument in measure function calls. `eval_method` takes a function that returns
a tuple containing the supports and coefficients, which is will be used by
[`generate_measure_data`](@ref) to construct the measure data object.
A table of available `eval_method` options in our package is listed below.
Each method is limited on the dimension of parameter and/or the type of set
that it can apply for. `eval_method` can be overwritten by self-defined
functions as long as they preserve the same input and output argument types.

| `eval_method` | Uni/multi-variate? | Type of set |
| --- | --- | --- |
| [`mc_sampling`](@ref mc_sampling(::Number, ::Number, ::Int)) | Both | [`IntervalSet`](@ref) |
| [`mc_sampling`](@ref mc_sampling(::Distributions.UnivariateDistribution, ::InfiniteOpt.ParameterRef, ::Int)) | Both | [`DistributionSet`](@ref) |
| [`gauss_legendre`](@ref) | Univariate | Finite [`IntervalSet`](@ref) |
| [`gauss_laguerre`](@ref) | Univariate | Semi-infinite [`IntervalSet`](@ref) |
| [`gauss_hermite`](@ref) | Univariate | Infinite [`IntervalSet`](@ref) |

The default evaluation method is MC sampling.
If the integrated parameter is defined on a finite interval set
(using [`IntervalSet`](@ref)), the MC sampling method will draw samples from
a uniform distribution over that interval. For (semi-)infinite domains, the
MC sampling method will map the domain into a finite domain and do uniform
sampling on that finite domain. The mapping is based on a change of variable
that preserves the same integral value as follows, which is encoded in
the [`infinite_transform`](@ref) function:
```math
\int_{-\infty}^{\infty} f(x) dx = \int_{-1}^{1}f\left(\frac{t}{1-t^2}\right)\frac{1+t^2}{(1-t^2)^2}dt \\
\int_{a}^{\infty} f(x) dx = \int_{0}^{1}f\left(a + \frac{t}{1-t}\right)\frac{dt}{(1-t)^2} \\
\int_{-\infty}^{a} f(x) dx = \int_{0}^{1}f\left(a - \frac{1-t}{t}\right)\frac{dt}{t^2} \\
```
If the integrated infinite parameter is a random variable defined in a
[`DistributionSet`](@ref), then MC sampling will draw samples from the
underlying distribution.

For univariate infinite parameters, this package also supports basic Gaussian
quadrature schemes. Specifically, [`gauss_legendre`](@ref) generates supports
and coefficients based on Gauss-Legendre quadrature method for parameters in
finite interval. [`gauss_laguerre`](@ref) generates supports and coefficients
based on Gauss-Laguerre quadrature method for parameters in semi-infinite
interval. [`gauss_hermite`](@ref) generates supports and coefficients based
on Gauss-Hermite quadrature method for parameters infinite interval. For an
introduction to the mathematics of Gaussian quadrature methods, see
[here](https://en.wikipedia.org/wiki/Gaussian_quadrature).

## Expansion
In a model, each measure records the integrand expression and an evaluation
scheme that details the discretization scheme to approximate the integral.
The model will not expand the measures until the transcription stage, at which
a `JuMP.AbstractJuMPScalar` is created for each measure to represent how
the measure is modeled in a transcription model based on the stored
discretization scheme (see [Model Transcription](@ref transcription_docs) for
details on transcription). Additional point variables will be created in the
expansion process if the measure is evaluated at infinite parameter points that
do not have corresponding point variables yet.

Sometimes for extension purposes, one might want to expand a specific measure
before reaching the transcription stage. Alternatively, one might want to use
custom reformulation instead of the transcription encoded in this package, in
which expanding measures will also be useful. This can be done using the [`expand`](@ref)
function, which takes a [`MeasureRef`](@ref) object and returns a `JuMP.AbstractJuMPScalar`
based on the [`AbstractMeasureData`](@ref). For example, suppose we want to
integrate ``y^2`` in ``t``, with two supports ``t = 2.5`` and ``t = 7.5``.
We can set up and expand this measure as follows:
```jldoctest meas_basic
julia> tdata = DiscreteMeasureData(t, [5, 5], [2.5, 7.5])
DiscreteMeasureData(t, [5, 5], [2.5, 7.5], "measure", InfiniteOpt._w)

julia> mref3 = measure(y^2, tdata)
measure(y(t)²)

julia> expanded_measure = expand(mref3)
5 y(2.5)² + 5 y(7.5)²

julia> typeof(expanded_measure)
GenericQuadExpr{Float64,GeneralVariableRef}
```
In the expand call, two point variables, `y(2.5)` and `y(7.5)`, are created
because they are not defined in the model before the expand call. One can use
the [`expand_all_measures!`](@ref) function to expand all measures in a model,
which simply applies the [`expand`](@ref) to all measures stored in the model.

## Reduced Infinite Variables
Expanding partial integrals will introduce reduced infinite variables to the
model. To see what this means, suppose we have an infinite variable that is
parameterized by multiple infinite parameters defined as follows:
```jldoctest meas_basic
julia> @infinite_variable(model, T(x, t))
T(x, t)
```
Now say we want to integrate `T` over `t`. We can define a measure for the
integral similar to how we have defined other measures:
```jldoctest meas_basic
julia> mref4 = measure(T, tdata)
measure(T(x, t))
```
Now if we expand this measure, the measure data object `tdata` records the
supports for `t`, but no supports for `x` because `T` is not integrated over
`x` in this measure. Therefore, point variables cannot be defined in the
measure expansion.

Instead of point variables, each new variable in the measure expansion will be
represented using reduced infinite variables. Reduced infinite variables are
"reduced" from their original infinite variables in that they are parameterized
by less infinite parameters. In the example above, in the expansion each
reduced infinite variable for `T` should only be parameterized by `x` since
the value of `t` is fixed. The expanded measure now looks like this:
```jldoctest meas_basic
julia> expanded_measure = expand(mref4)
5 T(x, 2.5) + 5 T(x, 7.5)
```
where the expanded measure is a `JuMP.GenericAffExpr` that takes
[`ReducedInfiniteVariableRef`](@ref) in its terms. [`ReducedInfiniteVariableRef`](@ref)
refers to the information of the reduced infinite variable stored in its model,
in the type of [`ReducedInfiniteInfo`](@ref). Each [`ReducedInfiniteInfo`](@ref)
records a reference for its original infinite variable, and the value of the
fixed infinite parameter. One can query this information using
[`infinite_variable_ref`](@ref) and [`eval_supports`](@ref) function as follows:
```jldoctest meas_basic
julia> T1 = first(keys(expanded_measure.terms))
T(x, 2.5)

julia> typeof(T1)
ReducedInfiniteVariableRef

julia> infinite_variable_ref(T1)
T(x, t)

julia> eval_supports(T1)
Dict{Int64,Union{Number, SparseAxisArray{#s57,N,K} where K<:Tuple{Vararg{Any,N}} where N where #s57<:Number}} with 1 entry:
  2 => 2.5
```
All the `JuMP` functions extended for infinite variables are also extended for
reduced infinite variables, e.g. [`JuMP.lower_bound`](@ref JuMP.lower_bound(::ReducedInfiniteVariableRef)).

## Datatypes
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractMeasureData
DiscreteMeasureData
MultiDiscreteMeasureData
Measure
MeasureRef
AbstractReducedInfo
ReducedInfiniteInfo
ReducedInfiniteVariableRef
```

## Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:function]
```
```@docs
expect
support_sum
measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing})
measure_defaults
set_measure_defaults
DiscreteMeasureData(::ParameterRef, ::Vector{<:Number}, ::Vector{<:Number})
DiscreteMeasureData(::AbstractArray{<:ParameterRef}, ::Vector{<:Number}, ::Vector{<:AbstractArray})
measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)
add_measure
measure_function
measure_data
used_by_constraint(::MeasureRef)
used_by_measure(::MeasureRef)
used_by_objective(::MeasureRef)
is_used(::MeasureRef)
JuMP.is_valid(::InfiniteModel, ::MeasureRef)
JuMP.delete(::InfiniteModel, ::MeasureRef)
JuMP.name(::MeasureRef)
JuMP.set_name(::MeasureRef, ::String)
expand
expand_all_measures!
measure_name
parameter_refs(::AbstractMeasureData)
supports(::AbstractMeasureData)
measure_data_in_hold_bounds
expand_measure
make_point_variable_ref
make_reduced_variable_ref(::InfiniteVariableRef, ::Int, ::Union{Number,JuMP.Containers.SparseAxisArray{<:Number}})
make_reduced_variable_ref(::InfiniteVariableRef, ::Dict)
infinite_variable_ref(::ReducedInfiniteVariableRef)
eval_supports(::ReducedInfiniteVariableRef)
parameter_refs(::ReducedInfiniteVariableRef)
JuMP.name(::ReducedInfiniteVariableRef)
JuMP.has_lower_bound(::ReducedInfiniteVariableRef)
JuMP.lower_bound(::ReducedInfiniteVariableRef)
JuMP.LowerBoundRef(::ReducedInfiniteVariableRef)
JuMP.has_upper_bound(::ReducedInfiniteVariableRef)
JuMP.upper_bound(::ReducedInfiniteVariableRef)
JuMP.UpperBoundRef(::ReducedInfiniteVariableRef)
JuMP.is_fixed(::ReducedInfiniteVariableRef)
JuMP.fix_value(::ReducedInfiniteVariableRef)
JuMP.FixRef(::ReducedInfiniteVariableRef)
JuMP.start_value(::ReducedInfiniteVariableRef)
JuMP.is_binary(::ReducedInfiniteVariableRef)
JuMP.BinaryRef(::ReducedInfiniteVariableRef)
JuMP.is_integer(::ReducedInfiniteVariableRef)
JuMP.IntegerRef(::ReducedInfiniteVariableRef)
used_by_constraint(::ReducedInfiniteVariableRef)
used_by_measure(::ReducedInfiniteVariableRef)
JuMP.is_valid(::InfiniteModel, ::ReducedInfiniteVariableRef)
JuMP.delete(::InfiniteModel, ::ReducedInfiniteVariableRef)
```

## MeasureEvalMethods Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt.MeasureEvalMethods]
Order   = [:function]
```

```@docs
InfiniteOpt.MeasureEvalMethods.generate_measure_data
InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs
InfiniteOpt.MeasureEvalMethods.eval_method_registry
InfiniteOpt.MeasureEvalMethods.register_eval_method
InfiniteOpt.MeasureEvalMethods.mc_sampling(::Number, ::Number, ::Int)
InfiniteOpt.MeasureEvalMethods.mc_sampling(::Distributions.UnivariateDistribution, ::InfiniteOpt.ParameterRef, ::Int)
InfiniteOpt.MeasureEvalMethods.gauss_legendre
InfiniteOpt.MeasureEvalMethods.gauss_hermite
InfiniteOpt.MeasureEvalMethods.gauss_laguerre
InfiniteOpt.MeasureEvalMethods.infinite_transform
```
