# Extensions
Here we provide guidance to various ways `InfinitOpt` can be extended.

## Overview
Extendibility is one of the core ideas of `InfinitOpt` so that it can serve as a
convenient tool for those developing and implementing advanced techniques for
infinite dimensional optimization problems. Thus, `InfiniteOpt` is developed in
a modular manner to readily accommodate user-defined functionality and/or to
serve as useful base in writing a `JuMP` extension. Admittedly, this modularity
might be imperfect and comments/suggestions are welcomed to help us improve this.

## Infinite Sets
Infinite sets are used to characterize the behavior of infinite parameters and
used to govern the behavior of supports in `InfiniteOpt`. Here we walk through
how user-defined sets can be added to various degrees of functionality. A
template is provided in `./InfiniteOpt/test/extensions/infinite_set.jl`. The
extension steps employed are:
1. Define the new `struct` infinite set type (only thing required as bare minimum)
2. Extend [`InfiniteOpt.supports_in_set`](@ref) (enables error checking of supports)
3. Extend [`InfiniteOpt.generate_support_values`](@ref) (enables support generation via `num_supports` keyword arguments)
4. Extend [`InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs`](@ref) and register the new set type via [`register_eval_method`](@ref) (enables measure evaluation methods)
5. If a lower bound and upper bound can be reported, extend `JuMP` lower bound and upper bound methods

As an example, let's create a disjoint interval set as an infinite set type.
This corresponds to the set ``[lb_1, ub_1] \cup [lb_2, ub_2]`` where
``ub_1 \leq lb_2``. First, we need to create the `DataType` with inheritance from [`AbstractInfiniteSet`](@ref):
```jldoctest set_ext; output = false
using InfiniteOpt, JuMP

struct DisjointSet <: InfiniteOpt.AbstractInfiniteSet
    lb1::Float64
    ub1::Float64
    lb2::Float64
    ub2::Float64
    # constructor
    function DisjointSet(lb1::Number, ub1::Number, lb2::Number, ub2::Number)
        if lb1 > ub1 || lb2 > ub2 || ub1 > lb2
            error("Invalid bounds")
        end
        return new(convert(Float64, lb1), convert(Float64, ub1),
                   convert(Float64, lb2), convert(Float64, ub2))
    end
end

# output


```
Notice that we also define the constructor function to error check and convert as
needed (this is recommended, but not required). For basic functionality this is
all we have to do to add a set in `InfiniteOpt`.

We can now define infinite parameters using this set via
[`@infinite_parameter`](@ref) both anonymously and explicitly:
```jldoctest set_ext
julia> model = InfiniteModel();

julia> t = @infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), base_name = "t")
t

julia> @infinite_parameter(model, t in DisjointSet(0, 1, 3, 4))
t
```
Once defined (without further extension), these parameters can be used as normal
with the following limitations:
- Supports must be specified manually (`num_supports` is not enabled)
- Supports will not be checked if they are in the domain of the infinite set
- The [`DiscreteMeasureData`](@ref) must be provided explicitly to evaluate measures
- Set bounds cannot be queried.
However, all of these limitations can be eliminated by extending a few functions
as outlined above.

To enable support domain checking which is useful to avoid strange bugs, we will
extend [`InfiniteOpt.supports_in_set`](@ref). This returns a `Bool` to indicate
if a vector of supports are in the set's domain:
```jldoctest set_ext; output = false
function InfiniteOpt.supports_in_set(supports::Union{Number, Vector{<:Number}},
                                     set::DisjointSet)::Bool
    return all((set.lb1 .<= supports .<= set.ub1) .| (set.lb2 .<= supports .<= set.ub2))
end

# output


```
Now the checks are enabled so, the following would yield an error because the
support is not in the set domain:
```jldoctest set_ext
julia> @infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), supports = 2)
ERROR: In `@infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), supports = 2)`: Supports violate the set domain bounds.
```

To enable automatic support generation via the `num_supports` keyword and with
functions such as [`fill_in_supports!`](@ref), we will extend
[`InfiniteOpt.generate_support_values`](@ref):
```jldoctest set_ext; output = false
function InfiniteOpt.generate_support_values(set::DisjointSet;
                                             num_supports::Int = 50,
                                             sig_fig::Int = 5)::Array
    length_ratio = (set.ub1 - set.lb1) / (set.ub1 - set.lb1 + set.ub2 - set.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1 = collect(range(set.lb1, stop = set.ub1, length = num_supports1))
    supports2 = collect(range(set.lb2, stop = set.ub2, length = num_supports2))
    return round.([supports1; supports2], sigdigits = sig_fig)
end

# output


```
Now automatic support generation is enabled, for example:
```jldoctest set_ext
julia> par = @infinite_parameter(model, set = DisjointSet(0, 2, 3, 4), num_supports = 10)
noname

julia> supports(par)
10-element Array{Float64,1}:
 0.0
 0.33333
 0.66667
 1.0
 1.3333
 1.6667
 2.0
 3.0
 3.5
 4.0
```
Note that in some cases [`generate_support_values`](@ref) might not be extendable
since returning supports for a single infinite parameter may not be appropriate.
This is typically true when a single set is used for multi-dimensional infinite
parameters. In these cases, [`generate_and_add_supports!`](@ref) will need to
be extended. This is exemplified with multivariate distribution sets that
correspond to an array of parameters where `generate_and_add_supports!` is
extended to generate the supports and then add the appropriate slices to each
dependent infinite parameter.

Next we can enable typical measure definition by extending
[`InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs`](@ref) for our
new set type. This function creates support values and corresponding coefficients
for a measure using the `method` argument if appropriate. Also, if the `method`
argument is used then we'll need to register valid method calls for our new set
via [`register_eval_method`](@ref). If we wish to ignore the `method` then we'll
need to set `check_method = false` when calling [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}))
or using [`set_measure_defaults`](@ref). Continuing our example, we obtain:
```jldoctest set_ext; output = false
const JuMPC = JuMP.Containers

function InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs(
    set::DisjointSet,
    params::Union{InfiniteOpt.ParameterRef, AbstractArray{<:InfiniteOpt.ParameterRef}},
    num_supports::Int,
    lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
    ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
    method::Function)::Tuple
    length_ratio = (set.ub1 - set.lb1) / (set.ub1 - set.lb1 + set.ub2 - set.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1, coeffs1 = method(set.lb1, set.ub1, num_supports1)
    supports2, coeffs2 = method(set.lb2, set.ub2, num_supports2)
    return ([supports1; supports2], [coeffs1; coeffs2])
end

register_eval_method(model, DisjointSet, [mc_sampling, gauss_legendre])

# output


```
Now we can define measures as normal:
```jldoctest set_ext
julia> mref = measure(t^2, t)
measure(t²)
```

Finally, we can extend the appropriate `JuMP` upper and lower bound functions
if desired which are:
- [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::AbstractInfiniteSet))
- [`JuMP.lower_bound`](@ref JuMP.lower_bound(::AbstractInfiniteSet))
- [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::AbstractInfiniteSet, ::Number))
- [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::AbstractInfiniteSet))
- [`JuMP.upper_bound`](@ref JuMP.upper_bound(::AbstractInfiniteSet))
- [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::AbstractInfiniteSet, ::Number))
However, if we want `has_lower_bound = false` and `has_upper_bound = false` then
no extension is needed. For our current example we won't do this since lower
and upper bounds aren't exactly clear for a disjoint interval. Please refer to
the template in `./InfiniteOpt/test/extensions/infinite_set.jl` to see how this
is done.

## Measure Evaluation Techniques

## Measure Data
Measures are used to evaluate over infinite domains. Users may wish to employ
measure abstractions that cannot be readily represented with coefficients and
discretized supports, and thus may wish to extend `InfiniteOpt`'s
measure framework to accommodate other paradigms. This can be accomplished my
implementing a user-defined measure data structure that inherits from
[`AbstractMeasureData`](@ref). A template for how such an extension is
accomplished is provided in `./InfiniteOpt/test/extensions/measure_data.jl`. The
extension steps employed are:
1. Define the new data struct inheriting from [`AbstractMeasureData`](@ref) (required)
2. Extend [`InfiniteOpt.parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) (required)
3. Extend [`InfiniteOpt.expand_measure`](@ref) (required)
4. Extend [`InfiniteOpt.supports`](@ref supports(::AbstractMeasureData)) (required if parameter supports are employed in any way)
5. Extend [`InfiniteOpt.measure_data_in_hold_bounds`](@ref) (enables hold variable bound checking with measures)
6. Extend [`InfiniteOpt.measure_name`](@ref) (enables meaningful measure naming)
7. Make simple measure constructor wrapper of [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)) to ease definition.

To illustrate how this process can be done, let's consider extending `InfiniteOpt`
to include measure for assessing the variance of random expressions. The
variance of an expression ``f(x, \xi)`` where ``x \in \mathbb{R}^n`` are hold
variables and ``\xi \in \mathbb{R}^m`` are random infinite parameters is defined:
```math
\mathbb{V}[f(x, \xi)] = \mathbb{E}\left[(f(x, \xi) - \mathbb{E}[f(x, \xi)])^2 \right].
```
Note, we could just accomplish this by nested use of [`expect`](@ref), but we
implement this example to illustrate the mechanics of extension.

First, let's define our new `struct` inheriting from `AbstractMeasureData`:
```jldoctest measure_data; output = false
using InfiniteOpt, JuMP, Distributions

const JuMPC = JuMP.Containers

struct DiscreteVarianceData <: AbstractMeasureData
    parameter_refs::Union{ParameterRef, JuMPC.SparseAxisArray{<:ParameterRef}}
    supports::Vector
    name::String
    # constructor
    function DiscreteVarianceData(
        parameter_refs::Union{ParameterRef, AbstractArray{<:ParameterRef}},
        supports::Vector, name::String = "Var")
        # convert input as necessary to proper array format
        if parameter_refs isa AbstractArray
            parameter_refs = convert(JuMPC.SparseAxisArray, parameter_refs)
            supports = [convert(JuMPC.SparseAxisArray, arr) for arr in supports]
        end
        return new(parameter_refs, supports, name)
    end
end

# output


```

We have defined our data type, so let's extend the measure data query
methods to enable its definition. These include
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)),
[`supports`](@ref supports(::AbstractMeasureData)), and
[`measure_name`](@ref) (optional):
```jldoctest measure_data; output = false
function InfiniteOpt.parameter_refs(data::DiscreteVarianceData)
    return data.parameter_refs
end

function InfiniteOpt.supports(data::DiscreteVarianceData)::Vector
    return data.supports
end

function InfiniteOpt.measure_name(data::DiscreteVarianceData)::String
    return data.name
end

# output


```
Note that extending `supports` is not needed for abstractions that don't involve
discretization of the infinite parameter(s), such as the case for outer certain
outer approximation techniques. Our extension is now sufficiently constructor to
allow us to define out new variance measure via
[`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)). For
example,
```jldoctest measure_data; setup = :(using Random; Random.seed!(42))
# Setup the infinite model
model = InfiniteModel()
@infinite_parameter(model, xi in Normal(), num_supports = 2) # few for simplicity
@infinite_variable(model, y(xi))
@hold_variable(model, z)

# Define out new variance measure
data = DiscreteVarianceData(xi, supports(xi))
mref = measure(2y + z, data)

# output
Var(2 y(xi) + z)
```
Thus, we can define measure references that employ this our new data type.

We can define variance measures now, but now let's extend
[`expand_measure`](@ref) so that they can be expanded into finite expressions:
```jldoctest measure_data; output = false
function InfiniteOpt.expand_measure(expr::JuMP.AbstractJuMPScalar,
                                    data::DiscreteVarianceData,
                                    write_model::JuMP.AbstractModel,
                                    point_mapper::Function)::JuMP.AbstractJuMPScalar
    # define the expectation data
    expect_data = DiscreteMeasureData(
                      data.parameter_refs,
                      1 / length(data.supports) * ones(length(data.supports)),
                      data.supports, name = data.name)
    # define the mean
    mean = measure(expr, expect_data)
    # return the expansion of the variance using the data mean
    return expand_measure((copy(expr) - mean)^2, expect_data, write_model, point_mapper)
end

# output


```
Notice that we reformulated our abstraction in terms of measures with
[`DiscreteMeasureData`](@ref) so that we could leverage the existing
[`expand_measure`](@ref) library. Now, new measure type can be expanded and
moreover infinite models using this new type can be optimized. Let's try
expanding the measure we already defined:
```jldoctest measure_data
julia> expand(mref)
2 y(-0.55603)² + 2 y(-0.44438)² + 2 z*y(-0.55603) + 2 z*y(-0.44438) - 4 y(-0.55603)² - 4 y(-0.44438)*y(-0.55603) - 4 z*y(-0.55603) + 0 z² - 2 y(-0.55603)*z - 2 y(-0.44438)*z + y(-0.55603)² + y(-0.44438)²
```

Finally, as per recommendation let's make a wrapper method to make defining
variance measures more convenient:
```jldoctest measure_data; output = false
function variance(expr::Union{JuMP.GenericAffExpr, GeneralVariableRef},
                  params::Union{ParameterRef, AbstractArray{<:ParameterRef}};
                  name::String = "Var", num_supports::Int = 10,
                  use_existing::Bool = false)::MeasureRef
    # get the supports
    if use_existing
        supps = supports.(params)
    else
        supps = generate_support_values(infinite_set(first(params)),
                                           num_supports = num_supports)
    end
    # make the data
    data = DiscreteVarianceData(params, supps, name)
    # built the measure
    return measure(expr, data)
end

# output

variance (generic function with 1 method)
```
Notice in this case that we only permit linear expressions for `expr` since
it will be squared by our new measure and we currently only support quadratic
expressions. (This could be overcome by defining place a place holder variable
for `expr`.

Now let's use our constructor to repeat the above measure example:
```jldoctest measure_data
julia> expand(variance(2y + z, xi, use_existing = true))
2 y(-0.55603)² + 2 y(-0.44438)² + 2 z*y(-0.55603) + 2 z*y(-0.44438) - 4 y(-0.55603)² - 4 y(-0.44438)*y(-0.55603) - 4 z*y(-0.55603) + 0 z² - 2 y(-0.55603)*z - 2 y(-0.44438)*z + y(-0.55603)² + y(-0.44438)²
```

We have done it! Now go and extend away!

## [Optimizer Models] (@id extend_optimizer_model)
`InfiniteOpt` provides a convenient interface and abstraction for modeling
infinite dimensional optimization problems. By default, `InfiniteModel`s are
reformulated into a solvable `JuMP.Model` (referred to as an optimizer model)
via `TranscriptionOpt` which discretizes the model in accordance with the
infinite parameter supports. However, uses my wish to employ some other
reformulation method to produce the optimizer model. This section will explain
how this can be done in `InfiniteOpt`. A template for implementing this
extension is provided in `./InfiniteOpt/test/extensions/optimizer_model.jl`.
Our default sub-package `InfiniteOpt.TranscriptionOpt` also serves as a good
example. 

## Wrapper Packages
