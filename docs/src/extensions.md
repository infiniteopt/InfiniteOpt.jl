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
no extension is needed. For our current, example we won't do this since lower
and upper bounds aren't exactly clear for a disjoint interval. Please refer to
the template in `./InfiniteOpt/test/extensions/infinite_set.jl` to see how this
is done.

## Measure Evaluation Techniques

## Measure Data
Measures are used to evaluate over infinite domains. Users may wish to employ
measures that are strictly integrals and thus may wish to extend `InfiniteOpt`'s
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
7. Make simple measure constructor wrapper of [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)) to ease definition

TODO: Provide an extension example using CVaR.

## [Optimizer Models] (@id extend_optimizer_model)

## Wrapper Packages
