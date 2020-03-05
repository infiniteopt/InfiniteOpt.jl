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
following steps employed are:
1. Define the new `struct` infinite set type (only thing required as bare minimum)
2. Extend [`InfiniteOpt.supports_in_set`](@ref) (enables error checking of supports)
3. Extend [`InfiniteOpt.generate_support_values`](@ref) (enables support generation via `num_supports` keyword arguments)
4. Extend [`InfiniteOpt.MeasureEvalMethods.measure_dispatch`](@ref) (enables measure evaluation methods)
5. If a lower bound and upper bound can be reported, extend `JuMP` lower bound and upper bound methods

As an example, let's create a disjoint interval set as an infinite set type.
This corresponds to the set ``[lb_1, ub_1] \cup [lb_2, ub_2]`` where
``ub_1 \leq lb_2``. First, we need to create the `DataType` with inheritance from [`AbstractInfiniteSet`](@ref):
```jldoctest set_ext
using InfiniteOpt, JuMP

struct DisjointSet <: InfiniteOpt.AbstractInfiniteSet
    lb1::Float64
    ub1::Float64
    lb2::Float64
    ub2::Float64
    # constructor
    function DisjointSet(lb1, ub1, lb2, ub2)
        if lb1 > ub1 || lb2 > ub2 || ub1 > lb2
            error("Invalid bounds")
        end
        return new(convert(Float64, lb1), convert(Float64, ub1),
                   convert(Float64, lb2), convert(Float64, ub2))
    end
end
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

## Measure Evaluation Techniques

## Measure Data

## [Optimizer Models] (@id extend_optimizer_model)

## Wrapper Packages
