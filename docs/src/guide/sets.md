# [Infinite Sets] (@id infinite_sets_normal)
A guide and manual to the definition and use of infinite sets in
`InfiniteOpt`. The Datatypes and Methods sections at the end comprise the manual,
and the above sections comprise the guide.  

## Basic Usage
Interval sets are what characterize the behavior of infinite parameters in
`InfiniteOpt`, since they comprise the domains of infinite parameters. However,
most users will not need to work directly with infinite sets and can instead
focus on the use of infinite parameters as defined via [`@infinite_parameter`](@ref)
and as discussed on the [Infinite Parameters](@ref inf_par_page) page.

However, for convenience below we summary the infinite sets associated with
`InfiniteOpt`:

| Set Type        | Domain              |
|:---------------:|:-------------------:|
| IntervalSet     | ``[lb, ub]``        |
| UniDistributionSet | ``\sim \mathcal{D}``|
| MultiDistributionSet | ``\sim \mathcal{D}``|

## IntervalSets
The domain of a given infinite parameter is described by an infinite set
inherited from [`AbstractInfiniteSet`](@ref). `InfiniteOpt` natively supports
two such sets. The first is [`IntervalSet`](@ref) which describes a continuous
interval from some lower bound up to some upper bound. Typically, this range
is inclusive of the bounds. Such sets often arise for parameters that pertain to
time and/or spatial position. For example, to define a position interval
``[-2, 2]`` we would use:
```jldoctest; setup = :(using InfiniteOpt)
julia> set = IntervalSet(-2, 2)
[-2, 2]
```
Note that (semi-)infinite bounds are acceptable, as shown in the following example:
```jldoctest
julia> infinite_set = IntervalSet(-Inf, Inf)
[-Inf, Inf]
```

## DistributionSets
The second kind of set is that of [`UniDistributionSet`](@ref) which contains
univariate distribution from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl)
that characterizes an infinite parameter. Similarly, [`MultiDistributionSet`](@ref) contains 
multivariate distributions, wchih are only valid when
defining an array of parameters with equivalent dimensions. For example, let's
make a `UniDistributionSet` that depends on a Beta distribution:
```jldoctest; setup = :(using InfiniteOpt, Distributions)
julia> set = UniDistributionSet(Beta(2,2))
Beta{Float64}(α=2.0, β=2.0)
```
`MultiDistributionSet`s can be created in a similar way. For example, we can
make a `MultiDistributionSet` that depends on a 2-D normal distribution as follows:
```jldoctest; setup = :(using InfiniteOpt, Distributions)
julia> dist = MvNormal([0., 0.], [1. 0.; 0. 2.]);

julia> set = MultiDistributionSet(dist)
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)
```
User-defined distributions are also permissible so long as they are created in
accordance with `Distributions.jl`.

## Bound Query/Modification Methods for Infinite Sets
Once an infinite set is created, one can query the lower bound and upper bound of the set 
similar to how one queries the bounds of a `JuMP` variable. Thus, the functions 
[`JuMP.has_lower_bound`](@ref), [`JuMP.has_upper_bound`](@ref), [`JuMP.lower_bound`](@ref), [`JuMP.upper_bound`](@ref) 
are all applicable to infinite sets mentioned above. For example, for an [`IntervalSet`](@set)
`[-2, 2]` we can query the bound information as follows:
```jldoctest; setup = :(using InfiniteOpt)
julia> set = IntervalSet(-2, 2);

julia> has_lower_bound(set)
true

julia> has_upper_bound(set)
true

julia> lower_bound(set)
-2.0

julia> upper_bound(set)
2.0
```
In addition, we can also apply [`JuMP.set_lower_bound`](@ref) and [`JuMP.set_upper_bound`](@ref) 
to [`IntervalSet`](@set)s to generate a new set with updated bounds. Note that this will not modify the
original set. For example, we can change the bounds of the set `[-2, 2]` as follows:
```jldoctest; setup = :(using InfiniteOpt; set = IntervalSet(-2, 2))
julia> set_lower_bound(set, -1)
[-1, 2]

julia> set_upper_bound(set, 1)
[-2, 1]
```

## Support Generation for Infinite Sets
`InfiniteOpt` provides a systematic interface to generate support points for infinite sets.
This is crucial as support generation decides how each infinite-dimensional parameter, which is subject
to certain infinite set, is discretized later in the transcription stage. The interface will allow users
to automatically generate support points using our default methods. Later we will also show that users can 
also input support points manually for an infinite parameter.

In `InfiniteOpt` supports can be generated via [`generate_supports`](@ref) function. For example, let's 
generate 5 equidistant support points for the `IntervalSet` [-2, 2]:
```jldoctest; setup = :(using InfiniteOpt; set = IntervalSet(-2, 2))
julia> supps, label = generate_supports(set, num_supports = 5)
([-2.0, -1.0, 0.0, 1.0, 2.0], :uniform_grid)

julia> supps
5-element Array{Float64,1}:
 -2.0
 -1.0
  0.0
  1.0
  2.0

julia> label
:uniform_grid
```
Note that the number of supports generated is specified via
`num_supports` keyword argument, which will take a default value of 10 if not specified. 
The function `generate_supports` returns a vector of the supports generated, and a label that symbolizes
the underlying method. In this case the label returned is `UniformGrid`, which is the default 
support generation method for `IntervalSet`s. Another support generation method implemented for `IntervalSet`s
is `MCSample`, which is to sample from a uniform distribution over the interval. To use this mehtod, users
need to specify a second positional argument, as shown in the following example:
```jldoctest; setup = :(using InfiniteOpt, Random; Random.seed!(0); set = IntervalSet(-2, 2))
julia> generate_supports(set, MCSample, num_supports = 5)
([1.29459003191, 1.64142615171, -1.34173680747, -1.29068461413, -0.884479562675], :mc_sample)
```
In this case, the returned label is `MCSample`, instead of `UniformGrid`.

`generate_supports` can also be applied to `DistributionSets`. The default (and currently only) method
implemented for `DistributionSets` is `WeightedSample`, which generates Monte Carlo samples that are 
weighted based on the underlying probability density function of the distribution. 
For example, a set of support points for a 2-D normal distribution can be generated as follows:
```setup = :(using InfiniteOpt, Random; Random.seed!(0))
julia> dist = MvNormal([0., 0.], [1. 0.;0. 2.]);

julia> set = MultiDistributionSet(dist);

julia> supps, label = generate_supports(set, num_supports = 3)
([0.679107426036 -0.353007400301 0.586617074633; 1.17155358277 -0.190712174623 0.420496392851], :weighted_sample)

julia> supps
2×3 Array{Float64,2}:
 0.679107  -0.353007  0.586617
 1.17155   -0.190712  0.420496
```

For those who are interested in coding up their own support generation functions, [`generate_supports`](@ref) is
an inteface that calls the proper [`generate_support_values`](@ref) function based on the type of set and value of mehtod.
Therefore, to use custom support generation methods, users can implement extensions for [`generate_support_values`](@ref) 
with a different method label from the existing methods. See [Extensions](@ref) for full details.

## User Defined Sets
Furthermore, custom infinite sets that inherit `AbstractInfiniteSet` can also
be defined. See [Extensions](@ref) for more information.

## Datatypes
```@index
Pages   = ["sets.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractInfiniteSet
InfiniteScalarSet
IntervalSet
UniDistributionSet
InfiniteArraySet
MultiDistributionSet
CollectionSet
```

## Methods
```@index
Pages   = ["sets.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
collection_sets
JuMP.has_lower_bound(::AbstractInfiniteSet)
JuMP.lower_bound(::AbstractInfiniteSet)
JuMP.set_lower_bound(::AbstractInfiniteSet, ::Real)
JuMP.has_upper_bound(::AbstractInfiniteSet)
JuMP.upper_bound(::AbstractInfiniteSet)
JuMP.set_upper_bound(::AbstractInfiniteSet, ::Real)
supports_in_set
generate_supports
InfiniteOpt.generate_support_values
```
