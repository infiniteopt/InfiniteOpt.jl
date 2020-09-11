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
Note that infinite bounds are acceptable.

## UniDistributionSets
The second kind of set is that of [`UniDistributionSet`](@ref) which contains
distribution from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl)
that characterizes an infinite parameter. Univariate and multivariate
distributions can be used, but multivariate distributions are only valid when
defining an array of parameters with equivalent dimensions. For example, let's
make a `UniDistributionSet` that depends on a Beta distribution:
```jldoctest; setup = :(using InfiniteOpt, Distributions)
julia> set = UniDistributionSet(Beta(2,2))
Beta{Float64}(α=2.0, β=2.0)
```
User-defined distributions are also permissible so long as they are created in
accordance with `Distributions.jl`.

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
