# [Infinite Domains](@id infinite_domains_docs)
A guide for infinite domains in `InfiniteOpt`. See the respective 
[technical manual](@ref infinite_domains_manual) for more details.

!!! note 
    Previous versions of `InfiniteOpt` referred to infinite domains as infinite 
    sets. Hence, all the methods and datatypes have been updated accordingly. 

## Basic Usage
Interval domains are what characterize the behavior of infinite parameters in 
`InfiniteOpt`, since they comprise the domains of infinite parameters. However, 
most users will not need to work directly with infinite domains and can instead 
focus on the use of infinite parameters as defined via [`@infinite_parameter`](@ref) 
and as discussed on the [Infinite Parameters](@ref inf_par_docs) page.

However, for convenience below we summarize the infinite domains associated with 
`InfiniteOpt`:

| Domain Type                       | Domain                                     | Type                        |
|:------------------------------:|:------------------------------------------:|:---------------------------:|
| [`IntervalDomain`](@ref)          | ``[lb, ub]``                               | [`InfiniteScalarDomain`](@ref) |
| [`UniDistributionDomain`](@ref)   | ``\sim \mathcal{D} \subseteq \mathbb{R}``  | [`InfiniteScalarDomain`](@ref) |
| [`MultiDistributionDomain`](@ref) | ``\sim \mathcal{D} \subseteq \mathbb{R}^n``| [`InfiniteArrayDomain`](@ref)  |
| [`CollectionDomain`](@ref)        | Combination of Univariate Domains          | [`InfiniteArrayDomain`](@ref)  |

## Infinite Domain Classes
The domain of a given infinite parameter(s) is described by an infinite domain (domain) 
inherited from [`AbstractInfiniteDomain`](@ref). `InfiniteOpt` natively supports 
two domain sub-groups, namely [`InfiniteScalarDomain`](@ref)s and [`InfiniteArrayDomain`](@ref)s. 
These correspond to a single independent infinite parameter and a dependent multi-dimensional 
group of infinite parameters, respectively. We describe each group's natively 
supported domains below.

### Univariate Domains
Univariate infinite domains (i.e., [`InfiniteScalarDomain`](@ref)s) are one-dimensional 
domains (``\subseteq \mathbb{R}``) that describe the behavior of one single independent 
infinite parameter (i.e., infinite parameters made using `independent = true`). The 
two natively supported concrete types are [`IntervalDomain`](@ref)s and [`UniDistributionDomain`](@ref)s.

[`IntervalDomain`](@ref)s describe a continuous interval from some lower bound up to 
some upper bound. Where the range is inclusive of the bounds. Such domains often 
arise for parameters that pertain to time and/or spatial position. For example, 
to define a position interval ``[-2, 2]`` we would use:
```jldoctest; setup = :(using InfiniteOpt)
julia> domain = IntervalDomain(-2, 2)
[-2, 2]
```
Note that (semi-)infinite bounds are acceptable, as shown in the following example:
```jldoctest; setup = :(using InfiniteOpt)
julia> infinite_domain = IntervalDomain(-Inf, Inf)
[-Inf, Inf]
```

[`UniDistributionDomain`](@ref)s pertain to the co-domain of a univariate distribution. 
In other words, these correspond to the underlying distributions that characterize 
uncertain scalar parameters. These domains are compatible with any univariate 
distribution native to [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl).
For example, let's make a `UniDistributionDomain` that depends on a Beta distribution:
```jldoctest; setup = :(using InfiniteOpt)
julia> using Distributions;

julia> domain = UniDistributionDomain(Beta(2,2))
Beta{Float64}(α=2.0, β=2.0)
```
User-defined distributions are also permissible so long as they are created in
accordance with `Distributions.jl`.

### Multivariate Domains
Multivariate infinite domains (i.e., [`InfiniteArrayDomain`])(@ref)s are multi-dimensional 
domains that define the behavior of a group of dependent infinite parameters 
(i.e., an array of infinite parameters where `independent = false`). This is a 
unique feature to `InfiniteOpt` that enables a much richer set of possibilities 
for modeling infinite domain. Natively two domain types are supported: 
[`MultiDistributionDomain`](@ref)s and [`CollectionDomain`](@ref)s. 

[`MultiDistributionDomain`](@ref)s correspond to the co-domain of a multi-variate 
(or matrix-variate) distribution which characterizes the behavior of multi-dimensional 
uncertain parameters. Again, these correspond to any appropriate distribution 
defined in `Distributions.jl`. For example, we can make a `MultiDistributionDomain` 
that depends on a 2D normal distribution as follows:
```jldoctest; setup = :(using InfiniteOpt)
julia> using Distributions;

julia> dist = MvNormal([0., 0.], [1. 0.; 0. 2.]);

julia> domain = MultiDistributionDomain(dist)
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)
```

!!! note 
    The dimensions (shape) of a chosen distribution used in an `MultiDistriubtionDomain` 
    must match those of the corresponding infinite parameter array.

Finally, [`CollectionDomain`](@ref)s are a dependent collection of [`InfiniteScalarDomain`](@ref)s
that correspond to a group of infinite parameters that are treated dependently. 
This can be useful when the user wishes to have complete control over how the 
supports are generated for a group independent parameters where the default 
combinatorial approach in not wanted. For example, let's make a set of `IntervalDomain`s:
```jldoctest; setup = :(using InfiniteOpt)
julia> domain = CollectionDomain([IntervalDomain(-2, 2), IntervalDomain(-1, 4)])
CollectionDomain with 2 domains:
 [-2, 2]
 [-1, 4]
```
Now we could use this domain in define a two-dimensional infinite parameter of which 
we can have the freedom to define a non-combinatorial support grid.

## Bound Query/Modification Methods for Infinite Domains
Once an infinite domain is created, one can query the lower bound and upper bound of the domain 
similar to how one queries the bounds of a `JuMP` variable. Thus, the functions 
[`JuMP.has_lower_bound`](@ref), [`JuMP.has_upper_bound`](@ref), [`JuMP.lower_bound`](@ref), [`JuMP.upper_bound`](@ref) 
are all applicable to infinite domains mentioned above. For example, for an `IntervalDomain`
`[-2, 2]` we can query the bound information as follows:
```jldoctest; setup = :(using InfiniteOpt)
julia> domain = IntervalDomain(-2, 2);

julia> has_lower_bound(domain)
true

julia> has_upper_bound(domain)
true

julia> lower_bound(domain)
-2.0

julia> upper_bound(domain)
2.0
```
In addition, we can also apply [`JuMP.set_lower_bound`](@ref) and [`JuMP.set_upper_bound`](@ref) 
to `IntervalDomain`s to generate a new domain with updated bounds. Note that this will not modify the
original domain. For example, we can change the bounds of the set `[-2, 2]` as follows:
```jldoctest; setup = :(using InfiniteOpt; domain = IntervalDomain(-2, 2))
julia> set_lower_bound(domain, -1)
[-1, 2]

julia> set_upper_bound(domain, 1)
[-2, 1]
```

## Support Generation for Infinite Domains
`InfiniteOpt` provides a systematic interface to generate support points for 
infinite domains. This is crucial as support generation decides how each 
infinite-dimensional parameter, which is subject to certain infinite domain, is 
discretized later in the transcription stage. The interface will allow users to 
automatically generate support points using our default methods. Later we will 
also show that users can also input support points manually for an infinite 
parameter. Please note that these methods are called by the 
[`@infinite_parameter`](@ref) macro when the `num_supports` keyword is used. 
Thus, users typically will not need to use this interface directly.

In `InfiniteOpt` supports can be generated via [`generate_supports`](@ref) 
function. For example, let's generate 5 equidistant support points for the 
`IntervalDomain` [-2, 2]:
```jldoctest; setup = :(using InfiniteOpt; domain = IntervalDomain(-2, 2))
julia> supps, label = generate_supports(domain, num_supports = 5)
([-2.0, -1.0, 0.0, 1.0, 2.0], UniformGrid)
```
Note that the number of supports generated is specified via `num_supports` 
keyword argument, which will take a default value of 10 if not specified. The 
function `generate_supports` returns a vector of the supports generated, and a 
label that denotes the underlying method. In this case the label returned is 
`UniformGrid`, which is the default support generation method for 
`IntervalDomain`s. Another support generation method implemented for 
`IntervalDomain`s is `MCSample`, which is to sample from a uniform distribution 
over the interval. To use this method, users need to specify a second positional 
argument, as shown in the following example:
```jldoctest; setup = :(using InfiniteOpt, Random; Random.seed!(0); domain = IntervalDomain(-2, 2))
julia> generate_supports(domain, MCSample, num_supports = 5, sig_digits = 5)
([1.2946, 1.6414, -1.3417, -1.2907, -0.88448], MCSample)
```
In this case, the returned label is `MCSample`, instead of `UniformGrid`.

`generate_supports` can also be applied to `DistributionDomains`. The default 
(and currently only) method implemented for `DistributionDomains` is 
`WeightedSample`, which generates Monte Carlo samples that are weighted based on 
the underlying probability density function of the distribution. For example, a 
domain of support points for a 2D normal distribution can be generated as 
follows:
```setup = :(using InfiniteOpt, Random; Random.seed!(0))
julia> dist = MvNormal([0., 0.], [1. 0.;0. 2.]);

julia> domain = MultiDistributionDomain(dist);

julia> supps, label = generate_supports(domain, num_supports = 3)
([0.679107426036 -0.353007400301 0.586617074633; 1.17155358277 -0.190712174623 0.420496392851], WeightedSample)
```

For those who are interested in coding up their own support generation functions, 
[`generate_supports`](@ref) is an interface that calls the proper 
[`generate_support_values`](@ref) function based on the type of domain and value 
of method. Therefore, to use custom support generation methods, users can 
implement extensions for [`generate_support_values`](@ref) with a different 
method label from the existing methods. See [Extensions](@ref) for full details.

## User Defined Domains
Furthermore, custom infinite domains that inherit `AbstractInfiniteDomain` can 
also be defined. See [Extensions](@ref) for more information.
