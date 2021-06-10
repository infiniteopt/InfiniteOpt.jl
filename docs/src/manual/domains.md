# [Infinite Domains](@id infinite_domains_manual)
A technical manual for infinite domains in `InfiniteOpt`. See the respective 
[guide](@ref infinite_domains_docs) for more information.

## Domain Types
```@docs
AbstractInfiniteDomain
InfiniteScalarDomain
IntervalDomain
UniDistributionDomain
InfiniteArrayDomain
MultiDistributionDomain
CollectionDomain
```

## Domain Methods
```@docs
collection_domains
JuMP.has_lower_bound(::AbstractInfiniteDomain)
JuMP.lower_bound(::AbstractInfiniteDomain)
JuMP.set_lower_bound(::AbstractInfiniteDomain, ::Real)
JuMP.has_upper_bound(::AbstractInfiniteDomain)
JuMP.upper_bound(::AbstractInfiniteDomain)
JuMP.set_upper_bound(::AbstractInfiniteDomain, ::Real)
```

## Support Point Labels
```@docs
AbstractSupportLabel
All
PublicLabel
UserDefined
UniformGrid
SampleLabel
MCSample
WeightedSample
Mixture
UniqueMeasure
InfiniteOpt.generate_unique_label
MeasureBound
InternalLabel
InfiniteOpt.MeasureToolbox.InternalGaussLobatto
```

## Support Point Methods
```@docs
supports_in_domain
generate_supports
InfiniteOpt.generate_support_values
```
