################################################################################
#                              ACCESSER METHODS
################################################################################
"""
    collection_domains(domain::AbstractInfiniteDomain)

Return the array of domains associated with a `CollectionDomain`. Error if the input
domain is not a `CollectionDomain`.
"""
function collection_domains(domain::AbstractInfiniteDomain)
    type = typeof(domain)
    error("`collection_domains` not defined for infinite domain of type $type.")
end
function collection_domains(domain::CollectionDomain{S}
                         )::Vector{S} where {S <: InfiniteScalarDomain}
    return domain.domains
end

################################################################################
#                               BASIC EXTENSIONS
################################################################################
Base.length(s::CollectionDomain)::Int = length(collection_domains(s))
Base.length(s::MultiDistributionDomain)::Int = length(s.distribution)
Base.length(s::InfiniteScalarDomain)::Int = 1

################################################################################
#                           SUPPORT VALIDITY METHODS
################################################################################
"""
    supports_in_domain(supports::Union{Real, Vector{<:Real}, Array{<:Real, 2}},
                    domain::AbstractInfiniteDomain)::Bool

Used to check if `supports` are in the domain of `domain`. Returns `true` if
`supports` are in domain of `domain` and returns `false` otherwise.
This is primarily an internal method for performing checks but can be extended
for user-defined domain types. Extending this is optional, but recommended where
possible. Note by fallback, this returns `true` for unrecognized domain types such
that an error won't be thrown.
"""
function supports_in_domain end

# IntervalDomain
function supports_in_domain(supports::Union{Real, Vector{<:Real}},
                         domain::IntervalDomain)::Bool
    min_support = minimum(supports)
    max_support = maximum(supports)
    if min_support < JuMP.lower_bound(domain) || max_support > JuMP.upper_bound(domain)
        return false
    end
    return true
end

# UnivariateDistribution domain
function supports_in_domain(supports::Union{Real, Vector{<:Real}},
                         domain::UniDistributionDomain)::Bool
    return all(Distributions.insupport(domain.distribution, supports))
end

# MultivariateDistribution domain
function supports_in_domain(
    supports::Array{<:Real},
    domain::MultiDistributionDomain{<:Distributions.MultivariateDistribution{S}}
    )::Bool where {S <: Distributions.ValueSupport}
    if length(domain.distribution) != size(supports, 1)
        error("Support dimensions does not match distribution dimensions.")
    end
    return all(Distributions.insupport(domain.distribution, supports))
end

# CollectionDomain
function supports_in_domain(supports::Array{<:Real, 2}, domain::CollectionDomain)::Bool
    domains = collection_domains(domain)
    if length(domains) != size(supports, 1)
        error("Support dimensions does not match CollectionDomain dimensions.")
    end
    for i in eachindex(domains)
        if !supports_in_domain(supports[i, :], domains[i])
            return false
        end
    end
    return true
end

# Fallback
function supports_in_domain(supports, domain::AbstractInfiniteDomain)::Bool
    return true
end

################################################################################
#                              LOWER BOUND METHODS
################################################################################
"""
    JuMP.has_lower_bound(domain::AbstractInfiniteDomain)::Bool

Return `Bool` indicating if `domain` has a lower bound that can be determined. This
should be extended for user-defined infinite domains. It defaults to `false`
for unrecognized domain types.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> has_lower_bound(domain)
true
```
"""
function JuMP.has_lower_bound(domain::AbstractInfiniteDomain)::Bool # fallback
    return false
end

# IntervalDomain
JuMP.has_lower_bound(domain::IntervalDomain)::Bool = true

# DistributionDomain (Univariate)
JuMP.has_lower_bound(domain::UniDistributionDomain)::Bool = true

# CollectionDomain
function JuMP.has_lower_bound(domain::CollectionDomain)::Bool
    for s in collection_domains(domain)
        if !JuMP.has_lower_bound(s)
            return false
        end
    end
    return true
end

"""
    JuMP.lower_bound(domain::AbstractInfiniteDomain)::Union{Real, Vector{<:Real}}

Return the lower bound of `domain` if one exists. This should be extended for
user-defined infinite domains if appropriate. Errors if `JuMP.has_lower_bound`
returns `false`. Extensions are enabled by `JuMP.has_lower_bound(domain)` and
`JuMP.lower_bound(domain)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> lower_bound(domain)
0.0
```
"""
function JuMP.lower_bound(domain::AbstractInfiniteDomain) # fallback
    type = typeof(domain)
    error("`JuMP.lower_bound` not defined for infinite domain of type $type.")
end

# IntervalDomain
JuMP.lower_bound(domain::IntervalDomain)::Float64 = domain.lower_bound

# DistributionDomain (Univariate)
function JuMP.lower_bound(domain::UniDistributionDomain)::Real
    return minimum(domain.distribution)
end

# CollectionDomain
function JuMP.lower_bound(domain::CollectionDomain)::Vector{<:Real}
    return [JuMP.lower_bound(i) for i in collection_domains(domain)]
end

"""
    JuMP.set_lower_bound(domain::AbstractInfiniteDomain,
                         lower::Union{Real, Vector{<:Real}})::AbstractInfiniteDomain

Set and return the lower bound of `domain` if such an operation makes sense. Errors
if the type of `domain` does not support this operation or has not been extended.
User-defined domain types should extend this if appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> set_lower_bound(domain, 0.5)
[0.5, 1]
```
"""
function JuMP.set_lower_bound(domain::AbstractInfiniteDomain,
                              lower::Union{Real, Vector{<:Real}}) # fallback
    type = typeof(domain)
    error("`JuMP.set_lower_bound` not defined for infinite domain of type $type.")
end

# IntervalDomain
function JuMP.set_lower_bound(domain::IntervalDomain, lower::Real)::IntervalDomain
    return IntervalDomain(lower, domain.upper_bound)
end

# DistributionDomain
function JuMP.set_lower_bound(domain::Union{UniDistributionDomain,
                                         MultiDistributionDomain},
                              lower::Union{Real, Vector{<:Real}})
    error("Cannot set the lower bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

# CollectionDomain
function JuMP.set_lower_bound(domain::CollectionDomain, lower::Vector{<:Real})::CollectionDomain
    domains = collection_domains(domain)
    new_domains = [JuMP.set_lower_bound(domains[i], lower[i]) for i in eachindex(domains)]
    return CollectionDomain(new_domains)
end

################################################################################
#                              UPPER BOUND METHODS
################################################################################
"""
    JuMP.has_upper_bound(domain::AbstractInfiniteDomain)::Bool

Return `Bool` indicating if `domain` has a upper bound that can be determined. This
should be extended for user-defined infinite domains. It defaults to `false`
for unrecognized domain types.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> has_upper_bound(domain)
true
```
"""
function JuMP.has_upper_bound(domain::AbstractInfiniteDomain)::Bool # fallback
    return false
end

# IntervalDomain
JuMP.has_upper_bound(domain::IntervalDomain)::Bool = true

# DistributionDomain (Univariate)
JuMP.has_upper_bound(domain::UniDistributionDomain)::Bool = true

# CollectionDomain
function JuMP.has_upper_bound(domain::CollectionDomain)::Bool
    for i in collection_domains(domain)
        if !JuMP.has_upper_bound(i)
            return false
        end
    end
    return true
end

"""
    JuMP.upper_bound(domain::AbstractInfiniteDomain)::Union{Real, Vector{<:Real}}

Return the upper bound of `domain` if one exists. This should be extended for
user-defined infinite domains if appropriate. Errors if `JuMP.has_upper_bound`
returns `false`. Extensions are enabled by `JuMP.has_upper_bound(domain)` and
`JuMP.upper_bound(domain)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> upper_bound(domain)
1.0
```
"""
function JuMP.upper_bound(domain::AbstractInfiniteDomain) # fallback
    type = typeof(domain)
    error("`JuMP.upper_bound` not defined for infinite domain of type $type.")
end

# IntervalDomain
JuMP.upper_bound(domain::IntervalDomain)::Float64 = domain.upper_bound

# DistributionDomain (Univariate)
function JuMP.upper_bound(domain::UniDistributionDomain)::Real
    return maximum(domain.distribution)
end

# CollectionDomain
function JuMP.upper_bound(domain::CollectionDomain)::Vector{<:Real}
    return [JuMP.upper_bound(i) for i in collection_domains(domain)]
end


"""
    JuMP.set_upper_bound(domain::AbstractInfiniteDomain,
                         upper::Real)::AbstractInfiniteDomain

Set and return the upper bound of `domain` if such an aoperation makes sense. Errors
if the type of `domain` does not support this operation or has not been extended.
User-defined domain types should extend this if appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> domain = InfiniteDomain(0, 1);

julia> set_upper_bound(domain, 0.5)
[0, 0.5]
```
"""
function JuMP.set_upper_bound(domain::AbstractInfiniteDomain, upper::Real) # fallback
    type = typeof(domain)
    error("`JuMP.set_lower_bound` not defined for infinite domain of type $type.")
end

# IntervalDomain
function JuMP.set_upper_bound(domain::IntervalDomain, upper::Real)::IntervalDomain
    return IntervalDomain(domain.lower_bound, upper)
end

# DistributionDomain
function JuMP.set_upper_bound(domain::Union{UniDistributionDomain,
                                         MultiDistributionDomain}, lower::Real)
    error("Cannot set the upper bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

# CollectionDomain
function JuMP.set_upper_bound(domain::CollectionDomain, lower::Vector{<:Real})::CollectionDomain
    domains = collection_domains(domain)
    new_domains = [JuMP.set_upper_bound(domains[i], lower[i]) for i in eachindex(domains)]
    return CollectionDomain(new_domains)
end
