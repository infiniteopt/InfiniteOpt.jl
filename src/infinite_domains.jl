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

################################################################################
#                        SUPPORT AND LABEL GENERATION
################################################################################
"""
    AbstractSupportLabel

An abstract type for support label types. These are used to distinguish different 
kinds of supports that are added to infinite parameters.
"""
abstract type AbstractSupportLabel end 

"""
    All <: AbstractSupportLabel

This support label is unique in that it isn't associated with a particular set of 
supports, but rather is used used to indicate that all supports should be used.
"""
struct All <: AbstractSupportLabel end

# Filler label for NoGenerativeSupports
struct _NoLabel <: AbstractSupportLabel end

"""
    PublicLabel <: AbstractSupportLabel

An abstract type used to denote that labels that should be given to the user by 
default.
"""
abstract type PublicLabel <: AbstractSupportLabel end 

"""
    UserDefined <: PublicLabel

A support label for supports that are supplied by the user directly to an infinite 
parameter. 
"""
struct UserDefined <: PublicLabel end 

"""
    UniformGrid <: PublicLabel

A support label for supports that are generated uniformly accross a given interval.
"""
struct UniformGrid <: PublicLabel end

"""
    SampleLabel <: PublicLabel

An abstract type for labels of supports that are generated via some sampling technique.
"""
abstract type SampleLabel <: PublicLabel end

"""
    MCSample <: SampleLabel

A support label for supports that are generated via Monte Carlo Sampling.
"""
struct MCSample <: SampleLabel end 

"""
    WeightedSample <: SampleLabel

A support label for supports that are generated by sampling from a statistical 
distribution.
"""
struct WeightedSample <: SampleLabel end 

"""
    Mixture <: PublicLabel

A support label for multi-dimensional supports that are generated from a variety 
of methods.
"""
struct Mixture <: PublicLabel end

"""
    UniqueMeasure{S::Symbol} <: PublicLabel

A support label for supports that are provided from the `DiscreteMeasureData` 
associated with a measure where a unique label is generated to distinguish those 
supports. This is done by invoking [`generate_unique_label`](@ref).
"""
struct UniqueMeasure{S} <: PublicLabel end

"""
    MeasureBound <: PublicLabel
    
A support label for supports that are generated using the upper and lower bounds
for `FunctionalDiscreteMeasureData`.
"""
struct MeasureBound <: PublicLabel end

"""
    InternalLabel <: AbstractSupportLabel

An abstract type for support labels that are associated with supports that should 
not be reported to the user by default.
"""
abstract type InternalLabel <: AbstractSupportLabel end 

"""
    generate_unique_label()::Type{UniqueMeasure}

Generate and return a unique support label for measures.
"""
function generate_unique_label()::DataType
    return UniqueMeasure{gensym()}
end

# Define default values of sig_digits and num_supports keywords
const DefaultSigDigits = 12
const DefaultNumSupports = 10

# a user interface of generate_support_values
"""
    generate_supports(domain::AbstractInfiniteDomain
                      [method::Type{<:AbstractSupportLabel}];
                      [num_supports::Int = DefaultNumSupports,
                      sig_digits::Int = DefaultSigDigits]
                      )::Tuple{Array{<:Real}, DataType}

Generate `num_supports` support values with `sig_digits` significant digits in
accordance with `domain` and return them along with the correct generation label(s).
`IntervalDomain`s generate supports uniformly with label `UniformGrid` and
distribution domains generate them randomly accordingly to the
underlying distribution. Moreover, `method` indicates the generation method that
should be used. These `methods` correspond to parameter support labels. Current
labels that can be used as generation methods include (but may not be defined
for certain domain types):
- [`MCSample`](@ref): Uniformly distributed Monte Carlo samples.
- [`WeightedSample`](@ref): Monte Carlo samples that are weighted by an underlying PDF.
- [`UniformGrid`](@ref): Samples that are generated uniformly over the domain.

Extensions that employ user-defined infinite domain types and/or methods
should extend [`generate_support_values`](@ref) to enable this. Errors if the
`domain` type and /or methods are unrecognized. This is intended as an internal
method to be used by methods such as [`generate_and_add_supports!`](@ref).
"""
function generate_supports(domain::AbstractInfiniteDomain;
                           num_supports::Int = DefaultNumSupports,
                           sig_digits::Int = DefaultSigDigits
                           )::Tuple
    return generate_support_values(domain, num_supports = num_supports,
                                   sig_digits = sig_digits)
end

# 2 arguments
function generate_supports(domain::AbstractInfiniteDomain,
                           method::Type{<:AbstractSupportLabel};
                           num_supports::Int = DefaultNumSupports,
                           sig_digits::Int = DefaultSigDigits
                           )::Tuple
    return generate_support_values(domain, method,
                                   num_supports = num_supports,
                                   sig_digits = sig_digits)
end

"""
    generate_support_values(domain::AbstractInfiniteDomain,
                            [method::Type{MyMethod} = MyMethod];
                            [num_supports::Int = DefaultNumSupports,
                            sig_digits::Int = DefaultSigDigits]
                            )::Tuple{Array{<:Real}, Symbol}

A multiple dispatch method for [`generate_supports`](@ref). This will return
a tuple where the first element are the supports and the second is their
label. This can be extended for user-defined infinite domains and/or generation
methods. When defining a new domain type the default method dispatch should
make `method` an optional argument (making it the default). Otherwise, other
method dispatches for a given domain must ensure that `method` is positional
argument without a default value (contrary to the definition above). Note that the 
`method` must be a subtype of either [`PublicLabel`](@ref) or [`InternalLabel`](@ref).
"""
function generate_support_values(domain::AbstractInfiniteDomain,
                                 args...; kwargs...)
    if isempty(args)
        error("`generate_support_values` has not been extended for infinite domains " * 
              "of type `$(typeof(domain))`. This automatic support generation is not " * 
              "implemented.")
    else
        error("`generate_support_values` has not been extended for infinite domains " * 
              "of type `$(typeof(domain))` with the generation method `$(args[1])`. " * 
              "This automatic support generation is not implemented.")
    end
end

# IntervalDomain and UniformGrid
function generate_support_values(domain::IntervalDomain,
                                 method::Type{UniformGrid} = UniformGrid;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, DataType}
    lb = JuMP.lower_bound(domain)
    ub = JuMP.upper_bound(domain)
    new_supports = round.(range(lb, stop = ub, length = num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# IntervalDomain and MCSample
function generate_support_values(domain::IntervalDomain,
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, DataType}
    lb = JuMP.lower_bound(domain)
    ub = JuMP.upper_bound(domain)
    dist = Distributions.Uniform(lb, ub)
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# UniDistributionDomain and MultiDistributionDomain (with multivariate only)
function generate_support_values(
    domain::Union{UniDistributionDomain, MultiDistributionDomain{<:Distributions.MultivariateDistribution}},
    method::Type{WeightedSample} = WeightedSample;
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{<:Real}, DataType}
    dist = domain.distribution
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# UniDistributionDomain and MCSample 
function generate_support_values(
    domain::UniDistributionDomain,
    method::Type{MCSample};
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Vector{Float64}, DataType}
    return generate_support_values(domain, WeightedSample; num_supports = num_supports, 
                                   sig_digits = sig_digits)[1], method # TODO use an unwieghted sample...
end

# MultiDistributionDomain (matrix-variate distribution)
function generate_support_values(
    domain::MultiDistributionDomain{<:Distributions.MatrixDistribution},
    method::Type{WeightedSample} = WeightedSample;
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{Float64, 2}, DataType}
    dist = domain.distribution
    raw_supports = Distributions.rand(dist, num_supports)
    new_supports = Array{Float64}(undef, length(dist), num_supports)
    for i in 1:size(new_supports, 2)
        new_supports[:, i] = round.(reduce(vcat, raw_supports[i]),
                                    sigdigits = sig_digits)
    end
    return new_supports, method
end

# Generate the supports for a collection domain
function _generate_collection_supports(domain::CollectionDomain, num_supports::Int,
                                       sig_digits::Int)::Array{Float64, 2}
    domains = collection_domains(domain)
    # build the support array transpose to fill in column order (leverage locality)
    trans_supports = Array{Float64, 2}(undef, num_supports, length(domains))
    for i in eachindex(domains)
        @inbounds trans_supports[:, i] = generate_support_values(domains[i],
                                                   num_supports = num_supports,
                                                   sig_digits = sig_digits)[1]
    end
    return permutedims(trans_supports)
end

function _generate_collection_supports(domain::CollectionDomain,
                                       method::Type{<:AbstractSupportLabel},
                                       num_supports::Int,
                                       sig_digits::Int)::Array{Float64, 2}
    domains = collection_domains(domain)
    # build the support array transpose to fill in column order (leverage locality)
    trans_supports = Array{Float64, 2}(undef, num_supports, length(domains))
    for i in eachindex(domains)
        @inbounds trans_supports[:, i] = generate_support_values(domains[i],
                                                   method,
                                                   num_supports = num_supports,
                                                   sig_digits = sig_digits)[1]
    end
    return permutedims(trans_supports)
end

# CollectionDomain (IntervalDomains)
function generate_support_values(domain::CollectionDomain{IntervalDomain},
                                 method::Type{UniformGrid} = UniformGrid;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(domain, num_supports, sig_digits)
    return new_supports, method
end

function generate_support_values(domain::CollectionDomain{IntervalDomain},
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(domain, method, num_supports, sig_digits)
    return new_supports, method
end

# CollectionDomain (UniDistributionDomains)
function generate_support_values(domain::CollectionDomain{<:UniDistributionDomain},
                                 method::Type{WeightedSample} = WeightedSample;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(domain, num_supports, sig_digits)
    return new_supports, method
end

# CollectionDomain (InfiniteScalarDomains)
function generate_support_values(domain::CollectionDomain,
                                 method::Type{Mixture} = Mixture;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(domain, num_supports, sig_digits)
    return new_supports, method
end

# CollectionDomain (InfiniteScalarDomains) using purely MC sampling
# this is useful for measure support generation
function generate_support_values(domain::CollectionDomain,
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(domain, method, num_supports, sig_digits)
    return new_supports, method
end

# For label All: dispatch to default methods
function generate_support_values(domain::AbstractInfiniteDomain, ::Type{All};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits)
    return generate_support_values(domain, num_supports = num_supports,
                                   sig_digits = sig_digits)
end
