################################################################################
#                              ACCESSER METHODS
################################################################################
"""
    collection_sets(set::AbstractInfiniteSet)

Return the array of sets associated with a `CollectionSet`. Error if the input
set is not a `CollectionSet`.
"""
function collection_sets(set::AbstractInfiniteSet)
    type = typeof(set)
    error("`collection_sets` not defined for infinite set of type $type.")
end
function collection_sets(set::CollectionSet{S}
                         )::Vector{S} where {S <: InfiniteScalarSet}
    return set.sets
end

################################################################################
#                               BASIC EXTENSIONS
################################################################################
Base.length(s::CollectionSet)::Int = length(collection_sets(s))
Base.length(s::MultiDistributionSet)::Int = length(s.distribution)
Base.length(s::InfiniteScalarSet)::Int = 1

################################################################################
#                           SUPPORT VALIDITY METHODS
################################################################################
"""
    supports_in_set(supports::Union{Real, Vector{<:Real}, Array{<:Real, 2}},
                    set::AbstractInfiniteSet)::Bool

Used to check if `supports` are in the domain of `set`. Returns `true` if
`supports` are in domain of `set` and returns `false` otherwise.
This is primarily an internal method for performing checks but can be extended
for user-defined set types. Extending this is optional, but recommended where
possible. Note by fallback, this returns `true` for unrecognized set types such
that an error won't be thrown.
"""
function supports_in_set end

# IntervalSet
function supports_in_set(supports::Union{Real, Vector{<:Real}},
                         set::IntervalSet)::Bool
    min_support = minimum(supports)
    max_support = maximum(supports)
    if min_support < JuMP.lower_bound(set) || max_support > JuMP.upper_bound(set)
        return false
    end
    return true
end

# UnivariateDistribution set
function supports_in_set(supports::Union{Real, Vector{<:Real}},
                         set::UniDistributionSet)::Bool
    return all(Distributions.insupport(set.distribution, supports))
end

# MultivariateDistribution set
function supports_in_set(
    supports::Array{<:Real},
    set::MultiDistributionSet{<:Distributions.MultivariateDistribution{S}}
    )::Bool where {S <: Distributions.ValueSupport}
    if length(set.distribution) != size(supports, 1)
        error("Support dimensions does not match distribution dimensions.")
    end
    return all(Distributions.insupport(set.distribution, supports))
end

# CollectionSet
function supports_in_set(supports::Array{<:Real, 2}, set::CollectionSet)::Bool
    sets = collection_sets(set)
    if length(sets) != size(supports, 1)
        error("Support dimensions does not match CollectionSet dimensions.")
    end
    for i in eachindex(sets)
        if !supports_in_set(supports[i, :], sets[i])
            return false
        end
    end
    return true
end

# Fallback
function supports_in_set(supports, set::AbstractInfiniteSet)::Bool
    return true
end

################################################################################
#                              LOWER BOUND METHODS
################################################################################
"""
    JuMP.has_lower_bound(set::AbstractInfiniteSet)::Bool

Return `Bool` indicating if `set` has a lower bound that can be determined. This
should be extended for user-defined infinite sets. It defaults to `false`
for unrecognized set types.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> has_lower_bound(set)
true
```
"""
function JuMP.has_lower_bound(set::AbstractInfiniteSet)::Bool # fallback
    return false
end

# IntervalSet
JuMP.has_lower_bound(set::IntervalSet)::Bool = true

# DistributionSet (Univariate)
JuMP.has_lower_bound(set::UniDistributionSet)::Bool = true

# CollectionSet
function JuMP.has_lower_bound(set::CollectionSet)::Bool
    for s in collection_sets(set)
        if !JuMP.has_lower_bound(s)
            return false
        end
    end
    return true
end

"""
    JuMP.lower_bound(set::AbstractInfiniteSet)::Union{Real, Vector{<:Real}}

Return the lower bound of `set` if one exists. This should be extended for
user-defined infinite sets if appropriate. Errors if `JuMP.has_lower_bound`
returns `false`. Extensions are enabled by `JuMP.has_lower_bound(set)` and
`JuMP.lower_bound(set)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> lower_bound(set)
0.0
```
"""
function JuMP.lower_bound(set::AbstractInfiniteSet) # fallback
    type = typeof(set)
    error("`JuMP.lower_bound` not defined for infinite set of type $type.")
end

# IntervalSet
JuMP.lower_bound(set::IntervalSet)::Float64 = set.lower_bound

# DistributionSet (Univariate)
function JuMP.lower_bound(set::UniDistributionSet)::Real
    return minimum(set.distribution)
end

# CollectionSet
function JuMP.lower_bound(set::CollectionSet)::Vector{<:Real}
    return [JuMP.lower_bound(i) for i in collection_sets(set)]
end

"""
    JuMP.set_lower_bound(set::AbstractInfiniteSet,
                         lower::Union{Real, Vector{<:Real}})::AbstractInfiniteSet

Set and return the lower bound of `set` if such an aoperation makes sense. Errors
if the type of `set` does not support this operation or has not been extended.
User-defined set types should extend this if appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> set_lower_bound(set, 0.5)
[0.5, 1]
```
"""
function JuMP.set_lower_bound(set::AbstractInfiniteSet,
                              lower::Union{Real, Vector{<:Real}}) # fallback
    type = typeof(set)
    error("`JuMP.set_lower_bound` not defined for infinite set of type $type.")
end

# IntervalSet
function JuMP.set_lower_bound(set::IntervalSet, lower::Real)::IntervalSet
    return IntervalSet(lower, set.upper_bound)
end

# DistributionSet
function JuMP.set_lower_bound(set::Union{UniDistributionSet,
                                         MultiDistributionSet},
                              lower::Union{Real, Vector{<:Real}})
    error("Cannot set the lower bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

# CollectionSet
function JuMP.set_lower_bound(set::CollectionSet, lower::Vector{<:Real})::CollectionSet
    sets = collection_sets(set)
    new_sets = [JuMP.set_lower_bound(sets[i], lower[i]) for i in eachindex(sets)]
    return CollectionSet(new_sets)
end

################################################################################
#                              UPPER BOUND METHODS
################################################################################
"""
    JuMP.has_upper_bound(set::AbstractInfiniteSet)::Bool

Return `Bool` indicating if `set` has a upper bound that can be determined. This
should be extended for user-defined infinite sets. It defaults to `false`
for unrecognized set types.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> has_upper_bound(set)
true
```
"""
function JuMP.has_upper_bound(set::AbstractInfiniteSet)::Bool # fallback
    return false
end

# IntervalSet
JuMP.has_upper_bound(set::IntervalSet)::Bool = true

# DistributionSet (Univariate)
JuMP.has_upper_bound(set::UniDistributionSet)::Bool = true

# CollectionSet
function JuMP.has_upper_bound(set::CollectionSet)::Bool
    for i in collection_sets(set)
        if !JuMP.has_upper_bound(i)
            return false
        end
    end
    return true
end

"""
    JuMP.upper_bound(set::AbstractInfiniteSet)::Union{Real, Vector{<:Real}}

Return the upper bound of `set` if one exists. This should be extended for
user-defined infinite sets if appropriate. Errors if `JuMP.has_upper_bound`
returns `false`. Extensions are enabled by `JuMP.has_upper_bound(set)` and
`JuMP.upper_bound(set)`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> upper_bound(set)
1.0
```
"""
function JuMP.upper_bound(set::AbstractInfiniteSet) # fallback
    type = typeof(set)
    error("`JuMP.upper_bound` not defined for infinite set of type $type.")
end

# IntervalSet
JuMP.upper_bound(set::IntervalSet)::Float64 = set.upper_bound

# DistributionSet (Univariate)
function JuMP.upper_bound(set::UniDistributionSet)::Real
    return maximum(set.distribution)
end

# CollectionSet
function JuMP.upper_bound(set::CollectionSet)::Vector{<:Real}
    return [JuMP.upper_bound(i) for i in collection_sets(set)]
end


"""
    JuMP.set_upper_bound(set::AbstractInfiniteSet,
                         upper::Real)::AbstractInfiniteSet

Set and return the upper bound of `set` if such an aoperation makes sense. Errors
if the type of `set` does not support this operation or has not been extended.
User-defined set types should extend this if appropriate.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP)
julia> set = InfiniteSet(0, 1);

julia> set_upper_bound(set, 0.5)
[0, 0.5]
```
"""
function JuMP.set_upper_bound(set::AbstractInfiniteSet, upper::Real) # fallback
    type = typeof(set)
    error("`JuMP.set_lower_bound` not defined for infinite set of type $type.")
end

# IntervalSet
function JuMP.set_upper_bound(set::IntervalSet, upper::Real)::IntervalSet
    return IntervalSet(set.lower_bound, upper)
end

# DistributionSet
function JuMP.set_upper_bound(set::Union{UniDistributionSet,
                                         MultiDistributionSet}, lower::Real)
    error("Cannot set the upper bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

# CollectionSet
function JuMP.set_upper_bound(set::CollectionSet, lower::Vector{<:Real})::CollectionSet
    sets = collection_sets(set)
    new_sets = [JuMP.set_upper_bound(sets[i], lower[i]) for i in eachindex(sets)]
    return CollectionSet(new_sets)
end

################################################################################
#                        SUPPORT AND LABEL GENERATION
################################################################################
# Define generation labels to be stored with supports
const All = :all
const UserDefined = :user_defined
const MCSample = :mc_sample
const WeightedSample = :weighted_sample
const UniformGrid = :uniform_grid
const Mixture = :mixture

# Define default values of sig_digits and num_supports keywords
const DefaultSigDigits = 12
const DefaultNumSupports = 10

"""
    generate_supports(set::AbstractInfiniteSet
                      method::Union{Symbol, Nothing} = Nothing;
                      [num_supports::Int = DefaultNumSupports,
                      sig_digits::Int = DefaultSigDigits]
                      )::Tuple{Array{<:Real}, Symbol}

Generate `num_supports` support values with `sig_digits` significant digits in
accordance with `set` and return them along with the correct generation label(s).
`IntervalSet`s generate supports uniformly with label `Grid` and
distribution sets generate them randomly accordingly to the
underlying distribution. Extensions that employ user-defined infinite set types
should extend according to the new set type. Errors if the `set` is a type that
has not been explicitly extended. This is intended as an internal method to be
used by [`generate_and_add_supports!`](@ref) and [`build_parameter`](@ref).
"""
# a user interface of generate_support_values
function generate_supports(set::AbstractInfiniteSet,
                           method::Union{Symbol, Nothing} = nothing;
                           num_supports::Int = DefaultNumSupports,
                           sig_digits::Int = DefaultSigDigits)::Tuple
    if method === nothing
        return generate_support_values(set,
                                       num_supports = num_supports,
                                       sig_digits = sig_digits)
    else
        return generate_support_values(set, Val(method),
                                       num_supports = num_supports,
                                       sig_digits = sig_digits)
    end
end

# Fallback
function generate_support_values(set::AbstractInfiniteSet,
                                 method::Val = Val(:nothing); kwargs...)
     error("Method `$(typeof(method).parameters[1])` is not supported for " *
           "set of type $(typeof(set))")
end

# IntervalSet and UniformGrid
function generate_support_values(set::IntervalSet,
                                 ::Val{UniformGrid} = Val(UniformGrid);
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, Symbol}
    lb = JuMP.lower_bound(set)
    ub = JuMP.upper_bound(set)
    new_supports = round.(range(lb, stop = ub, length = num_supports),
                          sigdigits = sig_digits)
    return new_supports, UniformGrid
end

# IntervalSet and MCSample
function generate_support_values(set::IntervalSet,
                                 ::Val{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, Symbol}
    lb = JuMP.lower_bound(set)
    ub = JuMP.upper_bound(set)
    dist = Distributions.Uniform(lb, ub)
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, MCSample
end

# UniDistributionSet and MultiDistributionSet (with multivariate only)
function generate_support_values(
    set::Union{UniDistributionSet, MultiDistributionSet{<:Distributions.MultivariateDistribution}},
    ::Val{WeightedSample} = Val(WeightedSample);
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{<:Real}, Symbol}
    dist = set.distribution
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, WeightedSample
end

# MultiDistributionSet (matrix-variate distribution)
function generate_support_values(
    set::MultiDistributionSet{<:Distributions.MatrixDistribution},
    ::Val{WeightedSample} = Val(WeightedSample);
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{Float64, 2}, Symbol}
    dist = set.distribution
    raw_supports = Distributions.rand(dist, num_supports)
    new_supports = Array{Float64}(undef, length(dist), num_supports)
    for i in 1:size(new_supports, 2)
        new_supports[:, i] = round.(reduce(vcat, raw_supports[i]),
                                    sigdigits = sig_digits)
    end
    return new_supports, WeightedSample
end

# Generate the supports for a collection set
function _generate_collection_supports(set::CollectionSet, num_supports::Int,
                                       sig_digits::Int)::Array{Float64, 2}
    sets = collection_sets(set)
    # build the support array transpose to fill in column order (leverage locality)
    trans_supports = Array{Float64, 2}(undef, num_supports, length(sets))
    for i in eachindex(sets)
        @inbounds trans_supports[:, i] = generate_support_values(sets[i],
                                                   num_supports = num_supports,
                                                   sig_digits = sig_digits)[1]
    end
    return permutedims(trans_supports)
end

function _generate_collection_supports(set::CollectionSet,
                                       method::Symbol,
                                       num_supports::Int,
                                       sig_digits::Int)::Array{Float64, 2}
    sets = collection_sets(set)
    # build the support array transpose to fill in column order (leverage locality)
    trans_supports = Array{Float64, 2}(undef, num_supports, length(sets))
    for i in eachindex(sets)
        @inbounds trans_supports[:, i] = generate_support_values(sets[i],
                                                   Val(method),
                                                   num_supports = num_supports,
                                                   sig_digits = sig_digits)[1]
    end
    return permutedims(trans_supports)
end

# CollectionSet (IntervalSets)
function generate_support_values(set::CollectionSet{IntervalSet},
                                 ::Val{UniformGrid} = Val(UniformGrid);
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, Symbol}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, UniformGrid
end

function generate_support_values(set::CollectionSet{IntervalSet},
                                 ::Val{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, Symbol}
    new_supports = _generate_collection_supports(set, MCSample, num_supports, sig_digits)
    return new_supports, MCSample
end

# CollectionSet (UniDistributionSets)
function generate_support_values(set::CollectionSet{<:UniDistributionSet},
                                 ::Val{WeightedSample} = Val(WeightedSample);
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, Symbol}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, WeightedSample
end

# CollectionSet (InfiniteScalarSets)
function generate_support_values(set::CollectionSet,
                                 ::Val{Mixture} = Val(Mixture);
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, Symbol}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, Mixture
end

# CollectionSet (InfiniteScalarSets) using purely MC sampling
# this is useful for measure support generation
function generate_support_values(set::CollectionSet,
                                 ::Val{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, Symbol}
    new_supports = _generate_collection_supports(set, MCSample, num_supports, sig_digits)
    return new_supports, MCSample
end
