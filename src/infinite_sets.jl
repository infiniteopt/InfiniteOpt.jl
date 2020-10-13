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

Set and return the lower bound of `set` if such an operation makes sense. Errors
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
    UniqueMeasure{S <: Val{<:Symbol}} <: PublicLabel

A support label for supports that are provided from the `DiscreteMeasureData` 
associated with a measure where a unique label is generated to distinguish those 
supports. This is done by invoking [`generate_unique_label`](@ref).
"""
struct UniqueMeasure{S} <: PublicLabel end

"""
    InternalLabel <: AbstractSupportLabel

An abstract type for support labels that are associated with supports that should 
not be reported to the user by default.
"""
abstract type InternalLabel <: AbstractSupportLabel end 

"""
    OrthogonalCollocationNode <: AbstractSupportLabel

An support label for additional supports generated by the orthogonal collocation 
method (i.e. nodes between public supports of an infinite parameter).
"""
struct OrthogonalCollocationNode <: InternalLabel end

"""
    generate_unique_label()::Type{UniqueMeasure}

Generate and return a unique support label for measures.
"""
function generate_unique_label()::DataType
    return UniqueMeasure{Val{gensym()}}
end

# Define default values of sig_digits and num_supports keywords
const DefaultSigDigits = 12
const DefaultNumSupports = 10

# a user interface of generate_support_values
"""
    generate_supports(set::AbstractInfiniteSet
                      [method::Type{<:AbstractSupportLabel}];
                      [num_supports::Int = DefaultNumSupports,
                      sig_digits::Int = DefaultSigDigits]
                      )::Tuple{Array{<:Real}, DataType}

Generate `num_supports` support values with `sig_digits` significant digits in
accordance with `set` and return them along with the correct generation label(s).
`IntervalSet`s generate supports uniformly with label `UniformGrid` and
distribution sets generate them randomly accordingly to the
underlying distribution. Moreover, `method` indicates the generation method that
should be used. These `methods` correspond to parameter support labels. Current
labels that can be used as generation methods include (but may not be defined
for certain set types):
- [`MCSample`](@ref): Uniformly distributed Monte Carlo samples.
- [`WeightedSample`](@ref): Monte Carlo samples that are weighted by an underlying PDF.
- [`UniformGrid`](@ref): Samples that are generated uniformly over the set domain.

Extensions that employ user-defined infinite set types and/or methods
should extend [`generate_support_values`](@ref) to enable this. Errors if the
`set` type and /or methods are unrecognized. This is intended as an internal
method to be used by methods such as [`generate_and_add_supports!`](@ref).
"""
function generate_supports(set::AbstractInfiniteSet;
                           num_supports::Int = DefaultNumSupports,
                           sig_digits::Int = DefaultSigDigits
                           )::Tuple
    return generate_support_values(set, num_supports = num_supports,
                                   sig_digits = sig_digits)
end

# 2 arguments
function generate_supports(set::AbstractInfiniteSet,
                           method::Type{<:AbstractSupportLabel};
                           num_supports::Int = DefaultNumSupports,
                           sig_digits::Int = DefaultSigDigits
                           )::Tuple
    return generate_support_values(set, method,
                                   num_supports = num_supports,
                                   sig_digits = sig_digits)
end

"""
    generate_support_values(set::AbstractInfiniteSet,
                            [method::Type{MyMethod} = MyMethod];
                            [num_supports::Int = DefaultNumSupports,
                            sig_digits::Int = DefaultSigDigits]
                            )::Tuple{Array{<:Real}, Symbol}

A multiple dispatch method for [`generate_supports`](@ref). This will return
a tuple where the first element are the supports and the second is their
label. This can be extended for user-defined infinite sets and/or generation
methods. When defining a new set type the default method dispatch should
make `method` an optional argument (making it the default). Otherwise, other
method dispatches for a given set must ensure that `method` is positional
argument without a default value (contrary to the definition above). Note that the 
`method` must be a subtype of either [`PublicLabel`](@ref) or [`InternalLabel`](@ref).
"""
function generate_support_values(set::AbstractInfiniteSet,
                                 args...; kwargs...)
    if isempty(args)
        error("`generate_support_values` has not been extended for infinite sets " * 
              "of type `$(typeof(set))`. This automatic support generation is not " * 
              "implemented.")
    else
        error("`generate_support_values` has not been extended for infinite sets " * 
              "of type `$(typeof(set))` with the generation method `$(args[1])`. " * 
              "This automatic support generation is not implemented.")
    end
end

# IntervalSet and UniformGrid
function generate_support_values(set::IntervalSet,
                                 method::Type{UniformGrid} = UniformGrid;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, DataType}
    lb = JuMP.lower_bound(set)
    ub = JuMP.upper_bound(set)
    new_supports = round.(range(lb, stop = ub, length = num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# IntervalSet and MCSample
function generate_support_values(set::IntervalSet,
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits,
                                 )::Tuple{Vector{<:Real}, DataType}
    lb = JuMP.lower_bound(set)
    ub = JuMP.upper_bound(set)
    dist = Distributions.Uniform(lb, ub)
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# UniDistributionSet and MultiDistributionSet (with multivariate only)
function generate_support_values(
    set::Union{UniDistributionSet, MultiDistributionSet{<:Distributions.MultivariateDistribution}},
    method::Type{WeightedSample} = WeightedSample;
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{<:Real}, DataType}
    dist = set.distribution
    new_supports = round.(Distributions.rand(dist, num_supports),
                          sigdigits = sig_digits)
    return new_supports, method
end

# UniDistributionSet and MCSample 
function generate_support_values(
    set::UniDistributionSet,
    method::Type{MCSample};
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Vector{Float64}, DataType}
    return generate_support_values(set, WeightedSample; num_supports = num_supports, 
                                   sig_digits = sig_digits)[1], method # TODO use an unwieghted sample...
end

# MultiDistributionSet (matrix-variate distribution)
function generate_support_values(
    set::MultiDistributionSet{<:Distributions.MatrixDistribution},
    method::Type{WeightedSample} = WeightedSample;
    num_supports::Int = DefaultNumSupports,
    sig_digits::Int = DefaultSigDigits
    )::Tuple{Array{Float64, 2}, DataType}
    dist = set.distribution
    raw_supports = Distributions.rand(dist, num_supports)
    new_supports = Array{Float64}(undef, length(dist), num_supports)
    for i in 1:size(new_supports, 2)
        new_supports[:, i] = round.(reduce(vcat, raw_supports[i]),
                                    sigdigits = sig_digits)
    end
    return new_supports, method
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
                                       method::Type{<:AbstractSupportLabel},
                                       num_supports::Int,
                                       sig_digits::Int)::Array{Float64, 2}
    sets = collection_sets(set)
    # build the support array transpose to fill in column order (leverage locality)
    trans_supports = Array{Float64, 2}(undef, num_supports, length(sets))
    for i in eachindex(sets)
        @inbounds trans_supports[:, i] = generate_support_values(sets[i],
                                                   method,
                                                   num_supports = num_supports,
                                                   sig_digits = sig_digits)[1]
    end
    return permutedims(trans_supports)
end

# CollectionSet (IntervalSets)
function generate_support_values(set::CollectionSet{IntervalSet},
                                 method::Type{UniformGrid} = UniformGrid;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, method
end

function generate_support_values(set::CollectionSet{IntervalSet},
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(set, method, num_supports, sig_digits)
    return new_supports, method
end

# CollectionSet (UniDistributionSets)
function generate_support_values(set::CollectionSet{<:UniDistributionSet},
                                 method::Type{WeightedSample} = WeightedSample;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, method
end

# CollectionSet (InfiniteScalarSets)
function generate_support_values(set::CollectionSet,
                                 method::Type{Mixture} = Mixture;
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(set, num_supports, sig_digits)
    return new_supports, method
end

# CollectionSet (InfiniteScalarSets) using purely MC sampling
# this is useful for measure support generation
function generate_support_values(set::CollectionSet,
                                 method::Type{MCSample};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits
                                 )::Tuple{Array{<:Real}, DataType}
    new_supports = _generate_collection_supports(set, method, num_supports, sig_digits)
    return new_supports, method
end

# For label All: dispatch to default methods
function generate_support_values(set::AbstractInfiniteSet, ::Type{All};
                                 num_supports::Int = DefaultNumSupports,
                                 sig_digits::Int = DefaultSigDigits)
    return generate_support_values(set, num_supports = num_supports,
                                   sig_digits = sig_digits)
end
