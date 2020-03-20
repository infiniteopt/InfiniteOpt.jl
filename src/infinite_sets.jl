"""
    supports_in_set(supports::Union{Number, Vector{<:Number}},
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
function supports_in_set(supports::Union{Number, Vector{<:Number}},
                         set::IntervalSet)::Bool
    min_support = minimum(supports)
    max_support = maximum(supports)
    if min_support < set.lower_bound || max_support > set.upper_bound
        return false
    end
    return true
end

# UnivariateDistribution set
function supports_in_set(supports::Union{Number, Vector{<:Number}},
                         set::DistributionSet{<:Distributions.UnivariateDistribution})::Bool
    min_support = minimum(supports)
    max_support = maximum(supports)
    check1 = min_support < minimum(set.distribution)
    check2 = max_support > maximum(set.distribution)
    if check1 || check2
        return false
    end
    return true
end

# Fallback
function supports_in_set(supports::Union{Number, Vector{<:Number}},
                         set::AbstractInfiniteSet)::Bool
    return true
end

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
JuMP.has_lower_bound(set::DistributionSet{<:Distributions.UnivariateDistribution})::Bool = true

"""
    JuMP.lower_bound(set::AbstractInfiniteSet)::Number

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
JuMP.lower_bound(set::IntervalSet)::Number = set.lower_bound

# DistributionSet (Univariate)
function JuMP.lower_bound(set::DistributionSet{<:Distributions.UnivariateDistribution})::Number
    return minimum(set.distribution)
end

"""
    JuMP.set_lower_bound(set::AbstractInfiniteSet,
                         lower::Number)::AbstractInfiniteSet

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
function JuMP.set_lower_bound(set::AbstractInfiniteSet, lower::Number) # fallback
    type = typeof(set)
    error("`JuMP.set_lower_bound` not defined for infinite set of type $type.")
end

# IntervalSet
function JuMP.set_lower_bound(set::IntervalSet, lower::Number)::IntervalSet
    return IntervalSet(lower, set.upper_bound)
end

# DistributionSet
function JuMP.set_lower_bound(set::DistributionSet, lower::Number)
    error("Cannot set the lower bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

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
JuMP.has_upper_bound(set::DistributionSet{<:Distributions.UnivariateDistribution})::Bool = true

"""
    JuMP.upper_bound(set::AbstractInfiniteSet)::Number

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
JuMP.upper_bound(set::IntervalSet)::Number = set.upper_bound

# DistributionSet (Univariate)
function JuMP.upper_bound(set::DistributionSet{<:Distributions.UnivariateDistribution})::Number
    return maximum(set.distribution)
end

"""
    JuMP.set_upper_bound(set::AbstractInfiniteSet,
                         upper::Number)::AbstractInfiniteSet

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
function JuMP.set_upper_bound(set::AbstractInfiniteSet, upper::Number) # fallback
    type = typeof(set)
    error("`JuMP.set_lower_bound` not defined for infinite set of type $type.")
end

# IntervalSet
function JuMP.set_upper_bound(set::IntervalSet, upper::Number)
    return IntervalSet(set.lower_bound, upper)
end

# DistributionSet
function JuMP.set_upper_bound(set::DistributionSet, lower::Number)
    error("Cannot set the upper bound of a distribution, try using " *
          "`Distributions.Truncated` instead.")
end

"""
    generate_support_values(set::AbstractInfiniteSet; [num_supports::Int = 10,
                            sig_fig::Int = 5])

Generate `num_supports` support values with `sig_figs` significant digits in
accordance with `set` and return them. `IntervalSet`s generate supports
uniformly and `DistributionSet`s generate them randomly accordingly to the
underlyign distribution. Extensions that employ user-defined infinite set types
should extend according to the new set type. Errors if the `set` is a type that
has not been explicitly extended. This is intended as an internal method to be
used by [`generate_and_add_supports!`](@ref) and [`build_parameter`](@ref).
"""
function generate_support_values(set::AbstractInfiniteSet; num_supports::Int = 0,
                                 sig_fig::Int = 1)
    type = typeof(set)
    error("Unable to generate support values for unrecognized infinite set " *
          "type $type")
end

# IntervalSet
function generate_support_values(set::IntervalSet; num_supports::Int = 10,
                         sig_fig::Int = 5)::Vector
    lb = set.lower_bound
    ub = set.upper_bound
    new_supports = round.(collect(range(lb, stop = ub, length = num_supports)),
                          sigdigits = sig_fig)
    return new_supports
end

# DistributionSet
function generate_support_values(set::DistributionSet; num_supports::Int = 10,
                         sig_fig::Int = 5)::Array
    new_supports = round.(Distributions.rand(set.distribution, num_supports),
                          sigdigits = sig_fig)
    return new_supports
end
