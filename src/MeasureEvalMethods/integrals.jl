################################################################################
#                          EVALUATION METHOD DATATYPES
################################################################################
"""
    AbstractIntegralMethod

An abstract type for integral evaluation methods use in combination with
`integral` and `generate_integral_data`.
"""
abstract type AbstractIntegralMethod end

"""
    Automatic <: AbstractIntegralMethod

An integral evaluation type for automically selecting an appropriate integral
evaluation method.
"""
struct Automatic <: AbstractIntegralMethod end

"""
    AbstractUnivariateMethod <: AbstractIntegralMethod

An abstract type for integral evaluation methods for 1-dimensional integrals.
"""
abstract type AbstractUnivariateMethod <: AbstractIntegralMethod end

"""
    UniTrapezoid <: AbstractUnivariateMethod

An integral evalution method that uses the trapezoid rule to in combination
with all parameter supports available when the integral is expanded and/or when
the infinite model is optimized, whichever comes first. Note this method will
ignore the `num_supports` keyword argument. Note this is valid only for finite
integral domains.
"""
struct UniTrapezoid <: AbstractUnivariateMethod end

"""
    UniMCSampling <: AbstractUnivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral. This variant will add more supports to the model as
needed to satisfy `num_supports` and it will include all supports with the
`MCSample` label up till the integral is expanded and/or when
the infinite model is optimized, whichever comes first. Note this is valid only
for finite integral domains.
"""
struct UniMCSampling <: AbstractUnivariateMethod end

"""
    UniIndepMCSampling <: AbstractUnivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral similar to [`UniMCSampling`](@ref MeasureEvalMethods.UniMCSampling).
However, this variant will generate its own set of supports and ignore all other
supports with the `MCSample` label. Note this is valid only for finite integral
domains.
"""
struct UniIndepMCSampling <: AbstractUnivariateMethod end

"""
    Quadrature <: AbstractUnivariateMethod

A general integral evaluation method that will automatically select the
appropriate quadrature method to approximate the integral. Please note that this
will generate a unique set of parameter supports and will ignore existing supports
when the integral is evaluated and thus should be used with caution. However,
this method is able to handle infinite and semi-infinite integral domains.
"""
struct Quadrature <: AbstractUnivariateMethod end

"""
    GaussHermite <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Hermite quadrature to
evaluate integrals. This is valid for infinite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
"""
struct GaussHermite <: AbstractUnivariateMethod end

"""
    GaussLegendre <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Legendre quadrature to
evaluate integrals. This is valid for finite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
"""
struct GaussLegendre <: AbstractUnivariateMethod end

"""
    GaussLaguerre <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Laguerre quadrature to
evaluate integrals. This is valid for semi-infinite integral domains. Note this
will generate its own set of supports and will ignore other parameter supports.
"""
struct GaussLaguerre <: AbstractUnivariateMethod end

"""
    AbstractMultivariateMethod <: AbstractIntegralMethod

An abstract type for integral evaluation methods for multi-dimensional integrals.
"""
abstract type AbstractMultivariateMethod <: AbstractIntegralMethod end

"""
    MultiMCSampling <: AbstractMultivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral. This variant will add more supports to the model as
needed to satisfy `num_supports` and it will include all supports with the
`MCSample` label up till the integral is expanded and/or when
the infinite model is optimized, whichever comes first. Note this is valid only
for finite integral domains.
"""
struct MultiMCSampling <: AbstractMultivariateMethod end

"""
    MultiIndepMCSampling <: AbstractMultivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral similar to [`MultiMCSampling`](@ref MeasureEvalMethods.MultiMCSampling).
However, this variant will generate its own set of supports and ignore all other
supports with the `MCSample` label. Note this is valid only for finite integral
domains.
"""
struct MultiIndepMCSampling <: AbstractMultivariateMethod end

################################################################################
#                       DATA GENERATION METHODS (UNIVARIATE)
################################################################################
# TODO update these as needed for below
function _truncated_set(set::InfiniteOpt.IntervalSet,
                        lb::Real, ub::Real)::InfiniteOpt.IntervalSet
    new_lb = maximum(set.lower_bound, lb)
    new_ub = minimum(set.upper_bound, ub)
    return InfiniteOpt.IntervalSet(new_lb, new_ub)
end

"""
    generate_integral_data(
        prefs::Union{InfiniteOpt.GeneralVariableRef, Vector{InfiniteOpt.GeneralVariableRef}},
        lower_bounds::Union{Real, Vector{<:Real}},
        upper_bounds::Union{Real, Vector{<:Real}},
        method::Type{V}; [num_supports::Int = InfiniteOpt.DefaultNumSupports,
        weight_func::Function = InfiniteOpt.default_weight,
        extra_kwargs...]
        )::InfiniteOpt.AbstractMeasureData where {V <: AbstractIntegralMethod}

Generate the appropriate concrete realization of `AbstractMeasureData` using
`method`. Here `prefs`, `lower_bounds`, and `upper_bounds` will always have a
1 to 1 correspondence when this is called from `integral`. Please refer
to the method docstrings for an explanation of each one.

User-defined method extensions should first
define a concrete `method` type inheriting from `AbstractUnivariateMethod` or
`AbstractMultivariateMethod` as appropriate and then implement extend this
method using that type for `method`.
"""
function generate_integral_data end

# General fallback
function generate_integral_data(prefs, lb, ub, method; kwargs...)
    error("`generate_integral_data` is not defined for method type `$(method)` " *
          "in accordance with the arugments given.")
end

# Single pref with Automatic
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
    lower_bound::Real,
    upper_bound::Real,
    method::Type{Automatic};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )::InfiniteOpt.AbstractMeasureData
    # TODO finish this is the forgiving method of users
    # If finite should dispatch to trapezoid and ignore num_supports
    # if infinite or semi-infinite shold use appropriate quadrature (error if pref is dependent)
    if InfiniteOpt._index_type(pref) == DependentParameterIndex
        error("Univariate integral data cannot be generated for single dependent parameter.")
    end
    inf_bound_num = (lower_bound == -Inf) && (upper_bound == Inf)
    if inf_bound_num == 0 # finite interval
        method = UniTrapezoid
    elseif inf_bound_num == 1 # semi-infinite interval
        method = GaussLaguerre
    else # infinite interval
        method = GaussHermite
    end
    return generate_integral_data(pref, lower_bound, upper_bound, method,
                                  num_supports = num_supports,
                                  weight_func = weight_func)
end

function _trapezoid_coeff(supps::AbstractArray{<:Real, 1})::Vector{<:Real}
    len = length(supps)
    coeffs = zero(length)
    for i in eachindex(supps)
        if i == 1
            coeffs[i] = 1 / 2 * (supps[2] - supps[1])
        elseif i == len
            coeffs[i] = 1 / 2 * (supps[len] - supps[len-1])
        else
            coeffs[i] = 1 / 2 * (supps[i+1] - supps[i-1])
        end
    end
    return coeffs
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{UniTrapezoid};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    return FunctionalDiscreteMeasureData(pref, _trapezoid_coeff, num_supports,
                                         some_label, weight_func, lower_bound,
                                         upper_bound, false)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{Quadrature};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)

    inf_bound_num = (lb == -Inf) + (ub == Inf)
    if inf_bound_num == 0 # finite interval
        method = GaussLegendre
    elseif inf_bound_num == 1 # semi-infinite interval
        method = GaussLaguerre
    else # infinite interval
        method = GaussHermite
    end
    return generate_integral_data(pref, lower_bound, upper_bound, method,
                                  num_supports = num_supports, weight_func = weight_func)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{GaussLegendre};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    if lower_bound == -Inf || upper_bound == Inf
        error("Gauss Legendre quadrature can only be applied on finite intervals.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (upper_bound - lower_bound) / 2 * supports .+ (upper_bound + lower_bound) / 2
    coeffs = (upper_bound - lower_bound) / 2 * coeffs
    return DiscreteMeasureData(pref, coeffs, supports, gensym(), weight_func,
                               lower_bound, upper_bound, false)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{GaussLaguerre};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    if upper_bound == Inf
        if lower_bound == -Inf
            error("The range is infinite. Use other measure evaluation methods.")
        end
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = supports .+ lower_bound
    elseif lb == -Inf
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = -supports .+ upper_bound
    else
        error("The interval is finite. Use other measure evaluation methods.")
    end

    return DiscreteMeasureData(pref, coeffs, supports, gensym(), weight_func,
                               lower_bound, upper_bound, false)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{GaussHermite};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    if lower_bound != -Inf || upper_bound != Inf
        error("The range is not infinite. Use other measure evaluation methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return DiscreteMeasureData(pref, coeffs, supports, gensym(),
                               weight_func, lower_bound, upper_bound, false)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{UniMCSampling};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    if lower_bound == -Inf || upper_bound == Inf
        error("Univariate MC sampling is not applicable to (semi-)infinite intervals.")
    end

    new_lb = maximum(lower_bound(pref), lower_bound)
    new_ub = minimum(upper_bound(pref), upper_bound)

    supports = supports(pref, label = MCSample)
    filter!(x->(x >= new_lb && x <= new_ub), supports)

    num_samples = length(supports)
    if num_samples < num_supports
        new_supports, _ = generate_support_values(new_set, Val(MCSample),
                                                      num_supports = num_supports - num_samples)
        add_supports(pref, new_supports, label = MCSample) # this ensures new supports have MCSample label
    end
    function coeff_func(supps::AbstractArray{<:Real, 1})::Vector{Float64}
        return ones(num_supports) / num_supports * (new_ub - new_lb)
    end
    return FunctionalDiscreteMeasureData(pref, coeff_func, num_supports, MCSample,
                                         lower_bound, upper_bound, false)
end

function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Type{UniIndepMCSampling};
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight)
    if lower_bound == -Inf || upper_bound == Inf
        error("Univariate MC sampling is not applicable to (semi-)infinite intervals.")
    end

    new_lb = maximum(lower_bound(pref), lower_bound)
    new_ub = minimum(upper_bound(pref), upper_bound)

    supports, _ = generate_support_values(new_set, Val(MCSample),
                                          num_supports = num_supports - num_samples)
    coeffs = ones(num_supports) / num_supports * (new_ub - new_lb)
    return DiscreteMeasureData(pref, coeffs, supports, gensym(),
                               weight_func, lower_bound, upper_bound, false)
end

################################################################################
#                       DATA GENERATION METHODS (MULTIVARIATE)
################################################################################

# Vector prefs with Automatic
function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real},
    method::Type{Automatic};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )::InfiniteOpt.AbstractMeasureData
    # TODO finish this is the forgiving method of users
    # this should dispatch to use MultiMCSampling
    return generate_integral_data(prefs, lower_bound, upper_bound, MultiMCSampling,
                                  num_supports = num_supports,
                                  weight_func = weight_func)
end


function _truncate_collection_set_with_vol(
    prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real})::Tuple{CollectionSet{IntervalSet}, Float64}
    set = IntervalSet[]
    vol = 1.
    for i in eachindex(prefs)
        new_lb = maximum(lower_bound(prefs[i]), lower_bound[i])
        new_ub = minimum(upper_bound(prefs[i]), upper_bound[i])
        push!(set, IntervalSet(new_lb[i], new_ub[i]))
        vol *= (new_ub - new_lb)
    end
    set = CollectionSet(set)
    return (set, vol)
end

function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real},
    method::Type{MultiMCSampling};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )::InfiniteOpt.AbstractMeasureData

    for i in eachindex(lower_bound)
        if lower_bound[i] == -Inf || upper_bound[i] == Inf
            error("MC Sampling is not applicable to (semi-)infinite intervals.")
        end
    end

    set, vol = _truncate_collection_set_with_vol(prefs, lower_bound, upper_bound)
    supports = supports(prefs, label = MCSample)
    valid_supports_idx = [supports_in_set(i, set) for i in eachcol(supports)]
    supports = supports[valid_supports_idx]

    num_samples = size(supports)[2]
    if num_samples < num_supports
        new_supports, _ = generate_support_values(set, Val(MCSample),
                                                  num_supports = num_supports - num_samples)
        add_supports(pref, new_supports, label = MCSample) # this ensures new supports have MCSample label
    end
    function coeff_func(supps::AbstractArray{<:Real, 2})::Vector{Float64}
        return ones(num_samples) / num_samples * vol
    end
    return FunctionalDiscreteMeasureData(pref, coeff_func, num_supports, MCSample,
                                         lower_bound, upper_bound, false)
end

function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real},
    method::Type{MultiIndepMCSampling};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )::InfiniteOpt.AbstractMeasureData

    for i in eachindex(lower_bound)
        if lower_bound[i] == -Inf || upper_bound[i] == Inf
            error("MC Sampling is not applicable to (semi-)infinite intervals.")
        end
    end
    set, vol = _truncate_collection_set_with_vol(prefs, lower_bound, upper_bound)
    supports, _ = generate_support_values(set, Val(MCSample), num_supports = num_supports)
    coeffs = 1 / vol * num_supports
    return DiscreteMeasureData(prefs, coeffs, supports, gensym(), weight_func,
                               lower_bound, upper_bound, false)
end

################################################################################
#                               HELPER METHODS
################################################################################
#=
"""
    infinite_transform(lb::Number, ub::Number, num_supports::Int;
                       [sub_method::Function = mc_sampling,
                       transform_x::Function = _default_x,
                       transform_dx::Function = _default_dx,
                       t_lb::Number = -convert(Number, lb == -Inf && ub == Inf),
                       t_ub::Number = 1., kwargs...])::Tuple

Returns a tuple that contains supports and coefficients generated for a
parameter in an infinite or semi-infinite interval. It works by transforming
the original unbounded interval to a finite interval, on which a support
generation method for finite intervals is applied. Then, the generated supports
are transformed back to the original interval. The user is allowed to specify
the support generation method for finite intevals to use, as well as the
transform function. The default transform function is
``t \\in [-\\infty, \\infty] \\rightarrow x \\in [-1, 1]: t(x) = \\frac{t}{1-t^2}``
``t \\in [a, \\infty] \\rightarrow x \\in [0, 1]: t(x) = a + \\frac{t}{1-t}``
``t \\in [-\\infty, a] \\rightarrow x \\in [0, 1]: t(x) = a - \\frac{1-t}{t}``

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> (supps, coeffs) = infinite_transform(-Inf, Inf, 5, sub_method = gauss_legendre)
([-5.06704059565454, -0.7583532171678754, 0.0, 0.7583532171678754, 5.06704059565454], [13.490960583398396, 1.2245949721571516, 0.5688888888888889, 1.2245949721571516, 13.490960583398396])
```
"""
function infinite_transform(set::InfiniteOpt.IntervalSet,
                            params::Union{InfiniteOpt.ParameterRef,
                            AbstractArray{<:InfiniteOpt.ParameterRef}},
                            num_supports::Int,
                            lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                            ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                            sub_method::Val = Val(mc_sampling);
                            transform_x::Function = _default_x,
                            transform_dx::Function = _default_dx,
                            t_lb::Number = -convert(Number, lb == -Inf && ub == Inf),
                            t_ub::Number = 1.)::Tuple
    # transform (semi-)infinite domain to finite domain
    if lb != -Inf && ub != Inf
        error("The range is not (semi-)infinite. Use evaluation methods for " *
              "bounded domains.")
    end
    (t_supports, t_coeffs) = generate_integral_data(set, params, num_supports, t_lb, t_ub, sub_method)
    supports = transform_x.(t_supports, lb, ub)
    coeffs = t_coeffs .* transform_dx.(t_supports, lb, ub)
    return (supports, coeffs)
end

function _default_x(t::Number, lb::Number, ub::Number)::Number
    if lb > -Inf
        return lb + t / (1 - t)
    elseif ub < Inf
        return ub - (1 - t) / t
    else
        return t / (1 - t^2)
    end
end

function _default_dx(t::Number, lb::Number, ub::Number)::Number
    if lb > -Inf
        return 1 / (1 - t)^2
    elseif ub < Inf
        return 1 / t^2
    else
        return (1 + t^2) / (1 - t^2)^2
    end
end
=#
# TODO: consider truncated distribution
# TODO: consider adding uniform grids
