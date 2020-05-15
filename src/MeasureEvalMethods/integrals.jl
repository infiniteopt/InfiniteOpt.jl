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
#                             DATA GENERATION METHODS
################################################################################
# TODO update these as needed for below
function _truncated_set(set::InfiniteOpt.IntervalSet,
                        lb::Union{Real, Nothing},
                        ub::Union{Real, Nothing})::InfiniteOpt.IntervalSet
    isa(lb, Nothing) ? new_lb = set.lower_bound : new_lb = lb
    isa(ub, Nothing) ? new_ub = set.lower_bound : new_lb = ub
    return InfiniteOpt.IntervalSet(new_lb, new_ub)
end

function _truncated_set(set::InfiniteOpt.UniDistributionSet,
                        lb::Union{Real, Nothing},
                        ub::Union{Real, Nothing})::InfiniteOpt.UniDistributionSet
    isa(lb, Nothing) ? new_lb = -Inf : new_lb = lb
    isa(ub, Nothing) ? new_ub = Inf : new_lb = ub
    if new_lb == -Inf && new_ub == Inf
        return set
    else
        return InfiniteOpt.UniDistributionSet(Distributions.truncated(dist, lb, ub))
    end
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
    return
end

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
    return
end

#=
# MC Sampling methods
# IntervalSets & UniDistributionSet
function generate_integral_data(set::InfiniteScalarSet,
                                      param::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Real, Nothing},
                                      ub::Union{Real, Nothing},
                                      method::Val{mc_sampling})::Tuple
    _nothing_test(method, lb, ub)
    # use the correct set for generating supports
    new_set = _truncated_set(set, lb, ub)
    supports = supports(param, label = MCSample)
    num_samples = length(supports)
    if num_samples >= num_supports
        coeffs = ones(num_samples) / num_samples * (new_ub - new_lb)
    else
        new_supports, _ = generate_support_values(new_set, Val(MCSample),
                                                      num_supports = num_supports - num_samples)
        union!(supports, new_supports)
        coeffs = ones(num_supports) / num_supports * (new_ub - new_lb)
    end
    add_supports(param, supports, label = MCSample) # this ensures new supports have MCSample label
    return (supports, coeffs, label)
end

# MultiDistributionSet
function generate_integral_data(set::InfiniteOpt.MultiDistributionSet,
                                      params::AbstractArray{<:InfiniteOpt.GeneralVariableRef},
                                      num_supports::Int,
                                      lb::Union{AbstractArray{<:Real}, Nothing},
                                      ub::Union{AbstractArray{<:Real}, Nothing},
                                      method::Val{mc_sampling})::Tuple
    # use the correct set for generating supports
    if !isa(lb, Nothing) || !isa(ub, Nothing)
        @warn("MC sampling for truncated multivariate distribution is not supported. " *
              "Lower and upper bounds are thus ignored.")
    end
    supports = supports(params, label = MCSample)
    num_samples = size(supports)[2]
    if num_samples >= num_supports
        coeffs = ones(num_samples) / num_samples
    else
        new_supports, _ = generate_support_values(new_set, Val(MCSample),
                                                  num_supports = num_supports - num_samples)
        union!(supports, new_supports)
        coeffs = ones(num_supports) / num_supports
    end
    add_supports(params, supports, label = MCSample) # this ensures new supports have MCSample label
    return (supports, coeffs, MCSample)
end

# CollectionSet
function generate_integral_data(set::InfiniteOpt.CollectionSet,
                                      params::AbstractArray{<:InfiniteOpt.GeneralVariableRef},
                                      num_supports::Int,
                                      lb::Union{AbstractArray{<:Real}, Nothing},
                                      ub::Union{AbstractArray{<:Real}, Nothing},
                                      method::Val{mc_sampling})::Tuple
    # Truncate each dimension accordingly
    if !isa(lb, Nothing) || !isa(ub, Nothing)
        if isa(lb, Nothing)
            lb = -Inf * ones(length(params))
        end
        if isa(ub, Nothing)
            ub = Inf * ones(length(params))
        end
        new_set = CollectionSet([_truncated_set(set.sets[i], lb[i], ub[i]) for i in eachindex(set.sets)])
    else
        new_set = set
    end
    supports = supports(params, label = MCSample)
    num_samples = size(supports)[2]
    if num_samples >= num_supports
        coeffs = ones(num_samples) / num_samples
    else
        new_supports, _ = generate_support_values(set, Val(MCSample),
                                                  num_supports = num_supports - num_samples)
        supports = hcat(supports, new_supports)
        coeffs = ones(curr_num_supps) / num_supports
    end
    add_supports(params, supports, label = MCSample) # this ensures new supports have MCSample label
    return (supports, coeffs, MCSample)
end

# TODO incorporate the below into the appropriate functions
# if eval_method == quadrature || eval_method == gauss_legendre ||
#    eval_method == gauss_hermite || eval_method == gauss_laguerre
#     if num_params > 1
#         error("Quadrature method is not supported for multivariate measures.")
#     end
#     inf_bound_num = (lb == -Inf) + (ub == Inf)
#     if inf_bound_num == 0
#         if eval_method != quadrature || eval_method != gauss_legendre
#             @warn("$(eval_method) is not compatible with finite intervals. " *
#                   "The integral will use Gauss-Legendre quadrature instead.")
#         end
#         kwargs[:eval_method] = gauss_legendre
#     elseif inf_bound_num == 1
#         if eval_method != quadrature || eval_method != gauss_laguerre
#             @warn("$(eval_method) is not compatible with semi-infinite intervals. " *
#                   "The integral will use Gauss-Laguerre quadrature instead.")
#         end
#         kwargs[:eval_method] = gauss_laguerre
#     else
#         if eval_method != quadrature || eval_method != gauss_hermite
#             @warn("$(eval_method) is not compatible with infinite intervals. " *
#                   "The integral will use Gauss-Hermite quadrature instead.")
#         end
#         kwargs[:eval_method] = gauss_hermite
#     end
# end

# Gaussian Quadrature Methods for IntervalSet
# Gauss-Legendre for finite intervals
function generate_integral_data(set::InfiniteScalarSet,
                                      param::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Real, Nothing},
                                      ub::Union{Real, Nothing},
                                      method::Val{gauss_legendre})::Tuple
    _nothing_test(method, lb, ub)
    if lb == -Inf || ub == Inf
        error("Gauss Legendre quadrature can only be applied on finite intervals.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs, gauss_legendre)
end


function generate_integral_data(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Real, Nothing},
                                      ub::Union{Real, Nothing},
                                      method::Val{gauss_hermite})::Tuple
    _nothing_test(method, lb, ub)
    if lb != -Inf || ub != Inf
        error("The range is not finite. Use other measure evaluation methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs, gauss_hermite)
end

function generate_integral_data(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{gauss_laguerre})::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    # default Gaussian quadrature method for semi-infinite domain
    if ub == Inf
        if lb == -Inf
            error("The range is infinite. Use other measure evaluation methods.")
        end
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = supports .+ lb
        return (supports, coeffs, gauss_laguerre)
    elseif lb == -Inf
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = -supports .+ ub
        return (supports, coeffs, gauss_laguerre)
    else
        error("The range is bounded. Use other measure evaluation methods.")
    end
end


function generate_integral_data(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{trapezoid})::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    if trapz
    supports, coeffs = _trapezoid_univ(lb, ub, num_supports)
    return (supports, coeffs, trapezoid)
end

function _trapezoid_univ(lb::Number, ub::Number, num_supports::Int)::Tuple
    increment = (ub - lb) / (num_supports - 1)
    supps = [(i - 1) * increment + lb for i in 1:num_supports]
    coeffs = increment * ones(num_supports)
    coeffs[1] = coeffs[1] / 2
    coeffs[num_supports] = coeffs[num_supports] / 2
    return (supps, coeffs)
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

function generate_integral_data(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.GeneralVariableRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{trapezoid};
                                      is_functional = false)::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    if is_functional
        supports = []
        coeffs = _trapezoid_coeff
    else
        supports, coeffs = _trapezoid_univ(lb, ub, num_supports)
    end
    return (supports, coeffs, trapezoid)
end

function _nothing_test(method::Val,
                       lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                       ub::Union{Number, JuMPC.SparseAxisArray, Nothing})
    if isa(lb, Nothing) || isa(ub, Nothing)
        error("Method " * string(method)[6:length(string(method))-3] *
              " cannot be nothing.")
    end
end

# test if the lower/upper bound is univariate
function _univariate_bounds(method::Val,
                            lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                            ub::Union{Number, JuMPC.SparseAxisArray, Nothing})
    if isa(lb, JuMPC.SparseAxisArray) || isa(ub, JuMPC.SparseAxisArray)
        error("Method " * string(method)[6:length(string(method))-3] *
              " is not applicable to multivariate infinite parameters.")
    end
end

################################################################################
#                               HELPER METHODS
################################################################################
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
