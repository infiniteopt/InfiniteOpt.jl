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
evaluation method. Contains no fields.
"""
struct Automatic <: AbstractIntegralMethod end

"""
    AbstractUnivariateMethod <: AbstractIntegralMethod

An abstract type for integral evaluation methods for 1-dimensional integrals.
"""
abstract type AbstractUnivariateMethod <: AbstractIntegralMethod end

"""
<<<<<<< HEAD
Abstract type for Finite 1-D integral evaluation using a variety of Gaussian Quadrature methods over the entire domain.
"""
abstract type FiniteGaussQuad <: AbstractUnivariateMethod end

"""
=======
>>>>>>> f36e98edf94dd41f108f13ce112e07afc3151f70
    UniTrapezoid <: AbstractUnivariateMethod

An integral evalution method that uses the trapezoid rule to in combination
with all parameter supports available when the integral is expanded and/or when
the infinite model is optimized, whichever comes first. Note this method will
ignore the `num_supports` keyword argument. Note this is valid only for finite
integral domains. Contains no fields.
"""
struct UniTrapezoid <: AbstractUnivariateMethod end

"""
    UniMCSampling <: AbstractUnivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral. This variant will add more supports to the model as
needed to satisfy `num_supports` and it will include all supports with the
`MCSample` label up till the integral is expanded and/or when
the infinite model is optimized, whichever comes first. Note this is valid only
for finite integral domains. Contains no fields.
"""
struct UniMCSampling <: AbstractUnivariateMethod end

"""
    UniIndepMCSampling <: AbstractUnivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral similar to [`UniMCSampling`](@ref MeasureToolbox.UniMCSampling).
However, this variant will generate its own set of supports and ignore all other
supports with the `MCSample` label. Note this is valid only for finite integral
domains. This is not compatible with individual dependent parameters. Contains no fields.
"""
struct UniIndepMCSampling <: AbstractUnivariateMethod end

"""
    Quadrature <: AbstractUnivariateMethod

A general integral evaluation method that will automatically select the
appropriate quadrature method to approximate the integral. Please note that this
will generate a unique set of parameter supports and will ignore existing supports
when the integral is evaluated and thus should be used with caution. However,
this method is able to handle infinite and semi-infinite integral domains.
This is not compatible with individual dependent parameters. Contains no fields.
"""
struct Quadrature <: AbstractUnivariateMethod end

"""
<<<<<<< HEAD
Integral evaluation method that allows for the user to specify supports to be included in
quadrature evaluation. This method uses Gauss Lobatto quadrature to decompose the overall Integral
into smaller integrals that span the user defined supports as follows:

``\\int_{x_1}^{x_3} f(x) dx = \\int_{x_1}^{x_2} f(x) + \\int_{x_2}^{x_3} f(x)``

where the integrals are evaluated using Gauss Lobatto quadrature:

``\\int f(x) dx \\approx \\sum_{i=1}^{n} \\alpha_i f(x_i)``

=======
    GaussHermite <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Hermite quadrature to
evaluate integrals. This is valid for infinite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters. Contains no fields.
>>>>>>> f36e98edf94dd41f108f13ce112e07afc3151f70
"""
struct FEGaussLobatto <: AbstractUnivariateMethod end

"""
<<<<<<< HEAD
    GaussLegendre <: FiniteGaussQuad
=======
    GaussLegendre <: AbstractUnivariateMethod
>>>>>>> f36e98edf94dd41f108f13ce112e07afc3151f70

An integral evaulation method that uses Gauss-Legendre quadrature to
evaluate integrals. This is valid for finite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters. Contains no fields.
"""
struct GaussLegendre <: FiniteGaussQuad end

"""
    GaussRadau <: FiniteGaussQuad

An integral evaulation method that uses Gauss-Radau quadrature to
evaluate integrals. This is valid for finite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters. Contains no fields.
"""
struct GaussRadau <: FiniteGaussQuad end

"""
    GaussLobatto <: FiniteGaussQuad

An integral evaulation method that uses Gauss-Lobatto quadrature to
evaluate integrals. This is valid for finite integral domains. Note this will
generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters. Contains no fields.
"""
struct GaussLobatto <: FiniteGaussQuad end

"""
    GaussJacobi <: FiniteGaussQuad

An integral evaulation method that uses Gauss-Jacobi quadrature to
evaluate integrals. It will take this form:

``\\int_{-1}^{1} f(x) (1-x)^\\alpha (1+x)^\\beta dx \\approx \\sum_{i=1}^{n} \\alpha_i f(x_i)``

Where, 

`` (1-x)^\\alpha (1+x)^\\beta`` 

is the weight function.
This is valid for finite integral domains. This requires the user
to input the alpha and beta shape parameters for their function. This will then
generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters.
If α or β is < -1, an error will be returned.  

**Fields**
- `α::Float64`: Shape parameter that must be > -1
- `β::Float64`: Shape parameter that must be > -1
"""
struct GaussJacobi <: FiniteGaussQuad
    α::Float64
    β::Float64
    function GaussJacobi(α::Float64, β::Float64)
        if α < -1 || β < -1
            error("α and β must be greater than -1")
        end
        return new(α, β)
    end
end

"""
    GaussChebyshev <: FiniteGaussQuad

An integral evaulation method that uses Gauss-Chebyshev quadrature to
evaluate integrals. This is valid for finite integral domains. This requires the user
to input the order of Guass-Chebyshev Quadrature they want to use. 
If the order is not between 1 and 4 an error will be returned. 
The integral evaluated is as follows:

``\\int_{-1}^{1} f(x) w(x) \\approx \\sum_{i=1}^{n} \\alpha_i f(x_i)``

The weight functions are as follows: 

- 1st order: ``w(x)  =  \\frac{1}{\\sqrt{1-x^2}}``
- 2nd order: ``w(x) = {\\sqrt{1-x^2}}``
- 3rd order: ``w(x) = \\sqrt{(1+x)/(1-x)}``
- 4th order: ``w(x) = \\sqrt{(1-x)/(1+x)}``

This will then generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters.

**Fields**
- `order::Int64`: Specifies the order of Gauss-Chebyshev Quadrature. Must be between 1 and 4.
"""
struct GaussChebyshev <: FiniteGaussQuad
    order::Int64
    function GaussChebyshev(order::Int64)
        if order < 1 || order > 4
            error("Order must be between 1 and 4")
        end
        return new(order)
    end
end

"""
    GaussLaguerre <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Laguerre quadrature to
evaluate integrals. This is valid for semi-infinite integral domains. 

This method evaluates the following integral:

``\\int_{0}^{+∞} f(x) e^{-x} \\approx \\sum_{i=1}^{n} \\alpha_i f(x_i)``

Using the weight function:

``w(x) = e^{-x}``

Note this will generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters.
"""
struct GaussLaguerre <: AbstractUnivariateMethod end

"""
    GaussHermite <: AbstractUnivariateMethod

An integral evaulation method that uses Gauss-Hermite quadrature to
evaluate integrals. This is valid for infinite integral domains. 
It will take this form:

``\\int_{-∞}^{∞} f(x) e^{-x^2} \\approx \\sum_{i=1}^{n} \\alpha_i f(x_i)``

Using the weight function:
``w(x)`` = ``e^{-x^2}``

Note this will generate its own set of supports and will ignore other parameter supports.
This is not compatible with individual dependent parameters.
"""
struct GaussHermite <: AbstractUnivariateMethod end

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
for finite integral domains. If an array of independent infinite parameters
is specified, they must use the same amount of supports. Contains no fields.
"""
struct MultiMCSampling <: AbstractMultivariateMethod end

"""
    MultiIndepMCSampling <: AbstractMultivariateMethod

An integral evaluation method that uses uniform Monte Carlo sampling to
approximate the integral similar to [`MultiMCSampling`](@ref MeasureToolbox.MultiMCSampling).
However, this variant will generate its own set of supports and ignore all other
supports with the `MCSample` label. Note this is valid only for finite integral
domains. Contains no fields.
"""
struct MultiIndepMCSampling <: AbstractMultivariateMethod end



################################################################################
#                             INTERNAL SUPPORT LABELS
################################################################################
"""
    InternalGaussLobatto <: InfiniteOpt.InternalLabel

A support label Gauss Lobatto points that are used as generative supports.
"""
struct InternalGaussLobatto <: InfiniteOpt.InternalLabel end

################################################################################
#                       DATA GENERATION METHODS (UNIVARIATE)
################################################################################
"""
    generate_integral_data(
        prefs::Union{InfiniteOpt.GeneralVariableRef, Vector{InfiniteOpt.GeneralVariableRef}},
        lower_bounds::Union{Real, Vector{<:Real}},
        upper_bounds::Union{Real, Vector{<:Real}},
        method::V; [num_supports::Int = InfiniteOpt.DefaultNumSupports,
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
          "in accordance with the arguments given.")
end

# Single pref with Automatic
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
    lower_bound::Real,
    upper_bound::Real,
    method::Automatic;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight,
    kwargs...)::InfiniteOpt.AbstractMeasureData
    is_depend = InfiniteOpt._index_type(pref) == DependentParameterIndex
    inf_bound_num = (lower_bound == -Inf) + (upper_bound == Inf)
    if inf_bound_num == 0 # finite interval
        method = UniTrapezoid()
    elseif inf_bound_num == 1 && !is_depend # semi-infinite interval
        method = GaussLaguerre()
    elseif inf_bound_num == 2 && !is_depend # infinite interval
        method = GaussHermite()
    else
        error("Cannot generate measure data for individual dependent parameters " *
              "with infinite or semi-infinite domains.")
    end
    return generate_integral_data(pref, lower_bound, upper_bound, method,
                                  num_supports = num_supports,
                                  weight_func = weight_func)
end

#Function for generating coefficients for FE Gauss Lobatto Quadrature
function _lobatto_coeff(points::Array{Float64}, num_nodes::Int = 10)
    sort!(points)
    (_, coeffs) = FastGaussQuadrature.gausslobatto(num_nodes)
    coeff_vector = zeros((length(points) - 1)*(num_nodes - 1) + 1)
    idx = 1:num_nodes
    for i in 1:length(points)-1
        coeffs_modified = ((points[i+1] - points[i]) / 2 * coeffs)
        coeff_vector[idx] += coeffs_modified
        idx = idx .+ (num_nodes - 1)
    end
    return coeff_vector
end

#Single Pref Finite Gauss Lobatto Quadrature
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
        lower_bound::Real,
        upper_bound::Real,
        method::FEGaussLobatto;
        num_nodes::Int = 3,
        weight_func::Function = InfiniteOpt.default_weight,
        kwargs...)
    if num_nodes < 2
        error("Num supports must be greater than 2")
    end
    (nodes, _) = FastGaussQuadrature.gausslobatto(num_nodes)
    info = UniformGenerativeInfo(nodes[2:end-1], InternalGaussLobatto, -1, 1)
    coeff_func(supps) = _lobatto_coeff(supps, num_nodes)
    return InfiniteOpt.FunctionalDiscreteMeasureData(pref, coeff_func, 0,
                        InfiniteOpt.All, info, weight_func,
                        lower_bound, upper_bound, false)
end

# Univariate trapezoid coefficient function
function _trapezoid_coeff(supps::Vector{<:Real})::Vector{Float64}
    len = length(supps)
    len >= 2 || error("Cannot invoke trapezoid rule on integral if there are " * 
                      "less than 2 supports added to its infinite parameter. Ensure " * 
                      "all infinite parameters have supports.")
    coeffs = Vector{Float64}(undef, len)
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

# Single pref trapezoid rule
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::UniTrapezoid;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    return InfiniteOpt.FunctionalDiscreteMeasureData(pref, _trapezoid_coeff, 0,
                                                     InfiniteOpt.All, 
                                                     InfiniteOpt.NoGenerativeSupports(), 
                                                     weight_func, lower_bound, 
                                                     upper_bound, false)
end

# Useful error function for dependent parameters
function _ensure_independent_param(pref::InfiniteOpt.GeneralVariableRef,
                                   method)::Nothing
    if InfiniteOpt._index_type(pref) == InfiniteOpt.DependentParameterIndex
        error("Cannot generate measure data for individual dependent parameters " *
              "via `$(method)`. Try using `UniTrapezoid` for finite domains.")
    end
    return
end

# Single pref general quadrature dispatch
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::Quadrature;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    _ensure_independent_param(pref, method)
    inf_bound_num = (lower_bound == -Inf) + (upper_bound == Inf)
    if inf_bound_num == 0 # finite interval
        method = GaussLegendre()
    elseif inf_bound_num == 1 # semi-infinite interval
        method = GaussLaguerre()
    else # infinite interval
        method = GaussHermite()
    end
    return generate_integral_data(pref, lower_bound, upper_bound, method,
                                  num_supports = num_supports, weight_func = weight_func)
end

#Single Pref Gauss-Legendre
function _make_nodes_weights(method::GaussLegendre, num_nodes::Int)
    return FastGaussQuadrature.gausslegendre(num_nodes)
end

#Single Pref Gauss-Lobatto
function _make_nodes_weights(method::GaussLobatto, num_nodes::Int)
    return FastGaussQuadrature.gausslobatto(num_nodes)
end

#Single Pref Gauss-Radau
function _make_nodes_weights(method::GaussRadau, num_nodes::Int)
    return FastGaussQuadrature.gaussradau(num_nodes)
end

#Single Pref Gauss-Chebyshev
function _make_nodes_weights(method::GaussChebyshev, num_nodes::Int)
    return FastGaussQuadrature.gausschebyshev(num_nodes, method.order)
end

#Single Pref Gauss-Jacobi
function _make_nodes_weights(method::GaussJacobi, num_nodes::Int)
    return FastGaussQuadrature.gaussjacobi(num_nodes, method.α, method.β)
end

# Single pref Finite Gauss Quadrature
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::FiniteGaussQuad;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    _ensure_independent_param(pref, method)
    if lower_bound == -Inf || upper_bound == Inf
        @warn("The `$(typeof(method))` method can only be applied on finite intervals, " *
              "switching to an appropriate method.")
        return generate_integral_data(pref, lower_bound, upper_bound, Quadrature(),
                                      num_supports = num_supports,
                                      weight_func = weight_func)
    end
    # (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    (supports, coeffs) = _make_nodes_weights(method, num_supports)
    supports = (upper_bound - lower_bound) / 2 * supports .+ (upper_bound + lower_bound) / 2
    coeffs = (upper_bound - lower_bound) / 2 * coeffs
    return InfiniteOpt.DiscreteMeasureData(pref, coeffs, supports, 
                                           InfiniteOpt.generate_unique_label(),
                                           weight_func, lower_bound, upper_bound,
                                           false)
end


# Single pref Gauss-Laguerre
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::GaussLaguerre;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    _ensure_independent_param(pref, method)
    if upper_bound == Inf
        if lower_bound == -Inf
            @warn("Gauss Laguerre quadrature can only be applied on semi-infinite intervals, " *
                  "switching to an appropriate method.")
            return generate_integral_data(pref, lower_bound, upper_bound, GaussHermite(),
                                          num_supports = num_supports,
                                          weight_func = weight_func)
        end
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = supports .+ lower_bound
    elseif lower_bound == -Inf
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = -supports .+ upper_bound
    else
        @warn("Gauss Laguerre quadrature can only be applied on semi-infinite intervals, " *
              "switching to an appropriate method.")
        return generate_integral_data(pref, lower_bound, upper_bound, GaussLegendre(),
                                      num_supports = num_supports,
                                      weight_func = weight_func)
    end

    return InfiniteOpt.DiscreteMeasureData(pref, coeffs, supports, 
                                           InfiniteOpt.generate_unique_label(),
                                           weight_func, lower_bound, upper_bound,
                                           false)
end

# Single pref Gauss-Hermite
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::GaussHermite;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    _ensure_independent_param(pref, method)
    if lower_bound != -Inf || upper_bound != Inf
        @warn("Gauss Hermite quadrature can only be applied on infinite intervals, " *
              "switching to an appropriate method.")
        return generate_integral_data(pref, lower_bound, upper_bound, Quadrature(),
                                      num_supports = num_supports,
                                      weight_func = weight_func)
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return InfiniteOpt.DiscreteMeasureData(pref, coeffs, supports, 
                                           InfiniteOpt.generate_unique_label(),
                                           weight_func, lower_bound, upper_bound,
                                           false)
end

# Single pref uniform Monte Carlo sampling
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::UniMCSampling;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    # check and process the arguments
    if lower_bound == -Inf || upper_bound == Inf
        error("Univariate MC sampling is not supported for (semi-)infinite intervals.")
    elseif InfiniteOpt._index_type(pref) == InfiniteOpt.DependentParameterIndex
        if !(num_supports in [0, InfiniteOpt.DefaultNumSupports])
            @warn("Cannot specify a nonzero minimum supports for an individual " *
                  "dependent parameter. Setting `num_supports = 0`.")
        end
        num_supports = 0
    end
    # make the coefficient function
    function _coeffs(supps::Vector{<:Real})::Vector{Float64}
        len = length(supps)
        return (upper_bound - lower_bound) * ones(len) / len
    end
    # prepare the data
    return InfiniteOpt.FunctionalDiscreteMeasureData(pref, _coeffs, num_supports,
                                                     InfiniteOpt.MCSample,
                                                     InfiniteOpt.NoGenerativeSupports(),
                                                     weight_func,
                                                     lower_bound, upper_bound,
                                                     false)
end

# Single pref independent Monte Carlo sampling
function generate_integral_data(pref::InfiniteOpt.GeneralVariableRef,
                                lower_bound::Real,
                                upper_bound::Real,
                                method::UniIndepMCSampling;
                                num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                weight_func::Function = InfiniteOpt.default_weight,
                                kwargs...)
    _ensure_independent_param(pref, method)
    if lower_bound == -Inf || upper_bound == Inf
        error("Univariate MC sampling is not applicable to (semi-)infinite intervals.")
    end
    set = InfiniteOpt.IntervalSet(lower_bound, upper_bound)
    supports, _ = InfiniteOpt.generate_support_values(set, InfiniteOpt.MCSample,
                                                      num_supports = num_supports)
    coeffs = (upper_bound - lower_bound) * ones(num_supports) / num_supports
    return InfiniteOpt.DiscreteMeasureData(pref, coeffs, supports,
                                           InfiniteOpt.MCSample, weight_func,
                                           lower_bound, upper_bound, false)
end

################################################################################
#                       DATA GENERATION METHODS (MULTIVARIATE)
################################################################################
# Vector of prefs with Automatic
function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bounds::Vector{<:Real},
    upper_bounds::Vector{<:Real},
    method::Automatic;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight,
    kwargs...)::InfiniteOpt.AbstractMeasureData
    return generate_integral_data(prefs, lower_bounds, upper_bounds,
                                  MultiMCSampling(), num_supports = num_supports,
                                  weight_func = weight_func)
end

# Vector of prefs using Monte Carlo samples
function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bounds::Vector{<:Real},
    upper_bounds::Vector{<:Real},
    method::MultiMCSampling;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight,
    kwargs...)::InfiniteOpt.AbstractMeasureData
    # check bounds
    if any(lb == -Inf for lb in lower_bounds) || any(ub == Inf for ub in upper_bounds)
        error("MC Sampling is not applicable to (semi-)infinite intervals.")
    end
    # make the coefficients function
    vol = prod(upper_bounds .- lower_bounds)
    function _coeffs(supps::Array{<:Real, 2})::Vector{Float64}
        num_samples = size(supps, 2)
        return vol * ones(num_samples) / num_samples
    end
    # make the data
    return InfiniteOpt.FunctionalDiscreteMeasureData(prefs, _coeffs,
                                                     num_supports,
                                                     InfiniteOpt.MCSample,
                                                     weight_func,
                                                     lower_bounds, upper_bounds,
                                                     false)
end

## Generate multivariate Monte Carlo samples in accordance with the parameter type
# DependentParameterRefs
function _make_multi_mc_supports(prefs::Vector{InfiniteOpt.DependentParameterRef},
                                 lbs::Vector{<:Real}, ubs::Vector{<:Real},
                                 num_supps::Int)::Matrix{Float64}
    sets = [InfiniteOpt.IntervalSet(lbs[i], ubs[i]) for i in eachindex(lbs)]
    set = InfiniteOpt.CollectionSet(sets)
    supports, _ = InfiniteOpt.generate_support_values(set, InfiniteOpt.MCSample,
                                                      num_supports = num_supps)
    return supports
end

# Vector of prefs with independent Monte Carlo supports
function generate_integral_data(prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bounds::Vector{<:Real},
    upper_bounds::Vector{<:Real},
    method::MultiIndepMCSampling;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight,
    kwargs...)::InfiniteOpt.AbstractMeasureData
    # check bounds
    if any(lb == -Inf for lb in lower_bounds) || any(ub == Inf for ub in upper_bounds)
        error("MC Sampling is not applicable to (semi-)infinite intervals.")
    end
    # prepare the supports
    sets = [InfiniteOpt.IntervalSet(lower_bounds[i], upper_bounds[i])
            for i in eachindex(lower_bounds)]
    set = InfiniteOpt.CollectionSet(sets)
    supports, _ = InfiniteOpt.generate_support_values(set, InfiniteOpt.MCSample,
                                                      num_supports = num_supports)
    # prepare the coefficients
    coeffs = prod(upper_bounds .- lower_bounds) * ones(num_supports) / num_supports
    # prepare the data
    return InfiniteOpt.DiscreteMeasureData(prefs, coeffs, supports,
                                           InfiniteOpt.MCSample, weight_func,
                                           lower_bounds, upper_bounds, false)
end

################################################################################
#                           UNIVARIATE INTEGRALS
################################################################################
# Define default keyword arguments for 1-D integrals
const UniIntegralDefaults = Dict(:eval_method => Automatic(),
                                 :num_supports => InfiniteOpt.DefaultNumSupports,
                                 :num_nodes => 3,
                                 :weight_func => InfiniteOpt.default_weight)

"""
    uni_integral_defaults()::Dict{Symbol, Any}

Get the default keyword argument values for defining one-dimensional integrals.

```julia-repl
julia> uni_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports          => 10
  :eval_method           => Automatic()
  :weight_func           => default_weight
```
"""
uni_integral_defaults()::Dict{Symbol, Any} = UniIntegralDefaults

"""
    set_uni_integral_defaults(; kwargs...)::Nothing

Set the default keyword argument settings for one-dimensional integrals.
The keyword arguments of this function will be recorded in the default keyword
argument dictionary. These will determine the default keyword argument values
when calling [`integral`](@ref MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::InfiniteOpt.GeneralVariableRef,::Real, ::Real))
with a single infinite parameter.

**Example**
```julia-repl
julia> uni_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports          => 10
  :eval_method           => Automatic()
  :weight_func           => default_weight

julia> set_uni_integral_defaults(num_supports = 5, eval_method = Quadrature(),
                                 new_kwarg = true)

julia> uni_integral_defaults()
Dict{Symbol,Any} with 4 entries:
  :new_kwarg             => true
  :num_supports          => 5
  :eval_method           => Quadrature()
  :weight_func           => default_weight
```
"""
function set_uni_integral_defaults(; kwargs...)::Nothing
    merge!(UniIntegralDefaults, kwargs)
    return
end

"""
    integral(expr::JuMP.AbstractJuMPScalar,
             pref::GeneralVariableRef,
             [lower_bound::Real = lower_bound(pref),
             upper_bound::Real = upper_bound(pref);
             kwargs...])::GeneralVariableRef

Returns a measure reference that evaluates the integral of `expr` with respect
to infinite parameter `pref` from `lower_bound` to `upper_bound`. This thus
considers integrals of the form: ``\\int_{p \\in P} expr(p) w(p) dp`` where ``p``
is an infinite parameter and ``w`` is the weight function is 1 by default. This
function provides a high-level interface that ultimately constructs an appropriate
concrete form of [`AbstractMeasureData`](@ref) via [`generate_integral_data`](@ref)
in accordance with the keyword arugment `eval_method` that is then used with
[`measure`](@ref). Note that it is preferred to call [`@integral`](@ref) when
`expr` is not just a single variable reference. Errors for bad bound input.

The keyword arguments are as follows:
- `eval_method::AbstractUnivariateMethod`: Used to determine the
    numerical evaluation scheme. Possible choices include:
    - [`Automatic`](@ref MeasureToolbox.Automatic)
    - [`UniTrapezoid`](@ref MeasureToolbox.UniTrapezoid)
    - [`UniMCSampling`](@ref MeasureToolbox.UniMCSampling)
    - [`UniIndepMCSampling`](@ref MeasureToolbox.UniIndepMCSampling)
    - [`Quadrature`](@ref MeasureToolbox.Quadrature)
    - [`GaussHermite`](@ref MeasureToolbox.GaussHermite)
    - [`GaussLegendre`](@ref MeasureToolbox.GaussLegendre)
    - [`FEGaussLobatto`](@ref MeasureToolbox.FEGaussLobatto)
    - [`GaussLageurre`](@ref MeasureToolbox.GaussLaguerre)
    - [`GaussLobatto`](@ref MeasureToolbox.GaussLobatto)
    - [`GaussChebyshev`](@ref MeasureToolbox.GaussChebyshev)
    - [`GaussRadau`](@ref MeasureToolbox.GaussRadau)
    - [`GaussJacobi`](@ref MeasureToolbox.GaussJacobi)
- `num_supports`: The minimum number of supports to be generated (if used by
    `eval_method`)
- `weight_func`: ``w(p)`` above with parameter value inputs and scalar output

See [`set_uni_integral_defaults`](@ref) to update the default keyword argument
values for all one-dimensional integral calls.

**Example**
```julia-repl
julia> @infinite_parameter(model, x in [0, 1])
x

julia> @infinite_variable(model, f(x))
f(x)

julia> int = integral(f, x)
∫{x ∈ [0, 1]}[f(x)]

julia> expand(int)
0.2 f(0.8236475079774124) + 0.2 f(0.9103565379264364) + 0.2 f(0.16456579813368521) + 0.2 f(0.17732884646626457) + 0.2 f(0.278880109331201)
```
"""
function integral(expr::JuMP.AbstractJuMPScalar,
                  pref::InfiniteOpt.GeneralVariableRef,
                  lower_bound::Real = NaN,
                  upper_bound::Real = NaN;
                  kwargs...)::InfiniteOpt.GeneralVariableRef
    # check parameter formatting
    InfiniteOpt._check_params(pref)
    # fill in bounds if needed
    set = infinite_set(pref)
    if isnan(lower_bound) && JuMP.has_lower_bound(set)
        lower_bound = JuMP.lower_bound(pref)
    end
    if isnan(upper_bound) && JuMP.has_upper_bound(set)
        upper_bound = JuMP.upper_bound(pref)
    end
    # ensure valid bounds
    if lower_bound >= upper_bound
        error("Invalid integral bounds, ensure that lower_bound < upper_bound.")
    elseif !InfiniteOpt.supports_in_set([lower_bound, upper_bound], set)
        error("Integral bounds violate the infinite domain.")
    end
    # prepare the keyword arguments and make the measure data
    processed_kwargs = merge(uni_integral_defaults(), kwargs)
    eval_method = pop!(processed_kwargs, :eval_method)
    data = generate_integral_data(pref, lower_bound, upper_bound,
                                  eval_method; processed_kwargs...)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "integral")
end

"""
    ∫(expr::JuMP.AbstractJuMPScalar,
      pref::GeneralVariableRef,
      [lower_bound::Real = NaN,
      upper_bound::Real = NaN;
      kwargs...])::GeneralVariableRef

A convenient wrapper for [`integral`](@ref). The `∫` unicode symbol is produced 
via `\\int`.
"""
function ∫(expr::JuMP.AbstractJuMPScalar,
           pref::InfiniteOpt.GeneralVariableRef,
           lower_bound::Real = NaN,
           upper_bound::Real = NaN;
           kwargs...)::InfiniteOpt.GeneralVariableRef
    return integral(expr, pref, lower_bound, upper_bound; kwargs...)
end

################################################################################
#                            MULTIVARIATE INTEGRALS
################################################################################
# Define default keyword arguments for multi-D integrals
const MultiIntegralDefaults = Dict(:eval_method => Automatic(),
                                   :num_supports => InfiniteOpt.DefaultNumSupports,
                                   :weight_func => InfiniteOpt.default_weight)

"""
    multi_integral_defaults()::Dict{Symbol, Any}

Get the default keyword argument values for defining multi-dimensional integrals.

```julia-repl
julia> multi_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports          => 10
  :eval_method           => Automatic()
  :weight_func           => default_weight
```
"""
multi_integral_defaults()::Dict{Symbol, Any} = MultiIntegralDefaults

"""
    set_multi_integral_defaults(; kwargs...)::Nothing

Set the default keyword argument settings for multi-dimesnional integrals.
The keyword arguments of this function will be recorded in the default keyword
argument dictionary. These will determine the default keyword argument values
when calling [`integral`](@ref MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::AbstractArray{GeneralVariableRef},::Union{Real, AbstractArray{<:Real}}, ::Union{Real, AbstractArray{<:Real}}))
with an array of infinite parameters.

**Example**
```julia-repl
julia> multi_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports          => 10
  :eval_method           => Automatic()
  :weight_func           => default_weight

julia> set_multi_integral_defaults(num_supports = 5, new_kwarg = true)

julia> multi_integral_defaults()
Dict{Symbol,Any} with 4 entries:
  :new_kwarg             => true
  :num_supports          => 5
  :eval_method           => Automatic()
  :weight_func           => default_weight
```
"""
function set_multi_integral_defaults(; kwargs...)::Nothing
    merge!(MultiIntegralDefaults, kwargs)
    return
end

Base.isnan(arr::AbstractArray{<:Real})::Bool = all(Base.isnan.(arr))
"""
    integral(expr::JuMP.AbstractJuMPScalar,
             prefs::AbstractArray{GeneralVariableRef},
             [lower_bounds::Union{Real, AbstractArray{<:Real}} = [lower_bound(pref)...],
             upper_bounds::Union{Real, AbstractArray{<:Real}} = [upper_bound(pref)...];
             kwargs...])::GeneralVariableRef

Returns a measure reference that evaluates the integral of `expr` with respect
to infinite parameters `prefs` from `lower_bounds` to `upper_bounds`. This thus
considers integrals of the form: ``\\int_{p \\in P} expr(p) w(p) dp`` where ``p``
is an infinite parameter and ``w`` is the weight function is 1 by default. This
function provides a high-level interface that ultimately constructs an appropriate
concrete form of [`AbstractMeasureData`](@ref) via [`generate_integral_data`](@ref)
in accordance with the keyword arugment `eval_method` that is then used with
[`measure`](@ref). Note that it is preferred to call [`@integral`](@ref) when
`expr` is not just a single variable reference. Errors when the container types
and dimensions do not match or the bounds are invalid.

The keyword arguments are as follows:
- `eval_method::AbstractMultivariateMethod`: Used to determine the
    numerical evaluation scheme. Possible choices include:
    - [`Automatic`](@ref MeasureToolbox.Automatic)
    - [`MultiMCSampling`](@ref MeasureToolbox.MultiMCSampling)
    - [`MultiIndepMCSampling`](@ref MeasureToolbox.MultiIndepMCSampling)
- `num_supports`: The minimum number of supports to be generated (if used by
    `eval_method`)
- `weight_func`: ``w(p)`` above with parameter value inputs and scalar output

See [`set_multi_integral_defaults`](@ref) to update the default keyword argument
values for all multi-dimensional integral calls.

**Example**
```julia-repl
julia> @infinite_parameter(model, x[1:2] in [0, 1], independent = true);

julia> @infinite_variable(model, f(x));

julia> int = integral(f, x)
∫{x ∈ [0, 1]^2}[f(x)]
```
"""
function integral(expr::JuMP.AbstractJuMPScalar,
                  prefs::AbstractArray{InfiniteOpt.GeneralVariableRef},
                  lower_bounds::Union{Real, AbstractArray{<:Real}} = NaN,
                  upper_bounds::Union{Real, AbstractArray{<:Real}} = NaN;
                  kwargs...)::InfiniteOpt.GeneralVariableRef
    # check parameter formatting
    InfiniteOpt._check_params(prefs)
    # fill in the lower bounds if needed
    if isnan(lower_bounds) && JuMP.has_lower_bound(first(prefs))
        lbs = map(p -> JuMP.lower_bound(p), prefs)
    elseif !isnan(lower_bounds) && lower_bounds isa Real
        lbs = map(p -> lower_bounds, prefs)
    else
        lbs = lower_bounds
    end
    # fill in the upper bounds if needed
    if isnan(upper_bounds) && JuMP.has_upper_bound(first(prefs))
        ubs = map(p -> JuMP.upper_bound(p), prefs)
    elseif !isnan(upper_bounds) && upper_bounds isa Real
        ubs = map(p -> upper_bounds, prefs)
    else
        ubs = upper_bounds
    end
    # do initial bound checks
    if InfiniteOpt._keys(prefs) != InfiniteOpt._keys(lbs) ||
       InfiniteOpt._keys(prefs) != InfiniteOpt._keys(ubs)
        error("Array keys of infinite parameters and bounds do not match.")
    end
    # process the arrays
    vector_prefs = InfiniteOpt._make_ordered_vector(prefs)
    vector_lbs = InfiniteOpt._make_ordered_vector(lbs)
    vector_ubs = InfiniteOpt._make_ordered_vector(ubs)
    # ensure valid bounds
    if any(vector_lbs[i] >= vector_ubs[i] for i in eachindex(vector_lbs))
        error("Invalid integral bounds, ensure that lower_bounds < upper_bounds.")
    end
    dprefs = map(p -> dispatch_variable_ref(p), vector_prefs)
    InfiniteOpt._check_bounds_in_set(dprefs, vector_lbs, vector_ubs)
    # prepare the keyword arguments and make the measure data
    processed_kwargs = merge(multi_integral_defaults(), kwargs)
    eval_method = pop!(processed_kwargs, :eval_method)
    data = generate_integral_data(vector_prefs, vector_lbs, vector_ubs,
                                  eval_method; processed_kwargs...)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "integral")
end

"""
    ∫(expr::JuMP.AbstractJuMPScalar,
      prefs::AbstractArray{GeneralVariableRef},
      [lower_bounds::Union{Real, AbstractArray{<:Real}} = NaN,
      upper_bounds::Union{Real, AbstractArray{<:Real}} = NaN;
      kwargs...])::GeneralVariableRef

A convenient wrapper for [`integral`](@ref). The unicode symbol `∫` is produced 
via `\\int`.
"""
function ∫(expr::JuMP.AbstractJuMPScalar,
           prefs::AbstractArray{InfiniteOpt.GeneralVariableRef},
           lower_bounds::Union{Real, AbstractArray{<:Real}} = NaN,
           upper_bounds::Union{Real, AbstractArray{<:Real}} = NaN;
           kwargs...)::InfiniteOpt.GeneralVariableRef
    return integral(expr, prefs, lower_bounds, upper_bounds; kwargs...)
end

"""
    @integral(expr::JuMP.AbstractJuMPScalar,
              prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}},
              [lower_bounds::Union{Real, AbstractArray{<:Real}} = default_bounds,
              upper_bounds::Union{Real, AbstractArray{<:Real}} = default_bounds;
              kwargs...])::GeneralVariableRef

An efficient wrapper for [`integral`](@ref integral(::JuMP.AbstractJuMPScalar, ::InfiniteOpt.GeneralVariableRef, ::Real, ::Real))
and [`integral`](@ref integral(::JuMP.AbstractJuMPScalar, ::AbstractArray{InfiniteOpt.GeneralVariableRef}, ::Union{Real, AbstractArray{<:Real}}, ::Union{Real, AbstractArray{<:Real}})).
Please see the above doc strings for more information.
"""
macro integral(expr, prefs, args...)
    _error(str...) = InfiniteOpt._macro_error(:integral, (expr, prefs, args...), 
                                              str...)
    extra, kw_args, requestedcontainer = InfiniteOpt._extract_kw_args(args)
    if length(extra) != 0 && length(extra) != 2
        _error("Incorrect number of positional arguments for @integral. " *
               "Must provide both bounds or no bounds.")
    end
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( integral($expression, $prefs, $(extra...); ($(kw_args...))) )
    return esc(mref)
end

"""
    @∫(expr::JuMP.AbstractJuMPScalar,
       prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}},
       [lower_bounds::Union{Real, AbstractArray{<:Real}} = default_bounds,
       upper_bounds::Union{Real, AbstractArray{<:Real}} = default_bounds;
       kwargs...])::GeneralVariableRef

A convenient wrapper for [`@integral`](@ref). The unicode symbol `∫` is produced 
via `\\int`.
"""
macro ∫(expr, prefs, args...)
    return esc(:( @integral($expr, $prefs, $(args...)) ))
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
julia> (supps, coeffs) = infinite_transform(-Inf, Inf, 5, sub_method = gauss_legendre())
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
