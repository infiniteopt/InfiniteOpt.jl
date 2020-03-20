const mc_sampling = :mc_sampling
const gauss_hermite = :gauss_hermite
const gauss_legendre = :gauss_legendre
const gauss_laguerre = :gauss_laguerre
const trapezoid = :trapezoid

"""
    generate_measure_data(params::Union{InfiniteOpt.ParameterRef,
                          AbstractArray{<:InfiniteOpt.ParameterRef}},
                          num_supports::Int,
                          lb::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing,
                          ub::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing;
                          eval_method::Function = mc_sampling, name::String = "",
                          weight_func::Function = InfiniteOpt._w,
                          check_method::Bool = true, kwargs...
                          )::InfiniteOpt.AbstractMeasureData

Generate an [`AbstractMeasureData`](@ref) object that automatically generate
supports based on a set of given information. The information is required to
include parameters and number of supports. Other optional information to input
includes lower and upper bounds, measure name, support generation method, and
weight functions. The users could supply extra keyword arguments if necessary
for their custom support generation methods.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(m, x in [0., 1.]);

julia> measure_data = generate_measure_data(x, 3, 0.3, 0.7, method = gauss_legendre)
DiscreteMeasureData(x, [0.1111111111111111, 0.17777777777777776, 0.1111111111111111], [0.3450806661517033, 0.5, 0.6549193338482967], "", InfiniteOpt._w)
```
"""
function generate_measure_data(params::Union{InfiniteOpt.ParameterRef,
                               AbstractArray{<:InfiniteOpt.ParameterRef}},
                               num_supports::Int,
                               lb::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing,
                               ub::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing;
                               eval_method::Symbol = mc_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w,
                               check_method::Bool = true, kwargs...
                               )::InfiniteOpt.AbstractMeasureData
    if isa(params, InfiniteOpt.ParameterRef)
        set = InfiniteOpt._parameter_set(params)
    else
        params = convert(JuMPC.SparseAxisArray, params)
        set = InfiniteOpt._parameter_set(first(params))
        sets = InfiniteOpt._parameter_set.(params)
        if any(typeof(set) != typeof(s) for s in sets)
            error("Automatic measure data generation for multi-dimensional " *
                  "parameters with mixed set types is not supported.")
        end
    end
    (supports, coeffs) = generate_supports_and_coeffs(set, params, num_supports, lb, ub, Val(eval_method); kwargs...)
    return InfiniteOpt.DiscreteMeasureData(params, coeffs, supports,
                                           name = name, weight_function = weight_func)
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.AbstractInfiniteSet,
                                 params::Union{InfiniteOpt.ParameterRef,
                                 AbstractArray{<:InfiniteOpt.ParameterRef}},
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{eval_method}; [kwargs...])::Tuple

Generate supports and coefficients using the arguments and
`kwargs` as appropriate. This will dispatch to `eval_method` in accordance with the
type of `set`. This is intended as an internal method for
[`generate_measure_data`](@ref InfiniteOpt.MeasureEvalMethods.generate_measure_data)
and will need to be extended for user-defined
infinite set types. The output tuple should consists of two data objects, the
first being the supports and the second being the coefficients. Supports must
be a `Vector` of `Number` or of `AbstractArray{:Number}`, and coefficients
must be a `Vector` of `Number`. Extensions will be needed to implement their
own version of function `generate_supports_and_coeffs(set::AbstractInfiniteSet, params, num_supports, lb, ub, method::Val{my_method};kwargs...)`.
Refer to the other dispatches of `generate_supports_and_coeffs` for details of
specific evaluation method implementations.
"""
function generate_supports_and_coeffs end

# fallback
function generate_supports_and_coeffs(set::InfiniteOpt.AbstractInfiniteSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method; kwargs...)::Tuple
    error("`generate_supports_and_coeffs` is not extended for parameters in sets " *
          "of type $(typeof(set)) with method " * string(method)[6:length(string(method))-3] * ".")
end

# MC Sampling for IntervalSet
"""
    generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                 params::Union{InfiniteOpt.ParameterRef,
                                 AbstractArray{<:InfiniteOpt.ParameterRef}},
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{mc_sampling})::Tuple

Return a tuple that contains supports and coefficients generated by Monte Carlo
sampling from a uniform distribution between the lower and upper bounds provided.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(seed = true))
julia> @infinite_parameter(model, t in [0., 1.]);

julia> (supps, coeffs) = generate_supports_and_coeffs(IntervalSet(0,1), t, 5, 0., 1., Val(mc_sampling))
([0.8236475079774124, 0.9103565379264364, 0.16456579813368521, 0.17732884646626457, 0.278880109331201], [0.2, 0.2, 0.2, 0.2, 0.2])

julia> supps
5-element Array{Float64,1}:
 0.8236475079774124
 0.9103565379264364
 0.16456579813368521
 0.17732884646626457
 0.278880109331201
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{mc_sampling})::Tuple
    _nothing_test(method, lb, ub)
    if isa(params, InfiniteOpt.ParameterRef)
        return _mc_sampling(set, params, num_supports, lb, ub)
    else
        samples_dict = Dict()
        for i in eachindex(lb)
            (samples_dict[i], _) = _mc_sampling(set, params, num_supports, lb[i], ub[i])
        end
        samples = Array{JuMPC.SparseAxisArray, 1}(undef, num_supports)
        for j in 1:num_supports
            samples[j] = JuMPC.SparseAxisArray(Dict(k => samples_dict[k][j]
                                                    for k in eachindex(lb)))
        end
        return (samples, ones(num_supports) / num_supports * prod(ub .- lb))
    end
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.DistributionSet,
                                 params::Union{InfiniteOpt.ParameterRef,
                                 AbstractArray{<:InfiniteOpt.ParameterRef}},
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{mc_sampling})::Tuple

Return a tuple that contains supports and coefficients generated by Monte Carlo
sampling from a given distribution.

**Example**
```jldoctest; setup = :(using InfiniteOpt, Distributions, Random; model = InfiniteModel(seed = true))
julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> @infinite_parameter(model, x in dist)
x

julia> (supps, coeffs) = generate_supports_and_coeffs(DistributionSet(dist), x, 10, nothing, nothing, Val(mc_sampling))
([0.6791074260357777, 0.8284134829000359, -0.3530074003005963, -0.13485387193052173, 0.5866170746331097, 0.29733585084941616, 0.06494754854834232, -0.10901738508171745, -0.514210390833322, 1.5743302021369892], [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.DistributionSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{mc_sampling})::Tuple
    dist = set.distribution
    # create truncated distribution if necessary
    if !isa(lb, Nothing) || !isa(ub, Nothing)
        if isa(dist, Distributions.MultivariateDistribution)
            @warn("Truncated distribution for multivariate distribution is " *
                  "not supported. Lower bounds and upper bounds are ignored.")
        else
            isa(lb, Number) ? temp_lb = lb : temp_lb = -Inf
            isa(ub, Number) ? temp_ub = ub : temp_ub = Inf
            dist = Distributions.Truncated(dist, temp_lb, temp_ub)
        end
    end
    return _mc_sampling(dist, params, num_supports)
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                 params::InfiniteOpt.ParameterRef,
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{gauss_legendre})::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Legendre quadrature method. This is useful for univariate parameter in a
finite interval.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 1])
t

julia> (supps, coeffs) = generate_supports_and_coeffs(IntervalSet(0,1), t, 5, 0., 1., Val(gauss_legendre))
([0.04691007703066802, 0.23076534494715845, 0.5, 0.7692346550528415, 0.9530899229693319], [0.11846344252809454, 0.23931433524968324, 0.28444444444444444, 0.23931433524968324, 0.11846344252809454])

julia> supps
5-element Array{Float64,1}:
 0.04691007703066802
 0.23076534494715845
 0.5
 0.7692346550528415
 0.9530899229693319
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.ParameterRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{gauss_legendre})::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs)
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                 params::InfiniteOpt.ParameterRef,
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{gauss_hermite})::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Hermite quadrature method. This is useful for univariate parameter in an
infinite interval.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [-Inf, Inf])
t

julia> (supps, coeffs) = generate_supports_and_coeffs(IntervalSet(-Inf,Inf), t, 5, -Inf, Inf, Val(gauss_hermite))
([-2.0201828704560856, -0.9585724646138196, -8.881784197001252e-16, 0.9585724646138196, 2.0201828704560856], [1.1814886255359844, 0.986580996751429, 0.9453087204829428, 0.986580996751429, 1.1814886255359844])

julia> supps
5-element Array{Float64,1}:
 -2.0201828704560856
 -0.9585724646138196
 -8.881784197001252e-16
  0.9585724646138196
  2.0201828704560856
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.ParameterRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{gauss_hermite})::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    if lb != -Inf || ub != Inf
        error("Lower/upper bound is not infinity. Use other measure evaluation " *
              "methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs)
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                 params::InfiniteOpt.ParameterRef,
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{gauss_laguerre})::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Laguerre quadrature method. This is useful for univariate parameter in a
semi-infinite interval.


**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [-Inf, 0])
t

julia> (supps, coeffs) = generate_supports_and_coeffs(IntervalSet(-Inf,0), t, 5, -Inf, 0, Val(gauss_laguerre))
([-0.2635603197181408, -1.413403059106515, -3.596425771040715, -7.08581000585883, -12.640800844275773], [0.6790940422077494, 1.638487873602747, 2.7694432423708255, 4.3156569009208585, 7.219186354354335])

julia> supps
5-element Array{Float64,1}:
  -0.2635603197181408
  -1.413403059106515
  -3.596425771040715
  -7.08581000585883
 -12.640800844275773
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.ParameterRef,
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
        return (supports, coeffs)
    elseif lb == -Inf
        (supports, coeffs) = FastGaussQuadrature.gausslaguerre(num_supports)
        coeffs = coeffs .* exp.(supports)
        supports = -supports .+ ub
        return (supports, coeffs)
    else
        error("The range is bounded. Use other measure evaluation methods.")
    end
end

"""
    generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                 params::InfiniteOpt.ParameterRef,
                                 num_supports::Int,
                                 lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                 method::Val{trapezoid})::Tuple

Return a tuple that contains supports and coefficients generated a uniform
trapezoid rule. This is useful for both univariate parameter in a bounded
interval. For multivariate parameters, `num_supports` will be the number of
grid points on one dimension.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 1])
t

julia> (supps, coeffs) = generate_supports_and_coeffs(IntervalSet(0,1), t, 3, 0, 1, Val(trapezoid))
([0.0, 0.5, 1.0], [0.25, 0.5, 0.25])

julia> supps
3-element Array{Float64,1}:
 0.0
 0.5
 1.0
```
"""
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::InfiniteOpt.ParameterRef,
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Val{trapezoid})::Tuple
    _nothing_test(method, lb, ub)
    _univariate_bounds(method, lb, ub)
    return _trapezoid_univ(lb, ub, num_supports)
end

function _trapezoid_univ(lb::Number, ub::Number, num_supports::Int)::Tuple
    increment = (ub - lb) / (num_supports - 1)
    supps = [(i - 1) * increment + lb for i in 1:num_supports]
    coeffs = increment * ones(num_supports)
    coeffs[1] = coeffs[1] / 2
    coeffs[num_supports] = coeffs[num_supports] / 2
    return (supps, coeffs)
end

# MC sampling for univariate IntervalSet
function _mc_sampling(set::InfiniteOpt.IntervalSet,
                      params::Union{InfiniteOpt.ParameterRef,
                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                      num_supports::Int,
                      lb::Number, ub::Number)::Tuple
    # MC sampling from uniform distribution over the interval [lb, ub]
    if lb == -Inf || ub == Inf
        return infinite_transform(set, params, num_supports, lb, ub, Val(mc_sampling))
    else
        samples = rand(num_supports) .* (ub - lb) .+ lb
        return (samples, ones(num_supports) / num_supports * (ub - lb))
    end
end

# MC sampling for univariate DistributionSet
function _mc_sampling(dist::Distributions.UnivariateDistribution,
                     param::InfiniteOpt.ParameterRef,
                     num_supports::Int)::Tuple
    # MC sampling - Distribution Set
    samples = rand(dist, num_supports)
    return (samples, ones(num_supports) / num_supports)
end

# MC sampling for multivariate DistributionSet
function _mc_sampling(dist::Distributions.MultivariateDistribution,
                     params::AbstractArray{<:InfiniteOpt.ParameterRef},
                     num_supports::Int)::Tuple
    samples_matrix = rand(dist, num_supports)
    ordered_pairs = sort(collect(params.data), by=x->x.second.index)
    samples_dict = Dict()
    for i in eachindex(ordered_pairs)
        samples_dict[ordered_pairs[i].first] = samples_matrix[i, :]
    end
    samples = Array{JuMPC.SparseAxisArray, 1}()
    for j in 1:num_supports
        append!(samples, [JuMP.Containers.SparseAxisArray(
                          Dict(k.first => samples_dict[k.first][j]
                          for k in ordered_pairs))])
    end
    return (samples, ones(num_supports) / num_supports)
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
    (t_supports, t_coeffs) = generate_supports_and_coeffs(set, params, num_supports, t_lb, t_ub, sub_method)
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

# TODO: consider truncated distribution
# TODO: consider adding uniform grids
