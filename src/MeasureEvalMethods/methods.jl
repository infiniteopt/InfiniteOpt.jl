"""
    eval_method_registry(model::InfiniteOpt.InfiniteModel
                         )::Dict{Type, Set{Function}}
Return the registry of the model that records which measure evaluation methods
are valid for different [`AbstractInfiniteSet`](@ref).
"""
function eval_method_registry(model::InfiniteOpt.InfiniteModel
                              )::Dict{Type, Set{Function}}
    return model.meas_method_registry
end

"""
    register_eval_method(model::InfiniteOpt.InfiniteModel,
                         set_type::Type,
                         methods::Union{Function, Array{Function, 1}})
Set the registry of the model that records which measure evaluation methods
are valid for different [`AbstractInfiniteSet`](@ref). This function allows for
addition of new set types and modification of acceptable methods for existing
types.
"""
function register_eval_method(model::InfiniteOpt.InfiniteModel,
                              set_type::Type,
                              methods::Union{Function, Array{Function, 1}})
    if !(set_type <: InfiniteOpt.AbstractInfiniteSet)
        error("Set type must be a subtype of AbstractInfiniteSet.")
    end
    if methods isa Function
        methods = [methods]
    end
    if set_type in keys(model.meas_method_registry)
        union!(model.meas_method_registry[set_type], Set{Function}(methods))
    else
        model.meas_method_registry[set_type] = Set{Function}(methods)
    end
    return
end

# check if a method is valid for a set
function _set_method_check(model::InfiniteOpt.InfiniteModel,
                           set::InfiniteOpt.AbstractInfiniteSet,
                           method::Function)
    registry = eval_method_registry(model)
    set_type = _set_type(model, set)
    if !(set_type in keys(registry))
        error("The parameter set type $(typeof(set)) does not have valid " *
              "measure evaluation methods.")
    end
    if !(method in registry[set_type])
        error("Method $(method) is not valid for set type $(typeof(set)).")
    end
    return
end

# return parameterized set type without parameters
# this is needed because typeof() returns parameterized set type with parameters
function _set_type(model::InfiniteOpt.InfiniteModel,
                   set::InfiniteOpt.AbstractInfiniteSet)::Union{Type, Nothing}
    for i in keys(model.meas_method_registry)
        if isa(set, i)
            return i
        end
    end
    return nothing
end

"""
    generate_measure_data(params::Union{InfiniteOpt.ParameterRef,
                          AbstractArray{<:InfiniteOpt.ParameterRef}},
                          num_supports::Int,
                          lb::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing,
                          ub::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing;
                          eval_method::Function = mc_sampling, name::String = "",
                          weight_func::Function = InfiniteOpt._w, kwargs...
                          )::InfiniteOpt.AbstractMeasureData

Generate an [`AbstractMeasureData`](@ref) object that automatically generate
supports based on a set of given information. The information is required to
include parameters and number of supports. Other optional information to input
includes lower and upper bounds, measure name, support generation method, and
weight functions. The users could supply extra keyword arguments if necessary
for their custom support generation methods.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(m, x in [0., 1.])
x

julia> measure_data = generate_measure_data(x, 3, 0.3, 0.7, method = gauss_legendre)
DiscreteMeasureData(x, [0.1111111111111111, 0.17777777777777776, 0.1111111111111111], [0.3450806661517033, 0.5, 0.6549193338482967], "", InfiniteOpt._w)
```
"""
function generate_measure_data(params::Union{InfiniteOpt.ParameterRef,
                               AbstractArray{<:InfiniteOpt.ParameterRef}},
                               num_supports::Int,
                               lb::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing,
                               ub::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing;
                               eval_method::Function = mc_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w, kwargs...
                               )::InfiniteOpt.AbstractMeasureData
    if isa(params, InfiniteOpt.ParameterRef)
        set = InfiniteOpt._parameter_set(params)
        model = params.model
    else
        params = convert(JuMPC.SparseAxisArray, params)
        set = InfiniteOpt._parameter_set(first(params))
        model = first(params).model
    end
    _set_method_check(model, set, eval_method)
    (supports, coeffs) = generate_supports_and_coeffs(set, params, num_supports, lb, ub, eval_method; kwargs...)
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
                                 method::Function; [kwargs...])::Tuple

Call `method` to generate supports and coefficients using the arguments and
`kwargs` as appropriate. This will dispatch to `method` in accordance with the
type of `set`. This is intended as an internal method for
[`generate_measure_data`](@ref InfiniteOpt.MeasureEvalMethods.generate_measure_data)
and will need to be extended for user-defined
infinite set types. Extensions will also need to consider whether there are
appropriate methods for `method` and extend those as needed.
"""
function generate_supports_and_coeffs(set::InfiniteOpt.AbstractInfiniteSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Function; kwargs...)::Tuple
    error("`generate_supports_and_coeffs` is not extended for parameters in sets " *
          "of type $(typeof(set)).")
end

# IntervalSet
function generate_supports_and_coeffs(set::InfiniteOpt.IntervalSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Function; kwargs...)::Tuple
    return method(lb, ub, num_supports; kwargs...)
end

# DistributionSet
function generate_supports_and_coeffs(set::InfiniteOpt.DistributionSet,
                                      params::Union{InfiniteOpt.ParameterRef,
                                      AbstractArray{<:InfiniteOpt.ParameterRef}},
                                      num_supports::Int,
                                      lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                                      method::Function; kwargs...)::Tuple
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
    return method(dist, params, num_supports; kwargs...)
end

"""
    mc_sampling(lb::Union{JuMPC.SparseAxisArray, Number},
                ub::Union{JuMPC.SparseAxisArray, Number},
                num_supports::Int)::Tuple

Return a tuple that contains supports and coefficients generated by Monte Carlo
sampling from a uniform distribution between the lower and upper bounds provided.

**Example**
```jldoctest; setup = :(using InfiniteOpt, Random; Random.seed!(0))
julia> (supps, coeffs) = mc_sampling(0., 1., 5)
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
function mc_sampling(lb::Number, ub::Number, num_supports::Int)::Tuple
    # MC sampling from uniform distribution over the interval [lb, ub]
    if lb == -Inf || ub == Inf
        return infinite_transform(lb, ub, num_supports)
    else
        samples = rand(num_supports) .* (ub - lb) .+ lb
        return (samples, ones(num_supports) / num_supports * (ub - lb))
    end
end

# MC sampling - multi-dim version
function mc_sampling(lb::JuMPC.SparseAxisArray, ub::JuMPC.SparseAxisArray,
                     num_supports::Int)::Tuple
    samples_dict = Dict()
    for i in eachindex(lb)
        (samples_dict[i], _) = mc_sampling(lb[i], ub[i], num_supports)
    end
    samples = Array{JuMPC.SparseAxisArray, 1}(undef, num_supports)
    for j in 1:num_supports
        samples[j] = JuMP.Containers.SparseAxisArray(Dict(k => samples_dict[k][j]
                                                        for k in eachindex(lb)))
    end
    return (samples, ones(num_supports) / num_supports * prod(ub .- lb))
end

"""
    mc_sampling(dist::Distributions.NonMatrixDistribution,
                param::InfiniteOpt.ParameterRef,
                num_supports::Int)::Tuple

Return a tuple that contains supports and coefficients generated by Monte Carlo
sampling from a given distribution.

**Example**
```jldoctest; setup = :(using InfiniteOpt, Distributions, Random; m = InfiniteModel(seed = true))
julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> @infinite_parameter(m, x in dist)
x

julia> mc_sampling(dist, x, 10)
([0.6791074260357777, 0.8284134829000359, -0.3530074003005963, -0.13485387193052173, 0.5866170746331097, 0.29733585084941616, 0.06494754854834232, -0.10901738508171745, -0.514210390833322, 1.5743302021369892], [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
```
"""
function mc_sampling(dist::Distributions.UnivariateDistribution,
                     param::InfiniteOpt.ParameterRef,
                     num_supports::Int)::Tuple
    # MC sampling - Distribution Set
    samples = rand(dist, num_supports)
    return (samples, ones(num_supports) / num_supports)
end

function mc_sampling(dist::Distributions.MultivariateDistribution,
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

"""
    gauss_legendre(lb::Number, ub::Number, num_supports::Int)::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Legendre quadrature method. This is useful for univariate parameter in a
finite interval.

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> (supps, coeffs) = gauss_legendre(0., 1., 5)
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
function gauss_legendre(lb::Number, ub::Number, num_supports::Int)::Tuple
    # default Gaussian quadrature method for infinite domain
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs)
end

"""
    gauss_hermite(lb::Number, ub::Number, num_supports::Int)::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Hermite quadrature method. This is useful for univariate parameter in an
infinite interval.

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> (supps, coeffs) = gauss_hermite(-Inf, Inf, 5)
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
function gauss_hermite(lb::Number, ub::Number, num_supports::Int)::Tuple
    # default Gaussian quadrature method for infinite domain
    if lb != -Inf || ub != Inf
        error("Lower/upper bound is not infinity. Use other measure evaluation " *
              "methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs)
end

"""
    gauss_laguerre(lb::Number, ub::Number, num_supports::Int)::Tuple

Return a tuple that contains supports and coefficients generated using
Gauss-Laguerre quadrature method. This is useful for univariate parameter in a
semi-infinite interval.

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> (supps, coeffs) = gauss_laguerre(-Inf, 0., 5)
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
function gauss_laguerre(lb::Number, ub::Number, num_supports::Int)::Tuple
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
function infinite_transform(lb::Number, ub::Number, num_supports::Int;
                            sub_method::Function = mc_sampling,
                            transform_x::Function = _default_x,
                            transform_dx::Function = _default_dx,
                            t_lb::Number = -convert(Number, lb == -Inf && ub == Inf),
                            t_ub::Number = 1.)::Tuple
    # transform (semi-)infinite domain to finite domain
    if lb != -Inf && ub != Inf
        error("The range is not (semi-)infinite. Use evaluation methods for " *
              "bounded domains.")
    end
    (t_supports, t_coeffs) = sub_method(t_lb, t_ub, num_supports)
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

# Default method registration
const default_set_types = [InfiniteOpt.IntervalSet, InfiniteOpt.DistributionSet]
const default_methods = [[mc_sampling, gauss_legendre, gauss_laguerre, gauss_hermite],
                         [mc_sampling]]

# TODO: consider truncated distribution
# TODO: consider adding uniform grids
