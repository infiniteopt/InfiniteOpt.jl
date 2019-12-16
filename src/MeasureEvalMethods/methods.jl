function generate_measure_data(params::Union{InfiniteOpt.ParameterRef,
                               AbstractArray{<:InfiniteOpt.ParameterRef}},
                               num_supports::Int,
                               lb::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing,
                               ub::Union{Number, JuMPC.SparseAxisArray, Nothing} = nothing;
                               method::Function = MC_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w, kwargs...
                               )::InfiniteOpt.AbstractMeasureData
    if isa(params, InfiniteOpt.ParameterRef)
        set = InfiniteOpt._parameter_set(params)
    else
        params = convert(JuMPC.SparseAxisArray, params)
        set = InfiniteOpt._parameter_set(first(params))
    end
    (supports, coeffs) = measure_dispatch(set, params, num_supports, lb, ub, method; kwargs...)
    return InfiniteOpt.DiscreteMeasureData(params, coeffs, supports,
                                           name = name, weight_function = weight_func)
end

function measure_dispatch(set::InfiniteOpt.IntervalSet,
                          params::Union{InfiniteOpt.ParameterRef,
                          AbstractArray{<:InfiniteOpt.ParameterRef}},
                          num_supports::Int,
                          lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                          ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                          method::Function; kwargs...)::Tuple
    return method(lb, ub, num_supports; kwargs...)
end

function measure_dispatch(set::InfiniteOpt.DistributionSet,
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

# fallback
function measure_dispatch(set::InfiniteOpt.AbstractInfiniteSet,
                          params::Union{InfiniteOpt.ParameterRef,
                          AbstractArray{<:InfiniteOpt.ParameterRef}},
                          num_supports::Int,
                          lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
                          ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
                          method::Function; kwargs...)::Tuple
    error("Measure dispatch function is not extended for parameters in sets " *
          "of type $(typeof(set)).")
    return
end

# MC sampling from uniform distribution over the interval [lb, ub]
function MC_sampling(lb::Number, ub::Number, num_supports::Int; kwargs...)::Tuple
    if lb == -Inf || ub == Inf
        return infinite_transform(lb, ub, num_supports)
    else
        samples = rand(num_supports) .* (ub - lb) .+ lb
        return (samples, ones(num_supports) / num_supports * (ub - lb))
    end
end

# MC sampling - multi-dim version
function MC_sampling(lb::JuMPC.SparseAxisArray, ub::JuMPC.SparseAxisArray,
                     num_supports::Int; kwargs...)::Tuple
    samples_dict = Dict()
    for i in eachindex(lb)
        (samples_dict[i], _) = MC_sampling(lb[i], ub[i], num_supports)
    end
    samples = Array{JuMPC.SparseAxisArray, 1}(undef, num_supports)
    for j in 1:num_supports
        samples[j] = JuMP.Containers.SparseAxisArray(Dict(k => samples_dict[k][j]
                                                        for k in eachindex(lb)))
    end
    return (samples, ones(num_supports) / num_supports * prod(ub .- lb))
end

# MC sampling - Distribution Set
function MC_sampling(dist::Distributions.UnivariateDistribution,
                     param::InfiniteOpt.ParameterRef,
                     num_supports::Int; kwargs...)::Tuple
    samples = rand(dist, num_supports)
    return (samples, ones(num_supports) / num_supports)
end

function MC_sampling(dist::Distributions.MultivariateDistribution,
                     params::AbstractArray{<:InfiniteOpt.ParameterRef},
                     num_supports::Int; kwargs...)::Tuple
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
# default Gaussian quadrature method for bounded domain
function Gauss_Legendre(lb::Number, ub::Number, num_supports::Int; kwargs...)::Tuple
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs)
end

# default Gaussian quadrature method for infinite domain
function Gauss_Hermite(lb::Number, ub::Number, num_supports::Int; kwargs...)::Tuple
    if lb != -Inf || ub != Inf
        error("Lower/upper bound is not infinity. Use other measure evaluation " *
              "methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs)
end

# default Gaussian quadrature method for semi-infinite domain
function Gauss_Laguerre(lb::Number, ub::Number, num_supports::Int; kwargs...)::Tuple
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

# transform (semi-)infinite domain to finite domain
function infinite_transform(lb::Number, ub::Number, num_supports::Int;
                            sub_method::Function = MC_sampling,
                            transform_x::Function = _default_x,
                            transform_dx::Function = _default_dx,
                            t_lb::Number = -convert(Number, lb == -Inf && ub == Inf),
                            t_ub::Number = 1., kwargs...)::Tuple
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

# TODO: consider truncated distribution
# TODO: consider adding uniform grids
