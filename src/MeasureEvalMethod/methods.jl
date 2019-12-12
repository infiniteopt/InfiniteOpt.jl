function generate_measure_data(param::InfiniteOpt.ParameterRef,
                               lb::Union{Number, Nothing},
                               ub::Union{Number, Nothing}, num_supports::Int;
                               method::Function = MC_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w
                               )::InfiniteOpt.DiscreteMeasureData
    set = InfiniteOpt._parameter_set(param)
    if isa(set, DistributionSet)
        (supports, coeffs) = method(set.distribution, num_supports)
    else
        (supports, coeffs) = method(lb, ub, num_supports)
    end
    return InfiniteOpt.DiscreteMeasureData(param, coeffs, supports,
                                           name = name, weight_function = weight_func)
end

#function generate_measure_data(params::JuMPC.SparseAxisArray{<:InfiniteOpt.ParameterRef},
function generate_measure_data(params::AbstractArray{<:InfiniteOpt.ParameterRef},
                               lb::Union{JuMPC.SparseAxisArray, Nothing},
                               ub::Union{JuMPC.SparseAxisArray, Nothing},
                               num_supports::Int;
                               method::Function = MC_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w
                               )::InfiniteOpt.MultiDiscreteMeasureData
    params = convert(JuMPC.SparseAxisArray, params)
    set = 1;
    for i in params
        set = InfiniteOpt._parameter_set(i)
        break
    end
    if isa(set, DistributionSet)
        (supports, coeffs) = method(set.distribution, params, num_supports)
    else
        (supports, coeffs) = method(lb, ub, num_supports)
    end
    return InfiniteOpt.DiscreteMeasureData(params, coeffs, supports,
                                           name = name, weight_function = weight_func)
end

# MC sampling from uniform distribution over the interval [lb, ub]
function MC_sampling(lb::Number, ub::Number, num_supports::Int)::Tuple
    if lb == -Inf || ub == Inf
        return infinite_transform(lb, ub, num_supports)
    else
        samples = rand(num_supports) .* (ub - lb) .+ lb
        return (samples, ones(num_supports) / num_supports * (ub - lb))
    end
end

# MC sampling - multi-dim version
function MC_sampling(lb::JuMPC.SparseAxisArray, ub::JuMPC.SparseAxisArray,
                     num_supports::Int)::Tuple
    samples_dict = Dict()
    for i in eachindex(lb)
        (samples_dict[i], _) = MC_sampling(lb[i], ub[i], num_supports)
    end
    samples = Array{JuMPC.SparseAxisArray, 1}()
    for j in 1:num_supports
        append!(samples, [JuMP.Containers.SparseAxisArray(Dict(k => samples_dict[k][j]
                          for k in eachindex(lb)))])
    end
    return (samples, ones(num_supports) / num_supports * prod(ub .- lb))

end

# MC sampling - Distribution Set
function MC_sampling(dist::Distributions.UnivariateDistribution,
                     num_supports::Int)::Tuple
    samples = rand(dist, num_supports)
    return (samples, ones(num_supports) / num_supports)
end

function MC_sampling(dist::Distributions.MultivariateDistribution,
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
# default Gaussian quadrature method for bounded domain
function Gauss_Legendre(lb::Number, ub::Number, num_supports::Int)::Tuple
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs)
end

# default Gaussian quadrature method for infinite domain
function Gauss_Hermite(lb::Number, ub::Number, num_supports::Int)::Tuple
    if lb != -Inf || ub != Inf
        error("Lower/upper bound is not infinity. Use other measure evaluation " *
              "methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs)
end

# default Gaussian quadrature method for semi-infinite domain
function Gauss_Laguerre(lb::Number, ub::Number, num_supports::Int)::Tuple
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
                            t_ub::Number = 1.)::Tuple
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
