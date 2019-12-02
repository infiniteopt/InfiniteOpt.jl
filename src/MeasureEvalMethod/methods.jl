function generate_measure_data(param::InfiniteOpt.ParameterRef,
                               lb::Float64, ub::Float64, num_supports::Int;
                               method::Function = MC_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w
                               )::InfiniteOpt.DiscreteMeasureData
    (supports, coeffs) = method(lb, ub, num_supports)
    return InfiniteOpt.DiscreteMeasureData(param, coeffs, supports, name, weight_func)
end

function generate_measure_data(param::Vector{InfiniteOpt.ParameterRef},
                               lb::Vector{Float64}, ub::Vector{Float64},
                               num_supports::Int;
                               method::Function = MC_sampling, name::String = "",
                               weight_func::Function = InfiniteOpt._w
                               )::InfiniteOpt.MultiDiscreteMeasureData
    (supports, coeffs) = method(lb, ub, num_supports)
    return InfiniteOpt.MultiDiscreteMeasureData(param, coeffs, supports, name, weight_func)
end

# MC sampling from uniform distribution over the interval [lb, ub]
function MC_sampling(lb::Float64, ub::Float64, num_supports::Int)::Tuple
    if lb == -Inf || ub == Inf
        (samples, _) = infinite_transform(lb, ub, num_supports)
    else
        samples = rand(num_supports) .* (ub - lb) .+ lb
    end
    return (samples, ones(num_supports) / num_supports * (ub - lb))
end

# MC sampling - multi-dim version
function MC_sampling(lb::Vector{Float64}, ub::Vector{Float64},
                     num_supports::Int)::Tuple
    samples = zeros(length(lb), num_supports)
    for i in eachindex(lb)
        (samples[i, :], _) = MC_sampling(lb[i], ub[i], num_supports)
    end
    return (samples, ones(num_supports) / num_supports * prod(ub .- lb))
end

# default Gaussian quadrature method for bounded domain
function Gauss_Legendre(lb::Float64, ub::Float64, num_supports::Int)::Tuple
    (supports, coeffs) = FastGaussQuadrature.gausslegendre(num_supports)
    supports = (ub - lb) / 2 * supports .+ (ub + lb) / 2
    coeffs = (ub - lb) / 2 * coeffs
    return (supports, coeffs)
end

# default Gaussian quadrature method for infinite domain
function Gauss_Hermite(lb::Float64, ub::Float64, num_supports::Int)::Tuple
    if lb != -Inf || ub != Inf
        error("Lower/upper bound is not infinity. Use other measure evaluation " *
              "methods.")
    end
    (supports, coeffs) = FastGaussQuadrature.gausshermite(num_supports)
    coeffs = coeffs .* exp.(supports.^2)
    return (supports, coeffs)
end

# default Gaussian quadrature method for semi-infinite domain
function Gauss_Laguerre(lb::Float64, ub::Float64, num_supports::Int)::Tuple
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
function infinite_transform(lb::Float64, ub::Float64, num_supports::Int;
                            sub_method::Function = MC_sampling,
                            transform_x::Function = _default_x,
                            transform_dx::Function = _default_dx,
                            t_lb::Float64 = -convert(Float64, lb == -Inf && ub == Inf),
                            t_ub::Float64 = 1)::Tuple
    if lb != -Inf || ub != Inf
        error("The range is not (semi-)infinite. Use evaluation methods for " *
              "bounded domains.")
    end
    (t_supports, t_coeffs) = sub_method(t_lb, t_ub, num_supports)
    supports = transform_x.(t_supports, lb, ub)
    coeffs = t_coeffs .* transform_dx.(t_supports, lb, ub)
    return (supports, coeffs)
end

function _default_x(t::Float64, lb::Float64, ub::Float64)::Float64
    if lb > -Inf
        return lb + t / (1 - t)
    elseif ub < Inf
        return ub - (1 - t) / t
    else
        return t / (1 - t^2)
    end
end

function _default_dx(t::Float64, lb::Float64, ub::Float64)::Float64
    if lb > -Inf
        return 1 / (1 - t)^2
    elseif ub < Inf
        return 1 / t^2
    else
        return (1 + t^2) / (1 - t^2)^2
    end
end

# TODO: consider truncated distribution
