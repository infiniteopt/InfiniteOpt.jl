## Define the new measure evaluation method
# Make alias for our new method
struct NewUniEvalMethod <: InfiniteOpt.MeasureToolbox.AbstractUnivariateMethod end
struct NewMultiEvalMethod <: InfiniteOpt.MeasureToolbox.AbstractMultivariateMethod end

# Extend generate_supports_and_coeffs for scalar params
function InfiniteOpt.MeasureToolbox.generate_integral_data(
    pref::InfiniteOpt.GeneralVariableRef,
    lower_bound::Real,
    upper_bound::Real,
    method::Type{NewUniEvalMethod};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function
    )::InfiniteOpt.AbstractMeasureData # REPLACE WITH ACTUAL ALIAS
    # REPLACE WITH ACTUAL FUNCTIONALITY
    # ADD CHECKS IF NECESSARY
    increment = (upper_bound - lower_bound) / (num_supports - 1)
    supports = [lower_bound + (i - 1) * increment for i in 1:num_supports]
    coeffs = ones(num_supports) / num_supports * (upper_bound - lower_bound)
    # MAKE SURE IT RETURNS APPROPRIATE MEASURE DATA
    return InfiniteOpt.DiscreteMeasureData(
        pref, coeffs, supports,
        weight_function = weight_func,
        lower_bound = lower_bound, 
        upper_bound = upper_bound)
end

# Extend generate_supports_and_coeffs for multi-dimensional params
function InfiniteOpt.MeasureToolbox.generate_integral_data(
    prefs::Vector{InfiniteOpt.GeneralVariableRef},
    lower_bounds::Vector{<:Real},
    upper_bounds::Vector{<:Real},
    method::Type{NewMultiEvalMethod};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function
    )::InfiniteOpt.AbstractMeasureData # REPLACE WITH ACTUAL ARGUMENTS
    # REPLACE WITH ACTUAL FUNCTIONALITY
    # ADD CHECKS IF NECESSARY
    if !all(i -> (dispatch_variable_ref(i) isa InfiniteOpt.IndependentParameterRef), prefs)
        @warn("The method is implemented for independent multivariate parameters.")
    end
    supports_dict = Dict()
    for i in eachindex(lower_bounds)
        increment = (upper_bounds[i] - lower_bounds[i]) / (num_supports - 1)
        supports_dict[i] = [lower_bounds[i] + (j - 1) * increment for j in 1:num_supports]
    end
    supports = Array{Vector{Float64}, 1}(undef, num_supports)
    for j in 1:num_supports
        supports[j] = [supports_dict[k][j] for k in keys(lower_bounds)]
    end
    coeffs = ones(num_supports) / num_supports * prod(upper_bounds .- lower_bounds)
    # MAKE SURE IT RETURNS APPROPRIATE MEASURE DATA
    return InfiniteOpt.DiscreteMeasureData(
        prefs, coeffs, supports,
        weight_function = weight_func,
        lower_bounds = lower_bounds, 
        upper_bounds = upper_bounds)
end
