# Define the new measure evaluation method
# myNewEval for univariate IntervalSet
const NewEvalMethod = :NewEvalMethod
function InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs(
    set::InfiniteOpt.IntervalSet,
    params::Union{InfiniteOpt.ParameterRef,
    AbstractArray{<:InfiniteOpt.ParameterRef}},
    num_supports::Int,
    lb::Number,
    ub::Number,
    method::Val{NewEvalMethod})::Tuple
    # REPLACE WITH ACTUAL FUNCTIONALITY
    increment = (ub - lb) / (num_supports - 1)
    supports = [lb + (i - 1) * increment for i in 1:num_supports]
    # MAKE SURE IT RETURNS APPROPRIATE DATA IN THE TUPLE
    return (supports, ones(num_supports) / num_supports * (ub - lb))
end

function InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs(
    set::InfiniteOpt.IntervalSet,
    params::Union{InfiniteOpt.ParameterRef,
    AbstractArray{<:InfiniteOpt.ParameterRef}},
    num_supports::Int,
    lb::JuMPC.SparseAxisArray,
    ub::JuMPC.SparseAxisArray,
    method::Val{NewEvalMethod}; independent::Bool = true)::Tuple # REPLACE WITH ACTUAL ARGUMENTS
    # REPLACE WITH ACTUAL FUNCTIONALITY
    if !independent
        @warn("The method is implemented for independent multivariate parameters.")
    end
    supports_dict = Dict()
    for i in eachindex(lb)
        (supports_dict[i], _) = generate_supports_and_coeffs(set, params,
                                 num_supports, lb[i], ub[i], Val(NewEvalMethod))
    end
    supports = Array{JuMPC.SparseAxisArray, 1}(undef, num_supports)
    for j in 1:num_supports
        supports[j] = JuMP.Containers.SparseAxisArray(Dict(k => supports_dict[k][j]
                                                        for k in eachindex(lb)))
    end
    # MAKE SURE IT RETURNS APPROPRIATE DATA IN THE TUPLE
    return (supports, ones(num_supports) / num_supports * prod(ub .- lb))
end
