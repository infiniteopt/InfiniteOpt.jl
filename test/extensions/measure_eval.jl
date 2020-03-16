# Define the new measure evaluation method
# myNewEval for univariate IntervalSet
function NewEvalMethod(lb::Number, ub::Number, num_supports::Int)::Tuple # REPLACE WITH ACTUAL ARGUMENTS
    # REPLACE WITH ACTUAL FUNCTIONALITY
    increment = (ub - lb) / (num_supports - 1)
    supports = [lb + (i - 1) * increment for i in 1:num_supports]
    # MAKE SURE IT RETURNS APPROPRIATE DATA IN THE TUPLE
    return (supports, ones(num_supports) / num_supports * (ub - lb))
end

# myNewEval multivariate IntervalSet
function NewEvalMethod(lb::JuMPC.SparseAxisArray, ub::JuMPC.SparseAxisArray,
                       num_supports::Int; independent::Bool = true)::Tuple # REPLACE WITH ACTUAL ARGUMENTS
    # REPLACE WITH ACTUAL FUNCTIONALITY
    if !independent
        @warn("The method is implemented for independent multivariate parameters.")
    end
    supports_dict = Dict()
    for i in eachindex(lb)
        (supports_dict[i], _) = NewEvalMethod(lb[i], ub[i], num_supports)
    end
    supports = Array{JuMPC.SparseAxisArray, 1}(undef, num_supports)
    for j in 1:num_supports
        supports[j] = JuMP.Containers.SparseAxisArray(Dict(k => supports_dict[k][j]
                                                        for k in eachindex(lb)))
    end
    # MAKE SURE IT RETURNS APPROPRIATE DATA IN THE TUPLE
    return (supports, ones(num_supports) / num_supports * prod(ub .- lb))
end
