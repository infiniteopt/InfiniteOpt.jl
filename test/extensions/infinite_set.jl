# Define the new set structure
struct MyNewSet <: InfiniteOpt.AbstractInfiniteSet
    attr1::Float64 # REPLACE WITH DESIRED TYPE
    attr2::Float64 # REPLACE WITH DESIRED TYPE
    # ADD MORE ATTRIBUTES AS NEEDED
    # constructor
    function MyNewSet(attr1, attr2)
        # PUT APPROPRIATE CHECKS HERE AND OTHER NEEDED OPERATIONS
        return new(convert(Float64, attr1), convert(Float64, attr2))
    end
end

# Extend supports_in_set
function InfiniteOpt.supports_in_set(supports::Union{Number, Vector{<:Number}},
                                     set::MyNewSet)::Bool
    # DETERMINE IF SUPPORTS ARE IN THE DOMAIN OF `set`
    in_set = set.attr1 <= minimum(supports) && maximum(supports) <= set.attr2 # REPLACE WITH ACTUAL CHECK
    return in_set
end

# Extend generate_support_values if possible
function InfiniteOpt.generate_support_values(set::MyNewSet;
                                            num_supports::Int = 50,
                                            sig_fig::Int = 5)::Array
    # REPLACE BELOW WITH METHODS TO GENERATE `num_samples` with `sig_fig`
    supports = collect(range(set.attr1, stop = set.attr2, length = num_supports))
    return round.(supports, sigdigits = sig_fig)
end

# Extend generate_and_add_supports! if supports are not generated for parameters individually
# Here we'll assume generate_support_values is sufficient

# Extend measure_dispatch in enable measure() with methods
function InfiniteOpt.MeasureEvalMethods.measure_dispatch(
    set::MyNewSet,
    params::Union{InfiniteOpt.ParameterRef, AbstractArray{<:InfiniteOpt.ParameterRef}},
    num_supports::Int,
    lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
    ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
    method::Function; kwargs...
    )::Tuple
    # ADD CHECKS IF NEEDED
    return method(lb, ub, num_supports; kwargs...) # REPLACE WITH PROPER DISPATCH
end

# Extend JuMP.has_lower_bound (optional if the answer is always false)
function JuMP.has_lower_bound(set::MyNewSet)::Bool
    # INSERT NECESSARY CHECKS IF NEEDED
    has_bound = true # REPLACE WITH ACTUAL RESULT
    return has_bound
end

# Extend JuMP.lower_bound if appropriate
function JuMP.lower_bound(set::MyNewSet)
    return set.attr1 # REPLACE WITH ACTUAL
end

# Extend JuMP.set_lower_bound if appropriate
function JuMP.set_lower_bound(set::MyNewSet, lower::Number)::MyNewSet
    return MyNewSet(lower, set.attr2) # REPLACE WITH ACTUAL CONSTRUCTOR
end

# Extend JuMP.has_supper_bound (optional if the answer is always false)
function JuMP.has_upper_bound(set::MyNewSet)::Bool
    # INSERT NECESSARY CHECKS IF NEEDED
    has_bound = true # REPLACE WITH ACTUAL RESULT
    return has_bound
end

# Extend JuMP.lower_bound if appropriate
function JuMP.upper_bound(set::MyNewSet)
    return set.attr2 # REPLACE WITH ACTUAL
end

# Extend JuMP.set_upper_bound if appropriate
function JuMP.set_upper_bound(set::MyNewSet, upper::Number)::MyNewSet
    return MyNewSet(set.attr1, upper) # REPLACE WITH ACTUAL CONSTRUCTOR
end
