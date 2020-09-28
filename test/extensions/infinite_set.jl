# Define the new set structure
struct MyNewSet <: InfiniteOpt.InfiniteScalarSet
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
function InfiniteOpt.supports_in_set(supports::Union{Real, Vector{<:Real}},
                                     set::MyNewSet)::Bool
    # DETERMINE IF SUPPORTS ARE IN THE DOMAIN OF `set`
    in_set = all(set.attr1 .<= supports .<= set.attr2) # REPLACE WITH ACTUAL CHECK
    return in_set
end

# Extend generate_support_values if possible
function InfiniteOpt.generate_support_values(set::MyNewSet;
                                             num_supports::Int = 10,
                                             sig_digits::Int = 5)::Tuple{Vector{<:Real}, DataType}
    # REPLACE BELOW WITH METHODS TO GENERATE `num_samples` with `sig_fig` 
    supports = collect(range(set.attr1, stop = set.attr2, length = num_supports))
    return round.(supports, sigdigits = sig_digits), UniformGrid
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
function JuMP.set_lower_bound(set::MyNewSet, lower::Real)::MyNewSet
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
function JuMP.set_upper_bound(set::MyNewSet, upper::Real)::MyNewSet
    return MyNewSet(set.attr1, upper) # REPLACE WITH ACTUAL CONSTRUCTOR
end
