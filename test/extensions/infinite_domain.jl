# Define the new domain structure
struct MyNewDomain <: InfiniteOpt.InfiniteScalarDomain
    attr1::Float64 # REPLACE WITH DESIRED TYPE
    attr2::Float64 # REPLACE WITH DESIRED TYPE
    # ADD MORE ATTRIBUTES AS NEEDED
    # constructor
    function MyNewDomain(attr1, attr2)
        # PUT APPROPRIATE CHECKS HERE AND OTHER NEEDED OPERATIONS
        return new(convert(Float64, attr1), convert(Float64, attr2))
    end
end

# Extend supports_in_domain
function InfiniteOpt.supports_in_domain(
    supports::Union{Real, Vector{<:Real}},
    domain::MyNewDomain
    )::Bool
    # DETERMINE IF SUPPORTS ARE IN THE DOMAIN OF `domain`
    in_domain = all(domain.attr1 .<= supports .<= domain.attr2) # REPLACE WITH ACTUAL CHECK
    return in_domain
end

# Extend generate_support_values if possible
function InfiniteOpt.generate_support_values(
    domain::MyNewDomain;
    num_supports::Int = 10,
    sig_digits::Int = 5
    )::Tuple{Vector{Float64}, DataType}
    # REPLACE BELOW WITH METHODS TO GENERATE `num_samples` with `sig_fig` 
    supports = collect(range(domain.attr1, stop = domain.attr2, length = num_supports))
    return round.(supports, sigdigits = sig_digits), UniformGrid
end

# Extend JuMP.has_lower_bound (optional if the answer is always false)
function JuMP.has_lower_bound(domain::MyNewDomain)::Bool
    # INSERT NECESSARY CHECKS IF NEEDED
    has_bound = true # REPLACE WITH ACTUAL RESULT
    return has_bound
end

# Extend JuMP.lower_bound if appropriate
function JuMP.lower_bound(domain::MyNewDomain)
    return domain.attr1 # REPLACE WITH ACTUAL
end

# Extend JuMP.set_lower_bound if appropriate
function JuMP.set_lower_bound(domain::MyNewDomain, lower::Real)::MyNewDomain
    return MyNewDomain(lower, domain.attr2) # REPLACE WITH ACTUAL CONSTRUCTOR
end

# Extend JuMP.has_upper_bound (optional if the answer is always false)
function JuMP.has_upper_bound(domain::MyNewDomain)::Bool
    # INSERT NECESSARY CHECKS IF NEEDED
    has_bound = true # REPLACE WITH ACTUAL RESULT
    return has_bound
end

# Extend JuMP.lower_bound if appropriate
function JuMP.upper_bound(domain::MyNewDomain)
    return domain.attr2 # REPLACE WITH ACTUAL
end

# Extend JuMP.set_upper_bound if appropriate
function JuMP.set_upper_bound(domain::MyNewDomain, upper::Real)::MyNewDomain
    return MyNewDomain(domain.attr1, upper) # REPLACE WITH ACTUAL CONSTRUCTOR
end

# Extend InfiniteOpt.MeasureToolbox.generate_expect_data if wanted
function InfiniteOpt.MeasureToolbox.generate_expect_data(
    domain::MyNewDomain, 
    pref::GeneralVariableRef, 
    num_supports::Int; 
    kwargs...
    )
    new_domain = IntervalDomain(domain.attr1, domain.attr2)
    return generate_expect_data(new_domain, pref, num_supports; kwargs...) # REPLACE WITH ACTUAL
end
