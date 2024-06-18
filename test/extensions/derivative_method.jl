## A template for defining a new derivative evaluation method 

# Define the method struct inheriting from either `GenerativeDerivativeMethod` or `NonGenerativeDerivativeMethod`
# Here we'll use `GenerativeDerivativeMethod` assuming our method will add extra supports
# This should contain any information that is needed to implement the method
struct MyDerivMethod <: GenerativeDerivativeMethod
    my_attr::Float64 # REPLACE WITH ACTUAL ATTRIBUTE
    # ADD ANY MORE INFORMATION THAT IS NEEDED
end

# Extend `allows_high_order_derivatives`
# Return a `Bool` on whether this method will explictly support derivatives 
# with an order greater than 1 (e.g., 2nd derivatives)
# If we return `false`, the InfiniteOpt will automatically reformulate higher order 
# derivatives into 1st derivatives that the method can handle
InfiniteOpt.allows_high_order_derivatives(method::MyDerivMethod) = false # TODO replace with desired output

# Extend `generative_support_info` (only needed for generative methods)
function InfiniteOpt.generative_support_info(method::MyDerivMethod)
    info = UniformGenerativeInfo([method.my_attr], InternalLabel)
    return info # REPLACE WITH NEEDED AbstractGenerativeInfo THAT IS NEEDED TO MAKE THE GENERATIVE SUPPORTS
end

# Extend `derivative_expr_data` (see docs for details)
function InfiniteOpt.derivative_expr_data(
    dref::GeneralVariableRef, # the dervative
    order::Int, # the derivative order
    supps::Vector{Float64}, # the ordered supports of infinite parameter
    method::MyDerivMethod
    )
    # generate the support indices to be used for each call of `make_indexed_derivative_expr`
    idxs = 1:length(supps)-1 # TODO REPLACE WITH ACTUAL IMPLEMENTATION
    # generate any additional data that is needed as interators
    supp_diffs = (supps[i + 1] - supps[i] for i in idxs) # TODO REPLACE WITH ACTUAL IMPLEMENTATION
    # return the indexes and the other iterators
    return idxs, supp_diffs # TODO MODIFY AS NEEDED BASED ON ABOVE
end

# Extend `make_indexed_derivative_expr` (see docs for details)
function InfiniteOpt.make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::JuMP.AbstractModel,
    ::MyDerivMethod,
    supp_diff
    )
    # generate the derivative expression h corresponding to the equation of 
    # the form h = 0 (REPLACE BELOW WITH ACTUAL IMPLEMENTATION)
    d = InfiniteOpt.make_reduced_expr(dref, pref, supps[idx], write_model)
    v1 = InfiniteOpt.make_reduced_expr(vref, pref, supps[idx], write_model)
    v2 = InfiniteOpt.make_reduced_expr(vref, pref, supps[idx + 1], write_model)
    return JuMP.@expression(write_model, -supp_diff * d + v2 - v1)
end
