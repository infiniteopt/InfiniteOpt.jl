## A template for defining a new derivative evaluation method 

# Define the method struct inheriting from either `GenerativeDerivativeMethod` or `NonGenerativeDerivativeMethod`
# Here we'll use `GenerativeDerivativeMethod` assuming our method will add extra supports
# This should contain any information that is needed to implement the method
struct MyDerivMethod <: GenerativeDerivativeMethod
    my_attr::Bool # REPLACE WITH ACTUAL ATTRIBUTE
    # ADD ANY MORE INFORMATION THAT IS NEEDED
end

# Extend `support_label` (only needed for generative methods)
function InfiniteOpt.support_label(method::MyDerivMethod)::DataType 
    return InternalLabel # REPLACE WITH DESIRED SUPPORT LABEL FOR ANY EXTRA SUPPORTS THAT ARE ADDED
end

# Extend `generate_derivative_supports` (only needed for generative methods)
function InfiniteOpt.generate_derivative_supports(
    pref::IndependentParameterRef,
    method::MyDerivMethod
    )::Vector{Float64}
    # collect the supports already added to `pref`
    curr_supps = supports(pref, label = All)
    # generate the extra supports and return them (REPLACE BELOW WITH ACTUAL)
    min_s = minimum(curr_supps)
    max_s = maximum(curr_supps)
    return rand(2) * (max_s - min_s) .+ min_s
end

# Extend `evaluate_derivative`
# Generative methods must include calls of `add_derivative_supports`
# It will likely also be convenient to use `make_reduced_expr`
function InfiniteOpt.evaluate_derivative(
    dref::GeneralVariableRef, 
    method::MyDerivMethod,
    write_model::JuMP.AbstractModel
    )::Vector{JuMP.AbstractJuMPScalar}
    # get the basic derivative information 
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    # make sure internal supports are added to the model
    InfiniteOpt.add_derivative_supports(pref)
    # generate the derivative expressions h_i corresponding to the equations of 
    # the form h_i = 0 (REPLACE BELOW WITH ACTUAL IMPLEMENTATION)
    supps = supports(pref, label = All)
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(supps) - 1)
    for i in eachindex(exprs)
        d = InfiniteOpt.make_reduced_expr(dref, pref, supps[i], write_model)
        v1 = InfiniteOpt.make_reduced_expr(vref, pref, supps[i], write_model)
        v2 = InfiniteOpt.make_reduced_expr(vref, pref, supps[i + 1], write_model)
        change = supps[i + 1] - supps[i]
        exprs[i] = JuMP.@expression(write_model, -change * d + v2 - v1)
    end
    return exprs
end
