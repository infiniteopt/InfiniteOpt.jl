## A template for defining a new derivative evaluation method 

# Define the method struct inheriting from either `GenerativeDerivativeMethod` or `NonGenerativeDerivativeMethod`
# Here we'll use `GenerativeDerivativeMethod` assuming our method will add extra supports
# This should contain any information that is needed to implement the method
struct MyDerivMethod <: GenerativeDerivativeMethod
    my_attr::Float64 # REPLACE WITH ACTUAL ATTRIBUTE
    # ADD ANY MORE INFORMATION THAT IS NEEDED
end

# Extend `generative_support_info` (only needed for generative methods)
function InfiniteOpt.generative_support_info(method::MyDerivMethod)
    info = UniformGenerativeInfo([method.my_attr], InternalLabel)
    return info # REPLACE WITH NEEDED AbstractGenerativeInfo THAT IS NEEDED TO MAKE THE GENERATIVE SUPPORTS
end

# Extend `evaluate_derivative`
# Generative methods must include calling `add_generative_supports`
# It will likely also be convenient to use `make_reduced_expr`
function InfiniteOpt.evaluate_derivative(
    dref::GeneralVariableRef, 
    method::MyDerivMethod,
    write_model::JuMP.AbstractModel
    )::Vector{JuMP.AbstractJuMPScalar}
    # get the basic derivative information 
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    # make sure generative supports are added to the model
    InfiniteOpt.add_generative_supports(pref)
    # generate the derivative expressions h_i corresponding to equations of 
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
