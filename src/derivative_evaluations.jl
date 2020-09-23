################################################################################
#                        EVALUATE_DERIVIATIVE DEFINITIONS
################################################################################
"""
    evaluate_derivative(dref::GeneralVariableRef, data::AbstractDerivativeMethod,
                        write_model::JuMP.AbstractModel)::PUT_FORMAT_HERE

Add text here. (NOTE maybe something different than `dref` could be given)
"""
function evaluate_derivative end

# TODO add realizations with a fallback too 
# NOTE These should use `make_point_variable_ref` and `make_reduced_variable_ref` as needed

# TODO consider the handling and possible addition of supports for collocation methods

################################################################################
#                              EVALUATION METHODS
################################################################################
"""
    evaluate(dref::DerivativeRef)::PUT_FORMAT_HERE

Add text here.
"""
function evaluate(dref::DerivativeRef) # TODO place output format here 
    # TODO utilize `evaluate_derivative` where `write_model` is the `InfiniteModel` 
    # and the evaluation constraints are added.
end

"""
    evaluate_all_derivatives!(model::InfiniteModel)::Nothing

Add text here.
"""
function evaluate_all_derivatives!(model::InfiniteModel)::Nothing
    # TODO utilize `evaluate_derivative` where `write_model` is the `InfiniteModel` 
    # and the evaluation constraints are added for all the derivatives.
end