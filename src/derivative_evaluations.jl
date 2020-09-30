################################################################################
#                                HELPER METHODS
################################################################################
"""
    generate_derivative_supports(pref::IndependentParameterRef, 
                                 method::GenerativeDerivativeMethod)::Vector{Float64}

Generate and return a vector any additional supports needed by `method`. This is 
intended as an internal method and will need to be extended for user-defined 
derivative methods that are generative.
"""
function generate_derivative_supports(
    pref::IndependentParameterRef, 
    method::AbstractDerivativeMethod
    )
    error("`generate_derivative_supports` not extended for derivative method of " * 
          "type $(typeof(method)).")
end

# NonGenerativeDerivativeMethod
function generate_derivative_supports(
    pref::IndependentParameterRef, 
    method::NonGenerativeDerivativeMethod
    )::Vector{Float64}
    return Float64[]
end

# OrthogonalCollocation
function generate_derivative_supports(
    pref::IndependentParameterRef, 
    method::OrthogonalCollocation
    )::Vector{Float64}
    # TODO GENERATE SUPPORTS AND RETURN THEM
end

"""
    add_derivative_supports(pref::Union{IndependentParameterRef, DependentParameterRef})::Nothing

Add any supports `pref` that are needed for derivative evaluation. This is intended 
as a helper method for derivative evaluation and depends [`generate_derivative_supports`](@ref InfiniteOpt.generate_derivative_supports) 
which will need to be extended for user-defined derivative methods that generate supports. 
In such cases, it is necessary to also extend 
[`support_label`](@ref InfiniteOpt.support_label(::AbstractDerivativeMethod)) Errors if 
such is not defined for the current derivative method associated with `pref`. 
"""
function add_derivative_supports(pref::IndependentParameterRef)::Nothing 
    if !has_derivative_supports(pref)
        method = derivative_method(pref)
        supps = generate_derivative_supports(pref, method)
        if !isempty(supps)
            add_supports(pref, supps, label = support_label(method))
            set_has_derivative_supports(pref, true)
        end
    end
    return
end

# Define for DependentParameterRef
function add_derivative_supports(pref::DependentParameterRef)::Nothing
    return 
end

"""
    make_reduced_expr(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                      support::Float64, write_model::Union{InfiniteModel, JuMP.Model})

Given the argument variable `vref` and the operator parameter `pref` from a 
derivative, build and return the reduced expression in accordance to the support 
`support` with respect to `pref`. New point/reduced variables will be written to 
`write_model`. This is solely intended as a helper function for derivative 
evaluation.
"""
function make_reduced_expr(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                           support::Float64, write_model::JuMP.AbstractModel)
    return make_reduced_expr(vref, _index_type(vref), pref, support, write_model)
end

# MeasureIndex
function make_reduced_expr(mref, ::Type{MeasureIndex}, pref, support, write_model)
    data = DiscreteMeasureData(pref, [1], [support], InternalLabel, # NOTE the label will not be used
                               default_weight, NaN, NaN, false)
    return expand_measure(mref, data, write_model)
end

# InfiniteVariableIndex
function make_reduced_expr(vref, 
    ::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}}, 
    pref, support, write_model
    )::GeneralVariableRef
    prefs = parameter_list(vref)
    if length(prefs) == 1
        return make_point_variable_ref(write_model, vref, [support])
    else 
        pindex = findfirst(isequal(pref), prefs)
        return make_reduced_variable_ref(write_model, vref, [pindex], [support])
    end
end

# ReducedVariableIndex
function make_reduced_expr(vref, ::Type{ReducedVariableIndex}, pref, support, 
                           write_model)::GeneralVariableRef
    dvref = dispatch_variable_ref(vref)
    ivref = infinite_variable_ref(vref)
    var_prefs = parameter_list(dvref)
    orig_prefs = parameter_list(ivref)
    eval_supps = eval_supports(dvref)
    pindex = findfirst(isequal(pref), orig_prefs)
    if length(var_prefs) == 1
        processed_support = _make_point_support(orig_prefs, eval_supps, pindex, support)
        return make_point_variable_ref(write_model, ivref, processed_support)
    else
        indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], indices)
        push!(indices, pindex)
        push!(vals, support)
        return make_reduced_variable_ref(write_model, ivref, indices, vals)
    end
end

################################################################################
#                        EVALUATE_DERIVIATIVE DEFINITIONS
################################################################################
"""
    evaluate_derivative(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                        method::AbstractDerivativeMethod,
                        write_model::JuMP.AbstractModel)::PUT_FORMAT_HERE

Build expression for derivative of `vref` with respect to `pref` evaluated at each
support of `pref` using finite difference scheme.
"""
function evaluate_derivative(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                             method::AbstractDerivativeMethod,
                             write_model::JuMP.AbstractModel)
    error("`evaluate_derivative` not defined for derivative method of type " * 
          "$(typeof(method)).")
end

"""
    evaluate_derivative(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                        method::FiniteDifference,
                        write_model::JuMP.AbstractModel)::PUT_FORMAT_HERE
    
Build expression for derivative of `vref` with respect to `pref` evaluated at each
support of `pref` using finite difference scheme.
"""
function evaluate_derivative(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                             method::FiniteDifference, 
                             write_model::JuMP.AbstractModel)

    n_supps = num_supports(pref)
    if n_supps <= 1
        error("$(pref) does not have enough supports to apply any finite difference methods.")
    end
    ordered_supps = sort(supports(pref))
    expr_dict = Dict{Float64, JuMP.AbstractJuMPScalar}()
    for i in eachindex(ordered_supps)
        curr_value = ordered_supps[i]
        if i == 1
            expr_dict[curr_value] = _make_difference_expr(vref, pref, i, ordered_supps, write_model, FDForward)
        elseif i == n_supps
            expr_dict[curr_value] = _make_difference_expr(vref, pref, i, ordered_supps, write_model, FDBackward)
        else
            expr_dict[curr_value] = _make_difference_expr(vref, pref, i, ordered_supps, write_model, method.technique)
        end
    end
    return expr_dict
end

function _make_difference_expr(vref::GeneralVariableRef, pref::GeneralVariableRef,
                               index::Int, 
                               ordered_supps::Vector{Float64},
                               write_model::JuMP.AbstractModel, 
                               type::Type{FDForward})::JuMP.AbstractJuMPScalar
    
    curr_value = ordered_supps[index]
    next_value = ordered_supps[index+1]
    return JuMP.@expression(InfiniteOpt._Model, (make_reduced_expr(vref, pref, next_value, write_model) 
                                            - make_reduced_expr(vref, pref, curr_value, write_model))
                                            / (next_value - curr_value) )
end

function _make_difference_expr(vref::GeneralVariableRef, pref::GeneralVariableRef,
                               index::Int, 
                               ordered_supps::Vector{Float64}, 
                               write_model::JuMP.AbstractModel, 
                               type::Type{FDCentral})::JuMP.AbstractJuMPScalar
    prev_value = ordered_supps[index-1]
    next_value = ordered_supps[index+1]
    return JuMP.@expression(InfiniteOpt._Model, (make_reduced_expr(vref, pref, next_value, write_model) 
                                            - make_reduced_expr(vref, pref, prev_value, write_model)) 
                                            / (next_value - prev_value) )
end

function _make_difference_expr(vref::GeneralVariableRef, pref::GeneralVariableRef,
                               index::Int, 
                               ordered_supps::Vector{Float64}, 
                               write_model::JuMP.AbstractModel,
                               type::Type{FDBackward})::JuMP.AbstractJuMPScalar
    prev_value = ordered_supps[index-1]
    curr_value = ordered_supps[index]
    return JuMP.@expression(InfiniteOpt._Model, (make_reduced_expr(vref, pref, curr_value, write_model) 
                                            - make_reduced_expr(vref, pref, prev_value, write_model)) 
                                            / (curr_value - prev_value) )
end

function evaluate_derivative(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                             method::OrthogonalCollocation, 
                             write_model::JuMP.AbstractModel)
# NOTE These should use `add_derivative_supports` and `make_reduced_expr`.
    
end


################################################################################
#                              EVALUATION METHODS
################################################################################
"""
    evaluate(dref::DerivativeRef)::Nothing

Numerically evaluate `dref` by computing its auxiliary derivative constraints 
(e.g., collocation equations) and add them to the model. For normal usage, it is 
recommended that this method not be called directly and instead have TranscriptionOpt 
handle these equations, since preemptive evaluation can lead to invalid relations 
if the support structure is modified. Errors if `evaluate_derivative` is not 
defined for the derivative method employed.

**Example**
```julia-repl 
TODO ADD EXAMPLE
```
"""
function evaluate(dref::DerivativeRef)::Nothing
    # collect the basic info
    vref = argument_variable(dref)
    pref = operator_parameter(dref)
    method = derivative_method(pref)
    model = JuMP.owner_model(dref)
    # get the expressions 
    exprs = evaluate_derivative(vref, pref, method, model)
    # add the constraints

    # TODO ADD DERIVATIVE AUXILIARY CONSTRAINTS HERE

    # update the parameter status 
    _data_object(pref).has_deriv_constrs = true
    return
end

"""
    evaluate_all_derivatives!(model::InfiniteModel)::Nothing

Evaluate all the derivatives in `model` by adding the corresponding auxiliary 
equations to `model`. See [`evaluate`](@ref) for more information.

**Example**
```julia-repl 
TODO ADD EXAMPLE
```
"""
function evaluate_all_derivatives!(model::InfiniteModel)::Nothing
    for dref in all_derivatives(model)
        evaluate(dref)
    end
    return
end