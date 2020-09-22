################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::InfiniteDerivativeIndex
                               )::InfiniteDerivativeRef
    return InfiniteDerivativeRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::DerivativeData{<:InfiniteDerivative}
                          )::InfiniteDerivativeIndex
    return MOIUC.add_item(model.infinite_derivs, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{InfiniteDerivative}
    )::MOIUC.CleverDict{InfiniteDerivativeIndex, DerivativeData{<:InfiniteDerivative{<:GeneralDerivativeRef}}}
    return model.infinite_derivs
end

# Extend _data_dictionary (reference based)
function _data_dictionary(dref::InfiniteDerivativeRef
    )::MOIUC.CleverDict{InfiniteDerivativeIndex, DerivativeData{<:InfiniteDerivative{<:GeneralDerivativeRef}}}
    return JuMP.owner_model(dref).infinite_derivs
end

# Extend _data_object
function _data_object(dref::InfiniteDerivativeRef
    )::DerivativeData{<:InfiniteDerivative{<:GeneralDerivativeRef}}
    object = get(_data_dictionary(dref), JuMP.index(dref), nothing)
    object === nothing && error("Invalid infinite derivative reference, cannot find " *
    "corresponding derivative in the model. This is likely " *
    "caused by using the reference of a deleted derivative.")
    return object
end

# Extend _core_variable_object
function _core_variable_object(dref::InfiniteDerivativeRef
    )::InfiniteDerivative{<:GeneralDerivativeRef}
    return _data_object(dref).derivative
end

# Define getter function for deriv.is_vector_start
function _is_vector_start(dref::InfiniteDerivativeRef)::Bool
    return _core_derivative_object(dref).is_vector_start
end

# Define getter function for deriv.eval_method
function _eval_method(dref::InfiniteDerivativeRef)::AbstractDerivativeEvalData
    return _core_derivative_object(dref).eval_method
end

################################################################################
#                          DEFINTION METHODS
################################################################################
const DefaultEvalMethod = Integral(10, MeasureToolbox.Automatic)

# Define build_derivative
function build_derivative(_error::Function, info::JuMP.VariableInfo, 
                          numerator::GeneralVariableRef, 
                          denominator::GeneralVariableRef; 
                          eval_method::AbstractDerivativeMethod = DefaultEvalMethod
                          )::InfiniteDerivative
    # check the derivative numerator and denominator 
    if !(_index_type(denominator) <: InfiniteParameterIndex)
        _error("Derivatives must be with respect to an infinite parameter, but " * 
               "$denominator was given.")
    elseif !(_parameter_number(denominator) in _parameter_numbers(numerator))
        _error("Derivative numerators must dependent on the infinite parameter " * 
               "it is with respect to. Here $numerator does not depend on " * 
               "$denominator.")
    end
    # check and format the info correctly
    prefs = raw_parameter_refs(numerator)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # make the derivative and return it 
    return InfiniteDerivative(new_info, is_vect_func, numerator, denominator, 
                              eval_method)
end

# Define add_derivative 
function add_derivative(model::InfiniteModel, d::InfiniteDerivative, 
                        name::String = "")::GeneralVariableRef
    # check the validity 
    vref = dispatch_variable_ref(d.variable_ref)
    pref = dispatch_variable_ref(d.parameter_ref)
    JuMP.check_belongs_to_model(vref, model)
    JuMP.check_belongs_to_model(pref, model)
    # add it to the model and make the reference
    data_object = DerivativeData(d, name)
    dindex = _add_data_object(model, data_object)
    dref = InfiniteDerivativeRef(model, dindex)
    # update the mappings 
    push!(_derivative_dependencies(vref), dindex)
    push!(_derivative_dependencies(pref), dindex)
    # add the info constraints
    gvref = _make_variable_ref(model, dindex)
    _set_info_constraints(d.info, gvref, dref)
    return gvref
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################

