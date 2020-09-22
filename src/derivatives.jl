################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel, index::DerivativeIndex)::DerivativeRef
    return DerivativeRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                            object::VariableData{<:Derivative} )::DerivativeIndex
    return MOIUC.add_item(model.derivatives, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{Derivative})
    return model.derivatives
end

# Extend _data_dictionary (reference based)
function _data_dictionary(dref::DerivativeRef)
    return JuMP.owner_model(dref).derivatives
end

# Extend _data_object
function _data_object(dref::DerivativeRef
    )::VariableData{<:Derivative{<:GeneralVariableRef}}
    object = get(_data_dictionary(dref), JuMP.index(dref), nothing)
    object === nothing && error("Invalid derivative reference, cannot find " *
    "corresponding derivative in the model. This is likely " *
    "caused by using the reference of a deleted derivative.")
    return object
end

# Extend _core_variable_object
function _core_variable_object(dref::DerivativeRef
    )::Derivative{<:GeneralVariableRef}
    return _data_object(dref).derivative
end

# Define getter function for deriv.is_vector_start
function _is_vector_start(dref::DerivativeRef)::Bool
    return _core_variable_object(dref).is_vector_start
end

"""

"""
function derivative_argument(dref::DerivativeRef)::GeneralVariableRef
    return _core_variable_object(dref).variable_ref
end

"""

"""
function operator_parameter(dref::DerivativeRef)::GeneralVariableRef
    return _core_variable_object(dref).parameter_ref
end

"""

"""
function eval_method(dref::DerivativeRef)::AbstractDerivativeMethod
    return _core_variable_object(dref).eval_method
end

# Extend _object_numbers
function _object_numbers(dref::DerivativeRef)::Vector{Int}
    return _object_numbers(derivative_argument(dref))
end

# Extend _parameter_numbers
function _parameter_numbers(dref::DerivativeRef)::Vector{Int}
    return _parameter_numbers(derivative_argument(dref))
end

################################################################################
#                          DEFINTION METHODS
################################################################################
const DefaultEvalMethod = Integral(10, MeasureToolbox.Automatic)

"""

"""
function build_derivative(_error::Function, info::JuMP.VariableInfo, 
                          numerator::GeneralVariableRef, 
                          denominator::GeneralVariableRef; 
                          eval_method::AbstractDerivativeMethod = DefaultEvalMethod
                          )::Derivative
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
    return Derivative(new_info, is_vect_func, numerator, denominator, 
                              eval_method)
end

"""

"""
function add_derivative(model::InfiniteModel, d::Derivative, 
                        name::String = "")::GeneralVariableRef
    # check the validity 
    vref = dispatch_variable_ref(d.variable_ref)
    pref = dispatch_variable_ref(d.parameter_ref)
    JuMP.check_belongs_to_model(vref, model)
    JuMP.check_belongs_to_model(pref, model)
    # add it to the model and make the reference
    data_object = DerivativeData(d, name)
    dindex = _add_data_object(model, data_object)
    dref = DerivativeRef(model, dindex)
    # update the mappings 
    push!(_derivative_dependencies(vref), dindex)
    push!(_derivative_dependencies(pref), dindex)
    # add the info constraints
    gvref = _make_variable_ref(model, dindex)
    _set_info_constraints(d.info, gvref, dref)
    return gvref
end

# TODO make derivative function for use in expressions

################################################################################
#                           PARAMETER REFERENCES
################################################################################
"""
    raw_parameter_refs(dref::DerivativeRef)::VectorTuple{GeneralVariableRef}

Return the raw [`VectorTuple`](@ref) of the parameter references that `dref`
depends on. This is primarily an internal method where
[`parameter_refs`](@ref parameter_refs(::DerivativeRef))
is intended as the preferred user function.
"""
function raw_parameter_refs(dref::DerivativeRef)::VectorTuple{GeneralVariableRef}
    return raw_parameter_refs(derivative_argument(dref))
end

"""
    parameter_refs(dref::DerivativeRef)::Tuple

Return the parameter references associated with the infinite derivative `dref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `dref`.

**Example**
```julia-repl
julia> parameter_refs(deriv)
(t,)
```
"""
function parameter_refs(dref::DerivativeRef)::Tuple
    return parameter_refs(derivative_argument(dref))
end

"""
    parameter_list(dref::DerivativeRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `dref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(dref::DerivativeRef))
is intended as the preferred user function.
"""
function parameter_list(dref::DerivativeRef)::Vector{GeneralVariableRef}
    return raw_parameter_refs(dref).values
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for infinite derivatives
function _update_variable_info(dref::DerivativeRef,
                               info::JuMP.VariableInfo)::Nothing
    new_info = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                                 info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 info.has_start, info.start, info.binary, info.integer)
    numer = derivative_argument(dref)
    denom = operator_parameter(dref)
    method = eval_method(dref)
    is_vect_func = _is_vector_start(dref)
    new_deriv = Derivative(new_info, is_vect_func, numer, denom, method)
    _set_core_variable_object(dref, new_deriv)
    return
end

"""
    set_start_value_function(dref::DerivativeRef,
                             start::Union{Real, Function})::Nothing

Set the start value function of `dref`. If `start::Real` then a function is
generated to such that the start value will be `start` for the entire infinite
domain. If `start::Function` then this function should map to a scalar start value
given a support value arguments matching the format of the parameter elements in
`parameter_refs(dref)`.

**Example**
```julia-repl
julia> set_start_value_function(dref, 1) # all start values will be 1

julia> set_start_value_function(dref, my_func) # each value will be made via my_func
```
"""
function set_start_value_function(dref::DerivativeRef,
                                  start::Union{Real, Function})::Nothing
    info = _variable_info(dref)
    set_optimizer_model_ready(JuMP.owner_model(dref), false)
    prefs = raw_parameter_refs(dref)
    temp_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 true, start, info.binary, info.integer)
    new_info, is_vect_func = _check_and_format_infinite_info(error, temp_info, prefs)
    numer = derivative_argument(dref)
    denom = operator_parameter(dref)
    method = eval_method(dref)
    new_deriv = InfiniteDreivative(new_info, is_vect_func, numer, denom, method)
    _set_core_variable_object(dref, new_deriv)
    return
end

"""
    reset_start_value_function(dref::DerivativeRef)::Nothing

Remove the existing start value function and return to the default. Generally,
this is triggered by deleting an infinite parameter that `dref` depends on.

**Example**
```julia-repl
julia> reset_start_value_function(dref)
```
"""
function reset_start_value_function(dref::DerivativeRef)::Nothing
    info = _variable_info(dref)
    set_optimizer_model_ready(JuMP.owner_model(dref), false)
    start_func = (s::Vector{<:Real}) -> NaN
    new_info = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                                 info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 false, start_func, info.binary, info.integer)
    numer = derivative_argument(dref)
    denom = operator_parameter(dref)
    method = eval_method(dref)
    new_deriv = InfiniteDreivative(new_info, true, numer, denom, method)
    _set_core_variable_object(dref, new_deriv)
    return
end

# TODO maybe throw errors for binary and integer?

################################################################################
#                              MODEL QUERIES
################################################################################
# derivative_by_name

# num_derivatives

# all_derivatives

################################################################################
#                                 DELETION
################################################################################
