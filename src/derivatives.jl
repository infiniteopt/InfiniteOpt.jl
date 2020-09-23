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
    return _data_object(dref).variable
end

# Define getter function for deriv.is_vector_start
function _is_vector_start(dref::DerivativeRef)::Bool
    return _core_variable_object(dref).is_vector_start
end

"""
    derivative_argument(dref::DerivativeRef)::GeneralVariableRef

Returns the infinite variable/derivative reference that is the input the
differential operator (i.e., the dependent variable of the derivative).

**Example**
```julia-repl
julia> derivative_argument(dref) 
x(t)
```
"""
function derivative_argument(dref::DerivativeRef)::GeneralVariableRef
    return _core_variable_object(dref).variable_ref
end

"""
    operator_parameter(dref::DerivativeRef)::GeneralVariableRef

Returns the infinite parameter reference that is what the differential operator 
is operating with respect to (i.e., the independent  variable of the derivative).

**Example**
```julia-repl
julia> operator_parameter(dref) 
t
```
"""
function operator_parameter(dref::DerivativeRef)::GeneralVariableRef
    return _core_variable_object(dref).parameter_ref
end

"""
    eval_method(dref::DerivativeRef)::AbstractDerivativeMethod

Returns the evaluation method employed by `dref` that determines the numerical 
computation scheme that will be used to evaluate the derivative.

**Example**
```julia-repl
julia> eval_method(dref) 
Integral(10, Automatic)
```
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
    build_derivative(_error::Function, info::JuMP.VariableInfo, 
                     argument_ref::GeneralVariableRef, 
                     parameter_ref::GeneralVariableRef; 
                     [eval_method::AbstractDerivativeMethod = DefaultEvalMethod]
                     )::Derivative

Constructs and returns a [`Derivative`](@ref) with a differential operator that 
depends on `parameter_ref` and operates on `argument_ref`. Variable `info` can also 
be provided to associate this derivative with bounds and a starting value function 
like that of infinite variables. Here `eval_method` refers to the 
`AbstractDerivativeMethod` use to numerically evaulate the derivative. Errors 
when `argument_ref` is not an infinite/reduced variable or derivative that depends 
on `parameter_ref`.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1]); @infinite_variable(m, x(t));

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> build_derivative(error, info, x, t)
Derivative{GeneralVariableRef, Integral}(VariableInfo{Float64,Float64,Float64,Function}(false, 0.0, false, 0.0, false, 0.0, false, start_func, false, false), true, x(t), t, Integral(10, Automatic))
````
"""
function build_derivative(_error::Function, info::JuMP.VariableInfo, 
                          argument_ref::GeneralVariableRef, 
                          parameter_ref::GeneralVariableRef; 
                          eval_method::AbstractDerivativeMethod = DefaultEvalMethod
                          )::Derivative
    # check the derivative numerator and denominator 
    if !(_index_type(parameter_ref) <: InfiniteParameterIndex)
        _error("Derivatives must be with respect to an infinite parameter, but " * 
               "$parameter_ref was given.")
    elseif !(_parameter_number(parameter_ref) in _parameter_numbers(argument_ref))
        _error("Derivative arguments must dependent on the infinite parameter " * 
               "it is with respect to. Here $argument_ref does not depend on " * 
               "$parameter_ref.")
    elseif argument_ref == parameter_ref
        _error("Cannot specify the operator parameter as the derivative argument.")
    end
    # check and format the info correctly
    prefs = raw_parameter_refs(argument_ref)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # make the derivative and return it 
    return Derivative(new_info, is_vect_func, argument_ref, parameter_ref, eval_method)
end

"""
    add_derivative(model::InfiniteModel, d::Derivative, 
                   [name::String = ""])::GeneralVariableRef

Adds a derivative `d` to `model` and returns a `GeneralVariableRef` that points 
to it. Errors if the derivative dependencies do not belong to `model`. 

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1]); @infinite_variable(m, x(t));

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> d = build_derivative(error, info, x, t);

julia> dref = add_derivative(m, d)
∂/∂t[x(t)]
```
"""
function add_derivative(model::InfiniteModel, d::Derivative, 
                        name::String = "")::GeneralVariableRef
    # check the validity 
    vref = dispatch_variable_ref(d.variable_ref)
    pref = dispatch_variable_ref(d.parameter_ref)
    JuMP.check_belongs_to_model(vref, model)
    JuMP.check_belongs_to_model(pref, model)
    # add it to the model and make the reference
    data_object = VariableData(d, name)
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

# Default derivative info 
const _DefaultInfo = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                        false, NaN, false, NaN, false, NaN, false, s -> NaN, 
                        false, false)

## Define function that build derivative expressions following rules of calculus
# GeneralVariableRef
function _build_deriv_expr(vref::GeneralVariableRef, pref, eval_method
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    if vref == pref 
        return 1.0
    elseif _parameter_number(pref) in _parameter_numbers(vref)
        d = Derivative(_DefaultInfo, true, vref, pref, eval_method)
        return add_derivative(JuMP.owner_model(vref), d)
    else 
        return 0.0
    end
end

# AffExpr
function _build_deriv_expr(aff::JuMP.GenericAffExpr, pref, eval_method
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    return JuMP.@expression(_Model, sum(c * _build_deriv_expr(v, pref, eval_method) 
                            for (c, v) in JuMP.linear_terms(aff)))
end

# Quad Expr (implements product rule)
function _build_deriv_expr(quad::JuMP.GenericQuadExpr, pref, eval_method
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    return JuMP.@expression(_Model, sum(c * (_build_deriv_expr(v1, pref, eval_method) * v2 + 
                                        v1 * _build_deriv_expr(v2, pref, eval_method)) 
                                        for (c, v1, v2) in JuMP.quad_terms(quad)) + 
                                        _build_deriv_expr(quad.aff, pref, eval_method))
end

# Real number
function _build_deriv_expr(expr::Real, pref, eval_method)::Float64
    return 0.0
end

# Fallback 
function _build_deriv_expr(expr, pref, eval_method)
    error("Unsupported expression type $(typeof(expr)) for derivative definition.")
end

# Define recursive function to call _build_deriv_expr 
function _recursive_deriv_build(expr, prefs, eval_method
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    if isempty(prefs)
        return expr 
    else
        new_expr = _build_deriv_expr(expr, first(prefs), eval_method)
        return _recursive_deriv_build(new_expr, prefs[2:end], eval_method)
    end
end

"""
    deriv(expr::JuMP.AbstractJuMPScalar, pref1::GeneralVariableRef[, ....]; 
          [eval_method::AbstractDerivativeMethod = DefaultEvalMethod]
          )::Union{JuMP.AbstractJuMPScalar, Float64}

Apply appropriate calculus methods to define and return the derivative expression of `expr` 
with respect to the infinite parameter(s) `pref1`, pref2`, etc. in that respective 
order. This will implicilty build and add individual [`Derivative`](@ref)s as 
appropriate. Errors if no infinite parameter is given or if the parameters are 
not infinite. Here `eval_method` refers to the concrete `AbstractDerivativeMethod` 
that will be used to evaluate the individual derivatives numerically.

**Example**
```julia-repl 
deriv(x^2 + z, t, t)
# TODO ADD RESULT ONCE PRINTING IS FINALIZED
```
"""
function deriv(expr, prefs::GeneralVariableRef...; 
               eval_method::AbstractDerivativeMethod = DefaultEvalMethod
               )::Union{JuMP.AbstractJuMPScalar, Float64}
    # Check inputs 
    if !all(_index_type(pref) <: InfiniteParameterIndex for pref in prefs)
        error("Can only take derivative with respect to infinite parameters.")
    elseif length(prefs) == 0
        error("Must specify at least one infinite parameter with which to define " * 
              "the derivative operator.")
    end
    # Build the derivative expression
    return _recursive_deriv_build(expr, prefs, eval_method)
end

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
primarily an internal method where [`parameter_refs`](@ref parameter_refs(::DerivativeRef))
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
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    method = eval_method(dref)
    new_deriv = Derivative(new_info, true, vref, pref, method)
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
