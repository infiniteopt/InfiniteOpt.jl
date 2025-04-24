################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel, index::DerivativeIndex)
    return DerivativeRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:Derivative}
    )
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
function _data_object(dref::DerivativeRef)
    object = get(_data_dictionary(dref), JuMP.index(dref), nothing)
    if isnothing(object) 
        error("Invalid derivative reference, cannot find ",
              "corresponding derivative in the model. This is likely ",
              "caused by using the reference of a deleted derivative.")
    end
    return object
end

"""
    core_object(dref::DerivativeRef)::Derivative

Retrieve the underlying core [`Derivative`](@ref) object for `dref`. 
This is intended as an advanced method for developers.
"""
function core_object(dref::DerivativeRef)
    return _data_object(dref).variable
end

# Define getter function for deriv.is_vector_start
function _is_vector_start(dref::DerivativeRef)
    return core_object(dref).is_vector_start
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
function derivative_argument(dref::DerivativeRef)
    return core_object(dref).variable_ref
end

"""
    derivative_order(dref::DerivativeRef)::Int

Return the what order the derivative is.

**Example**
```julia-repl
julia> derivative_order(dref) 
1
```
"""
function derivative_order(dref::DerivativeRef)
    return core_object(dref).order
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
function operator_parameter(dref::DerivativeRef)
    return core_object(dref).parameter_ref
end

"""
    derivative_method(dref::DerivativeRef)::AbstractDerivativeMethod

Returns the evaluation method employed by `dref` that determines the numerical 
computation scheme that will be used to evaluate the derivative. Note that this 
is set on by the infinite parameter with respect to which the derivative is 
defined.

**Example**
```julia-repl
julia> derivative_method(dref) 
FiniteDifference(Backward, true)
```
"""
function derivative_method(dref::DerivativeRef)
    return derivative_method(operator_parameter(dref))
end

"""
    parameter_group_int_indices(dref::DerivativeRef)::Vector{Int}

Return the list of infinite parameter group integer indices used by `dref`.
"""
function parameter_group_int_indices(dref::DerivativeRef)
    return parameter_group_int_indices(derivative_argument(dref))
end

# Extend _parameter_numbers
function _parameter_numbers(dref::DerivativeRef)
    return _parameter_numbers(derivative_argument(dref))
end

# Get derivative index for variable-parameter pair if one exists 
function _existing_derivative_index(
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef,
    order::Int
    )
    model = JuMP.owner_model(vref)
    return get(model.deriv_lookup, (vref, pref, order), nothing)
end

# Get access the derivative constraint indices 
function _derivative_constraint_dependencies(dref::DerivativeRef)
    return _data_object(dref).deriv_constr_indices
end

"""
    has_derivative_constraints(dref::DerivativeRef)::Bool

Return a `Bool` whether `dref` has been evaluated within the `InfiniteModel` and 
has derivative constraints that have been added to the `InfiniteModel`. Note this 
does not indicate if such constraints have been added to the transformation backend. Thus, 
with normal usage (i.e., not using `evaluate`) this should always return `false`.
"""
function has_derivative_constraints(dref::DerivativeRef)
    return !isempty(_derivative_constraint_dependencies(dref))
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(
    ::DerivativeRef, 
    d::D, 
    data::VariableData{D}
    ) where {D <: Derivative}
    data.variable = d
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(
    dref::DerivativeRef, 
    d::D1, 
    data::VariableData{D2}
    ) where {D1, D2}
    new_data = VariableData(d, data.name, data.lower_bound_index,
                            data.upper_bound_index, data.fix_index,
                            data.zero_one_index, data.integrality_index,
                            data.measure_indices, data.constraint_indices,
                            data.in_objective, data.point_var_indices,
                            data.semi_infinite_var_indices, 
                            data.derivative_indices, data.deriv_constr_indices)
    _data_dictionary(dref)[JuMP.index(dref)] = new_data
    return
end

# Extend _set_core_object for DerivativeRefs
function _set_core_object(dref::DerivativeRef, d::Derivative)
    _adaptive_data_update(dref, d, _data_object(dref))
    return
end

################################################################################
#                          DEFINTION METHODS
################################################################################
# If possible and applicable, unnest a derivative of a derivative
# For instance if we take the 1st derivative of a 2nd derivative with respect to 
# the same infinite parameter, this will ensure we produce a single 3rd order derivative
function _unnest_derivative(vref, pref, order)
    if vref.index_type == DerivativeIndex && operator_parameter(vref) == pref
        order += derivative_order(vref)
        vref = derivative_argument(vref)
    end
    return vref, order
end

"""
    build_derivative(_error::Function, info::JuMP.VariableInfo, 
                     argument_ref::GeneralVariableRef, 
                     parameter_ref::GeneralVariableRef,
                     order::Int = 1
                     )::Derivative

Constructs and returns a [`Derivative`](@ref) with a differential operator that 
depends on `parameter_ref` and operates on `argument_ref`. The order of the derivative 
is dictated by `order`. Variable `info` can also 
be provided to associate this derivative with bounds and a starting value function 
like that of infinite variables. Errors when `argument_ref` is not an 
infinite/semi-infinite variable or derivative that depends on `parameter_ref`.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1]); @variable(m, y, Infinite(t));

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> build_derivative(error, info, y, t)
Derivative{GeneralVariableRef}(VariableInfo{Float64,Float64,Float64,Function}(false, 0.0, false, 0.0, false, 0.0, false, start_func, false, false), true, y(t), t, 1)
````
"""
function build_derivative(
    _error::Function, 
    info::JuMP.VariableInfo, 
    argument_ref::GeneralVariableRef, 
    parameter_ref::GeneralVariableRef,
    order::Int = 1
    )
    # check the derivative numerator and denominator 
    if !(_index_type(parameter_ref) <: InfiniteParameterIndex)
        _error("Derivatives must be with respect to an infinite parameter, but " * 
               "$parameter_ref was given.")
    elseif !(_parameter_number(parameter_ref) in _parameter_numbers(argument_ref))
        _error("Derivative arguments must dependent on the infinite parameter " * 
               "it is with respect to. Here $argument_ref does not depend on " * 
               "$parameter_ref.")
    elseif isequal(argument_ref, parameter_ref)
        _error("Cannot specify the operator parameter as the derivative argument.")
    elseif order <= 0 
        _error("The specified derivative order $order is not a positive integer.")
    end
    # check and format the info correctly
    prefs = raw_parameter_refs(argument_ref)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # make the derivative and return it 
    vref, order = _unnest_derivative(argument_ref, parameter_ref, order)
    return Derivative(new_info, is_vect_func, vref, parameter_ref, order)
end

"""
    Deriv{V, P} <: InfOptVariableType

A `DataType` to assist in making derivative variables. This can be passed as an 
extra argument to `@variable` to make such a variable: 
```julia 
@variable(model, var_expr, Deriv(inf_var, inf_par [, order]), kwargs...)
```
Here, `inf_var` is the infinite variable that is being operated on and `inf_par` 
is the infinite parameter that the derivative is defined with respect to. `order` 
is an optional argument and defaults to 1.

**Fields**
- `argument::V`: The infinite variable being operated on.
- `operator_parameter::P`: The infinite parameter that determines the derivative.
- `order::Int`: The derivative order (defaults to 1 if not given).
"""
struct Deriv{V, P} <: InfOptVariableType 
    argument::V 
    operator_parameter::P
    order::Int
    function Deriv(arg::V, pref::P, order::Int = 1) where {V, P}
        return new{V, P}(arg, pref, order)
    end
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, 
                        var_type::Deriv)::InfiniteVariable{GeneralVariableRef}

Build and return a first order derivative based on `info` and `var_type`. Errors 
if the information in `var_type` is invalid. See [`Deriv`](@ref) for more 
information.

**Example**
```julia-repl
julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> deriv_var = build_variable(error, info, Deriv(y, t));
```
"""
function JuMP.build_variable(
    _error::Function, 
    info::JuMP.VariableInfo, 
    var_type::Deriv;
    extra_kw_args...
    )
    # Error extra keyword arguments
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with derivatives.")
    end
    # check for valid inputs
    ivref = var_type.argument
    if !(ivref isa GeneralVariableRef)
        _error("Expected an infinite variable reference dependency ",
               "of type `GeneralVariableRef`, but got an argument of type ",
               "$(typeof(ivref)).")
    end
    pref = var_type.operator_parameter
    if !(pref isa GeneralVariableRef)
        _error("Expected an infinite parameter reference dependency ",
               "of type `GeneralVariableRef`, but got an argument of type ",
               "$(typeof(pref)).")
    end
    return build_derivative(_error, info, ivref, pref, var_type.order)
end

"""
    add_derivative(model::InfiniteModel, d::Derivative, 
                   [name::String = ""])::GeneralVariableRef

Adds a derivative `d` to `model` and returns a `GeneralVariableRef` that points 
to it. Errors if the derivative dependencies do not belong to `model`. Note that 
`d` should be built using [`build_derivative`](@ref) to avoid nuance internal 
errors.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1]); @variable(m, y, Infinite(t));

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> d = build_derivative(error, info, y, t, 2);

julia> dref = add_derivative(m, d)
∂²/∂t²[y(t)]
```
"""
function add_derivative(model::InfiniteModel, d::Derivative, name::String = "")
    # check the validity 
    vref = dispatch_variable_ref(d.variable_ref)
    pref = dispatch_variable_ref(d.parameter_ref)
    JuMP.check_belongs_to_model(vref, model)
    JuMP.check_belongs_to_model(pref, model)
    # check if we already have a derivative 
    existing_index = _existing_derivative_index(d.variable_ref, d.parameter_ref, d.order)
    if isnothing(existing_index)
        # add it to the model and make the reference
        data_object = VariableData(d, name)
        dindex = _add_data_object(model, data_object)
        dref = DerivativeRef(model, dindex)
        # update the mappings 
        push!(_derivative_dependencies(vref), dindex)
        push!(_derivative_dependencies(pref), dindex)
        # update the derivative lookup dict
        model.deriv_lookup[(d.variable_ref, d.parameter_ref, d.order)] = dindex
        # add the info constraints
        gvref = GeneralVariableRef(model, dindex)
        _set_info_constraints(d.info, gvref, dref)
    else
        dref = DerivativeRef(model, existing_index)
        gvref = GeneralVariableRef(model, existing_index)
        old_info = _variable_info(dref)
        if old_info.has_lb || old_info.has_ub || old_info.has_fix || old_info.has_start
            @warn "Overwriting $dref, any previous properties (e.g., lower bound " * 
                  "or start value) will be lost/changed."
        end
        _update_info_constraints(d.info, gvref, dref)
        _set_core_object(dref, d)
        if !isempty(name)
            set_name(dref, name)
        end
    end
    return gvref
end

# Extend JuMP.add_variable to enable macro definition 
function JuMP.add_variable(model::InfiniteModel, d::Derivative, name::String = "")
    return add_derivative(model, d, name)
end

# Helper method to to quickly build and add derivatives internally
function _build_add_derivative(vref, pref, order)
    vref, order = _unnest_derivative(vref, pref, order)
    dindex = _existing_derivative_index(vref, pref, order)
    model = JuMP.owner_model(vref)
    if isnothing(dindex)
        info = VariableInfo(false, NaN, false, NaN, false, NaN, false, 
                            s -> NaN, false, false)
        d = Derivative(info, true, vref, pref, order)
        return add_derivative(model, d)
    else 
        return GeneralVariableRef(model, dindex)
    end
end

## Define function that build derivative expressions following rules of calculus
# GeneralVariableRef
function _build_deriv_expr(vref::GeneralVariableRef, pref, order)
    if iszero(order) # needed for the product rule
        return vref
    elseif isequal(vref, pref)
        return isone(order) ? 1.0 : 0.0
    elseif _parameter_number(pref) in _parameter_numbers(vref)
        return _build_add_derivative(vref, pref, order)
    else 
        return 0.0
    end
end

# AffExpr
function _build_deriv_expr(aff::JuMP.GenericAffExpr, pref, order)
    return @_expr(sum(c * _build_deriv_expr(v, pref, order) 
                        for (c, v) in JuMP.linear_terms(aff); init = 0.0))
end

# Quad Expr (implements product rule)
function _build_deriv_expr(quad::JuMP.GenericQuadExpr, pref, order)
    return @_expr(sum(c * sum(binomial(order, k) * _build_deriv_expr(v1, pref, order-k) * _build_deriv_expr(v2, pref, k) for k = 0:order) 
                        for (c, v1, v2) in JuMP.quad_terms(quad)) + 
                        _build_deriv_expr(quad.aff, pref, order))
end

# Real number
function _build_deriv_expr(expr::Real, pref, order)
    return 0.0
end

# Fallback 
function _build_deriv_expr(expr, pref, order)
    error("Unsupported expression type $(typeof(expr)) for derivative definition.")
end

# Define recursive function to call _build_deriv_expr 
function _recursive_deriv_build(expr, prefs, orders)
    if isempty(prefs)
        return expr 
    else
        new_expr = _build_deriv_expr(expr, first(prefs), first(orders))
        return _recursive_deriv_build(new_expr, prefs[2:end], orders[2:end])
    end
end

# Given a tuple parameter references, return tuple without repeated values and another 
# tuple that reports the number of repeats of each in a row (the order)
function _detect_repeated(pref_tup)
    length(pref_tup) == 1 && return pref_tup, (1,)
    counts = Int[]
    prefs = [first(pref_tup)]
    count = 0
    for pref in pref_tup
        if pref == prefs[end]
            count += 1
        else
            push!(counts, count)
            push!(prefs, pref)
            count = 1
        end
    end
    push!(counts, count)
    return Tuple(prefs), Tuple(counts)
end

"""
    deriv(expr::JuMP.AbstractJuMPScalar, pref1::GeneralVariableRef[, ....]
          )::Union{JuMP.AbstractJuMPScalar, Float64}

Apply appropriate calculus methods to define and return the derivative expression of `expr` 
with respect to the infinite parameter(s) `pref1`, pref2`, etc. in that respective 
order. This will implicilty build and add individual [`Derivative`](@ref)s as 
appropriate. Errors if no infinite parameter is given or if the parameters are 
not infinite.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1])
t

julia> @variable(m, y, Infinite(t))
y(t)

julia> @variable(m, z)
z

julia> deriv_expr = deriv(y^2 + z, t, t)
2 ∂²/∂t²[y(t)]*y(t) + 2 ∂/∂t[y(t)]²
```
"""
function deriv(expr, prefs::GeneralVariableRef...)
    # Check inputs 
    if !all(_index_type(pref) <: InfiniteParameterIndex for pref in prefs)
        error("Can only take derivative with respect to infinite parameters.")
    elseif length(prefs) == 0
        error("Must specify at least one infinite parameter with which to define " * 
              "the derivative operator.")
    end
    # Build the derivative expression
    processed_prefs, pref_orders = _detect_repeated(prefs)
    return _recursive_deriv_build(expr, processed_prefs, pref_orders)
end

"""
    @deriv(expr, pref_expr1[, ...])::Union{JuMP.AbstractJuMPScalar, Float64}

The macro variant of [`deriv`](@ref) that is more efficient for expression building 
and enables symbolic differential operator parameter defintions via `pref_expr`s. 
Like `deriv` expr can be any InfiniteOpt expression and the appropriate calculus 
rules will applied to `expr` to take its derivative with respect to the indicated 
infinite parameters detailed by the `pref_expr`s. The resulting derivative 
expression will contain individual derivatives that were created and added to the 
InfiniteModel as needed. Here each `pref_expr` arugment can be of the form:
- `pref::GeneralVariableRef`: An indiviudal infinite parameter reference
- `(pref::GeneralVariableRef)^(p::Int)`: An infinite parameter applied `p` times.
Thus, the syntax `@deriv(expr, pref^2)` is equivalent to `@deriv(expr, pref, pref)`. 

This will error if `pref_expr` is an unrecongnized syntax, no infinite parameter 
is given, or if any of the specified parameters are not infinite.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1])
t

julia> @variable(m, y, Infinite(t))
y(t)

julia> @variable(m, z)
z

julia> deriv_expr = @deriv(y^2 + z, t^2)
2 ∂²/∂t²[y(t)]*y(t) + 2 ∂/∂t[y(t)]²
```
"""
macro deriv(expr, args...)
    # make an error function
    error_fn = JuMPC.build_error_fn(:deriv, (expr, args...), __source__)

    # process the inputs
    extra, _ = JuMPC.parse_macro_arguments(
        error_fn, 
        args, 
        valid_kwargs = Symbol[]
    )

    # expand the parameter references as needed with powers
    pref_exprs = []
    for p in extra
        if isexpr(p, :call) && p.args[1] == :^
            append!(pref_exprs, repeat([esc(p.args[2])], p.args[3]))
        else
            push!(pref_exprs, esc(p))
        end
    end
    # prepare the code to call deriv
    expression = _MA.rewrite_and_return(expr)
    return :( deriv($expression, $(pref_exprs...)) ) # TODO call _recursive_deriv_build directly
end

"""
    ∂(expr::JuMP.AbstractJuMPScalar, pref1::GeneralVariableRef[, ....]
      )::Union{JuMP.AbstractJuMPScalar, Float64}

This serves as a convenient unicode wrapper for [`deriv`](@ref). The `∂` is 
produced via `\\partial`.
"""
function ∂(
    expr, 
    prefs::GeneralVariableRef...
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    return deriv(expr, prefs...)
end

"""
    @∂(expr, pref_expr1[, ...])::Union{JuMP.AbstractJuMPScalar, Float64}

This serves as a convenient unicode wrapper for [`@deriv`](@ref). The `∂` is 
produced via `\\partial`.
"""
macro ∂(expr, args...)
    return esc(:( @deriv($expr, $(args...)) ))
end

################################################################################
#                           PARAMETER REFERENCES
################################################################################
"""
    raw_parameter_refs(dref::DerivativeRef)::VectorTuple

Return the raw [`VectorTuple`](@ref InfiniteOpt.Collections.VectorTuple) of the 
parameter references that `dref` depends on. This is primarily an internal method 
where [`parameter_refs`](@ref parameter_refs(::DerivativeRef)) 
is intended as the preferred user function.
"""
function raw_parameter_refs(dref::DerivativeRef)
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
function parameter_refs(dref::DerivativeRef)
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
function _update_variable_info(
    dref::DerivativeRef,
    info::JuMP.VariableInfo
    )::Nothing
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    is_vect_func = _is_vector_start(dref)
    new_deriv = Derivative(_format_infinite_info(info), is_vect_func, vref, pref, order)
    _set_core_object(dref, new_deriv)
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
function set_start_value_function(
    dref::DerivativeRef,
    start::Union{Real, Function}
    )::Nothing
    info = _variable_info(dref)
    set_transformation_backend_ready(JuMP.owner_model(dref), false)
    prefs = raw_parameter_refs(dref)
    temp_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 true, start, info.binary, info.integer)
    new_info, is_vect_func = _check_and_format_infinite_info(error, temp_info, prefs)
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    new_deriv = Derivative(new_info, is_vect_func, vref, pref, order)
    _set_core_object(dref, new_deriv)
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
    set_transformation_backend_ready(JuMP.owner_model(dref), false)
    new_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 false, s -> NaN, info.binary, info.integer)
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    new_deriv = Derivative(new_info, true, vref, pref, order)
    _set_core_object(dref, new_deriv)
    return
end

# TODO maybe throw errors for binary and integer?

################################################################################
#                              MODEL QUERIES
################################################################################
"""
    num_derivatives(model::InfiniteModel)::Int

Returns the number of derivatives that have been defined in `model`. Note that 
nested derivatives will be counted in accordance with their components (e.g., 
``\\frac{d^2 x(t)}{dt^2} = ``\\frac{d}{dt}\\left(\\frac{d x(t)}{dt} \\right)`` 
will count as 2 derivatives.)

**Example**
```julia-repl
julia> num_derivatives(model)
12
```
"""
function num_derivatives(model::InfiniteModel)::Int 
    return length(_data_dictionary(model, Derivative))
end

"""
    all_derivatives(model::InfiniteModel)::Vector{GeneralVariableRef}

Returns a list of all the individual derivatives stored in `model`. 

**Example**
```julia-repl
julia> all_derivatives(model)
3-element Array{GeneralVariableRef,1}:
 ∂/∂t[T(x, t)]
 ∂/∂x[T(x, t)]
 ∂/∂x[∂/∂x[T(x, t)]]
```
"""
function all_derivatives(model::InfiniteModel)::Vector{GeneralVariableRef}
    vrefs_list = Vector{GeneralVariableRef}(undef, num_derivatives(model))
    for (i, (index, _)) in enumerate(_data_dictionary(model, Derivative))
        vrefs_list[i] = GeneralVariableRef(model, index)
    end
    return vrefs_list
end

################################################################################
#                      EVALUATION METHOD MODIFICATIONS/QUERIES
################################################################################
# Fallback (needs to be done via the infinite parameter)
function set_derivative_method(::DerivativeRef, ::AbstractDerivativeMethod)
    error("Cannot specify the derivative method in terms of the derivative, must " * 
          " do so in terms of the operator parameter by calling `set_derivative_method` " * 
          " using the operator parameter.")
end

"""
    derivative_constraints(dref::DerivativeRef)::Vector{InfOptConstraintRef}

Return a list of the derivative evaluation constraints for `dref` that have been 
added directly to the `InfiniteModel` associated with `dref`. An empty vector is 
returned is there are no such constraints.
"""
function derivative_constraints(dref::DerivativeRef)::Vector{InfOptConstraintRef}
    return [InfOptConstraintRef(JuMP.owner_model(dref), idx) 
            for idx in _derivative_constraint_dependencies(dref)]
end

"""
    delete_derivative_constraints(dref::DerivativeRef)::Nothing

Delete any derivative constraints of `dref` that have been directly added to the 
`InfiniteModel`.
"""
function delete_derivative_constraints(dref::DerivativeRef)::Nothing 
    model = JuMP.owner_model(dref)
    for cref in derivative_constraints(dref)
        JuMP.delete(model, cref)
    end
    empty!(_derivative_constraint_dependencies(dref))
    pref = dispatch_variable_ref(operator_parameter(dref))
    if all(!has_derivative_constraints(DerivativeRef(model, idx)) 
            for idx in _derivative_dependencies(pref))
        _set_has_derivative_constraints(pref, false)
    end
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(dref::DerivativeRef)::Nothing
    # remove variable info constraints associated with dref
    _delete_info_constraints(dref)
    # update variable and parameter mapping
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    filter!(e -> e != JuMP.index(dref), _derivative_dependencies(vref))
    filter!(e -> e != JuMP.index(dref), _derivative_dependencies(pref))
    model = JuMP.owner_model(dref)
    # delete the variable-parameter pair from the anti-duplication dictionary 
    delete!(model.deriv_lookup, (vref, pref, order))
    # delete any derivative constraints associated with this derivative 
    delete_derivative_constraints(dref)
    # delete associated point variables and mapping
    for index in copy(_point_variable_dependencies(dref))
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated semi-infinite variables and mapping
    for index in copy(_semi_infinite_variable_dependencies(dref))
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated derivative variables and mapping 
    for index in copy(_derivative_dependencies(dref))
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    return
end
