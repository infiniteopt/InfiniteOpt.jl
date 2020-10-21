################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel, index::DerivativeIndex)::DerivativeRef
    return DerivativeRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::VariableData{<:Derivative})::DerivativeIndex
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
function derivative_method(dref::DerivativeRef)::AbstractDerivativeMethod
    return derivative_method(operator_parameter(dref))
end

# Extend _object_numbers
function _object_numbers(dref::DerivativeRef)::Vector{Int}
    return _object_numbers(derivative_argument(dref))
end

# Extend _parameter_numbers
function _parameter_numbers(dref::DerivativeRef)::Vector{Int}
    return _parameter_numbers(derivative_argument(dref))
end

# Get derivative index for variable-parameter pair if one exists 
function _existing_derivative_index(
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef
    )::Union{DerivativeIndex, Nothing}
    model = JuMP.owner_model(vref)
    return get(model.deriv_lookup, (vref, pref), nothing)
end

# Get access the derivative constraint indices 
function _derivative_constraint_dependencies(dref::DerivativeRef
    )::Vector{ConstraintIndex}
    return _data_object(dref).deriv_constr_indices
end

"""
    has_derivative_constraints(dref::DerivativeRef)::Bool

Return a `Bool` whether `dref` has been evaluated within the `InfiniteModel` and 
has derivative constraints that have been added to the `InfiniteModel`. Note this 
does not indicate if such constraints have been added to the optimizer model. Thus, 
with normal usage (i.e., not using `evaluate`) this should always return `false`.
"""
function has_derivative_constraints(dref::DerivativeRef)::Bool
    return !isempty(_derivative_constraint_dependencies(dref))
end

################################################################################
#                          DEFINTION METHODS
################################################################################
# Set a dictionary to store the eval_method defaults
const _DerivativeDefaults = Dict{Symbol, Any}(
    :info => JuMP.VariableInfo{Float64, Float64, Float64, Function}(           
                false, NaN, false, NaN, false, NaN, false, s -> NaN, false, false),
    :is_vect_start => true
    )

# Get the derivative defaults
function _derivative_defaults()::Dict{Symbol, Any}
    return _DerivativeDefaults
end

# Specify the derivative defaults
function _set_derivative_default(key::Symbol, value)::Nothing 
    _derivative_defaults()[key] = value
    return 
end

"""
    build_derivative(_error::Function, info::JuMP.VariableInfo, 
                     argument_ref::GeneralVariableRef, 
                     parameter_ref::GeneralVariableRef
                     )::Derivative

Constructs and returns a [`Derivative`](@ref) with a differential operator that 
depends on `parameter_ref` and operates on `argument_ref`. Variable `info` can also 
be provided to associate this derivative with bounds and a starting value function 
like that of infinite variables. Errors when `argument_ref` is not an 
infinite/reduced variable or derivative that depends on `parameter_ref`. Also, 
errors if such a derivative was already added to the model.

**Example**
```julia-repl 
julia> @infinite_parameter(m, t in [0, 1]); @infinite_variable(m, x(t));

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> build_derivative(error, info, x, t)
Derivative{GeneralVariableRef}(VariableInfo{Float64,Float64,Float64,Function}(false, 0.0, false, 0.0, false, 0.0, false, start_func, false, false), true, x(t), t)
````
"""
function build_derivative(_error::Function, info::JuMP.VariableInfo, 
                          argument_ref::GeneralVariableRef, 
                          parameter_ref::GeneralVariableRef
                          )::Derivative{GeneralVariableRef}
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
    elseif _existing_derivative_index(argument_ref, parameter_ref) !== nothing
        _error("Cannot double specify derivatives, a derivative operating on " * 
               "$(argument_ref) with $parameter_ref was already added.")
    end
    # check and format the info correctly
    prefs = raw_parameter_refs(argument_ref)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # make the derivative and return it 
    return Derivative(new_info, is_vect_func, argument_ref, parameter_ref)
end

# Extend _make_variable to enable JuMP.build_variable for macro definition 
function _make_variable(_error::Function, info::JuMP.VariableInfo, ::Type{Deriv};
    argument::Union{GeneralVariableRef, Nothing} = nothing,
    operator_parameter::Union{GeneralVariableRef, Nothing} = nothing,
    extra_kw_args...
    )::Derivative{GeneralVariableRef}
    for (kwarg, _) in extra_kw_args
        _error("Keyword argument $kwarg is not for use with derivatives.")
    end
    # check that we have been given parameter references
    if argument === nothing || operator_parameter === nothing
        _error("Must specify the derivative argument variable and the operator " *
               "infinite parameter that the derivative is with respect to.")
    end
    return build_derivative(_error, info, argument, operator_parameter)
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
    # update the derivative lookup dict (assumes that this is anew index) --> checked by build_derivative
    model.deriv_lookup[(d.variable_ref, d.parameter_ref)] = dindex
    # add the info constraints
    gvref = _make_variable_ref(model, dindex)
    _set_info_constraints(d.info, gvref, dref)
    return gvref
end

# Extend JuMP.add_variable to enable macro definition 
function JuMP.add_variable(model::InfiniteModel, d::Derivative, 
                           name::String = "")::GeneralVariableRef
    return add_derivative(model, d, name)
end

## Define function that build derivative expressions following rules of calculus
# GeneralVariableRef
function _build_deriv_expr(vref::GeneralVariableRef, pref
    )::Union{GeneralVariableRef, Float64}
    if vref == pref 
        return 1.0
    elseif _parameter_number(pref) in _parameter_numbers(vref)
        dindex = _existing_derivative_index(vref, pref)
        model = JuMP.owner_model(vref)
        if dindex === nothing
            defaults = _derivative_defaults()
            d = Derivative(defaults[:info], defaults[:is_vect_start], vref, pref)
            return add_derivative(model, d)
        else 
            return _make_variable_ref(model, dindex)
        end
    else 
        return 0.0
    end
end

# AffExpr
function _build_deriv_expr(aff::JuMP.GenericAffExpr, pref)
    return JuMP.@expression(_Model, sum(c * _build_deriv_expr(v, pref) 
                            for (c, v) in JuMP.linear_terms(aff)))
end

# Quad Expr (implements product rule)
function _build_deriv_expr(quad::JuMP.GenericQuadExpr, pref)
    return JuMP.@expression(_Model, sum(c * (_build_deriv_expr(v1, pref) * v2 + 
                                        v1 * _build_deriv_expr(v2, pref)) 
                                        for (c, v1, v2) in JuMP.quad_terms(quad)) + 
                                        _build_deriv_expr(quad.aff, pref))
end

# Real number
function _build_deriv_expr(expr::Real, pref)::Float64
    return 0.0
end

# Fallback 
function _build_deriv_expr(expr, pref)
    error("Unsupported expression type $(typeof(expr)) for derivative definition.")
end

# Define recursive function to call _build_deriv_expr 
function _recursive_deriv_build(expr, prefs
    )::Union{JuMP.AbstractJuMPScalar, Float64}
    if isempty(prefs)
        return expr 
    else
        new_expr = _build_deriv_expr(expr, first(prefs))
        return _recursive_deriv_build(new_expr, prefs[2:end])
    end
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

julia> @infinite_variable(m, x(t))
x(t)

julia> @hold_variable(m, z)
z

julia> deriv_expr = deriv(x^2 + z, t, t)
2 ∂/∂t[∂/∂t[x(t)]]*x(t) + 2 ∂/∂t[x(t)]²
```
"""
function deriv(expr, prefs::GeneralVariableRef...
               )::Union{JuMP.AbstractJuMPScalar, Float64}
    # Check inputs 
    if !all(_index_type(pref) <: InfiniteParameterIndex for pref in prefs)
        error("Can only take derivative with respect to infinite parameters.")
    elseif length(prefs) == 0
        error("Must specify at least one infinite parameter with which to define " * 
              "the derivative operator.")
    end
    # Build the derivative expression
    return _recursive_deriv_build(expr, prefs)
end

"""
    @deriv(expr, pref_expr1[, ...]
           )::Union{JuMP.AbstractJuMPScalar, Float64}

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

julia> @infinite_variable(m, x(t))
x(t)

julia> @hold_variable(m, z)
z

julia> deriv_expr = @deriv(x^2 + z, t^2)
2 ∂/∂t[∂/∂t[x(t)]]*x(t) + 2 ∂/∂t[x(t)]²
```
"""
macro deriv(expr, args...)
    # process the arugments
    extra, kw_args, requestedcontainer = _extract_kw_args(args)
    # expand the parameter references as needed wiwht powers
    pref_exprs = []
    for p in extra
        if isexpr(p, :call) && p.args[1] == :^
            append!(pref_exprs, repeat([p.args[2]], p.args[3]))
        else
            push!(pref_exprs, p)
        end
    end
    # prepare the code to call deriv
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    code = :( deriv($expression, $(pref_exprs...); ($(kw_args...))) ) # TODO throw error if has kw_args?
    return esc(code)
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
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    is_vect_func = _is_vector_start(dref)
    new_deriv = Derivative(new_info, is_vect_func, vref, pref)
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
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    new_deriv = Derivative(new_info, is_vect_func, vref, pref)
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
    defaults = _derivative_defaults()
    start_func = defaults[:info].start 
    has_start = defaults[:info].has_start 
    new_info = JuMP.VariableInfo{Float64, Float64, Float64, Function}(
                                 info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 has_start, start_func, info.binary, info.integer)
    vref = derivative_argument(dref)
    pref = operator_parameter(dref)
    new_deriv = Derivative(new_info, defaults[:is_vect_start], vref, pref)
    _set_core_variable_object(dref, new_deriv)
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
        vrefs_list[i] = _make_variable_ref(model, index)
    end
    return vrefs_list
end

################################################################################
#                      EVALUATION METHOD MODIFICATIONS/QUERIES
################################################################################
# Fallback (needs to be done via the infinite parameter)
function set_derivative_method(dref::DerivativeRef, method::AbstractDerivativeMethod)
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
    filter!(e -> e != JuMP.index(dref), _derivative_dependencies(vref))
    filter!(e -> e != JuMP.index(dref), _derivative_dependencies(pref))
    model = JuMP.owner_model(dref)
    # delete the variable-parameter pair from the anti-duplication dictionary 
    delete!(model.deriv_lookup, (vref, pref))
    # delete any derivative constraints associated with this derivative 
    delete_derivative_constraints(dref)
    # delete associated point variables and mapping
    for index in _point_variable_dependencies(dref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated reduced variables and mapping
    for index in _reduced_variable_dependencies(dref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated derivative variables and mapping 
    for index in _derivative_dependencies(dref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    return
end
