################################################################################
#                       INFINITE PARAMETER FUNCTION METHODS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::ParameterFunctionIndex
                               )::ParameterFunctionRef
    return ParameterFunctionRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::ParameterFunctionData
                          )::ParameterFunctionIndex
    return MOIUC.add_item(model.param_functions, object)
end

# Extend _data_dictionary (reference based)
function _data_dictionary(fref::ParameterFunctionRef
    )::MOIUC.CleverDict{ParameterFunctionIndex, ParameterFunctionData{ParameterFunction{GeneralVariableRef}}}
    return JuMP.owner_model(fref).param_functions
end

# Extend _data_object
function _data_object(fref::ParameterFunctionRef
    )::ParameterFunctionData{ParameterFunction{GeneralVariableRef}}
  object = get(_data_dictionary(fref), JuMP.index(fref), nothing)
  object === nothing && error("Invalid parameter function reference, cannot find " *
                        "corresponding object in the model. This is likely " *
                        "caused by using the reference of a deleted function.")
  return object
end

# Extend _core_variable_object
function _core_variable_object(fref::ParameterFunctionRef
    )::ParameterFunction{GeneralVariableRef}
    return _data_object(fref).func
end

# Extend _object_numbers
function _object_numbers(fref::ParameterFunctionRef)::Vector{Int}
    return _core_variable_object(fref).object_nums
end

# Extend _parameter_numbers
function _parameter_numbers(fref::ParameterFunctionRef)::Vector{Int}
    return _core_variable_object(fref).parameter_nums
end

"""
    build_parameter_function(_error::Function, func::Function, 
        parameter_refs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}, Tuple}
        )::ParameterFunction{GeneralVariableRef}

Build an [`ParameterFunction`](@ref) object that employs a infinite 
parameter function `func` that takes instances of the infinite parameter(s) as 
input. This can ultimately by incorporated into expressions to enable nonlinear 
infinite parameter behavior and/or incorporate data over infinite domains. Errors 
if the infinite parameter tuple is formatted incorrectly. The allowed format 
follows that of infinite variables. Also errors if the function doesn't accept 
a support realization of the `parameter_refs` as input.

**Example**
```julia-repl 
julia> f = build_parameter_function(error, sin, t)
ParameterFunction{GeneralVariableRef}(sin, (t,), [1], [1])
```
"""
function build_parameter_function(
    _error::Function, 
    func::Function, 
    parameter_refs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}, Tuple}
    )::ParameterFunction{GeneralVariableRef}
    # process the parameter reference inputs
    prefs = VectorTuple(parameter_refs)
    # check the arguments
    _check_parameter_tuple(_error, prefs)
    _check_valid_function(_error, func, prefs)
    # get the parameter object numbers
    object_nums = Int[]
    for pref in prefs 
        union!(object_nums, _object_number(pref))
    end
    # get the parameter numbers 
    param_nums = [_parameter_number(pref) for pref in prefs]
    # make the variable and return
    return ParameterFunction(func, prefs, object_nums, param_nums)
end 

# Used to update parameter-parameter function mappings
function _update_param_var_mapping(fref::ParameterFunctionRef,
                                   prefs::VectorTuple)::Nothing
    for pref in prefs
        dependency_list = _parameter_function_dependencies(pref)
        if !(JuMP.index(fref) in dependency_list)
            push!(dependency_list, JuMP.index(fref))
        end
    end
    return
end

"""
    add_parameter_function(model::InfiniteModel, pfunc::ParameterFunction, 
                           [name::String = ""])::GeneralVariableRef

Add an [`ParameterFunction`](@ref) `pfunc` to the `model` using `name` for 
printing and return a `GeneralVariableRef` such that it can be embedded in 
expressions. Errors if the infinite parameters `pfunc` points to do not belong to 
`model`.

**Example**
```julia-repl
julia> f = build_parameter_function(error, sin, t);

julia> fref = add_parameter_function(model, f, "sin")
sin(t)
```
"""
function add_parameter_function(
    model::InfiniteModel, 
    pfunc::ParameterFunction, 
    name::String = ""
    )::GeneralVariableRef
    _check_parameters_valid(model, pfunc.parameter_refs)
    if isempty(name) 
        name = String(first(methods(pfunc.func)).name)
    end
    data_object = ParameterFunctionData(pfunc, name)
    findex = _add_data_object(model, data_object)
    fref = ParameterFunctionRef(model, findex)
    _update_param_var_mapping(fref, pfunc.parameter_refs)
    return _make_variable_ref(model, findex)
end

"""
    JuMP.name(fref::ParameterFunctionRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to return the base name of
`fref`.

**Example**
```julia-repl
julia> name(fref)
"func_name"
```
"""
function JuMP.name(fref::ParameterFunctionRef)::String 
    object = get(_data_dictionary(fref), JuMP.index(fref), nothing)
    return object === nothing ? "" : object.name
end

"""
    JuMP.set_name(fref::ParameterFunctionRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of infinite parameter functions.

**Example**
```julia-repl
julia> set_name(fref, "func_name")

julia> name(fref)
"func_name"
```
"""
function JuMP.set_name(fref::ParameterFunctionRef, name::String)::Nothing
    _data_object(fref).name = name
    return
end

"""
    parameter_function(func::Function, 
                       pref_inputs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}, Tuple}; 
                       [name::String = ""])::GeneralVariableRef

Make an infinite parameter function and return a `GeneralVariableRef` that can be 
embedded in InfiniteOpt expressions. This serves as a convenient wrapper for 
[`build_parameter_function`](@ref) and [`add_parameter_function`](@ref). Here 
`func` denotes the function that will take a support of infinite parameters as 
input (formatted like `pref_inputs`) and will return a scalar value. Errors if 
`func` will not take such as support as input or if `pref_inputs` follow an invalid 
input format.

**Example**
```julia-repl 
julia> p_func = parameter_function(sin, t)
sin(t)
```
"""
function parameter_function(func::Function, pref_inputs; 
                            name::String = "")::GeneralVariableRef
    f = build_parameter_function(error, func, pref_inputs)
    model = JuMP.owner_model(first(f.parameter_refs))
    return add_parameter_function(model, f, name)
end

"""
    raw_parameter_refs(fref::ParameterFunctionRef)::VectorTuple{GeneralVariableRef}

Return the raw [`VectorTuple`](@ref) of the parameter references that `fref`
depends on. This is primarily an internal method where
[`parameter_refs`](@ref parameter_refs(::ParameterFunctionRef))
is intended as the preferred user function.
"""
function raw_parameter_refs(fref::ParameterFunctionRef)::VectorTuple{GeneralVariableRef}
    return _core_variable_object(fref).parameter_refs 
end

"""
    parameter_refs(fref::ParameterFunctionRef)::Tuple

Return the parameter references associated with `fref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `fref`.

**Example**
```julia-repl
julia> parameter_refs(p_func)
(t,)
```
"""
function parameter_refs(fref::ParameterFunctionRef)::Tuple
    return Tuple(raw_parameter_refs(fref))
end

"""
    parameter_list(fref::ParameterFunctionRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `fref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(::ParameterFunctionRef))
is intended as the preferred user function.
"""
function parameter_list(fref::ParameterFunctionRef)::Vector{GeneralVariableRef}
    return raw_parameter_refs(fref).values
end

"""
    raw_function(fref::ParameterFunctionRef)::Function

Returns the raw function behind `fref` that takes a particular support of `fref`'s 
infinite parameters as input. 
"""
function raw_function(fref::ParameterFunctionRef)::Function
    return _core_variable_object(fref).func
end

# Extend _semi_infinite_variable_dependencies
function _semi_infinite_variable_dependencies(
    fref::ParameterFunctionRef
     )::Vector{SemiInfiniteVariableIndex}
    return _data_object(fref).semi_infinite_var_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(
    fref::ParameterFunctionRef
     )::Vector{DerivativeIndex}
    return _data_object(fref).derivative_indices
end

# Extend _measure_dependencies
function _measure_dependencies(fref::ParameterFunctionRef)::Vector{MeasureIndex}
    return _data_object(fref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(fref::ParameterFunctionRef)::Vector{ConstraintIndex}
    return _data_object(fref).constraint_indices
end

"""
    used_by_semi_infinite_variable(fref::ParameterFunctionRef)::Bool

Return a `Bool` indicating if `fref` is used by a semi-infinite infinite variable.

**Example**
```julia-repl
julia> used_by_semi_infinite_variable(fref)
false
```
"""
function used_by_semi_infinite_variable(fref::ParameterFunctionRef)::Bool
    return !isempty(_semi_infinite_variable_dependencies(fref))
end

"""
    used_by_derivative(fref::ParameterFunctionRef)::Bool

Return a `Bool` indicating if `fref` is used by a derivative.

**Example**
```julia-repl
julia> used_by_derivative(vref)
true
```
"""
function used_by_derivative(fref::ParameterFunctionRef)::Bool
    return !isempty(_derivative_dependencies(fref))
end

"""
    used_by_measure(fref::ParameterFunctionRef)::Bool

Return a `Bool` indicating if `fref` is used by a measure.

**Example**
```julia-repl
julia> used_by_measure(fref)
true
```
"""
function used_by_measure(fref::ParameterFunctionRef)::Bool
    return !isempty(_measure_dependencies(fref))
end

"""
    used_by_constraint(fref::ParameterFunctionRef)::Bool

Return a `Bool` indicating if `fref` is used by a constraint.

**Example**
```julia-repl
julia> used_by_constraint(fref)
false
```
"""
function used_by_constraint(fref::ParameterFunctionRef)::Bool
    return !isempty(_constraint_dependencies(fref))
end

"""
    is_used(fref::ParameterFunctionRef)::Bool

Return a `Bool` indicating if `fref` is used in the model.

**Example**
```julia-repl
julia> is_used(fref)
true
```
"""
function is_used(fref::ParameterFunctionRef)::Bool
    return used_by_measure(fref) || used_by_constraint(fref) || 
           used_by_semi_infinite_variable(fref) || used_by_derivative(fref)
end

"""
    JuMP.delete(model::InfiniteModel, fref::ParameterFunctionRef)::Nothing

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)) to delete
infinite parameter functions and their dependencies. Errors if `fref` is invalid,
meaning it has already been deleted or it belongs to another model.
"""
function JuMP.delete(model::InfiniteModel, fref::ParameterFunctionRef)::Nothing 
    @assert JuMP.is_valid(model, fref) "Infinite parameter function is invalid."
    # update the optimizer model status
    if is_used(fref)
        set_optimizer_model_ready(model, false)
    end
    # update parameter mapping
    all_prefs = parameter_list(fref)
    for pref in all_prefs
        filter!(e -> e != JuMP.index(fref), _parameter_function_dependencies(pref))
    end
    gvref = _make_variable_ref(model, JuMP.index(fref))
    # delete associated semi-infinite variables and mapping
    for index in _semi_infinite_variable_dependencies(fref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated derivative variables and mapping 
    for index in _derivative_dependencies(fref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # remove from measures if used
    for mindex in _measure_dependencies(fref)
        mref = dispatch_variable_ref(model, mindex)
        func = measure_function(mref)
        data = measure_data(mref)
        if func isa GeneralVariableRef
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_meas = Measure(new_func, data, Int[], Int[], true)
        else
            _remove_variable(func, gvref)
            new_meas = build_measure(func, data)
        end
        _set_core_variable_object(mref, new_meas)
    end
    # remove from constraints if used
    for cindex in _constraint_dependencies(fref)
        cref = _temp_constraint_ref(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_constraint_object(cref, new_constr)
            empty!(_object_numbers(cref))
        else
            _remove_variable(func, gvref)
            _data_object(cref).object_nums = sort(_object_numbers(func))
        end
    end
    # delete the data object
    _delete_data_object(fref)
    return
end

################################################################################
#                              COMPARISON METHODS
################################################################################
## Extend for better comparisons than default
# GenericAffExpr
function Base.:(==)(aff1::JuMP.GenericAffExpr{C, V},
    aff2::JuMP.GenericAffExpr{C, V}
    )::Bool where {C, V <: GeneralVariableRef}
    return aff1.constant == aff2.constant && aff1.terms == aff2.terms
end

# GenericQuadExpr
function Base.:(==)(quad1::JuMP.GenericQuadExpr{C, V},
    quad2::JuMP.GenericQuadExpr{C, V}
    )::Bool where {C, V <: GeneralVariableRef}
    pairs1 = collect(pairs(quad1.terms))
    pairs2 = collect(pairs(quad2.terms))
    if length(pairs1) != length(pairs2)
        return false
    end
    for i in eachindex(pairs1)
        if pairs1[i][1].a != pairs2[i][1].a || pairs1[i][1].b != pairs2[i][1].b ||
           pairs1[i][2] != pairs2[i][2]
            return false
        end
    end
    return quad1.aff == quad2.aff
end

################################################################################
#                            VARIABLE LIST MAKING
################################################################################
## Determine which variables are present in a function
# GeneralVariableRef
function _all_function_variables(
    f::GeneralVariableRef
    )::Vector{GeneralVariableRef}
    return [f]
end

# GenericAffExpr
function _all_function_variables(
    f::JuMP.GenericAffExpr{C, V}
    )::Vector{V} where {C, V <: GeneralVariableRef}
    return collect(keys(f.terms))
end

# GenericQuadExpr
function _all_function_variables(
    f::JuMP.GenericQuadExpr{C, V}
    )::Vector{V} where {C, V <: GeneralVariableRef}
    vref_set = Set(keys(f.aff.terms))
    for pair in keys(f.terms)
        push!(vref_set, pair.a)
        push!(vref_set, pair.b)
    end
    return collect(vref_set)
end

# Fallback
function _all_function_variables(f)
    error("`_all_function_variables` not defined for expression of type $(typeof(f)).")
end

################################################################################
#                            OBJECT NUMBER METHODS
################################################################################
## Return the unique set of object numbers in an expression
# Dispatch fallback (--> should be defined for each non-empty variable type)
_object_numbers(expr::DispatchVariableRef)::Vector{Int} = Int[]

# GeneralVariableRef
function _object_numbers(expr::GeneralVariableRef)::Vector{Int}
    return _object_numbers(dispatch_variable_ref(expr))
end

# GenericAffExpr
function _object_numbers(expr::JuMP.GenericAffExpr)::Vector{Int}
    obj_nums = Set{Int}()
    for vref in keys(expr.terms)
        union!(obj_nums, _object_numbers(vref))
    end
    return collect(obj_nums)
end

# GenericQuadExpr
function _object_numbers(expr::JuMP.GenericQuadExpr)::Vector{Int}
    obj_nums = Set{Int}()
    for vref in keys(expr.aff.terms)
        union!(obj_nums, _object_numbers(vref))
    end
    for pair in keys(expr.terms)
        union!(obj_nums, _object_numbers(pair.a))
        union!(obj_nums, _object_numbers(pair.b))
    end
    return collect(obj_nums)
end

# Variable list
function _object_numbers(vrefs::Vector{GeneralVariableRef})::Vector{Int}
    obj_nums = Set{Int}()
    for vref in vrefs
        union!(obj_nums, _object_numbers(vref))
    end
    return collect(obj_nums)
end

################################################################################
#                             PARAMETER NUMBER METHODS
################################################################################
## Return the unique set of parameter numbers in an expression
# Dispatch fallback (--> should be defined for each non-empty variable type)
_parameter_numbers(expr::DispatchVariableRef)::Vector{Int} = Int[]

# GeneralVariableRef
function _parameter_numbers(expr::GeneralVariableRef)::Vector{Int}
    return _parameter_numbers(dispatch_variable_ref(expr))
end

# GenericAffExpr
function _parameter_numbers(expr::JuMP.GenericAffExpr)::Vector{Int}
    par_nums = Set{Int}()
    for vref in keys(expr.terms)
        union!(par_nums, _parameter_numbers(vref))
    end
    return collect(par_nums)
end

# GenericQuadExpr
function _parameter_numbers(expr::JuMP.GenericQuadExpr)::Vector{Int}
    par_nums = Set{Int}()
    for vref in keys(expr.aff.terms)
        union!(par_nums, _parameter_numbers(vref))
    end
    for pair in keys(expr.terms)
        union!(par_nums, _parameter_numbers(pair.a))
        union!(par_nums, _parameter_numbers(pair.b))
    end
    return collect(par_nums)
end

################################################################################
#                             MODEL EXTRACTION METHODS
################################################################################
## Get the model from an expression
# GeneralVariableRef
function _model_from_expr(expr::GeneralVariableRef)::InfiniteModel
    return JuMP.owner_model(expr)
end

# AffExpr
function _model_from_expr(expr::JuMP.GenericAffExpr)::Union{InfiniteModel, Nothing}
    if isempty(expr.terms)
        return
    else
        return JuMP.owner_model(first(keys(expr.terms)))
    end
end

# QuadExpr
function _model_from_expr(expr::JuMP.GenericQuadExpr)::Union{InfiniteModel, Nothing}
    result = _model_from_expr(expr.aff)
    if result !== nothing
        return result
    elseif isempty(expr.terms)
        return
    else
        return JuMP.owner_model(first(keys(expr.terms)).a)
    end
end

# Vector{GeneralVariableRef}
function _model_from_expr(vrefs::Vector{GeneralVariableRef})::Union{InfiniteModel, Nothing}
    if isempty(vrefs)
        return
    else
        return JuMP.owner_model(first(vrefs))
    end
end

# Fallback
function _model_from_expr(expr)
    error("`_model_from_expr` not defined for expr of type $(typeof(expr)).")
end

################################################################################
#                            VARIABLE REMOVAL BOUNDS
################################################################################
## Delete variables from an expression
# GenericAffExpr
function _remove_variable(f::JuMP.GenericAffExpr,
                          vref::GeneralVariableRef)::Nothing
    if haskey(f.terms, vref)
        delete!(f.terms, vref)
    end
    return
end

# GenericQuadExpr
function _remove_variable(f::JuMP.GenericQuadExpr,
                          vref::GeneralVariableRef)::Nothing
    _remove_variable(f.aff, vref)
    for (pair, coef) in f.terms # TODO should we preserve non-vref in pair? --> this would differ from JuMP
        if pair.a == vref || pair.b == vref
            delete!(f.terms, pair)
        end
    end
    return
end

################################################################################
#                         COEFFICIENT MODIFICATION METHODS
################################################################################
## Modify linear coefficient of variable in expression
# GeneralVariableRef
function _set_variable_coefficient!(expr::GeneralVariableRef,
    var::GeneralVariableRef,
    coeff::Real
    )::JuMP.GenericAffExpr{Float64, GeneralVariableRef}
    # Determine if variable is that of the expression and change accordingly
    if expr == var
        return Float64(coeff) * var
    else
        return expr + Float64(coeff) * var
    end
end

# GenericAffExpr
function _set_variable_coefficient!(expr::JuMP.GenericAffExpr{C, V},
    var::V,
    coeff::Real
    )::JuMP.GenericAffExpr{C, V} where {C, V <: GeneralVariableRef}
    # Determine if variable is in the expression and change accordingly
    if haskey(expr.terms, var)
        expr.terms[var] = coeff
        return expr
    else
        return expr + coeff * var
    end
end

# GenericQuadExpr
function _set_variable_coefficient!(expr::JuMP.GenericQuadExpr{C, V},
    var::V,
    coeff::Real
    )::JuMP.GenericQuadExpr{C, V} where {C, V <: GeneralVariableRef}
    # Determine if variable is in the expression and change accordingly
    if haskey(expr.aff.terms, var)
        expr.aff.terms[var] = coeff
        return expr
    else
        return expr + coeff * var
    end
end

# Fallback
function _set_variable_coefficient!(expr, var::GeneralVariableRef, coeff::Real)
    error("Unsupported function type for coefficient modification.")
end

################################################################################
#                         PARAMETER REFERNCE METHODS
################################################################################
## Return an element of a parameter reference tuple given the model and index
# IndependentParameterIndex
function _make_param_tuple_element(model::InfiniteModel,
    idx::IndependentParameterIndex,
    )::GeneralVariableRef
    return _make_parameter_ref(model, idx)
end

# DependentParametersIndex
function _make_param_tuple_element(model::InfiniteModel,
    idx::DependentParametersIndex,
    )::Vector{GeneralVariableRef}
    dpref = DependentParameterRef(model, DependentParameterIndex(idx, 1))
    num_params = _num_parameters(dpref)
    return [GeneralVariableRef(model, idx.value, DependentParameterIndex, i)
            for i in 1:num_params]
end

"""
    parameter_refs(expr)::Tuple

Return the tuple of infinite parameter references that determine the infinite
dependencies of `expr`.

**Example**
```julia-repl
julia> parameter_refs(my_expr)
(t,)
```
"""
function parameter_refs(expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr})::Tuple
    model = _model_from_expr(expr)
    if model === nothing
        return ()
    else
        obj_nums = sort!(_object_numbers(expr))
        obj_indices = _param_object_indices(model)[obj_nums]
        return Tuple(_make_param_tuple_element(model, idx) for idx in obj_indices)
    end
end

################################################################################
#                             EXPRESSION SEARCHING
################################################################################
# Check expression for a particular variable type via a recursive search
#=
function _has_variable(vrefs::Vector{GeneralVariableRef},
                       vref::GeneralVariableRef; prior=GeneralVariableRef[]
                       )::Bool
    if vrefs[1] == vref
        return true
    elseif _index_type(vrefs[1]) == MeasureIndex
        dvref = dispatch_variable_ref(vrefs[1])
        if length(vrefs) > 1
            return _has_variable(_all_function_variables(measure_function(dvref)),
                                 vref, prior = append!(prior, vrefs[2:end]))
        else
            return _has_variable(_all_function_variables(measure_function(dvref)),
                                 vref, prior = prior)
        end
    elseif length(vrefs) > 1
        return _has_variable(vrefs[2:end], vref, prior = prior)
    elseif length(prior) > 0
        return _has_variable(prior, vref)
    else
        return false
    end
end
=#
