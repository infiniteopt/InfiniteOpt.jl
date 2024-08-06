################################################################################
#                            CORE OBJECT DISPATCHES
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::ParameterFunctionIndex
    )
    return ParameterFunctionRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel, object::ParameterFunctionData)
    return MOIUC.add_item(model.param_functions, object)
end

# Extend _data_dictionary (reference based)
function _data_dictionary(fref::ParameterFunctionRef)
    return JuMP.owner_model(fref).param_functions
end

# Extend _data_object
function _data_object(fref::ParameterFunctionRef)
    object = get(_data_dictionary(fref), JuMP.index(fref), nothing)
    if isnothing(object) 
        error("Invalid parameter function reference, cannot find ",
              "corresponding object in the model. This is likely ",
              "caused by using the reference of a deleted function.")
    end
    return object
end

"""
    core_object(fref::ParameterFunctionRef)::ParameterFunction

Retrieve the underlying core [`ParameterFucntion`](@ref) object for `fref`. 
This is intended as an advanced method for developers.
"""
function core_object(fref::ParameterFunctionRef)
    return _data_object(fref).func
end

"""
    parameter_group_int_indices(fref::ParameterFunctionRef)::Vector{Int}

Return the list of infinite parameter group integer indices used by `fref`.
"""
function parameter_group_int_indices(fref::ParameterFunctionRef)
    return core_object(fref).group_int_idxs
end

# Extend _parameter_numbers
function _parameter_numbers(fref::ParameterFunctionRef)
    return core_object(fref).parameter_nums
end

################################################################################
#                                 DEFINITION
################################################################################
# Helper method for checking functions of `ParameterFunction`s
function _check_param_func_method(_error::Function, func::Function, prefs)
    supp_type = typeof(Tuple(Vector{Float64}(undef, length(prefs)), prefs))
    if !hasmethod(func, supp_type)
        _error("Parameter function method `$(func)` is not defined for `$(func)(",
               join(Tuple(prefs), ", ") * ")`. Note that the infinite parameter ",
               "arguments `" * join(Tuple(prefs), ", ") * "` are checked via a ",
               "numeric support (each parameter is a `Float64`) of the same format. ",
               "If you get this message specifying a function for the ",
               "bounds/fix-value/start of a variable, then the function arguments ",
               "Do not properly accomodate the format of the variable's infinite ",
               "parameters.")
    end
    return
end

"""
    build_parameter_function(
        _error::Function, 
        func::Function, 
        parameter_refs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}, Tuple}
        )::ParameterFunction

Build an [`ParameterFunction`](@ref) object that employs a parameter function 
`func` that takes instances of the infinite parameter(s) as input. This can 
ultimately by incorporated into expressions to enable nonlinear infinite parameter 
behavior and/or incorporate data over infinite domains.

Here `func` should be of the form:
```
func(paramvals...)::Float64
```
where the formatting of `paramvals` is analagous to point variables (and will be 
based on the tuple of infinite parameter references given in `parameter_refs`). 

Errors if the infinite parameter tuple is formatted incorrectly. The allowed 
format follows that of infinite variables. Also errors if the function doesn't 
accept a support realization of the `parameter_refs` as input.

**Example**
```julia-repl 
julia> f = build_parameter_function(error, sin, t);
```
"""
function build_parameter_function(
    _error::Function, 
    func::Function, 
    parameter_refs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}, Tuple};
    extra_kwargs...
    )
    # check for unneeded keywords
    for (kwarg, _) in extra_kwargs
        _error("Keyword argument $kwarg is not for use with parameter functions.")
    end
    # process the parameter reference inputs and make checks
    prefs = Collections.VectorTuple(parameter_refs)
    _check_parameter_tuple(_error, prefs)
    _check_param_func_method(_error, func, prefs)
    # get the parameter group integer indices
    group_int_idxs = Int[]
    for pref in prefs 
        union!(group_int_idxs, parameter_group_int_index(pref))
    end
    # get the parameter numbers 
    param_nums = [_parameter_number(pref) for pref in prefs]
    # make the variable and return
    return ParameterFunction(func, prefs, param_nums, group_int_idxs)
end 

# Fallback for weird macro inputs
function build_parameter_function(_error::Function, func, prefs; kwargs...)
    _error("Invalid syntax causing unexpected parsing of the function and " * 
           "infinite parameter arguments.")
end 

# Used to update parameter-parameter function mappings
function _update_param_var_mapping(
    fref::ParameterFunctionRef,
    prefs::Collections.VectorTuple
    )
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
                           [name::String])::GeneralVariableRef

Add an [`ParameterFunction`](@ref) `pfunc` to the `model` using `name` for 
printing and return a `GeneralVariableRef` such that it can be embedded in 
expressions. Errors if the parameter function `pfunc` points to do not belong to 
`model`. Note that `pfunc` should be created using [`build_parameter_function`](@ref).

**Example**
```julia-repl
julia> f = build_parameter_function(error, sin, t);

julia> fref = add_parameter_function(model, f)
sin(t)
```
"""
function add_parameter_function(
    model::InfiniteModel, 
    pfunc::ParameterFunction, 
    name::String = String(nameof(pfunc.func))
    )
    _check_parameters_valid(model, pfunc.parameter_refs)
    data_object = ParameterFunctionData(pfunc, name)
    findex = _add_data_object(model, data_object)
    fref = ParameterFunctionRef(model, findex)
    _update_param_var_mapping(fref, pfunc.parameter_refs)
    return GeneralVariableRef(model, findex)
end

"""
    parameter_function(func::Function, 
                       pref_inputs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}, Tuple}; 
                       [name::String = [the name of `func`]]
                       )::GeneralVariableRef

Make a parameter function and return a `GeneralVariableRef` that can be 
embedded in InfiniteOpt expressions. This serves as a convenient wrapper for 
[`build_parameter_function`](@ref) and [`add_parameter_function`](@ref). For an 
even more convenient definition method see [`@parameter_function`](@ref).

Here `func` denotes the function that will take a support of infinite parameters as 
input (formatted like `pref_inputs`) and will return a scalar value. Specifically, 
`func` should be of the form:
```
func(paramvals...)::Float64
```
where the formatting of `paramvals` is analagous to point variables (and will be 
based on the tuple of infinite parameter references given in `parameter_refs`).
Moreover, `func` must be a function that returns a scalar numeric value. 

Errors if `func` will not take a support formatted like `pref_inputs` in 
combination with the `fargs` and `fkwargs` specified. Also errors if `pref_inputs` 
follow an invalid input format.

**Example**
```julia-repl 
julia> p_func = parameter_function(sin, t)
sin(t)

julia> p_func3 = parameter_function((t_supp) -> 2 * sin(2 * t_supp), t, name = "mysin")
mysin(t)

julia> p_func4 = parameter_function(t, name = "mysin") do t_supp
                    if t_supp <= 5
                        return sin(t_supp)
                    else 
                        return 2 * sin(2 * t_supp)
                    end
                 end

mysin(t)
```
"""
function parameter_function(
    func::Function, 
    pref_inputs; 
    name::String = String(nameof(func))
    )
    f = build_parameter_function(error, func, pref_inputs)
    model = JuMP.owner_model(first(f.parameter_refs))
    return add_parameter_function(model, f, name)
end

################################################################################
#                               GETTERS/SETTERS
################################################################################
"""
    JuMP.name(fref::ParameterFunctionRef)::String

Extend `JuMP.name` to return the base name of `fref`.

**Example**
```julia-repl
julia> name(fref)
"func_name"
```
"""
function JuMP.name(fref::ParameterFunctionRef)
    object = get(_data_dictionary(fref), JuMP.index(fref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(fref::ParameterFunctionRef, name::String)::Nothing

Extend `JuMP.set_name` to set the name of a parameter function.

**Example**
```julia-repl
julia> set_name(fref, "func_name")

julia> name(fref)
"func_name"
```
"""
function JuMP.set_name(fref::ParameterFunctionRef, name::String)
    _data_object(fref).name = name
    return
end

"""
    raw_parameter_refs(fref::ParameterFunctionRef)::VectorTuple

Return the raw [`VectorTuple`](@ref InfiniteOpt.Collections.VectorTuple) of the 
parameter references that `fref` depends on. This is primarily an internal method 
where [`parameter_refs`](@ref parameter_refs(::ParameterFunctionRef)) 
is intended as the preferred user function.
"""
function raw_parameter_refs(fref::ParameterFunctionRef)
    return core_object(fref).parameter_refs 
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
function parameter_refs(fref::ParameterFunctionRef)
    return Tuple(raw_parameter_refs(fref))
end

"""
    parameter_list(fref::ParameterFunctionRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `fref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(::ParameterFunctionRef))
is intended as the preferred user function.
"""
function parameter_list(fref::ParameterFunctionRef)
    return raw_parameter_refs(fref).values
end

"""
    raw_function(fref::ParameterFunctionRef)::Function

Returns the raw function behind `fref` that takes a particular support of `fref`'s 
infinite parameters as input. 
"""
function raw_function(fref::ParameterFunctionRef)
    return core_object(fref).func
end

"""
    call_function(fref::ParameterFunctionRef, support...)::Float64

Safely evaluates the [`raw_function`](@ref) of `fref` at a particular support `support`
point that matches the format of the infinite parameter tuple given when the `fref` 
was defined. This is essentially equivalent to `raw_function(fref)(supps...)`. 
"""
function call_function(fref::ParameterFunctionRef, supps...)
    pfunc = core_object(fref)
    return pfunc.func(supps...)
end

# Extend _semi_infinite_variable_dependencies
function _semi_infinite_variable_dependencies(fref::ParameterFunctionRef)
    return _data_object(fref).semi_infinite_var_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(fref::ParameterFunctionRef)
    return _data_object(fref).derivative_indices
end

# Extend _measure_dependencies
function _measure_dependencies(fref::ParameterFunctionRef)
    return _data_object(fref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(fref::ParameterFunctionRef)
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
function used_by_semi_infinite_variable(fref::ParameterFunctionRef)
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
function used_by_derivative(fref::ParameterFunctionRef)
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
function used_by_measure(fref::ParameterFunctionRef)
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
function used_by_constraint(fref::ParameterFunctionRef)
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
function is_used(fref::ParameterFunctionRef)
    return used_by_measure(fref) || used_by_constraint(fref) || 
           used_by_semi_infinite_variable(fref) || used_by_derivative(fref)
end

################################################################################
#                                   DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel, fref::ParameterFunctionRef)::Nothing

Extend `JuMP.delete` to delete parameter functions and their dependencies. Errors 
if `fref` is invalid, meaning it has already been deleted or it belongs to 
another model.
"""
function JuMP.delete(model::InfiniteModel, fref::ParameterFunctionRef)
    @assert JuMP.is_valid(model, fref) "Parameter function is invalid."
    # update the transformation backend status
    if is_used(fref)
        set_transformation_backend_ready(model, false)
    end
    # update parameter mapping
    all_prefs = parameter_list(fref)
    for pref in all_prefs
        filter!(e -> e != JuMP.index(fref), _parameter_function_dependencies(pref))
    end
    gvref = GeneralVariableRef(model, JuMP.index(fref))
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
        _set_core_object(mref, new_meas)
    end
    # remove from constraints if used
    for cindex in copy(_constraint_dependencies(fref))
        cref = InfOptConstraintRef(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_object(cref, new_constr)
            empty!(parameter_group_int_indices(cref))
        elseif func isa AbstractArray && any(isequal(gvref), func)
            JuMP.delete(model, cref)
        else
            _remove_variable(func, gvref)
            _data_object(cref).group_int_idxs = sort(parameter_group_int_indices(func))
        end
    end
    # delete the data object
    _delete_data_object(fref)
    return
end