################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::InfiniteVariableIndex
    )::InfiniteVariableRef
    return InfiniteVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:InfiniteVariable}
    )::InfiniteVariableIndex
    return MOIUC.add_item(model.infinite_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel, 
    ::Type{InfiniteVariable}
    )::MOIUC.CleverDict{InfiniteVariableIndex, VariableData{<:InfiniteVariable}}
    return model.infinite_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(
    vref::InfiniteVariableRef
    )::MOIUC.CleverDict{InfiniteVariableIndex, VariableData{<:InfiniteVariable}}
    return JuMP.owner_model(vref).infinite_vars
end

# Extend _data_object
function _data_object(vref::InfiniteVariableRef)
    object = Base.get(_data_dictionary(vref), JuMP.index(vref), nothing)
    if isnothing(object) 
        error("Invalid infinite variable reference, cannot find ",
              "corresponding variable in the model. This is likely ",
              "caused by using the reference of a deleted variable.")
    end
    return object
end

# Extend _core_variable_object
function _core_variable_object(vref::InfiniteVariableRef)
    return _data_object(vref).variable
end

# Extend _object_numbers
function _object_numbers(vref::InfiniteVariableRef)::Vector{Int}
    return _core_variable_object(vref).object_nums
end

# Extend _parameter_numbers
function _parameter_numbers(vref::InfiniteVariableRef)::Vector{Int}
    return _core_variable_object(vref).parameter_nums
end

# Define getter function for var.is_vector_start
function _is_vector_start(vref::InfiniteVariableRef)::Bool
    return _core_variable_object(vref).is_vector_start
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(
    vref::InfiniteVariableRef, 
    var::V, 
    data::VariableData{V}
    )::Nothing where {V <: InfiniteVariable}
    data.variable = var
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(
    vref::InfiniteVariableRef, 
    var::V1, 
    data::VariableData{V2}
    )::Nothing  where {V1, V2}
    new_data = VariableData(var, data.name, data.lower_bound_index,
                            data.upper_bound_index, data.fix_index,
                            data.zero_one_index, data.integrality_index,
                            data.measure_indices, data.constraint_indices,
                            data.in_objective, data.point_var_indices,
                            data.semi_infinite_var_indices, 
                            data.derivative_indices, data.deriv_constr_indices)
    _data_dictionary(vref)[JuMP.index(vref)] = new_data
    return
end

# Extend _set_core_variable_object for InfiniteVariableRefs
function _set_core_variable_object(
    vref::InfiniteVariableRef,
    var::InfiniteVariable
    )::Nothing
    _adaptive_data_update(vref, var, _data_object(vref))
    return
end

################################################################################
#                                 MODEL CHECKING
################################################################################
# TODO(odow): this is a work-around for a bug in Julia v1.9. I've no idea why it
# is necessary.
function JuMP.check_belongs_to_model(
    x::InfiniteVariableRef,
    model::InfiniteModel
    )
    if JuMP.owner_model(x) !== model
        throw(JuMP.VariableNotOwned(x))
    end
    return
end

################################################################################
#                          DEFINTION HELPER METHODS
################################################################################
"""
    Infinite{VT <: VectorTuple} <: InfOptVariableType

A `DataType` to assist in making infinite variables. This can be passed as an 
extra argument to `@variable` to make an infinite variable: 
```julia 
@variable(model, var_expr, Infinite(parameter_refs...), args..., kwargs...)
```
Here `parameter_refs` can be a single parameter reference, a single parameter 
array with parameters defined in the same macro call, or multiple arguments where 
each argument is either of the first two options listed.

**Fields**
- `parameter_refs::VT`: The infinite parameters the variable will depend on.
"""
struct Infinite{VT <: Collections.VectorTuple} <: InfOptVariableType
    parameter_refs::VT
    function Infinite(vt::VT) where {VT <: Collections.VectorTuple}
        return new{VT}(vt)
    end
end
function Infinite(args...)
    return Infinite(Collections.VectorTuple(args))
end

## Check that each parameter tuple element is formatted correctly
# IndependentParameterRefs
function _check_tuple_element(
    _error::Function,
    prefs::Vector{IndependentParameterRef}
    )::Nothing
    return
end

# DependentParameterRefs
function _check_tuple_element(
    _error::Function,
    prefs::Vector{DependentParameterRef}
    )::Nothing
    if length(prefs) != _num_parameters(first(prefs))
        _error("Infinite parameter tuple elements cannot depend on a subset ",
               "of dependent parameters.")
    end
    return
end

# Fallback
function _check_tuple_element(_error::Function, prefs)
    _error("Cannot have mixed parameter types in a tuple element and can only ",
           "specify infinite parameters.")
end

# Check parameter tuple, ensure all elements contain parameter references
function _check_parameter_tuple(
    _error::Function,
    raw_prefs::Collections.VectorTuple{GeneralVariableRef}
    )::Nothing
    allunique(raw_prefs) || _error("Cannot double specify infinite parameter ",
                                   "references.")
    for i in 1:size(raw_prefs, 1)
        prefs = dispatch_variable_ref.(raw_prefs[i, :])
        _check_tuple_element(_error, prefs)
    end
    return
end
function _check_parameter_tuple(
    _error::Function, 
    raw_prefs::Collections.VectorTuple{T}
    ) where {T}
    _error("Expected infinite parameter reference arguments, but got arguments ",
           "of type $T. This may be because `Infinite()` was given (i.e., no ",
           "infinite parameter references were specified as arguments).")
end

## Check and format the variable info considering functional start values
# Just a number given for the start value
function _check_and_format_infinite_info(
    _error::Function,
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Real},
    prefs::Collections.VectorTuple
    )
    # prepare the start value function and return the info
    start_func = (s::Vector{<:Real}) -> info.start
    return JuMP.VariableInfo(info.has_lb, convert(Float64, info.lower_bound), 
                             info.has_ub, convert(Float64, info.upper_bound),
                             info.has_fix, convert(Float64, info.fixed_value), 
                             !isnan(info.start), start_func, info.binary, 
                             info.integer), 
           true
end

# make function for checking start value function properties 
function _check_valid_function(
    _error::Function, 
    func::Function, 
    prefs::Collections.VectorTuple
    )::Nothing 
    input_format = typeof(Tuple(Vector{Float64}(undef, length(prefs)), prefs))
    if !hasmethod(func, input_format)
        _error("Specified function `$func` must be able to accept a `Float64` " * 
               "support realization of the infinite parameter tuple $(prefs).")
    end
    return 
end

# A function is given for the start value generation
function _check_and_format_infinite_info(
    _error::Function,
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, F},
    prefs::Collections.VectorTuple
    )::Tuple{JuMP.VariableInfo{Float64, Float64, Float64, F}, Bool} where {F <: Function}
    # check the function properties
    _check_valid_function(_error, info.start, prefs)
    # make the info and return
    return JuMP.VariableInfo(info.has_lb, convert(Float64, info.lower_bound), 
                             info.has_ub, convert(Float64, info.upper_bound),
                             info.has_fix, convert(Float64, info.fixed_value), 
                             info.has_start, info.start, info.binary, 
                             info.integer), 
           false
end

# Fallback
function _check_and_format_infinite_info(
    _error::Function,
    info::JuMP.VariableInfo,
    prefs::Collections.VectorTuple
    )
    _error("Unrecognized formatting for the variable information.")
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, 
                        var_type::Infinite)::InfiniteVariable{GeneralVariableRef}

Build and return an infinite variable based on `info` and `var_type`. Errors if 
the infinite parameter references included in `var_type` are invalid. See 
[`Infinite`](@ref) for more information.

**Example**
```julia-repl
julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite(t));
```
"""
function JuMP.build_variable(
    _error::Function, 
    info::JuMP.VariableInfo, 
    var_type::Infinite;
    extra_kwargs...
    )
    # check for unneeded keywords
    for (kwarg, _) in extra_kwargs
        _error("Keyword argument $kwarg is not for use with infinite variables.")
    end
    # get the parameter_refs
    prefs = var_type.parameter_refs
    # check the VectorTuple for validity and format
    _check_parameter_tuple(_error, prefs)
    # check and format the info (accounting for start value functions)
    new_info, is_vect_func = _check_and_format_infinite_info(_error, info, prefs)
    # get the parameter object numbers
    object_nums = Int[]
    for pref in prefs 
        union!(object_nums, _object_number(pref))
    end
    # make the variable and return
    return InfiniteVariable(new_info, prefs,
                            [_parameter_number(pref) for pref in prefs],
                            object_nums, is_vect_func)
end

# check the pref tuple contains only valid parameters
function _check_parameters_valid(
    model::InfiniteModel,
    prefs::Collections.VectorTuple
    )::Nothing
    for pref in prefs
        JuMP.check_belongs_to_model(pref, model)
    end
    return
end

# Used to update parameter-infinite variable mappings
function _update_param_var_mapping(
    vref::InfiniteVariableRef,
    prefs::Collections.VectorTuple
    )::Nothing
    for pref in prefs
        dependency_list = _infinite_variable_dependencies(pref)
        if !(JuMP.index(vref) in dependency_list)
            push!(dependency_list, JuMP.index(vref))
        end
    end
    return
end

"""
    JuMP.add_variable(model::InfiniteModel, var::InfiniteVariable,
                      [name::String = ""])::GeneralVariableRef

Extend the `JuMP.add_variable` function to accomodate infinite variable 
types. Adds a variable to an infinite model `model` and returns a 
[`GeneralVariableRef`](@ref). Primarily intended to be an internal function of 
the constructor macro `@variable`. However, it can be used in combination with
`JuMP.build_variable` to add  infinite variables to an infinite model object.
Errors if invalid parameter reference(s) are included in `var`.

**Example**
```julia-repl
julia> @infinite_parameter(m, t in [0, 10]);

julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite(t));

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)
```
"""
function JuMP.add_variable(
    model::InfiniteModel,
    v::InfiniteVariable,
    name::String = ""
    )::GeneralVariableRef 
    _check_parameters_valid(model, v.parameter_refs)
    data_object = VariableData(v, name)
    vindex = _add_data_object(model, data_object)
    vref = InfiniteVariableRef(model, vindex)
    _update_param_var_mapping(vref, v.parameter_refs)
    gvref = _make_variable_ref(model, vindex)
    _set_info_constraints(v.info, gvref, vref)
    model.name_to_var = nothing
    return gvref
end

################################################################################
#                              RESTRICTION METHODS
################################################################################
## Dispatch functions for functional syntax of making restricted variables
# Point variable
function _restrict_infinite_variable(
    ivref::GeneralVariableRef, 
    vt::Collections.VectorTuple{<:Real}
    )::GeneralVariableRef
    info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, 
                             false, false)
    new_var = JuMP.build_variable(error, info, Point(ivref, vt))
    return JuMP.add_variable(JuMP.owner_model(ivref), new_var, 
                             update_info = false)
end

# Infinite variable
function _restrict_infinite_variable(
    ivref::GeneralVariableRef, 
    vt::Collections.VectorTuple{<:GeneralVariableRef}
    )::GeneralVariableRef
    if parameter_list(ivref) != vt.values
        error("Unrecognized syntax for infinite variable restriction.")
    end
    @warn("Unnecessary use of functional infinite variable restriction syntax " *
          "that will cause performance degredations. This was probably caused " *
          "by using syntax like `y(t, x)` inside expressions. Instead just " *
          "use the infinite variable reference (e.g. `y`).")
    return ivref
end

# Semi-Infinite variable
function _restrict_infinite_variable(
    ivref::GeneralVariableRef, 
    vt::Collections.VectorTuple
    )::GeneralVariableRef
    info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, 
                             false, false)
    new_var = JuMP.build_variable(error, info, SemiInfinite(ivref, vt))
    return JuMP.add_variable(JuMP.owner_model(ivref), new_var)
end

"""
    restrict(ivref::GeneralVariableRef, supps...)::GeneralVariableRef

Restrict the input domain of an infinite variable/derivative `ivref` in 
accordance with the infinite parameters and/or values `supps`. Here `supps` must 
match the formatting of `ivref`'s infinite parameters. Here the following 
outputs are possible:
- Equivalent to `@variable(model, variable_type = Point(ivref, supps...)` if 
  `supps` are a complete support point
- Equivalent to `@variable(model, variable_type = SemiInfinite(ivref, supps...)` 
  if `supps` are a partial support point.

Conveniently, we can also invoke this method by calling `ivref(supps...)`.

Errors if ivref is not an infinite variable or derivative or the formatting of 
`supps` is incorrect. Will warn if `supps` only contain infinite parameters and 
will simply return `ivref`.

**Example**
```julia-repl
julia> restrict(y, 0, x)
y(0, [x[1], x[2]])

julia> restrict(y, 0, [0, 0])
y(0, [0, 0])

julia> y(0, x)
y(0, [x[1], x[2]])

julia> y(0, [0, 0])
y(0, [0, 0])
```
"""
function restrict(ivref::GeneralVariableRef, supps...)::GeneralVariableRef
    idx_type = _index_type(ivref)
    if idx_type != InfiniteVariableIndex && idx_type != DerivativeIndex
        error("The `vref(values..)` restriction syntax is only valid for ",
              "infinite variable and derivative references.")
    end
    return _restrict_infinite_variable(ivref, Collections.VectorTuple(supps))
end

# This enables functional calls (e.g., `y(0, x)` => semi-infinite variable)
function _functional_reference_call(
    ivref::GeneralVariableRef, 
    ::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}}, 
    supps...
    )::GeneralVariableRef
    return _restrict_infinite_variable(ivref, Collections.VectorTuple(supps))
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _semi_infinite_variable_dependencies
function _semi_infinite_variable_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Vector{SemiInfiniteVariableIndex}
    return _data_object(vref).semi_infinite_var_indices
end

# Extend _point_variable_dependencies
function _point_variable_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Vector{PointVariableIndex}
    return _data_object(vref).point_var_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Vector{DerivativeIndex}
    return _data_object(vref).derivative_indices
end

"""
    used_by_semi_infinite_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a semi-infinite variable.

**Example**
```julia-repl
julia> used_by_semi_infinite_variable(vref)
false
```
"""
function used_by_semi_infinite_variable(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Bool
    return !isempty(_semi_infinite_variable_dependencies(vref))
end

"""
    used_by_point_variable(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a point variable.

**Example**
```julia-repl
julia> used_by_point_variable(vref)
false
```
"""
function used_by_point_variable(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Bool
    return !isempty(_point_variable_dependencies(vref))
end

"""
    used_by_derivative(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used by a derivative.

**Example**
```julia-repl
julia> used_by_derivative(vref)
true
```
"""
function used_by_derivative(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    return !isempty(_derivative_dependencies(vref))
end

"""
    is_used(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia-repl
julia> is_used(vref)
false
```
"""
function is_used(vref::Union{InfiniteVariableRef, DerivativeRef})::Bool
    if used_by_measure(vref) || used_by_constraint(vref)
        return true
    end
    for vindex in _point_variable_dependencies(vref)
        if is_used(PointVariableRef(JuMP.owner_model(vref), vindex))
            return true
        end
    end
    for vindex in _semi_infinite_variable_dependencies(vref)
        if is_used(SemiInfiniteVariableRef(JuMP.owner_model(vref), vindex))
            return true
        end
    end
    for dindex in _derivative_dependencies(vref)
        if is_used(DerivativeRef(JuMP.owner_model(vref), dindex))
            return true
        end
    end
    return false
end

################################################################################
#                           PARAMETER REFERENCES
################################################################################
"""
    raw_parameter_refs(vref::InfiniteVariableRef)::VectorTuple

Return the raw [`VectorTuple`](@ref InfiniteOpt.Collections.VectorTuple) of the 
parameter references that `vref` depends on. This is primarily an internal method 
where [`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef)) 
is intended as the preferred user function.
"""
function raw_parameter_refs(vref::InfiniteVariableRef)
    return _core_variable_object(vref).parameter_refs
end

"""
    parameter_refs(vref::InfiniteVariableRef)::Tuple

Return the parameter references associated with the infinite variable `vref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `vref`.

**Example**
```julia-repl
julia> @variable(model, T, Infinite(t))
T(t)

julia> parameter_refs(T)
(t,)
```
"""
function parameter_refs(vref::InfiniteVariableRef)
    return Tuple(raw_parameter_refs(vref))
end

"""
    parameter_list(vref::InfiniteVariableRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `vref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(vref::InfiniteVariableRef))
is intended as the preferred user function.
"""
function parameter_list(vref::InfiniteVariableRef)::Vector{GeneralVariableRef}
    return raw_parameter_refs(vref).values
end

# get parameter list from raw VectorTuple
function parameter_list(
    prefs::Collections.VectorTuple{GeneralVariableRef}
    )::Vector{GeneralVariableRef}
    return prefs.values
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
## format infinite info 
# Good to go
function _format_infinite_info(
    info::I
    )::I where {I <: JuMP.VariableInfo{Float64, Float64, Float64, <:Function}}
    return info
end

# Convert as needed 
function _format_infinite_info(
    info::JuMP.VariableInfo{V, W, T, F}
    )::JuMP.VariableInfo{Float64, Float64, Float64, F} where {V, W, T, F <: Function}
    return JuMP.VariableInfo(info.has_lb, convert(Float64, info.lower_bound), 
                             info.has_ub, convert(Float64, info.upper_bound), 
                             info.has_fix, convert(Float64, info.fixed_value),
                             info.has_start, info.start, info.binary, 
                             info.integer)
end

# Set info for infinite variables
function _update_variable_info(
    vref::InfiniteVariableRef,
    info::JuMP.VariableInfo
    )::Nothing
    prefs = raw_parameter_refs(vref)
    param_nums = _parameter_numbers(vref)
    obj_nums = _object_numbers(vref)
    is_vect_func = _is_vector_start(vref)
    new_var = InfiniteVariable(_format_infinite_info(info), prefs, param_nums, 
                               obj_nums, is_vect_func)
    _set_core_variable_object(vref, new_var)
    return
end

# Specify start_value fallback for infinite variables
function JuMP.start_value(vref::Union{InfiniteVariableRef, DerivativeRef})
    error("`start_value` not defined for infinite variables, consider calling " *
          "`start_value_function` instead.")
end

# Specify set_start_value fallback for infinite variables
function JuMP.set_start_value(
    vref::Union{InfiniteVariableRef, DerivativeRef}, 
    value::Real
    )
    error("`set_start_value` not defined for infinite variables, consider calling " *
          "`set_start_value_function` instead.")
end

"""
    start_value_function(vref::Union{InfiniteVariableRef, DerivativeRef})::Union{Nothing, Function}

Return the function that is used to generate the start values of `vref` for
particular support values. Returns `nothing` if no start behavior has been
specified.

**Example**
```julia-repl
julia> start_value_function(vref)
my_start_func
```
"""
function start_value_function(
    vref::Union{InfiniteVariableRef, DerivativeRef}
    )::Union{Nothing, Function}
    if _variable_info(vref).has_start
        return _variable_info(vref).start
    else
        return
    end
end

"""
    set_start_value_function(vref::InfiniteVariableRef,
                             start::Union{Real, Function})::Nothing

Set the start value function of `vref`. If `start::Real` then a function is
generated to such that the start value will be `start` for the entire infinite
domain. If `start::Function` then this function should map to a scalar start value
given a support value arguments matching the format of the parameter elements in
`parameter_refs(vref)`.

**Example**
```julia-repl
julia> set_start_value_function(vref, 1) # all start values will be 1

julia> set_start_value_function(vref, my_func) # each value will be made via my_func
```
"""
function set_start_value_function(
    vref::InfiniteVariableRef,
    start::Union{Real, Function}
    )::Nothing
    info = _variable_info(vref)
    set_optimizer_model_ready(JuMP.owner_model(vref), false)
    prefs = raw_parameter_refs(vref)
    temp_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                  info.upper_bound, info.has_fix, info.fixed_value,
                                  true, start, info.binary, info.integer)
    new_info, is_vect_func = _check_and_format_infinite_info(error, temp_info, prefs)
    obj_nums = _object_numbers(vref)
    param_nums = _parameter_numbers(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums, is_vect_func)
    _set_core_variable_object(vref, new_var)
    # TODO update point variable start values as appropriate
    return
end

"""
    reset_start_value_function(vref::InfiniteVariableRef)::Nothing

Remove the existing start value function and return to the default. Generally,
this is triggered by deleting an infinite parameter that `vref` depends on.

**Example**
```julia-repl
julia> reset_start_value_function(vref)
```
"""
function reset_start_value_function(vref::InfiniteVariableRef)::Nothing
    info = _variable_info(vref)
    set_optimizer_model_ready(JuMP.owner_model(vref), false)
    start_func = (s::Vector{<:Real}) -> NaN
    new_info = JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                                 info.upper_bound, info.has_fix, info.fixed_value,
                                 false, start_func, info.binary, info.integer)
    prefs = raw_parameter_refs(vref)
    obj_nums = _object_numbers(vref)
    param_nums = _parameter_numbers(vref)
    new_var = InfiniteVariable(new_info, prefs, param_nums, obj_nums, true)
    _set_core_variable_object(vref, new_var)
    # TODO update point variable start values as appropriate
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::InfiniteVariableRef)::Nothing
    # remove variable info constraints associated with vref
    _delete_info_constraints(vref)
    # update parameter mapping
    all_prefs = parameter_list(vref)
    for pref in all_prefs
        filter!(e -> e != JuMP.index(vref), _infinite_variable_dependencies(pref))
    end
    model = JuMP.owner_model(vref)
    # delete associated point variables and mapping
    for index in _point_variable_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated semi-infinite variables and mapping
    for index in _semi_infinite_variable_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # delete associated derivative variables and mapping 
    for index in _derivative_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    return
end
