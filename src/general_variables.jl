################################################################################
#                               BASE EXTENSIONS
################################################################################
# Extend Base.copy for new variable types
Base.copy(v::GeneralVariableRef) = v
Base.copy(v::DispatchVariableRef) = v

# Extend Base.:(==) for GeneralVariableRef
function Base.:(==)(v::T, w::T)::Bool where {T <: GeneralVariableRef}
    return v.model === w.model && v.raw_index == w.raw_index &&
           v.index_type == w.index_type && v.param_index == w.param_index
end

# Extend Base.:(==) for DispatchVariableRef
function Base.:(==)(v::T, w::U)::Bool where {T <: DispatchVariableRef,
                                             U <: DispatchVariableRef}
    return v.model === w.model && v.index == w.index
end

# Extend Base.broadcastable
Base.broadcastable(v::GeneralVariableRef) = Ref(v)
Base.broadcastable(v::DispatchVariableRef) = Ref(v)

# Extend Base.length
Base.length(v::GeneralVariableRef)::Int = 1
Base.length(v::DispatchVariableRef)::Int = 1

# Extend JuMP functions
JuMP.isequal_canonical(v::GeneralVariableRef, w::GeneralVariableRef)::Bool = v == w
JuMP.isequal_canonical(v::DispatchVariableRef, w::DispatchVariableRef)::Bool = v == w

# Extract the root name of a variable reference (removes the bracketed container indices)
function _remove_name_index(vref::GeneralVariableRef)::String
    name = JuMP.name(vref)
    first_bracket = findfirst(isequal('['), name)
    if first_bracket === nothing
        return name
    else
        # Hacky fix to handle invalid Unicode
        try
            return name[1:first_bracket-1]
        catch
            return name[1:first_bracket-2]
        end
    end
end

# Define basic attribute getters
_index_type(vref::GeneralVariableRef)::DataType = vref.index_type
_raw_index(vref::GeneralVariableRef)::Int = vref.raw_index
_param_index(vref::GeneralVariableRef)::Int = vref.param_index

################################################################################
#                          BASIC REFERENCE ACCESSERS
################################################################################
"""
    JuMP.index(vref::GeneralVariableRef)::AbstractInfOptIndex

Extend [`JuMP.index`](@ref JuMP.index(::JuMP.VariableRef)) to return the
appropriate index of `vref`.

**Example**
```julia-repl
julia> index(vref)
FiniteVariableIndex(1)
```
"""
function JuMP.index(vref::GeneralVariableRef)::AbstractInfOptIndex
    index_type = _index_type(vref)
    if index_type == DependentParameterIndex
        return index_type(DependentParametersIndex(_raw_index(vref)),
                           _param_index(vref))
    else
        return index_type(_raw_index(vref))
    end
end

"""
    JuMP.index(vref::DispatchVariableRef)::AbstractInfOptIndex

Extend [`JuMP.index`](@ref JuMP.index(::JuMP.VariableRef)) to return the
appropriate index of `vref`.
"""
function JuMP.index(vref::DispatchVariableRef)::AbstractInfOptIndex
    return vref.index
end

"""
    JuMP.owner_model(vref::GeneralVariableRef)::InfiniteModel

Extend [`JuMP.owner_model`](@ref JuMP.owner_model(::JuMP.AbstractVariableRef)) to
return the model where `vref` is stored.

**Example**
```julia-repl
julia> owner_model(vref)
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variable: 1
Derivatives: 0
Measures: 0
`FiniteVariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
`FiniteVariableRef`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Names registered in the model: vref
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
function JuMP.owner_model(vref::GeneralVariableRef)::InfiniteModel
    return vref.model
end

"""
    JuMP.owner_model(vref::DispatchVariableRef)::InfiniteModel

Extend [`JuMP.owner_model`](@ref JuMP.owner_model(::JuMP.AbstractVariableRef)) to
return the model where `vref` is stored.
"""
function JuMP.owner_model(vref::DispatchVariableRef)::InfiniteModel
    return vref.model
end

################################################################################
#                          DISPATCH VARIABLE MAKERS
################################################################################
"""
    dispatch_variable_ref(model::InfiniteModel, index::AbstractInfOptIndex)

Return the variable reference associated the type of `index`. This needs to be
defined for each variable reference type.
"""
function dispatch_variable_ref end

"""
    dispatch_variable_ref(vef::GeneralVariableRef)::DispatchVariableRef

Return the concrete [`DispatchVariableRef`](@ref) this associated with `vref`.
This relies on `dispatch_variable_ref` being extended for the index type,
otherwise an `MethodError` is thrown.
"""
function dispatch_variable_ref(vref::GeneralVariableRef)::DispatchVariableRef
    model = JuMP.owner_model(vref)
    idx = JuMP.index(vref)
    return dispatch_variable_ref(model, idx)
end

################################################################################
#                          CORE DATA METHODS
################################################################################
"""
    _add_data_object(model::InfiniteModel, object::AbstractDataObject)::ObjectIndex

Add `object` to the appropriate `CleverDict` in `model` and return the its
index. This needs to be defined for the type of `object`. These definitions
need to use `MOIUC.add_item` to add the object to the `CleverDict`.
"""
function _add_data_object end

"""
    _data_dictionary(vref::DispatchVariableRef)::MOIUC.CleverDict

Return the `CleverDict` that stores data objects for the type of `vref`. This
needs to be defined for the type of `vref`.
"""
function _data_dictionary end

"""
    _data_dictionary(vref::GeneralVariableRef)::MOIUC.CleverDict

Return the `CleverDict` that stores data objects for the type of `vref`. It
relies on `_data_dictionary` being defined
for the underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _data_dictionary(vref::GeneralVariableRef)::MOIUC.CleverDict
    return _data_dictionary(dispatch_variable_ref(vref))
end

"""
    _data_object(vref::DispatchVariableRef)::AbstractDataObject

Return the data object associated with `vref`, in other words the object its
index points to in the `InfiniteModel`. This needs to be defined for the type
of `vref`. This should use `_data_dictionary` to access the `CleverDict` that
the object is stored in.
"""
function _data_object end

"""
    _data_object(vref::GeneralVariableRef)::AbstractDataObject

Return the data object associated with `vref`, in other words the object its
index points to in the `InfiniteModel`. It relies on `_data_object` being defined
for the underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _data_object(vref::GeneralVariableRef)::AbstractDataObject
    return _data_object(dispatch_variable_ref(vref))
end

"""
    _delete_data_object(vref::DispatchVariableRef)::Nothing

Delete the concrete `AbstractDataObject` associated with `vref`.
"""
function _delete_data_object(vref::DispatchVariableRef)::Nothing
    delete!(_data_dictionary(vref), JuMP.index(vref))
    return
end

################################################################################
#                                 NAME METHODS
################################################################################
# Dispatch fallback
function JuMP.name(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.name` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.name(vref::GeneralVariableRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to
return the name of `vref`. It relies on `JuMP.name` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.name(vref::GeneralVariableRef)::String
    return JuMP.name(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_name(vref::DispatchVariableRef, name::String)
    throw(ArgumentError("`JuMP.set_name` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_name(vref::GeneralVariableRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.name(::JuMP.VariableRef, ::String)) to
set the name of `vref`. It relies on `JuMP.set_name` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_name(vref::GeneralVariableRef, name::String)::Nothing
    return JuMP.set_name(dispatch_variable_ref(vref), name)
end

################################################################################
#                             VALIDITY METHODS
################################################################################
"""
    JuMP.is_valid(model::InfiniteModel, vref::DispatchVariableRef)::Bool

Extend [`JuMP.is_valid`](@ref JuMP.is_valid(::JuMP.Model, ::JuMP.VariableRef)) to
return `Bool` if `vref` is a valid reference.
"""
function JuMP.is_valid(model::InfiniteModel, vref::DispatchVariableRef)::Bool
    return model === JuMP.owner_model(vref) &&
           haskey(_data_dictionary(vref), JuMP.index(vref))
end

# DependentParameterRef
function JuMP.is_valid(model::InfiniteModel, vref::DependentParameterRef)::Bool
    return model === JuMP.owner_model(vref) &&
           haskey(_data_dictionary(vref), JuMP.index(vref).object_index)
end

"""
    JuMP.is_valid(model::InfiniteModel, vref::GeneralVariableRef)::Bool

Extend [`JuMP.is_valid`](@ref JuMP.is_valid(::JuMP.Model, ::JuMP.VariableRef)) to
return `Bool` if `vref` is a valid reference.

**Example**
```julia-repl
julia> is_valid(model, vref)
true
```
"""
function JuMP.is_valid(model::InfiniteModel, vref::GeneralVariableRef)::Bool
    return JuMP.is_valid(model, dispatch_variable_ref(vref))
end

################################################################################
#                             CORE OBJECT METHODS
################################################################################
"""
    _core_variable_object(vref::DispatchVariableRef)::Union{InfOptParameter, InfOptVariable, Measure}

Return the core object that `vref` points to. This needs to be extended for type
of `vref`. This should use `_data_object` to access the data object where the
variable object is stored.
"""
function _core_variable_object end

"""
    _core_variable_object(vref::GeneralVariableRef)::Union{InfOptParameter, InfOptVariable, Measure}

Return the core object that `vref` points to. This is enabled
with appropriate definitions of `_core_variable_object` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _core_variable_object(vref::GeneralVariableRef)
    return _core_variable_object(dispatch_variable_ref(vref))
end

"""
    _set_core_variable_object(vref::DispatchVariableRef, object)::Nothing

Sets the core object that `vref` points to `object`. This needs to be extended
for types of `vref` and `object`. This should use `_data_object` to access the
data object where the variable object is stored.
"""
function _set_core_variable_object end

################################################################################
#                               DEPENDENCY METHODS
################################################################################
# Define wrappers for internal usage methods and their templates
for op = (:_infinite_variable_dependencies, :_semi_infinite_variable_dependencies,
          :_point_variable_dependencies, :_measure_dependencies,
          :_constraint_dependencies, :_derivative_dependencies, 
          :_derivative_constraint_dependencies, :_parameter_function_dependencies,
          :_generative_measures)
    @eval begin
        # define the api template
        func = $op
        """
            $func(vref::DispatchVariableRef)::Vector{AbstractInfOptIndex}

        Return the indices of these entities that depend on `vref`. This needs to
        be extended for type of `vref`. This should use `_data_object` to access the
        data object where the name is stored if appropriate.
        """
        function $(op) end
        # define the dispatch version
        """
            $func(vref::GeneralVariableRef)::Vector{AbstractInfOptIndex}

        Return the indices of these entities that depend on `vref`. This is enabled
        with appropriate definitions of `$func` for the
        underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
        """
        function $op(vref::GeneralVariableRef)
            return $op(dispatch_variable_ref(vref))
        end
    end
end

################################################################################
#                              USED BY METHODS
################################################################################
# Define the usage method wrappers and their fallbacks
for op = (:used_by_infinite_variable, :used_by_semi_infinite_variable,
          :used_by_point_variable, :used_by_measure, :used_by_constraint,
          :used_by_objective, :used_by_derivative, :is_used, :has_derivative_constraints,
          :used_by_parameter_function)
    @eval begin
        # define the fallback method
        func = $op
        function $op(vref::DispatchVariableRef)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(vref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            $func(vref::GeneralVariableRef)::Bool

        Define `$func` for general variable references. It relies on `$func`
        being defined for the underlying `DispatchVariableRef`, otherwise an
        `ArugmentError` is thrown. See the underlying docstrings for more
        information.
        """
        function $op(vref::GeneralVariableRef)::Bool
            return $op(dispatch_variable_ref(vref))
        end
    end
end

################################################################################
#                              DELETE METHODS
################################################################################
# Dispatch fallback
function JuMP.delete(model::InfiniteModel, vref)
    throw(ArgumentError("`JuMP.delete` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.delete(model::InfiniteModel, vref::GeneralVariableRef)::Nothing

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)) to
delete `vref` and its dependencies. It relies on `JuMP.delete`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.delete(model::InfiniteModel, vref::GeneralVariableRef)::Nothing
    return JuMP.delete(model, dispatch_variable_ref(vref))
end

"""
    JuMP.delete(model::InfiniteModel,
                prefs::AbstractArray{<:GeneralVariableRef})::Nothing

Extend `JuMP.delete` to delete a group of dependent infinite parameters and
their dependencies. An `ArugmentError` is thrown if `prefs` are not dependent
infinite parameters.
"""
function JuMP.delete(model::InfiniteModel,
                     prefs::AbstractArray{<:GeneralVariableRef})::Nothing
    return JuMP.delete(model, dispatch_variable_ref.(prefs))
end

################################################################################
#                              PARAMETER METHODS
################################################################################
"""
    _parameter_number(pref::DispatchVariableRef)::Int

Return the parameter creation number for `pref` assuming it is an infinite
parameter. This needs to be defined for the type of `pref`. This should use
the `_data_object` to get the number.
"""
function _parameter_number end

"""
    _parameter_number(pref::GeneralVariableRef)::Int

Return the parameter creation number for `pref` assuming it is an infinite
parameter. It relies on `_parameter_number` being properly defined for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _parameter_number(pref::GeneralVariableRef)::Int
    return _parameter_number(dispatch_variable_ref(pref))
end

"""
    _object_number(pref::DispatchVariableRef)::Int

Return the object number for `pref` assuming it is an infinite
parameter. This needs to be defined for the type of `pref`. This should use
the `_data_object` to get the number.
"""
function _object_number end

"""
    _object_number(pref::GeneralVariableRef)::Int

Return the object number for `pref` assuming it is an infinite
parameter. It relies on `_object_number` being properly defined for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _object_number(pref::GeneralVariableRef)::Int
    return _object_number(dispatch_variable_ref(pref))
end

# Define 1 argument user method wrappers and their fallbacks
for op = (:infinite_set, :num_supports, :significant_digits, :has_supports,
          :supports, :delete_supports, :fill_in_supports!, :parameter_value,
          :derivative_method, :has_generative_supports, :has_internal_supports,
          :add_generative_supports, :raw_function, :generative_support_info)
    @eval begin
        # define the fallback method
        func = $op
        function $op(pref; kwargs...)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(pref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            $func(prefs; [kwargs...])

        Define `$func` for general variable references. It relies on `$func`
        being defined for the underlying `DispatchVariableRef`, otherwise an
        `ArugmentError` is thrown. See the underlying docstrings for more
        information. Note that this is a auto generated wrapper and the underlying
        method may or may not use `kwargs`.
        """
        function $op(prefs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}; kwargs...)
            return $op(dispatch_variable_ref.(prefs); kwargs...)
        end
    end
end

# Dispatch fallback
function set_infinite_set(pref, set::AbstractInfiniteSet)
    throw(ArgumentError("`set_infinite_set` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    set_infinite_set(pref::GeneralVariableRef, set::InfiniteScalarSet)::Nothing

Specify the scalar infinite set of the infinite parameter `pref` to `set`. Note
this will reset/delete all the supports contained in the
underlying parameter object. Also, errors if `pref` is used
by a measure. An `ArgumentError` is thrown if `pref` is not an infinite parameter.
"""
function set_infinite_set(pref::GeneralVariableRef,
                          set::InfiniteScalarSet)::Nothing
    return set_infinite_set(dispatch_variable_ref(pref), set)
end

"""
    set_infinite_set(prefs::AbstractArray{<:GeneralVariableRef},
                     set::InfiniteArraySet)::Nothing

Specify the multi-dimensional infinite set of the dependent infinite parameters
`prefs` to `set`. Note this will reset/delete all the supports contained in the
underlying [`DependentParameters`](@ref) object. This will error if the not all
of the dependent infinite parameters are included or if any of them are used by
measures. An `ArgumentError` is thrown if `prefs` are not dependent infinite
parameters.
"""
function set_infinite_set(prefs::AbstractArray{<:GeneralVariableRef},
                          set::InfiniteArraySet)::Nothing
    return set_infinite_set(dispatch_variable_ref.(prefs), set)
end

# Better fallbacks for supports
function supports(pref::DispatchVariableRef; kwargs...)
    throw(ArgumentError("`supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end
function supports(prefs::AbstractArray; kwargs...)
    throw(ArgumentError("`supports` not defined for variable reference type(s) " *
                        "`$(typeof(prefs))`."))
end

# Dispatch fallback
function set_supports(pref, supports; kwargs...)
    throw(ArgumentError("`set_supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    set_supports(pref::GeneralVariableRef, supports::Union{Real, Vector{<:Real}};
                 [force::Bool = false])::Nothing

Set the support points associated with a single infinite
parameter `pref`. An `ArgumentError` is thrown if `pref` is not an independent
infinite parameter.
"""
function set_supports(pref::GeneralVariableRef,
                      supports::Union{Real, Vector{<:Real}};
                      force::Bool = false, 
                      label::Type{<:AbstractSupportLabel} = UserDefined
                      )::Nothing
    return set_supports(dispatch_variable_ref(pref), supports,
                        force = force, label = label)
end

"""
    set_supports(
        prefs::Union{Vector{GeneralVariableRef}, AbstractArray{<:GeneralVariableRef}},
        supports::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}};
        [force::Bool = false]
        )::Nothing

Set the support points associated with dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function set_supports(prefs::AbstractArray{<:GeneralVariableRef},
                      supports::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}};
                      label::Type{<:AbstractSupportLabel} = UserDefined, 
                      force::Bool = false
                      )::Nothing
    return set_supports(dispatch_variable_ref.(prefs), supports, label = label,
                        force = force)
end

# Dispatch fallback
function add_supports(pref, supports; kwargs...)
    throw(ArgumentError("`add_supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    add_supports(pref::GeneralVariableRef,
                 supports::Union{Real, Vector{<:Real}})::Nothing

Add the support points `supports` to a single infinite
parameter `pref`. An `ArgumentError` is thrown if `pref` is not an independent
infinite parameter.
"""
function add_supports(pref::GeneralVariableRef,
                      supports::Union{Real, Vector{<:Real}};
                      check::Bool = true, 
                      label::Type{<:AbstractSupportLabel} = UserDefined
                      )::Nothing
    return add_supports(dispatch_variable_ref(pref), supports,
                        check = check, label = label)
end

"""
    add_supports(
        prefs::Union{Vector{GeneralVariableRef}, AbstractArray{<:GeneralVariableRef}},
        supports::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}}
        )::Nothing

Add the support points `supports` to the dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function add_supports(prefs::AbstractArray{<:GeneralVariableRef},
                      supports::Union{Array{<:Real, 2}, Vector{<:AbstractArray{<:Real}}};
                      label::Type{<:AbstractSupportLabel} = UserDefined, 
                      check::Bool = true
                      )::Nothing
    return add_supports(dispatch_variable_ref.(prefs), supports, label = label,
                        check = check)
end

# Fallback
function JuMP.set_value(vref::DispatchVariableRef, value::Real)
    throw(ArgumentError("`JuMP.set_value` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_value(vref::DispatchVariableRef, value::Real)::Nothing

Extend `JuMP.set_value` to
set the value of `vref`. It relies on `JuMP.set_value`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.set_value(vref::GeneralVariableRef, value::Real)::Nothing
    return JuMP.set_value(dispatch_variable_ref(vref), value)
end

# Dispatch fallback
function set_derivative_method(pref::DispatchVariableRef, method)
    throw(ArgumentError("`set_derivative_method` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    set_derivative_method(pref::GeneralVariableRef,
                          method::AbstractDerivativeMethod
                          )::Nothing

Specify the numerical derivative evaluation technique associated with `pref`.
An `ArgumentError` is thrown if `pref` is not an infinite parameter.
"""
function set_derivative_method(pref::GeneralVariableRef,
                               method::AbstractDerivativeMethod
                               )::Nothing
    return set_derivative_method(dispatch_variable_ref(pref), method)
end

# Define parameter status setters 
for op = (:_set_has_generative_supports, :_set_has_internal_supports,
          :_set_has_derivative_constraints)
    @eval begin
        # define the fallback method
        func = $op
        function $op(vref::DispatchVariableRef, status)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(vref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        function $op(vref::GeneralVariableRef, status::Bool)::Nothing
            return $op(dispatch_variable_ref(vref), status)
        end
    end
end

# Dispatch fallback
function call_function(fref::DispatchVariableRef, support...)
    throw(ArgumentError("`call_function` not defined for variable reference type(s) " *
                        "`$(typeof(fref))`."))
end

"""
    call_function(fref::GeneralVariableRef, support...)::Float64

Call the parameter function of `fref` at `support`.
An `ArgumentError` is thrown if `fref` is not a parameter function.
"""
function call_function(fref::GeneralVariableRef, support...)::Float64
    return call_function(dispatch_variable_ref(fref), support...)
end

################################################################################
#                              VARIABLE METHODS
################################################################################
# Define single argument variable method wrappers and their fallbacks
for op = (:raw_parameter_refs, :parameter_refs, :parameter_list,
          :start_value_function, :reset_start_value_function,
          :infinite_variable_ref, :eval_supports, :raw_parameter_values,
          :parameter_values, :parameter_bounds, :delete_parameter_bounds)
    @eval begin
        # define the fallback method
        func = $op
        function $op(vref::DispatchVariableRef)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(vref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            $func(vref::GeneralVariableRef)

        Define `$func` for general variable references. It relies on `$func`
        being defined for the underlying `DispatchVariableRef`, otherwise an
        `ArugmentError` is thrown. See the underlying docstrings for more
        information.
        """
        function $op(vref::GeneralVariableRef)
            return $op(dispatch_variable_ref(vref))
        end
    end
end

# Dispatch fallback
function set_start_value_function(vref::DispatchVariableRef, start)
    throw(ArgumentError("`set_start_value_function` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    set_start_value_function(vref::GeneralVariableRef, start::Union{Real, Function})::Nothing

Set the start value function of `vref`. It relies on `set_start_value_function`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function set_start_value_function(vref::GeneralVariableRef, start)::Nothing
    return set_start_value_function(dispatch_variable_ref(vref), start)
end

"""
    has_parameter_bounds(vref::GeneralVariableRef)::Bool

Return a `Bool` indicating if `vref` is limited to a sub-domain as defined
by parameter bound.
"""
function has_parameter_bounds(vref::GeneralVariableRef)::Bool
    return has_parameter_bounds(dispatch_variable_ref(vref))
end

# Dispatch fallback
function set_parameter_bounds(vref::DispatchVariableRef, bounds; kwargs...)
    throw(ArgumentError("`set_parameter_bounds` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    set_parameter_bounds(vref::GeneralVariableRef,
                         bounds::ParameterBounds{GeneralVariableRef};
                         [force::Bool = false])::Nothing

Specify a new set of parameter bounds for a finite variable `vref`.
An `ArgumentError` is thrown if `vref` is not a finite variable.
"""
function set_parameter_bounds(vref::GeneralVariableRef,
                              bounds::ParameterBounds{GeneralVariableRef};
                              force::Bool = false, _error::Function = error
                              )::Nothing
    return set_parameter_bounds(dispatch_variable_ref(vref), bounds,
                                force = force, _error = _error)
end

# Dispatch fallback
function add_parameter_bounds(vref::DispatchVariableRef, bounds; kwargs...)
    throw(ArgumentError("`add_parameter_bounds` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    add_parameter_bounds(vref::GeneralVariableRef,
                         bounds::ParameterBounds{GeneralVariableRef}
                         )::Nothing

Specify more parameter bounds for a finite variable `vref`.
An `ArgumentError` is thrown if `vref` is not a finite variable.
"""
function add_parameter_bounds(vref::GeneralVariableRef,
                              bounds::ParameterBounds{GeneralVariableRef};
                              _error::Function = error
                              )::Nothing
    return add_parameter_bounds(dispatch_variable_ref(vref), bounds,
                                _error = _error)
end

################################################################################
#                               MEASURE METHODS
################################################################################
# Define measure queries and their fallbacks
for op = (:measure_function, :measure_data, :is_analytic, :expand)
    @eval begin
        # define the fallback method
        func = $op
        function $op(mref::DispatchVariableRef)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(mref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            $func(mref::GeneralVariableRef)

        Define `$func` for general variable references. Errors if `mref` does
        not correspond to a `MeasureRef`. See the underlying docstrings for more
        information.
        """
        function $op(mref::GeneralVariableRef)
            return $op(dispatch_variable_ref(mref))
        end
    end
end

################################################################################
#                               DERIVATIVE METHODS
################################################################################
# Define measure queries and their fallbacks
for op = (:derivative_argument, :operator_parameter, :evaluate, 
          :derivative_constraints, :delete_derivative_constraints)
    @eval begin
        # define the fallback method
        func = $op
        function $op(dref::DispatchVariableRef)
            str = string("`", func, "` not defined for variable reference type " *
                         "`$(typeof(dref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            $func(dref::GeneralVariableRef)

        Define `$func` for general variable references. Errors if `dref` does
        not correspond to a `DerivativeRef`. See the underlying docstrings for more
        information.
        """
        function $op(dref::GeneralVariableRef)
            return $op(dispatch_variable_ref(dref))
        end
    end
end

################################################################################
#                            VARIABLE INFO METHODS
################################################################################
# Define the 1 argument JuMP variable info methods
for op = (:has_lower_bound, :has_upper_bound, :is_fixed, :is_binary, :is_integer,
          :lower_bound, :upper_bound, :fix_value, :start_value,
          :set_binary, :set_integer,
          :LowerBoundRef, :UpperBoundRef, :FixRef, :BinaryRef, :IntegerRef,
          :delete_lower_bound, :delete_upper_bound, :unfix, :unset_binary,
          :unset_integer)
    @eval begin
        # define the fallback method
        func = JuMP.$op
        function JuMP.$op(vref::DispatchVariableRef)
            str = string("`JuMP.", func, "` not defined for variable reference type " *
                         "`$(typeof(vref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            JuMP.$func(vref::GeneralVariableRef)

        Define `JuMP.$func` for general variable references. It relies on `JuMP.$func`
        being defined for the underlying `DispatchVariableRef`, otherwise an
        `ArugmentError` is thrown. See the underlying docstrings for more
        information.
        """
        function JuMP.$op(vref::GeneralVariableRef)
            return JuMP.$op(dispatch_variable_ref(vref))
        end
    end
end

# Define the 2 argument setting methods (except for fix)
for op = (:set_lower_bound, :set_upper_bound, :set_start_value)
    @eval begin
        # define the fallback method
        func = JuMP.$op
        function JuMP.$op(vref::DispatchVariableRef, value)
            str = string("`JuMP.", func, "` not defined for variable reference type " *
                         "`$(typeof(vref))`.")
            throw(ArgumentError(str))
        end
        # define the dispatch version
        """
            JuMP.$func(vref::GeneralVariableRef, value::Real)::Nothing

        Define `JuMP.$func` for general variable references. It relies on `JuMP.$func`
        being defined for the underlying `DispatchVariableRef`, otherwise an
        `ArugmentError` is thrown. See the underlying docstrings for more
        information.
        """
        function JuMP.$op(vref::GeneralVariableRef, value::Real)::Nothing
            return JuMP.$op(dispatch_variable_ref(vref), value)
        end
    end
end

# Dispatch fallback for JuMP.fix
function JuMP.fix(vref::DispatchVariableRef, value::Real; force::Bool = false)
    throw(ArgumentError("`JuMP.fix` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.fix(vref::GeneralVariableRef, value::Real; force::Bool = false)::Nothing

Define `JuMP.fix` for general variable references. It relies on `JuMP.fix`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown. See the underlying docstrings for more
information.
"""
function JuMP.fix(vref::GeneralVariableRef, value::Real;
                  force::Bool = false)::Nothing
    return JuMP.fix(dispatch_variable_ref(vref), value, force = force)
end
