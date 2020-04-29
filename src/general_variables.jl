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
JuMP.variable_type(model::InfiniteModel)::DataType = GeneralVariableRef

# Extract the root name of a variable reference (removes the bracketed container indices)
function _remove_name_index(vref::GeneralVariableRef)::String
    name = JuMP.name(vref)
    first_bracket = findfirst(isequal('['), name)
    if first_bracket == nothing
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
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref))
julia> index(vref)
HoldVariableIndex(1)
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
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref <= 1))
julia> owner_model(vref)
An InfiniteOpt Model
Feasibility problem with:
Variable: 1
`HoldVariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
`HoldVariableRef`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
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

Delete the concrete `AbstractDataObject` associated with `vref`. Will error if
`vref` is a `DependentParameterRef`.
"""
function _delete_data_object(vref::DispatchVariableRef)::Nothing
    delete!(_data_dictionary(vref), JuMP.index(vref))
    return
end

# Error for DependentParameterRefs
function _delete_data_object(vref::DependentParameterRef)
    error("Cannot delete the data object only given a single parameter $vref.")
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
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @hold_variable(model, vref))
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
function _core_variable_object(vref::GeneralVariableRef)::Union{InfOptParameter, InfOptVariable, Measure}
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
#                             DEPENDENCY METHODS
################################################################################
"""
    _infinite_variable_dependencies(vref::DispatchVariableRef)::Vector{InfiniteVariableIndex}

Return the indices of infinite variables that depend on `vref`. This needs to
be extended for type of `vref`. This should use `_data_object` to access the
data object where the name is stored if appropriate.
"""
function _infinite_variable_dependencies end

"""
    _infinite_variable_dependencies(vref::GeneralVariableRef)::Vector{InfiniteVariableIndex}

Return the indices of infinite variables that depend on `vref`. This is enabled
with appropriate definitions of `_infinite_variable_dependencies` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _infinite_variable_dependencies(vref::GeneralVariableRef)::Vector{InfiniteVariableIndex}
    return _infinite_variable_dependencies(dispatch_variable_ref(vref))
end

"""
    _reduced_variable_dependencies(vref::DispatchVariableRef)::Vector{ReducedInfiniteVariableIndex}

Return the indices of reduced variables that depend on `vref`. This needs to
be extended for type of `vref`. This should use `_data_object` to access the
data object where the name is stored if appropriate.
"""
function _reduced_variable_dependencies end

"""
    _reduced_variable_dependencies(vref::GeneralVariableRef)::Vector{ReducedInfiniteVariableIndex}

Return the indices of reduced variables that depend on `vref`. This is enabled
with appropriate definitions of `_reduced_variable_dependencies` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _reduced_variable_dependencies(vref::GeneralVariableRef)::Vector{ReducedInfiniteVariableIndex}
    return _reduced_variable_dependencies(dispatch_variable_ref(vref))
end

"""
    _point_variable_dependencies(vref::DispatchVariableRef)::Vector{PointVariableIndex}

Return the indices of point variables that depend on `vref`. This needs to
be extended for type of `vref`. This should use `_data_object` to access the
data object where the name is stored if appropriate.
"""
function _point_variable_dependencies end

"""
    _point_variable_dependencies(vref::GeneralVariableRef)::Vector{PointVariableIndex}

Return the indices of point variables that depend on `vref`. This is enabled
with appropriate definitions of `_point_variable_dependencies` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _point_variable_dependencies(vref::GeneralVariableRef)::Vector{PointVariableIndex}
    return _point_variable_dependencies(dispatch_variable_ref(vref))
end

"""
    _measure_dependencies(vref::DispatchVariableRef)::Vector{MeasureIndex}

Return the indices of measures that depend on `vref`. This needs to
be extended for type of `vref`. This should use `_data_object` to access the
data object where the name is stored if appropriate.
"""
function _measure_dependencies end

"""
    _measure_dependencies(vref::GeneralVariableRef)::Vector{MeasureIndex}

Return the indices of measures that depend on `vref`. This is enabled
with appropriate definitions of `_measure_dependencies` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _measure_dependencies(vref::GeneralVariableRef)::Vector{MeasureIndex}
    return _measure_dependencies(dispatch_variable_ref(vref))
end

"""
    _constraint_dependencies(vref::DispatchVariableRef)::Vector{ConstraintIndex}

Return the indices of constraints that depend on `vref`. This needs to
be extended for type of `vref`. This should use `_data_object` to access the
data object where the name is stored if appropriate.
"""
function _constraint_dependencies end

"""
    _constraint_dependencies(vref::GeneralVariableRef)::Vector{ConstraintIndex}

Return the indices of constraints that depend on `vref`. This is enabled
with appropriate definitions of `_constraint_dependencies` for the
underlying `DispatchVariableRef`, otherwise an `MethodError` is thrown.
"""
function _constraint_dependencies(vref::GeneralVariableRef)::Vector{ConstraintIndex}
    return _constraint_dependencies(dispatch_variable_ref(vref))
end

################################################################################
#                              USED BY METHODS
################################################################################
# Dispatch fallback
function used_by_infinite_variable(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_infinite_variable` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_infinite_variable(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by an infinite variable. This requires that
`used_by_infinite_variable` be defined for the underlying dispatch variable.
"""
function used_by_infinite_variable(vref::GeneralVariableRef)::Bool
    return used_by_infinite_variable(dispatch_variable_ref(vref))
end

# Dispatch fallback
function used_by_reduced_variable(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_reduced_variable` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_reduced_variable(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by a reduced variable. This requires that
`used_by_reduced_variable` be defined for the underlying dispatch variable.
"""
function used_by_reduced_variable(vref::GeneralVariableRef)::Bool
    return used_by_reduced_variable(dispatch_variable_ref(vref))
end

# Dispatch fallback
function used_by_point_variable(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_point_variable` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_point_variable(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by a point variable. This requires that
`used_by_point_variable` be defined for the underlying dispatch variable.
"""
function used_by_point_variable(vref::GeneralVariableRef)::Bool
    return used_by_point_variable(dispatch_variable_ref(vref))
end

# Dispatch fallback
function used_by_measure(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_measure` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_measure(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by a measure. This requires that
`used_by_measure` be defined for the underlying dispatch variable.
"""
function used_by_measure(vref::GeneralVariableRef)::Bool
    return used_by_measure(dispatch_variable_ref(vref))
end

# Dispatch fallback
function used_by_objective(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_objective` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_objective(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by an objective. This requires that
`used_by_objective` be defined for the underlying dispatch variable.
"""
function used_by_objective(vref::GeneralVariableRef)::Bool
    return used_by_objective(dispatch_variable_ref(vref))
end

# Dispatch fallback
function used_by_constraint(vref::DispatchVariableRef)
    throw(ArgumentError("`used_by_constraint` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    used_by_constraint(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used by a constraint. This requires that
`used_by_constraint` be defined for the underlying dispatch variable.
"""
function used_by_constraint(vref::GeneralVariableRef)::Bool
    return used_by_constraint(dispatch_variable_ref(vref))
end

# Dispatch fallback
function is_used(vref::DispatchVariableRef)
    throw(ArgumentError("`is_used` not defined for " *
                        "variable reference type `$(typeof(vref))`."))
end

"""
    is_used(vref::GeneralVariableRef)::Bool

Return `Bool` whether `vref` is used. This requires that
`is_used` be defined for the underlying dispatch variable.
"""
function is_used(vref::GeneralVariableRef)::Bool
    return is_used(dispatch_variable_ref(vref))
end

################################################################################
#                              DELETE METHODS
################################################################################
# Dispatch fallback
function JuMP.delete(model::InfiniteModel, vref::DispatchVariableRef)
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

# Dispatch fallback
function infinite_set(pref)
    throw(ArgumentError("`infinite_set` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    infinite_set(pref::GeneralVariableRef)::InfiniteScalarSet

Return the infinite set associated with `pref`. Errors if `pref` is not an
infinite parameter or if an [`InfiniteScalarSet`](@ref) is not defined for `pref`.
"""
function infinite_set(pref::GeneralVariableRef)::InfiniteScalarSet
    return infinite_set(dispatch_variable_ref(pref))
end

"""
    infinite_set(prefs::AbstractArray{<:GeneralVariableRef})::InfiniteArraySet

Return the multi-dimensional infinite set associated with `prefs`. Errors if
`prefs` are not an dependent infinite parameters of if they are not the full set
of dependent parameters.
"""
function infinite_set(prefs::AbstractArray{<:GeneralVariableRef})::InfiniteArraySet
    return infinite_set(dispatch_variable_ref.(prefs))
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

# Dispatch fallback
function num_supports(pref; kwargs...)
    throw(ArgumentError("`num_supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    num_supports(pref::GeneralVariableRef; [label::Symbol = All])::Int

Return the number of support points associated with a single infinite
parameter `pref`. Specify a subset of supports via `label` to only count the
supports with `label`. An `ArgumentError` is thrown if `pref` is not an infinite
parameter.
"""
function num_supports(pref::GeneralVariableRef; label::Symbol = All)::Int
    return num_supports(dispatch_variable_ref(pref), label = label)
end

"""
    num_supports(prefs::AbstractArray{<:GeneralVariableRef};
                 [label::Symbol = All])::Int

Return the number of support points associated with dependent infinite
parameters `prefs`. Specify a subset of supports via `label` to only count the
supports with `label`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function num_supports(prefs::AbstractArray{<:GeneralVariableRef};
                      label::Symbol = All)::Int
    return num_supports(dispatch_variable_ref.(prefs), label = label)
end

# Dispatch fallback
function has_supports(pref)
    throw(ArgumentError("`has_supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    has_supports(pref::GeneralVariableRef)::Bool

Return `Bool` if there are support points associated with a single infinite
parameter `pref`. An `ArgumentError` is thrown if `pref` is not an infinite
parameter.
"""
function has_supports(pref::GeneralVariableRef)::Bool
    return has_supports(dispatch_variable_ref(pref))
end

"""
    has_supports(prefs::AbstractArray{<:GeneralVariableRef})::Bool

Return `Bool` if there are support points associated with dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
```
"""
function has_supports(prefs::AbstractArray{<:GeneralVariableRef})::Bool
    return has_supports(dispatch_variable_ref.(prefs))
end

# Dispatch fallback
function supports(pref; kwargs...)
    throw(ArgumentError("`supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    supports(pref::GeneralVariableRef; [label::Symbol = All])::Vector{Float64}

Return the support points associated with a single infinite
parameter `pref`. Specify a subset of supports via `label` to only count the
supports with `label`. An `ArgumentError` is thrown if `pref` is not an infinite
parameter.
"""
function supports(pref::GeneralVariableRef; label::Symbol = All)::Vector{Float64}
    return supports(dispatch_variable_ref(pref), label = label)
end

"""
    supports(prefs::AbstractArray{<:GeneralVariableRef};
             [label::Symbol = All]
             )::Union{AbstractArray{<:Vector{<:Float64}}, Array{Float64, 2}}

Return the support points associated with dependent infinite
parameters `prefs`. Specify a subset of supports via `label` to only count the
supports with `label`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function supports(prefs::AbstractArray{<:GeneralVariableRef};
                  label::Symbol = All
                  )::Union{AbstractArray{<:Vector{<:Float64}}, Array{Float64, 2}}
    return supports(dispatch_variable_ref.(prefs), label = label)
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
                      force::Bool = false, label::Symbol = UserDefined
                      )::Nothing
    return set_supports(dispatch_variable_ref(pref), supports,
                        force = force, label = label)
end

"""
    set_supports(
        prefs::Union{Vector{GeneralVariableRef}, AbstractArray{<:GeneralVariableRef}},
        supports::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}};
        [force::Bool = false]
        )::Nothing

Set the support points associated with dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function set_supports(prefs::AbstractArray{<:GeneralVariableRef},
                      supports::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}};
                      label::Symbol = UserDefined, force::Bool = false
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
                      check::Bool = true, label::Symbol = UserDefined
                      )::Nothing
    return add_supports(dispatch_variable_ref(pref), supports,
                        check = check, label = label)
end

"""
    add_supports(
        prefs::Union{Vector{GeneralVariableRef}, AbstractArray{<:GeneralVariableRef}},
        supports::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}}
        )::Nothing

Add the support points `supports` to the dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
"""
function add_supports(prefs::AbstractArray{<:GeneralVariableRef},
                      supports::Union{Array{<:Real, 2}, AbstractArray{<:Vector{<:Real}}};
                      label::Symbol = UserDefined, check::Bool = true
                      )::Nothing
    return add_supports(dispatch_variable_ref.(prefs), supports, label = label,
                        check = check)
end

# Dispatch fallback
function delete_supports(pref)
    throw(ArgumentError("`delete_supports` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    delete_supports(pref::GeneralVariableRef)::Nothing

Delete the support points associated with a single infinite
parameter `pref`. An `ArgumentError` is thrown if `pref` is not an infinite
parameter.
"""
function delete_supports(pref::GeneralVariableRef)::Nothing
    return delete_supports(dispatch_variable_ref(pref))
end

"""
    delete_supports(prefs::AbstractArray{<:GeneralVariableRef})::Nothing

Delete the support points associated with dependent infinite
parameters `prefs`. An `ArgumentError` is thrown if `prefs` is are not
dependent infinite parameters.
```
"""
function delete_supports(prefs::AbstractArray{<:GeneralVariableRef})::Nothing
    return delete_supports(dispatch_variable_ref.(prefs))
end

# Dispatch fallback
function fill_in_supports!(pref; kwargs...)
    throw(ArgumentError("`fill_in_supports!` not defined for variable reference type(s) " *
                        "`$(typeof(pref))`."))
end

"""
    fill_in_supports!(pref::GeneralVariableRef; [num_supports::Int = 10,
                      sig_figs::Int = 5])::Nothing

Fill in the support points associated with a single infinite
parameter `pref` up to `num_supports`. An `ArgumentError` is thrown if `pref`
is not an infinite parameter.
"""
function fill_in_supports!(pref::GeneralVariableRef; num_supports::Int = 10,
                           sig_figs::Int = 5)::Nothing
    return fill_in_supports!(dispatch_variable_ref(pref),
                             num_supports = num_supports, sig_figs = sig_figs)
end

"""
    fill_in_supports!(prefs::AbstractArray{<:GeneralVariableRef};
                      [num_supports::Int = 10, sig_figs::Int = 5,
                       modify::Bool = true])::Nothing

Fill in the support points associated with dependent infinite
parameters `prefs` up to `num_supports`. An `ArgumentError` is thrown if `prefs`
is are not dependent infinite parameters.
```
"""
function fill_in_supports!(prefs::AbstractArray{<:GeneralVariableRef};
                           num_supports::Int = 10, sig_figs::Int = 5,
                           modify::Bool = true)::Nothing
    return fill_in_supports!(dispatch_variable_ref.(prefs),
                             num_supports = num_supports, sig_figs = sig_figs,
                             modify = modify)
end

# Fallback
function JuMP.value(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.value` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.value(vref::GeneralVariableRef)::Float64

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.VariableRef)) to
return the value of `vref`. It relies on `JuMP.value`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.value(vref::GeneralVariableRef)::Float64
    return JuMP.value(dispatch_variable_ref(vref))
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

################################################################################
#                              VARIABLE METHODS
################################################################################
# Dispatch fallback
function raw_parameter_refs(vref::DispatchVariableRef)
    throw(ArgumentError("`raw_parameter_refs` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    raw_parameter_refs(vref::GeneralVariableRef)::VectorTuple{GeneralVariableRef}

Return the raw [`VectorTuple`](@ref) of the parameter references that `vref`
depends on. An `ArgumentError` is thrown if `vref` is not an infinite/reduced
variable.
"""
function raw_parameter_refs(vref::GeneralVariableRef)::VectorTuple{GeneralVariableRef}
    return raw_parameter_refs(dispatch_variable_ref(vref))
end

# Dispatch fallback
function parameter_refs(vref::DispatchVariableRef)
    throw(ArgumentError("`parameter_refs` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    parameter_refs(vref::GeneralVariableRef)::Tuple

Return the parameter references associated with the infinite variable `vref`. This
is formatted as a Tuple of containing the parameter references as they inputted
to define `vref`. An `ArgumentError` is thrown if `vref` is not an infinite/reduced
variable.
"""
function parameter_refs(vref::GeneralVariableRef)::Tuple
    return parameter_refs(dispatch_variable_ref(vref))
end

# Dispatch fallback
function parameter_list(vref::DispatchVariableRef)
    throw(ArgumentError("`parameter_list` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    parameter_list(vref::GeneralVariableRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `vref` depends on.
An `ArgumentError` is thrown if `vref` is not an infinite/reduced variable.
"""
function parameter_list(vref::GeneralVariableRef)::Vector{GeneralVariableRef}
    return parameter_list(dispatch_variable_ref(vref))
end

# Dispatch fallback
function set_parameter_refs(vref::DispatchVariableRef, prefs::Tuple)
    throw(ArgumentError("`set_parameter_refs` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    set_parameter_refs(vref::GeneralVariableRef, prefs::Tuple)::Nothing

Specify a new parameter reference tuple `prefs` for the infinite variable `vref`.
An `ArgumentError` is thrown if `vref` is not an infinite variable.
"""
function set_parameter_refs(vref::GeneralVariableRef, prefs::Tuple)::Nothing
    return set_parameter_refs(dispatch_variable_ref(vref), prefs)
end

# Dispatch fallback
function add_parameter_ref(vref::DispatchVariableRef,
    pref::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
    )
    throw(ArgumentError("`add_parameter_ref` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    add_parameter_ref(vref::GeneralVariableRef,
        pref::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
        )::Nothing

Add additional parameter reference or group of parameter references to be
associated with the infinite variable `vref`. An `ArgumentError` is thrown if
`vref` is not an infinite variable.
"""
function add_parameter_ref(vref::GeneralVariableRef,
    pref::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
    )::Nothing
    return add_parameter_ref(dispatch_variable_ref(vref), pref)
end

# Dispatch fallback
function infinite_variable_ref(vref::DispatchVariableRef)
    throw(ArgumentError("`infinite_variable_ref` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    infinite_variable_ref(vref::GeneralVariableRef)::GeneralVariableRef

Return the infinite variable reference associated with `vref`. An
`ArgumentError` is thrown if `vref` is not a reduced/point variable.
"""
function infinite_variable_ref(vref::GeneralVariableRef)::GeneralVariableRef
    return infinite_variable_ref(dispatch_variable_ref(vref))
end

# Dispatch fallback
function eval_supports(vref::DispatchVariableRef)
    throw(ArgumentError("`eval_supports` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    eval_supports(vref::GeneralVariableRef)::Dict{Int, Float64}

Return the evaluation supports associated with the reduced infinite variable
`vref`. An `ArgumentError` is thrown if `vref` is not a reduced variable.
"""
function eval_supports(vref::GeneralVariableRef)::Dict{Int, Float64}
    return eval_supports(dispatch_variable_ref(vref))
end

# Dispatch fallback
function raw_parameter_values(vref::DispatchVariableRef)
    throw(ArgumentError("`raw_parameter_values` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    raw_parameter_values(vref::GeneralVariableRef)::Vector{Float64}

Return the raw support point values associated with the point variable `vref`.
An `ArgumentError` is thrown if `vref` is not a point variable.
"""
function raw_parameter_values(vref::GeneralVariableRef)::Vector{Float64}
    return raw_parameter_values(dispatch_variable_ref(vref))
end

# Dispatch fallback
function parameter_values(vref::DispatchVariableRef)
    throw(ArgumentError("`parameter_values` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    parameter_values(vref::GeneralVariableRef)::Tuple

Return the support point associated with the point variable `vref`.
An `ArgumentError` is thrown if `vref` is not a point variable.
"""
function parameter_values(vref::GeneralVariableRef)::Tuple
    return parameter_values(dispatch_variable_ref(vref))
end

# Dispatch fallback
function parameter_bounds(vref::DispatchVariableRef)
    throw(ArgumentError("`parameter_bounds` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    parameter_bounds(vref::GeneralVariableRef)::ParameterBounds{GeneralVariableRef}

Return the [`ParameterBounds`](@ref) object associated with the hold variable
`vref`. An `ArgumentError` is thrown if `vref` is not a hold variable.
"""
function parameter_bounds(vref::GeneralVariableRef
                          )::ParameterBounds{GeneralVariableRef}
    return parameter_bounds(dispatch_variable_ref(vref))
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

Specify a new set of parameter bounds for a hold variable `vref`.
An `ArgumentError` is thrown if `vref` is not a hold variable.
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

Specify more parameter bounds for a hold variable `vref`.
An `ArgumentError` is thrown if `vref` is not a hold variable.
"""
function add_parameter_bounds(vref::GeneralVariableRef,
                              bounds::ParameterBounds{GeneralVariableRef};
                              _error::Function = error
                              )::Nothing
    return add_parameter_bounds(dispatch_variable_ref(vref), bounds,
                                _error = _error)
end

# Dispatch fallback
function delete_parameter_bounds(vref::DispatchVariableRef)
    throw(ArgumentError("`delete_parameter_bounds` not defined for variable reference type(s) " *
                        "`$(typeof(vref))`."))
end

"""
    delete_parameter_bounds(vref::GeneralVariableRef)::Nothing

Delete all the parameter bounds of the hold variable `vref`. An `ArgumentError`
is thrown if `vref` is not a hold variable.
"""
function delete_parameter_bounds(vref::GeneralVariableRef)::Nothing
    return delete_parameter_bounds(dispatch_variable_ref(vref))
end

################################################################################
#                            LOWER BOUND METHODS
################################################################################
# Dispatch fallback
function JuMP.has_lower_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.has_lower_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.has_lower_bound(vref::GeneralVariableRef)::Bool

Extend [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::JuMP.VariableRef)) to
return `Bool` if `vref` has lower bound. It relies on `JuMP.has_lower_bound`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.has_lower_bound(vref::GeneralVariableRef)::Bool
    return JuMP.has_lower_bound(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.lower_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.lower_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.lower_bound(vref::GeneralVariableRef)::Float64

Extend [`JuMP.lower_bound`](@ref JuMP.lower_bound(::JuMP.VariableRef)) to
return the lower bound of `vref`. It relies on `JuMP.lower_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.lower_bound(vref::GeneralVariableRef)::Float64
    return JuMP.lower_bound(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_lower_bound(vref::DispatchVariableRef, value::Real)
    throw(ArgumentError("`JuMP.set_lower_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_lower_bound(vref::GeneralVariableRef, value::Real)::Nothing

Extend [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::JuMP.VariableRef, ::Number)) to
set the lower bound of `vref`. It relies on `JuMP.set_lower_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_lower_bound(vref::GeneralVariableRef, value::Real)::Nothing
    return JuMP.set_lower_bound(dispatch_variable_ref(vref), value)
end

# Dispatch fallback
function JuMP.LowerBoundRef(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.LowerBoundRef` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.LowerBoundRef(vref::GeneralVariableRef)::InfOptConstraintRef

Extend [`JuMP.LowerBoundRef`](@ref JuMP.LowerBoundRef(::JuMP.VariableRef)) to
return the lower bound constraint of `vref`. It relies on `JuMP.LowerBoundRef` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.LowerBoundRef(vref::GeneralVariableRef)::InfOptConstraintRef
    return JuMP.LowerBoundRef(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.delete_lower_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.delete_lower_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.delete_lower_bound(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.delete_lower_bound`](@ref JuMP.delete_lower_bound(::JuMP.VariableRef)) to
delete the lower bound of `vref`. It relies on `JuMP.delete_lower_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.delete_lower_bound(vref::GeneralVariableRef)::Nothing
    return JuMP.delete_lower_bound(dispatch_variable_ref(vref))
end

################################################################################
#                            UPPER BOUND METHODS
################################################################################
# Dispatch fallback
function JuMP.has_upper_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.has_upper_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.has_upper_bound(vref::GeneralVariableRef)::Bool

Extend [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::JuMP.VariableRef)) to
return `Bool` if `vref` has upper bound. It relies on `JuMP.has_upper_bound`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.has_upper_bound(vref::GeneralVariableRef)::Bool
    return JuMP.has_upper_bound(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.upper_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.upper_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.upper_bound(vref::GeneralVariableRef)::Float64

Extend [`JuMP.upper_bound`](@ref JuMP.upper_bound(::JuMP.VariableRef)) to
return the upper bound of `vref`. It relies on `JuMP.upper_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.upper_bound(vref::GeneralVariableRef)::Float64
    return JuMP.upper_bound(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_upper_bound(vref::DispatchVariableRef, value::Real)
    throw(ArgumentError("`JuMP.set_upper_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_upper_bound(vref::GeneralVariableRef, value::Real)::Nothing

Extend [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::JuMP.VariableRef, ::Number)) to
set the upper bound of `vref`. It relies on `JuMP.set_upper_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_upper_bound(vref::GeneralVariableRef, value::Real)::Nothing
    return JuMP.set_upper_bound(dispatch_variable_ref(vref), value)
end

# Dispatch fallback
function JuMP.UpperBoundRef(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.UpperBoundRef` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.UpperBoundRef(vref::GeneralVariableRef)::InfOptConstraintRef

Extend [`JuMP.UpperBoundRef`](@ref JuMP.UpperBoundRef(::JuMP.VariableRef)) to
return the upper bound constraint of `vref`. It relies on `JuMP.UpperBoundRef` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.UpperBoundRef(vref::GeneralVariableRef)::InfOptConstraintRef
    return JuMP.UpperBoundRef(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.delete_upper_bound(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.delete_upper_bound` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.delete_upper_bound(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.delete_upper_bound`](@ref JuMP.delete_upper_bound(::JuMP.VariableRef)) to
delete the upper bound of `vref`. It relies on `JuMP.delete_upper_bound` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.delete_upper_bound(vref::GeneralVariableRef)::Nothing
    return JuMP.delete_upper_bound(dispatch_variable_ref(vref))
end

################################################################################
#                                FIXING METHODS
################################################################################
# Dispatch fallback
function JuMP.is_fixed(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.is_fixed` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.is_fixed(vref::GeneralVariableRef)::Bool

Extend [`JuMP.is_fixed`](@ref JuMP.is_fixed(::JuMP.VariableRef)) to
return `Bool` if `vref` is fixed. It relies on `JuMP.is_fixed`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.is_fixed(vref::GeneralVariableRef)::Bool
    return JuMP.is_fixed(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.fix_value(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.fix_value` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.fix_value(vref::GeneralVariableRef)::Float64

Extend [`JuMP.fix_value`](@ref JuMP.fix_value(::JuMP.VariableRef)) to
return the fixed value of `vref`. It relies on `JuMP.fix_value` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.fix_value(vref::GeneralVariableRef)::Float64
    return JuMP.fix_value(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.fix(vref::DispatchVariableRef, value::Real; force::Bool = false)
    throw(ArgumentError("`JuMP.fix` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.fix(vref::GeneralVariableRef, value::Real; force::Bool = false)::Nothing

Extend [`JuMP.fix`](@ref JuMP.fix(::JuMP.VariableRef, ::Number)) to
fix the value of `vref`. It relies on `JuMP.fix` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.fix(vref::GeneralVariableRef, value::Real;
                  force::Bool = false)::Nothing
    return JuMP.fix(dispatch_variable_ref(vref), value, force = force)
end

# Dispatch fallback
function JuMP.FixRef(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.FixRef` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.FixRef(vref::GeneralVariableRef)::InfOptConstraintRef

Extend [`JuMP.FixRef`](@ref JuMP.FixRef(::JuMP.VariableRef)) to
return the fix bound constraint of `vref`. It relies on `JuMP.FixRef` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.FixRef(vref::GeneralVariableRef)::InfOptConstraintRef
    return JuMP.FixRef(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.unfix(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.unfix` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.unfix(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.unfix`](@ref JuMP.unfix(::JuMP.VariableRef)) to
delete the upper bound of `vref`. It relies on `JuMP.unfix` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.unfix(vref::GeneralVariableRef)::Nothing
    return JuMP.unfix(dispatch_variable_ref(vref))
end

################################################################################
#                              START VALUE METHODS
################################################################################
# Dispatch fallback
function JuMP.start_value(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.start_value` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.start_value(vref::GeneralVariableRef)::Union{Nothing, Float64}

Extend [`JuMP.start_value`](@ref JuMP.start_value(::JuMP.VariableRef)) to
return the start value of `vref`. It relies on `JuMP.start_value` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.start_value(vref::GeneralVariableRef)::Union{Nothing, Float64}
    return JuMP.start_value(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_start_value(vref::DispatchVariableRef, value::Real)
    throw(ArgumentError("`JuMP.set_start_value` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_start_value(vref::GeneralVariableRef, value::Real)::Nothing

Extend [`JuMP.set_start_value`](@ref JuMP.set_start_value(::JuMP.VariableRef, ::Number)) to
set the upper bound of `vref`. It relies on `JuMP.set_start_value` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_start_value(vref::GeneralVariableRef, value::Real)::Nothing
    return JuMP.set_start_value(dispatch_variable_ref(vref), value)
end

################################################################################
#                                BINARY METHODS
################################################################################
# Dispatch fallback
function JuMP.is_binary(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.is_binary` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.is_binary(vref::GeneralVariableRef)::Bool

Extend [`JuMP.is_binary`](@ref JuMP.is_binary(::JuMP.VariableRef)) to
return `Bool` if `vref` is binary. It relies on `JuMP.is_binary`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.is_binary(vref::GeneralVariableRef)::Bool
    return JuMP.is_binary(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_binary(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.set_binary` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_binary(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.set_binary`](@ref JuMP.set_binary(::JuMP.VariableRef)) to
set `vref` as binary. It relies on `JuMP.set_binary` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_binary(vref::GeneralVariableRef)::Nothing
    return JuMP.set_binary(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.BinaryRef(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.BinaryRef` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.BinaryRef(vref::GeneralVariableRef)::InfOptConstraintRef

Extend [`JuMP.BinaryRef`](@ref JuMP.BinaryRef(::JuMP.VariableRef)) to
return the binary constraint of `vref`. It relies on `JuMP.BinaryRef` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.BinaryRef(vref::GeneralVariableRef)::InfOptConstraintRef
    return JuMP.BinaryRef(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.unset_binary(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.unset_binary` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.unset_binary(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.unset_binary`](@ref JuMP.unset_binary(::JuMP.VariableRef)) to
unset `vref` as binary. It relies on `JuMP.unset_binary` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.unset_binary(vref::GeneralVariableRef)::Nothing
    return JuMP.unset_binary(dispatch_variable_ref(vref))
end

################################################################################
#                               INTEGER METHODS
################################################################################
# Dispatch fallback
function JuMP.is_integer(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.is_integer` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.is_integer(vref::GeneralVariableRef)::Bool

Extend [`JuMP.is_integer`](@ref JuMP.is_integer(::JuMP.VariableRef)) to
return `Bool` if `vref` is integer. It relies on `JuMP.is_integer`
being defined for the underlying `DispatchVariableRef`, otherwise an
`ArugmentError` is thrown.
"""
function JuMP.is_integer(vref::GeneralVariableRef)::Bool
    return JuMP.is_integer(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.set_integer(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.set_integer` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.set_integer(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.set_integer`](@ref JuMP.set_integer(::JuMP.VariableRef)) to
set `vref` as integer. It relies on `JuMP.set_integer` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.set_integer(vref::GeneralVariableRef)::Nothing
    return JuMP.set_integer(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.IntegerRef(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.IntegerRef` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.IntegerRef(vref::GeneralVariableRef)::InfOptConstraintRef

Extend [`JuMP.IntegerRef`](@ref JuMP.IntegerRef(::JuMP.VariableRef)) to
return the integer constraint of `vref`. It relies on `JuMP.IntegerRef` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.IntegerRef(vref::GeneralVariableRef)::InfOptConstraintRef
    return JuMP.IntegerRef(dispatch_variable_ref(vref))
end

# Dispatch fallback
function JuMP.unset_integer(vref::DispatchVariableRef)
    throw(ArgumentError("`JuMP.unset_integer` not defined for variable reference type " *
                        "`$(typeof(vref))`."))
end

"""
    JuMP.unset_integer(vref::GeneralVariableRef)::Nothing

Extend [`JuMP.unset_integer`](@ref JuMP.unset_integer(::JuMP.VariableRef)) to
unset `vref` as integer. It relies on `JuMP.unset_integer` being defined
for the underlying `DispatchVariableRef`, otherwise an `ArugmentError` is thrown.
"""
function JuMP.unset_integer(vref::GeneralVariableRef)::Nothing
    return JuMP.unset_integer(dispatch_variable_ref(vref))
end
