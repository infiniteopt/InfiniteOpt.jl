################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::DependentParameterIndex
                               )::DependentParameterRef
    return DependentParameterRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::MultiParameterData{<:InfiniteArraySet}
                          )::DependentParametersIndex
    index = MOIUC.add_item(model.dependent_params, object)
    push!(model.param_object_indices, index)
    return index
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel,
    ::Type{DependentParameters})::MOIUC.CleverDict
    return model.dependent_params
end

# Extend _data_dictionary (ref based)
function _data_dictionary(pref::DependentParameterRef)::MOIUC.CleverDict
    return JuMP.owner_model(pref).dependent_params
end

# Extend _data_object
function _data_object(pref::DependentParameterRef)::MultiParameterData
    object = get(_data_dictionary(pref), JuMP.index(pref).object_index, nothing)
    object === nothing && error("Invalid dependent parameter reference, cannot find " *
                           "corresponding parameter in the model. This is likely " *
                           "caused by using the reference of a deleted parameter.")
    return object
end

# Extend _core_variable_object
function _core_variable_object(pref::DependentParameterRef)::DependentParameters
    return _data_object(pref).parameters
end

# Return the number of dependent parameters involved
function _num_parameters(pref::DependentParameterRef)::Int
    return length(_data_object(pref).names)
end

# Extend _delete_data_object
function _delete_data_object(vref::DependentParameterRef)::Nothing
    delete!(_data_dictionary(vref), JuMP.index(vref).object_index)
    return
end

################################################################################
#                             PARAMETER DEFINITION
################################################################################
# Store partially processed individual dependent parameters
struct _DependentParameter{S <: AbstractInfiniteSet}
    set::S
    supports::Vector{Float64}
    name::String
    function _DependentParameter(set::S,
        supports::Union{Vector{<:Real}, Real},
        name::String
        )::_DependentParameter{S} where {S <: AbstractInfiniteSet}
        if supports isa Real
            return new{S}(set, [supports], name)
        else
            return new{S}(set, supports, name)
        end
    end
end

## Use type dispatch to efficiently check the set(s)
# All InfiniteScalarSets
function _check_param_sets(_error::Function,
    params::AbstractArray{<:_DependentParameter{<:InfiniteScalarSet}}
    )::Nothing
    return
end

# All MultiDistributionSet with non SparseAxisArray
function _check_param_sets(_error::Function,
    params::AbstractArray{<:_DependentParameter{<:MultiDistributionSet}}
    )::Nothing
    dist = first(params).set.distribution
    set_type = typeof(first(params).set)
    if size(dist) != size(params)
        _error("The dimensions of the parameters and the multi-dimensional " *
               "distribution $dist.")
    end
    return
end

# All MultiDistributionSet with SparseAxisArray
function _check_param_sets(_error::Function,
    params::JuMPC.SparseAxisArray{<:_DependentParameter{<:MultiDistributionSet}}
    )
    _error("Cannot specify multiple-dimensional distribution set with a " *
           "`SparseAxisArray` of dependent infinite parameters.")
end

# All CollectionSet
function _check_param_sets(_error::Function,
    params::AbstractArray{<:_DependentParameter{<:CollectionSet}}
    )::Nothing
    sets = collection_sets(first(params).set)
    set_type = typeof(first(params).set)
    if length(sets) != length(params)
        _error("The dimensions of the parameters and the specified CollectionSet " *
               "do not match.")
    elseif params isa JuMPC.SparseAxisArray
        @warn("CollectionSet order may not match the given `SparseAxisArray` " *
              "of specified dependent infinite parameters, consider instead " *
              "specifying the `InfiniteScalarSet` for each parameter using " *
              "the `set` keyword and the appropriate indices.")
    end
    return
end

# All some InfiniteArraySet (for extensions)
function _check_param_sets(_error::Function,
    params::AbstractArray{<:_DependentParameter{<:InfiniteArraySet}}
    )::Nothing
    return
end

# Mixed sets
function _check_param_sets(_error::Function,
    params::AbstractArray{<:_DependentParameter}
    )::Nothing
    set_type = typeof(first(params).set)
    if !all(param.set isa InfiniteScalarSet for param in params)
        _error("Cannot specify multiple `InfiniteScalarSets` for one container " *
               "of infinite dependent parameters.")
    end
    return
end

# Fallback
function _check_param_sets(_error::Function, params)
    _error("Unrecognized infinite set input.")
end

## Use set type dispatch to make the proper InfiniteArraySet
# InfiniteArraySet
function _make_array_set(params::Vector{<:_DependentParameter{T}}
                         )::T where {T <: InfiniteArraySet}
    return first(params).set
end

# InfiniteScalarSets
function _make_array_set(params::Vector{<:_DependentParameter{T}}
                         )::CollectionSet{T} where {T <: InfiniteScalarSet}
    return CollectionSet([p.set for p in params])
end

# Build a DependentParameters object given an array of _DependentParameters
function _build_parameters(_error::Function,
                           params::AbstractArray{<:_DependentParameter};
                           num_supports::Int = 0,
                           sig_digits::Int = DefaultSigDigits,
                           extra_kw_args...)
    # error with extra keywords
    for (kwarg, _) in extra_kw_args
       _error("Unrecognized keyword argument $kwarg")
    end
    # check the formatting
    _check_param_sets(_error, params)
    # vectorize the parameter array
    indices = Collections._get_indices(params)
    ordered_params = Collections._make_ordered(params, indices)
    vector_params = _make_vector(ordered_params)
    # make the set
    set = _make_array_set(vector_params)
    # make the supports and labels
    lens = [length(p.supports) for p in vector_params]
    _allequal(lens) || _error("Inconsistent support dimensions.")
    # we have supports
    if first(lens) != 0
        # build the support array transpose to fill in column order (leverage locality)
        trans_supps = Array{Float64}(undef, first(lens), length(vector_params))
        for i in 1:size(trans_supps, 2)
            trans_supps[:, i] = vector_params[i].supports
        end
        supps = permutedims(trans_supps)
        supports_in_set(supps, set) || _error("Supports violate infinite set domain.")
        supps = round.(supps, sigdigits = sig_digits)
        label = UserDefined
        supp_dict = Dict{Vector{Float64}, Set{Symbol}}(@views supps[:, i] =>
                                         Set([label]) for i in 1:size(supps, 2))
    # we want to generate supports
    elseif num_supports != 0
        supps, label = generate_support_values(set, num_supports = num_supports,
                                               sig_digits = sig_digits)
        supp_dict = Dict{Vector{Float64}, Set{Symbol}}(@views supps[:, i] =>
                                         Set([label]) for i in 1:size(supps, 2))
    # no supports are specified
    else
        supp_dict = Dict{Vector{Float64}, Set{Symbol}}()
    end
    # make the parameter object
    names = [param.name for param in vector_params]
    return DependentParameters(set, supp_dict, sig_digits), names, indices
end

"""
    add_parameters(model::InfiniteModel,
                   params::DependentParameters,
                   names::Vector{String} = ["noname", "noname", ...],
                   indices = nothing
                   )::AbstractArray{<:GeneralVariableRef}

Add `params` to `model` and return an appropriate container of the dependent
infinite parameter references. This is intended as an internal method for use
with [`@dependent_parameters`](@ref). However, if desired users can use this
to create a container of infinite dependent parameter without the use of a
macro. `names` denote the name of each parameter and `indices` denote the
indices of the expected container as used by `Containers._make_array`
(implemented by `VectorTuple`s), by default a `Vector` is returned.

**Example**
```julia-repl
julia> using Distributions

julia> dist = MvNormal(ones(3)); # 3 dimensional

julia> set = MultiDistributionSet(dist); # 3 dimensional

julia> params = DependentParameters(set, Dict{Vector{Float64}, Set{Symbol}}(), 10);

julia> prefs = add_parameters(model, params, ["par1", "par2", "par3"])
3-element Array{GeneralVariableRef,1}:
 par1
 par2
 par3
```
"""
function add_parameters(model::InfiniteModel,
                        params::DependentParameters,
                        names::Vector{String} = String[],
                        indices = nothing
                        )::AbstractArray{<:GeneralVariableRef}
    # get the number of parameters
    num_params = length(params.set)
    # process the names
    if isempty(names)
        names = ["noname" for i in 1:num_params]
    end
    # process the indices
    if indices === nothing
        indices = CartesianIndices(1:num_params)
    end
    # make the parameter model object
    obj_num = length(_param_object_indices(model)) + 1
    first_param_num = model.last_param_num + 1
    last_param_num = model.last_param_num += num_params
    param_nums = first_param_num:last_param_num
    data_object = MultiParameterData(params, obj_num, param_nums, names)
    # add the data object to the model and make the references
    obj_index = _add_data_object(model, data_object)
    prefs = [GeneralVariableRef(model, obj_index.value, DependentParameterIndex, i)
             for i in 1:num_params]
    return Collections._make_array(prefs, indices)
end

# Construct an expression to build an infinite set(s) (use with @dependent_parameters)
function _construct_array_set(_error::Function, info::_ParameterInfoExpr)
    if (info.has_lb || info.has_ub) && !(info.has_lb && info.has_ub)
        _error("Must specify both an upper bound and a lower bound")
    elseif info.has_lb
        check = :(isa($(info.lower_bound), Real))
        return :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a real number."))
    elseif info.has_dist
        check = :(isa($(info.distribution), Distributions.UnivariateDistribution))
        return :($(check) ? UniDistributionSet($(info.distribution)) : MultiDistributionSet($(info.distribution)))
    elseif info.has_set
        check1 = :(isa($(info.set), AbstractInfiniteSet))
        check2 = :(isa($(info.set), Distributions.UnivariateDistribution))
        return :($(check1) ? $(info.set) : ($(check2) ? UniDistributionSet($(info.set)) : MultiDistributionSet($(info.set))))
    else
        _error("Must specify upper/lower bounds, a distribution, or a set")
    end
end

################################################################################
#                                  NAMING
################################################################################
# Get the parameter index in the DependentParameters object
_param_index(pref::DependentParameterRef)::Int = JuMP.index(pref).param_index

"""
    JuMP.name(pref::DependentParameterRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to return the names of
infinite dependent parameters.

**Example**
```julia-repl
julia> name(pref)
"par_name"
```
"""
function JuMP.name(pref::DependentParameterRef)::String
    object = get(_data_dictionary(pref), JuMP.index(pref).object_index, nothing)
    return object === nothing ? "" : object.names[_param_index(pref)]
end

"""
    JuMP.set_name(pref::DependentParameterRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of dependent infinite parameters.

**Example**
```julia-repl
julia> set_name(vref, "par_name")

julia> name(vref)
"para_name"
```
"""
function JuMP.set_name(pref::DependentParameterRef, name::String)::Nothing
    _data_object(pref).names[_param_index(pref)] = name
    JuMP.owner_model(pref).name_to_param = nothing
    return
end

################################################################################
#                           PARAMETER DEPENDENCIES
################################################################################
# Extend _infinite_variable_dependencies
function _infinite_variable_dependencies(pref::DependentParameterRef
    )::Vector{InfiniteVariableIndex}
    return _data_object(pref).infinite_var_indices
end

# Extend _measure_dependencies
function _measure_dependencies(pref::DependentParameterRef
    )::Vector{MeasureIndex}
    return _data_object(pref).measure_indices[_param_index(pref)]
end

# Extend _constraint_dependencies
function _constraint_dependencies(pref::DependentParameterRef
    )::Vector{ConstraintIndex}
    return _data_object(pref).constraint_indices[_param_index(pref)]
end

# Extend _derivative_dependencies
function _derivative_dependencies(pref::DependentParameterRef
    )::Vector{DerivativeIndex}
    return _data_object(pref).derivative_indices[_param_index(pref)]
end

"""
    used_by_infinite_variable(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
an infinite variable.

**Example**
```julia-repl
julia> used_by_objective(pref)
true
```
"""
function used_by_infinite_variable(pref::DependentParameterRef)::Bool
    return !isempty(_infinite_variable_dependencies(pref))
end

"""
    used_by_measure(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used
by a measure.

**Example**
```julia-repl
julia> used_by_measure(pref)
true
```
"""
function used_by_measure(pref::DependentParameterRef)::Bool
    return !isempty(_measure_dependencies(pref))
end

"""
    used_by_constraint(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
a constraint.

**Example**
```julia-repl
julia> used_by_constraint(pref)
false
```
"""
function used_by_constraint(pref::DependentParameterRef)::Bool
    return !isempty(_constraint_dependencies(pref))
end

"""
    used_by_derivative(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
a derivative.

**Example**
```julia-repl
julia> used_by_derivative(pref)
false
```
"""
function used_by_derivative(pref::DependentParameterRef)::Bool
    return !isempty(_derivative_dependencies(pref))
end

# Extend used by objective
used_by_objective(pref::DependentParameterRef)::Bool = false

"""
    is_used(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used in
the model.

**Example**
```julia-repl
julia> is_used(pref)
true
```
"""
function is_used(pref::DependentParameterRef)::Bool
    return used_by_measure(pref) || used_by_constraint(pref) ||
           used_by_infinite_variable(pref) || used_by_derivative(pref)
end

################################################################################
#                          PARAMETER OBJECT METHODS
################################################################################
# Extend _parameter_number
function _parameter_number(pref::DependentParameterRef)::Int
    return _data_object(pref).parameter_nums[_param_index(pref)]
end

# Extend _parameter_numbers
function _parameter_numbers(pref::DependentParameterRef)::Vector{Int}
    return [_parameter_number(pref)]
end

# Extend _object_number
function _object_number(pref::DependentParameterRef)::Int
    return _data_object(pref).object_num
end

# Extend _object_numbers
function _object_numbers(pref::DependentParameterRef)::Vector{Int}
    return [_object_number(pref)]
end

# Extend _set_core_variable_object
function _set_core_variable_object(pref::DependentParameterRef,
                                   params::DependentParameters)::Nothing
    _data_object(pref).parameters = params
    return
end

################################################################################
#                             INFINITE SET METHODS
################################################################################
## Get the individual infinite set if possible
# raw_set
function _parameter_set(pref::DependentParameterRef)::InfiniteArraySet
    return _core_variable_object(pref).set
end

# CollectionSet
function _parameter_set(set::CollectionSet{S},
                        pref::DependentParameterRef
                        )::S where {S <: InfiniteScalarSet}
    return collection_sets(set)[_param_index(pref)]
end

# InfiniteArraySet (Fallback)
function _parameter_set(set::InfiniteArraySet, pref::DependentParameterRef)
    error("An individual infinite set is not well-defined for $pref which " *
          "is part of a group of dependent infinite parameters that correspond " *
          "to an multi-dimensional infinite set of type `$(typeof(set))`.")
end

"""
    infinite_set(pref::DependentParameterRef)::InfiniteScalarSet

Return the infinite set associated with the particular infinite dependent
parameter `pref` if valid. Errors if the underlying [`DependentParameters`](@ref)
object does not use a [`CollectionSet`](@ref).

**Example**
```julia-repl
julia> infinite_set(x[1])
[-1, 1]
```
"""
function infinite_set(pref::DependentParameterRef)::InfiniteScalarSet
    return _parameter_set(_parameter_set(pref), pref)
end

# Check that prefs are complete
function _check_complete_param_array(
    prefs::AbstractArray{<:DependentParameterRef}
    )::Nothing
    if length(prefs) != _num_parameters(first(prefs))
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    end
    return
end

"""
    infinite_set(prefs::AbstractArray{<:DependentParameterRef})::InfiniteArraySet

Return the infinite set associated with the container of infinite dependent
parameters `prefs`. Errors if the container `prefs` is incomplete.

**Example**
```julia-repl
julia> infinite_set(x)
ZeroMeanDiagNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 1.0]
)
```
"""
function infinite_set(prefs::AbstractArray{<:DependentParameterRef}
                      )::InfiniteArraySet
    _check_complete_param_array(prefs)
    return _parameter_set(first(prefs))
end

# Update the underlying set and delete the supports
function _update_parameter_set(pref::DependentParameterRef,
                               new_set::InfiniteArraySet)::Nothing
    old_params = _core_variable_object(pref)
    new_supports = Dict{Vector{Float64}, Set{Symbol}}()
    new_params = DependentParameters(new_set, new_supports,
                                     significant_digits(pref))
    _set_core_variable_object(pref, new_params)
    if is_used(pref)
        set_optimizer_model_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    set_infinite_set(pref::DependentParameterRef,
                     set::InfiniteScalarSet)::Nothing

Specify the scalar infinite set of the dependent infinite parameter `pref` to
`set` if `pref` is part of a [`CollectionSet`](@ref), otherwise an error is
thrown. Note this will reset/delete all the supports contained in the
underlying [`DependentParameters`](@ref) object. Also, errors if `pref` is used
by a measure or if the new set type is different than the original.

**Example**
```julia-repl
julia> set_infinite_set(x[1], IntervalSet(0, 2))

julia> infinite_set(x[1])
[0, 2]
```
"""
function set_infinite_set(pref::DependentParameterRef,
                          set::InfiniteScalarSet)::Nothing
    old_set = _parameter_set(pref)
    if !(old_set isa CollectionSet)
        error("Cannot set the individual infinite set of $pref if the " *
              "underlying set is not a CollectionSet.")
    elseif used_by_measure(pref)
        error("Cannot override the infinite set of $pref since it is used by " *
              "a measure.")
    elseif !(CollectionSet{typeof(set)} <: typeof(old_set))
        error("Cannot change the underlying set type.")
    end
    param_idx = _param_index(pref)
    new_set = CollectionSet([i != param_idx ? collection_sets(old_set)[i] : set
                             for i in eachindex(collection_sets(old_set))])
    _update_parameter_set(pref, new_set)
    return
end

"""
    set_infinite_set(prefs::AbstractArray{<:DependentParameterRef},
                     set::InfiniteArraySet)::Nothing

Specify the multi-dimensional infinite set of the dependent infinite parameters
`prefs` to `set`. Note this will reset/delete all the supports contained in the
underlying [`DependentParameters`](@ref) object. This will error if the not all
of the dependent infinite parameters are included, if any of them are used by
measures, or if the new set type is different than the previous.

**Example**
```julia-repl
julia> set_infinite_set(x, CollectionSet([IntervalSet(0, 1), IntervalSet(0, 2)]))
```
"""
function set_infinite_set(prefs::AbstractArray{<:DependentParameterRef},
                          set::InfiniteArraySet)::Nothing
    if any(used_by_measure(pref) for pref in prefs)
        error("Cannot override the infinite set of $prefs since it is used by " *
              "a measure.")
    elseif !(typeof(set) <: typeof(infinite_set(prefs))) # this checks prefs
        error("Cannot change the underlying set type.")
    end
    _update_parameter_set(first(prefs), set)
    return
end

"""
    JuMP.has_lower_bound(pref::DependentParameterRef)::Bool

Extend the `JuMP.has_lower_bound` function to accomodate a single dependent
infinite parameter.
Return true if the set associated with `pref` has a defined lower bound or if a
lower bound can be found. Extensions with user-defined scalar infinite set types
should extend `JuMP.has_lower_bound(set::NewType)`.

**Example**
```julia-repl
julia> has_lower_bound(x[1])
true
```
"""
function JuMP.has_lower_bound(pref::DependentParameterRef)::Bool
    set = _parameter_set(pref)
    if set isa CollectionSet
        return JuMP.has_lower_bound(collection_sets(set)[_param_index(pref)])
    else
        return false
    end
end

"""
    JuMP.lower_bound(pref::DependentParameterRef)::Number

Extend the `JuMP.lower_bound` function to accomodate a single dependent infinite
parameter. Returns the lower bound associated with the infinite set. Errors if
such a bound is not well-defined.

**Example**
```julia-repl
julia> lower_bound(x[1])
0.0
```
"""
function JuMP.lower_bound(pref::DependentParameterRef)::Number
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    return JuMP.lower_bound(infinite_set(pref))
end

"""
    JuMP.set_lower_bound(pref::DependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_lower_bound` function to accomodate a single dependent
infinite parameter. Updates the infinite set lower bound if such an operation
is supported. Infinite scalar set extensions that seek to employ this should extend
`JuMP.set_lower_bound(set::NewType, lower::Number)`. This will call
[`set_infinite_set`](@ref) and will error if this is not well-defined. Note
that existing supports will be deleted.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_lower_bound(t, -1)

julia> lower_bound(t)
-1.0
```
"""
function JuMP.set_lower_bound(pref::DependentParameterRef, lower::Real)::Nothing
    set = infinite_set(pref)
    new_set = JuMP.set_lower_bound(set, lower)
    set_infinite_set(pref, new_set)
    return
end

"""
    JuMP.has_upper_bound(pref::DependentParameterRef)::Bool

Extend the `JuMP.has_upper_bound` function to accomodate a single dependent
infinite parameter.
Return true if the set associated with `pref` has a defined upper bound or if a
upper bound can be found. Extensions with user-defined scalar infinite set types
should extend `JuMP.has_upper_bound(set::NewType)`.

**Example**
```julia-repl
julia> has_upper_bound(x[1])
true
```
"""
function JuMP.has_upper_bound(pref::DependentParameterRef)::Bool
    set = _core_variable_object(pref).set
    if set isa CollectionSet
        return JuMP.has_upper_bound(collection_sets(set)[_param_index(pref)])
    else
        return false
    end
end

"""
    JuMP.upper_bound(pref::DependentParameterRef)::Number

Extend the `JuMP.upper_bound` function to accomodate a single dependent infinite
parameter. Returns the upper bound associated with the infinite set. Errors if
such a bound is not well-defined.

**Example**
```julia-repl
julia> upper_bound(x[1])
0.0
```
"""
function JuMP.upper_bound(pref::DependentParameterRef)::Number
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    return JuMP.upper_bound(infinite_set(pref))
end

"""
    JuMP.set_upper_bound(pref::DependentParameterRef, upper::Real)::Nothing

Extend the `JuMP.set_upper_bound` function to accomodate a single dependent
infinite parameter. Updates the infinite set upper bound if such an operation
is supported. Infinite scalar set extensions that seek to employ this should extend
`JuMP.set_upper_bound(set::NewType, upper::Number)`. This will call
[`set_infinite_set`](@ref) and will error if this is not well-defined. Note
that existing supports will be deleted.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1]))
julia> set_upper_bound(t, -1)

julia> upper_bound(t)
-1.0
```
"""
function JuMP.set_upper_bound(pref::DependentParameterRef, upper::Real)::Nothing
    set = infinite_set(pref)
    new_set = JuMP.set_upper_bound(set, upper)
    set_infinite_set(pref, new_set)
    return
end

################################################################################
#                              SUPPORT METHODS
################################################################################
# Get the raw supports
function _parameter_supports(pref::DependentParameterRef
                             )::Dict{Vector{Float64}, Set{Symbol}}
    return _core_variable_object(pref).supports
end

"""
    significant_digits(pref::DependentParameterRef)::Int

Return the number of significant digits enforced on the supports of `pref`.

**Example**
```julia-repl
julia> significant_digits(x[1])
12
```
"""
function significant_digits(pref::DependentParameterRef)::Int
    return _core_variable_object(pref).sig_digits
end

"""
    num_supports(pref::DependentParameterRef; [label::Symbol = All])::Int

Return the number of support points associated with a single dependent infinite
parameter `pref`. Specify a subset of supports via `label` to only count the
supports with `label`.

**Example**
```julia-repl
julia> num_supports(x[1])
2

julia> num_supports(x[1], label = MCSample)
0
```
"""
function num_supports(pref::DependentParameterRef; label::Symbol = All)::Int
    supp_dict = _parameter_supports(pref)
    if label == All
        return length(supp_dict)
    else
        return count(p -> label in p[2], supp_dict)
    end
end

"""
    num_supports(prefs::AbstractArray{<:DependentParameterRef};
                 [label::Symbol = All])::Int

Return the number of support points associated with dependent infinite
parameters `prefs`. Errors if not all from the same underlying object.
Specify a subset of supports via `label` to only count the supports with `label`.

**Example**
```julia-repl
julia> num_supports(x)
2
```
"""
function num_supports(prefs::AbstractArray{<:DependentParameterRef};
                      label::Symbol = All)::Int
    _check_complete_param_array(prefs)
    return num_supports(first(prefs), label = label)
end

"""
    has_supports(pref::DependentParameterRef)::Bool

Return true if `pref` has supports or false otherwise.

**Example**
```julia-repl
julia> has_supports(x[1])
true
```
"""
has_supports(pref::DependentParameterRef)::Bool = !isempty(_parameter_supports(pref))

"""
    has_supports(prefs::AbstractArray{<:DependentParameterRef})::Bool

Return true if `prefs` have supports or false otherwise. Errors if not all of the
infinite dependent parameters are from the same object.

**Example**
```julia-repl
julia> has_supports(x)
true
```
"""
function has_supports(prefs::AbstractArray{<:DependentParameterRef})::Bool
    _check_complete_param_array(prefs)
    return has_supports(first(prefs))
end

"""
    supports(pref::DependentParameterRef; [label::Symbol = All])::Vector{Float64}

Return the support points associated with `pref`. A subset of supports can be
returned via `label` to return just the supports associated with `label`.

**Example**
```julia-repl
julia> supports(x[1])
2-element Array{Float64,1}:
 0.0
 1.0
```
"""
function supports(pref::DependentParameterRef;
                  label::Symbol = All)::Vector{Float64}
    if label == All
        pindex = _param_index(pref)
        return Float64[supp[pindex] for supp in keys(_parameter_supports(pref))]
    else
        reduced_supps = findall(e -> label in e, _parameter_supports(pref))
        pindex = _param_index(pref)
        return Float64[supp[pindex] for supp in reduced_supps]
    end
end

"""
    supports(prefs::AbstractArray{<:DependentParameterRef};
             [label::Symbol = All]
             )::Union{AbstractArray{<:Vector{<:Float64}}, Array{Float64, 2}}

Return the support points associated with `prefs`. Errors if not all of the
infinite dependent parameters are from the same object. This will return a
matrix if `prefs` is `Vector`, otherwise an array of vectors is returned by
calling `supports.(prefs)`. A subset of supports can be
returned via `label` to return just the supports associated with `label`.

**Example**
```julia-repl
julia> supports(x) # columns are supports
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0
```
"""
function supports(prefs::AbstractArray{<:DependentParameterRef};
                  label::Symbol = All)::AbstractArray{<:Vector{<:Float64}}
    _check_complete_param_array(prefs)
    return supports.(prefs, label = label) # TODO make more efficient
end

# More efficient dispatch for Vectors
function supports(prefs::Vector{DependentParameterRef};
                  label::Symbol = All)::Array{Float64, 2}
    if !has_supports(prefs)
        return zeros(Float64, _num_parameters(first(prefs)), 0)
    elseif label == All
        raw_supps = keys(_parameter_supports(first(prefs)))
        if length(raw_supps) == 1
            return reduce(hcat, collect(raw_supps))
        else
            return reduce(hcat, raw_supps)
        end
    else
        raw_supps = findall(e -> label in e, _parameter_supports(first(prefs)))
        if isempty(raw_supps)
            return zeros(Float64, _num_parameters(first(prefs)), 0)
        else
            return reduce(hcat, raw_supps)
        end
    end
end

# Define method for overriding the current supports
function _update_parameter_supports(prefs::AbstractArray{<:DependentParameterRef},
                                    supports::Array{<:Real, 2},
                                    label::Symbol)::Nothing
    set = _parameter_set(first(prefs))
    new_supps = Dict{Vector{Float64}, Set{Symbol}}(@views supports[:, i] =>
                                      Set([label]) for i in 1:size(supports, 2))
    new_params = DependentParameters(set, new_supps, significant_digits(first(prefs)))
    _set_core_variable_object(first(prefs), new_params)
    if any(is_used(pref) for pref in prefs)
        set_optimizer_model_ready(JuMP.owner_model(first(prefs)), false)
    end
    return
end

# Process an array of vectors into a support matrix
function _make_support_matrix(prefs::AbstractArray{<:DependentParameterRef},
                              supports::AbstractArray{<:Vector{<:Real}}
                              )::Array{<:Real, 2}
    _keys(supports) == _keys(prefs) || error("Inconsistent support indices")
    lens = [length(supp) for supp in supports]
    _allequal(lens) || error("Inconsistent support dimensions.")
    trans_supps = Array{Float64}(undef, first(lens), length(prefs))
    for k in eachindex(prefs)
        trans_supps[:, _param_index(prefs[k])] = supports[k]
    end
    return permutedims(trans_supps)
end

"""
    set_supports(prefs::AbstractArray{<:DependentParameterRef},
                 supports::AbstractArray{<:Vector{<:Real}};
                 [force::Bool = false])::Nothing

Specify the support points for `prefs`. Errors if the supports violate the domain
of the infinite set, if the dimensions don't match up properly,
if `prefs` and `supports` have different indices, not all of the `prefs` are
from the same dependent infinite parameter container, there are existing
supports and `force = false`. Note that it is strongly preferred to use
`add_supports` if possible to avoid destroying measure dependencies.

```julia
    set_supports(prefs::Vector{DependentParameterRef},
                 supports::Array{<:Real, 2};
                 [force::Bool = false])::Nothing
```
Specify the supports for a vector `prefs` of dependent infinite parameters.
Here rows of `supports` correspond to `prefs` and the columns correspond to the
supports. This is more efficient than the above method and will error for the
same reasons.

**Example**
```julia-repl
julia> set_supports(y, [[0, 1], [0, 1]])

julia> set_supports(x, [0 1; 0 1])

julia> supports(x)
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0
```
"""
function set_supports(prefs::AbstractArray{<:DependentParameterRef},
                      supports::AbstractArray{<:Vector{<:Real}};
                      force::Bool = false,
                      label::Symbol = UserDefined
                      )::Nothing
    supps = _make_support_matrix(prefs, supports)
    set_supports(_make_vector(prefs), supps, force = force, label = label)
    return
end

# Efficient method for vector prefs and matrix of supports
function set_supports(prefs::Vector{DependentParameterRef},
                      supports::Array{<:Real, 2};
                      force::Bool = false,
                      label::Symbol = UserDefined
                      )::Nothing
    set = infinite_set(prefs) # this does a check on prefs
    if has_supports(prefs) && !force
        error("Unable set supports for $prefs since they already have supports." *
              " Consider using `add_supports` or use set `force = true` to " *
              "overwrite the existing supports.")
    elseif !supports_in_set(supports, set)
        error("Supports violate the domain of the infinite set.")
    end
    supports = round.(supports, sigdigits = significant_digits(first(prefs)))
    _update_parameter_supports(prefs, supports, label)
    return
end

# Error for single dependent parameters
function set_supports(pref::DependentParameterRef, supports; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    add_supports(prefs::AbstractArray{<:DependentParameterRef},
                 supports::AbstractArray{<:Vector{<:Real}})::Nothing

Add additional support points for `prefs`. Errors if the supports violate the domain
of the infinite set, if the dimensions don't match up properly,
if `prefs` and `supports` have different indices, or not all of the `prefs` are
from the same dependent infinite parameter container.

```julia
    add_supports(prefs::Vector{DependentParameterRef},
                 supports::Array{<:Real, 2})::Nothing
```
Specify the supports for a vector `prefs` of dependent infinite parameters.
Here rows of `supports` correspond to `prefs` and the columns correspond to the
supports. This is more efficient than the above method and will error for the
same reasons.

**Example**
```julia-repl
julia> add_supports(x, [[1], [1]])

julia> supports(x)
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0

julia> add_supports(x, ones(2, 1) * 0.5)

julia> supports(t)
2×3 Array{Float64,2}:
 0.0  1.0  0.5
 0.0  1.0  0.5
```
"""
function add_supports(prefs::AbstractArray{<:DependentParameterRef},
                      supports::AbstractArray{<:Vector{<:Real}};
                      label::Symbol = UserDefined, # interal keyword args
                      check::Bool = true)::Nothing
    supps = _make_support_matrix(prefs, supports)
    add_supports(_make_vector(prefs), supps, label = label, check = check)
    return
end

# More efficient version for supports in the correct format
function add_supports(prefs::Vector{DependentParameterRef},
                      supports::Array{<:Real, 2};
                      label::Symbol = UserDefined, # internal keyword args
                      check::Bool = true)::Nothing
    set = infinite_set(prefs) # this does a check on prefs
    if check && !supports_in_set(supports, set)
        error("Supports violate the domain of the infinite set.")
    end
    supports = round.(supports, sigdigits = significant_digits(first(prefs)))
    current_supports = _parameter_supports(first(prefs))
    for i in 1:size(supports, 2)
        s = @view(supports[:, i])
        if haskey(current_supports, s)
            push!(current_supports[s], label)
        else
            current_supports[s] = Set([label])
        end
    end
    if any(is_used(pref) for pref in prefs)
        set_optimizer_model_ready(JuMP.owner_model(first(prefs)), false)
    end
    return
end

# Error for single dependent parameters
function add_supports(pref::DependentParameterRef, supports; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    delete_supports(prefs::AbstractArray{<:DependentParameterRef})::Nothing

Delete the support points for `prefs`. Errors if any of the parameters are
used by a measure or if not all belong to the same set of dependent parameters.

**Example**
```julia-repl
julia> delete_supports(w)

```
"""
function delete_supports(prefs::AbstractArray{<:DependentParameterRef})::Nothing
    _check_complete_param_array(prefs)
    if any(used_by_measure(pref) for pref in prefs)
        error("Cannot delete supports with measure dependencies.")
    end
    empty!(_parameter_supports(first(prefs)))
    if any(is_used(pref) for pref in prefs)
        set_optimizer_model_ready(JuMP.owner_model(first(prefs)), false)
    end
    return
end

# Error for single dependent parameters
function delete_supports(pref::DependentParameterRef)
    error("Cannot delete the supports of a single dependent infinite parameter.")
end

# TODO resolve case that there are existing UniformGrid supports
"""
    generate_and_add_supports!(prefs::AbstractArray{<:DependentParameterRef},
                               set::InfiniteArraySet,
                               method::Union{Symbol, Nothing} = nothing;
                               [num_supports::Int = DefaultNumSupports])::Nothing

Generate supports for `prefs` via [`generate_support_values`](@ref) and add them
to `pref`. This is intended as an extendable internal method for
[`fill_in_supports!`](@ref fill_in_supports!(::AbstractArray{<:DependentParameterRef})).
Most extensions that employ user-defined infinite sets can typically enable this
by extending [`generate_support_values`](@ref). However, in some cases it may be
necessary to extend this when more complex operations need to take place then just
adding supports to a set of infinite parameters. Errors if the
infinite set type is not recognized.
"""
function generate_and_add_supports!(prefs::AbstractArray{<:DependentParameterRef},
                                    set::InfiniteArraySet,
                                    method::Union{Symbol, Nothing} = nothing;
                                    num_supports::Int = DefaultNumSupports)::Nothing
    new_supps, label = generate_supports(set, method,
                                         num_supports = num_supports,
                                    sig_digits = significant_digits(first(prefs)))
    add_supports(_make_vector(prefs), new_supps, check = false, label = label)
    return
end

"""
    fill_in_supports!(prefs::AbstractArray{<:DependentParameterRef};
                      [num_supports::Int = DefaultNumSupports,
                       modify::Bool = true])::Nothing

Automatically generate support points for a container of dependent infinite
parameters `prefs`. Generating up to `num_supports` for the parameters in accordance
with `generate_and_add_supports!`. Will add nothing if there are supports and
`modify = false`. Extensions that use user defined
set types should extend [`generate_and_add_supports!`](@ref) and/or
[`generate_support_values`](@ref) as needed. Errors if the infinite set type is
not recognized.

**Example**
```julia-repl
julia> fill_in_supports!(x, num_supports = 4)

julia> supports(x)
2×4 Array{Float64,2}:
 0.0  0.333  0.667  1.0
 0.0  0.333  0.667  1.0
```
"""
function fill_in_supports!(prefs::AbstractArray{<:DependentParameterRef};
                           num_supports::Int = DefaultNumSupports,
                           modify::Bool = true)::Nothing
    set = infinite_set(prefs) # does check for bad container
    current_amount = InfiniteOpt.num_supports(first(prefs))
    if (modify || current_amount == 0) && (current_amount < num_supports)
        generate_and_add_supports!(prefs, set,
                                   num_supports = num_supports - current_amount)
    end
    return
end

# Error for single dependent parameters
function fill_in_supports!(pref::DependentParameterRef; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    fill_in_supports!(model::InfiniteModel; [num_supports::Int = DefaultNumSupports,
                      modify::Bool = true])::Nothing

Automatically generate support points for all infinite parameters in model.
This calls `fill_in_supports!` for each parameter in the model.
See [`fill_in_supports!`](@ref)
for more information. Errors if one of the infinite set types is unrecognized.
Note that no supports will be added to a particular parameter if it already has
some and `modify = false`.

**Example**
```julia-repl
julia> fill_in_supports!(model, num_supports = 4)

julia> supports(t)
4-element Array{Float64,1}:
 0.0
 0.333
 0.667
 1.0
```
"""
function fill_in_supports!(model::InfiniteModel; num_supports::Int = DefaultNumSupports,
                           modify::Bool = true)::Nothing
    # fill in the the supports of each independent parameter
    for (key, data_object) in model.independent_params
        pref = dispatch_variable_ref(model, key)
        fill_in_supports!(pref, num_supports = num_supports, modify = modify)
    end
    # fill in the supports of each dependent parameter set
    for (key, data_object) in model.dependent_params
        prefs = [dispatch_variable_ref(model, DependentParameterIndex(key, i))
                 for i in 1:length(data_object.names)]
        fill_in_supports!(prefs, num_supports = num_supports, modify = modify)
    end
    return
end

################################################################################
#                          MODEL PARAMETER QUERIES
################################################################################
"""
    num_parameters(model::InfiniteModel,
                   [type::Type{InfOptParameter} = InfOptParameter])::Int

Return the number of `InfiniteOpt` parameters assigned to `model`. By default,
the total number of infinite and finite parameters is returned. The amount
of a particular type is obtained by specifying the concrete parameter type
of [`InfOptParameter`](@ref) via `type`. Type options include:
 - `InfOptParameter`: all parameters
 - `ScalarParameter`: all scalar parameters
 - `InfiniteParameter`: all infinite parameters
 - `FiniteParameter`: all finite parameters
 - `IndependentParameter`: all independent infinite parameters
 - `DependentParameters`: all dependent infinite parameters

**Example**
```julia-repl
julia> num_parameters(model)
3

julia> num_parameters(model, IndependentParameter)
2
```
"""
function num_parameters(model::InfiniteModel,
                        type::Type{InfOptParameter} = InfOptParameter
                        )::Int
    num_pars = num_parameters(model, IndependentParameter)
    num_pars += num_parameters(model, FiniteParameter)
    num_pars += num_parameters(model, DependentParameters)
    return num_pars
end

# Particular scalar parameter types
function num_parameters(model::InfiniteModel,
                        type::Type{C})::Int where {C <: ScalarParameter}
    return length(_data_dictionary(model, type))
end

# ScalarParameter
function num_parameters(model::InfiniteModel,
                        type::Type{ScalarParameter})::Int
    num_pars = num_parameters(model, FiniteParameter)
    num_pars += num_parameters(model, IndependentParameter)
    return num_pars
end

# DependentParameters
function num_parameters(model::InfiniteModel,
                        type::Type{DependentParameters})::Int
    num_pars = 0
    for (_, object) in _data_dictionary(model, type)
        num_pars += length(object.names)
    end
    return num_pars
end

# InfiniteParameter
function num_parameters(model::InfiniteModel,
                        type::Type{InfiniteParameter})::Int
    num_pars = num_parameters(model, IndependentParameter)
    num_pars += num_parameters(model, DependentParameters)
    return num_pars
end

"""
    all_parameters(model::InfiniteModel,
                   type::Type{InfOptParameter} = InfOptParameter
                   )::Vector{GeneralVariableRef}

Return a list of all the `InfiniteOpt` parameters assigned to `model`. By default,
all of the infinite and finite parameters is returned. The search is reduced to
a particular type is obtained by specifying the concrete parameter type
of [`InfOptParameter`](@ref) via `type`. Type options include:
- `InfOptParameter`: all parameters
- `ScalarParameter`: all scalar parameters
- `InfiniteParameter`: all infinite parameters
- `FiniteParameter`: all finite parameters
- `IndependentParameter`: all independent infinite parameters
- `DependentParameters`: all dependent infinite parameters

**Examples**
```julia-repl
julia> all_parameters(model)
4-element Array{GeneralVariableRef,1}:
 t
 x[1]
 x[2]
 alpha

julia> all_parameters(model, FiniteParameter)
1-element Array{GeneralVariableRef,1}:
 alpha
```
"""
function all_parameters(model::InfiniteModel,
                        type::Type{InfOptParameter} = InfOptParameter
                        )::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, IndependentParameter)
    append!(prefs_list, all_parameters(model, DependentParameters))
    append!(prefs_list, all_parameters(model, FiniteParameter))
    return prefs_list
end

# Particular scalar parameter types
function all_parameters(model::InfiniteModel,
                        type::Type{C}
                        )::Vector{GeneralVariableRef} where {C <: InfOptParameter}
    prefs_list = Vector{GeneralVariableRef}(undef, num_parameters(model, type))
    counter = 1
    for (index, _) in _data_dictionary(model, type)
        prefs_list[counter] = _make_parameter_ref(model, index)
        counter += 1
    end
    return prefs_list
end

# ScalarParameter
function all_parameters(model::InfiniteModel,
                        type::Type{ScalarParameter})::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, IndependentParameter)
    append!(prefs_list, all_parameters(model, FiniteParameter))
    return prefs_list
end

# DependentParameters
function all_parameters(model::InfiniteModel,
                        type::Type{DependentParameters}
                        )::Vector{GeneralVariableRef}
    prefs_list = Vector{GeneralVariableRef}(undef, num_parameters(model, type))
    counter = 1
    for (index, object) in _data_dictionary(model, type)
        for i in eachindex(object.names)
            dep_idx = DependentParameterIndex(index, i)
            prefs_list[counter] = _make_parameter_ref(model, dep_idx)
            counter += 1
        end
    end
    return prefs_list
end

# InfiniteParameter
function all_parameters(model::InfiniteModel,
                        type::Type{InfiniteParameter}
                        )::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, IndependentParameter)
    append!(prefs_list, all_parameters(model, DependentParameters))
    return prefs_list
end

################################################################################
#                                 DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel,
                prefs::AbstractArray{<:DependentParameterRef})::Nothing

Extend `JuMP.delete` to delete
dependent infinite parameters and their dependencies. All variables, constraints, and
measure functions that depend on `prefs` are updated to exclude them. Errors if the
parameters are contained in an `AbstractMeasureData` datatype that is employed by
a measure since the measure becomes invalid otherwise. Thus, measures that
contain this dependency must be deleted first. Note that
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) needs to be
extended to allow deletion of parameters when custom `AbstractMeasureData`
datatypes are used. Note that any dependent infinite variables will have their
start values reset to the default via [`reset_start_value_function`](@ref).

**Example**
```julia-repl
julia> print(model)
Min measure(g(t, x)*t + x) + z
Subject to
 z ≥ 0.0
 g(t, x) + z ≥ 42.0, ∀ t ∈ [0, 6], x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
 g(0.5, x) = 0, x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]

julia> delete(model, x)

julia> print(model)
Min measure(g(t)*t) + z
Subject to
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]
 g(0.5) = 0
```
"""
function JuMP.delete(model::InfiniteModel,
                     prefs::AbstractArray{<:DependentParameterRef})::Nothing
    @assert JuMP.is_valid(model, first(prefs)) "Parameter references are invalid."
    _check_complete_param_array(prefs)
    gvrefs = [_make_parameter_ref(model, JuMP.index(pref)) for pref in prefs]
    # ensure deletion is okay (prefs are not used by measure data)
    for pref in gvrefs
        for mindex in _measure_dependencies(pref)
            data = measure_data(dispatch_variable_ref(model, mindex))
            _check_param_in_data(pref, data)
        end
    end
    # update optimizer model status
    if any(is_used(pref) for pref in prefs)
        set_optimizer_model_ready(model, false)
    end
    # delete dependence of measures and constraints on prefs
    for pref in gvrefs
        _update_measures(model, pref)
        _update_constraints(model, pref)
    end
    # get the object and parameter numbers
    obj_num = _object_number(first(prefs))
    param_nums = collect(_data_object(first(prefs)).parameter_nums)
    # update infinite variables that depend on pref
    for vindex in _infinite_variable_dependencies(first(prefs))
        # remove the parameter dependences
        vref = InfiniteVariableRef(model, vindex)
        vprefs = raw_parameter_refs(vref)
        tup_index = findfirst(isequal(obj_num), _object_number.(vprefs[:, 1]))
        delete_indices = vprefs.ranges[tup_index]
        deleteat!(vprefs, tup_index, tuple_index = true)
        reset_start_value_function(vref)
        # update any point variables that depend on vref accordingly
        for pindex in _point_variable_dependencies(vref)
            pvref = PointVariableRef(model, pindex)
            deleteat!(raw_parameter_values(pvref), delete_indices)
        end
        # update any reduced variables that depend on vref accordingly
        for rindex in _reduced_variable_dependencies(vref)
            rvref = ReducedVariableRef(model, rindex)
            _update_reduced_variable(rvref, delete_indices)
        end
    end
    # delete parameter information stored in model
    _delete_data_object(first(prefs))
    # update the object numbers and parameter numbers
    _update_model_numbers(model, obj_num, param_nums)
    return
end
