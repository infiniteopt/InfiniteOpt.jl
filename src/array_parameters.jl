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
    return MOIUC.add_item(model.dependent_params, object)
end

# Extend _data_dictionary
function _data_dictionary(pref::DependentParameterRef)::MOIUC.CleverDict
    return JuMP.owner_model(pref).dependent_params
end

# Extend _data_object
function _data_object(pref::DependentParameterRef
                      )::MultiParameterData{<:InfiniteArraySet}
    return _data_dictionary(pref)[JuMP.index(pref).object_index]
end

# Extend _core_variable_object
function _core_variable_object(pref::DependentParameterRef)::DependentParameters
    return _data_object(pref).parameters
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
                           sig_figs::Int = 5,
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
        labels = [Set([UserDefined]) for i in 1:size(supps, 2)]
    # we want to generate supports
    elseif num_supports != 0
        supps, label = generate_support_values(set, num_supports = num_supports,
                                               sig_figs = sig_figs)
        labels = [label for i in 1:num_supports]
    # no supports are specified
    else
        supps = zeros(Float64, length(vector_params), 0)
        labels = Set{Symbol}[]
    end
    # make the parameter object
    names = [param.name for param in vector_params]
    return DependentParameters(set, supps, labels), names, indices
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

julia> params = DependentParameters(set, rand(dist, 10), [Set([McSample]) for i = 1:10]);

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
    num_params = size(params.supports, 1)
    # process the names
    if isempty(names)
        names = ["noname" for i in 1:num_params]
    end
    # process the indices
    if indices === nothing
        indices = CartesianIndices(1:num_params)
    end
    # make the parameter model object
    obj_num = model.last_object_num += 1
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
    return _data_object(pref).names[_param_index(pref)]
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
    return _data_object(pref).infinite_var_indices[_param_index(pref)]
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
    return used_by_measure(vref) || used_by_constraint(vref) ||
           used_by_infinite_variable(vref)
end

################################################################################
#                          PARAMETER OBJECT METHODS
################################################################################
# Extend _parameter_number
function _parameter_number(pref::DependentParameterRef)::Int
    return _data_object(pref).parameter_nums[_param_index(pref)]
end

# Extend _object_number
function _object_number(pref::DependentParameterRef)::Int
    return _data_object(pref).object_num
end

# Extend _set_core_variable_object
function _set_core_variable_object(pref::DependentParameterRef,
                                   params::DependentParameters)::Nothing
    _data_object(pref).parameters = object
    return
end

################################################################################
#                             INFINITE SET METHODS
################################################################################
## Get the individual infinite set if possible
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
    return _parameter_set(_core_variable_object(pref).set, pref)
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
    if length(prefs) == length(_data_object(first(prefs)).names)
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    end
    return _core_variable_object(first(prefs)).set
end

# Update the underlying set and delete the supports
function _update_parameter_set(pref::DependentParameterRef,
                               new_set::InfiniteArraySet)::Nothing
    old_params = _core_variable_object(pref)
    new_supports = zeros(size(old_params.supports, 1), 0)
    new_params = DependentParameters(new_set, new_supports,
                                     empty!(old_params.labels))
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
by a measure.

**Example**
```julia-repl
julia> set_infinite_set(x[1], IntervalSet(0, 2))

julia> infinite_set(x[1])
[0, 2]
```
"""
function set_infinite_set(pref::DependentParameterRef,
                          set::InfiniteScalarSet)::Nothing
    old_set = _core_variable_object(pref).set
    if old_set != CollectionSet
        error("Cannot set the individual infinite set of $pref if the " *
              "underlying set is not a CollectionSet.")
    elseif used_by_measure(pref)
        error("Cannot override the infinite set of $pref since it is used by " *
              "a measure.")
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
of the dependent infinite parameters are included or if any of them are used by
measures.

**Example**
```julia-repl
julia> set_infinite_set(x, CollectionSet([IntervalSet(0, 1), IntervalSet(0, 2)]))
```
"""
function set_infinite_set(prefs::AbstractArray{<:DependentParameterRef},
                          set::InfiniteArraySet)::Nothing
    if any(used_by_measure(pref) for pref in prefs)
        error("Cannot override the infinite set of $pref since it is used by " *
              "a measure.")
    end
    old_set = infinite_set(prefs) # this error checks
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
    set = _core_variable_object(pref).set
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
function _parameter_supports(pref::DependentParameterRef)::Array{Float64, 2}
    return _core_variable_object(pref).supports
end

# Get the raw support labels
function _parameter_support_labels(pref::DependentParameterRef)::Vector{Set{Symbol}}
    return _core_variable_object(pref).labels
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

julia> num_supports(x[1], label = McSample)
0
```
"""
function num_supports(pref::DependentParameterRef; label::Symbol = All)::Int
    if label == All
        return size(_parameter_supports(pref), 2)
    else
        return length(findall(x -> label in x, _parameter_support_labels(pref)))
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
    if length(prefs) == length(_data_object(first(prefs)).names)
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    end
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
has_supports(pref::DependentParameterRef)::Bool = num_supports(pref) > 0

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
    return num_supports(prefs) > 0
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
        return _parameter_supports(pref)[_param_index(pref), :]
    else
        indices = findall(x -> label in x, _parameter_support_labels(pref))
        return _parameter_supports(pref)[_param_index(pref), indices]
    end
end

"""
    supports(prefs::AbstractArray{<:DependentParameterRef};
             [label::Symbol = All])::AbstractArray{<:Float64}

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
                  label::Symbol = All)::AbstractArray{<:Float64}
    if length(prefs) == length(_data_object(first(prefs)).names)
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    end
    return supports.(prefs, label = label) # TODO make more efficient
end

# More efficient dispatch for Vectors
function supports(prefs::Vector{DependentParameterRef};
                  label::Symbol = All)::Array{Float64, 2}
    if length(prefs) == length(_data_object(first(prefs)).names)
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    end
    if label == All
        return _parameter_supports(pref)
    else
        indices = findall(x -> label in x, _parameter_support_labels(pref))
        return _parameter_supports(pref)[:, indices]
    end
end

# Define method for overriding the current supports
function _update_parameter_supports(prefs::AbstractArray{<:DependentParameterRef},
                                    supports::Array{<:Real, 2},
                                    labels::Vector{Set{Symbol}})::Nothing
    old_params = _core_variable_object(first(prefs))
    new_params = DependentParameters(old_params.set, supports, labels)
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
    keys(supports) == keys(prefs) || error("Inconsistent support indices")
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
    labels = [Set([label]) for i = 1:size(supports, 2)]
    _update_parameter_supports(prefs, supports, labels)
    return
end

# Error for single dependent parameters
function set_supports(pref::DependentParameterRef, supports; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    add_supports(prefs::AbstractArray{<:IndependentParameterRef},
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
function add_supports(prefs::AbstractArray{<:IndependentParameterRef},
                      supports::AbstractArray{<:Vector{<:Real}};
                      label::Symbol = UserDefined, # interal keyword args
                      check::Bool = true)::Nothing
    supps = _make_support_matrix(prefs, supports)
    add_supports(_make_vector(prefs), supps, label = label, check = check)
    return
end

# More efficient version for supports in the correct format
function add_supports(prefs::Vector{IndependentParameterRef},
                      supports::Array{<:Real, 2};
                      label::Symbol = UserDefined, # interal keyword args
                      check::Bool = true)::Nothing
    set = infinite_set(prefs) # this does a check on prefs
    if check && !supports_in_set(supports, set)
        error("Supports violate the domain of the infinite set.")
    end
    current_supports = _parameter_supports(first(prefs))
    new_supports = hcat(current_supports, supports)
    current_labels = _parameter_support_labels(first(prefs))
    labels = push!(current_labels, [Set([label]) for i = 1:size(supports, 2)])
    _update_parameter_supports(prefs, new_supports, labels)
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
    if length(prefs) == length(_data_object(first(prefs)).names)
        error("Dimensions of parameter container and the infinite set do not " *
              "match, ensure all related dependent parameters are included.")
    elseif any(used_by_measure(pref) for pref in prefs)
        error("Cannot delete supports with measure dependencies.")
    end
    _update_parameter_supports(prefs, zeros(length(prefs), 0), [Set{Symbol}()])
    return
end

# TODO add support generation methods

################################################################################
#                               OTHER METHODS
################################################################################
# Extract the root name of a parameter reference
function _root_name(pref::GeneralVariableRef)::String
    name = JuMP.name(pref)
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

################################################################################
#                                 DELETION
################################################################################
# TODO add deletion
