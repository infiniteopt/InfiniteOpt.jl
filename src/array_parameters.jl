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
                          object::MultiParameterData{<:DependentParameters}
                          )::DependentParametersIndex
    return MOIUC.add_item(model.dependent_params, object)
end

# Extend _data_dictionary
function _data_dictionary(pref::DependentParameterRef)::MOIUC.CleverDict
    return model.dependent_params
end

# Extend _data_object
function _data_object(pref::DependentParameterRef
                      )::MultiParameterData{<:DependentParameters}
    return _data_dictionary(pref)[JuMP.index(pref)]
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
    elseif !(params isa AbstractArray{<:_DependentParameter{<:set_type}})
        _error("Cannot specify multiple multi-dimensional distributions sets for one " *
               "set of dependent infinite parameters.")
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
    elseif !(params isa AbstractArray{<:_DependentParameter{<:set_type}})
        _error("Cannot specify multiple `CollectionSet`s for one group of " *
               "dependent infinite parameters.")
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
    set_type = typeof(first(params).set)
    if !(params isa AbstractArray{<:_DependentParameter{<:set_type}})
        _error("Cannot specify multiple `InfiniteArraySet`s for one group of " *
               "dependent infinite parameters.")
    end
    return
end

# Fallback
function _check_param_sets(_error::Function, params)
    _error("Unrecognized infinite set input.")
end

## Use set type dispatch to make the proper InfiniteArraySet
# InfiniteArraySet
function _make_array_set(params::Vector{<:_DependentParameter{T}
                         )::T where {T <: InfiniteArraySet}
    return first(params).set
end

# InfiniteScalarSets
function _make_array_set(params::Vector{<:_DependentParameter{T}
                         )::CollectionSet{T} where {T <: InfiniteScalarSet}
    return CollectionSet([p.set for p in params])
end

# Build a DependentParameters object given an array of _DependentParameters
function _build_parameters(_error::Function,
                           params::AbstractArray{<:_DependentParameter};
                           num_supports::Int = 0, sig_figs::Int = 5,
                           extra_kw_args...)::DependentParameters
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
        supps = collect(transpose(trans_supps))
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
infinite parameter references.
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
        check = :(isa($(info.lower_bound), Number))
        return :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a number."))
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
#                             OTHER STUFF
################################################################################
"""
    supports(prefs::AbstractArray{<:ParameterRef})::Vector

Return the support points associated with an array of `prefs` formatted as a
vector of supports following the format of the input array. If the
parameters are not independent then the supports of each parameter are simply
spliced together. Alternatively can call `supports.()` to more efficiently obtain
an array of the same input format whose parameter references have been replaced
with their supports (so long as the prefs are not independent). Errors if all
the parameter references do not have the same
group ID number (were intialized together as an array) or if the nonindependent
parameters have support vectors of different lengths. If the parameters are
independent then all the unique combinations are identified and returned as
supports. Warning this operation is computationally very expensive if there
exist a large number of combinations.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> x = @infinite_parameter(model, [i = 1:2], set = IntervalSet(-1, 1),
                               base_name = "x", independent = true)
2-element Array{ParameterRef,1}:
 x[1]
 x[2]

julia> set_supports.(x, [[-1, 1], [-1, 1]]);

julia> supports(x)
4-element Array{Array{Int64,1},1}:
 [-1, -1]
 [1, -1]
 [-1, 1]
 [1, 1]
```
"""
function supports(prefs::AbstractArray{<:ParameterRef})::Vector
    _allequal(group_id.(prefs)) || error("Array contains parameters from multiple" *
                                         " groups.")
    all(has_supports(prefs[k]) for k in keys(prefs)) || error("Not all parameters have supports.")
    lengths = [num_supports(prefs[k]) for k in keys(prefs)]
    indices = Collections._get_indices(prefs)
    prefs = Collections._make_ordered(prefs, indices)
    if !is_independent(first(prefs))
        _allequal(lengths) || error("Each nonindependent parameter must have " *
                                    "the same number of support points.")
        support_list = unique([[supports(prefs[k])[i] for k in keys(prefs)] for i in 1:lengths[1]])
    else
        all_supports = [supports(prefs[k]) for k in keys(prefs)]
        combos = Iterators.product(all_supports...)
        support_list = [[combo...] for combo in Iterators.take(combos, length(combos))]
    end
    range = 1:length(prefs)
    return [Collections._make_array(supp, range, indices) for supp in support_list]
end


"""
    delete_supports(pref::ParameterRef)

Delete the support points for `pref`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @infinite_parameter(model, t in [0, 1], supports = [0, 1]))
julia> delete_supports(t)

julia> supports(t)
ERROR: Parameter t does not have supports.
```
"""
function delete_supports(pref::ParameterRef)
    _update_parameter_supports(pref, Int[])
    return
end

# Multivariate distribution sets
function generate_and_add_supports!(pref::ParameterRef,
                                    set::DistributionSet{<:Distributions.MultivariateDistribution};
                                    num_supports::Int = 10, sig_fig::Int = 5)
    pref_group_id = group_id(pref)
    model = JuMP.owner_model(pref)
    associated_p_index = sort([i for i in 1:length(model.params)
                               if model.param_to_group_id[i] == pref_group_id])
    new_supports = generate_support_values(set, num_supports = num_supports, sig_fig = sig_fig)

    for i in 1:length(associated_p_index)
        pref_i = ParameterRef(model, associated_p_index[i])
        add_supports(pref_i, new_supports[i, :])
    end
    return
end

## Define functions to extract the names of parameters
# Extract the root name of a parameter reference
function _root_name(pref::ParameterRef)::String
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
