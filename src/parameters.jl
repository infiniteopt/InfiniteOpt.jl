# Define symbol inputs for different variable types
const Parameter = :Parameter

# Extend Base.copy for new variable types
function Base.copy(v::ParameterRef, new_model::InfiniteModel)
    return ParameterRef(new_model, v.index)
end

mutable struct _ParameterInfoExpr
    has_lb::Bool
    lower_bound::Any
    has_ub::Bool
    upper_bound::Any
    has_dist::Bool
    distribution::Any
    has_set::Bool
    set::Any
end

function _is_set_keyword(kw::Expr)
    return kw.args[1] in [:set, :lower_bound, :upper_bound, :distribution]
end

function _ParameterInfoExpr(; lower_bound = NaN, upper_bound = NaN,
                            distribution = NaN, set = NaN)
    # isnan(::Expr) is not defined so we need to do !== NaN
    _ParameterInfoExpr(lower_bound !== NaN, lower_bound, upper_bound !== NaN,
                       upper_bound, distribution !== NaN, distribution,
                       set !== NaN, set)
end

function JuMP._set_lower_bound_or_error(_error::Function,
                                        info::_ParameterInfoExpr, lower)
    info.has_lb && _error("Cannot specify parameter lower_bound twice")
    info.has_dist && _error("Cannot specify parameter lower_bound and " *
                            "distribution")
    info.has_set && _error("Cannot specify parameter lower_bound and set")
    info.has_lb = true
    info.lower_bound = lower
end

function JuMP._set_upper_bound_or_error(_error::Function,
                                        info::_ParameterInfoExpr, upper)
    info.has_ub && _error("Cannot specify parameter upper_bound twice")
    info.has_dist && _error("Cannot specify parameter upper_bound and " *
                            "distribution")
    info.has_set && _error("Cannot specify parameter upper_bound and set")
    info.has_ub = true
    info.upper_bound = upper
end

function _dist_or_error(_error::Function, info::_ParameterInfoExpr, dist)
    info.has_dist && _error("Cannot specify parameter distribution twice")
    (info.has_lb || info.has_ub) && _error("Cannot specify parameter " *
                                           "distribution and upper/lower bounds")
    info.has_set && _error("Cannot specify parameter distribution and set")
    info.has_dist = true
    info.distribution = dist
end

function _set_or_error(_error::Function, info::_ParameterInfoExpr, set)
    info.has_set && _error("Cannot specify variable fixed value twice")
    (info.has_lb || info.has_ub) && _error("Cannot specify parameter set and " *
                                           "upper/lower bounds")
    info.has_dist && _error("Cannot specify parameter set and distribution")
    info.has_set = true
    info.set = set
end

function _constructor_set(_error::Function, info::_ParameterInfoExpr)
    if (info.has_lb || info.has_ub) && !(info.has_lb && info.has_ub)
        _error("Must specify both an upper bound and a lower bound")
    elseif info.has_lb
        check = :(isa($(info.lower_bound), Number))
        return :($(check) ? IntervalSet(convert(Float64, $(info.lower_bound)), convert(Float64, $(info.upper_bound))) : error("Bounds must be a number."))
    elseif info.has_dist
        check = :(isa($(info.distribution), Distributions.NonMatrixDistribution))
        return :($(check) ? DistributionSet($(info.distribution)) : error("Distribution must be a subtype of Distributions.NonMatrixDistribution."))
    elseif info.has_set
        check = :(isa($(info.set), AbstractInfiniteSet))
        return :($(check) ? $(info.set) : error("Set must be a subtype of AbstractInfiniteSet."))
    else
        _error("Must specify upper/lower bounds, a distribution, or a set")
    end
end

function _check_supports_in_bounds(_error::Function,
                                   supports::Union{Number, Vector{<:Number}},
                                   set::AbstractInfiniteSet)
    min_support = minimum(supports)
    max_support = maximum(supports)
    if isa(set, IntervalSet)
       if min_support < set.lower_bound || max_support > set.upper_bound
           _error("Support points violate the interval set bounds.")
       end
    elseif isa(set, DistributionSet)
       if typeof(set.distribution) <: Distributions.UnivariateDistribution
           if min_support < minimum(set.distribution) || max_support > maximum(set.distribution)
               _error("Support points violate the distribution set bounds.")
           end
       end
    end
    return
end

"""
build_parameter(_error::Function, set::AbstractInfiniteSet, num_params::Int;
                supports::Union{Number, Vector{<:Number}} = Number[],
                independent::Bool = false, extra_kw_args...)
Build an infinite parameter to the model in a manner similar to `JuMP.build_variable`.
"""
function build_parameter(_error::Function, set::AbstractInfiniteSet, num_params::Int;
                         supports::Union{Number, Vector{<:Number}} = Number[],
                         independent::Bool = false, extra_kw_args...)
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    if length(supports) != 0
        _check_supports_in_bounds(_error, supports, set)
    end
    if isa(set, DistributionSet)
        if isa(set.distribution, Distributions.MultivariateDistribution)
            if num_params != length(set.distribution)
                _error("Multivariate distribution must match dimension of parameter.")
            end
        end
    end
    unique_supports = sort(unique(supports))
    if length(unique_supports) != length(supports)
        @warn("Support points are not unique, eliminating redundant points.")
    end
    return InfOptParameter(set, unique_supports, independent)
end

"""
    add_parameter(model::InfiniteModel, v::InfOptParameter, name::String="")
Add an infinite parameter to the model in a manner similar to `JuMP.add_variable`.
"""
function add_parameter(model::InfiniteModel, v::InfOptParameter, name::String="")
    model.next_param_index += 1
    pref = ParameterRef(model, model.next_param_index)
    model.params[JuMP.index(pref)] = v
    model.param_to_group_id[JuMP.index(pref)] = model.next_param_id
    JuMP.set_name(pref, name)
    return pref
end

"""
    used_by_constraint(pref::ParameterRef)::Bool
Return a Boolean indicating if `pref` is used by a constraint.
"""
function used_by_constraint(pref::ParameterRef)::Bool
    return haskey(JuMP.owner_model(pref).param_to_constrs, JuMP.index(pref))
end

"""
    used_by_measure(pref::ParameterRef)::Bool
Return a Boolean indicating if `pref` is used by a measure.
"""
function used_by_measure(pref::ParameterRef)::Bool
    return haskey(JuMP.owner_model(pref).param_to_meas, JuMP.index(pref))
end

"""
    used_by_variable(pref::ParameterRef)::Bool
Return a Boolean indicating if `pref` is used by an infinite variable.
"""
function used_by_variable(pref::ParameterRef)::Bool
    return haskey(JuMP.owner_model(pref).param_to_vars, JuMP.index(pref))
end

"""
    is_used(pref::ParameterRef)::Bool
Return a Boolean indicating if `pref` is used in the model.
"""
function is_used(pref::ParameterRef)::Bool
    return used_by_measure(pref) || used_by_constraint(pref) || used_by_variable(pref)
end

# Define internal functions for deleting elements of parameter tuples
function _remove_parameter(pref::ParameterRef, delete_pref::ParameterRef)
    if pref == delete_pref
        return
    else
        return pref
    end
end

function _remove_parameter(arr::JuMP.Containers.SparseAxisArray{<:ParameterRef},
                           delete_pref::ParameterRef)
    data = filter(v -> v.second != delete_pref, arr.data)
    arr = JuMP.Containers.SparseAxisArray(data)
    return length(arr) != 0 ? arr : nothing
end

function _remove_parameter(prefs::Tuple, delete_pref::ParameterRef)
    pref_list = Any[_remove_parameter(pref, delete_pref) for pref in prefs]
    filter!(x -> x != nothing, pref_list)
    pref_tuple = ()
    for pref in pref_list
        pref_tuple = (pref_tuple..., pref)
    end
    return pref_tuple
end

function _root_name(vref::InfiniteVariableRef)
    name = JuMP.name(vref)
    return name[1:findfirst(isequal('('), name)-1]
end

"""
    JuMP.delete(model::InfiniteModel, pref::ParameterRef)
Extend the `JuMP.delete` function to accomodate infinite parameters
"""
function JuMP.delete(model::InfiniteModel, pref::ParameterRef)
    @assert JuMP.is_valid(model, pref)
    if is_used(pref)
        set_optimizer_model_status(model, false)
    end
    if used_by_variable(pref)
        for vindex in model.param_to_vars[JuMP.index(pref)]
            prefs = _remove_parameter(model.vars[vindex].parameter_refs, pref)
            vref = InfiniteVariableRef(model, vindex)
            _update_variable_param_refs(vref, prefs)
            JuMP.set_name(vref, _root_name(vref))
            if used_by_measure(vref)
                for mindex in model.var_to_meas[JuMP.index(vref)]
                    measure = model.measures[mindex]
                    mref = MeasureRef(model, mindex)
                    JuMP.set_name(mref, _make_meas_name(measure))
                end
            end
        end
        delete!(model.param_to_vars, JuMP.index(pref))
    end
    if used_by_measure(pref)
        for mindex in model.param_to_meas[JuMP.index(pref)]
            if isa(model.measures[mindex].func, ParameterRef)
                data = model.measures[mindex].data
                model.measures[mindex] = Measure(zero(JuMP.AffExpr), data)
            else
                _remove_variable(model.measures[mindex].func, pref)
            end
            measure = model.measures[mindex]
            mref = MeasureRef(model, mindex)
            JuMP.set_name(mref, _make_meas_name(measure))
        end
        delete!(model.param_to_meas, JuMP.index(pref))
    end
    if used_by_constraint(pref)
        for cindex in model.param_to_constrs[JuMP.index(pref)]
            if isa(model.constrs[cindex].func, ParameterRef)
                set = model.constrs[cindex].set
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr), set)
            else
                _remove_variable(model.constrs[cindex].func, pref)
            end
        end
        delete!(model.param_to_constrs, JuMP.index(pref))
    end
    delete!(model.params, JuMP.index(pref))
    delete!(model.param_to_name, JuMP.index(pref))
    delete!(model.param_to_group_id, JuMP.index(pref))
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, pref::ParameterRef)
Extend the `JuMP.is_valid` function to accomodate infinite parameters.
"""
function JuMP.is_valid(model::InfiniteModel, pref::ParameterRef)
    return (model === JuMP.owner_model(pref) && JuMP.index(pref) in keys(model.params))
end

"""
    JuMP.name(pref::ParameterRef)
Extend the `JuMP.name` function to accomodate infinite parameters
"""
JuMP.name(pref::ParameterRef) = JuMP.owner_model(pref).param_to_name[JuMP.index(pref)]

"""
    JuMP.set_name(pref::ParameterRef, name::String)
Extend the `JuMP.set_name` function to accomodate infinite parameters.
"""
function JuMP.set_name(pref::ParameterRef, name::String)
    JuMP.owner_model(pref).param_to_name[JuMP.index(pref)] = name
    JuMP.owner_model(pref).name_to_param = nothing
    return
end

"""
    num_parameters(model::InfiniteModel)
Return the number of infinite parameters.
"""
num_parameters(model::InfiniteModel) = length(model.params)

# Internal functions
_parameter_set(pref::ParameterRef) = JuMP.owner_model(pref).params[JuMP.index(pref)].set
_parameter_supports(pref::ParameterRef) = JuMP.owner_model(pref).params[JuMP.index(pref)].supports
function _update_parameter_set(pref::ParameterRef, set::AbstractInfiniteSet)
    supports = JuMP.owner_model(pref).params[JuMP.index(pref)].supports
    independent = JuMP.owner_model(pref).params[JuMP.index(pref)].independent
    JuMP.owner_model(pref).params[JuMP.index(pref)] = InfOptParameter(set, supports, independent)
    if is_used(pref)
        set_optimizer_model_status(JuMP.owner_model(pref), false)
    end
    return
end
function _update_parameter_supports(pref::ParameterRef, supports::Vector{<:Number})
    set = JuMP.owner_model(pref).params[JuMP.index(pref)].set
    independent = JuMP.owner_model(pref).params[JuMP.index(pref)].independent
    JuMP.owner_model(pref).params[JuMP.index(pref)] = InfOptParameter(set, supports, independent)
    if is_used(pref)
        set_optimizer_model_status(JuMP.owner_model(pref), false)
    end
    return
end

"""
    infinite_set(pref::ParameterRef)::AbstractInfiniteSet
Return the infinite set of `pref`.
"""
function infinite_set(pref::ParameterRef)::AbstractInfiniteSet
    return _parameter_set(pref)
end

"""
    set_infinite_set(pref::ParameterRef, set::AbstractInfiniteSet)
Specify the infinite set of `pref`.
"""
function set_infinite_set(pref::ParameterRef, set::AbstractInfiniteSet)
    _update_parameter_set(pref, set)
    return
end

"""
    JuMP.has_lower_bound(pref::ParameterRef)
Extend the `JuMP.has_lower_bound` function to accomodate infinite parameters.
"""
function JuMP.has_lower_bound(pref::ParameterRef)
    set = _parameter_set(pref)
    if isa(set, IntervalSet)
        return true
    elseif isa(set, DistributionSet)
        if typeof(set.distribution) <: Distributions.UnivariateDistribution
            return true
        else
            false
        end
    else
        type = typeof(set)
        error("Undefined infinite set type $type for lower bound checking.")
    end
end

"""
    JuMP.lower_bound(pref::ParameterRef)::Number
Extend the `JuMP.lower_bound` function to accomodate infinite parameters.
"""
function JuMP.lower_bound(pref::ParameterRef)::Number
    set = _parameter_set(pref)
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    if isa(set, IntervalSet)
        return set.lower_bound
    else isa(set, DistributionSet)
        if typeof(set.distribution) <: Distributions.UnivariateDistribution
            return minimum(set.distribution)
        end
    end
end

"""
    JuMP.set_lower_bound(pref::ParameterRef, lower::Number)
Extend the `JuMP.set_lower_bound` function to accomodate infinite parameters.
"""
function JuMP.set_lower_bound(pref::ParameterRef, lower::Number)
    set = _parameter_set(pref)
    if isa(set, DistributionSet)
        error("Cannot set the lower bound of a distribution, try using `Distributions.Truncated` instead.")
    elseif !isa(set, IntervalSet)
        error("Parameter $(pref) is not an interval set.")
    end
    _update_parameter_set(pref, IntervalSet(Float64(lower), set.upper_bound))
    return
end

"""
    JuMP.has_upper_bound(pref::ParameterRef)
Extend the `JuMP.has_upper_bound` function to accomodate infinite parameters.
"""
function JuMP.has_upper_bound(pref::ParameterRef)
    set = _parameter_set(pref)
    if isa(set, IntervalSet)
        return true
    elseif isa(set, DistributionSet)
        if typeof(set.distribution) <: Distributions.UnivariateDistribution
            return true
        else
            false
        end
    else
        type = typeof(set)
        error("Undefined infinite set type $type for lower upper checking.")
    end
end

"""
    JuMP.upper_bound(pref::ParameterRef)::Number
Extend the `JuMP.upper_bound` function to accomodate infinite parameters.
"""
function JuMP.upper_bound(pref::ParameterRef)::Number
    set = _parameter_set(pref)
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    if isa(set, IntervalSet)
        return set.upper_bound
    else isa(set, DistributionSet)
        if typeof(set.distribution) <: Distributions.UnivariateDistribution
            return maximum(set.distribution)
        end
    end
end

"""
    JuMP.set_upper_bound(pref::ParameterRef, lower::Number)
Extend the `JuMP.set_upper_bound` function to accomodate infinite parameters.
"""
function JuMP.set_upper_bound(pref::ParameterRef, upper::Number)
    set = _parameter_set(pref)
    if isa(set, DistributionSet)
        error("Cannot set the upper bound of a distribution, try using `Distributions.Truncated` instead.")
    elseif !isa(set, IntervalSet)
        error("Parameter $(pref) is not an interval set.")
    end
    _update_parameter_set(pref, IntervalSet(set.lower_bound, Float64(upper)))
    return
end

"""
    num_supports(pref::ParameterRef)
Return the number of support points associated with `pref`.
"""
num_supports(pref::ParameterRef) = length(_parameter_supports(pref))

"""
    has_supports(pref::ParameterRef)
Return true if `pref` has supports or false otherwise.
"""
has_supports(pref::ParameterRef) = num_supports(pref) > 0

"""
    supports(pref::ParameterRef)
Return the support points associated with `pref`.
"""
function supports(pref::ParameterRef)
    !has_supports(pref) && error("Parameter $pref does not have supports.")
    return _parameter_supports(pref)
end

"""
    supports(prefs::AbstractArray{<:ParameterRef})
Return the support points associated with an array of `prefs`. Error if not from
same group.
"""
function supports(prefs::AbstractArray{<:ParameterRef})
    prefs = convert(JuMP.Containers.SparseAxisArray, prefs)
    !_only_one_group(prefs) && error("Array contains parameters from multiple groups.")
    for (k, pref) in prefs.data
        !has_supports(pref) && error("Parameter $pref does not have supports.")
    end
    lengths = [num_supports(pref) for (k, pref) in prefs.data]
    if !is_independent(collect(values(prefs.data))[1])
        length(unique(lengths)) != 1 && error("Each nonindependent parameter must " *
                                              "have the same number of support " *
                                              "points.")
        support_list = Vector{JuMP.Containers.SparseAxisArray}(undef, lengths[1])
        for i = 1:length(support_list)
            support_list[i] = JuMP.Containers.SparseAxisArray(Dict(k => supports(pref)[i] for (k, pref) in prefs.data))
        end
    else
        all_keys = collect(keys(prefs))
        all_supports = [supports(pref) for (k, pref) in prefs.data]
        support_list = Vector{JuMP.Containers.SparseAxisArray}(undef, prod(lengths))
        counter = 1
        for combo in Iterators.product(all_supports...)
            support_list[counter] = JuMP.Containers.SparseAxisArray(Dict(all_keys[i] => combo[i] for i = 1:length(combo)))
            counter += 1
        end
    end
    return support_list
end

"""
    set_supports(pref::ParameterRef, supports::Vector{<:Number})
Specify the support points for `pref`.
"""
function set_supports(pref::ParameterRef, supports::Vector{<:Number})
    set = _parameter_set(pref)
    _check_supports_in_bounds(error, supports, set)
    unique_supports = sort(unique(supports))
    if length(unique_supports) != length(supports)
        @warn("Support points are not unique, eliminating redundant points.")
    end
    _update_parameter_supports(pref, unique_supports)
    return
end

"""
    add_supports(pref::ParameterRef, supports::Union{Number, Vector{<:Number}})
Add additional support points for `pref`.
"""
function add_supports(pref::ParameterRef, supports::Union{Number, Vector{<:Number}})
    set = _parameter_set(pref)
    _check_supports_in_bounds(error, supports, set)
    current_supports = _parameter_supports(pref)
    _update_parameter_supports(pref, sort(unique([current_supports; supports])))
    return
end

"""
    delete_supports(pref::ParameterRef)
Delete the support points for `pref`.
"""
function delete_supports(pref::ParameterRef)
    _update_parameter_supports(pref, Number[])
    return
end

"""
    group_id(pref::ParameterRef)::Int
Return the group ID number for `pref`.
"""
function group_id(pref::ParameterRef)::Int
    return JuMP.owner_model(pref).param_to_group_id[JuMP.index(pref)]
end

"""
    group_id(prefs::AbstractArray{<:ParameterRef})::Int
Return the group ID number for a group of `prefs`. Error if contains multiple groups.
"""
function group_id(prefs::AbstractArray{<:ParameterRef})
    prefs = convert(JuMP.Containers.SparseAxisArray, prefs)
    groups = _groups(prefs)
    length(unique(groups)) != 1 && error("Array contains parameters from multiple groups.")
    return groups[1]
end

"""
    is_independent(pref::ParameterRef)::Bool
Returns true for `pref` if it is independent or false otherwise.
"""
function is_independent(pref::ParameterRef)::Bool
    return JuMP.owner_model(pref).params[JuMP.index(pref)].independent
end

"""
    parameter_by_name(model::InfiniteModel, name::String)
Return the parameter reference assoociated with a parameter.
"""
function parameter_by_name(model::InfiniteModel, name::String)
    if model.name_to_param === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_param = Dict{String, Int}()
        for (param, param_name) in model.param_to_name
            if haskey(model.name_to_param, param_name)
                # -1 is a special value that means this string does not map to
                # a unique variable name.
                model.name_to_param[param_name] = -1
            else
                model.name_to_param[param_name] = param
            end
        end
    end
    index = get(model.name_to_param, name, nothing)
    if index isa Nothing
        return nothing
    elseif index == -1
        error("Multiple parameters have the name $name.")
    else
        return ParameterRef(model, index)
    end
    return
end

"""
    all_parameters(model::InfiniteModel)
Return all of the infinite parameters as a vector of type `ParameterRef`.
"""
function all_parameters(model::InfiniteModel)
    pref_list = Vector{ParameterRef}(undef, num_parameters(model))
    indexes = sort([index for index in keys(model.params)])
    counter = 1
    for index in indexes
        pref_list[counter] = ParameterRef(model, index)
        counter += 1
    end
    return pref_list
end

# Define functions to extract the names of parameters
function _names(arr::JuMP.Containers.SparseAxisArray{<:ParameterRef})
    return [JuMP.name(arr[k]) for k in keys(arr.data)]
end

function _root_name(pref::ParameterRef)
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

function _root_names(prefs::Tuple)
    root_names = Vector{String}(undef, length(prefs))
    for i = 1:length(root_names)
        if isa(prefs[i], ParameterRef)
            root_names[i] = _root_name(prefs[i])
        else
            root_names[i] = _root_name(collect(values(prefs[i].data))[1])
        end
    end
    return root_names
end

function _only_one_name(arr::JuMP.Containers.SparseAxisArray{<:ParameterRef})
    names = _names(arr)
    root_names = Vector{String}(undef, length(names))
    for i = 1:length(root_names)
        first_bracket = findfirst(isequal('['), names[i])
        root_names[i] = names[i][1:first_bracket-1]
    end
    return length(unique(root_names)) == 1
end

_only_one_name(pref::ParameterRef) = true

# Internal functions for group checking
function _groups(arr::JuMP.Containers.SparseAxisArray{<:ParameterRef})
    return [group_id(arr[k]) for k in keys(arr.data)]
end

function _groups(prefs::Tuple)
    groups = Vector{Int}(undef, length(prefs))
    for i = 1:length(groups)
        if isa(prefs[i], ParameterRef)
            groups[i] = group_id(prefs[i])
        else
            groups[i] = group_id(prefs[i])
        end
    end
    return groups
end

function _only_one_group(arr::JuMP.Containers.SparseAxisArray{<:ParameterRef})
    groups = _groups(arr)
    return length(unique(groups)) == 1
end

_only_one_group(pref::ParameterRef) = true

# Internal function for extracting parameter references
function _list_parameter_refs(prefs::Tuple)
    list = ParameterRef[]
    for pref in prefs
        if isa(pref, ParameterRef)
            push!(list, pref)
        else
            for k in keys(pref.data)
                push!(list, pref[k])
            end
        end
    end
    return list
end
