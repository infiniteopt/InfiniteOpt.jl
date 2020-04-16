################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################


################################################################################
#                             PARAMETER DEFINITION
################################################################################
function build_parameter(_error::Function, set::InfiniteArraySet;
                         label::Set{Symbol} = Set{Symbol}(),
                         num_supports::Int = 0,
                         supports::Union{Number, Vector{<:Number}} = Number[],
                         sig_fig::Int = 5,
                         extra_kw_args...)::DependentParameters
    return nothing
end

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
