# Define symbol input
const Parameter = :Parameter

# Extend Base.copy for new variable types
Base.copy(v::ParameterRef, new_model::InfiniteModel) = ParameterRef(new_model, v.index)

"""
    JuMP.add_variable(model::InfiniteModel, v::InfOptParameter, name::String="")
Extend the `JuMP.add_variable` function to accomodate our infinite parameters.
"""
function JuMP.add_variable(model::InfiniteModel, v::InfOptParameter, name::String="")
    model.next_param_index += 1
    pref = ParameterRef(model, model.next_param_index)
    model.params[pref.index] = v
    JuMP.set_name(pref, name)
    return pref
end

"""
    JuMP.delete(model::InfiniteModel, pref::ParameterRef)
Extend the `JuMP.delete` function to accomodate infinite parameters
"""
function JuMP.delete(model::InfiniteModel, pref::ParameterRef)
    @assert JuMP.is_valid(model, pref)
    delete!(model.params, JuMP.index(pref))
    delete!(model.param_to_name, JuMP.index(pref))
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
    return
end

# Define functions to extract the names of parameters
function _get_names(arr::AbstractArray{<:ParameterRef})
    if isa(arr, JuMP.Containers.SparseAxisArray)
        return [JuMP.name(arr[k]) for k in keys(arr.data)]
    else
        return [JuMP.name(arr[k]) for k in CartesianIndices(arr)]
    end
end

function _get_names(arr::Array{<:AbstractArray{<:ParameterRef}})
    names = String[]
    for a in arr
        if isa(a, JuMP.Containers.SparseAxisArray)
            names = [names; [JuMP.name(arr[k]) for k in keys(arr.data)]]
        else
            names = [names; [JuMP.name(arr[k]) for k in CartesianIndices(arr)]]
        end
    end
    return names
end

function _get_root_names(arr::Union{AbstractArray{<:ParameterRef}, Array{<:AbstractArray{<:ParameterRef}}})
    names = _get_names(arr)
    root_names = Vector{String}(undef, length(names))
    for i = 1:length(names)
        first_bracket = findfirst(isequal('['), names[i])
        if first_bracket == nothing
            root_names[i] = names[i]
        else
            root_names[i] = names[i][1:first_bracket-1]
        end
    end
    return unique(root_names)
end
