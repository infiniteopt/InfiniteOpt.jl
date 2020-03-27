"""
    VectorTuple{T}


"""
struct VectorTuple{T}
    values::Vector{T}
    starts::Vector{Int64}
    ends::Vector{Int64}
    indices::Vector
end

## Parse the proper indices for each container type
# Array
_get_indices(arr::Array) = CartesianIndices(arr)

# DenseAxisArray
_get_indices(arr::JuMPC.DenseAxisArray) = (CartesianIndices(arr), axes(arr))

# SparseAxisArray
_get_indices(arr::JuMPC.SparseAxisArray) = sort(collect(keys(arr)))

# Fallback
_get_indices(arr) = nothing

## Order arrays by index if necessary (to take care of SparseAxisArrays)
# SparseAxisArray
function _make_ordered(arr::JuMPC.SparseAxisArray, indices::Vector)
    return [arr[i] for i in indices]
end

# other
_make_ordered(arr, indices) = arr

# Constructor from Tuple
function VectorTuple(tuple::Tuple)::VectorTuple
    num_elements = length(tuple)
    element_lengths = [length(tuple[i]) for i in eachindex(tuple)]
    starts = Vector{Int64}(undef, num_elements)
    ends = Vector{Int64}(undef, num_elements)
    for i in eachindex(starts)
        if i == 1
            starts[i] = 1
            ends[i] = element_lengths[i]
        else
            starts[i] = starts[i - 1] + element_lengths[i - 1]
            ends[i] = starts[i] + element_lengths[i] - 1
        end
    end
    indices = [_get_indices(tuple[i]) for i in eachindex(tuple)]
    if any(k isa JuMPC.SparseAxisArray for k in tuple)
        tuple = Tuple(_make_ordered(tuple[i], indices[i]) for i in eachindex(tuple))
    end
    return VectorTuple([(tuple...)...], starts, ends, indices)
end

# Default Constructor
function VectorTuple()::VectorTuple
    return VectorTuple(Any[], Int64[], Int64[], Any[])
end

# Determine the number of unqiue tuple elements
function tuple_length(vt::VectorTuple)::Int
    return length(vt.starts)
end

## Extend Base.getindex
# linear index along vector
Base.getindex(vt::VectorTuple, i) = Base.getindex(vt.values, i)

# Tuple indexing
function Base.getindex(vt::VectorTuple, i, j)
    return Base.getindex(vt.values[vt.starts[i]:vt.ends[i]], j)
end

## Extend Base.setindex!
# linear index along vector
function Base.setindex!(vt::VectorTuple, v, i)
    return Base.setindex!(vt.values, v, i)
end

# Tuple indexing
function Base.setindex!(vt::VectorTuple, v, i, j) # TODO make work with colons
    return Base.setindex!(vt.values, v, vt.starts[i] + j - 1)
end

# Extend Base.length
Base.length(vt::VectorTuple) = Base.length(vt.values)

# Extend Base.isempty
Base.isempty(vt::VectorTuple) = Base.isempty(vt.values)

# Extend Base.empty!
function Base.empty!(vt::VectorTuple)
    empty!(vt.values)
    empty!(vt.starts)
    empty!(vt.ends)
    empty!(vt.indices)
    return vt
end

# Extend Base.in
function Base.in(item::T, vt::VectorTuple{T}) where {T}
    return Base.in(item, vt.values)
end

# Extend Base.eachindex
Base.eachindex(vt::VectorTuple) = Base.eachindex(vt::VectorTuple)

# Extend Base.keys
Base.keys(vt::VectorTuple) = Base.keys(vt::VectorTuple)

## Extend Base.push!
# Single item
function Base.push!(vt::VectorTuple{T}, item::T) where {T}
    push!(vt.values, item)
    push!(vt.starts, vt.ends[end] + 1)
    push!(vt.ends, vt.ends[end] + 1)
    push!(vt.indices, nothing)
    return vt
end

# AbstractArray
function Base.push!(vt::VectorTuple{T}, items::AbstractArray{<:T}) where {T}
    push!(vt.values, items...)
    push!(vt.starts, vt.ends[end] + 1)
    push!(vt.ends, length(vt.values))
    push!(vt.indices, _get_indices(items))
    return vt
end

## Update the indices if an element is deleted
# SparseAxisArray
function _update_indices(inds::Vector{<:Tuple}, i::Int)
    return deleteat!(inds, i)
end

# DenseAxisArray
function _update_indices(inds::Tuple, i::Int)
    new_inds = [k.I for k in JuMPC.DenseAxisArrayKeys(Base.Iterators.product(inds[2]...))]
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Array
function _update_indices(inds::CartesianIndices, i::Int)
    new_inds = [Tuple(inds[k]) for k in 1:length(inds)]
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Extend Base.deleteat!
function Base.deleteat!(vt::VectorTuple, i::Int)
    deleteat!(vt.values, i)
    group = findfirst(vt.starts >= i)
    if vt.starts[group] == vt.ends[group]
        deleteat!(vt.starts, group)
        deleteat!(vt.ends, group)
        deleteat!(vt.indices, group)
    else
        vt.indices[group] = _update_indices(vt.indices[group],
                                            i - vt.starts[group] + 1)
    end
    vt.starts[vt.starts > i] .-= 1
    vt.ends[vt.ends >= i] .-= 1
    return vt
end

# Extend Base.findfirst
Base.findfirst(pred::Function, vt::VectorTuple) = Base.findfirst(pred, vt.values)

## Make the original array for a particular tuple element
# Array
function _make_array(values::Vector, first::Int, last::Int,
                     indices::CartesianIndices)
    return [values[first:last][i] for i in LinearIndices(indices)]
end

# DenseAxisArray
function _make_array(values::Vector, first::Int, last::Int, indices::Tuple)
    array = _make_array(values, first, last, indices[1])
    return JuMPC.DenseAxisArray(array, indices[2]...)
end

# SparseAxisArray
function _make_array(values::Vector, first::Int, last::Int, indices::Vector)
    data = Dict(indices[i] => values[first:last][i] for i in 1:length(indices))
    return JuMPC.SparseAxisArray(data)
end

# Other
function _make_array(values::Vector, first::Int, last::Int, indices::Nothing)
    return values[first:last]
end

# Define Tuple conversion
function Base.Tuple(vt::VectorTuple; use_indices = true)::Tuple
    if use_indices
        return Tuple(_make_array(vt.values, vt.starts[i], vt.ends[i],
                                 vt.indices[i]) for i = 1:tuple_length(vt))
    else
        return Tuple(vt.values[vt.starts[i]:vt.ends[i]] for i = 1:tuple_length(vt))
    end
end

## Enable pretty printing
# Extend Base.string
function Base.string(vt::VectorTuple)::String
    return string(Tuple(vt, use_indices = false))
end

# Show VectorTuples in REPLMode
function Base.show(io::IO, vt::VectorTuple)
    print(io, string(vt))
    return
end

# Show VectorTuples in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", vt::VectorTuple)
    print(io, string(vt))
end
