"""
    VectorTuple{T}

A collection DataType for storing a `Tuple` of singular elements of type `T`
and/or `AbstractArray{<:T}`s in a convenient vectorized form that utilizes
linear indexing. Moreover `VectorTuple`s can be modified using standard vector
operations such as `empty!`, `push!`, and `deleteat!`. `VectorTuple`s should
be defined from an original tuple via `VectorTuple(tuple)` or by listing the
tuple elements `VectorTuple(items...)`. Note this is still an experimental
type and is primarily intended to store infinite parameter reference tuples and
point variable support value tuples. Some of the notable capabilities are
exemplified below.

**Example**
```julia-repl
julia> tuple = (3, [-2, 4], ones(2, 2))
(3, [-2, 4], [1.0 1.0; 1.0 1.0])

julia> vt = VectorTuple(tuple) # make by listing items (notice everything is a vector)
([3.0], [3.0, 4.0], [1.0, 1.0, 1.0, 1.0])

julia> vt[2] # linear indexing
-2.0

julia> vt[2, 2] # tuple indexing (note the second index is treated linearly)
4.0

julia> vt[6:end] # linear slicing
2-element Array{Float64,1}:
 1.0
 1.0

julia> vt[2:3, :] # tuple slicing
2-element Array{Array{Float64,1},1}:
 [-2.0, 4.0]
 [1.0, 1.0, 1.0, 1.0]

julia> tuple2 = Tuple(vt) # rebuild original Tuple with original indices
([3.0], [-2.0, 4.0], [1.0 1.0; 1.0 1.0])

julia> push!(vt, [42., 42]) # add new tuple element
([3.0], [-2.0, 4.0], [1.0, 1.0, 1.0, 1.0], [42.0, 42.0])

julia> deleteat!(vt, 4) # delete an element via linear indexing
([3.0], [-2.0, 4.0], [1.0, 1.0, 1.0])

julia> Tuple(vt) # The 3rd element becomes a SparseAxisArray because of deletion
([3.0], [-2.0, 4.0],   [1, 2]  =  1.0
  [2, 2]  =  1.0
  [2, 1]  =  1.0)

julia> deleteat!(vt, 3, tuple_index = true) # delete a whole tuple element
([3.0], [-2.0, 4.0])
```
"""
struct VectorTuple{T}
    values::Vector{T}
    starts::Vector{Int64}
    ends::Vector{Int64}
    indices::Vector
end

## Parse the proper indices for each container type
# Array
_get_indices(arr::Array)::CartesianIndices = CartesianIndices(arr)

# DenseAxisArray
_get_indices(arr::JuMPC.DenseAxisArray)::Tuple = (CartesianIndices(arr), axes(arr))

# SparseAxisArray
_get_indices(arr::JuMPC.SparseAxisArray)::Vector = sort(collect(keys(arr)))

# Fallback
_get_indices(arr) = nothing

## Order arrays by index if necessary (to take care of SparseAxisArrays)
# SparseAxisArray
function _make_ordered(arr::JuMPC.SparseAxisArray, indices::Vector)::Vector
    return [arr[i] for i in indices]
end

# other
_make_ordered(arr, indices) = arr

# Constructor from a Tuple
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

# Constructor from various arguments (like splatting the tuple)
function VectorTuple(items...)::VectorTuple
    return VectorTuple(items)
end

# Default Constructor
function VectorTuple()::VectorTuple
    return VectorTuple(Any[], Int64[], Int64[], Any[])
end

# Extend Base.:(==)
function Base.:(==)(vt1::VectorTuple, vt2::VectorTuple)::Bool
    return all(getproperty(vt1, f) == getproperty(vt2, f) for f in fieldnames(VectorTuple))
end

# Determine the number of unqiue tuple elements
function tuple_length(vt::VectorTuple)::Int
    return length(vt.starts)
end

## Extend first and last indices
# Base.firstindex (linear)
Base.firstindex(vt::VectorTuple)::Int = Base.firstindex(vt.values)

# Base.firstindex (tuple)
Base.firstindex(vt::VectorTuple, d::Int)::Int = 1

# Base.lastindex (linear)
Base.lastindex(vt::VectorTuple)::Int = Base.lastindex(vt.values)

# Base.lastindex (tuple)
function Base.lastindex(vt::VectorTuple, d::Int)::Int
    return d == 1 ? tuple_length(vt) : error("lastindex(::VectorTuple, " *
                                             "d::Int > 1) not defined.")
end

## Extend Base.getindex
# linear index along vector
Base.getindex(vt::VectorTuple, i) = Base.getindex(vt.values, i)

# Tuple indexing (single position)
function Base.getindex(vt::VectorTuple, i::Int, j)
    return Base.getindex(vt.values[vt.starts[i]:vt.ends[i]], j)
end

# Tuple indexing (with a colon)
function Base.getindex(vt::VectorTuple, i::Colon, j)
    return [Base.getindex(vt, i, j) for i = 1:tuple_length(vt)]
end

# Tuple indexing (with an iterable) --> assumed backup
function Base.getindex(vt::VectorTuple, is, j)
    return [Base.getindex(vt, i, j) for i in is]
end

## Extend Base.setindex!
# linear index along vector
function Base.setindex!(vt::VectorTuple, v, i)
    return Base.setindex!(vt.values, v, i)
end

# Tuple indexing (with integer i and colon j)
function Base.setindex!(vt::VectorTuple, v, i::Int, j::Colon)
    linear_index = vt.starts[i]:vt.ends[i]
    return Base.setindex!(vt.values, v, linear_index)
end

# Tuple indexing (with integer i and any non colon j)
function Base.setindex!(vt::VectorTuple, v, i::Int, j)
    linear_index = (vt.starts[i] - 1) .+ j
    if any(linear_index .> vt.ends[i])
        throw(BoundsError(vt, [i, j]))
    end
    return Base.setindex!(vt.values, v, linear_index)
end

# Tuple indexing (with colon i and integer j)
function Base.setindex!(vt::VectorTuple, v, i::Colon, j::Int)
    linear_index = (vt.starts[i] .- 1) .+ j
    return Base.setindex!(vt.values, v, linear_index)
end

# Tuple indexing (with colon i and noninteger j)
function Base.setindex!(vt::VectorTuple, v, i::Colon, j)
    return Base.setindex!(vt, v, 1:tuple_length(vt), j)
end

# Tuple indexing (with array i and any j)
function Base.setindex!(vt::VectorTuple, v, i, j)
    @assert length(v) == length(i)
    return [Base.setindex!(vt, v[k], i[k], j) for k in 1:length(i)]
end

# Extend Base.length
Base.length(vt::VectorTuple)::Int = Base.length(vt.values)

# Extend Base.isempty
Base.isempty(vt::VectorTuple)::Bool = Base.isempty(vt.values)

# Extend Base.eachindex
Base.eachindex(vt::VectorTuple) = Base.eachindex(vt.values)

# Extend Base.keys
Base.keys(vt::VectorTuple) = Base.keys(vt.values)

# Extend Base.findfirst
Base.findfirst(pred::Function, vt::VectorTuple) = Base.findfirst(pred, vt.values)

# Extend Base.findall
Base.findall(pred::Function, vt::VectorTuple) = Base.findall(pred, vt.values)

# Extend Base.in
function Base.in(item, vt::VectorTuple)::Bool
    return Base.in(item, vt.values)
end

# Extend Base.empty!
function Base.empty!(vt::VectorTuple)::VectorTuple
    empty!(vt.values)
    empty!(vt.starts)
    empty!(vt.ends)
    empty!(vt.indices)
    return vt
end

## Extend Base.push!
# Single item
function Base.push!(vt::VectorTuple{T}, item::T)::VectorTuple where {T}
    push!(vt.values, item)
    push!(vt.starts, vt.ends[end] + 1)
    push!(vt.ends, vt.ends[end] + 1)
    push!(vt.indices, nothing)
    return vt
end

# AbstractArray
function Base.push!(vt::VectorTuple{T},
                    items::AbstractArray{<:T})::VectorTuple where {T}
    indices = _get_indices(items)
    push!(vt.values, _make_ordered(items, indices)...)
    push!(vt.starts, vt.ends[end] + 1)
    push!(vt.ends, length(vt.values))
    push!(vt.indices, indices)
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
function Base.deleteat!(vt::VectorTuple, i::Int; tuple_index = false)
    if tuple_index
        len = length(vt.starts[i]:vt.ends[i])
        start = vt.starts[i]
        deleteat!(vt.values, vt.starts[i]:vt.ends[i])
        deleteat!(vt.starts, i)
        deleteat!(vt.ends, i)
        deleteat!(vt.indices, i)
        vt.starts[vt.starts .> start] .-= len
        vt.ends[vt.ends .>= start] .-= len
    else
        group = findfirst(vt.ends .>= i)
        deleteat!(vt.values, i)
        if vt.starts[group] == vt.ends[group]
            deleteat!(vt.starts, group)
            deleteat!(vt.ends, group)
            deleteat!(vt.indices, group)
        else
            vt.indices[group] = _update_indices(vt.indices[group],
                                                i - vt.starts[group] + 1)
        end
        vt.starts[vt.starts .> i] .-= 1
        vt.ends[vt.ends .>= i] .-= 1
    end
    return vt
end

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
    return first == last ? values[first] : values[first:last]
end

# Define Tuple construction from a VectorTuple
function Base.Tuple(vt::VectorTuple; use_indices = true)::Tuple
    if use_indices
        return Tuple(_make_array(vt.values, vt.starts[i], vt.ends[i],
                                 vt.indices[i]) for i = 1:tuple_length(vt))
    else
        return Tuple(vt.starts[i] == vt.ends[i] ? vt.values[vt.starts[i]] : vt.values[vt.starts[i]:vt.ends[i]] for i = 1:tuple_length(vt))
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
