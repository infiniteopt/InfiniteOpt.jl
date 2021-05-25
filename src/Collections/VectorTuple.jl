"""
    VectorTuple{T}

A collection DataType for storing a `Tuple` of singular elements of type `T`
and/or `AbstractArray{<:T}`s in a convenient vector form that utilizes
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
    ranges::Vector{UnitRange{Int}}
    indices::Vector{Any}
end

## Parse the proper indices for each container type
# Array
function _get_indices(
    arr::Array{T, N}
    )::CartesianIndices{N, NTuple{N, Base.OneTo{Int}}} where {T, N}
    return CartesianIndices(arr)
end

# DenseAxisArray
function _get_indices(arr::JuMPC.DenseAxisArray{T, N, Ax}
    )::Tuple{CartesianIndices{N, NTuple{N, Base.OneTo{Int}}}, Ax} where {T, N, Ax}
    return CartesianIndices(arr), axes(arr)
end

# SparseAxisArray
function _get_indices(arr::JuMPC.SparseAxisArray{T, N, K})::Vector{K} where {T, N, K}
    return sort!(collect(InfiniteOpt._keys(arr)))
end

# Fallback
_get_indices(arr)::Nothing = nothing

## Order arrays by index if necessary (to take care of SparseAxisArrays)
# SparseAxisArray
function _make_ordered(arr::JuMPC.SparseAxisArray{T},
                       indices::Vector)::Vector{T} where {T}
    return [arr[i] for i in indices]
end

# other
_make_ordered(arr, indices) = arr

# internal function for constructor
function _make_vt_components(tuple::Tuple)
    # form the range vector
    ranges = Vector{UnitRange{Int}}(undef, length(tuple))
    for i in eachindex(ranges)
        element_length = length(tuple[i])
        if i === 1
            ranges[i] = UnitRange(1, element_length)
        else
            ranges[i] = UnitRange(ranges[i-1].stop + 1, ranges[i-1].stop + element_length)
        end
    end
    # get the indices and reorder the tuple if needed
    indices = Any[_get_indices(tuple[i]) for i in eachindex(tuple)]
    if any(k isa JuMPC.SparseAxisArray for k in tuple)
        tuple = Tuple(_make_ordered(tuple[i], indices[i]) for i in eachindex(tuple))
    end
    # expand the elements into a single vector
    vals = [j for i in tuple for j in i]
    return vals, ranges, indices
end

# Constructor from a Tuple with type
function VectorTuple{T}(tuple::Tuple)::VectorTuple{T} where {T}
    if isempty(tuple)
        return VectorTuple{T}()
    end
    vals, ranges, indices = _make_vt_components(tuple)
    return VectorTuple{T}(vals, ranges, indices)
end

# Constructor from a Tuple without type
function VectorTuple(tuple::Tuple)::VectorTuple
    if isempty(tuple)
        return VectorTuple()
    end
    vals, ranges, indices = _make_vt_components(tuple)
    return VectorTuple(vals, ranges, indices)
end

# Constructor from various arguments (like splatting the tuple) with type
function VectorTuple{T}(items...)::VectorTuple{T} where {T}
    return VectorTuple{T}(items)
end

# Constructor from various arguments (like splatting the tuple) without type
function VectorTuple(items...)::VectorTuple
    return VectorTuple(items)
end

# Default Constructor with type
function VectorTuple{T}()::VectorTuple{T} where {T}
    return VectorTuple(T[], UnitRange{Int}[], Any[])
end

# Default Constructor without type
function VectorTuple()::VectorTuple
    return VectorTuple(Any[], UnitRange{Int}[], Any[])
end

# Extend Base.:(==)
function Base.:(==)(vt1::VectorTuple, vt2::VectorTuple)::Bool
    return all(getproperty(vt1, f) == getproperty(vt2, f) for f in fieldnames(VectorTuple))
end

# Extend Base.size
function Base.size(vt::VectorTuple, d::Int)::Int
    d > 1 && throw(ArgumentError("invalid VectorTuple dimension $d"))
    return length(vt.ranges)
end

# Determine if 2 VectorTuples have the same structure
function same_structure(vt1::VectorTuple, vt2::VectorTuple;
                        use_indices = true)::Bool
    same_size = vt1.ranges == vt2.ranges
    return same_size && (!use_indices || (vt1.indices == vt2.indices))
end

# Get the element type
param_type(vt::VectorTuple{T}) where {T} = T

# Extend Base.copy
function Base.copy(vt::VectorTuple{T})::VectorTuple{T} where {T}
    return VectorTuple(copy(vt.values), copy(vt.ranges), copy(vt.indices))
end

## Extend first and last indices
# Base.firstindex (linear)
Base.firstindex(vt::VectorTuple)::Int = firstindex(vt.values)

# Base.firstindex (tuple)
Base.firstindex(vt::VectorTuple, d::Int)::Int = 1

# Base.lastindex (linear)
Base.lastindex(vt::VectorTuple)::Int = lastindex(vt.values)

# Base.lastindex (tuple)
function Base.lastindex(vt::VectorTuple, d::Int)::Int
    return d == 1 ? size(vt, 1) : throw(ArgumentError("lastindex(::VectorTuple, " *
                                                      "d::Int > 1) not defined."))
end

## Extend Base.getindex
# linear index along vector
Base.getindex(vt::VectorTuple, i) = getindex(vt.values, i)

# Tuple indexing (single position)
function Base.getindex(vt::VectorTuple, i::Int, j)
    return getindex(vt.values[vt.ranges[i]], j)
end

# Tuple indexing (with a colon)
function Base.getindex(vt::VectorTuple, i::Colon, j)
    return [getindex(vt, i, j) for i = 1:size(vt, 1)]
end

# Tuple indexing (with an iterable) --> assumed backup
function Base.getindex(vt::VectorTuple, is, j)
    return [getindex(vt, i, j) for i in is]
end

## Extend Base.setindex!
# linear index along vector
function Base.setindex!(vt::VectorTuple, v, i)
    return setindex!(vt.values, v, i)
end

# Tuple indexing (with integer i and colon j)
function Base.setindex!(vt::VectorTuple, v, i::Int, j::Colon)
    return setindex!(vt.values, v, vt.ranges[i])
end

# Tuple indexing (with integer i and any non colon j)
function Base.setindex!(vt::VectorTuple, v, i::Int, j)
    linear_index = (vt.ranges[i].start - 1) .+ j
    if any(linear_index .> vt.ranges[i].stop)
        throw(BoundsError(vt, [i, j]))
    end
    return setindex!(vt.values, v, linear_index)
end

# Tuple indexing (with colon i and integer j)
function Base.setindex!(vt::VectorTuple, v, i::Colon, j::Int)
    linear_indexes = [vt.ranges[k].start - 1 + j for k in eachindex(vt.ranges)]
    return setindex!(vt.values, v, linear_indexes)
end

# Tuple indexing (with colon i and noninteger j)
function Base.setindex!(vt::VectorTuple, v, i::Colon, j)
    return setindex!(vt, v, 1:size(vt, 1), j)
end

# Tuple indexing (with array i and any j)
function Base.setindex!(vt::VectorTuple, v, i, j)
    @assert length(v) == length(i)
    return [setindex!(vt, v[k], i[k], j) for k in 1:length(i)]
end

# Extend simple Base redirect methods with VectorTuple arguments
for op = (:length, :isempty, :eachindex, :keys)
    @eval Base.$op(vt::VectorTuple) = $op(vt.values)
end

# Extend simple Base 2 argument methods
for op = (:findfirst, :findall)
    @eval Base.$op(f::Function, vt::VectorTuple) = $op(f, vt.values)
end

# Extend Base.in
Base.in(i, vt::VectorTuple)::Bool = in(i, vt.values)

# Extend Base.iterate
Base.iterate(vt::VectorTuple) = iterate(vt.values)
Base.iterate(vt::VectorTuple, i) = iterate(vt.values, i)

# Extend Base.empty!
function Base.empty!(vt::VectorTuple{T})::VectorTuple{T} where {T}
    empty!(vt.values)
    empty!(vt.ranges)
    empty!(vt.indices)
    return vt
end

## Extend Base.push!
# Single item
function Base.push!(vt::VectorTuple{T}, item::T)::VectorTuple{T} where {T}
    push!(vt.values, item)
    push!(vt.ranges, UnitRange(vt.ranges[end].stop + 1, vt.ranges[end].stop + 1))
    push!(vt.indices, nothing)
    return vt
end

# AbstractArray
function Base.push!(vt::VectorTuple{T},
                    items::AbstractArray{<:T})::VectorTuple{T} where {T}
    indices = _get_indices(items)
    push!(vt.values, _make_ordered(items, indices)...)
    push!(vt.ranges, UnitRange(vt.ranges[end].stop + 1, length(vt.values)))
    push!(vt.indices, indices)
    return vt
end

## Update the indices if an element is deleted
# SparseAxisArray
function _update_indices(inds::Vector{<:Tuple}, i)
    new_inds = deleteat!(inds, i)
    return length(new_inds) >= 2 ? new_inds : nothing
end

# DenseAxisArray
function _update_indices(inds::Tuple, i)
    new_inds = reduce(vcat, Base.Iterators.product(inds[2]...))
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Array
function _update_indices(inds::CartesianIndices, i)
    new_inds = [Tuple(inds[k]) for k in 1:length(inds)]
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Extend Base.deleteat!
function Base.deleteat!(vt::VectorTuple{T}, i::Int;
                        tuple_index::Bool = false)::VectorTuple{T} where {T}
    if tuple_index
        len = length(vt.ranges[i])
        start = vt.ranges[i].start
        deleteat!(vt.values, vt.ranges[i])
        deleteat!(vt.ranges, i)
        deleteat!(vt.indices, i)
        for j in eachindex(vt.ranges)
            if vt.ranges[j].start > start
                vt.ranges[j] = vt.ranges[j] .- len
            end
        end
    else
        group = findfirst(x->in(i, x), vt.ranges)
        deleteat!(vt.values, i)
        if length(vt.ranges[group]) == 1
            deleteat!(vt.ranges, group)
            deleteat!(vt.indices, group)
        else
            vt.indices[group] = _update_indices(vt.indices[group],
                                                i - vt.ranges[group].start + 1)
        end
        for j in eachindex(vt.ranges)
            if vt.ranges[j].start > i
                vt.ranges[j] = vt.ranges[j] .- 1
            elseif vt.ranges[j].stop >= i
                vt.ranges[j] = UnitRange(vt.ranges[j].start, vt.ranges[j].stop - 1)
            end
        end
    end
    return vt
end

# Extend Base.deleteat! (boolean list)
function Base.deleteat!(vt::VectorTuple{T}, inds::AbstractVector{<:Bool};
                        tuple_index::Bool = false)::VectorTuple{T} where {T}
    if tuple_index
        deleteat!(vt.values, [j for i in vt.ranges[inds] for j in i])
        deleteat!(vt.indices, inds)
        prev_sum = 0
        for i in eachindex(vt.ranges)
            if inds[i]
                prev_sum += length(vt.ranges[i])
            else
                vt.ranges[i] = vt.ranges[i] .- prev_sum
            end
        end
        deleteat!(vt.ranges, inds)
    else
        deleteat!(vt.values, inds)
        prev_sum = 0
        delete_inds = []
        for i in eachindex(vt.ranges)
            delete_sum = sum(inds[vt.ranges[i]])
            if length(vt.ranges[i]) == delete_sum
                push!(delete_inds, i)
            else
                group_inds = vt.ranges[i][inds[vt.ranges[i]]] .+ (1 - vt.ranges[i].start)
                if length(group_inds) != 0
                    vt.indices[i] = _update_indices(vt.indices[i], group_inds)
                end
                vt.ranges[i] = UnitRange(vt.ranges[i].start - prev_sum,
                                         vt.ranges[i].stop - prev_sum - delete_sum)
            end
            prev_sum += delete_sum
        end
        deleteat!(vt.ranges, delete_inds)
        deleteat!(vt.indices, delete_inds)
    end
    return vt
end

# TODO add deleteat! for sequential indices

## Make the original array for a particular tuple element
# Array
function _make_array(values::Vector{T},
                     indices::CartesianIndices{N})::Array{T, N} where {T, N}
    return [values[i] for i in LinearIndices(indices)]
end

# DenseAxisArray
function _make_array(values::Vector, indices::Tuple)
    array = _make_array(values, indices[1])
    return JuMPC.DenseAxisArray(array, indices[2]...)
end

# SparseAxisArray
function _make_array(values::Vector, indices::Vector)
    data = Dict(indices[i] => values[i] for i in eachindex(indices))
    return JuMPC.SparseAxisArray(data)
end

# Other
function _make_array(values::Vector, indices::Nothing)
    return length(values) == 1 ? first(values) : values
end

# Define Tuple construction from a VectorTuple
function Base.Tuple(vt::VectorTuple; use_indices = true)::Tuple
    if use_indices
        return Tuple(_make_array(vt.values[vt.ranges[i]], vt.indices[i])
                     for i in eachindex(vt.ranges))
    else
        return Tuple(length(vt.ranges[i]) == 1 ? vt.values[vt.ranges[i].start] : vt.values[vt.ranges[i]]
                     for i in eachindex(vt.ranges))
    end
end

# Define Tuple using construction from a values vector in combination with a VT
function Base.Tuple(values::Vector, vt::VectorTuple; use_indices = true)::Tuple
    @assert length(values) == length(vt)
    if use_indices
        return Tuple(_make_array(values[vt.ranges[i]], vt.indices[i])
                     for i in eachindex(vt.ranges))
    else
        return Tuple(length(vt.ranges[i]) == 1 ? values[vt.ranges[i].start] : values[vt.ranges[i]]
                     for i in eachindex(vt.ranges))
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
