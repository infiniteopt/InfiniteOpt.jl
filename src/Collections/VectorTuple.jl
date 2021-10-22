"""
    VectorTuple{T, I <: Tuple}

A collection DataType for storing a `Tuple` of singular elements of type `T`
and/or `AbstractArray{<:T}`s in a convenient vector form that utilizes
linear indexing. Here `I` is denotes the type of a `Tuple` that stores the 
indices of each tuple element as given by `indices`. `VectorTuple`s should be 
defined from an original tuple via `VectorTuple(tuple)` or by listing the tuple 
elements `VectorTuple(items...)`. Note this is still an experimental type and is 
primarily intended to store infinite parameter reference tuples and point 
variable support value tuples. Some of the notable capabilities are exemplified 
below.

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

julia> inds = [true, true, true, false, true, true, true];

julia> restricted_copy(vt, delete_locs) # make a copy with deleted elements
([3.0], [-2.0, 4.0], [1.0, 1.0, 1.0])

julia> Tuple(vt) # The 3rd element becomes a SparseAxisArray because of deletion
([3.0], [-2.0, 4.0],   [1, 2]  =  1.0
  [2, 2]  =  1.0
  [2, 1]  =  1.0)
```
"""
struct VectorTuple{T, I <: Tuple}
    values::Vector{T}
    ranges::Vector{UnitRange{Int}}
    indices::I
end

## Order arrays by index if necessary (to take care of SparseAxisArrays)
# SparseAxisArray
function _make_ordered(arr::JuMPC.SparseAxisArray{T}, inds)::Vector{T} where {T}
    return vectorize(arr, inds)
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
    # get the indices
    inds = Tuple(indices(t) for t in tuple)
    # prepare the values vector
    vals = [j for i in eachindex(tuple) for j in _make_ordered(tuple[i], inds[i])]
    # return everything
    return vals, ranges, inds
end

# Constructor from a Tuple without type
function VectorTuple(tuple::Tuple)
    vals, ranges, indices = _make_vt_components(tuple)
    return VectorTuple(vals, ranges, indices)
end

# Constructor from various arguments (like splatting the tuple) without type
function VectorTuple(items...)
    return VectorTuple(items)
end

# Extend Base.isequal
function Base.isequal(vt1::VectorTuple, vt2::VectorTuple)::Bool
    return all(isequal(getproperty(vt1, f), getproperty(vt2, f)) for f in fieldnames(VectorTuple))
end

# Extend Base.:(==)
function Base.:(==)(vt1::VectorTuple, vt2::VectorTuple)::Bool
    return isequal(vt1, vt2)
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

# Extend Base.copy
function Base.copy(vt::VectorTuple)
    return VectorTuple(copy(vt.values), copy(vt.ranges), vt.indices)
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

## Update the indices if an element is deleted
# SparseAxisArray
function _update_indices(inds::ContainerIndices{1, <:Vector}, i)
    new_axes = deleteat!(copy(inds.axes), i)
    cart_inds = CartesianIndices(new_axes)
    return length(new_axes) >= 2 ? ContainerIndices(cart_inds, new_axes) : nothing
end

# DenseAxisArray
function _update_indices(inds::ContainerIndices, i)
    new_axes = reduce(vcat, Base.Iterators.product(inds.axes...))
    new_inds = ContainerIndices(CartesianIndices(new_axes), new_axes)
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Array
function _update_indices(inds::ContainerIndices{N, Nothing}, i) where {N}
    new_axes = [Tuple(inds.indices[k]) for k in 1:length(inds.indices)]
    new_inds = ContainerIndices(CartesianIndices(new_axes), new_axes)
    return _update_indices(new_inds, i) # call the SparseAxisArray version
end

# Make a restricted copy (a copy with deleted elements)
function restricted_copy(vt::VectorTuple, inds::AbstractVector{<:Bool})
    # get the new values
    new_values = vt.values[inds]
    # process the new ranges and indices
    inv_inds = .!inds
    prev_sum = 0
    delete_inds = Int[]
    new_ranges = copy(vt.ranges)
    new_indices = collect(vt.indices)
    for i in eachindex(new_ranges)
        delete_sum = sum(inv_inds[new_ranges[i]])
        if length(new_ranges[i]) == delete_sum
            push!(delete_inds, i)
        else
            group_inds = new_ranges[i][inv_inds[new_ranges[i]]] .+ (1 - new_ranges[i].start)
            if length(group_inds) != 0
                new_indices[i] = _update_indices(new_indices[i], group_inds)
            end
            new_ranges[i] = UnitRange(new_ranges[i].start - prev_sum,
                                      new_ranges[i].stop - prev_sum - delete_sum)
        end
        prev_sum += delete_sum
    end
    deleteat!(new_ranges, delete_inds)
    deleteat!(new_indices, delete_inds)
    # make the new VectorTuple and return 
    return VectorTuple(new_values, new_ranges, Tuple(i for i in new_indices))
end

## Make the original array for a particular tuple element
# AbstractArray
function _make_array(vals::Vector, inds::ContainerIndices)
    return unvectorize(vals, inds)
end

# Other
function _make_array(vals::Vector, inds::Nothing)
    return length(vals) == 1 ? first(vals) : vals
end

# Define Tuple construction from a VectorTuple
function Base.Tuple(vt::VectorTuple; use_indices = true)
    if use_indices
        return Tuple(_make_array(vt.values[vt.ranges[i]], vt.indices[i])
                     for i in eachindex(vt.ranges))
    else
        return Tuple(length(vt.ranges[i]) == 1 ? vt.values[vt.ranges[i].start] : vt.values[vt.ranges[i]]
                     for i in eachindex(vt.ranges))
    end
end

# Define Tuple using construction from a values vector in combination with a VT
function Base.Tuple(vals::Vector, vt::VectorTuple; use_indices = true)
    @assert length(vals) == length(vt)
    if use_indices
        return Tuple(_make_array(vals[vt.ranges[i]], vt.indices[i])
                     for i in eachindex(vt.ranges))
    else
        return Tuple(length(vt.ranges[i]) == 1 ? vals[vt.ranges[i].start] : vals[vt.ranges[i]]
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
