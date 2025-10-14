"""
    VectorTuple{T}

A collection DataType for storing a `Tuple` of singular elements of type `T`
and/or `Array{T, N}`s in a convenient vector form that utilizes
linear indexing. Here `I` is denotes the type of a `Tuple` that stores the 
indices of each tuple element as given by `indices`. `VectorTuple`s should be 
defined from an original tuple via `VectorTuple(tuple)` or by listing the tuple 
elements `VectorTuple(items...)`. Note this is still an experimental type and is 
primarily intended to store infinite parameter reference tuples and point 
variable support value tuples. Some of the notable capabilities are exemplified 
below.

**Example**
```julia-repl
julia> julia> tuple = (3., [-2., 4.], ones(2, 2))
(3.0, [-2.0, 4.0], [1.0 1.0; 1.0 1.0])

julia> vt = VectorTuple(tuple) # make by listing items (notice everything is a vector)
(3.0, [-2.0, 4.0], [1.0 1.0; 1.0 1.0])

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
(3.0, [-2.0, 4.0], [1.0 1.0; 1.0 1.0])

julia> restricted_copy(vt, [true, false, true]) # make a copy with deleted elements
(3.0, [1.0 1.0; 1.0 1.0])
```
"""
struct VectorTuple{T}
    values::Vector{T}
    ranges::Vector{UnitRange{Int}}
    dimensions::Vector{Int}
    num_columns::Vector{Int}
end

# Extract the dimension of a tuple element
_get_dimension(x) = 0
_get_dimension(::Vector) = 1
_get_dimension(::Matrix) = 2
function _get_dimension(::AbstractArray)
    error("VectorTuple only supports tuples with scalars, ",
          "vectors, and/matrices.")
end

# Extract number of columns (defaults to 1 for anything not a matrix)
_get_numcol(m::Matrix) = size(m)[2]
_get_numcol(x) = 1

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
    # get the dimensions
    dims = Int[_get_dimension(t) for t in tuple]
    numcols = Int[_get_numcol(t) for t in tuple]
    # prepare the values vector
    vals = [j for i in eachindex(tuple) for j in tuple[i]]
    # return everything
    return vals, ranges, dims, numcols
end

# Constructor from a Tuple without type
function VectorTuple(tuple::Tuple)
    vals, ranges, dims, numcols = _make_vt_components(tuple)
    return VectorTuple(vals, ranges, dims, numcols)
end

# Constructor from various arguments (like splatting the tuple) without type
function VectorTuple(items...)
    return VectorTuple(items)
end

# Extend Base.isequal
function Base.isequal(vt1::VectorTuple, vt2::VectorTuple)
    return all(isequal(getproperty(vt1, f), getproperty(vt2, f)) for f in fieldnames(VectorTuple))
end

# Extend Base.:(==)
function Base.:(==)(vt1::VectorTuple, vt2::VectorTuple)
    return isequal(vt1, vt2)
end

# Extend Base.size
function Base.size(vt::VectorTuple, d::Int)
    d > 1 && throw(ArgumentError("invalid VectorTuple dimension $d"))
    return length(vt.ranges)
end

# Determine if 2 VectorTuples have the same structure
function same_structure(vt1::VectorTuple, vt2::VectorTuple)
    same_size = vt1.ranges == vt2.ranges
    same_dims = vt1.dimensions == vt2.dimensions
    same_numcols = vt1.num_columns == vt2.num_columns
    return same_size && same_dims && same_numcols
end

# Extend Base.copy
function Base.copy(vt::VectorTuple)
    return VectorTuple(
        copy(vt.values),
        copy(vt.ranges),
        copy(vt.dimensions), 
        copy(vt.num_columns)
    )
end

## Extend first and last indices
# Base.firstindex (linear)
Base.firstindex(vt::VectorTuple) = firstindex(vt.values)

# Base.firstindex (tuple)
Base.firstindex(vt::VectorTuple, d::Int) = 1

# Base.lastindex (linear)
Base.lastindex(vt::VectorTuple) = lastindex(vt.values)

# Base.lastindex (tuple)
function Base.lastindex(vt::VectorTuple, d::Int)
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
function Base.getindex(vt::VectorTuple, ::Colon, j)
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
function Base.setindex!(vt::VectorTuple, v, i::Int, ::Colon)
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
function Base.setindex!(vt::VectorTuple, v, ::Colon, j::Int)
    linear_indexes = [vt.ranges[k].start - 1 + j for k in eachindex(vt.ranges)]
    return setindex!(vt.values, v, linear_indexes)
end

# Tuple indexing (with colon i and noninteger j)
function Base.setindex!(vt::VectorTuple, v, ::Colon, j)
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

# Make a restricted copy (a copy with deleted tuple elements) # TODO fix this
function restricted_copy(vt::VectorTuple{T}, inds::AbstractVector{<:Bool}) where {T}
    new_dims = vt.dimensions[inds]
    new_numcols = vt.num_columns[inds]
    sliced_ranges = vt.ranges[inds]
    new_values = T[vt.values[i] for is in sliced_ranges for i in is]
    new_ranges = Vector{UnitRange{Int}}(undef, length(new_dims))
    last_idx = 0
    for i in eachindex(new_ranges)
        @inbounds new_ranges[i] = UnitRange(last_idx + 1, length(sliced_ranges[i]) + last_idx)
        last_idx = new_ranges[i].stop
    end
    return VectorTuple{T}(new_values, new_ranges, new_dims, new_numcols)
end

# Define Tuple construction from a VectorTuple
function Base.Tuple(vals::Vector, vt::VectorTuple)
    @assert length(vals) == length(vt)
    return Tuple(begin
        dim = vt.dimensions[i]
        if iszero(dim)
            only(vals[vt.ranges[i]])
        elseif isone(dim)
            vals[vt.ranges[i]]
        else
            numcols = vt.num_columns
            reshape(vals[vt.ranges[i]], :, vt.num_columns[i])
        end 
    end
    for i in eachindex(vt.ranges)
    )
end

# Define Tuple using construction from a values vector in combination with a VT
function Base.Tuple(vt::VectorTuple)
    return Base.Tuple(vt.values, vt)
end

## Enable pretty printing
# Extend Base.string
function Base.string(vt::VectorTuple)
    return string(Tuple(vt))
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
