"""
    ContainerIndices{N, Ax}

A `DataType` for storing the indexing information of arbitrary Julia arrays 
and JuMP containers. This is to be used to enable 
[`vectorize`](@ref InfiniteOpt.Collections.vectorize) and 
[`unvectorize`](@ref InfiniteOpt.Collections.unvectorize). `ContainerIndices` 
should be constructed via [`indices`](@ref InfiniteOpt.Collections.indices).
"""
struct ContainerIndices{N, Ax}
    indices::CartesianIndices{N, NTuple{N, Base.OneTo{Int}}}
    axes::Ax
end

## Make constructors for different array types
# Array
function ContainerIndices(
    arr::Array{T, N}
    )::ContainerIndices{N, Nothing} where {T, N}
    return ContainerIndices(CartesianIndices(arr), nothing)
end

# DenseAxisArray
function ContainerIndices(
    arr::JuMPC.DenseAxisArray{T, N, Ax}
    )::ContainerIndices{N, Ax} where {T, N, Ax}
    return ContainerIndices(CartesianIndices(arr), axes(arr))
end

# SparseAxisArray
function ContainerIndices(
    arr::JuMPC.SparseAxisArray{T, N, K}
    )::ContainerIndices{1, Vector{K}} where {T, N, K}
    axs = sort!(collect(keys(arr.data)))
    return ContainerIndices(CartesianIndices(axs), axs)
end

# Extend Base.:(==)
function Base.:(==)(i1::ContainerIndices, i2::ContainerIndices)::Bool
    return i1.indices == i2.indices && i1.axes == i2.axes
end

# Extend Base.size
function Base.size(inds::ContainerIndices)
    return size(inds.indices) # This should be avoided for SparseAxisArrays
end

# Extend Base.length
function Base.length(inds::ContainerIndices)::Int
    return length(inds.indices)
end

"""
    indices(arr)::Union{ContainerIndices, Nothing}

Return the indices of `arr` to aid with calls made to 
[`vectorize`](@ref InfiniteOpt.Collections.vectorize) and 
[`unvectorize`](@ref InfiniteOpt.Collections.unvectorize).
"""
function indices(arr)::Nothing
    return
end

# Array
function indices(arr::Array{T, N})::ContainerIndices{N, Nothing} where {T, N}
    return ContainerIndices(arr)
end

# DenseAxisArray
function indices(
    arr::JuMPC.DenseAxisArray{T, N, Ax}
    )::ContainerIndices{N, Ax} where {T, N, Ax}
    return ContainerIndices(arr)
end

# SparseAxisArray
function indices(
    arr::JuMPC.SparseAxisArray{T, N, K}
    )::ContainerIndices{1, Vector{K}} where {T, N, K}
    return ContainerIndices(arr)
end

"""
    vectorize(arr, [inds::Union{Nothing, ContainerIndices}])::Vector 

Return the vectorized version of `arr`. If `arr` is not an `AbstractArray` then 
it is returned unchanged. Optionally, the indices `inds` can be given to enhance 
the efficiency if they were already made extracted via the 
[`indices`](@ref InfiniteOpt.Collections.indices) method.
"""
function vectorize(arr::T)::T where {T}
    return arr
end

# Single argument Vector
function vectorize(arr::Vector{T})::Vector{T} where {T}
    return arr
end

# Single argument Array/DenseAxisArray
function vectorize(arr::AbstractArray{T})::Vector{T} where {T}
    if isempty(arr)
        return T[]
    else
        return reduce(vcat, arr)
    end
end

# Single argument SparseAxisArray
function vectorize(arr::JuMPC.SparseAxisArray{T})::Vector{T} where {T}
    ks = sort!(collect(keys(arr.data)))
    return [arr[i] for i in ks]
end

# Double argument singleton
function vectorize(arr::T, inds::Nothing)::T where {T}
    return vectorize(arr)
end

# Double argument Vector/Array/DenseAxisArray
function vectorize(
    arr::AbstractArray{T}, 
    inds::ContainerIndices
    )::Vector{T} where {T}
    return vectorize(arr)
end

# Double argument SparseAxisArray
function vectorize(
    arr::JuMPC.SparseAxisArray{T}, 
    inds::ContainerIndices
    )::Vector{T} where {T}
    return [arr[i] for i in inds.axes]
end

"""
    unvectorize(vect, inds::Union{Nothing, ContainerIndices})

Return the unvectorized version of `vect` as dictated by `inds`. Here `vect` 
comes from a call of [`vectorize`](@ref InfiniteOpt.Collections.vectorize) and 
`inds` comes from a call of [`indices`](@ref InfiniteOpt.Collections.indices).
"""
function unvectorize(vect::T, inds::Nothing)::T where {T}
    return vect
end

# Vector
function unvectorize(
    vect::Vector{T}, 
    inds::ContainerIndices{1, Nothing}
    )::Vector{T} where {T}
    return vect
end

# Array
function unvectorize(
    vect::Vector{T}, 
    inds::ContainerIndices{N, Nothing}
    )::Array{T, N} where {T, N}
    return [vect[i] for i in LinearIndices(inds.indices)]
end

# DenseAxisArray (vector)
function unvectorize(vect::Vector, inds::ContainerIndices{1})
    return JuMPC.DenseAxisArray(vect, inds.axes...)
end

# DenseAxisArray (general)
function unvectorize(vect::Vector, inds::ContainerIndices)
    arr = [vect[i] for i in LinearIndices(inds.indices)]
    return JuMPC.DenseAxisArray(arr, inds.axes...)
end

# SparseAxisArray
function unvectorize(
    vect::Vector{T}, 
    inds::ContainerIndices{1, Vector{K}}
    )::JuMPC.SparseAxisArray{T, N, K} where {T, N, K <: NTuple{N}}
    dict = Dict(inds.axes[i] => vect[i] for i in eachindex(inds.axes))
    return JuMPC.SparseAxisArray(dict)
end
