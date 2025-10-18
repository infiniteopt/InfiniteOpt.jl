"""
    vectorize(arr)::Tuple{Vector, Tuple{Int, Int}} 

Return the vectorized version of `arr`. If `arr` is not an `Array` then 
it is returned unchanged. The information needed to unvectorize `arr` is
also provided
"""
vectorize(arr) = arr, (0, 1)
vectorize(arr::Vector) = arr, (1, 1)
function vectorize(arr::Matrix{T}) where {T}
    if isempty(arr) 
        return T[], (2, 0)
    else
        reduce(vcat, arr), (2, size(arr)[2])
    end
end
vectorize(::T) where {T <: AbstractArray} = error("Unsupported container type `$(T)`.")

"""
    unvectorize(vect, info::Tuple{Int, Int})

Return the unvectorized version of `vect` as dictated by `info`. Here `vect` and
`info` come from a call of [`vectorize`](@ref InfiniteOpt.Collections.vectorize).
"""
function unvectorize(vect, info::Tuple{Int, Int})
    if first(info) < 2 
        return vect 
    else
        return reshape(vect, : , last(info))
    end
end
