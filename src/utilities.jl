## Extend convert to handle JuMP containers
# Array -> SparseAxisArray
function Base.convert(::Type{JuMPC.SparseAxisArray}, arr::Array)
    data = Dict(Tuple(k) => arr[k] for k in CartesianIndices(arr))
    return JuMPC.SparseAxisArray(data)
end

# DenseAxisArray -> SparseAxisArray
function Base.convert(::Type{JuMPC.SparseAxisArray},
                      arr::JuMPC.DenseAxisArray)
    data = Dict(k.I => arr[k] for k in keys(arr))
    return JuMPC.SparseAxisArray(data)
end

# Hack to make the keys function work for sparse arrays
Base.keys(d::JuMPC.SparseAxisArray) = keys(d.data)

# Hacky fix to compare SparseAxisArrays
function Base.isapprox(a::JuMPC.SparseAxisArray, b::JuMPC.SparseAxisArray)::Bool
    return all(isapprox.(a, b))
end

## Define functions to convert a JuMP array into a vector (need for @BDconstraint)
# AbstractArray
function _make_vector(arr::AbstractArray)
    return [arr[i] for i in keys(arr)]
end

# Array (do nothing)
function _make_vector(arr::Array)
    return arr
end

# Something else
function _make_vector(arr)
    return arr
end

## Define efficient function to check if all elements in array are equal
# method found at https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal/47578613
@inline function _allequal(x::AbstractArray)
    length(x) < 2 && return true
    e1 = first(x)
    for value in x
        value == e1 || return false
    end
    return true
end
