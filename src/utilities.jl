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

# Convert numbers to jump objects 
function Base.convert(::Type{JuMP.AbstractJuMPScalar}, c::Number) 
    return zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) + c
end

# Make workaround for keys of containers
_keys(a::JuMPC.SparseAxisArray) = keys(a.data)
_keys(a::AbstractArray) = keys(a)

## Define efficient function to check if all elements in array are equal
# method found at https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal/47578613
@inline function _allequal(x::AbstractArray)::Bool
    length(x) < 2 && return true
    e1 = first(x)
    for value in x
        value == e1 || return false
    end
    return true
end

## Convert JuMP variable info to only use Float64
# Change needed
function _make_float_info(info::JuMP.VariableInfo
    )::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    return JuMP.VariableInfo{Float64, Float64, Float64, Float64}(
        info.has_lb, info.lower_bound, info.has_ub, info.upper_bound,
        info.has_fix, info.fixed_value, info.has_start, info.start,
        info.binary, info.integer)
end

# No change needed
function _make_float_info(
    info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    )::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    return info
end

# Make a string of NamedTuple arguments 
function _kwargs2string(nt::NamedTuple)::String
    if length(nt) == 1
        return string(nt)[2:end-2]
    else
        return string(nt)[2:end-1]
    end
end
