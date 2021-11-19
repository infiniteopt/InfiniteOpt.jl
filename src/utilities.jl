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
        isequal(value, e1) || return false
    end
    return true
end

# Extend comparison for JuMP.VariableInfo 
function Base.:(==)(info1::JuMP.VariableInfo, info2::JuMP.VariableInfo)::Bool 
    return info1.has_lb == info2.has_lb && 
           (!info1.has_lb || info1.lower_bound == info2.lower_bound) && 
           info1.has_ub == info2.has_ub && 
           (!info1.has_ub || info1.upper_bound == info2.upper_bound) && 
           info1.has_fix == info2.has_fix && 
           (!info1.has_fix || info1.fixed_value == info2.fixed_value) && 
           info1.has_start == info2.has_start && 
           (!info1.has_start || info1.start == info2.start) && 
           info1.binary == info2.binary && info1.integer == info2.integer
end

## Convert JuMP variable info to only use Float64
# Change needed
function _make_float_info(info::JuMP.VariableInfo)
    return JuMP.VariableInfo(info.has_lb, Float64(info.lower_bound), info.has_ub, 
                             Float64(info.upper_bound), info.has_fix, 
                             Float64(info.fixed_value), info.has_start, 
                             Float64(info.start), info.binary, info.integer)
end

# No change needed
function _make_float_info(
    info::JuMP.VariableInfo{Float64, Float64, Float64, Float64}
    )
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
