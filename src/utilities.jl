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
function _kwargs2string(nt::NamedTuple)
    if length(nt) == 1
        return string(nt)[2:end-2]
    else
        return string(nt)[2:end-1]
    end
end
