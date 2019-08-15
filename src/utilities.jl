## Extend convert to handle JuMP expression types
# GenericAffExpr
function Base.convert(::Type{JuMP.GenericAffExpr{C,V}},
                      x::JuMP.GenericAffExpr) where {C,V}
    return JuMP.GenericAffExpr{C,V}(x.constant, x.terms)
end

# GenericQuadExpr
function Base.convert(::Type{JuMP.GenericQuadExpr{C,V}},
                      x::JuMP.GenericQuadExpr) where {C,V}
    return JuMP.GenericQuadExpr{C,V}(x.aff, x.terms)
end

# UnorderedPair
function Base.convert(::Type{JuMP.UnorderedPair{T}},
                      x::JuMP.UnorderedPair) where {T}
    return JuMP.UnorderedPair{T}(x.a, x.b)
end

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

# function Base.convert(::Type{Array}, arr::JuMPC.SparseAxisArray)
#
# end
#
# function Base.convert(::Type{JuMPC.DenseAxisArray},
#                       arr::JuMPC.SparseAxisArray)
#
# end

# Hack to make the keys function work for sparse arrays
Base.keys(d::JuMPC.SparseAxisArray) = keys(d.data)
# Hacky fix to compare SparseAxisArrays
function Base.isapprox(a::JuMPC.SparseAxisArray, b::JuMPC.SparseAxisArray)::Bool
    return all(isapprox.(a, b))
end

# Attempt to convert variable type of GenericAffExpr if possible
function _possible_convert(type::DataType,
                           aff::JuMP.GenericAffExpr{C, V}) where {C, V}
    valids = [k isa type for k in keys(aff.terms)]
    if sum(valids) == length(keys(aff.terms))
        return JuMP.GenericAffExpr{C, type}(aff.constant, aff.terms)
    else
        return aff
    end
end

# Attempt to convert variable type of GenericQuadExpr if possible
function _possible_convert(type::DataType,
                           quad::JuMP.GenericQuadExpr{C, V}) where {C, V}
    valids_a = [k.a isa type for k in keys(quad.terms)]
    valids_b = [k.b isa type for k in keys(quad.terms)]
    len_terms = length(keys(quad.terms))
    aff = _possible_convert(type, quad.aff)
    check1 = sum(valids_a) == len_terms
    check2 = sum(valids_b) == len_terms
    check3 = isa(aff, JuMP.GenericAffExpr{C, type})
    if check1 && check2 && check3
        return JuMP.GenericQuadExpr{C, type}(aff, quad.terms)
    else
        return quad
    end
end
