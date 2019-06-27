# Extend convert to handle JuMP expression types
function Base.convert(::Type{JuMP.GenericAffExpr{C,V}}, x::JuMP.GenericAffExpr) where {C,V}
    return JuMP.GenericAffExpr{C,V}(x.constant, x.terms)
end
function Base.convert(::Type{JuMP.GenericQuadExpr{C,V}}, x::JuMP.GenericQuadExpr) where {C,V}
    return JuMP.GenericQuadExpr{C,V}(x.aff, x.terms)
end
function Base.convert(::Type{JuMP.UnorderedPair{T}}, x::JuMP.UnorderedPair) where {T}
    return JuMP.UnorderedPair{T}(x.a, x.b)
end

# Extend convert to handle JuMP containers
function Base.convert(::Type{JuMP.Containers.SparseAxisArray}, arr::Array)
    data = Dict(Tuple(k) => arr[k] for k in CartesianIndices(arr))
    return JuMP.Containers.SparseAxisArray(data)
end

function Base.convert(::Type{JuMP.Containers.SparseAxisArray}, arr::JuMP.Containers.DenseAxisArray)
    data = Dict(k.I => arr[k] for k in keys(arr))
    return JuMP.Containers.SparseAxisArray(data)
end

# function Base.convert(::Type{Array}, arr::JuMP.Containers.SparseAxisArray)
#
# end
#
# function Base.convert(::Type{JuMP.Containers.DenseAxisArray}, arr::JuMP.Containers.SparseAxisArray)
#
# end

# Hack to make the keys function work for sparse arrays
Base.keys(d::JuMP.Containers.SparseAxisArray) = keys(d.data)

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
    if sum(valids_a) == len_terms && sum(valids_b) == len_terms && isa(aff, JuMP.GenericAffExpr{C, type})
        return JuMP.GenericQuadExpr{C, type}(aff, quad.terms)
    else
        return quad
    end
end
