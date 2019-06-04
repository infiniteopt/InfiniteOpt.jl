const InfiniteAffExpr = JuMP.GenericAffExpr{Float64, Union{InfiniteVariableRef, FiniteVariableRef}}
const FiniteAffExpr = JuMP.GenericAffExpr{Float64, FiniteVariableRef}

# function _expr_type_parser(type::DataType)
#     if isa(type, Union{InfiniteVariableRef, PointVariableRef})
#         return
# end

# function JuMP.GenericAffExpr(constant::V, kv1::Pair{K, V}, kv2::Pair{W, V}) where {K, W, V}
#     return JuMP.GenericAffExpr{V, Union{K, W}}(constant, JuMP._new_ordered_dict(Union{K, W}, V, kv1, kv2))
# end
