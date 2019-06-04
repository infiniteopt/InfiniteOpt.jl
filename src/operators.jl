# AbstractVariableRef--AbstractVariableRef
# Base.:+(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef} = JuMP.GenericAffExpr(0.0, lhs => 1.0, rhs => 1.0)
# Base.:-(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef} = JuMP.GenericAffExpr(0.0, lhs => 1.0, rhs => -1.0)
# function Base.:*(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
#     JuMP.GenericQuadExpr(JuMP.GenericAffExpr{Float64, Union{V, W}}(), JuMP.UnorderedPair{Union{V, W}}(lhs, rhs) => 1.0)
# end

# GenericAffExpr--AbstractVariableRef
# function Base.:+(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
#     return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, Union{V, W}}(lhs.constant, lhs.terms), one(C), rhs)
# end
# function Base.:-(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
#     return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, Union{V, W}}(lhs.constant, lhs.terms), -one(C), rhs)
# end
# function Base.:*(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
#     if !iszero(lhs.constant)
#         result = JuMP.GenericQuadExpr{C, Union{V, W}}(lhs.constant * rhs)
#     else
#         result = zero(JuMP.GenericQuadExpr{C, Union{V, W}})
#     end
#     for (coef, var) in JuMP.linear_terms(lhs)
#         JuMP.add_to_expression!(result, coef, var, rhs)
#     end
#     return result
# end
