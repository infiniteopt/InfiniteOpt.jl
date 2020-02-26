# TODO change operations to avoid nonlinear expressions obtained from operations done on measures

## GeneralVariableRef--GeneralVariableRef extensions for different types
# var1 + var2
function Base.:+(lhs::V, rhs::W)::JuMP.GenericAffExpr where {V <: GeneralVariableRef,
                                                             W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericAffExpr(0.0, Pair{type, Float64}(lhs, 1.0),
                               Pair{type, Float64}(rhs, 1.0))
end

# var1 - var2
function Base.:-(lhs::V, rhs::W)::JuMP.GenericAffExpr where {V <: GeneralVariableRef,
                                                             W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    if lhs == rhs
        return zero(JuMP.GenericAffExpr{Float64, type})
    else
        return JuMP.GenericAffExpr(0.0, Pair{type, Float64}(lhs, 1.0),
                                   Pair{type, Float64}(rhs, -1.0))
    end
end

# var1 * var2
function Base.:*(lhs::V, rhs::W)::JuMP.GenericQuadExpr where {V <: GeneralVariableRef,
                                                              W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(JuMP.GenericAffExpr{Float64, type}(),
                                JuMP.UnorderedPair{type}(lhs, rhs) => 1.0)
end

# GeneralVariableRef--GenericAffExpr
# var + aff
function Base.:+(lhs::V,rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericAffExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    # For the variables to have the proper order in the result,
    # we need to add the lhs first.
    type = _var_type_parser(V, W)
    result = zero(JuMP.GenericAffExpr{C, type})
    result.constant = rhs.constant
    sizehint!(result, length(JuMP.linear_terms(rhs)) + 1)
    JuMP.add_to_expression!(result, one(C), lhs)
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, coef, var)
    end
    return result
end

# var - aff
function Base.:-(lhs::V, rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericAffExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    # For the variables to have the proper order in the result,
    # we need to add the lhs first.
    type = _var_type_parser(V, W)
    result = zero(JuMP.GenericAffExpr{C, type})
    result.constant = -rhs.constant
    sizehint!(result, length(JuMP.linear_terms(rhs)) + 1)
    JuMP.add_to_expression!(result, one(C), lhs)
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var)
    end
    return result
end

# var * aff
function Base.:*(lhs::V, rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    if !iszero(rhs.constant)
        result = JuMP.GenericQuadExpr{C, type}(convert(JuMP.GenericAffExpr{C, type},
                                                       lhs * rhs.constant))
    else
        result = zero(JuMP.GenericQuadExpr{C, type})
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, coef, lhs, var)
    end
    return result
end

# GeneralVariableRef--GenericQuadExpr
# var + quad
function Base.:+(v::V, q::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    new_q = convert(JuMP.GenericQuadExpr{C, type}, q)
    return JuMP.GenericQuadExpr(v + new_q.aff, copy(new_q.terms))
 end

 # var - quad
 function Base.:-(v::V, q::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                  V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = convert(JuMP.GenericQuadExpr{C, type}, -q)
    # This makes an unnecessary copy of aff, but it's important for v to appear
    # first.
    result.aff = v + result.aff
    return result
end

## Extend operators for GenericAffExpr--GeneralVariableRef of different types
# aff + var
function Base.:+(lhs::JuMP.GenericAffExpr{C, V}, rhs::W)::JuMP.GenericAffExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, type}(lhs.constant,
                                   lhs.terms), one(C), rhs)
end

# aff - var
function Base.:-(lhs::JuMP.GenericAffExpr{C, V}, rhs::W)::JuMP.GenericAffExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, type}(lhs.constant,
                                   lhs.terms), -one(C), rhs)
end

# aff * var
function Base.:*(lhs::JuMP.GenericAffExpr{C, V}, rhs::W)::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    if !iszero(lhs.constant)
        result = JuMP.GenericQuadExpr{C, type}(convert(JuMP.GenericAffExpr{C,
                                                     type}, lhs.constant * rhs))
    else
        result = zero(JuMP.GenericQuadExpr{C, type})
    end
    for (coef, var) in JuMP.linear_terms(lhs)
        JuMP.add_to_expression!(result, coef, var, rhs)
    end
    return result
end

## Extend operators for AffExpr--AffExpr for different types
# aff + aff
function Base.:+(lhs::JuMP.GenericAffExpr{C, V},
                 rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericAffExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    if length(JuMP.linear_terms(lhs)) > 50 || length(JuMP.linear_terms(rhs)) > 50
        if length(JuMP.linear_terms(lhs)) > 1
            JuMP.operator_warn(JuMP.owner_model(first(JuMP.linear_terms(lhs))[2]))
        end
    end
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(copy(convert(JuMP.GenericAffExpr{C, type}, lhs)),
                                   convert(JuMP.GenericAffExpr{C, type}, rhs))
end

# aff - aff
function Base.:-(lhs::JuMP.GenericAffExpr{C, V},
                 rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericAffExpr  where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = copy(convert(JuMP.GenericAffExpr{C, type}, lhs))
    result.constant -= rhs.constant
    sizehint!(result, length(JuMP.linear_terms(lhs)) +
                      length(JuMP.linear_terms(rhs)))
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var)
    end
    return result
end

# aff * aff
function Base.:*(lhs::JuMP.GenericAffExpr{C, V},
                 rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = zero(JuMP.GenericQuadExpr{C, type})
    JuMP.add_to_expression!(result, convert(JuMP.GenericAffExpr{C, type}, lhs),
                            convert(JuMP.GenericAffExpr{C, type}, rhs))
    return result
end

## Extend for operators on GenericAffExpr--GenericQuadExpr of different types
# aff + quad
function Base.:+(lhs::JuMP.GenericAffExpr{C, V},
                 rhs::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs + rhs.aff,
                                copy(convert(JuMP.GenericQuadExpr{C, type}, rhs).terms))
end

# aff - quad
function Base.:-(lhs::JuMP.GenericAffExpr{C, V},
                 rhs::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = -convert(JuMP.GenericQuadExpr{C, type}, rhs)
    # This makes an unnecessary copy of aff, but it's important for a to appear
    # first.
    result.aff = lhs + result.aff
    return result
end

## Extend operators for GenericQuadExpr--GeneralVariableRef with different types
# quad + var
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V}, rhs::W)::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff + rhs,
                                copy(convert(JuMP.GenericQuadExpr{C, type}, lhs).terms))
end

# quad - var
function Base.:-(lhs::JuMP.GenericQuadExpr{C, V}, rhs::W)::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff - rhs,
                                copy(convert(JuMP.GenericQuadExpr{C, type}, lhs).terms))
end

# Extend operators for GenericQuadExpr--GenericAffExpr with different types
# quad + aff
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V},
                 rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff + rhs,
                                copy(convert(JuMP.GenericQuadExpr{C, type}, lhs).terms))
end

# quad - aff
function Base.:-(lhs::JuMP.GenericQuadExpr{C, V},
                 rhs::JuMP.GenericAffExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff - rhs,
                                copy(convert(JuMP.GenericQuadExpr{C, type}, lhs).terms))
end

## Extend operators for GenericQuadExpr--GenericQuadExpr with differnt types
# quad + quad
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V},
                 rhs::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = copy(convert(JuMP.GenericQuadExpr{C, type}, lhs))
    for (coef, var1, var2) in JuMP.quad_terms(rhs)
        JuMP.add_to_expression!(result, coef, var1, var2)
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, coef, var)
    end
    result.aff.constant += rhs.aff.constant
    return result
end

# quad - quad
function Base.:-(lhs::JuMP.GenericQuadExpr{C, V},
                 rhs::JuMP.GenericQuadExpr{C, W})::JuMP.GenericQuadExpr where {C,
                 V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    result = copy(convert(JuMP.GenericQuadExpr{C, type}, lhs))
    for (coef, var1, var2) in JuMP.quad_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var1, var2)
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var)
    end
    result.aff.constant -= rhs.aff.constant
    return result
end
