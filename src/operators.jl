# AbstractVariableRef--AbstractVariableRef
function Base.:+(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericAffExpr(0.0, Pair{type, Float64}(lhs, 1.0), Pair{type, Float64}(rhs, 1.0))
end

function Base.:-(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericAffExpr(0.0, Pair{type, Float64}(lhs, 1.0), Pair{type, Float64}(rhs, -1.0))
end

function Base.:*(lhs::V, rhs::W) where {V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(JuMP.GenericAffExpr{Float64, type}(), JuMP.UnorderedPair{type}(lhs, rhs) => 1.0)
end

# AbstractVariableRef--GenericAffExpr
function Base.:+(lhs::V, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    # For the variables to have the proper order in the result, we need to add the lhs first.
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

function Base.:-(lhs::V, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    # For the variables to have the proper order in the result, we need to add the lhs first.
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

function Base.:*(lhs::V, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    if !iszero(rhs.constant)
        result = GenericQuadExpr{C, type}(lhs*rhs.constant)
    else
        result = zero(JuMP.GenericQuadExpr{C, type})
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, coef, lhs, var)
    end
    return result
end

# GenericAffExpr--AbstractVariableRef
function Base.:+(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, type}(lhs.constant, lhs.terms), one(C), rhs)
end

function Base.:-(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(JuMP.GenericAffExpr{C, type}(lhs.constant, lhs.terms), -one(C), rhs)
end

function Base.:*(lhs::JuMP.GenericAffExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    if !iszero(lhs.constant)
        result = JuMP.GenericQuadExpr{C, type}(lhs.constant * rhs)
    else
        result = zero(JuMP.GenericQuadExpr{C, type})
    end
    for (coef, var) in JuMP.linear_terms(lhs)
        JuMP.add_to_expression!(result, coef, var, rhs)
    end
    return result
end

# AffExpr--AffExpr
function Base.:+(lhs::JuMP.GenericAffExpr{C, V}, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP._JuMPTypes, W <: JuMP._JuMPTypes}
    if length(JuMP.linear_terms(lhs)) > 50 || length(JuMP.linear_terms(rhs)) > 50
        if length(linear_terms(lhs)) > 1
            operator_warn(JuMP.owner_model(first(JuMP.linear_terms(lhs))[2]))
        end
    end
    type = _var_type_parser(V, W)
    return JuMP.add_to_expression!(convert(JuMP.GenericAffExpr{C, type}, lhs), convert(JuMP.GenericAffExpr{C, type}, rhs))
end

function Base.:-(lhs::JuMP.GenericAffExpr{C, V}, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    result = convert(JuMP.GenericAffExpr{C, type}, lhs)
    result.constant -= rhs.constant
    sizehint!(result, length(JuMP.linear_terms(lhs)) + length(JuMP.linear_terms(rhs)))
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var)
    end
    return result
end

function Base.:*(lhs::JuMP.GenericAffExpr{C, V}, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    result = zero(JuMP.GenericQuadExpr{C, type})
    JuMP.add_to_expression!(result, convert(JuMP.GenericAffExpr{C, type}, lhs), convert(JuMP.GenericAffExpr{C, type}, rhs))
    return result
end

# GenericAffExpr--GenericQuadExpr
function Base.:+(lhs::JuMP.GenericAffExpr{C, V}, rhs::JuMP.GenericQuadExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs + rhs.aff, convert(JuMP.GenericQuadExpr{C, type}, rhs).terms)
end

function Base.:-(lhs::JuMP.GenericAffExpr{C, V}, rhs::JuMP.GenericQuadExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    result = -convert(JuMP.GenericQuadExpr{C, type}, rhs)
    # This makes an unnecessary copy of aff, but it's important for a to appear
    # first.
    result.aff = lhs + result.aff
    return result
end

# GenericQuadExpr--AbstractVariableRef
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff + rhs, convert(JuMP.GenericQuadExpr{C, type}, lhs).terms)
end

function Base.:-(lhs::JuMP.GenericQuadExpr{C, V}, rhs::W) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff - rhs, convert(JuMP.GenericQuadExpr{C, type}, lhs).terms)
end

# GenericQuadExpr--GenericAffExpr
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V}, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff + rhs, convert(JuMP.GenericQuadExpr{C, type}, lhs).terms)
end

function Base.:-(lhs::JuMP.GenericQuadExpr{C, V}, rhs::JuMP.GenericAffExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    return JuMP.GenericQuadExpr(lhs.aff - rhs, convert(JuMP.GenericQuadExpr{C, type}, lhs).terms)
end

# GenericQuadExpr--GenericQuadExpr
function Base.:+(lhs::JuMP.GenericQuadExpr{C, V}, rhs::JuMP.GenericQuadExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    result = convert(JuMP.GenericQuadExpr{C, type}, lhs)
    for (coef, var1, var2) in JuMP.quad_terms(rhs)
        JuMP.add_to_expression!(result, coef, var1, var2)
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, coef, var)
    end
    result.aff.constant += rhs.aff.constant
    return result
end

function Base.:-(lhs::JuMP.GenericQuadExpr{C, V}, rhs::JuMP.GenericQuadExpr{C, W}) where {C, V <: JuMP.AbstractVariableRef, W <: JuMP.AbstractVariableRef}
    type = _var_type_parser(V, W)
    result = convert(JuMP.GenericQuadExpr{C, type}, lhs)
    for (coef, var1, var2) in JuMP.quad_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var1, var2)
    end
    for (coef, var) in JuMP.linear_terms(rhs)
        JuMP.add_to_expression!(result, -coef, var)
    end
    result.aff.constant -= rhs.aff.constant
    return result
end
