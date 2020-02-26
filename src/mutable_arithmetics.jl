# Abstract expression argument
const _ArgExpr{C, V} = Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}, V}

## Return the promoted datastructure for the indicated variable type
# GenericAffExpr
function _arg_promote(type::Type{<:GeneralVariableRef},
                      expr::JuMP.GenericAffExpr{C, V}) where {C,
                      V <: GeneralVariableRef}
    return convert(JuMP.GenericAffExpr{C, type}, expr)
end

# GenericQuadExpr
function _arg_promote(type::Type{<:GeneralVariableRef},
                      expr::JuMP.GenericQuadExpr{C, V}) where {C,
                      V <: GeneralVariableRef}
    return convert(JuMP.GenericQuadExpr{C, type}, expr)
end

# Fallback (variables, constants, etc.)
function _arg_promote(type::Type{<:GeneralVariableRef}, var::GeneralVariableRef)
    return JuMP.GenericAffExpr(0.0, Pair{type, Float64}(var, 1.0))
end

## Extend MutableArithmetics.mutable_operate to handle mixed variable types
# 2 argument adding
function _MA.mutable_operate!(::typeof(+),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, x)
end

# 2 argument adding/multiplying
function _MA.mutable_operate!(::typeof(_MA.add_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, x)
end

# 2 argument substracting
function _MA.mutable_operate!(::typeof(-),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, -1, x)
end

# 2 argument substracting/multiplying
function _MA.mutable_operate!(::typeof(_MA.sub_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, -1, x)
end

# 3 argument adding/multiplying (no constants)
function _MA.mutable_operate!(::typeof(_MA.add_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W},
                              y::_ArgExpr{C, Q}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef,
                              Q <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    type = _var_type_parser(type, Q)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    y = _arg_promote(type, y)
    return JuMP.add_to_expression!(expr, x, y)
end

# 3 argument adding/multiplying (2nd constant)
function _MA.mutable_operate!(::typeof(_MA.add_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W},
                              y::JuMP._Constant) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, x, y)
end

# 3 argument adding/multiplying (1st constant)
function _MA.mutable_operate!(::typeof(_MA.add_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::JuMP._Constant,
                              y::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    y = _arg_promote(type, y)
    return JuMP.add_to_expression!(expr, x, y)
end

# 3 argument substracting/multiplying (no constants)
function _MA.mutable_operate!(::typeof(_MA.sub_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W},
                              y::_ArgExpr{C, Q}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef,
                              Q <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    type = _var_type_parser(type, Q)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    y = _arg_promote(type, y)
    return JuMP.add_to_expression!(expr, -x, y)
end

# 3 argument substracting/multiplying (1st constant)
function _MA.mutable_operate!(::typeof(_MA.sub_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::JuMP._Constant,
                              y::_ArgExpr{C, W}) where {C,
                              V <: GeneralVariableRef, W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    y = _arg_promote(type, y)
    return JuMP.add_to_expression!(expr, -x, y)
end

# 3 argument substracting/multiplying (terminal constant)
function _MA.mutable_operate!(::typeof(_MA.sub_mul),
                              expr::JuMP._GenericAffOrQuadExpr{C, V},
                              x::_ArgExpr{C, W},
                              y::JuMP._Constant) where {C, V <: GeneralVariableRef,
                                                        W <: GeneralVariableRef}
    type = _var_type_parser(V, W)
    expr = _arg_promote(type, expr)
    x = _arg_promote(type, x)
    return JuMP.add_to_expression!(expr, x, -y)
end
