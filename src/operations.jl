function Base.promote_rule(V::Type{<:GeneralVariableRef}, W::Type{<:GeneralVariableRef})
    type = _var_type_parser(V, W)
    return JuMP.GenericAffExpr{Float64, type}
end
