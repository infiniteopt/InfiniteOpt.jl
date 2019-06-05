# Define infinite expressions
const InfiniteAffExpr = JuMP.GenericAffExpr{Float64, Union{InfiniteVariableRef, GeneralVariableRef}}
const InfiniteQuadExpr = JuMP.GenericQuadExpr{Float64, Union{InfiniteVariableRef, GeneralVariableRef}}

# Define type hierchical parser
function _var_type_parser(V::Type{<: JuMP.AbstractVariableRef}, W::Type{<: JuMP.AbstractVariableRef})
    if V == W
        return V
    elseif V isa Type{<: FiniteVariableRef} && W isa Type{<: FiniteVariableRef}
        return FiniteVariableRef
    elseif V isa Type{<: GeneralVariableRef} && W isa Type{<: GeneralVariableRef}
        return GeneralVariableRef
    else
        return JuMP.AbstractVariableRef
    end
end

# Extend handle mixed variable input
function JuMP.add_to_expression!(quad::JuMP.GenericQuadExpr{C,Z}, new_coef::C,
                            new_var1::V, new_var2::W) where {C,Z<:JuMP.AbstractVariableRef,V,W}
    type = _var_type_parser(Z, W)
    key = JuMP.UnorderedPair{type}(new_var1, new_var2)
    JuMP._add_or_set!(quad.terms, key, new_coef)
    return quad
end

# Extend convert to handle JuMP expression types
Base.convert(::Type{JuMP.GenericAffExpr{C,V}}, x::JuMP.GenericAffExpr) where {C,V} = JuMP.GenericAffExpr{C,V}(x.constant, x.terms)
Base.convert(::Type{JuMP.GenericQuadExpr{C,V}}, x::JuMP.GenericQuadExpr) where {C,V} = JuMP.GenericQuadExpr{C,V}(x.aff, x.terms)
Base.convert(::Type{JuMP.UnorderedPair{T}}, x::JuMP.UnorderedPair) where {T} = JuMP.UnorderedPair{T}(x.a, x.b)
