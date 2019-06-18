# Define type hierchical parser
function _var_type_parser(V::Type{<: JuMP.AbstractVariableRef}, W::Type{<: JuMP.AbstractVariableRef})
    if V == W
        return V
    elseif V isa Type{<: FiniteVariableRef} && W isa Type{<: FiniteVariableRef}
        return FiniteVariableRef
    elseif V isa Type{<: MeasureFiniteVariableRef} && W isa Type{<: MeasureFiniteVariableRef}
        return MeasureFiniteVariableRef
    elseif V isa Type{<: GeneralVariableRef} && W isa Type{<: GeneralVariableRef}
        return GeneralVariableRef
    else
        return JuMP.AbstractVariableRef
    end
end

# Extend handle mixed variable input
function JuMP.add_to_expression!(quad::JuMP.GenericQuadExpr{C, Z}, new_coef::C,
                                 new_var1::V, new_var2::W) where {C, Z <: JuMP.AbstractVariableRef, V, W}
    type = _var_type_parser(Z, W)
    key = JuMP.UnorderedPair{type}(new_var1, new_var2)
    JuMP._add_or_set!(quad.terms, key, new_coef)
    return quad
end

# Extend convert to handle JuMP expression types
Base.convert(::Type{JuMP.GenericAffExpr{C,V}}, x::JuMP.GenericAffExpr) where {C,V} = JuMP.GenericAffExpr{C,V}(x.constant, x.terms)
Base.convert(::Type{JuMP.GenericQuadExpr{C,V}}, x::JuMP.GenericQuadExpr) where {C,V} = JuMP.GenericQuadExpr{C,V}(x.aff, x.terms)
Base.convert(::Type{JuMP.UnorderedPair{T}}, x::JuMP.UnorderedPair) where {T} = JuMP.UnorderedPair{T}(x.a, x.b)

# Determine which variables are present in a function
_all_function_variables(f::GeneralVariableRef) = [f]
_all_function_variables(f::JuMP.GenericAffExpr) = [v for v in keys(f.terms)]
function _all_function_variables(f::JuMP.GenericQuadExpr)
    aff_vars = _all_function_variables(f.aff)
    var_pairs = [k for k in keys(f.terms)]
    a_vars = [pair.a for pair in var_pairs]
    b_vars = [pair.b for pair in var_pairs]
    return unique([aff_vars; a_vars; b_vars])
end

# delete variables from an expression
function _remove_variable(f::JuMP.GenericAffExpr, vref::GeneralVariableRef)
    if haskey(f.terms, vref)
        delete!(f.terms, vref)
    end
    return
end
function _remove_variable(f::JuMP.GenericQuadExpr, vref::GeneralVariableRef)
    _remove_variable(f.aff, vref)
    var_pairs = [k for k in keys(f.terms)]
    for i = 1:length(var_pairs)
        if var_pairs[i].a == vref
            delete!(f.terms, var_pairs[i])
        elseif var_pairs[i].b == vref
            delete!(f.terms, var_pairs[i])
        end
    end
    return
end
