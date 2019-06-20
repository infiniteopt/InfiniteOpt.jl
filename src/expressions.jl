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
_all_function_variables(f::JuMP.GenericAffExpr) = [vref for vref in keys(f.terms)]
function _all_function_variables(f::JuMP.GenericQuadExpr)
    aff_vrefs = _all_function_variables(f.aff)
    vref_pairs = [k for k in keys(f.terms)]
    a_vrefs = [pair.a for pair in vref_pairs]
    b_vrefs = [pair.b for pair in vref_pairs]
    return unique([aff_vrefs; a_vrefs; b_vrefs])
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
    vref_pairs = [k for k in keys(f.terms)]
    for i = 1:length(vref_pairs)
        if vref_pairs[i].a == vref
            delete!(f.terms, vref_pairs[i])
        elseif vref_pairs[i].b == vref
            delete!(f.terms, vref_pairs[i])
        end
    end
    return
end

# Check expression for a particular variable type
# function _has_variable(f::JuMP.AbstractJuMPScalar, vref::GeneralVariableRef)
#     vrefs = _all_function_variables(f)
#     if vref in vrefs
#         return true
#     end
#     # search recursively along nested measures
#     filter!(x -> isa(x, MeasureRef), vrefs)
#     if length(vrefs) != 0
#         for mref in vrefs
#             print(mref)
#             return _has_variable(measure_function(mref), vref)
#         end
#         return false
#     else
#         return false
#     end
# end

# TODO make recursion work
# Check expression for a particular variable type
function _has_variable(f::JuMP.AbstractJuMPScalar, vref::GeneralVariableRef)
    vrefs = _all_function_variables(f)
    if ref == vref
        return true
    end
    if isa(ref, MeasureRef)
        print(ref)
        return _has_variable(measure_function(ref), vref)
    end
end
