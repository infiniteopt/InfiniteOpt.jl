## Extend for better comparisons than default
# GenericAffExpr
function Base.:(==)(aff1::JuMP.GenericAffExpr{C, V},
                    aff2::JuMP.GenericAffExpr{C, W}) where {C, V <: GeneralVariableRef,
                                                            W <: GeneralVariableRef}
    return aff1.constant == aff2.constant && collect(pairs(aff1.terms)) == collect(pairs(aff2.terms))
end

# GenericQuadExpr
function Base.:(==)(quad1::JuMP.GenericQuadExpr{C, V},
                    quad2::JuMP.GenericQuadExpr{C, W}) where {C, V <: GeneralVariableRef,
                                                              W <: GeneralVariableRef}
    pairs1 = collect(pairs(quad1.terms))
    pairs2 = collect(pairs(quad2.terms))
    if length(pairs1) != length(pairs2)
        return false
    end
    for i in eachindex(pairs1)
        if pairs1[i][1].a != pairs2[i][1].a || pairs1[i][1].b != pairs2[i][1].b || pairs1[i][2] != pairs2[i][2]
            return false
        end
    end
    return quad1.aff == quad2.aff
end

# TODO make variable iterator by extending Base.iterate with a custom type
# TODO see https://docs.julialang.org/en/latest/manual/interfaces/#man-interface-iteration-1
## Determine which variables are present in a function
# GeneralVariableRef
function _all_function_variables(f::GeneralVariableRef)::Vector{<:GeneralVariableRef}
    return [f]
end

# GenericAffExpr
function _all_function_variables(f::JuMP.GenericAffExpr)::Vector{<:GeneralVariableRef}
    return GeneralVariableRef[vref for vref in keys(f.terms)]
end

# GenericQuadExpr
function _all_function_variables(f::JuMP.GenericQuadExpr)::Vector{<:GeneralVariableRef}
    aff_vrefs = _all_function_variables(f.aff)
    vref_pairs = [k for k in keys(f.terms)]
    a_vrefs = GeneralVariableRef[pair.a for pair in vref_pairs]
    b_vrefs = GeneralVariableRef[pair.b for pair in vref_pairs]
    return unique([aff_vrefs; a_vrefs; b_vrefs])
end

# Fallback
function _all_function_variables(f)
    error("Can only use InfiniteOpt variables and expressions.")
end

## Return a tuple of the parameter references in an expr
# FiniteVariableRef
_all_parameter_refs(expr::FiniteVariableRef)::Tuple = ()

# InfiniteVariableRef
_all_parameter_refs(expr::InfiniteVariableRef)::Tuple = Tuple(raw_parameter_refs(expr), use_indices = false)

# ParameterRef
_all_parameter_refs(expr::ParameterRef)::Tuple = (expr, )

# ReducedInfiniteVariableRef
_all_parameter_refs(expr::ReducedInfiniteVariableRef)::Tuple = Tuple(raw_parameter_refs(expr), use_indices = false)

# GenericAffExpr
function _all_parameter_refs(expr::JuMP.GenericAffExpr{C,
                              <:GeneralVariableRef})::Tuple where {C}
    pref_list = []
    for var in keys(expr.terms)
        push!(pref_list, _all_parameter_refs(var)...)
    end
    groups = _group.(pref_list)
    unique_groups = unique(groups)
    return Tuple(pref_list[findfirst(isequal(unique_groups[i]), groups)]
                 for i in eachindex(unique_groups))
end

# GenericQuadExpr
function _all_parameter_refs(expr::JuMP.GenericQuadExpr{C,
                             <:GeneralVariableRef})::Tuple where {C}
    pref_list = Any[i for i in _all_parameter_refs(expr.aff)]
    for pair in keys(expr.terms)
        push!(pref_list, _all_parameter_refs(pair.a)...)
        push!(pref_list, _all_parameter_refs(pair.b)...)
    end
    groups = _group.(pref_list)
    unique_groups = unique(groups)
    return Tuple(pref_list[findfirst(isequal(unique_groups[i]), groups)]
                 for i in eachindex(unique_groups))
end

## Delete variables from an expression
# GenericAffExpr
function _remove_variable(f::JuMP.GenericAffExpr, vref::GeneralVariableRef)
    if haskey(f.terms, vref)
        delete!(f.terms, vref)
    end
    return
end

# GenericQuadExpr
function _remove_variable(f::JuMP.GenericQuadExpr, vref::GeneralVariableRef)
    _remove_variable(f.aff, vref)
    vref_pairs = [k for k in keys(f.terms)]
    for i in eachindex(vref_pairs) # TODO retain the good variable in pair
        if vref_pairs[i].a == vref
            delete!(f.terms, vref_pairs[i])
        elseif vref_pairs[i].b == vref
            delete!(f.terms, vref_pairs[i])
        end
    end
    return
end

## Modify linear coefficient of variable in expression
# GeneralVariableRef
function _set_variable_coefficient!(expr::GeneralVariableRef,
                                    var::GeneralVariableRef,
                                    coeff::Real)::JuMP.GenericAffExpr
    # Determine if variable is that of the expression and change accordingly
    if expr == var
        return coeff * var
    else
        return expr + coeff * var
    end
end

# GenericAffExpr
function _set_variable_coefficient!(expr::JuMP.GenericAffExpr,
                                    var::GeneralVariableRef,
                                    coeff::Real)::JuMP.GenericAffExpr
    # Determine if variable is in the expression and change accordingly
    if haskey(expr.terms, var)
        expr.terms[var] = coeff
        return expr
    else
        return expr + coeff * var
    end
end

# GenericQuadExpr
function _set_variable_coefficient!(expr::JuMP.GenericQuadExpr,
                                    var::GeneralVariableRef,
                                    coeff::Real)::JuMP.GenericQuadExpr
    # Determine if variable is in the expression and change accordingly
    if haskey(expr.aff.terms, var)
        expr.aff.terms[var] = coeff
        return expr
    else
        return expr + coeff * var
    end
end

# Fallback
function _set_variable_coefficient!(expr, var::GeneralVariableRef, coeff::Real)
    error("Unsupported function type for coefficient modification.")
end

# Check expression for a particular variable type via a recursive search
# This is tested in test/measures.jl
function _has_variable(vrefs::Vector{<:GeneralVariableRef},
                       vref::GeneralVariableRef; prior=[])
    if vrefs[1] == vref
        return true
    elseif isa(vrefs[1], MeasureRef)
        if length(vrefs) > 1
            return _has_variable(_all_function_variables(measure_function(vrefs[1])),
                          vref, prior = GeneralVariableRef[prior; vrefs[2:end]])
        else
            return _has_variable(_all_function_variables(measure_function(vrefs[1])),
                                 vref, prior = prior)
        end
    elseif length(vrefs) > 1
        return _has_variable(vrefs[2:end], vref, prior = prior)
    elseif length(prior) > 0
        return _has_variable(prior, vref)
    else
        return false
    end
end
