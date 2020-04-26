## Extend for better comparisons than default
# GenericAffExpr
function Base.:(==)(aff1::JuMP.GenericAffExpr{C, V},
                    aff2::JuMP.GenericAffExpr{C, V}) where {C, V <: GeneralVariableRef}
    return aff1.constant == aff2.constant && collect(pairs(aff1.terms)) == collect(pairs(aff2.terms))
end

# GenericQuadExpr
function Base.:(==)(quad1::JuMP.GenericQuadExpr{C, V},
                    quad2::JuMP.GenericQuadExpr{C, V}) where {C, V <: GeneralVariableRef}
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

## Determine which variables are present in a function
# GeneralVariableRef
function _all_function_variables(f::GeneralVariableRef)::Vector{GeneralVariableRef}
    return [f]
end

# GenericAffExpr
function _all_function_variables(f::JuMP.GenericAffExpr)::Vector{GeneralVariableRef}
    return collect(keys(f.terms))
end

# GenericQuadExpr
function _all_function_variables(f::JuMP.GenericQuadExpr)::Vector{GeneralVariableRef}
    vref_set = Set(keys(f.aff.terms))
    for pair in keys(f.terms)
        union!(vref_set, pair.a)
        union!(vref_set, pair.b)
    end
    return collect(vref_set)
end

# Fallback
function _all_function_variables(f)
    error("`_all_function_variables` not defined for expression of type $(typeof(f)).")
end

## Return the unique set of object numbers in an expression
# Dispatch fallback (--> should be defined for each non-empty variable type)
_object_numbers(expr::DispatchVariableRef)::Vector{Int} = Int[]

# GeneralVariableRef
function _object_numbers(expr::GeneralVariableRef)::Vector{Int}
    return _object_numbers(dispatch_variable_ref(expr))
end

# GenericAffExpr
function _object_numbers(expr::JuMP.GenericAffExpr)::Vector{Int}
    obj_nums = Set{Int}()
    for vref in keys(expr.terms)
        union!(obj_nums, _object_numbers(vref))
    end
    return collect(obj_nums)
end

# GenericQuadExpr
function _object_numbers(expr::JuMP.GenericQuadExpr)::Vector{Int}
    obj_nums = Set(_object_numbers(expr.aff))
    for pair in keys(expr.terms)
        union!(obj_nums, _object_numbers(pair.a))
        union!(obj_nums, _object_numbers(pair.b))
    end
    return collect(obj_nums)
end

# Variable list
function _object_numbers(vrefs::Vector{GeneralVariableRef})::Vector{Int}
    obj_nums = Set{Int}()
    for vref in vrefs
        union!(obj_nums, _object_numbers(vref))
    end
    return collect(obj_nums)
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
