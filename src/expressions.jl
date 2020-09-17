################################################################################
#                              COMPARISON METHODS
################################################################################
## Extend for better comparisons than default
# GenericAffExpr
function Base.:(==)(aff1::JuMP.GenericAffExpr{C, V},
    aff2::JuMP.GenericAffExpr{C, V}
    )::Bool where {C, V <: GeneralVariableRef}
    return aff1.constant == aff2.constant && aff1.terms == aff2.terms
end

# GenericQuadExpr
function Base.:(==)(quad1::JuMP.GenericQuadExpr{C, V},
    quad2::JuMP.GenericQuadExpr{C, V}
    )::Bool where {C, V <: GeneralVariableRef}
    pairs1 = collect(pairs(quad1.terms))
    pairs2 = collect(pairs(quad2.terms))
    if length(pairs1) != length(pairs2)
        return false
    end
    for i in eachindex(pairs1)
        if pairs1[i][1].a != pairs2[i][1].a || pairs1[i][1].b != pairs2[i][1].b ||
           pairs1[i][2] != pairs2[i][2]
            return false
        end
    end
    return quad1.aff == quad2.aff
end

################################################################################
#                            VARIABLE LIST MAKING
################################################################################
## Determine which variables are present in a function
# GeneralVariableRef
function _all_function_variables(
    f::GeneralVariableRef
    )::Vector{GeneralVariableRef}
    return [f]
end

# GenericAffExpr
function _all_function_variables(
    f::JuMP.GenericAffExpr{C, V}
    )::Vector{V} where {C, V <: GeneralVariableRef}
    return collect(keys(f.terms))
end

# GenericQuadExpr
function _all_function_variables(
    f::JuMP.GenericQuadExpr{C, V}
    )::Vector{V} where {C, V <: GeneralVariableRef}
    vref_set = Set(keys(f.aff.terms))
    for pair in keys(f.terms)
        push!(vref_set, pair.a)
        push!(vref_set, pair.b)
    end
    return collect(vref_set)
end

# Fallback
function _all_function_variables(f)
    error("`_all_function_variables` not defined for expression of type $(typeof(f)).")
end

################################################################################
#                            OBJECT NUMBER METHODS
################################################################################
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
    obj_nums = Set{Int}()
    for vref in keys(expr.aff.terms)
        union!(obj_nums, _object_numbers(vref))
    end
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

################################################################################
#                             PARAMETER NUMBER METHODS
################################################################################
## Return the unique set of parameter numbers in an expression
# Dispatch fallback (--> should be defined for each non-empty variable type)
_parameter_numbers(expr::DispatchVariableRef)::Vector{Int} = Int[]

# GeneralVariableRef
function _parameter_numbers(expr::GeneralVariableRef)::Vector{Int}
    return _parameter_numbers(dispatch_variable_ref(expr))
end

# GenericAffExpr
function _parameter_numbers(expr::JuMP.GenericAffExpr)::Vector{Int}
    par_nums = Set{Int}()
    for vref in keys(expr.terms)
        union!(par_nums, _parameter_numbers(vref))
    end
    return collect(par_nums)
end

# GenericQuadExpr
function _parameter_numbers(expr::JuMP.GenericQuadExpr)::Vector{Int}
    par_nums = Set{Int}()
    for vref in keys(expr.aff.terms)
        union!(par_nums, _parameter_numbers(vref))
    end
    for pair in keys(expr.terms)
        union!(par_nums, _parameter_numbers(pair.a))
        union!(par_nums, _parameter_numbers(pair.b))
    end
    return collect(par_nums)
end

################################################################################
#                             MODEL EXTRACTION METHODS
################################################################################
## Get the model from an expression
# GeneralVariableRef
function _model_from_expr(expr::GeneralVariableRef)::InfiniteModel
    return JuMP.owner_model(expr)
end

# AffExpr
function _model_from_expr(expr::JuMP.GenericAffExpr)::Union{InfiniteModel, Nothing}
    if isempty(expr.terms)
        return
    else
        return JuMP.owner_model(first(keys(expr.terms)))
    end
end

# QuadExpr
function _model_from_expr(expr::JuMP.GenericQuadExpr)::Union{InfiniteModel, Nothing}
    result = _model_from_expr(expr.aff)
    if result !== nothing
        return result
    elseif isempty(expr.terms)
        return
    else
        return JuMP.owner_model(first(keys(expr.terms)).a)
    end
end

# Vector{GeneralVariableRef}
function _model_from_expr(vrefs::Vector{GeneralVariableRef})::Union{InfiniteModel, Nothing}
    if isempty(vrefs)
        return
    else
        return JuMP.owner_model(first(vrefs))
    end
end

# Fallback
function _model_from_expr(expr)
    error("`_model_from_expr` not defined for expr of type $(typeof(expr)).")
end

################################################################################
#                            VARIABLE REMOVAL BOUNDS
################################################################################
## Delete variables from an expression
# GenericAffExpr
function _remove_variable(f::JuMP.GenericAffExpr,
                          vref::GeneralVariableRef)::Nothing
    if haskey(f.terms, vref)
        delete!(f.terms, vref)
    end
    return
end

# GenericQuadExpr
function _remove_variable(f::JuMP.GenericQuadExpr,
                          vref::GeneralVariableRef)::Nothing
    _remove_variable(f.aff, vref)
    for (pair, coef) in f.terms # TODO should we preserve non-vref in pair? --> this would differ from JuMP
        if pair.a == vref || pair.b == vref
            delete!(f.terms, pair)
        end
    end
    return
end

################################################################################
#                         COEFFICIENT MODIFICATION METHODS
################################################################################
## Modify linear coefficient of variable in expression
# GeneralVariableRef
function _set_variable_coefficient!(expr::GeneralVariableRef,
    var::GeneralVariableRef,
    coeff::Real
    )::JuMP.GenericAffExpr{Float64, GeneralVariableRef}
    # Determine if variable is that of the expression and change accordingly
    if expr == var
        return Float64(coeff) * var
    else
        return expr + Float64(coeff) * var
    end
end

# GenericAffExpr
function _set_variable_coefficient!(expr::JuMP.GenericAffExpr{C, V},
    var::V,
    coeff::Real
    )::JuMP.GenericAffExpr{C, V} where {C, V <: GeneralVariableRef}
    # Determine if variable is in the expression and change accordingly
    if haskey(expr.terms, var)
        expr.terms[var] = coeff
        return expr
    else
        return expr + coeff * var
    end
end

# GenericQuadExpr
function _set_variable_coefficient!(expr::JuMP.GenericQuadExpr{C, V},
    var::V,
    coeff::Real
    )::JuMP.GenericQuadExpr{C, V} where {C, V <: GeneralVariableRef}
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

################################################################################
#                         PARAMETER REFERNCE METHODS
################################################################################
## Return an element of a parameter reference tuple given the model and index
# IndependentParameterIndex
function _make_param_tuple_element(model::InfiniteModel,
    idx::IndependentParameterIndex,
    )::GeneralVariableRef
    return _make_parameter_ref(model, idx)
end

# DependentParametersIndex
function _make_param_tuple_element(model::InfiniteModel,
    idx::DependentParametersIndex,
    )::Vector{GeneralVariableRef}
    dpref = DependentParameterRef(model, DependentParameterIndex(idx, 1))
    num_params = _num_parameters(dpref)
    return [GeneralVariableRef(model, idx.value, DependentParameterIndex, i)
            for i in 1:num_params]
end

"""
    parameter_refs(expr)::Tuple

Return the tuple of infinite parameter references that determine the infinite
dependencies of `expr`.

**Example**
```julia-repl
julia> parameter_refs(my_expr)
(t,)
```
"""
function parameter_refs(expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr})::Tuple
    model = _model_from_expr(expr)
    if model === nothing
        return ()
    else
        obj_nums = sort!(_object_numbers(expr))
        obj_indices = _param_object_indices(model)[obj_nums]
        return Tuple(_make_param_tuple_element(model, idx) for idx in obj_indices)
    end
end

################################################################################
#                             EXPRESSION SEARCHING
################################################################################
# Check expression for a particular variable type via a recursive search
#=
function _has_variable(vrefs::Vector{GeneralVariableRef},
                       vref::GeneralVariableRef; prior=GeneralVariableRef[]
                       )::Bool
    if vrefs[1] == vref
        return true
    elseif _index_type(vrefs[1]) == MeasureIndex
        dvref = dispatch_variable_ref(vrefs[1])
        if length(vrefs) > 1
            return _has_variable(_all_function_variables(measure_function(dvref)),
                                 vref, prior = append!(prior, vrefs[2:end]))
        else
            return _has_variable(_all_function_variables(measure_function(dvref)),
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
=#
