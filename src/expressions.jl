################################################################################
#                              VARIABLE ITERATION
################################################################################
## Create helper methods to interrogate the variables of an expr w/ a function
# Constant
function _interrogate_variables(interrogator::Function, c)
    return
end

# GeneralVariableRef
function _interrogate_variables(
    interrogator::Function, 
    v::JuMP.AbstractVariableRef
    )
    interrogator(v)
    return
end

# AffExpr
function _interrogate_variables(
    interrogator::Function, 
    aff::JuMP.GenericAffExpr
    )
    for (v, _) in aff.terms
        interrogator(v)
    end
    return
end

# QuadExpr
function _interrogate_variables(
    interrogator::Function, 
    quad::JuMP.GenericQuadExpr
    )
    for (p, _) in quad.terms
        interrogator(p.a)
        interrogator(p.b)
    end
    _interrogate_variables(interrogator, quad.aff)
    return
end

# NonlinearExpr (avoid recursion to handle deeply nested expressions)
function _interrogate_variables(interrogator::Function, nlp::JuMP.GenericNonlinearExpr)
    stack = Vector{Any}[nlp.args]
    while !isempty(stack)
        args = pop!(stack)
        for arg in args
            if arg isa JuMP.GenericNonlinearExpr
                push!(stack, arg.args)
            else
                _interrogate_variables(interrogator, arg)
            end
        end
    end
    return
end

# AbstractArray
function _interrogate_variables(interrogator::Function, arr::AbstractArray)
    for ex in arr
        _interrogate_variables(interrogator, ex)
    end
    return
end

################################################################################
#                            VARIABLE LIST MAKING
################################################################################
"""
    all_expression_variables(expr::JuMP.AbstractJuMPScalar)::Vector

Returns a vector of all the variable references contained in `expr`.

**Example**
```julia-repl
julia> all_expr_variables(y^2 + z - t)
3-element Array{GeneralVariableRef,1}:
 y(t)
 z
 t
```
"""
function all_expression_variables(f)
    error("`all_expression_variables` not defined for expression of type $(typeof(f)).")
end

# GeneralVariableRef
function all_expression_variables(f::GeneralVariableRef)
    return [f]
end

# GenericAffExpr
function all_expression_variables(f::JuMP.GenericAffExpr)
    return collect(keys(f.terms))
end

# GenericQuadExpr
function all_expression_variables(f::JuMP.GenericQuadExpr)
    vref_set = Set(keys(f.aff.terms))
    for (pair, _) in f.terms
        push!(vref_set, pair.a, pair.b)
    end
    return collect(vref_set)
end

# NonlinearExpr or array of expressions
function all_expression_variables(f::Union{JuMP.GenericNonlinearExpr, AbstractArray})
    vref_set = Set{GeneralVariableRef}()
    _interrogate_variables(v -> push!(vref_set, v), f)
    return collect(vref_set)
end

################################################################################
#                        GROUP/PARAMETER NUMBER METHODS
################################################################################
## Return the unique set of parameter group integer indices in an expression
# Dispatch fallback (--> should be defined for each non-empty variable type)
parameter_group_int_indices(v::DispatchVariableRef) = Int[]

"""
    parameter_group_int_indices(vref::GeneralVariableRef)::Vector{Int}

Return the list of infinite parameter group integer indices used by `vref`.
"""
function parameter_group_int_indices(v::GeneralVariableRef)
    return parameter_group_int_indices(dispatch_variable_ref(v))
end

"""
    parameter_group_int_indices(expr::JuMP.AbstractJuMPScalar)::Vector{Int}

Return the list of infinite parameter group integer indices used by `expr`.
"""
function parameter_group_int_indices(expr)
    group_int_idxs = Set{Int}()
    _interrogate_variables(v -> union!(group_int_idxs, parameter_group_int_indices(v)), expr)
    return collect(group_int_idxs)
end

################################################################################
#                            VARIABLE REMOVAL BOUNDS
################################################################################
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
    for (pair, _) in f.terms
        if isequal(pair.a, vref) || isequal(pair.b, vref)
            delete!(f.terms, pair)
        end
    end
    return
end

# Helper functions for nonlinear variable deletion
function _remove_variable_from_leaf(c, vref)
    return c
end
function _remove_variable_from_leaf(n_vref::GeneralVariableRef, vref)
    return isequal(n_vref, vref) ? 0.0 : n_vref
end
function _remove_variable_from_leaf(
    ex::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}, 
    vref
    )
    _remove_variable(ex, vref)
    return ex
end

# Nonlinear (avoid recursion to handle deeply nested expressions)
function _remove_variable(f::JuMP.GenericNonlinearExpr, vref::GeneralVariableRef)
    stack = Tuple{Vector{Any}, Int}[]
    for i in eachindex(f.args) # should be reverse, but order doesn't matter
        push!(stack, (f.args, i))
    end
    while !isempty(stack)
        arr, idx = pop!(stack)
        expr = arr[idx]
        if expr isa JuMP.GenericNonlinearExpr
            for i in eachindex(expr.args) # should be reverse, but order doesn't matter
                push!(stack, (expr.args, i))
            end
        else
            arr[idx] = _remove_variable_from_leaf(expr, vref)
        end
    end
    return
end

# AbstractArray 
function _remove_variable(arr::AbstractArray, vref::GeneralVariableRef)
    for f in arr
        _remove_variable(f, vref)
    end
    return
end

################################################################################
#                                MAPPING METHODS
################################################################################
"""
    map_expression(
        transform::Function, 
        expr::JuMP.AbstractJuMPScalar,
        [write_model::Union{InfiniteModel, AbstractTransformationBackend}]
    )::JuMP.AbstractJuMPScalar

Map and return a new expression of `expr` where each variable is transformed 
via `transform`. This can be helpful for writing user extensions. For this 
to work with new `AbstractTransformationBackend`s, 
[`JuMP.variable_ref_type`](@ref) needs to be defined for the backend. If 
`write_model` is not provided, the `JuMP.owner_model(expr)` is used.
`JuMP.owner_model(expr)` is used.
"""
function map_expression(
    transform::Function,
    v::JuMP.AbstractVariableRef,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    return transform(v)
end

# Constant
function map_expression(
    transform::Function,
    c::Number,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    return c
end

# AffExpr
function map_expression(
    transform::Function,
    aff::JuMP.GenericAffExpr,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    # flatten needed in case expression becomes nonlinear with transform
    return JuMP.flatten!(@_expr(sum(c * transform(v) 
                        for (c, v) in JuMP.linear_terms(aff)) + 
                        JuMP.constant(aff)))
end

# QuadExpr
function map_expression(
    transform::Function,
    quad::JuMP.GenericQuadExpr,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    # flatten needed in case expression becomes nonlinear with transform
    return JuMP.flatten!(@_expr(sum(c * transform(v1) * transform(v2) 
                        for (c, v1, v2) in JuMP.quad_terms(quad)) + 
                        map_expression(transform, quad.aff, write_model)))
end

# NonlinearExpr 
function map_expression(
    transform::Function,
    nlp::JuMP.GenericNonlinearExpr,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    stack = Tuple{Vector{Any}, Vector{Any}}[]
    NewVrefType = JuMP.variable_ref_type(write_model)
    new_nlp = JuMP.GenericNonlinearExpr{NewVrefType}(nlp.head, Any[])
    push!(stack, (nlp.args, new_nlp.args))
    while !isempty(stack)
        args, cloned = pop!(stack)
        for arg in args
            if arg isa JuMP.GenericNonlinearExpr
                new_expr = JuMP.GenericNonlinearExpr{NewVrefType}(arg.head, Any[])
                push!(stack, (arg.args, new_expr.args))
            else
                new_expr = map_expression(transform, arg, write_model)
            end
            push!(cloned, new_expr)
        end
    end
    return new_nlp
end

# Two arguments
function map_expression(transform::Function, expr::JuMP.AbstractJuMPScalar)
    model = JuMP.owner_model(expr)
    isnothing(model) && return expr
    return map_expression(transform, expr, model)
end
function map_expression(transform::Function, c::Number)
    return c
end

################################################################################
#                            EXPRESSION RESTICTIONS
################################################################################
"""
    restrict(expr::JuMP.AbstractJuMPScalar, supps...)::JuMP.AbstractJuMPScalar

Restrict an infinite expression `expr` to be enforced over infinite parameter 
supports `supps`. This is limited to expressions only contain infinite variables 
with the same kind of infinite parameter dependencies. Note that more conveniently 
the expression can be treated as a function for the syntax `expr(supps...)`. 

**Example**
```julia-repl
julia> ex = @expression(model, 3y - 2)
3 y(t) - 2

julia> restrict(ex, 0)
3 y(0) - 2

julia> ex(0)
3 y(0) - 2
```
"""
function restrict(expr::JuMP.AbstractJuMPScalar, supps...)
    # check to make sure all variables depend on the same infinite parameters
    pref_tuples = Set()
    _interrogate_variables(expr) do v
        if v.index_type <: InfiniteParameterIndex
            error("Restrictions on expressions with infinite parameters are not supported.")
        elseif !(v.index_type <: Union{FiniteVariableIndex, PointVariableIndex, FiniteParameterIndex})
            push!(pref_tuples, raw_parameter_refs(v))
        end
    end
    if length(pref_tuples) > 1
        error("Unable to restrict expression `$expr` with supports `$supps`. " *
              "Not all the variables use the same infinite parameters in the " *
              "same format.")
    end
    # restrict the expression using supps and return
    return map_expression(expr) do v
        if isempty(parameter_group_int_indices(v))
            return v
        else
            return restrict(v, supps...)
        end
    end
end

## Make expressions callable for restrictions
# AffExprs
function (aff::JuMP.GenericAffExpr{Float64, GeneralVariableRef})(supps...) 
    return restrict(aff, supps...)
end

# QuadExprs
function (quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef})(supps...) 
    return restrict(quad, supps...)
end

# NonlinearExprs
function (nl::JuMP.GenericNonlinearExpr{GeneralVariableRef})(supps...) 
    return restrict(nl, supps...)
end

################################################################################
#                             COEFFICIENT METHODS
################################################################################
## Modify linear coefficient of variable in expression
# GeneralVariableRef
function _set_variable_coefficient!(expr::GeneralVariableRef,
    var::GeneralVariableRef,
    coeff::Real
    )
    # Determine if variable is that of the expression and change accordingly
    if isequal(expr, var)
        return Float64(coeff) * var
    else
        return expr + Float64(coeff) * var
    end
end

# GenericAffExpr
function _set_variable_coefficient!(expr::JuMP.GenericAffExpr,
    var::GeneralVariableRef,
    coeff::Real
    )
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
    coeff::Real
    )
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

## Query the coefficient of a variable 
# GeneralVariableRef 
function _affine_coefficient(
    func::GeneralVariableRef, 
    var::GeneralVariableRef
    )
    return isequal(func, var) ? 1.0 : 0.0
end

# GenericAffExpr
function _affine_coefficient(
    func::GenericAffExpr, 
    var::GeneralVariableRef
    )
    return get(func.terms, var, 0.0)
end

# GenericQuadExpr
function _affine_coefficient(
    func::GenericQuadExpr, 
    var::GeneralVariableRef
    )
    return get(func.aff.terms, var, 0.0)
end

# Fallback 
function _affine_coefficient(func, var::GeneralVariableRef)
    error("Unsupported function type for coefficient queries.")
end

################################################################################
#                         PARAMETER REFERENCE METHODS
################################################################################
"""
    parameter_refs(expr)::Tuple

Return the tuple of parameter references that determine the infinite
dependencies of `expr`.

**Example**
```julia-repl
julia> parameter_refs(my_expr)
(t,)
```
"""
function parameter_refs(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr}
    )
    model = JuMP.owner_model(expr)
    isnothing(model) && return ()
    prefs = parameter_refs(model)
    group_int_idxs = parameter_group_int_indices(expr)
    length(prefs) == length(group_int_idxs) && return prefs
    return Tuple(prefs[i] for i in eachindex(prefs) if i in group_int_idxs)
end
