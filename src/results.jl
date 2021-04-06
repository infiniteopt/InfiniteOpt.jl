################################################################################
#                                  STATUS QUERIES
################################################################################
"""
    JuMP.termination_status(model::InfiniteModel)

Extend [`termination_status`](@ref JuMP.termination_status(::JuMP.Model)) to
return the `MOI.TerminationStatus` in accordance with the optimizer model.

**Example**
```julia-repl
julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4
```
"""
function JuMP.termination_status(model::InfiniteModel)
    return JuMP.termination_status(optimizer_model(model))
end

"""
    JuMP.raw_status(model::InfiniteModel)

Extend [`raw_status`](@ref JuMP.raw_status(::JuMP.Model)) to return the status
reported by the solver in accordance with the optimizer model.

**Example**
```julia-repl
julia> raw_status(model) # Ipopt
"Solve_Succeeded"
```
"""
function JuMP.raw_status(model::InfiniteModel)
    return JuMP.raw_status(optimizer_model(model))
end

"""
     JuMP.primal_status(model::InfiniteModel; [result::Int = 1])

Extend [`primal_status`](@ref JuMP.primal_status(::JuMP.Model)) to return the
`MOI.PrimalStatus` reported in accordance with the optimizer model and the
result index `result` of the most recent solution obtained.

**Example**
```julia-repl
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
"""
function  JuMP.primal_status(model::InfiniteModel; result::Int = 1)
    return JuMP.primal_status(optimizer_model(model); result = result)
end

"""
    JuMP.dual_status(model::InfiniteModel; [result::Int = 1])

Extend [`dual_status`](@ref JuMP.dual_status(::JuMP.Model)) to return the
`MOI.DualStatus` reported in accordance with the optimizer model and the
result index `result` of the most recent solution obtained.

**Example**
```julia-repl
julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
"""
function JuMP.dual_status(model::InfiniteModel; result::Int = 1)
    return JuMP.dual_status(optimizer_model(model); result = result)
end

"""
    JuMP.solve_time(model::InfiniteModel)

Extend [`solve_time`](@ref JuMP.solve_time(::JuMP.Model)) to return the
time used by the solver to terminate reported in accordance with the optimizer
model. This will error if not supported by the solver.

**Example**
```julia-repl
julia> solve_time(model)
0.004999876022338867
```
"""
function JuMP.solve_time(model::InfiniteModel)
    return JuMP.solve_time(optimizer_model(model))
end

"""
    JuMP.has_values(model::InfiniteModel; [result::Int = 1])::Bool

Extend [`has_values`](@ref JuMP.has_values(::JuMP.Model)) to return a `Bool`
whether variable values are available in accordance with the optimizer model and
the result index `result` of the most recent solution obtained.

**Example**
```julia-repl
julia> has_values(model)
true
```
"""
function JuMP.has_values(model::InfiniteModel; result::Int = 1)::Bool
    return JuMP.primal_status(model; result = result) != MOI.NO_SOLUTION
end

"""
    JuMP.has_duals(model::InfiniteModel; [result::Int = 1])::Bool

Extend [`has_duals`](@ref JuMP.has_duals(::JuMP.Model)) to return a `Bool`
whether constraint duals are available in accordance with the optimizer model and
the result index `result` of the most recent solution obtained.

**Example**
```julia-repl
julia> has_duals(model)
true
```
"""
function JuMP.has_duals(model::InfiniteModel; result::Int = 1)::Bool
    return JuMP.dual_status(model; result = result) != MOI.NO_SOLUTION
end

"""
    simplex_iterations(model::InfiniteModel)

Gets the cumulative number of simplex iterations during the most-recent optimization.

Solvers must implement `MOI.SimplexIterations()` to use this function.
"""
function JuMP.simplex_iterations(model::InfiniteModel)
    return JuMP.simplex_iterations(optimizer_model(model))
end

"""
    barrier_iterations(model::InfiniteModel)

Gets the cumulative number of barrier iterations during the most-recent optimization.

Solvers must implement `MOI.BarrierIterations()` to use this function.
"""
function JuMP.barrier_iterations(model::InfiniteModel)
    return JuMP.barrier_iterations(optimizer_model(model))
end

"""
    node_count(model::InfiniteModel)

Gets the total number of branch-and-bound nodes explored during the most recent
optimization in a Mixed Integer Program.

Solvers must implement `MOI.NodeCount()` to use this function.
"""
function JuMP.node_count(model::InfiniteModel)
    return JuMP.node_count(optimizer_model(model))
end

################################################################################
#                              OBJECTIVE QUERIES
################################################################################
"""
    JuMP.objective_bound(model::InfiniteModel)::Float64

Extend [`objective_bound`](@ref JuMP.objective_bound(::JuMP.Model)) to return
the objective bound in accordance with the optimizer model.

**Example**
```julia-repl
julia> objective_bound(model)
42.0
```
"""
function JuMP.objective_bound(model::InfiniteModel)::Float64
    return JuMP.objective_bound(optimizer_model(model))
end

"""
    JuMP.objective_value(model::InfiniteModel; [result::Int = 1])::Float64

Extend [`objective_value`](@ref JuMP.objective_value(::JuMP.Model)) to return
the objective value in accordance with the optimizer model and the
result index `result` of the most recent solution obtained.

**Example**
```julia-repl
julia> objective_value(model)
42.0
```
"""
function JuMP.objective_value(model::InfiniteModel; result::Int = 1)::Float64
    return JuMP.objective_value(optimizer_model(model); result = result)
end

"""
    JuMP.dual_objective_value(model::InfiniteModel; [result::Int = 1])::Float64

Extend [`dual_objective_value`](@ref JuMP.dual_objective_value(::JuMP.Model)) to
return the dual objective value in accordance with the optimizer model and the
result index `result` of the most recent solution obtained. Errors
if the solver does not support this.

**Example**
```julia-repl
julia> dual_objective_value(model)
42.0
```
"""
function JuMP.dual_objective_value(model::InfiniteModel; result::Int = 1)::Float64
    return JuMP.dual_objective_value(optimizer_model(model); result = result)
end

################################################################################
#                                 VALUE QUERIES
################################################################################
"""
    map_value([ref/expr], key::Val{ext_key_name}, result::Int; kwargs...)

Map the value(s) of `ref` to its counterpart in the optimizer model type that is
distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable`, `optimizer_model_expression`, and/or
`optimizer_model_constraint`. Such as is the case with reformuations that do not
have a direct mapping between variables and/or constraints in the original
infinite form. Otherwise, `optimizer_model_variable`,
`optimizer_model_expression`, and `optimizer_model_constraint` are used to make
these mappings by default where `kwargs` are passed on these functions. Here 
`result` is the result index used in `value`.
"""
function map_value end

# Default method that depends on optimizer_model_variable --> making extensions easier
function map_value(vref::GeneralVariableRef, key, result::Int; kwargs...)
    opt_vref = optimizer_model_variable(vref, key; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> JuMP.value(v; result = result), opt_vref)
    else
        return JuMP.value(opt_vref; result = result)
    end
end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_value(cref::InfOptConstraintRef, key, result::Int; kwargs...)
    opt_cref = optimizer_model_constraint(cref, key; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.value(c; result = result), opt_cref)
    else
        return JuMP.value(opt_cref; result = result)
    end
end

# Default method that depends on optimizer_model_expression --> making extensions easier
function map_value(expr::JuMP.AbstractJuMPScalar, key, result::Int; kwargs...)
    opt_expr = optimizer_model_expression(expr, key; kwargs...)
    if opt_expr isa AbstractArray
        return map(e -> JuMP.value(e; result = result), opt_expr)
    else
        return JuMP.value(opt_expr; result = result)
    end
end

## Define dispatch methods to collect value of parameters 
# InfiniteParameter 
function _get_value(pref, ::Type{<:InfiniteParameterIndex}, result; kwargs...)
    label = get(kwargs, :label, PublicLabel)
    return supports(pref, label = label)
end

# FiniteParameter 
function _get_value(pref, ::Type{<:FiniteParameterIndex}, result; kwargs...)
    return parameter_value(pref)
end

# Others 
function _get_value(vref, index_type, result; kwargs...)
    return map_value(vref, Val(optimizer_model_key(JuMP.owner_model(vref))), 
                     result; kwargs...)
end

"""
    JuMP.value(vref::GeneralVariableRef; [result::Int = 1, 
               label::Type{<:AbstractSupportLabel} = PublicLabel,
               ndarray::Bool = false, kwargs...])

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.VariableRef)) to return the value(s)
of `vref` in accordance with its reformulation variable(s) stored in the optimizer
model and the result index `result` of the most recent solution obtained. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
if a result exists before asking for values. 
    
The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the values associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the values of infinite 
variables are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the variable has multiple 
infinite parameter dependencies.

To provide context for the
results it may be helpful to also query the variable's `parameter_refs` and
`supports` which will have a one-to-one correspondence with the value(s).
It may also be helpful to query via [`optimizer_model_variable`](@ref) to
retrieve the variables(s) that these values are based on. These functions should 
all be called with the same keyword arugments for consistency.

For extensions, this only works if
[`optimizer_model_variable`](@ref) has been extended correctly and/or
[`map_value`](@ref) has been extended for variables.

**Example**
```julia-repl
julia> value(z)
42.0
```
"""
function JuMP.value(vref::GeneralVariableRef; result::Int = 1, kwargs...)
    return _get_value(vref, _index_type(vref), result; kwargs...)
end

"""
    JuMP.value(cref::InfOptConstraintRef; [result::Int = 1,
               label::Type{<:AbstractSupportLabel} = PublicLabel,
               ndarray::Bool = false, kwargs...])

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the value(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model and the result index `result` of the most recent
solution obtained. Use [`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel))
to check if a result exists before asking for values. 
    
The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the values associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the values of infinite 
constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

To provide context for
the results it may be helpful to also query the constraint's `parameter_refs`
and `supports` which will have a one-to-one correspondence with the value(s).
It may also be helpful to query via [`optimizer_model_constraint`](@ref) to
retrieve the constraint(s) that these values are based on. By default, only the 
values corresponding to public supports are returned. These functions should 
all be called with the same keyword arugments for consistency.

For extensions, this only
works if [`optimizer_model_constraint`](@ref) has been extended correctly and/or
[`map_value`](@ref) has been extended for constraints. 

**Example**
```julia-repl
julia> value(c1)
4-element Array{Float64,1}:
 -0.0
 20.9
 20.9
 20.9
```
"""
function JuMP.value(cref::InfOptConstraintRef; result::Int = 1, 
                    kwargs...)
    return map_value(cref, Val(optimizer_model_key(JuMP.owner_model(cref))), 
                     result; kwargs...)
end

"""
    JuMP.value(expr::JuMP.AbstractJuMPScalar; [result::Int = 1, 
               label::Type{<:AbstractSupportLabel} = PublicLabel,
               ndarray::Bool = false, kwargs...])

Return the value(s) of `expr` in accordance with the optimized variable values
the result index `result` of the most recent solution obtained. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check if a result
exists before asking for values. 
    
The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the values associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the values of infinite 
expressions are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the expression has multiple 
infinite parameter dependencies.
    
To provide context for the results it may be
helpful to also query the expression's `parameter_refs` and `supports` which
will have a one-to-one correspondence with the value(s). It may also be helpful
to query via [`optimizer_model_expression`](@ref) to retrieve the expression(s)
that these values are based on. These should use the same keyword arguments for 
consistency.

For extensions, this only works if [`optimizer_model_expression`](@ref) has been
extended correctly and/or [`map_value`](@ref) has been extended for expressions.

**Example**
```julia-repl
julia> value(my_finite_expr)
23.34

julia> value(my_infinite_expr)
4-element Array{Float64,1}:
 -0.0
 20.9
 20.9
 20.9
```
"""
function JuMP.value(expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}};
    result::Int = 1,
    kwargs...
    ) where {C, V <: GeneralVariableRef}
    # get the model
    model = _model_from_expr(expr)
    # if no model then the expression only contains a constant
    if model === nothing
        return JuMP.constant(expr)
    # otherwise let's call map_value
    else
        key = optimizer_model_key(model)
        return map_value(expr, Val(key), result; kwargs...)
    end
end

################################################################################
#                                 REDUCED COST
################################################################################
"""
    map_reduced_cost(vref::GeneralVariableRef, key::Val{ext_key_name}, 
                      result::Int; kwargs...)

Map the reduced cost(s) of `vref` to its counterpart in the optimizer model type that is
distininguished by its extension key `key` as type `Val{ext_key_name}`.
This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable`. Such as is the case with reformulations 
that do not have a direct mapping between variables in the original
infinite form. Otherwise, `optimizer_model_variable`, is used to make
these mappings by default where `kwargs` are passed on these functions. Here 
`result` is the result index used in `value`.
"""
function map_reduced_cost end

# Default definition for when optimizer_model_variable is defined
function map_reduced_cost(vref::GeneralVariableRef, key; kwargs...)
    opt_vref = optimizer_model_variable(vref, key; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> JuMP.reduced_cost(v), opt_vref)
    else
        return JuMP.reduced_cost(opt_vref)
    end
end

"""
    JuMP.reduced_cost(vref::GeneralVariableRef)

Extending [`JuMP.reduced_cost`](@ref JuMP.reduced_cost(::JuMP.VariableRef)). This 
returns the reduced cost(s) of a variable. This will be a vector of scalar values 
for an infinite variable or will be a scalar value for finite variables. 

**Example**
```julia-repl
julia> reduced_cost(x)
12.81
```
"""
function JuMP.reduced_cost(vref::GeneralVariableRef; kwargs...)
    return map_reduced_cost(vref, Val(optimizer_model_key(JuMP.owner_model(vref))); 
                            kwargs...)
end

################################################################################
#                             OPTIMIZER INDEX QUERIES
################################################################################
"""
    map_optimizer_index(ref, key::Val{ext_key_name}; kwargs...)

Map the `MathOptInterface` index(es) of `ref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default where `kwargs` are passed on as well.
"""
function map_optimizer_index end

# Default method that depends on optimizer_model_variable --> making extensions easier
function map_optimizer_index(vref::GeneralVariableRef, key; kwargs...)
    opt_vref = optimizer_model_variable(vref, key; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> JuMP.optimizer_index(v), opt_vref)
    else
        return JuMP.optimizer_index(opt_vref)
    end
end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_optimizer_index(cref::InfOptConstraintRef, key; kwargs...)
    opt_cref = optimizer_model_constraint(cref, key; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.optimizer_index(c), opt_cref)
    else
        return JuMP.optimizer_index(opt_cref)
    end
end

"""
    JuMP.optimizer_index(vref::GeneralVariableRef; 
                         [label::Type{<:AbstractSupportLabel} = PublicLabel,
                         ndarray::Bool = false, kwargs...])

Extend [`JuMP.optimizer_index`](@ref JuMP.optimizer_index(::JuMP.VariableRef)) to
return the `MathOptInterface` index(es) of `vref` in accordance with its
reformulation variable(s) stored in the optimizer model.

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the optimizer indices associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the indices of infinite 
variables are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the variable has multiple 
infinite parameter dependencies.

It may also be helpful to query via [`optimizer_model_variable`](@ref) to
retrieve the variables(s) that these indices are based on. These should use the 
same keyword arguments for consistency.

For extensions, this
only works if [`optimizer_model_variable`](@ref) has been extended correctly
and/or [`map_optimizer_index`](@ref) has been extended for variables.

**Example**
```julia-repl
julia> optimizer_index(x)
4-element Array{MathOptInterface.VariableIndex,1}:
 MathOptInterface.VariableIndex(2)
 MathOptInterface.VariableIndex(3)
 MathOptInterface.VariableIndex(4)
 MathOptInterface.VariableIndex(5)
```
"""
function JuMP.optimizer_index(vref::GeneralVariableRef; kwargs...)
    return map_optimizer_index(vref, Val(optimizer_model_key(JuMP.owner_model(vref))); 
                               kwargs...)
end

"""
    JuMP.optimizer_index(cref::InfOptConstraintRef; 
                         [label::Type{<:AbstractSupportLabel} = PublicLabel,
                         ndarray::Bool = false, kwargs...])

Extend [`JuMP.optimizer_index`](@ref JuMP.optimizer_index(::JuMP.ConstraintRef{JuMP.Model}))
to return the `MathOptInterface` index(es) of `cref` in accordance with its
reformulation constraints(s) stored in the optimizer model. 

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the optimizer indices associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the indices of infinite 
constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

It may also be helpful to query via [`optimizer_model_constraint`](@ref) to
retrieve the constraints(s) that these indices are based on. The same keyword 
arguments should be used for consistency.

For extensions, this
only works if [`optimizer_model_constraint`](@ref) has been extended correctly
and/or [`map_optimizer_index`](@ref) has been extended for constraints.

**Example**
```julia-repl
julia> optimizer_index(c1)
4-element Array{MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}},1}:
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(1)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(2)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(3)
 MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.GreaterThan{Float64}}(4)
```
"""
function JuMP.optimizer_index(cref::InfOptConstraintRef; kwargs...)
    return map_optimizer_index(cref, Val(optimizer_model_key(JuMP.owner_model(cref))); 
                               kwargs...)
end

################################################################################
#                                  DUAL QUERIES
################################################################################
"""
    map_dual(cref::InfOptConstraintRef, key::Val{ext_key_name}, result::Int; 
             kwargs...)

Map the dual(s) of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default where `kwargs` are also pass on to. Here `result` is 
the result index that is used in `dual`. 
"""
function map_dual end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_dual(cref::InfOptConstraintRef, key, result::Int; kwargs...)
    opt_cref = optimizer_model_constraint(cref, key; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.dual(c; result = result), opt_cref)
    else
        return JuMP.dual(opt_cref; result = result)
    end
end

"""
    JuMP.dual(cref::InfOptConstraintRef; [result::Int = 1, 
              label::Type{<:AbstractSupportLabel} = PublicLabel,
              ndarray::Bool = false, kwargs...])

Extend [`JuMP.dual`](@ref JuMP.dual(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the dual(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model and the result index `result` of the most recent
solution obtained. Use [`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel))
to check if a result exists before asking for duals. 

The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the duals associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the duals of infinite 
constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

It may also be helpful to
query via [`optimizer_model_constraint`](@ref) to retrieve the constraint(s)
that these duals are based on. Calling `parameter_refs` and `supports` may also
be insightful. Be sure to use the same keyword arguments for consistency.

For extensions, this only
works if [`optimizer_model_constraint`](@ref) has been extended correctly and/or
[`map_dual`](@ref) has been extended for constraints.

**Example**
```julia-repl
julia> dual(c1)
4-element Array{Float64,1}:
 -42.0
 -42.0
 32.3
 0.0
```
"""
function JuMP.dual(cref::InfOptConstraintRef; result::Int = 1, kwargs...)
    return map_dual(cref, Val(optimizer_model_key(JuMP.owner_model(cref))),
                    result; kwargs...)
end

# Error redriect for variable call
function JuMP.dual(vref::GeneralVariableRef; kwargs...)
    return JuMP.dual(JuMP.VariableRef(JuMP.Model(), MOI.VariableIndex(1)))
end

################################################################################
#                              SHADOW PRICE QUERIES
################################################################################
"""
    map_shadow_price(cref::InfOptConstraintRef, key::Val{ext_key_name}; 
                     kwargs...)

Map the shadow price(s) of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default where `kwargs` are passed on to.
"""
function map_shadow_price end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_shadow_price(cref::InfOptConstraintRef, key; kwargs...)
    opt_cref = optimizer_model_constraint(cref, key; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.shadow_price(c), opt_cref)
    else
        return JuMP.shadow_price(opt_cref)
    end
end

"""
    JuMP.shadow_price(cref::InfOptConstraintRef; 
                      [label::Type{<:AbstractSupportLabel} = PublicLabel,
                      ndarray::Bool = false, kwargs...])

Extend [`JuMP.shadow_price`](@ref JuMP.shadow_price(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the shadow price(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model. Use [`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel))
to check if a result exists before asking for duals. 
    
The keyword arugments `label` and `ndarray` are what `TranscriptionOpt` employ 
and `kwargs` denote extra ones that user extensions may employ.

By default only the shadow prices associated with public supports are returned, the 
full set can be accessed via `label = All`. Moreover, the prices of infinite 
constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

It may also be
helpful to query via [`optimizer_model_constraint`](@ref) to retrieve the
constraint(s) that these shadow prices are based on. Calling `parameter_refs` and
`supports` may also be insightful. Be sure to use the same keyword arguments for 
consistency.

For extensions, this only
works if [`optimizer_model_constraint`](@ref) has been extended correctly and/or
[`map_shadow_price`](@ref) has been extended for constraints. 

**Example**
```julia-repl
julia> shadow_price(c1)
4-element Array{Float64,1}:
 42.0
 42.0
 -32.3
 -0.0
```
"""
function JuMP.shadow_price(cref::InfOptConstraintRef; kwargs...)
    return map_shadow_price(cref, Val(optimizer_model_key(JuMP.owner_model(cref))); 
                            kwargs...)
end

################################################################################
#                           LP SENSITIVITY ANALYSIS
################################################################################
"""
    InfOptSensitivityReport

A wrapper `DataType` for [`JuMP.SensitivityReport`](@ref)s in `InfiniteOpt`. 
These are generated based on the optimizer model and should be made via the use of 
[`lp_sensitivity_report`](@ref JuMP.lp_sensitivity_report(::InfiniteModel)). Once 
made these can be indexed to get the sensitivies with respect to variables and/or 
constraints. The indexing syntax for these is: 
```julia
report[ref::[GeneralVariableRef/InfOptConstraintRef]; 
       [label::Type{<:AbstractSupportLabel} = PublicLabel,
       ndarray::Bool = false, kwargs...]]
```

This is enabled in user-defined optimizer model extensions by appropriately 
extending [`optimizer_model_variable`](@ref) and [`optimizer_model_constraint`](@ref).

**Fields**
- `opt_report::JuMP.SensitivityReport`: The LP sensitivity captured from the optimizer model.
"""
struct InfOptSensitivityReport 
    opt_report::JuMP.SensitivityReport
end

# Extend Base.getindex for variables on InfOptSensitivityReport
function Base.getindex(s::InfOptSensitivityReport, v::GeneralVariableRef; kwargs...)
    key = Val(optimizer_model_key(JuMP.owner_model(v)))
    opt_vref = optimizer_model_variable(v, key; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> s.opt_report[v], opt_vref)
    else
        return s.opt_report[opt_vref]
    end
end

# Extend Base.getindex for constraints on InfOptSensitivityReport
function Base.getindex(s::InfOptSensitivityReport, c::InfOptConstraintRef; kwargs...)
    key = Val(optimizer_model_key(JuMP.owner_model(c)))
    opt_cref = optimizer_model_constraint(c, key; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> s.opt_report[c], opt_cref)
    else
        return s.opt_report[opt_cref]
    end
end

"""
    JuMP.lp_sensitivity_report(model::InfiniteModel; 
                               [atol::Float64 = 1e-8])::InfOptSensitivityReport

Extends [`JuMP.lp_sensitivity_report`](@ref JuMP.lp_sensitivity_report(::JuMP.Model))
to generate and return an LP sensitivity report in accordance with the optimizer 
model. See [`InfOptSensitivityReport`](@ref) for syntax details on how to query 
it. `atol` denotes the optimality tolerance and should match that used by the 
solver to compute the basis. Please refer to `JuMP`'s documentation for more 
technical information on interpretting the output of the report.

**Example**
```julia-repl 
julia> report = lp_sensitivity_report(model);

julia> report[x]
(0.0, 0.5)
```
"""
function JuMP.lp_sensitivity_report(
    model::InfiniteModel; 
    atol::Float64 = 1e-8
    )::InfOptSensitivityReport
    opt_model = optimizer_model(model)
    opt_report = JuMP.lp_sensitivity_report(opt_model, atol = atol)
    return InfOptSensitivityReport(opt_report)
end
