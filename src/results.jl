################################################################################
#                                  MODEL QUERIES
################################################################################
# Simple model queries
for op in (:termination_status, :raw_status, :solve_time, :simplex_iterations,
           :barrier_iterations, :node_count, :objective_bound, :relative_gap)
    @eval begin 
        @doc """
            JuMP.$($op)(model::InfiniteModel)

        Extend [`JuMP.$($op)`](https://jump.dev/JuMP.jl/v0.21.8/reference/solutions/#JuMP.$($op)) 
        for `InfiniteModel`s in accordance with that reported by its optimizer 
        model. Errors if such a query is not supported or if the optimizer model 
        hasn't be solved.
        """
        function JuMP.$op(model::InfiniteModel)
            return JuMP.$op(optimizer_model(model))
        end
    end
end

# Simple result dependent model queries
for op in (:primal_status, :dual_status, :has_values, :has_duals, 
           :objective_value, :dual_objective_value)
    @eval begin 
        @doc """
            JuMP.$($op)(model::InfiniteModel; [result::Int = 1])

        Extend [`JuMP.$($op)`](https://jump.dev/JuMP.jl/v0.21.8/reference/solutions/#JuMP.$($op)) 
        for `InfiniteModel`s in accordance with that reported by its optimizer 
        model and the result index `result` of the most recent solution obtained. 
        Errors if such a query is not supported or if the optimizer model hasn't 
        be solved.
        """
        function JuMP.$op(model::InfiniteModel; result::Int = 1)
            return JuMP.$op(optimizer_model(model); result = result)
        end
    end
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
function _get_value(pref, ::Type{FiniteParameterIndex}, result; kwargs...)
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

Extend `JuMP.value` to return the value(s) of `vref` in accordance with its 
reformulation variable(s) stored in the optimizer model and the result index 
`result` of the most recent solution obtained. Use
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

Extend `JuMP.value` to return the value(s) of `cref` in accordance with its 
reformulation constraint(s) stored in the optimizer model and the result index 
`result` of the most recent solution obtained. Use 
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check if a result 
exists before asking for values. 
    
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

Extend `JuMP.reduced_cost`. This returns the reduced cost(s) of a variable. This 
will be a vector of scalar values for an infinite variable or will be a scalar 
value for finite variables. 

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

Extend `JuMP.optimizer_index` to return the `MathOptInterface` index(es) of 
`vref` in accordance with its reformulation variable(s) stored in the optimizer 
model.

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

Extend `JuMP.optimizer_index` to return the `MathOptInterface` index(es) of 
`cref` in accordance with its reformulation constraints(s) stored in the 
optimizer model. 

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

Extend `JuMP.dual` to return the dual(s) of `cref` in accordance with its 
reformulation constraint(s) stored in the optimizer model and the result index 
`result` of the most recent solution obtained. Use 
[`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel)) to check if a result 
exists before asking for duals. 

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

Extend `JuMP.shadow_price` to return the shadow price(s) of `cref` in accordance 
with its reformulation constraint(s) stored in the optimizer model. Use 
[`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel)) to check if a result 
exists before asking for duals. 
    
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

A wrapper `DataType` for `JuMP.SensitivityReport`s in `InfiniteOpt`. 
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

Extends `JuMP.lp_sensitivity_report` to generate and return an LP sensitivity 
report in accordance with the optimizer model. See 
[`InfOptSensitivityReport`](@ref) for syntax details on how to query it. `atol` 
denotes the optimality tolerance and should match that used by the solver to 
compute the basis. Please refer to `JuMP`'s documentation for more technical 
information on interpretting the output of the report.

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
