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

"""
    map_value(ref, key::Val{ext_key_name}, result::Int)

Map the value(s) of `ref` to its counterpart in the optimizer model type that is
distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default. Here `result` is the result index used in `value`.
"""
function map_value end

# Default method that depends on optimizer_model_variable --> making extensions easier
function map_value(vref::GeneralVariableRef, key, result::Int)
    opt_vref = optimizer_model_variable(vref, key)
    if opt_vref isa AbstractArray
        return JuMP.value.(opt_vref; result = result)
    else
        return JuMP.value(opt_vref; result = result)
    end
end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_value(cref::GeneralConstraintRef, key, result::Int)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.value.(opt_cref; result = result)
    else
        return JuMP.value(opt_cref; result = result)
    end
end

"""
    JuMP.value(vref::InfOptVariableRef; [result::Int = 1])

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.VariableRef)) to return the value(s)
of `vref` in accordance with its reformulation variable(s) stored in the optimizer
model and the result index `result` of the most recent solution obtained. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
if a result exists before asking for values. For extensions, this only works if
[`optimizer_model_variable`](@ref) has been extended correctly and/or
[`map_value`](@ref) has been extended for variables.

**Example**
```julia-repl
julia> value(z)
42.0
```
"""
function JuMP.value(vref::InfOptVariableRef; result::Int = 1)
    return map_value(vref, Val(optimizer_model_key(JuMP.owner_model(vref))), result)
end

"""
    JuMP.value(cref::GeneralConstraintRef; [result::Int = 1])

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the value(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model and the result index `result` of the most recent
solution obtained. Use [`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel))
to check if a result exists before asking for values. For extensions, this only
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
function JuMP.value(cref::GeneralConstraintRef; result::Int = 1)
    return map_value(cref, Val(optimizer_model_key(JuMP.owner_model(cref))), result)
end

"""
    map_optimizer_index(ref, key::Val{ext_key_name})

Map the `MathOptInterface` index(es) of `ref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default.
"""
function map_optimizer_index end

# Default method that depends on optimizer_model_variable --> making extensions easier
function map_optimizer_index(vref::GeneralVariableRef, key)
    opt_vref = optimizer_model_variable(vref, key)
    if opt_vref isa AbstractArray
        return JuMP.optimizer_index.(opt_vref)
    else
        return JuMP.optimizer_index(opt_vref)
    end
end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_optimizer_index(cref::GeneralConstraintRef, key)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.optimizer_index.(opt_cref)
    else
        return JuMP.optimizer_index(opt_cref)
    end
end

"""
    JuMP.optimizer_index(vref::InfOptVariableRef)

Extend [`JuMP.optimizer_index`](@ref JuMP.optimizer_index(::JuMP.VariableRef)) to
return the `MathOptInterface` index(es) of `vref` in accordance with its
reformulation variable(s) stored in the optimizer model. For extensions, this
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
function JuMP.optimizer_index(vref::InfOptVariableRef)
    return map_optimizer_index(vref, Val(optimizer_model_key(JuMP.owner_model(vref))))
end

"""
    JuMP.optimizer_index(cref::GeneralConstraintRef)

Extend [`JuMP.optimizer_index`](@ref JuMP.optimizer_index(::JuMP.ConstraintRef{JuMP.Model}))
to return the `MathOptInterface` index(es) of `cref` in accordance with its
reformulation constraints(s) stored in the optimizer model. For extensions, this
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
function JuMP.optimizer_index(cref::GeneralConstraintRef)
    return map_optimizer_index(cref, Val(optimizer_model_key(JuMP.owner_model(cref))))
end

"""
    map_dual(cref::GeneralConstraintRef, key::Val{ext_key_name}, result::Int)

Map the dual(s) of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default. Here `result` is the result index that is used in `dual`.
"""
function map_dual end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_dual(cref::GeneralConstraintRef, key, result::Int)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.dual.(opt_cref; result = result)
    else
        return JuMP.dual(opt_cref; result = result)
    end
end

"""
    JuMP.dual(cref::GeneralConstraintRef; [result::Int = 1])

Extend [`JuMP.dual`](@ref JuMP.dual(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the dual(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model and the result index `result` of the most recent
solution obtained. Use [`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel))
to check if a result exists before asking for duals. For extensions, this only
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
function JuMP.dual(cref::GeneralConstraintRef; result::Int = 1)
    return map_dual(cref, Val(optimizer_model_key(JuMP.owner_model(cref))),
                    result)
end

# Error redriect dor variable call
function JuMP.dual(vref::GeneralVariableRef; result::Int = 1)
    return JuMP.dual(JuMP.VariableRef(JuMP.Model(), MOI.VariableIndex(1)))
end

"""
    map_shadow_price(cref::GeneralConstraintRef, key::Val{ext_key_name})

Map the shadow price(s) of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default.
"""
function map_shadow_price end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_shadow_price(cref::GeneralConstraintRef, key)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.shadow_price.(opt_cref)
    else
        return JuMP.shadow_price(opt_cref)
    end
end

"""
    JuMP.shadow_price(cref::GeneralConstraintRef)

Extend [`JuMP.shadow_price`](@ref JuMP.shadow_price(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the shadow price(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model. Use [`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel))
to check if a result exists before asking for duals. For extensions, this only
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
function JuMP.shadow_price(cref::GeneralConstraintRef)
    return map_shadow_price(cref, Val(optimizer_model_key(JuMP.owner_model(cref))))
end

"""
    map_lp_rhs_perturbation_range(cref::GeneralConstraintRef,
                                  key::Val{ext_key_name}, toler::Float64)

Map the RHS perturbation range of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `cref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_constraint` is used to make these mappings by default. Here
`toler` corresponds to the `feasibility_tolerance` used by `lp_rhs_perturbation_range`.
"""
function map_lp_rhs_perturbation_range end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_lp_rhs_perturbation_range(cref::GeneralConstraintRef, key,
                                       toler::Float64)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.lp_rhs_perturbation_range.(opt_cref;
                                               feasibility_tolerance = toler)
    else
        return JuMP.lp_rhs_perturbation_range(opt_cref;
                                              feasibility_tolerance = toler)
    end
end

"""
    JuMP.lp_rhs_perturbation_range(cref::GeneralConstraintRef;
                                   [feasibility_tolerance::Float64 = 1e-8])

Extend [`JuMP.lp_rhs_perturbation_range`](@ref JuMP.lp_rhs_perturbation_range(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the range(s) of the RHS of `cref` for which the shadow price(s) are valid
in accordance with its reformulation constraint(s)
stored in the optimizer model. For extensions, this only
works if [`optimizer_model_constraint`](@ref) has been extended correctly and/or
[`map_lp_rhs_perturbation_range`](@ref) has been implemented.

**Example**
```julia-repl
julia> lp_rhs_perturbation_range(c1)
4-element Array{Tuple{Float64,Float64},1}:
 (-42.0, Inf)
 (-Inf, 42.0)
 (-Inf, 42.0)
 (-Inf, 42.0)
```
"""
function JuMP.lp_rhs_perturbation_range(cref::GeneralConstraintRef;
                                        feasibility_tolerance::Float64 = 1e-8)
    return map_lp_rhs_perturbation_range(cref, Val(optimizer_model_key(JuMP.owner_model(cref))),
                                         feasibility_tolerance)
end

"""
    map_lp_objective_perturbation_range(vref::InfOptVariableRef,
                                        key::Val{ext_key_name}, toler::Float64)

Map the reduced cost range(s) of `vref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `vref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` is used to make these mappings by default. Here
`toler` corresponds to the `optimality_tolerance` used by
`lp_objective_perturbation_range`.
"""
function map_lp_objective_perturbation_range end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_lp_objective_perturbation_range(vref::GeneralVariableRef, key,
                                             toler::Float64)
    opt_vref = optimizer_model_variable(vref)
    if opt_vref isa AbstractArray
        return JuMP.lp_objective_perturbation_range.(opt_vref;
                                                     optimality_tolerance = toler)
    else
        return JuMP.lp_objective_perturbation_range(opt_vref;
                                                    optimality_tolerance = toler)
    end
end

"""
    JuMP.lp_objective_perturbation_range(vref::InfOptVariableRef;
                                         [optimality_tolerance::Float64 = 1e-8])

Extend [`JuMP.lp_objective_perturbation_range`](@ref JuMP.lp_objective_perturbation_range(::JuMP.VariableRef))
to return the range(s) that the reduced cost(s) of `vref` remain valid
in accordance with its reformulation variables(s)
stored in the optimizer model. For extensions, this only
works if [`optimizer_model_variable`](@ref) has been extended correctly and/or
if [`map_lp_objective_perturbation_range`](@ref) has been implemented.

**Example**
```julia-repl
julia> lp_objective_perturbation_range(z)
(-2.0, Inf)
```
"""
function JuMP.lp_objective_perturbation_range(vref::GeneralVariableRef;
                                              optimality_tolerance::Float64 = 1e-8)
    return map_lp_objective_perturbation_range(vref, Val(optimizer_model_key(JuMP.owner_model(vref))),
                                               optimality_tolerance)
end
