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
     JuMP.primal_status(model::InfiniteModel)

Extend [`primal_status`](@ref JuMP.primal_status(::JuMP.Model)) to return the
`MOI.PrimalStatus` reported in accordance with the optimizer model.

**Example**
```julia-repl
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
"""
function  JuMP.primal_status(model::InfiniteModel)
    return JuMP.primal_status(optimizer_model(model))
end

"""
    JuMP.dual_status(model::InfiniteModel)

Extend [`dual_status`](@ref JuMP.dual_status(::JuMP.Model)) to return the
`MOI.DualStatus` reported in accordance with the optimizer model.

**Example**
```julia-repl
julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
"""
function JuMP.dual_status(model::InfiniteModel)
    return JuMP.dual_status(optimizer_model(model))
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
    JuMP.has_values(model::InfiniteModel)

Extend [`has_values`](@ref JuMP.has_values(::JuMP.Model)) to return a `Bool`
whether variable values are available in accordance with the optimizer model.

**Example**
```julia-repl
julia> has_values(model)
true
```
"""
JuMP.has_values(model::InfiniteModel) = JuMP.primal_status(model) != MOI.NO_SOLUTION

"""
    JuMP.has_duals(model::InfiniteModel)

Extend [`has_duals`](@ref JuMP.has_duals(::JuMP.Model)) to return a `Bool`
whether constraint duals are available in accordance with the optimizer model.

**Example**
```julia-repl
julia> has_duals(model)
true
```
"""
JuMP.has_duals(model::InfiniteModel) = JuMP.dual_status(model) != MOI.NO_SOLUTION

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
    JuMP.objective_value(model::InfiniteModel)::Float64

Extend [`objective_value`](@ref JuMP.objective_value(::JuMP.Model)) to return
the objective value in accordance with the optimizer model.

**Example**
```julia-repl
julia> objective_value(model)
42.0
```
"""
function JuMP.objective_value(model::InfiniteModel)::Float64
    return JuMP.objective_value(optimizer_model(model))
end

"""
    JuMP.dual_objective_value(model::InfiniteModel)::Float64

Extend [`dual_objective_value`](@ref JuMP.dual_objective_value(::JuMP.Model)) to
return the dual objective value in accordance with the optimizer model. Errors
if the solver does not support this.

**Example**
```julia-repl
julia> dual_objective_value(model)
42.0
```
"""
function JuMP.dual_objective_value(model::InfiniteModel)::Float64
    return JuMP.dual_objective_value(optimizer_model(model))
end

"""
    map_value(ref, key)

Map the value(s) of `ref` to its counterpart in the optimizer model type that is
distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default.
"""
function map_value end

# Default method that depends on optimizer_model_variable --> making extensions easier
function map_value(vref::GeneralVariableRef, key)
    opt_vref = optimizer_model_variable(vref, key)
    if opt_vref isa AbstractArray
        return JuMP.value.(opt_vref)
    else
        return JuMP.value(opt_vref)
    end
end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_value(cref::GeneralConstraintRef, key)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.value.(opt_cref)
    else
        return JuMP.value(opt_cref)
    end
end

"""
    JuMP.value(vref::InfOptVariableRef)

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.VariableRef)) to return the value(s)
of `vref` in accordance with its reformulation variable(s) stored in the optimizer
model. Use [`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
if a result exists before asking for values. For extensions, this only works if
[`optimizer_model_variable`](@ref) has been extended correctly and/or
[`map_value`](@ref) has been extended for variables.

**Example**
```julia-repl
julia> value(z)
42.0
```
"""
function JuMP.value(vref::InfOptVariableRef)
    return map_value(vref, Val(optimizer_model_key(JuMP.owner_model(vref))))
end

"""
    JuMP.value(cref::GeneralConstraintRef)

Extend [`JuMP.value`](@ref JuMP.value(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the value(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model. Use [`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel))
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
function JuMP.value(cref::GeneralConstraintRef)
    return map_value(cref, Val(optimizer_model_key(JuMP.owner_model(cref))))
end

"""
    map_optimizer_index(ref, key)

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
    map_dual(cref::GeneralConstraintRef, key)

Map the dual(s) of `cref` to its counterpart in the optimizer
model type that is distininguished by its extension key `key` as type `Val{ext_key_name}`.
Here `ref` need refer to methods for both variable references and constraint
references. This only needs to be defined for reformulation extensions that cannot
readily extend `optimizer_model_variable` and `optimizer_model_constraint`.
Such as is the case with reformuations that do not have a direct mapping between
variables and/or constraints in the original infinite form. Otherwise,
`optimizer_model_variable` and `optimizer_model_constraint` are used to make
these mappings by default.
"""
function map_dual end

# Default method that depends on optimizer_model_constraint --> making extensions easier
function map_dual(cref::GeneralConstraintRef, key)
    opt_cref = optimizer_model_constraint(cref)
    if opt_cref isa AbstractArray
        return JuMP.dual.(opt_cref)
    else
        return JuMP.dual(opt_cref)
    end
end

"""
    JuMP.dual(cref::GeneralConstraintRef)

Extend [`JuMP.dual`](@ref JuMP.dual(::JuMP.ConstraintRef{JuMP.Model, <:JuMP._MOICON}))
to return the dual(s) of `cref` in accordance with its reformulation constraint(s)
stored in the optimizer model. Use [`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel))
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
function JuMP.dual(cref::GeneralConstraintRef)
    return map_dual(cref, Val(optimizer_model_key(JuMP.owner_model(cref))))
end

# Error redriect dor variable call
function JuMP.dual(vref::GeneralVariableRef)
    return JuMP.dual(JuMP.VariableRef(JuMP.Model(), MOI.VariableIndex(1)))
end

"""
    map_shadow_price(cref::GeneralConstraintRef, key)

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
works if [`optimizer_model_constraint`](@ref) has been extended correctly.

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

# TODO add new functions
