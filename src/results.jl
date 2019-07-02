"""
    JuMP.termination_status(model::InfiniteModel)
Return the reason why the solver stopped (i.e., the MathOptInterface model
attribute `TerminationStatus`).
"""
function JuMP.termination_status(model::InfiniteModel)
    return JuMP.termination_status(optimizer_model(model))
end

"""
     JuMP.primal_status(model::InfiniteModel)
Return the status of the most recent primal solution of the solver (i.e., the
MathOptInterface model attribute `PrimalStatus`).
"""
function  JuMP.primal_status(model::InfiniteModel)
    return JuMP.primal_status(optimizer_model(model))
end

"""
    JuMP.dual_status(model::InfiniteModel)
Return the status of the most recent dual solution of the solver (i.e., the
MathOptInterface model attribute `DualStatus`).
"""
function JuMP.dual_status(model::InfiniteModel)
    return JuMP.dual_status(optimizer_model(model))
end

"""
    JuMP.solve_time(model::InfiniteModel)
If available, returns the solve time reported by the solver.
Returns "ArgumentError: ModelLike of type `Solver.Optimizer` does not support accessing
the attribute MathOptInterface.SolveTime()" if the attribute is
not implemented.
"""
# TODO include with next JuMP version
function solve_time(model::InfiniteModel)
    return MOI.get(optimizer_model(model), MOI.SolveTime())
end

"""
    has_values(model::InfiniteModel)
Return `true` if the solver has a primal solution available to query, otherwise
return `false`. See also [`value`](@ref).
"""
JuMP.has_values(model::InfiniteModel) = JuMP.primal_status(model) != MOI.NO_SOLUTION

"""
    JuMP.objective_bound(model::InfiniteModel)::Float64
Return the best known bound on the optimal objective value after a call to
`optimize!(model)`.
"""
function JuMP.objective_bound(model::InfiniteModel)::Float64
    return JuMP.objective_bound(optimizer_model(model))
end

"""
    JuMP.objective_value(model::InfiniteModel)
Return the objective value after a call to `optimize!(model)`.
"""
function JuMP.objective_value(model::InfiniteModel)::Float64
    return JuMP.objective_value(optimizer_model(model))
end

"""
    map_value(vref::GeneralVariableRef, key)
Map the value of `vref` to its counterpart in the optimizer model type distininguished
by its extension key `key` as type `Val{ext_key_name}`.
"""
function map_value end

"""
    JuMP.value(vref::GeneralVariableRef)
Get the value of this variable in the result returned by a solver. Use
[`JuMP.has_values`](@ref) to check if a result exists before asking for values.
"""
function JuMP.value(vref::GeneralVariableRef)
    return map_value(vref::GeneralVariableRef,
                     Val(optimizer_model_key(JuMP.owner_model(vref))))
end
