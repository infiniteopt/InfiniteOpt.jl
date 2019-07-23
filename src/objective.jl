"""
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                       f::JuMP.AbstractJuMPScalar)

Extend [`JuMP.set_objective`](@ref) to set the objective of infinite model
`model`. Errors if `f` contains infinite variables and/or parameters.

**Example**
```julia
julia> set_objective(model, MOI.MIN_SENSE, x + measure(g + 2, tdata))

julia> objective_function(model)
x + measure(g(t) + 2)
```
"""
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                            f::JuMP.AbstractJuMPScalar)
    if isa(f, InfiniteExpr) || isa(f, ParameterExpr)
        error("Objective function cannot contain infinite parameters/variables.")
    end
    model.objective_sense = sense
    model.objective_function = f
    for vindex in keys(model.var_in_objective)
        model.var_in_objective[vindex] = false
    end
    for mindex in keys(model.meas_in_objective)
        model.meas_in_objective[mindex] = false
    end
    vrefs = _all_function_variables(f)
    for vref in vrefs
        if isa(vref, InfOptVariableRef)
            model.var_in_objective[JuMP.index(vref)] = true
        elseif isa(vref, MeasureRef)
            model.meas_in_objective[JuMP.index(vref)] = true
        end
    end
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                       f::Real)

Extend [`JuMP.set_objective`](@ref) to set the objective of infinite model
`model` with a number.

**Example**
```julia
julia> set_objective(model, MOI.MIN_SENSE, 3)

julia> objective_function(model)
3
```
"""
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                            f::Real)
    model.objective_sense = sense
    model.objective_function = JuMP.GenericAffExpr{Float64, GlobalVariableRef}(f)
    for vindex in keys(model.var_in_objective)
        model.var_in_objective[vindex] = false
    end
    for mindex in keys(model.meas_in_objective)
        model.meas_in_objective[mindex] = false
    end
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.objective_sense(model::InfiniteModel)::MOI.OptimizationSense

Extend [`JuMP.objective_sense`](@ref) to return the objective sense of the
infinite model `model`.

**Example**
```julia
julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```
"""
function JuMP.objective_sense(model::InfiniteModel)::MOI.OptimizationSense
    return model.objective_sense
end

"""
    JuMP.set_objective_sense(model::InfiniteModel, sense::MOI.OptimizationSense)

Extend [`JuMP.set_objective_sense`](@ref) to set the objective sense of infinite
model `model`.

**Example**
```julia
julia> set_objective_sense(model, MOI.MIN_SENSE)

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```
"""
function JuMP.set_objective_sense(model::InfiniteModel,
                                  sense::MOI.OptimizationSense)
    model.objective_sense = sense
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.objective_function_type(model::InfiniteModel)::Type{<:JuMP.AbstractJuMPScalar}

Extend [`JuMP.objective_function_type`](@ref) to return the objective function
type of infinite model `model`.

**Example**
```julia
julia> objective_function_type(model)
GenericAffExpr{Float64, FiniteVariableRef}
```
"""
function JuMP.objective_function_type(model::InfiniteModel)::Type{<:JuMP.AbstractJuMPScalar}
    return typeof(model.objective_function)
end

"""
    JuMP.objective_function(model::InfiniteModel)::JuMP.AbstractJuMPScalar

Extend [`JuMP.objective_function`](@ref) to return the objective of infinite
model `model`.

**Example**
```julia
julia> objective_function(model)
x + measure(g(t))
```
"""
function JuMP.objective_function(model::InfiniteModel)::JuMP.AbstractJuMPScalar
    return model.objective_function
end

# Extend objective_function for infinite models
function JuMP.objective_function(model::InfiniteModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(model.objective_function isa FT)
        throw(InexactError(:objective_function, FT,
                           typeof(model.objective_function)))
    end
    return model.objective_function::FT
end
