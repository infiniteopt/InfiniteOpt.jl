"""
    JuMP.set_objective_function(model::InfiniteModel,
                                func::JuMP.AbstractJuMPScalar)

Extend [`JuMP.set_objective_function`](@ref) to set the objective expression of
infinite model `model`. Errors if `func` contains infinite variables and/or
parameters. Also errors if `func` contains invalid variables.

**Example**
```julia
julia> set_objective_function(model, x + measure(g + 2, tdata))

julia> objective_function(model)
x + measure(g(t) + 2)
```
"""
function JuMP.set_objective_function(model::InfiniteModel,
                                     func::JuMP.AbstractJuMPScalar)
    # check the function
    if isa(func, InfiniteExpr) || isa(func, ParameterExpr)
        error("Objective function cannot contain infinite parameters/variables.")
    end
    JuMP.check_belongs_to_model(func, model)
    # update the function
    model.objective_function = func
    # delete old mappings
    for vindex in keys(model.var_in_objective)
        model.var_in_objective[vindex] = false
    end
    for mindex in keys(model.meas_in_objective)
        model.meas_in_objective[mindex] = false
    end
    # update new mappings
    vrefs = _all_function_variables(func)
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
    JuMP.set_objective_function(model::InfiniteModel, func::Real)

Extend [`JuMP.set_objective_function`](@ref) to set the objective expression of
`model` with a number.

**Example**
```julia
julia> set_objective_function(model, 3)

julia> objective_function(model)
3
```
"""
function JuMP.set_objective_function(model::InfiniteModel, func::Real)
    # update function
    model.objective_function = JuMP.GenericAffExpr{Float64, GlobalVariableRef}(func)
    # delete old mappings
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
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                       func::Union{JuMP.AbstractJuMPScalar, Real})

Extend [`JuMP.set_objective`](@ref) to set the objective of infinite model
`model`. Errors if `func` contains infinite variables and/or parameters, or if
it does not belong to the model.

**Example**
```julia
julia> set_objective(model, MOI.MIN_SENSE, x + measure(g + 2, tdata))

julia> objective_function(model)
x + measure(g(t) + 2)
```
"""
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                            func::Union{JuMP.AbstractJuMPScalar, Real})
    JuMP.set_objective_sense(model, sense)
    JuMP.set_objective_function(model, func)
    return
end

# Fallback
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                            func)
    error("The objective function `$(func)` is not supported.")
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

"""
    JuMP.set_objective_coefficient(model::InfiniteModel,
                                   variable::GeneralVariableRef,
                                   coefficient::Real)

Set the linear objective coefficient associated with `variable` to `coefficient`.
Errors if the function type is unsupported.

**Example**
```julia
julia> objective_function(model)
x + y

julia> set_objective_coefficient(model, y, 2)

julia> objective_function(model)
x + 2 y
```
"""
function JuMP.set_objective_coefficient(model::InfiniteModel,
                                        variable::GeneralVariableRef,
                                        coeff::Real)
    new_expr = _set_variable_coefficient!(JuMP.objective_function(model),
                                          variable, coeff)
    JuMP.set_objective_function(model, new_expr)
    return
end
