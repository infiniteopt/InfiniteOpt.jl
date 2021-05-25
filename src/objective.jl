################################################################################
#                               ACCESSOR METHODS
################################################################################
"""
    JuMP.objective_sense(model::InfiniteModel)::MOI.OptimizationSense

Extend `JuMP.objective_sense` to return the objective sense of the infinite model 
`model`.

**Example**
```julia-repl
julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```
"""
function JuMP.objective_sense(model::InfiniteModel)::MOI.OptimizationSense
    return model.objective_sense
end

"""
    JuMP.objective_function(model::InfiniteModel)::JuMP.AbstractJuMPScalar

Extend `JuMP.objective_function` to return the objective of infinite model 
`model`.

**Example**
```julia-repl
julia> objective_function(model)
1
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
    JuMP.objective_function_type(model::InfiniteModel)::Type{<:JuMP.AbstractJuMPScalar}

Extend `JuMP.objective_function_type` to return the objective function type of 
infinite model `model`.

**Example**
```julia-repl
julia> objective_function_type(model)
GenericAffExpr{Float64,GeneralVariableRef}
```
"""
function JuMP.objective_function_type(model::InfiniteModel
    )::Type{<:JuMP.AbstractJuMPScalar}
    return typeof(JuMP.objective_function(model))
end

"""
    objective_has_measures(model::InfiniteModel)::Bool

Return `Bool` whether the objective function contains any measures.
"""
function objective_has_measures(model::InfiniteModel)::Bool
    return model.objective_has_measures
end

################################################################################
#                             DEFINITION METHODS
################################################################################
"""
    JuMP.set_objective_function(model::InfiniteModel,
                                func::JuMP.AbstractJuMPScalar)::Nothing

Extend `JuMP.set_objective_function` to set the objective expression of
infinite model `model`. Errors if `func` contains infinite variables and/or
parameters. Also errors if `func` contains invalid variables.

**Example**
```julia-repl
julia> set_objective_function(model, 2x + 1)

julia> objective_function(model)
2 x + 1
```
"""
function JuMP.set_objective_function(
    model::InfiniteModel,
    func::JuMP.AbstractJuMPScalar
    )::Nothing
    # gather the unique list of variable references for testing and mapping
    new_vrefs = _all_function_variables(func)
    # test in the model
    for vref in new_vrefs
        JuMP.check_belongs_to_model(vref, model)
    end
    if !isempty(_object_numbers(new_vrefs))
        error("Objective function cannot contain infinite parameters/variables.")
    end
    # delete old mappings
    old_vrefs = _all_function_variables(JuMP.objective_function(model))
    for vref in old_vrefs
        _data_object(vref).in_objective = false
    end
    model.objective_has_measures = false
    # update the function
    model.objective_function = func
    # update new mappings
    for vref in new_vrefs
        _data_object(vref).in_objective = true
        if _index_type(vref) == MeasureIndex
            model.objective_has_measures = true
        end
    end
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.set_objective_function(model::InfiniteModel, func::Real)::Nothing

Extend `JuMP.set_objective_function` to set the objective expression of
`model` with a number.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> set_objective_function(model, 3)

julia> objective_function(model)
3
```
"""
function JuMP.set_objective_function(model::InfiniteModel, func::Real)::Nothing
    # delete old mappings
    old_vrefs = _all_function_variables(JuMP.objective_function(model))
    for vref in old_vrefs
        _data_object(vref).in_objective = false
    end
    # update function
    model.objective_function = JuMP.GenericAffExpr{Float64, GeneralVariableRef}(func)
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.set_objective_sense(model::InfiniteModel,
                             sense::MOI.OptimizationSense)::Nothing

Extend `JuMP.set_objective_sense` to set the objective sense of infinite model 
`model`.

**Example**
```julia-repl
julia> set_objective_sense(model, MOI.MIN_SENSE)

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```
"""
function JuMP.set_objective_sense(
    model::InfiniteModel,
    sense::MOI.OptimizationSense
    )::Nothing
    model.objective_sense = sense
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                       func::Union{JuMP.AbstractJuMPScalar, Real})::Nothing

Extend `JuMP.set_objective` to set the objective of infinite model
`model`. Errors if `func` contains infinite variables and/or parameters, or if
it does not belong to the model.

**Example**
```julia-repl
julia> set_objective(model, MOI.MIN_SENSE, 2x + 1)

julia> objective_function(model)
2 x + 1
```
"""
function JuMP.set_objective(
    model::InfiniteModel, 
    sense::MOI.OptimizationSense,
    func::Union{JuMP.AbstractJuMPScalar, Real}
    )::Nothing
    JuMP.set_objective_sense(model, sense)
    JuMP.set_objective_function(model, func)
    return
end

# Fallback
function JuMP.set_objective(
    model::InfiniteModel, 
    sense::MOI.OptimizationSense,
    func
    )
    error("The objective function `$(func)` is not supported.")
end

"""
    JuMP.set_objective_coefficient(model::InfiniteModel,
                                   variable::GeneralVariableRef,
                                   coefficient::Real)::Nothing

Extend `JuMP.set_objective_coefficient` Set the linear objective coefficient 
associated with `variable` to `coefficient`. Errors if the function type is 
unsupported.

**Example**
```julia-repl
julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @objective(model, x + y)
x + y

julia> set_objective_coefficient(model, y, 2)

julia> objective_function(model)
x + 2 y
```
"""
function JuMP.set_objective_coefficient(
    model::InfiniteModel,
    variable::GeneralVariableRef,
    coeff::Real
    )::Nothing
    new_expr = _set_variable_coefficient!(JuMP.objective_function(model),
                                          variable, coeff)
    JuMP.set_objective_function(model, new_expr)
    return
end
