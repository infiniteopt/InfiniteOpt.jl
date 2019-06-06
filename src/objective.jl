"""
    JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense, f::JuMP.AbstractJuMPScalar)
Extend the `JuMP.set_objective` function to accomodate `InfiniteModel` objects.
"""
function JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense,
                            f::JuMP.AbstractJuMPScalar)
    if f isa InfiniteExpr
        error("Objective function cannot contain infinite variables.")
    end
    m.objective_sense = sense
    m.objective_function = f
    return
end

"""
    JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense, f::Real)
Extend the `JuMP.set_objective` function to accomodate `InfiniteModel` objects.
"""
function JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense, f::Real)
    m.objective_sense = sense
    m.objective_function = JuMP.GenericAffExpr{Float64, GlobalVariableRef}(f)
    return
end

"""
    JuMP.objective_sense(model::InfiniteModel)
Extend the `JuMP.objective_sense` function to accomodate `InfiniteModel` objects.
"""
JuMP.objective_sense(model::InfiniteModel) = model.objective_sense

"""
    JuMP.set_objective_sense(m::InfiniteModel, sense)
Extend the `JuMP.set_objective_sense` function to accomodate `InfiniteModel` objects.
"""
function JuMP.set_objective_sense(model::InfiniteModel, sense)
    model.objective_sense = sense
    return
end

"""
    JuMP.objective_function_type(model::InfiniteModel)
Extend the `JuMP.objective_function_type` function to accomodate `InfiniteModel` objects.
"""
JuMP.objective_function_type(model::InfiniteModel) = typeof(model.objective_function)

"""
    JuMP.objective_function(model::InfiniteModel)
Extend the `JuMP.objective_function` function to accomodate `InfiniteModel` objects.
"""
JuMP.objective_function(model::InfiniteModel) = model.objective_function

"""
    JuMP.objective_function(model::InfiniteModel, FT::Type)
Extend the `JuMP.objective_function` function to accomodate `InfiniteModel` objects.
"""
function JuMP.objective_function(model::InfiniteModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(model.objective_function isa FT)
        throw(InexactError(:objective_function, FT,
                           typeof(model.objective_function)))
    end
    return model.objective_function::FT
end
