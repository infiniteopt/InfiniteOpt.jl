function JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense,
                            f::JuMP.AbstractJuMPScalar)
    m.objective_sense = sense
    m.objective_function = f
end
function JuMP.set_objective(m::InfiniteModel, sense::MOI.OptimizationSense, f::Real)
    m.objective_sense = sense
    m.objective_function = JuMP.GenericAffExpr{Float64, MyVariableRef}(f)
end
JuMP.objective_sense(model::InfiniteModel) = model.objective_sense
function JuMP.set_objective_sense(model::InfiniteModel, sense)
    model.objective_sense = sense
end
JuMP.objective_function_type(model::InfiniteModel) = typeof(model.objective_function)
JuMP.objective_function(model::InfiniteModel) = model.objective_function
function JuMP.objective_function(model::InfiniteModel, FT::Type)
    # InexactError should be thrown, this is needed in `objective.jl`
    if !(model.objective_function isa FT)
        throw(InexactError(:objective_function, FT,
                           typeof(model.objective_function)))
    end
    return model.objective_function::FT
end
