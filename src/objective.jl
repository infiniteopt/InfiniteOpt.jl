

"""
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense, f::JuMP.AbstractJuMPScalar)
Extend the `JuMP.set_objective` function to accomodate `InfiniteModel` objects.
"""
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense,
                            f::JuMP.AbstractJuMPScalar)
    if f isa InfiniteExpr || f isa ParameterExpr
        error("Objective function cannot contain infinite parameters/variables.")
    end
    model.objective_sense = sense
    model.objective_function = f
    for vindex in keys(model.var_in_objective)
        model.var_in_objective[vindex] = false
    end
    vrefs = _all_function_variables(f)
    for vref in vrefs
        if isa(vref, InfOptVariableRef)
            model.var_in_objective[JuMP.index(vref)] = true
        end
    end
    return
end

"""
    JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense, f::Real)
Extend the `JuMP.set_objective` function to accomodate `InfiniteModel` objects.
"""
function JuMP.set_objective(model::InfiniteModel, sense::MOI.OptimizationSense, f::Real)
    model.objective_sense = sense
    model.objective_function = JuMP.GenericAffExpr{Float64, GlobalVariableRef}(f)
    for vindex in keys(model.var_in_objective)
        model.var_in_objective[vindex] = false
    end
    return
end

"""
    JuMP.objective_sense(model::InfiniteModel)
Extend the `JuMP.objective_sense` function to accomodate `InfiniteModel` objects.
"""
JuMP.objective_sense(model::InfiniteModel) = model.objective_sense

"""
    JuMP.set_objective_sense(model::InfiniteModel, sense)
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
