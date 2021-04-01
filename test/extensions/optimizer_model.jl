## A template for defining a new optimizer model

# Define a mutable struct for storing infinite model to optimizer model mappings
# plus other needed information
mutable struct NewReformData
    # Variable mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_mappings::Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}
    finvar_mappings::Dict{GeneralVariableRef, JuMP.VariableRef}

    # Map other variable info (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_to_supports::Dict{GeneralVariableRef, Vector{<:Tuple}}

    # Measure mappings and other info (ONLY NEEDED TO ENABLE MEASURE QUERIES)
    meas_mappings::Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}
    meas_to_supports::Dict{GeneralVariableRef, Vector{<:Tuple}}

    # Constraint mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    constr_mappings::Dict{InfOptConstraintRef, Vector{JuMP.ConstraintRef}}

    # Map other constraint info (REPLACE AND REFORMAT BELOW AS NEEDED)
    constr_to_supports::Dict{InfOptConstraintRef, Vector{<:Tuple}}

    # ADD OTHER MAPPING/MODEL INFORMATION HERE

    # default constructor
    function NewReformData() # REFORMAT BELOW IN ACCORDANCE WITH ABOVE ATTRIBUTES
        return new(Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{GeneralVariableRef, JuMP.VariableRef}(),
                   Dict{GeneralVariableRef, Vector{Tuple}}(),
                   Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{GeneralVariableRef, Vector{Tuple}}(),
                   Dict{InfOptConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfOptConstraintRef, Vector{Tuple}}()
                   )
    end
end

# Define the optimizer model key 
const OptKey = :ReformData # REPLACE WITH A DESIRED UNIQUE KEY

# Make a constructor for new optimizer model type (extension of JuMP.Model)
function NewReformModel(args...; kwargs...)::JuMP.Model # ADD EXPLICT ARGS AS NEEDED
    # initialize the JuMP Model
    model = JuMP.Model(args...; kwargs...)

    # ADD ADDITIONAL OPERATIONS IF NEEDED

    # add the extension data with a chosen optimizer model key and return
    model.ext[OptKey] = NewReformData()
    return model
end

# Make function for extracting the data from the model (optional)
function reform_data(model::JuMP.Model)::NewReformData
    # UPDATE THE NOMENCLATURE AS NEEDED
    haskey(model.ext, OptKey) || error("Model is not a NewReformModel.")
    return model.ext[OptKey]
end

# Extend build_optimizer_model! (enables build_optimizer_model! and optimize!)
function InfiniteOpt.build_optimizer_model!(model::InfiniteModel,
                                            key::Val{OptKey};
                                            my_kwarg::Bool = true) # ADD KEYWORD ARGUMENTS AS NEEDED
    # clear the model for a build/rebuild
    reform_model = clear_optimizer_model_build!(model)

    # IT MAY BE USEFUL TO CALL `expand_all_measures!` TO HANDLE MEASURES FIRST
    # otherwise can extend `add_measure_variable` and `delete_semi_infinite_variable` to
    # expand in place without modifying the infinite model

    # REPLACE LINES 73-89 WITH OPERATIONS TO BUILD A `NewReformModel` BASED ON `model`
    # these lines just generate artificial data, see `TranscriptionOpt` for a thorough example of implementation
    data = reform_data(reform_model)
    for vref in all_variables(model)
        if index(vref) isa Union{InfiniteVariableIndex, SemiInfiniteVariableIndex}
            data.infvar_mappings[vref] = [@variable(reform_model) for i = 1:2]
            data.infvar_to_supports[vref] = [(0.,), (1.,)]
        else
            data.finvar_mappings[vref] = @variable(reform_model)
        end
    end
    for mref in all_measures(model)
        data.meas_mappings[mref] = [@variable(reform_model) for i = 1:2]
        data.meas_to_supports[mref] = [(-1.,), (-2.,)]
    end
    for cref in all_constraints(model)
        data.constr_mappings[cref] = [@constraint(reform_model, @variable(reform_model) >= 0) for i = 1:2]
        data.constr_to_supports[cref] = [(2.,), (3.,)]
    end
    set_objective(reform_model, objective_sense(model), @variable(reform_model))
    # update the optimizer model status
    set_optimizer_model_ready(model, true)
    return
end

# Extend optimizer_model_variable if appropriate to enable variable related queries
function InfiniteOpt.optimizer_model_variable(vref::GeneralVariableRef,
                                              key::Val{OptKey};
                                              my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL VARIABLE(S)
    model = optimizer_model(JuMP.owner_model(vref))
    vindex = index(vref)
    if vindex isa Union{InfiniteVariableIndex, SemiInfiniteVariableIndex}
        map_dict = reform_data(model).infvar_mappings
    elseif vindex isa MeasureIndex
        map_dict = reform_data(model).meas_mappings
    else
        map_dict = reform_data(model).finvar_mappings
    end
    haskey(map_dict, vref) || error("Variable $vref not used in the optimizer model.")
    return map_dict[vref]
end

# Extend optimizer_model_expression if appropriate to enable expression related queries
function InfiniteOpt.optimizer_model_expression(expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
                                                key::Val{OptKey};
                                                my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO REFORMULATED EXPRESSIONS
    reform_expr = zero(AffExpr)
    return reform_expr
end

# Extend optimizer_model_constraint if appropriate to enable constraint related queries
function InfiniteOpt.optimizer_model_constraint(cref::InfOptConstraintRef,
                                                key::Val{OptKey};
                                                my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL CONSTRAINT(S)
    model = optimizer_model(JuMP.owner_model(cref))
    map_dict = reform_data(model).constr_mappings
    haskey(map_dict, cref) || error("Constraint $cref not used in the optimizer model.")
    return map_dict[cref]
end

# If appropriate extend variable_supports (enables support queries of infinite variables)
function InfiniteOpt.variable_supports(model::JuMP.Model,
                                       vref::Union{InfiniteVariableRef, SemiInfiniteVariableRef},
                                       key::Val{OptKey};
                                       my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = reform_data(model).infvar_to_supports
    gvref = InfiniteOpt._make_variable_ref(owner_model(vref), index(vref))
    haskey(map_dict, gvref) || error("Variable $gvref not used in the optimizer model.")
    return map_dict[gvref]
end

# If appropriate extend variable_supports for measures (enables support queries of measures)
function InfiniteOpt.variable_supports(model::JuMP.Model,
                                       mref::MeasureRef,
                                       key::Val{OptKey};
                                       my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = reform_data(model).meas_to_supports
    gvref = InfiniteOpt._make_variable_ref(owner_model(mref), index(mref))
    haskey(map_dict, gvref) || error("Variable $gvref not used in the optimizer model.")
    return map_dict[gvref]
end

# If appropriate extend expression_supports (enables support queries of expressions)
function InfiniteOpt.expression_supports(model::JuMP.Model,
                                         expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
                                         key::Val{OptKey};
                                         my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO SUPPORT(S) OF THE EXPRESSION(S)
    supps = [(-42.,), (1.,)]
    return supps
end

# If appropriate extend constraint_supports (enables support queries of constraints)
function InfiniteOpt.constraint_supports(model::JuMP.Model,
                                         cref::InfOptConstraintRef,
                                         key::Val{OptKey};
                                         my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = reform_data(model).constr_to_supports
    haskey(map_dict, cref) || error("Constraint $cref does not have supports in the optimizer model.")
    return length(map_dict[cref]) == 1 ? first(map_dict[cref]) : map_dict[cref]
end

# If it doesn't make sense to extend optimizer_model_variable, 
# optimizer_model_expression and/or optimizer_model_constraint then you'll need 
# to extend the following:
# - InfiniteOpt.map_value to enable JuMP.value
# - InfiniteOpt.map_optimizer_index to enable JuMP.optimizer_index
# - InfiniteOpt.map_dual to enable JuMP.dual
# - InfiniteOpt.map_shadow_price to enable JuMP.shadow_price
# - InfiniteOpt.map_lp_rhs_perturbation_range to enable JuMP.lp_rhs_perturbation_range
# - InfiniteOpt.map_lp_objective_perturbation_range to enable JuMP.lp_objective_perturbation_range
