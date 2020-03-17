## A template for defining a new optimizer model

# Define a mutable struct for storing infinite model to optimizer model mappings
# plus other needed information
mutable struct NewReformData
    # Variable mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_to_vars::Dict{InfiniteVariableRef, Vector{JuMP.VariableRef}}
    holdvar_to_var::Dict{HoldVariableRef, JuMP.VariableRef}
    ptvar_to_var::Dict{PointVariableRef, JuMP.VariableRef}

    # Map other variable info (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_to_supports::Dict{InfiniteVariableRef, Vector{<:Tuple}}

    # Constraint mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    infcon_to_constrs::Dict{InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}
    meascon_to_constrs::Dict{MeasureConstraintRef, Vector{JuMP.ConstraintRef}}
    fincon_to_constr::Dict{FiniteConstraintRef, JuMP.ConstraintRef}

    # Map other constraint info (REPLACE AND REFORMAT BELOW AS NEEDED)
    infcon_to_supports::Dict{InfiniteConstraintRef, Vector{<:Tuple}}
    meascon_to_supports::Dict{MeasureConstraintRef, Vector{<:Tuple}}
    infcon_to_params::Dict{InfiniteConstraintRef, Tuple}
    meascon_to_params::Dict{MeasureConstraintRef, Tuple}

    # ADD OTHER MAPPING/MODEL INFORMATION HERE
    # default constructor
    function NewReformData() # REFORMAT BELOW IN ACCORDANCE WITH ABOVE ATTRIBUTES
        return new(Dict{InfiniteVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{HoldVariableRef, JuMP.VariableRef}(),
                   Dict{PointVariableRef, JuMP.VariableRef}(),
                   Dict{InfiniteVariableRef, Vector{Tuple}}(),
                   Dict{InfiniteConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{MeasureConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{FiniteConstraintRef, JuMP.ConstraintRef}(),
                   Dict{InfiniteConstraintRef, Vector{Tuple}}(),
                   Dict{MeasureConstraintRef, Vector{Tuple}}(),
                   Dict{InfiniteConstraintRef, Tuple}(),
                   Dict{MeasureConstraintRef, Tuple}())
    end
end

# Make a constructor for new optimizer model type (extension of JuMP.Model)
function NewReformModel(args...; kwargs...)::JuMP.Model # ADD EXPLICT ARGS AS NEEDED
    # initialize the JuMP Model
    model = JuMP.Model(args...; kwargs...)
    # ADD ADDITIONAL OPERATIONS IF NEEDED
    # add the extension data with a chosen optimizer model key and return
    model.ext[:ReformData] = NewReformData() # REPLACE WITH ACTUAL KEY AND DATA
    return model
end

# Extend build_optimizer_model! (enables build_optimizer_model! and optimize!)
function InfiniteOpt.build_optimizer_model!(model::InfiniteModel,
                                            key::Val{:ReformData}; # INSERT NEW KEY
                                            my_kwarg::Bool = true) # ADD KEYWORD ARGUMENTS AS NEEDED
    # retrieve model mode
    mode = JuMP.backend(optimizer_model(model)).mode

    # IT MAY BE USEFUL TO CALL expand_all_measures! TO HANDLE MEASURES FIRST

    # REPLACE LINES 60-67 WITH OPERATIONS TO BUILD A NewReformModel BASED ON `model`
    reform_model = TranscriptionModel(model, caching_mode = mode)
    reform_data = NewReformData()
    for i in eachindex(fieldnames(TranscriptionData))
         value = copy(getfield(transcription_data(reform_model), fieldnames(TranscriptionData)[i]))
        setfield!(reform_data, fieldnames(NewReformData)[i], value)
    end
    delete!(reform_model.ext, :TransData)
    reform_model.ext[:ReformData] = reform_data

    # add the optimizer and its attributes
    add_infinite_model_optimizer(model, reform_model) # REPLACE reform_model WITH ACTUAL
    # add the the built optimizer model to the infinite model
    set_optimizer_model(model, reform_model) # REPLACE reform_model WITH ACTUAL
    set_optimizer_model_ready(model, true)
    return
end

# Make function for extracting the data from the model (optional)
function reform_data(model::JuMP.Model)::NewReformData
    # REPLACE THE KEY BELOW WITH ACTUAL OPTIMIZER MODEL KEY
    haskey(model.ext, :ReformData) || error("Model is not a NewReformModel.")
    return model.ext[:ReformData]
end

# Extend optimizer_model_variable if appropriate to enable variable related queries
function InfiniteOpt.optimizer_model_variable(vref::InfOptVariableRef,
                                              key::Val{:ReformData}; # REPLACE WITH KEY
                                              my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL VARIABLE(S)
    model = optimizer_model(JuMP.owner_model(vref))
    if vref isa InfiniteVariableRef
        map_dict = reform_data(model).infvar_to_vars
    elseif vref isa PointVariableRef
        map_dict = reform_data(model).ptvar_to_var
    else
        map_dict = reform_data(model).holdvar_to_var
    end
    haskey(map_dict, vref) || error("Variable $vref not used in the optimizer model.")
    return map_dict[vref]
end

# Extend optimizer_model_constraint if appropriate to enable constraint related queries
function InfiniteOpt.optimizer_model_constraint(cref::GeneralConstraintRef,
                                                key::Val{:ReformData}; # REPLACE WITH KEY
                                                my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL CONSTRAINT(S)
    model = optimizer_model(JuMP.owner_model(cref))
    if cref isa InfiniteConstraintRef
        map_dict = reform_data(model).infcon_to_constrs
    elseif cref isa MeasureConstraintRef
        map_dict = reform_data(model).meascon_to_constrs
    else
        map_dict = reform_data(model).fincon_to_constr
    end
    haskey(map_dict, cref) || error("Constraint $cref not used in the optimizer model.")
    return map_dict[cref]
end

# If appropriate extend variable_supports (enables support queries of infinite variables)
function InfiniteOpt.variable_supports(model::JuMP.Model,
                                       vref::InfiniteVariableRef,
                                       key::Val{:ReformData}; # REPLACE WITH KEY
                                       my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = reform_data(model).infvar_to_supports
    haskey(map_dict, vref) || error("Variable $vref not used in the optimizer model.")
    return map_dict[vref]
end

# If appropriate extend constraint_supports (enables support queries of infinite constraints)
function InfiniteOpt.constraint_supports(model::JuMP.Model,
                                         cref::GeneralConstraintRef,
                                         key::Val{:ReformData}; # REPLACE WITH KEY
                                         my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    if cref isa InfiniteConstraintRef
        map_dict = reform_data(model).infcon_to_supports
    else
        map_dict = reform_data(model).meascon_to_supports
    end
    haskey(map_dict, cref) || error("Constraint $cref does not have supports in the optimizer model.")
    return map_dict[cref]
end

# If appropriate extend constraint_parameter_refs (enables parameter reference queries of infinite constraints)
function InfiniteOpt.constraint_parameter_refs(model::JuMP.Model,
                                              cref::GeneralConstraintRef,
                                              key::Val{:ReformData}; # REPLACE WITH KEY
                                              my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED)
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    if cref isa InfiniteConstraintRef
        map_dict = reform_data(model).infcon_to_params
    else
        map_dict = reform_data(model).meascon_to_params
    end
    haskey(map_dict, cref) || error("Constraint $cref does not have parameter references in the optimizer model.")
    return map_dict[cref]
end

# If it doesn't make sense to extend optimizer_model_variable and/or optimizer_model_constraint
# then you'll need to extend the following:
# - InfiniteOpt.map_value to enable JuMP.value
# - InfiniteOpt.map_optimizer_index to enable JuMP.optimizer_index
# - InfiniteOpt.map_dual to enable JuMP.dual
# - InfiniteOpt.map_shadow_price to enable JuMP.shadow_price
# - InfiniteOpt.map_lp_rhs_perturbation_range to enable JuMP.lp_rhs_perturbation_range
# - InfiniteOpt.map_lp_objective_perturbation_range to enable JuMP.lp_objective_perturbation_range
