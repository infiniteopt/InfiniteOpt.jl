## A template for defining a new optimizer model

# Define a mutable struct for storing infinite model to optimizer model mappings
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

#
