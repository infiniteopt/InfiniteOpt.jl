"""
    InfiniteModel(args...; [kwargs...])
Return a infinite model object which extends a JuMP model object to contain
**Arguments**
-
```julia
julia> m = InfiniteModel(with_optimizer(Gurobi.Optimizer))

```
"""
mutable struct InfiniteModel <: JuMP.AbstractModel
    # Transcription data

    # Model data
    next_var_index::Int                                 # Next variable index is nextvaridx+1
    variables::Dict{Int, JuMP.ScalarVariable}       # Map varidx -> variable
    var_to_name::Dict{Int, String}                  # Map varidx -> name
    name_to_var::Union{Dict{String, Int}, Nothing}  # Map varidx -> name
    next_constr_index::Int                                 # Next constraint index is nextconidx+1
    constraints::Dict{Int, JuMP.AbstractConstraint}      # Map conidx -> variable
    constr_to_name::Dict{Int, String}      # Map conidx -> name
    name_to_constr::Union{Dict{String, Int}, Nothing} # Map name -> conidx
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    obj_dict::Dict{Symbol, Any}                     # Same that JuMP.Model's field `obj_dict`
    function InfiniteModel()
        new(0, Dict{Int, JuMP.AbstractVariable}(),
            Dict{Int, String}(), nothing,                        # Variables
            0, Dict{Int, JuMP.AbstractConstraint}(),
            Dict{Int, String}(), nothing,            # Constraints
            MOI.FEASIBILITY_SENSE,
            zero(JuMP.GenericAffExpr{Float64, GlobalVariableRef}),
            Dict{Symbol, Any}())
    end
end

Base.broadcastable(model::InfiniteModel) = Ref(model)

JuMP.object_dictionary(model::InfiniteModel) = model.obj_dict
