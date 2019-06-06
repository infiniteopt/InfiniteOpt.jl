# """
# InfOptVariable{S, T, U, V} <: JuMP.AbstractVariable
# A DataType for storing variable info
# **Fields**
# - `info::JuMP.VariableInfo{S, T, U, V}` Variable information.
# - `type::Symbol` Variable type (:Infinite, :Point, :Global).
# """
# struct InfOptVariable{S, T, U, V} <: JuMP.AbstractVariable
#     info::JuMP.VariableInfo{S, T, U, V}
#     type::Symbol
# end

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
    # Measure Data
    next_meas_index::Int
    measures::Dict{Int, Measure}
    meas_to_name::Dict{Int, String}

    # Variable data
    next_var_index::Int                             # Next variable index is nextvaridx+1
    vars::Dict{Int, InfOptVariable}                 # Map varidx -> variable
    var_to_name::Dict{Int, String}                  # Map varidx -> name
    name_to_var::Union{Dict{String, Int}, Nothing}  # Map varidx -> name

    # Constraint Data
    next_constr_index::Int                            # Next constraint index is nextconidx+1
    constrs::Dict{Int, JuMP.AbstractConstraint}       # Map conidx -> variable
    constr_to_name::Dict{Int, String}                 # Map conidx -> name
    name_to_constr::Union{Dict{String, Int}, Nothing} # Map name -> conidx

    # Objective Data
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar

    # Objects
    obj_dict::Dict{Symbol, Any} # Same that JuMP.Model's field `obj_dict`

    # Default constructor
    function InfiniteModel()
        new(0, Dict{Int, Measure}(), # Measures
            0, Dict{Int, JuMP.AbstractVariable}(),  # Variables
            Dict{Int, String}(), nothing,
            0, Dict{Int, JuMP.AbstractConstraint}(),
            Dict{Int, String}(), nothing,            # Constraints
            MOI.FEASIBILITY_SENSE,
            zero(JuMP.GenericAffExpr{Float64, FiniteVariableRef}),
            Dict{Symbol, Any}())
    end
end

Base.broadcastable(model::InfiniteModel) = Ref(model)

JuMP.object_dictionary(model::InfiniteModel) = model.obj_dict
