"""
    InfiniteVariableRef <: JuMP.AbstractVariableRef
A DataType for infinite dimensional decision variables.
**Fields**
- `model::Model` Infinite model.
- `index` Index of variable in model.
"""
struct InfiniteVariableRef <: JuMP.AbstractVariableRef
    model::Model
    index
end

"""
    InfiniteAffExpr <: JuMP.GenericAffExpr
A `GenericAffExpr` that contains infinite dimensional variables.
"""
const InfiniteAffExpr = JuMP.GenericAffExpr{JuMP.AffExpr, InfiniteVariableRef}
InfiniteAffExpr() = InfiniteAffExpr(InfiniteVariableRef[], JuMP.AffExpr[], JuMP.AffExpr())

"""
    Data
A DataType for storing the data necessary to manage the bookkeeping of the infinite variables `InfiniteVariableRef`,
the infinite constraints `InfiniteConstraints`, and the transcription data.
**Fields**
- `inf_constrs::Vector{JuMP.AbstractConstraint}` Constraints that involve infinite variables.
- `num_inf_vars::Int` The number of `InfiniteVariableRef` that have been added to the model.
- `var_names::Vector{AbstractString}` The symbolic name of each `InfiniteVariableRef`.
"""
mutable struct Data
    inf_constrs::Vector{JuMP.AbstractConstraint}

    # Recourse variable data
    num_inf_vars::Int
    var_names::Vector{AbstractString}

    # Default constructor
    function Data()
        new(InfiniteConstraint[],
            0,
            String[])
    end
end
