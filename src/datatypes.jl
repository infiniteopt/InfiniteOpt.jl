"""
    InfOptVariable{S, T, U, V} <: JuMP.AbstractVariable
A DataType for storing variable info
**Fields**
- `info::JuMP.VariableInfo{S, T, U, V}` Variable information.
- `type::Symbol` Variable type (:Infinite, :Point, :Global).
"""
struct InfOptVariable{S, T, U, V} <: JuMP.AbstractVariable
    info::JuMP.VariableInfo{S, T, U, V}
    type::Symbol
end

"""
    AbstractMeasureData
An abstract type to define data for measures as described by the `Measure` DataType.
"""
abstract type AbstractMeasureData end

"""
    MeasureData{T <: AbstractVector, V <: AbstractArray} <: AbstractMeasureData
A DataType for measure abstraction data where the measure abstraction is of the
form: ``measure = \\int_{\\tau \\in T} f(\\tau) w(\\tau) d\\tau \\approx \\frac{1}{N} \\sum_{i = 1}^N \\alpha_i f(\\tau_i) w(\\tau_i)``.
**Fields**
- `weight_func::Function` Weighting function ``w`` must map input of type `V` to a scalar value.
- `coeffs::T` Coefficients ``\\alpha_i`` for the above measure abstraction.
- `eval_points::V` Evaluation points ``\\tau_i`` for the above measure abstraction.
"""
struct DiscretizeData{T <: AbstractVector, V <: AbstractArray} <: AbstractMeasureData
    weight_func::Function
    coeffs::T
    eval_points::V
end

"""
    Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}
A DataType for measure abstractions.
**Fields**
- `name::String` Measure name
- `func::T` Infinite variable expression.
- `data::V` Data of the abstraction as described in a `AbstractMeasureData` subtype.
"""
struct Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}
    name::String
    func::T
    data::V
end

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
        new(0, Dict{Int, Measure}(), Dict{Int, String}(), # Measures
            0, Dict{Int, JuMP.AbstractVariable}(),  # Variables
            Dict{Int, String}(), nothing,
            0, Dict{Int, JuMP.AbstractConstraint}(),
            Dict{Int, String}(), nothing,            # Constraints
            MOI.FEASIBILITY_SENSE,
            zero(JuMP.GenericAffExpr{Float64, FiniteVariableRef}),
            Dict{Symbol, Any}())
    end
end

# Define basic InfiniteModel extensions
Base.broadcastable(model::InfiniteModel) = Ref(model)
JuMP.object_dictionary(model::InfiniteModel) = model.obj_dict

"""
    GeneralVariableRef <: JuMP.AbstractVariableRef
An abstract type to define new variable types.
"""
abstract type GeneralVariableRef <: JuMP.AbstractVariableRef end

"""
    MeasureFiniteVariableRef <: GeneralVariableRef
An abstract type to define finite variable and measure types.
"""
abstract type MeasureFiniteVariableRef <: GeneralVariableRef end

"""
    FiniteVariableRef <: GeneralVariableRef
An abstract type to define new finite variable types.
"""
abstract type FiniteVariableRef <: MeasureFiniteVariableRef end

"""
    GlobalVariableRef <: JuMP.FiniteVariableRef
A DataType for finite fixed variables (e.g., first stage variables,
steady-state variables).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct GlobalVariableRef <: FiniteVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

"""
    PointVariableRef <: JuMP.FiniteVariableRef
A DataType for variables defined at a transcipted point (e.g., second stage
variable at a particular scenario, dynamic variable at a discretized time point).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct PointVariableRef <: FiniteVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

"""
    InfiniteVariableRef <: JuMP.GeneralVariableRef
A DataType for untranscripted infinite dimensional variables (e.g., second stage
variables, dynamic variables).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct InfiniteVariableRef <: GeneralVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

# Define infinite expressions
const InfiniteAffExpr = JuMP.GenericAffExpr{Float64, Union{InfiniteVariableRef, GeneralVariableRef}}
const InfiniteQuadExpr = JuMP.GenericQuadExpr{Float64, Union{InfiniteVariableRef, GeneralVariableRef}}
const InfiniteExpr = Union{InfiniteAffExpr, InfiniteQuadExpr, InfiniteVariableRef}

"""
    MeasureRef <: JuMP.FiniteVariableRef
A DataType for referring to measure abstractions.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct MeasureRef <: MeasureFiniteVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

# Define finite measure expressions (note infinite expression take precedence)
const MeasureAffExpr = JuMP.GenericAffExpr{Float64, Union{MeasureFiniteVariableRef, MeasureRef}}
const MeasureQuadExpr = JuMP.GenericQuadExpr{Float64, Union{MeasureFiniteVariableRef, MeasureRef}}
const MeasureExpr = Union{MeasureAffExpr, MeasureQuadExpr, MeasureRef}

"""
    GeneralConstraintRef
An abstract type to define new variable types.
"""
abstract type GeneralConstraintRef end

"""
InfiniteConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain infinite variables.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct InfiniteConstraintRef <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::JuMP.AbstractShape
end

"""
    FiniteConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain finite variables.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct FiniteConstraintRef <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::JuMP.AbstractShape
end

"""
    MeasureConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain finite variables and measures.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct MeasureConstraintRef <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::JuMP.AbstractShape
end
