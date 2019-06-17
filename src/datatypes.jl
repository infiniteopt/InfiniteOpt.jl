"""
    AbstractInfiniteSet
An abstract type to define the sets describing infinite variable parameters.
"""
abstract type AbstractInfiniteSet end

"""
    InfOptParameter{T <: AbstractInfiniteSet} <: JuMP.AbstractVariable
A DataType for storing infinite parameter info
**Fields**
- `set::T` The set that contains the parameter.
"""
struct InfOptParameter{T <: AbstractInfiniteSet} <: JuMP.AbstractVariable
    set::T
end

"""
    AbstractMeasureData
An abstract type to define infinite, point, and global variables.
"""
abstract type InfOptVariable <: JuMP.AbstractVariable end

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
- 'name::String' Name of the measure that will be implemented.
- `weight_func::Function` Weighting function ``w`` must map input of type `V` to a scalar value.
- `coeffs::T` Coefficients ``\\alpha_i`` for the above measure abstraction.
- `eval_points::V` Evaluation points ``\\tau_i`` for the above measure abstraction.
"""
struct DiscretizeData{T <: AbstractVector, V <: AbstractArray} <: AbstractMeasureData
    name::String
    weight_func::Function
    coeffs::T
    eval_points::V
end

"""
    Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}
A DataType for measure abstractions.
**Fields**
- `func::T` Infinite variable expression.
- `data::V` Data of the abstraction as described in a `AbstractMeasureData` subtype.
"""
struct Measure{T <: JuMP.AbstractJuMPScalar, V <: AbstractMeasureData}
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

    # Parameter Data
    next_param_index::Int
    params::Dict{Int, InfOptParameter}
    param_to_name::Dict{Int, String}
    name_to_param::Union{Dict{String, Int}, Nothing}

    # Variable data
    next_var_index::Int                             # Next variable index is nextvaridx+1
    vars::Dict{Int, InfOptVariable}                 # Map varidx -> variable
    var_to_name::Dict{Int, String}                  # Map varidx -> name
    name_to_var::Union{Dict{String, Int}, Nothing}  # Map varidx -> name
    var_to_lower_bound::Dict{Int, Int}
    var_to_upper_bound::Dict{Int, Int}
    var_to_fix::Dict{Int, Int}
    var_to_zero_one::Dict{Int, Int}
    var_to_integrality::Dict{Int, Int}

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
            0, Dict{Int, InfOptParameter}(), Dict{Int, String}(), nothing, # Parameters
            0, Dict{Int, JuMP.AbstractVariable}(),  # Variables
            Dict{Int, String}(), nothing,
            Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}(),
            Dict{Int, Int}(), Dict{Int, Int}(),
            0, Dict{Int, JuMP.AbstractConstraint}(), # Constraints
            Dict{Int, String}(), nothing,
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

"""
    ParameterRef <: JuMP.GeneralVariableRef
A DataType for untranscripted infinite dimensional parameters that parameterize
the infinite variables.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct ParameterRef <: GeneralVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

"""
    InfiniteVariable{S, T, U, V} <: InfOptVariable
A DataType for storing infinite variable information
**Fields**
- `info::JuMP.VariableInfo{S, T, U, V}` Variable information.
- `parameter_refs::Tuple` The infinite parameters(s) with associated the variable.
"""
struct InfiniteVariable{S, T, U, V} <: InfOptVariable
    info::JuMP.VariableInfo{S, T, U, V}
    parameter_refs::Tuple
end

"""
    PointVariable{S, T, U, V} <: InfOptVariable
A DataType for storing point variable information
**Fields**
- `info::JuMP.VariableInfo{S, T, U, V}` Variable information.
- `infinite_variable_ref::InfiniteVariableRef`The infinite variable associated with the point variable.
- `parameter_values::Tuple` The infinite parameter evaluate values defining the point.
"""
struct PointVariable{S, T, U, V} <: InfOptVariable
    info::JuMP.VariableInfo{S, T, U, V}
    infinite_variable_ref::InfiniteVariableRef
    parameter_values::Tuple
end

"""
    GlobalVariable{S, T, U, V} <: InfOptVariable
A DataType for storing global variable information
**Fields**
- `info::JuMP.VariableInfo{S, T, U, V}` Variable information.
"""
struct GlobalVariable{S, T, U, V} <: InfOptVariable
    info::JuMP.VariableInfo{S, T, U, V}
end

# Define variable references without that aren't measures
const InfOptVariableRef = Union{InfiniteVariableRef, PointVariableRef, GlobalVariableRef}

# Define infinite expressions
const InfiniteExpr = Union{InfiniteVariableRef,
                           JuMP.GenericAffExpr{Float64, InfiniteVariableRef},
                           JuMP.GenericAffExpr{Float64, GeneralVariableRef},
                           JuMP.GenericQuadExpr{Float64, InfiniteVariableRef},
                           JuMP.GenericQuadExpr{Float64, GeneralVariableRef}}

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
const MeasureExpr = Union{MeasureRef,
                          JuMP.GenericAffExpr{Float64, MeasureRef},
                          JuMP.GenericAffExpr{Float64, MeasureFiniteVariableRef},
                          JuMP.GenericQuadExpr{Float64, MeasureRef},
                          JuMP.GenericQuadExpr{Float64, MeasureFiniteVariableRef}}

"""
    IntervalSet{T <: Number} <: AbstractInfiniteSet
A DataType that stores the lower and upper bounds for infinite parameters that are
rectangular and are associated with infinite variables.
**Fields**
- `lower_bound::T` Lower bound of the infinite parameter.
- `upper_bound::T` Upper bound of the infinite parameter.
"""
struct IntervalSet <: AbstractInfiniteSet
    lower_bound::Float64
    upper_bound::Float64
end

"""
    DistributionSet{T <: Distributions.NonMatrixDistribution} <: AbstractInfiniteSet
A DataType that stores the distribution for infinite parameters that are
random and are associated with infinite variables.
**Fields**
- `distribution::T` Distribution of the random parameter.
"""
struct DistributionSet{T <: Distributions.NonMatrixDistribution} <: AbstractInfiniteSet
    distribution::T
end

"""
    DiscreteSet{T <: Vector} <: AbstractInfiniteSet
A DataType that stores a dicrete set of points for infinite parameters that are
associated with infinite variables.
**Fields**
- `points::T` The discrete points that make up the set.
"""
struct DiscreteSet{T <: Vector} <: AbstractInfiniteSet
    points::T
end

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
struct InfiniteConstraintRef{S <: JuMP.AbstractShape} <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::S
end

"""
    FiniteConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain finite variables.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct FiniteConstraintRef{S <: JuMP.AbstractShape} <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::S
end

"""
    MeasureConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain finite variables and measures.
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct MeasureConstraintRef{S <: JuMP.AbstractShape} <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::S
end
