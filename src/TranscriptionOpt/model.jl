"""
    TranscriptionData

A DataType for storing the data mapping an `InfiniteModel` that has been
transcribed to a regular JuMP `Model` that contains the transcribed variables.

**Fields**
- `infinite_to_vars::Dict{InfiniteOpt.InfiniteVariableRef,
   Vector{JuMP.VariableRef}}`: Infinite variables to their transcribed variables.
- `global_to_var::Dict{InfiniteOpt.GlobalVariableRef, JuMP.VariableRef}`: Global
  variables to model variables.
- `point_to_var::Dict{InfiniteOpt.PointVariableRef, JuMP.VariableRef}`: Point
  variables to model variables.
- `infvar_to_supports::Dict{InfiniteOpt.InfiniteVariableRef, Dict}`: Infinite
  variables to transcribed supports indexed by their numeric aliases.
- `infinite_to_constrs::Dict{InfiniteOpt.InfiniteConstraintRef,
  Vector{JuMP.ConstraintRef}}`: Infinite constraints to their transcribed
  constraints.
- `measure_to_constrs::Dict{InfiniteOpt.MeasureConstraintRef,
  Vector{JuMP.ConstraintRef}}`: Measure constraints to model constraints.
- `finite_to_constr::Dict{InfiniteOpt.FiniteConstraintRef, JuMP.ConstraintRef}`:
  Finite constraints to model constraints.
- `infconstr_to_supports::Dict{InfiniteOpt.InfiniteConstraintRef, Dict}`: Infinite
  constraints to the transcribed supports indxed by their numeric aliases.
- `measconstr_to_supports::Dict{InfiniteOpt.MeasureConstraintRef, Dict}`:
  Measure constraints to the transcribed supports indxed by their numeric aliases.
- `infconstr_to_params::Dict{InfiniteOpt.InfiniteConstraintRef, Tuple}`: Infinite
  constraints to the parameter tuples associated with each transcribed support.
- `measconstr_to_params::Dict{InfiniteOpt.MeasureConstraintRef, Tuple}`: Measure
  constraints to the parameter tuples associated with each transcribed support.
"""
mutable struct TranscriptionData
    # Variable mapping
    infinite_to_vars::Dict{InfiniteOpt.InfiniteVariableRef,
                           Vector{JuMP.VariableRef}}
    global_to_var::Dict{InfiniteOpt.GlobalVariableRef, JuMP.VariableRef}
    point_to_var::Dict{InfiniteOpt.PointVariableRef, JuMP.VariableRef}

    # Variable support data
    infvar_to_supports::Dict{InfiniteOpt.InfiniteVariableRef, Vector{<:Tuple}}

    # Constraint mapping
    infinite_to_constrs::Dict{InfiniteOpt.InfiniteConstraintRef,
                             Vector{JuMP.ConstraintRef}}
    measure_to_constrs::Dict{InfiniteOpt.MeasureConstraintRef,
                            Vector{JuMP.ConstraintRef}}
    finite_to_constr::Dict{InfiniteOpt.FiniteConstraintRef, JuMP.ConstraintRef}

    # Constraint support data
    infconstr_to_supports::Dict{InfiniteOpt.InfiniteConstraintRef, Vector{<:Tuple}}
    measconstr_to_supports::Dict{InfiniteOpt.MeasureConstraintRef, Vector{<:Tuple}}
    infconstr_to_params::Dict{InfiniteOpt.InfiniteConstraintRef, Tuple}
    measconstr_to_params::Dict{InfiniteOpt.MeasureConstraintRef, Tuple}

    # Default constructor
    function TranscriptionData()
        return new(Dict{InfiniteOpt.InfiniteVariableRef,
                   Vector{JuMP.VariableRef}}(),
                   Dict{InfiniteOpt.GlobalVariableRef, JuMP.VariableRef}(),
                   Dict{InfiniteOpt.PointVariableRef, JuMP.VariableRef}(),
                   Dict{InfiniteOpt.InfiniteVariableRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef,
                        Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.MeasureConstraintRef,
                        Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.FiniteConstraintRef, JuMP.ConstraintRef}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.MeasureConstraintRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.InfiniteConstraintRef, Tuple}(),
                   Dict{InfiniteOpt.MeasureConstraintRef, Tuple}())
    end
end

"""
    TranscriptionModel(args...)::JuMP.Model

Return a JuMP `Model` with `TranscriptionData` included in the extension
data field. Accepts the same arguments as a typical JuMP `Model`.

**Example**
```julia
julia> TranscriptionModel()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
function TranscriptionModel(; kwargs...)::JuMP.Model
    model = JuMP.Model(; kwargs...)
    model.ext[:TransData] = TranscriptionData()
    return model
end
# Accept optimizer_factorys
function TranscriptionModel(optimizer_factory::JuMP.OptimizerFactory;
                            kwargs...)::JuMP.Model
    model = JuMP.Model(optimizer_factory; kwargs...)
    model.ext[:TransData] = TranscriptionData()
    return model
end

"""
    is_transcription_model(model::JuMP.Model)::Bool

Return true if `model` is a `TranscriptionModel` or false otherwise.

**Example**
```julia
julia> is_transcription_model(model)
true
```
"""
function is_transcription_model(model::JuMP.Model)::Bool
    return haskey(model.ext, :TransData)
end

"""
    transcription_data(model::JuMP.Model)::TranscriptionData

Return the `TranscriptionData` from a `TranscriptionModel`. Errors if it is not
a `TranscriptionModel`.
"""
function transcription_data(model::JuMP.Model)::TranscriptionData
    !is_transcription_model(model) && error("Model is not a transcription model.")
    return model.ext[:TransData]
end

"""
    transcription_variable(model::JuMP.Model,
                           vref::InfiniteOpt.InfOptVariableRef)

Return the transcribed variable reference(s) corresponding to `vref`. Errors
if no transcription variable is found.

**Example**
```julia
julia> transcription_variable(trans_model, infvar)
2-element Array{VariableRef,1}:
 infvar(support: 1)
 infvar(support: 2)

julia> transcription_variable(trans_model, gbvar)
gbvar
```
"""
function transcription_variable end

## Define the variable mapping functions
# GlobalVariableRef
function transcription_variable(model::JuMP.Model,
                                vref::InfiniteOpt.GlobalVariableRef)::JuMP.VariableRef
    !haskey(transcription_data(model).global_to_var, vref) && error("Variable " *
                             "reference $vref not used in transcription model.")
    return transcription_data(model).global_to_var[vref]
end
# InfiniteVariableRef
function transcription_variable(model::JuMP.Model,
                                vref::InfiniteOpt.InfiniteVariableRef)::Vector
    !haskey(transcription_data(model).infinite_to_vars, vref) && error("Variable" *
                             "reference $vref not used in transcription model.")
    return transcription_data(model).infinite_to_vars[vref]
end
# PointVariableRef
function transcription_variable(model::JuMP.Model,
                                vref::InfiniteOpt.PointVariableRef)::JuMP.VariableRef
    !haskey(transcription_data(model).point_to_var, vref) && error("Variable " *
                             "reference $vref not used in transcription model.")
    return transcription_data(model).point_to_var[vref]
end

"""
    InfiniteOpt.supports(model::JuMP.Model,
                         vref::InfiniteOpt.InfiniteVariableRef)::Vector

Return the support alias mapping associated with `vref` in the transcribed model.
Errors if `vref` does not have transcribed variables.
"""
function InfiniteOpt.supports(model::JuMP.Model,
                              vref::InfiniteOpt.InfiniteVariableRef)::Vector
    if !haskey(transcription_data(model).infvar_to_supports, vref)
        error("Variable reference $vref not used in transcription model.")
    end
    return transcription_data(model).infvar_to_supports[vref]
end

"""
    InfiniteOpt.supports(vref::InfiniteOpt.InfiniteVariableRef)::Vector

Return the support alias mapping associated with `vref` in the transcription
model. Errors if the infinite model does not contain a transcription model or if
`vref` is not transcribed.

**Example**
```julia
julia> supports(vref)
Dict{Int64,Tuple{Float64}} with 2 entries:
  2 => (1.0,)
  1 => (0.0,)
```
"""
function InfiniteOpt.supports(vref::InfiniteOpt.InfiniteVariableRef)::Vector
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(vref))
    return InfiniteOpt.supports(model, vref)
end

"""
    transcription_constraint(model::JuMP.Model,
                             cref::InfiniteOpt.GeneralConstraintRef)

Return the transcribed constraint reference(s) corresponding to `cref`. Errors
if `cref` has not been transcribed.

**Example**
```julia
julia> transcription_constraint(trans_model, fin_con)
fin_con : x(support: 1) - y <= 3.0
```
"""
function transcription_constraint end

## Define the cosntraint mapping functions
# InfiniteConstraintRef
function transcription_constraint(model::JuMP.Model,
                                  cref::InfiniteOpt.InfiniteConstraintRef)::Vector
    if !haskey(transcription_data(model).infinite_to_constrs, cref)
        error("Constraint reference $cref not used in transcription model.")
    end
    return transcription_data(model).infinite_to_constrs[cref]
end

# MeasureConstraintRef
function transcription_constraint(model::JuMP.Model,
                                  cref::InfiniteOpt.MeasureConstraintRef)::Vector
    if !haskey(transcription_data(model).measure_to_constrs, cref)
        error("Constraint reference $cref not used in transcription model.")
    end
    return transcription_data(model).measure_to_constrs[cref]
end

# FiniteConstraintRef
function transcription_constraint(model::JuMP.Model,
                                  cref::InfiniteOpt.FiniteConstraintRef)::JuMP.ConstraintRef
    if !haskey(transcription_data(model).finite_to_constr, cref)
        error("Constraint reference $cref not used in transcription model.")
    end
    return transcription_data(model).finite_to_constr[cref]
end

"""
    InfiniteOpt.supports(model::JuMP.Model,
                         cref::InfiniteOpt.GeneralConstraintRef)::Vector

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed.
"""
function InfiniteOpt.supports(model::JuMP.Model,
                              cref::InfiniteOpt.InfiniteConstraintRef)::Vector
    if !haskey(transcription_data(model).infconstr_to_supports, cref)
        error("Constraint reference $cref not used in transcription model.")
    end
    return transcription_data(model).infconstr_to_supports[cref]
end
# MeasureConstraintRef
function InfiniteOpt.supports(model::JuMP.Model,
                              cref::InfiniteOpt.MeasureConstraintRef)::Vector
    if !haskey(transcription_data(model).measconstr_to_supports, cref)
        error("Constraint reference $cref not used in transcription model " *
              "and/or is finite and doesn't have supports.")
    end
    return transcription_data(model).measconstr_to_supports[cref]
end

"""
    InfiniteOpt.supports(cref::InfiniteOpt.GeneralConstraintRef)::Vector

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed or if the infinite model does not have a transcription model.

**Example**
```julia
julia> supports(cref)
Dict{Int64,Tuple{Float64}} with 2 entries:
  2 => (1.0,)
  1 => (0.0,)
```
"""
function InfiniteOpt.supports(cref::InfiniteOpt.InfiniteConstraintRef)::Vector
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return InfiniteOpt.supports(model, cref)
end
# MeasureConstraintRef
function InfiniteOpt.supports(cref::InfiniteOpt.MeasureConstraintRef)::Vector
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return InfiniteOpt.supports(model, cref)
end

"""
    InfiniteOpt.parameter_refs(model::JuMP.Model,
                               cref::InfiniteOpt.GeneralConstraintRef)::Tuple

Return the a parameter reference tuple of all the parameters that parameterize
`cref` and correspond to the supports. Errors if `cref` has not been transcribed.
"""
function InfiniteOpt.parameter_refs(model::JuMP.Model,
                                    cref::InfiniteOpt.InfiniteConstraintRef)::Tuple
    if !haskey(transcription_data(model).infconstr_to_params, cref)
        error("Constraint reference $cref not used in transcription model " *
              "and/or is finite and doesn't have parameter references.")
    end
    return transcription_data(model).infconstr_to_params[cref]
end
# MeasureConstraintRef
function InfiniteOpt.parameter_refs(model::JuMP.Model,
                                    cref::InfiniteOpt.MeasureConstraintRef)::Tuple
    if !haskey(transcription_data(model).measconstr_to_params, cref)
        error("Constraint reference $cref not used in transcription model " *
              "and/or is finite and doesn't have parameter references.")
    end
    return transcription_data(model).measconstr_to_params[cref]
end

"""
    InfiniteOpt.parameter_refs(cref::InfiniteOpt.GeneralConstraintRef)::Tuple

Return the a parameter reference tuple of all the parameters that parameterize
`cref` and correspond to the supports. Errors if `cref` has not been transcribed
or if the infinite model does not have a transcription model associated with it.

**Example**
```julia
julia> parameter_refs(cref)
(t, x)
```
"""
function InfiniteOpt.parameter_refs(cref::InfiniteOpt.InfiniteConstraintRef)::Tuple
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return InfiniteOpt.parameter_refs(model, cref)
end
# MeasureConstraintRef
function InfiniteOpt.parameter_refs(cref::InfiniteOpt.MeasureConstraintRef)::Tuple
    model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return InfiniteOpt.parameter_refs(model, cref)
end
