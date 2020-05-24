################################################################################
#                              BASIC MODEL DEFINITION
################################################################################
"""
    TranscriptionData

A DataType for storing the data mapping an [`InfiniteOpt.InfiniteModel`](@ref)
that has been transcribed to a regular [`JuMP.Model`](@ref) that contains the
transcribed variables. This is stored in the `ext` field of a `JuMP.Model` to
make what is called a `TranscriptionModel` via the [`TranscriptionModel`](@ref)
constructor.

**Fields**
- `infvar_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}`:
   A lookup table of infinite variable transcriptions via support value.
- `infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.VariableRef}}`:
   Map infinite variables to their transcription variables.
- `infvar_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}`:
   Map infinite variables to their support values.
- `finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}`:
   Map finite variables to their transcription variables.
- `reduced_vars::Vector{InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}}`:
   Store the core reduced variable objects of reduced variables formed on transcription.
- `last_point_index::Int`: The last internal point variable index added.
- `measure_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}`:
   A lookup table of measure transcriptions via support value.
- `measure_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.AbstractJuMPScalar}}`:
   Map measures to transcription expressions.
- `measure_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}`:
   Map measures to their supports values (if the transcribed measure is still infinite).
- `constr_mappings::Dict{InfiniteOpt.InfOptConstraintRef, Vector{JuMP.ConstraintRef}}`:
   Map constraints to their transcriptions.
- `constr_supports::Dict{InfiniteOpt.InfOptConstraintRef, Vector{Tuple}}`:
   Map constraints to their support values.
- `supports::Tuple`: Store the collected parameter supports here.
"""
mutable struct TranscriptionData
    # Variable information
    infvar_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.VariableRef}}
    infvar_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}

    # Internal variables (created via internal measure expansions)
    reduced_vars::Vector{InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}}
    last_point_index::Int

    # Measure information
    measure_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}
    measure_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.AbstractJuMPScalar}}
    measure_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}

    # Constraint information
    constr_mappings::Dict{InfiniteOpt.InfOptConstraintRef{JuMP.ScalarShape},
                          Vector{JuMP.ConstraintRef}}
    constr_supports::Dict{InfiniteOpt.InfOptConstraintRef{JuMP.ScalarShape},
                          Vector{Tuple}}

    # Collected Supports
    supports::Tuple

    # Default constructor
    function TranscriptionData()
        return new( # variable info
                   Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}(),
                   # internal variables
                   Vector{InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}}(),
                   0,
                   # measure info
                   Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.AbstractJuMPScalar}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}(),
                   # constraint info
                   Dict{InfiniteOpt.InfOptConstraintRef{JuMP.ScalarShape},
                        Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.InfOptConstraintRef{JuMP.ScalarShape},
                        Vector{Vector{Float64}}}(),
                   # support storage
                   ())
    end
end

"""
    TranscriptionModel([optimizer_constructor;
                       caching_mode::MOIU.CachingOptimizerMode = MOIU.AUTOMATIC,
                       bridge_constraints::Bool = true])::JuMP.Model

Return a [`JuMP.Model`](@ref) with [`TranscriptionData`](@ref) included in the
`ext` data field. Accepts the same arguments as a typical JuMP `Model`.
More detailed variable and constraint naming can be enabled via `verbose_naming`.

**Example**
```julia-repl
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
# Accept optimizer constructors
function TranscriptionModel(optimizer_constructor;
                            kwargs...)::JuMP.Model
    model = JuMP.Model(optimizer_constructor; kwargs...)
    model.ext[:TransData] = TranscriptionData()
    return model
end

################################################################################
#                                BASIC QUERIES
################################################################################
"""
    is_transcription_model(model::JuMP.Model)::Bool

Return true if `model` is a `TranscriptionModel` or false otherwise.

**Example**
```julia-repl
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

################################################################################
#                              VARIABLE QUERIES
################################################################################
"""
    transcription_variable(model::JuMP.Model,
                           vref::InfiniteOpt.GeneralVariableRef)

Return the transcribed variable reference(s) corresponding to `vref`. Errors
if no transcription variable is found. Also can query via the syntax:
```julia
transcription_variable(vref::InfiniteOpt.GeneralVariableRef)
```
If the infinite model contains a built transcription model.

**Example**
```julia-repl
julia> transcription_variable(trans_model, infvar)
2-element Array{VariableRef,1}:
 infvar(support: 1)
 infvar(support: 2)

julia> transcription_variable(trans_model, hdvar)
hdvar

julia> transcription_variable(infvar)
2-element Array{VariableRef,1}:
 infvar(support: 1)
 infvar(support: 2)

julia> transcription_variable(hdvar)
hdvar
```
"""
function transcription_variable(model::JuMP.Model,
                                vref::InfiniteOpt.GeneralVariableRef)
    return transcription_variable(model, vref, InfiniteOpt._index_type(vref))
end

# define convenient aliases
const InfVarIndex = Union{InfiniteOpt.InfiniteVariableIndex,
                          InfiniteOpt.ReducedVariableIndex}
const FinVarIndex = Union{InfiniteOpt.HoldVariableIndex,
                          InfiniteOpt.PointVariableIndex}

## Define the variable mapping functions
# FinVarIndex
function transcription_variable(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V}
    )::JuMP.VariableRef where {V <: FinVarIndex}
    var = get(transcription_data(model).finvar_mappings, vref, nothing)
    if var === nothing
        error("Variable reference $vref not used in transcription model.")
    end
    return var
end

# InfVarIndex
function transcription_variable(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V}
    )::Vector{JuMP.VariableRef} where {V <: InfVarIndex}
    var = get(transcription_data(model).infvar_mappings, vref, nothing)
    if var === nothing
        error("Variable reference $vref not used in transcription model.")
    end
    return var
end

# Fallback
function transcription_variable(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type)
    error("`transcription_variable` not defined for variables with indices of " *
          "type $(index_type).")
end

# Dispatch for internal models
function transcription_variable(vref::InfiniteOpt.GeneralVariableRef)
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(vref))
    return transcription_variable(trans_model, vref)
end

"""
    InfiniteOpt.optimizer_model_variable(vref::InfiniteOpt.GeneralVariableRef,
                                         ::Val{:TransData})

Proper extension of [`InfiniteOpt.optimizer_model_variable`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_variable`](@ref).
"""
function InfiniteOpt.optimizer_model_variable(vref::InfiniteOpt.GeneralVariableRef,
                                              ::Val{:TransData})
    return transcription_variable(vref)
end

"""
    InfiniteOpt.variable_supports(model::JuMP.Model,
                                  vref::InfiniteOpt.DecisionVariableRef,
                                  key::Val{:TransData} = Val(:TransData))

Return the support alias mapping associated with `vref` in the transcription model.
Errors if `vref` does not have transcripted variables.
"""
function InfiniteOpt.variable_supports(model::JuMP.Model,
    dvref::Union{InfiniteOpt.InfiniteVariableRef, InfiniteOpt.ReducedVariableRef},
    key::Val{:TransData} = Val(:TransData))::Vector
    vref = InfiniteOpt._make_variable_ref(JuMP.owner_model(dvref), JuMP.index(dvref))
    if !haskey(transcription_data(model).infvar_mappings, vref)
        error("Variable reference $vref not used in transcription model.")
    elseif !haskey(transcription_data(model).infvar_supports, vref)
        prefs = InfiniteOpt.raw_parameter_refs(dvref)
        lookups = transcription_data(model).infvar_lookup[vref]
        type = typeof(Tuple(first(keys(lookups)), prefs))
        supps = Vector{type}(undef, length(lookups))
        for (s, i) in lookups
            supps[i] = Tuple(s, prefs)
        end
        transcription_data(model).infvar_supports[vref] = supps
    end
    return transcription_data(model).infvar_supports[vref]
end

"""
    lookup_by_support(model::JuMP.Model,
                      vref::InfiniteOpt.GeneralVariableRef,
                      support::Vector)

Return the transcription expression of `vref` defined at its `support`. This is
intended as a helper method for automated transcription.
"""
function lookup_by_support(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    support::Vector)
    return transcription_expression(model, vref, InfiniteOpt._index_type(vref),
                                    support)
end

# define error function for not being able to find a variable by its support
_supp_error() = error("Unable to locate transcription variable by support, consider " *
                      "rebuilding the infinite model with less significant digits.")

# InfiniteIndex
function lookup_by_support(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector
    )::JuMP.VariableRef where {V <: InfVarIndex}
    if !haskey(transcription_data(model).infvar_lookup, vref)
        error("Variable reference $vref not used in transcription model.")
    end
    idx = get(_supp_error, transcription_data(model).infvar_lookup[vref], support)
    return transcription_data(model).infvar_mappings[vref][idx]
end

# FiniteVariableIndex
function lookup_by_support(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector
    )::JuMP.VariableRef where {V <: FinVarIndex}
    if !haskey(transcription_data(model).finvar_mappings, vref)
        error("Variable reference $vref not used in transcription model.")
    end
    return transcription_data(model).finvar_mappings[vref]
end

"""
    InfiniteOpt.internal_reduced_variable(
        vref::InfiniteOpt.ReducedVariableRef,
        ::Val{:TransData}
        )::InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}

Return the internal reduced variable associated with `vref`, assuming it was
added internally during measure expansion at the transcription step. This
extends [`InfiniteOpt.internal_reduced_variable`](@ref) as described in its
docstring. Errors, if no such variable can be found.
"""
function InfiniteOpt.internal_reduced_variable(
    vref::InfiniteOpt.ReducedVariableRef,
    ::Val{:TransData}
    )::InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(vref))
    reduced_vars = transcription_data(trans_model).reduced_vars
    idx = -1 * JuMP.index(vref).value
    if idx in keys(reduced_vars)
        return reduced_vars[idx]
    else
        error("Invalid reduced variable reference, this likely is attributed " *
              "to its being deleted.")
    end
end

################################################################################
#                              MEASURE QUERIES
################################################################################
# MeasureIndex
function transcription_variable(model::JuMP.Model,
    mref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.MeasureRef}
    )
    exprs = get(transcription_data(model).measure_mappings, mref, nothing)
    if exprs === nothing
        error("Measure reference $mref not used in transcription model.")
    end
    return length(exprs) > 1 ? exprs : first(exprs)
end

# Extend transcription_expression
function lookup_by_support(model::JuMP.Model,
    mref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.MeasureIndex},
    support::Vector
    )::JuMP.AbstractJuMPScalar
    if !haskey(transcription_data(model).measure_lookup, mref)
        error("Measure reference $mref not used in transcription model.")
    end
    idx = get(_supp_error, transcription_data(model).measure_lookup[mref], support)
    return transcription_data(model).measure_mappings[mref][idx]
end

# Extend variable_supports
function InfiniteOpt.variable_supports(model::JuMP.Model,
    dmref::InfiniteOpt.MeasureRef,
    key::Val{:TransData} = Val(:TransData)
    )
    mref = InfiniteOpt._make_variable_ref(JuMP.owner_model(dmref), JuMP.index(dmref))
    if !haskey(transcription_data(model).measure_mappings, mref)
        error("Measure reference $mref not used in transcription model.")
    elseif !haskey(transcription_data(model).measure_supports, mref)
        lookups = transcription_data(model).measure_lookup[mref]
        prefs = InfiniteOpt.parameter_refs(dmref)
        vt_prefs = InfiniteOpt.Collections.VectorTuple(prefs)
        type = typeof(Tuple(first(keys(lookups)), vt_prefs))
        supps = Vector{type}(undef, length(lookups))
        for (supp, i) in lookups
            supps[i] = Tuple(supp, vt_prefs)
        end
        transcription_data(model).measure_supports[mref] = supps
    end
    supps = transcription_data(model).measure_supports[mref]
    return length(supps) > 1 ? supps : first(supps)
end

################################################################################
#                             EXPRESSION QUERIES
################################################################################
"""
    transcription_expression(model::JuMP.Model,
                             expr::JuMP.AbstractJuMPScalar)

Return the transcribed expression(s) corresponding to `expr`. Errors
if `expr` cannot be transcribed. Also can query via the syntax:
```julia
transcription_expression(expr::JuMP.AbstractJuMPScalar)
```
If the infinite model contains a built transcription model.

**Example**
```julia-repl
julia> transcription_expression(trans_model, my_expr)
fin_con : x(support: 1) - y <= 3.0

julia> transcription_expression(my_expr)
fin_con : x(support: 1) - y <= 3.0
```
"""
function transcription_expression(model::JuMP.Model,
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}
    )
    # get the object numbers of the expression and form the support iterator
    obj_nums = InfiniteOpt._object_numbers(expr)
    support_indices = support_index_iterator(model, obj_nums)
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(support_indices))
    counter = 1
    # iterate over the indices and compute the values
    for i in support_indices
        supp = index_to_support(model, i)
        @inbounds exprs[counter] = transcription_expression(model, expr, supp)
        counter += 1
    end
    # return the expressions
    return length(exprs) > 1 ? exprs : first(exprs)
end

# Define for variables
function transcription_expression(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef)
    return transcription_variable(model, vref)
end

# Dispatch for internal models
function transcription_expression(expr::JuMP.AbstractJuMPScalar)
    model = _model_from_expr(expr)
    if model === nothing
        return zero(JuMP.AffExpr) + JuMP.constant(expr)
    else
        trans_model = InfiniteOpt.optimizer_model(model)
    end
    return transcription_expression(trans_model, expr)
end

"""
    InfiniteOpt.optimizer_model_expression(expr::JuMP.AbstractJuMPScalar,
                                           ::Val{:TransData})

Proper extension of [`InfiniteOpt.optimizer_model_expression`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_expression`](@ref).
"""
function InfiniteOpt.optimizer_model_expression(
    expr::JuMP.AbstractJuMPScalar,
    ::Val{:TransData})
    return transcription_expression(expr)
end

"""
    InfiniteOpt.expression_supports(model::JuMP.Model,
                                    expr::JuMP.AbstractJuMPScalar,
                                    key::Val{:TransData} = Val(:TransData))

Return the support alias mappings associated with `expr`. Errors if `expr` cannot
be transcribed.
"""
function InfiniteOpt.expression_supports(model::JuMP.Model,
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
    key::Val{:TransData} = Val(:TransData))
    # get the object numbers of the expression and form the support iterator
    obj_nums = sort(InfiniteOpt._object_numbers(expr))
    support_indices = support_index_iterator(model, obj_nums)
    supps = Vector{Tuple}(undef, length(support_indices))
    param_supps = parameter_supports(model)
    counter = 1
    # iterate over the indices and compute the values
    for i in support_indices
        @inbounds supps[counter] = Tuple(param_supps[j][i[j]] for j in obj_nums)
        counter += 1
    end
    # return the supports
    return length(supps) > 1 ? supps : first(supps)
end

################################################################################
#                             CONSTRAINT QUERIES
################################################################################
"""
    transcription_constraint(model::JuMP.Model,
                             cref::InfiniteOpt.InfOptConstraintRef)

Return the transcribed constraint reference(s) corresponding to `cref`. Errors
if `cref` has not been transcribed. Also can query via the syntax:
```julia
transcription_constraint(cref::InfiniteOpt.InfOptConstraintRef)
```
If the infinite model contains a built transcription model.

**Example**
```julia-repl
julia> transcription_constraint(trans_model, fin_con)
fin_con : x(support: 1) - y <= 3.0

julia> transcription_constraint(fin_con)
fin_con : x(support: 1) - y <= 3.0
```
"""
function transcription_constraint(model::JuMP.Model,
                                  cref::InfiniteOpt.InfOptConstraintRef)
    constr = get(transcription_data(model).constr_mappings, cref, nothing)
    if constr === nothing
      error("Constraint reference $cref not used in transcription model.")
    end
    return length(constr) > 1 ? constr : first(constr)
end

# Dispatch for internal models
function transcription_constraint(cref::InfiniteOpt.InfOptConstraintRef)
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return transcription_constraint(trans_model, cref)
end

"""
    InfiniteOpt.optimizer_model_constraint(cref::InfiniteOpt.InfOptConstraintRef,
                                           ::Val{:TransData})

Proper extension of [`InfiniteOpt.optimizer_model_constraint`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_constraint`](@ref).
"""
function InfiniteOpt.optimizer_model_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    ::Val{:TransData})
    return transcription_constraint(cref)
end

"""
    InfiniteOpt.constraint_supports(model::JuMP.Model,
                                    cref::InfiniteOpt.InfOptConstraintRef,
                                    key::Val{:TransData} = Val(:TransData))

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed.
"""
function InfiniteOpt.constraint_supports(model::JuMP.Model,
                                         cref::InfiniteOpt.InfOptConstraintRef,
                                         key::Val{:TransData} = Val(:TransData))
    supps = get(transcription_data(model).constr_supports, cref, nothing)
    if supps === nothing
        error("Constraint reference $cref not used in transcription model.")
    end
    return length(supps) > 1 ? supps : first(supps)
end

################################################################################
#                             OTHER QUERIES
################################################################################
"""
    parameter_supports(model::JuMP.Model)::Tuple

Return the collected parameter support tuple that is stored in
`TranscriptionData.supports`.
"""
function parameter_supports(model::JuMP.Model)::Tuple
    return transcription_data(model).supports
end

"""
    support_index_iterator(model::JuMP.Model, [obj_nums::Vector{Int}])::CartesianIndices

Return the `CartesianIndices` that determine the indices of the unique combinations
of `TranscriptionData.supports` stored in `model`. If `obj_nums` is specified,
then the indices will only include the tuple elements uses indices are included
in the object numbers `obj_nums` and all others will be assigned the last index
which should correspond to an appropriately sized placeholder comprised of `NaN`s.
Note this method assumes that [`set_parameter_supports`](@ref) has already been
called and that the last elements of each support vector contains a placeholder
value.
"""
function support_index_iterator(model::JuMP.Model)::CartesianIndices
    raw_supps = parameter_supports(model)
    return CartesianIndices(ntuple(i -> 1:length(raw_supps[i]-1), length(raw_supps)))
end

# Generate for a subset of object numbers (use last index as placeholder --> support with NaNs)
function support_index_iterator(model::JuMP.Model,
                                obj_nums::Vector{Int})::CartesianIndices
    raw_supps = parameter_supports(model)
    lens = map(i -> length(i), raw_supps)
    # prepare the indices of each support combo
    # note that the actual supports are afrom 1:length-1 and the placeholders are at the ends
    return CartesianIndices(ntuple(i -> i in obj_nums ? (1:lens[i]-1) : (lens[i]:lens[i]),
                                   length(raw_supps)))
end

"""
    index_to_support(model::JuMP.Model, index::CartesianIndex)::Vector{Float64}

Given a particular support `index` generated via [`support_index_iterator`](@ref)
using `model`, return the corresponding support from `TranscriptionData.supports`
using placeholder `NaN`s as appropriate for tuple elements that are unneeded.
"""
function index_to_support(model::JuMP.Model,
                          index::CartesianIndex)::Vector{Float64}
    raw_supps = parameter_supports(model)
    return [j for i in eachindex(index.I) for j in raw_supps[i][index[i]]]
end
