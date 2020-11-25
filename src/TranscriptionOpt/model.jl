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
- `infvar_support_labels::Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}`: 
   Map the infinite variables to their support labels.
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
- `measure_support_labels::Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}`: 
   Map measures to their support labels if they have any.
- `constr_mappings::Dict{InfiniteOpt.InfOptConstraintRef, Vector{JuMP.ConstraintRef}}`:
   Map constraints to their transcriptions.
- `constr_supports::Dict{InfiniteOpt.InfOptConstraintRef, Vector{Tuple}}`:
   Map constraints to their support values.
- `constr_support_labels::Dict{InfiniteOpt.InfOptConstraintRef, Vector{Set{DataType}}}`: 
   Map constraints to their support labels.
- `supports::Tuple`: Store the collected parameter supports here.
- `support_labels::Tuple`: Store the collected parameter labels here.
- `has_internal_supports::Bool`: Where any internal supports collected?
"""
mutable struct TranscriptionData
    # Variable information
    infvar_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.VariableRef}}
    infvar_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}
    infvar_support_labels::Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}

    # Internal variables (created via internal measure expansions)
    reduced_vars::Vector{InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}}
    last_point_index::Int

    # Measure information
    measure_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}
    measure_mappings::Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.AbstractJuMPScalar}}
    measure_supports::Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}
    measure_support_labels::Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}

    # Constraint information
    constr_mappings::Dict{InfiniteOpt.InfOptConstraintRef,
                          Vector{JuMP.ConstraintRef}}
    constr_supports::Dict{InfiniteOpt.InfOptConstraintRef,
                          Vector{Tuple}}
    constr_support_labels::Dict{InfiniteOpt.InfOptConstraintRef,
                               Vector{Set{DataType}}}

    # Collected Supports
    supports::Tuple
    support_labels::Tuple
    has_internal_supports::Bool

    # Default constructor
    function TranscriptionData()
        return new( # variable info
                   Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}(),
                   # internal variables
                   Vector{InfiniteOpt.ReducedVariable{InfiniteOpt.GeneralVariableRef}}(),
                   0,
                   # measure info
                   Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{JuMP.AbstractJuMPScalar}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Tuple}}(),
                   Dict{InfiniteOpt.GeneralVariableRef, Vector{Set{DataType}}}(),
                   # constraint info
                   Dict{InfiniteOpt.InfOptConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfiniteOpt.InfOptConstraintRef, Vector{Vector{Float64}}}(),
                   Dict{InfiniteOpt.InfOptConstraintRef, Vector{Set{DataType}}}(),
                   # support storage
                   (), (), false)
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

"""
    has_internal_supports(model::JuMP.Model)::Bool

Return a `Bool` whether `model` has any internal supports that were collected.
"""
function has_internal_supports(model::JuMP.Model)::Bool
    return transcription_data(model).has_internal_supports
end

################################################################################
#                              VARIABLE QUERIES
################################################################################
# Define method for checking if the label needs to be accounted for 
function _ignore_label(model::JuMP.Model, label::Type{<:InfiniteOpt.AbstractSupportLabel})::Bool
    return label == InfiniteOpt.All || (!has_internal_supports(model) && label == InfiniteOpt.PublicLabel)
end

"""
    transcription_variable(model::JuMP.Model,
        vref::InfiniteOpt.GeneralVariableRef;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
         ndarray::Bool = false])

Return the transcribed variable reference(s) corresponding to `vref`. Errors
if no transcription variable is found. Also can query via the syntax:
```julia
transcription_variable(vref::InfiniteOpt.GeneralVariableRef; 
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
     ndarray::Bool = false])
```
If the infinite model contains a built transcription model. By default, this
method returns only transcribed variables associated with public supports. All the 
variables can be returned by setting `label = All`. 

If `vref` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of variable 
supports that use the label. This is defers from the default behavior which 
considers the union.

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
    vref::InfiniteOpt.GeneralVariableRef;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_variable(model, vref, InfiniteOpt._index_type(vref), label, ndarray)
end

# define convenient aliases
const InfVarIndex = Union{InfiniteOpt.InfiniteVariableIndex,
                          InfiniteOpt.ReducedVariableIndex,
                          InfiniteOpt.DerivativeIndex}
const FinVarIndex = Union{InfiniteOpt.HoldVariableIndex,
                          InfiniteOpt.PointVariableIndex}

## Define the variable mapping functions
# FinVarIndex
function transcription_variable(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
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
    index_type::Type{V},
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
    )::Array{JuMP.VariableRef} where {V <: InfVarIndex}
    vars = get(transcription_data(model).infvar_mappings, vref, nothing)
    if vars === nothing
        error("Variable reference $vref not used in transcription model.")
    end
    if ndarray 
        return make_ndarray(model, vref, vars, label)
    elseif _ignore_label(model, label)
        return vars
    else 
        labels = transcription_data(model).infvar_support_labels[vref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return vars[inds]
    end
end

# ParameterFunctionIndex
function transcription_variable(model::JuMP.Model,
    fref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.ParameterFunctionIndex},
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
    )
    # get the object numbers of the expression and form the support iterator
    obj_nums = InfiniteOpt._object_numbers(fref)
    support_indices = support_index_iterator(model, obj_nums)
    vals = Vector{Float64}(undef, length(support_indices))
    check_labels = length(vals) > 1 && !_ignore_label(model, label)
    label_inds = ones(Bool, length(vals))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        supp = index_to_support(model, idx)
        if check_labels && !any(l -> l <: label, index_to_labels(model, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds vals[i] = transcription_expression(model, fref, supp)
    end
    # return the values
    if ndarray
        return make_ndarray(model, fref, vals, label)
    else
        return vals[label_inds]
    end
end

# Fallback
function transcription_variable(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type,
    label,
    ndarray)
    error("`transcription_variable` not defined for variables with indices of " *
          "type $(index_type) and/or is not defined for labels of type $(label).")
end

# Dispatch for internal models
function transcription_variable(vref::InfiniteOpt.GeneralVariableRef; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel, 
    ndarray::Bool = false
    )
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(vref))
    return transcription_variable(trans_model, vref, label = label, 
                                  ndarray = ndarray)
end

"""
    InfiniteOpt.optimizer_model_variable(vref::InfiniteOpt.GeneralVariableRef,
        ::Val{:TransData};
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Proper extension of [`InfiniteOpt.optimizer_model_variable`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_variable`](@ref).
"""
function InfiniteOpt.optimizer_model_variable(vref::InfiniteOpt.GeneralVariableRef,
    ::Val{:TransData};
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false)
    return transcription_variable(vref, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.variable_supports(model::JuMP.Model,
        vref::InfiniteOpt.DecisionVariableRef,
        key::Val{:TransData} = Val(:TransData);
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the support alias mapping associated with `vref` in the transcription model.
Errors if `vref` does not have transcripted variables. See `transcription_variable` 
for an explanation of `ndarray`.
"""
function InfiniteOpt.variable_supports(model::JuMP.Model,
    dvref::Union{InfiniteOpt.InfiniteVariableRef, InfiniteOpt.ReducedVariableRef, 
                 InfiniteOpt.DerivativeRef},
    key::Val{:TransData} = Val(:TransData);
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )::Array
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
    supps = transcription_data(model).infvar_supports[vref]
    if ndarray 
        return make_ndarray(model, dvref, supps, label)
    elseif _ignore_label(model, label)
        return supps
    else 
        labels = transcription_data(model).infvar_support_labels[vref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return supps[inds]
    end
end

# ParameterFunctionRef 
function InfiniteOpt.variable_supports(model::JuMP.Model,
    dvref::InfiniteOpt.ParameterFunctionRef,
    key::Val{:TransData} = Val(:TransData);
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )::Array
    # get the object numbers of the expression and form the support iterator
    obj_nums = sort(InfiniteOpt._object_numbers(dvref))
    support_indices = support_index_iterator(model, obj_nums)
    supps = Vector{Tuple}(undef, length(support_indices))
    check_labels = length(supps) > 1 && !_ignore_label(model, label)
    param_supps = parameter_supports(model)
    label_inds = ones(Bool, length(supps))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        if check_labels && !any(l -> l <: label, index_to_labels(model, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds supps[i] = Tuple(param_supps[j][idx[j]] for j in obj_nums)
    end
    # return the supports
    if ndarray
        return make_ndarray(model, dvref, supps, label)
    else
        return supps[label_inds]
    end
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
    return lookup_by_support(model, vref, InfiniteOpt._index_type(vref), support)
end

# define error function for not being able to find a variable by its support
_supp_error() = error("Unable to locate transcription variable by support, consider " *
                      "rebuilding the infinite model with less significant digits. " *
                      "Note this might be due to partially evaluating dependent parameters " *
                      "which is not supported by TranscriptionOpt. Such is the case with " * 
                      "derivatives/measures that dependent on single dependent parameters.")

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

# ParameterFunctionIndex
function lookup_by_support(model::JuMP.Model,
    fref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.ParameterFunctionIndex},
    support::Vector
    )
    prefs = InfiniteOpt.raw_parameter_refs(fref)
    func = InfiniteOpt.raw_function(fref)
    return func(Tuple(support, prefs)...)
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
    index_type::Type{InfiniteOpt.MeasureIndex},
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool = false
    )
    exprs = get(transcription_data(model).measure_mappings, mref, nothing)
    if exprs === nothing
        error("Measure reference $mref not used in transcription model.")
    end
    if ndarray 
        return make_ndarray(model, mref, exprs, label)
    elseif length(exprs) > 1 && _ignore_label(model, label)
        return exprs
    elseif length(exprs) > 1
        labels = transcription_data(model).measure_support_labels[mref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return exprs[inds]
    else 
        return first(exprs)
    end
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
    key::Val{:TransData} = Val(:TransData);
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
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
    if ndarray
        return make_ndarray(model, dmref, supps, label)
    elseif length(supps) > 1 && _ignore_label(model, label)
        return supps
    elseif length(supps) > 1
        labels = transcription_data(model).measure_support_labels[mref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return supps[inds]
    else 
        return first(supps)
    end
end

################################################################################
#                             EXPRESSION QUERIES
################################################################################
"""
    transcription_expression(model::JuMP.Model,
        expr::JuMP.AbstractJuMPScalar;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the transcribed expression(s) corresponding to `expr`. Errors
if `expr` cannot be transcribed. Also can query via the syntax:
```julia
transcription_expression(expr::JuMP.AbstractJuMPScalar;
                         [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
                         ndarray::Bool = false])
```
If the infinite model contains a built transcription model. By default, this
method returns only transcribed expressions associated with public supports. All the 
expressions can be returned by setting `label = All`.

If `expr` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_expression(trans_model, my_expr)
x(support: 1) - y

julia> transcription_expression(my_expr)
x(support: 1) - y
```
"""
function transcription_expression(model::JuMP.Model,
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr};
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    # get the object numbers of the expression and form the support iterator
    obj_nums = InfiniteOpt._object_numbers(expr)
    support_indices = support_index_iterator(model, obj_nums)
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(support_indices))
    check_labels = length(exprs) > 1 && !_ignore_label(model, label)
    label_inds = ones(Bool, length(exprs))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        supp = index_to_support(model, idx)
        if check_labels && !any(l -> l <: label, index_to_labels(model, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds exprs[i] = transcription_expression(model, expr, supp)
    end
    # return the expressions
    if ndarray
        return make_ndarray(model, expr, exprs, label)
    else
        exprs = exprs[label_inds]
        return length(support_indices) > 1 ? exprs : first(exprs)
    end
end

# Define for variables
function transcription_expression(model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false)
    return transcription_variable(model, vref, label = label, ndarray = ndarray)
end

# Dispatch for internal models
function transcription_expression(expr::JuMP.AbstractJuMPScalar; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    model = InfiniteOpt._model_from_expr(expr)
    if model === nothing
        return zero(JuMP.AffExpr) + JuMP.constant(expr)
    else
        trans_model = InfiniteOpt.optimizer_model(model)
    end
    return transcription_expression(trans_model, expr, label = label, 
                                    ndarray = ndarray)
end

"""
    InfiniteOpt.optimizer_model_expression(expr::JuMP.AbstractJuMPScalar,
        ::Val{:TransData};
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Proper extension of [`InfiniteOpt.optimizer_model_expression`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_expression`](@ref).
"""
function InfiniteOpt.optimizer_model_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
    ::Val{:TransData};
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false)
    return transcription_expression(expr, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.expression_supports(model::JuMP.Model,
        expr::JuMP.AbstractJuMPScalar,
        key::Val{:TransData} = Val(:TransData);
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the support alias mappings associated with `expr`. Errors if `expr` cannot
be transcribed.
"""
function InfiniteOpt.expression_supports(model::JuMP.Model,
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
    key::Val{:TransData} = Val(:TransData);
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    # get the object numbers of the expression and form the support iterator
    obj_nums = sort(InfiniteOpt._object_numbers(expr))
    support_indices = support_index_iterator(model, obj_nums)
    supps = Vector{Tuple}(undef, length(support_indices))
    check_labels = length(supps) > 1 && !_ignore_label(model, label)
    param_supps = parameter_supports(model)
    label_inds = ones(Bool, length(supps))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        if check_labels && !any(l -> l <: label, index_to_labels(model, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds supps[i] = Tuple(param_supps[j][idx[j]] for j in obj_nums)
    end
    # return the supports
    if ndarray
        return make_ndarray(model, expr, supps, label)
    else
        supps = supps[label_inds]
        return length(support_indices) > 1 ? supps : first(supps)
    end
end

################################################################################
#                             CONSTRAINT QUERIES
################################################################################
"""
    transcription_constraint(model::JuMP.Model,
        cref::InfiniteOpt.InfOptConstraintRef;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the transcribed constraint reference(s) corresponding to `cref`. Errors
if `cref` has not been transcribed. Also can query via the syntax:
```julia
transcription_constraint(cref::InfiniteOpt.InfOptConstraintRef;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false])
```
If the infinite model contains a built transcription model. By default, this
method returns only transcribed constraints associated with public supports. All the 
constraints can be returned by setting `label = All`.

If `cref` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_constraint(trans_model, fin_con)
fin_con : x(support: 1) - y <= 3.0

julia> transcription_constraint(fin_con)
fin_con : x(support: 1) - y <= 3.0
```
"""
function transcription_constraint(model::JuMP.Model,
    cref::InfiniteOpt.InfOptConstraintRef;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    constr = get(transcription_data(model).constr_mappings, cref, nothing)
    if constr === nothing
      error("Constraint reference $cref not used in transcription model.")
    end
    if ndarray 
        return make_ndarray(model, cref, constr, label)
    elseif length(constr) > 1 && _ignore_label(model, label)
        return constr
    elseif length(constr) > 1
        labels = transcription_data(model).constr_support_labels[cref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return constr[inds]
    else 
        return first(constr)
    end
end

# Dispatch for internal models
function transcription_constraint(cref::InfiniteOpt.InfOptConstraintRef;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    trans_model = InfiniteOpt.optimizer_model(JuMP.owner_model(cref))
    return transcription_constraint(trans_model, cref, label = label, 
                                    ndarray = ndarray)
end

"""
    InfiniteOpt.optimizer_model_constraint(cref::InfiniteOpt.InfOptConstraintRef,
        ::Val{:TransData};
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
        ndarray::Bool = false])

Proper extension of [`InfiniteOpt.optimizer_model_constraint`](@ref) for
`TranscriptionModel`s. This simply dispatches to [`transcription_constraint`](@ref).
"""
function InfiniteOpt.optimizer_model_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    ::Val{:TransData};
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_constraint(cref, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.constraint_supports(model::JuMP.Model,
        cref::InfiniteOpt.InfOptConstraintRef,
        key::Val{:TransData} = Val(:TransData);
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed.
"""
function InfiniteOpt.constraint_supports(model::JuMP.Model,
    cref::InfiniteOpt.InfOptConstraintRef,
    key::Val{:TransData} = Val(:TransData);
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    supps = get(transcription_data(model).constr_supports, cref, nothing)
    if supps === nothing
        error("Constraint reference $cref not used in transcription model.")
    end
    if ndarray 
        return make_ndarray(model, cref, supps, label)
    elseif length(supps) > 1 && _ignore_label(model, label)
        return supps
    elseif length(supps) > 1
        labels = transcription_data(model).constr_support_labels[cref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return supps[inds]
    else 
        return first(supps)
    end
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
    return CartesianIndices(ntuple(i -> 1:length(raw_supps[i])-1, length(raw_supps)))
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

"""
    index_to_labels(model::JuMP.Model, index::CartesianIndex)::Set{DataType}

Given a particular support `index` generated via [`support_index_iterator`](@ref)
using `model`, return the corresponding support label set from `TranscriptionData.support_labels`.
"""
function index_to_labels(model::JuMP.Model,
                         index::CartesianIndex)::Set{DataType}
    raw_labels = transcription_data(model).support_labels
    labels = Set{DataType}()
    for (i, j) in enumerate(index.I)
        union!(labels, raw_labels[i][j])
    end
    return labels
end

################################################################################
#                               QUERY FORMATERS
################################################################################
# Helper function for getting the array type T 
function _get_array_type(array::Array{T, N}) where {T, N}
    return T
end

## Helper functions to consistently get object numbers 
# Fallback
function _get_object_numbers(ref)
    return InfiniteOpt._object_numbers(ref)
end 

# Expressions
function _get_object_numbers(expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr})
    return sort(InfiniteOpt._object_numbers(expr))
end 

"""
    make_narray(model::JuMP.Model, 
                ref::Union{JuMP.AbstractJuMPScalar, InfiniteOpt.InfOptConstraintRef}, 
                info::Vector, 
                label::Type{<:InfiniteOpt.AbstractSupportLabel})::Array 

Take the results`info` associated with `ref` and rearrange them into an 
n-dimensional array where the axes correspond to the infinite parameter dependencies 
in accordance with their creation. Note that this works by querying the object 
numbers. Thus, independent infinite parameters will each get their own dimension 
(even if they are defined at the same time in an array) and each dependent infinite 
parameter group will have its own dimension. 
"""
function make_ndarray(model::JuMP.Model, ref, info::Vector, label::DataType)
    # get the object numbers
    obj_nums = _get_object_numbers(ref)
    # return result if it is from a finite object
    if isempty(obj_nums)
        return info
    end
    # determine the dimensions of the new array
    raw_supps = parameter_supports(model)
    dims = Tuple(length(raw_supps[i]) - 1 for i in eachindex(raw_supps) if i in obj_nums)
    # check that the lengths match (otherwise we'll have some sparse set)
    # TODO add capability to avoid this problem (make reduced array by looking at the supports)
    if length(info) != prod(dims)
        error("Unable to make `ndarray`. This is likely due to the object being " * 
              "over a portion of the infinite-domain (e.g., bounded constraints and " * 
              "certain reduced infinite variables.")
    end
    # make and populate the array
    narray = Array{_get_array_type(info)}(undef, dims)
    for (i, idx) in enumerate(eachindex(narray))
        narray[idx] = info[i]
    end
    # rearrange the array as needed to match the object number order
    sorted_array = issorted(obj_nums) ? narray : permutedims(narray, sortperm(obj_nums)) 
    # consider the label specified (this will enforce the intersection of labels)
    if _ignore_label(model, label)
        return sorted_array
    else 
        labels = transcription_data(model).support_labels[obj_nums]
        inds = map(sets -> findall(s -> any(l -> l <: label, s), sets), labels)
        return sorted_array[inds...]
    end
end
