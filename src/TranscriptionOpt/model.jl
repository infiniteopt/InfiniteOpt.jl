################################################################################
#                             BASIC BACKEND DEFINITION
################################################################################
"""
    TranscriptionData

A DataType for storing the data mapping an [`InfiniteOpt.InfiniteModel`](@ref)
that has been transcribed to a regular `JuMP.Model` that contains the
transcribed variables. This is stored in the `data` field of 
[`InfiniteOpt.JuMPBackend`](@ref) to make what is called a `TranscriptionBackend` 
via the [`TranscriptionBackend`](@ref) constructor.

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
- `semi_infinite_vars::Vector{InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}}`:
   Store the core semi-infinite variable objects of semi-infinite variables formed on transcription.
- `semi_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Dict{Int, Float64}}, InfiniteOpt.GeneralVariableRef}`: 
  Lookup which semi-infinite variables have already been added.
- `last_point_index::Int`: The last internal point variable index added.
- `point_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Vector{Float64}}, InfiniteOpt.GeneralVariableRef}`: 
  Lookup which point variables have already been created internally.
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
    semi_infinite_vars::Vector{InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}}
    semi_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Dict{Int, Float64}}, InfiniteOpt.GeneralVariableRef}
    last_point_index::Int
    point_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Vector{Float64}}, InfiniteOpt.GeneralVariableRef}

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
                   Vector{InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}}(),
                   Dict{Tuple{InfiniteOpt.GeneralVariableRef, Dict{Int, Float64}}, InfiniteOpt.GeneralVariableRef}(),
                   0,
                   Dict{Tuple{InfiniteOpt.GeneralVariableRef, Vector{Float64}}, InfiniteOpt.GeneralVariableRef}(),
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
                   (), 
                   (), 
                   false
                   )
    end
end

# Extend Base.empty!
function Base.empty!(data::TranscriptionData)
    empty!(data.infvar_lookup)
    empty!(data.infvar_mappings)
    empty!(data.infvar_supports)
    empty!(data.infvar_support_labels)
    empty!(data.finvar_mappings)
    empty!(data.semi_infinite_vars)
    empty!(data.semi_lookup)
    data.last_point_index = 0
    empty!(data.point_lookup)
    empty!(data.measure_lookup)
    empty!(data.measure_mappings)
    empty!(data.measure_supports)
    empty!(data.measure_support_labels)
    empty!(data.constr_mappings)
    empty!(data.constr_supports)
    empty!(data.constr_support_labels)
    data.supports = ()
    data.support_labels = ()
    data.has_internal_supports = false
    return data
end

"""
    Transcription <: InfiniteOpt.AbstractJuMPTag

Dispatch tag needed for [`TranscriptionBackend`](@ref) to be based on 
[`InfiniteOpt.JuMPBackend`](@ref).
"""
struct Transcription <: InfiniteOpt.AbstractJuMPTag end 

"""
    TranscriptionBackend(
        [optimizer_constructor];
        [add_bridges::Bool = true]
        )::InfiniteOpt.JuMPBackend{Transcription}

Return an `InfiniteOpt.JuMPBackend` that uses [`TranscriptionData`](@ref) 
and the [`Transcription`](@ref) tag. Accepts the same arguments as a typical 
`JuMP.Model`. More detailed variable and constraint naming can be enabled 
via `verbose_naming`.

**Example**
```julia-repl
julia> backend = TranscriptionBackend();
```
"""
const TranscriptionBackend = InfiniteOpt.JuMPBackend{Transcription, Float64, TranscriptionData}

# Constructors
function TranscriptionBackend(; kwargs...)
    model = JuMP.Model(; kwargs...)
    return InfiniteOpt.JuMPBackend{Transcription}(model, TranscriptionData())
end
function TranscriptionBackend(optimizer_constructor; kwargs...)
    model = JuMP.Model(optimizer_constructor; kwargs...)
    return InfiniteOpt.JuMPBackend{Transcription}(model, TranscriptionData())
end

# Printing
function JuMP.show_backend_summary(
    io::IO,
    model::InfiniteOpt.InfiniteModel,
    backend::TranscriptionBackend
    )
    println(io, "  Backend type: TranscriptionBackend")
    # reformulation information
    data = transcription_data(backend)
    supp_tuple = data.supports
    prefs = InfiniteOpt.parameter_refs(model)
    for (i, supps) in enumerate(supp_tuple)
        # support info
        param_name = InfiniteOpt._get_param_group_name(prefs[i])
        println(io, "  `", param_name, "` transcribed over ", length(supps) - 1, " supports")
        # TODO add approximation method info (requires InfiniteOpt refactoring)
    end
    # solver name
    moi_summary = sprint(JuMP.show_backend_summary, backend.model)
    solver_str = filter(startswith("Solver"), split(moi_summary, "\n"))[1]
    println(io, "  ", solver_str)
    return
end

# Showing the backend
function Base.show(io::IO, backend::TranscriptionBackend)
    println(io, "A TranscriptionBackend that uses a")
    show(io, backend.model)
end

################################################################################
#                                BASIC QUERIES
################################################################################
"""
    transcription_data(backend::TranscriptionBackend)::TranscriptionData

Return the mapping data used by `backend`.
"""
function transcription_data(backend::TranscriptionBackend)
    return InfiniteOpt.transformation_data(backend)
end

"""
    has_internal_supports(backend::TranscriptionBackend)::Bool

Return a `Bool` whether `backend` has any internal supports that were collected.
"""
function has_internal_supports(backend::TranscriptionBackend)
    return transcription_data(backend).has_internal_supports
end

################################################################################
#                              VARIABLE QUERIES
################################################################################
# Define method for checking if the label needs to be accounted for 
function _ignore_label(
    backend::TranscriptionBackend, 
    label::Type{<:InfiniteOpt.AbstractSupportLabel}
    )
    return label == InfiniteOpt.All || 
           (!has_internal_supports(backend) && 
           label == InfiniteOpt.PublicLabel)
end

"""
    transcription_variable(
        vref::InfiniteOpt.GeneralVariableRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
         ndarray::Bool = false]
         )

Return the transcribed variable reference(s) corresponding to `vref`. Errors
if no transcription variable is found. Also can query via the syntax:
```julia
transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false]
    )
```
If the infinite model contains a built `TranscriptionBackend`. By default, this
method returns only transcribed variables associated with public supports. All the 
variables can be returned by setting `label = All`. 

If `vref` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of variable 
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_variable(infvar, trans_backend)
2-element Array{VariableRef,1}:
 infvar(support: 1)
 infvar(support: 2)

julia> transcription_variable(hdvar, trans_backend)
hdvar

julia> transcription_variable(infvar)
2-element Array{VariableRef,1}:
 infvar(support: 1)
 infvar(support: 2)

julia> transcription_variable(hdvar)
hdvar
```
"""
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_variable(
        vref,
        InfiniteOpt._index_type(vref),
        backend,
        label,
        ndarray
        )
end

# define convenient aliases
const InfVarIndex = Union{
    InfiniteOpt.InfiniteVariableIndex,
    InfiniteOpt.SemiInfiniteVariableIndex,
    InfiniteOpt.DerivativeIndex
    }
const FinVarIndex = Union{
    InfiniteOpt.FiniteVariableIndex,
    InfiniteOpt.PointVariableIndex
    }

## Define the variable mapping functions
# FinVarIndex
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
    ) where {V <: FinVarIndex}
    var = get(transcription_data(backend).finvar_mappings, vref, nothing)
    if isnothing(var)
        error("Variable reference $vref not used in transcription backend.")
    end
    return var
end

# InfVarIndex
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
    ) where {V <: InfVarIndex}
    vars = get(transcription_data(backend).infvar_mappings, vref, nothing)
    if isnothing(vars)
        error("Variable reference $vref not used in transcription backend.")
    end
    if ndarray 
        return make_ndarray(backend, vref, vars, label)
    elseif _ignore_label(backend, label)
        return vars
    else 
        labels = transcription_data(backend).infvar_support_labels[vref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return vars[inds]
    end
end

# ParameterFunctionIndex
function transcription_variable(
    fref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.ParameterFunctionIndex},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_int_idxs = InfiniteOpt.parameter_group_int_indices(fref)
    support_indices = support_index_iterator(backend, group_int_idxs)
    vals = Vector{Float64}(undef, length(support_indices))
    check_labels = length(vals) > 1 && !_ignore_label(backend, label)
    label_inds = ones(Bool, length(vals))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        supp = index_to_support(backend, idx)
        if check_labels && !any(l -> l <: label, index_to_labels(backend, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds vals[i] = transcription_expression(fref, backend, supp)
    end
    # return the values
    if ndarray
        return make_ndarray(backend, fref, vals, label)
    else
        return vals[label_inds]
    end
end

# Fallback
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type,
    backend::TranscriptionBackend,
    label,
    ndarray
    )
    error("`transcription_variable` not defined for variables with indices of " *
          "type $(index_type) and/or is not defined for labels of type $(label).")
end

# Dispatch for internal backends
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel, 
    ndarray::Bool = false
    )
    return transcription_variable(
        vref,
        JuMP.owner_model(vref).backend,
        label = label,
        ndarray = ndarray
        )
end

"""
    InfiniteOpt.transformation_variable(
        vref::InfiniteOpt.GeneralVariableRef,
        [backend::TranscriptionBackend];
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Proper extension of [`InfiniteOpt.transformation_variable`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_variable`](@ref).
"""
function InfiniteOpt.transformation_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_variable(vref, backend, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.variable_supports(
        vref::InfiniteOpt.DecisionVariableRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Return the support alias mapping associated with `vref` in the transcription backend.
Errors if `vref` does not have transcripted variables. See `transcription_variable` 
for an explanation of `ndarray`.
"""
function InfiniteOpt.variable_supports(
    dvref::Union{
        InfiniteOpt.InfiniteVariableRef,
        InfiniteOpt.SemiInfiniteVariableRef, 
        InfiniteOpt.DerivativeRef
        },
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    vref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(dvref), JuMP.index(dvref))
    if !haskey(transcription_data(backend).infvar_mappings, vref)
        error("Variable reference $vref not used in transcription backend.")
    elseif !haskey(transcription_data(backend).infvar_supports, vref)
        prefs = InfiniteOpt.raw_parameter_refs(dvref)
        lookups = transcription_data(backend).infvar_lookup[vref]
        type = typeof(Tuple(first(keys(lookups)), prefs))
        supps = Vector{type}(undef, length(lookups))
        for (s, i) in lookups
            supps[i] = Tuple(s, prefs)
        end
        transcription_data(backend).infvar_supports[vref] = supps
    end
    supps = transcription_data(backend).infvar_supports[vref]
    if ndarray 
        return make_ndarray(backend, dvref, supps, label)
    elseif _ignore_label(backend, label)
        return supps
    else 
        labels = transcription_data(backend).infvar_support_labels[vref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return supps[inds]
    end
end

# ParameterFunctionRef 
function InfiniteOpt.variable_supports(
    dvref::InfiniteOpt.ParameterFunctionRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_int_idxs = sort(InfiniteOpt.parameter_group_int_indices(dvref))
    support_indices = support_index_iterator(backend, group_int_idxs)
    supps = Vector{Tuple}(undef, length(support_indices))
    check_labels = length(supps) > 1 && !_ignore_label(backend, label)
    param_supps = parameter_supports(backend)
    label_inds = ones(Bool, length(supps))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        if check_labels && !any(l -> l <: label, index_to_labels(backend, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds supps[i] = Tuple(param_supps[j][idx[j]] for j in group_int_idxs)
    end
    # return the supports
    if ndarray
        return make_ndarray(backend, dvref, supps, label)
    else
        return supps[label_inds]
    end
end

"""
    lookup_by_support(
        vref::InfiniteOpt.GeneralVariableRef,
        backend::TranscriptionBackend,
        support::Vector
        )

Return the transcription expression of `vref` defined at its `support`. This is
intended as a helper method for automated transcription.
"""
function lookup_by_support(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend,
    support::Vector
    )
    return lookup_by_support(vref, InfiniteOpt._index_type(vref), backend, support)
end

# define error function for not being able to find a variable by its support
_supp_error() = error("""
    Unable to locate transcription variable by support, consider rebuilding the 
    infinite model with less significant digits. Note this might be due to partially
    evaluating dependent parameters which is not supported by TranscriptionOpt. Such 
    is the case with derivatives/measures that dependent on single dependent 
    parameters.
    """)

# InfiniteIndex
function lookup_by_support(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    support::Vector
    ) where {V <: InfVarIndex}
    if !haskey(transcription_data(backend).infvar_lookup, vref)
        error("Variable reference $vref not used in transcription backend.")
    end
    idx = get(_supp_error, transcription_data(backend).infvar_lookup[vref], support)
    return transcription_data(backend).infvar_mappings[vref][idx]
end

# ParameterFunctionIndex
function lookup_by_support(
    fref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.ParameterFunctionIndex},
    backend::TranscriptionBackend,
    support::Vector
    )
    prefs = InfiniteOpt.raw_parameter_refs(fref)
    func = InfiniteOpt.raw_function(fref)
    return InfiniteOpt.call_function(fref, Tuple(support, prefs)...)
end

# FiniteIndex
function lookup_by_support(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    support::Vector
    ) where {V <: FinVarIndex}
    if !haskey(transcription_data(backend).finvar_mappings, vref)
        error("Variable reference $vref not used in transcription backend.")
    end
    return transcription_data(backend).finvar_mappings[vref]
end

"""
    InfiniteOpt.internal_semi_infinite_variable(
        vref::InfiniteOpt.SemiInfiniteVariableRef,
        backend::TranscriptionBackend
        )::InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}

Return the internal semi-infinite variable associated with `vref`, assuming it was
added internally during measure expansion at the transcription step. This
extends [`InfiniteOpt.internal_semi_infinite_variable`](@ref) as described in its
docstring. Errors, if no such variable can be found.
"""
function InfiniteOpt.internal_semi_infinite_variable(
    vref::InfiniteOpt.SemiInfiniteVariableRef,
    backend::TranscriptionBackend
    )
    semi_infinite_vars = transcription_data(backend).semi_infinite_vars
    idx = -1 * JuMP.index(vref).value
    if idx in keys(semi_infinite_vars)
        return semi_infinite_vars[idx]
    else
        error("Invalid semi-infinite variable reference, this likely is attributed " *
              "to its being deleted.")
    end
end

################################################################################
#                              MEASURE QUERIES
################################################################################
# MeasureIndex
function transcription_variable(
    mref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.MeasureIndex},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel},
    ndarray::Bool = false
    )
    exprs = get(transcription_data(backend).measure_mappings, mref, nothing)
    if isnothing(exprs)
        error("Measure reference $mref not used in transcription backend.")
    end
    if ndarray 
        return make_ndarray(backend, mref, exprs, label)
    elseif length(exprs) > 1 && _ignore_label(backend, label)
        return exprs
    elseif length(exprs) > 1
        labels = transcription_data(backend).measure_support_labels[mref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return exprs[inds]
    else 
        return first(exprs)
    end
end

# Extend transcription_expression
function lookup_by_support(
    mref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.MeasureIndex},
    backend::TranscriptionBackend,
    support::Vector
    )
    if !haskey(transcription_data(backend).measure_lookup, mref)
        error("Measure reference $mref not used in transcription backend.")
    end
    idx = get(_supp_error, transcription_data(backend).measure_lookup[mref], support)
    return transcription_data(backend).measure_mappings[mref][idx]
end

# Extend variable_supports
function InfiniteOpt.variable_supports(
    dmref::InfiniteOpt.MeasureRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    mref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(dmref), JuMP.index(dmref))
    if !haskey(transcription_data(backend).measure_mappings, mref)
        error("Measure reference $mref not used in transcription backend.")
    elseif !haskey(transcription_data(backend).measure_supports, mref)
        lookups = transcription_data(backend).measure_lookup[mref]
        prefs = InfiniteOpt.parameter_refs(dmref)
        vt_prefs = InfiniteOpt.Collections.VectorTuple(prefs)
        type = typeof(Tuple(first(keys(lookups)), vt_prefs))
        supps = Vector{type}(undef, length(lookups))
        for (supp, i) in lookups
            supps[i] = Tuple(supp, vt_prefs)
        end
        transcription_data(backend).measure_supports[mref] = supps
    end
    supps = transcription_data(backend).measure_supports[mref]
    if ndarray
        return make_ndarray(backend, dmref, supps, label)
    elseif length(supps) > 1 && _ignore_label(backend, label)
        return supps
    elseif length(supps) > 1
        labels = transcription_data(backend).measure_support_labels[mref]
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
    transcription_expression(
        expr::JuMP.AbstractJuMPScalar,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Return the transcribed expression(s) corresponding to `expr`. Errors
if `expr` cannot be transcribed. Also can query via the syntax:
```julia
transcription_expression(
    expr::JuMP.AbstractJuMPScalar;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false]
    )
```
If the infinite model contains a built transcription backend. By default, this
method returns only transcribed expressions associated with public supports. All the 
expressions can be returned by setting `label = All`.

If `expr` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_expression(my_expr, backend)
x(support: 1) - y

julia> transcription_expression(my_expr)
x(support: 1) - y
```
"""
function transcription_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_int_idxs = InfiniteOpt.parameter_group_int_indices(expr)
    support_indices = support_index_iterator(backend, group_int_idxs)
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(support_indices))
    check_labels = length(exprs) > 1 && !_ignore_label(backend, label)
    label_inds = ones(Bool, length(exprs))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        supp = index_to_support(backend, idx)
        if check_labels && !any(l -> l <: label, index_to_labels(backend, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds exprs[i] = transcription_expression(expr, backend, supp)
    end
    # return the expressions
    if ndarray
        return make_ndarray(backend, expr, exprs, label)
    else
        exprs = exprs[label_inds]
        return length(support_indices) > 1 ? exprs : first(exprs)
    end
end

# Define for variables
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false)
    return transcription_variable(vref, backend, label = label, ndarray = ndarray)
end

# Dispatch for internal backends
function transcription_expression(
    expr::JuMP.AbstractJuMPScalar; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    model = JuMP.owner_model(expr)
    isnothing(model) && return zero(JuMP.AffExpr) + JuMP.constant(expr)
    return transcription_expression(
        expr,
        model.backend,
        label = label,
        ndarray = ndarray
        )
end

"""
    InfiniteOpt.transformation_expression(
        expr::JuMP.AbstractJuMPScalar,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Proper extension of [`InfiniteOpt.transformation_expression`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_expression`](@ref).
"""
function InfiniteOpt.transformation_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_expression(expr, backend, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.expression_supports(
        expr::JuMP.AbstractJuMPScalar,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Return the support alias mappings associated with `expr`. Errors if `expr` cannot
be transcribed.
"""
function InfiniteOpt.expression_supports(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_int_idxs = sort(InfiniteOpt.parameter_group_int_indices(expr))
    support_indices = support_index_iterator(backend, group_int_idxs)
    supps = Vector{Tuple}(undef, length(support_indices))
    check_labels = length(supps) > 1 && !_ignore_label(backend, label)
    param_supps = parameter_supports(backend)
    label_inds = ones(Bool, length(supps))
    # iterate over the indices and compute the values
    for (i, idx) in enumerate(support_indices)
        if check_labels && !any(l -> l <: label, index_to_labels(backend, idx))
            @inbounds label_inds[i] = false
        end
        @inbounds supps[i] = Tuple(param_supps[j][idx[j]] for j in group_int_idxs)
    end
    # return the supports
    if ndarray
        return make_ndarray(backend, expr, supps, label)
    else
        supps = supps[label_inds]
        return length(support_indices) > 1 ? supps : first(supps)
    end
end

################################################################################
#                             CONSTRAINT QUERIES
################################################################################
"""
    transcription_constraint(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false]
        )

Return the transcribed constraint reference(s) corresponding to `cref`. Errors
if `cref` has not been transcribed. Also can query via the syntax:
```julia
transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false]
    )
```
If the infinite model contains a built transcription backend. By default, this
method returns only transcribed constraints associated with public supports. All the 
constraints can be returned by setting `label = All`.

If `cref` is infinite and `ndarray = true` then an n-dimensional array will be 
returned in accordance with the infinite parameters that have unique object 
numbers. In this case, `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_constraint(fin_con, backend)
fin_con : x(support: 1) - y <= 3.0

julia> transcription_constraint(fin_con)
fin_con : x(support: 1) - y <= 3.0
```
"""
function transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    constr = get(transcription_data(backend).constr_mappings, cref, nothing)
    if isnothing(constr)
        error("Constraint reference $cref not used in transcription backend.")
    end
    if ndarray 
        return make_ndarray(backend, cref, constr, label)
    elseif length(constr) > 1 && _ignore_label(backend, label)
        return constr
    elseif length(constr) > 1
        labels = transcription_data(backend).constr_support_labels[cref]
        inds = map(s -> any(l -> l <: label, s), labels)
        return constr[inds]
    else 
        return first(constr)
    end
end

# Dispatch for internal backends
function transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_constraint(
        cref,
        JuMP.owner_model(cref).backend,
        label = label,
        ndarray = ndarray
        )
end

"""
    InfiniteOpt.transformation_constraint(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
        ndarray::Bool = false]
        )

Proper extension of [`InfiniteOpt.transformation_constraint`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_constraint`](@ref).
"""
function InfiniteOpt.transformation_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    return transcription_constraint(cref, backend, label = label, ndarray = ndarray)
end

"""
    InfiniteOpt.constraint_supports(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
        ndarray::Bool = false])

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed.
"""
function InfiniteOpt.constraint_supports(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel,
    ndarray::Bool = false
    )
    supps = get(transcription_data(backend).constr_supports, cref, nothing)
    if isnothing(supps)
        error("Constraint reference $cref not used in transcription backend.")
    end
    if ndarray 
        return make_ndarray(backend, cref, supps, label)
    elseif length(supps) > 1 && _ignore_label(backend, label)
        return supps
    elseif length(supps) > 1
        labels = transcription_data(backend).constr_support_labels[cref]
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
    parameter_supports(backend::TranscriptionBackend)::Tuple

Return the collected parameter support tuple that is stored in
`TranscriptionData.supports`.
"""
function parameter_supports(backend::TranscriptionBackend)
    return transcription_data(backend).supports
end

"""
    support_index_iterator(backend::TranscriptionBackend, [group_int_idxs::Vector{Int}])::CartesianIndices

Return the `CartesianIndices` that determine the indices of the unique combinations
of `TranscriptionData.supports` stored in `backend`. If `group_int_idxs` is specified,
then the indices will only include the tuple elements uses indices are included
in the parameter group integer indices `group_int_idxs` and all others will be assigned the last index
which should correspond to an appropriately sized placeholder comprised of `NaN`s.
Note this method assumes that [`set_parameter_supports`](@ref) has already been
called and that the last elements of each support vector contains a placeholder
value.
"""
function support_index_iterator(backend::TranscriptionBackend)
    raw_supps = parameter_supports(backend)
    return CartesianIndices(ntuple(i -> 1:length(raw_supps[i])-1, length(raw_supps)))
end

# Generate for a subset of parameter group integer indices (use last index as placeholder --> support with NaNs)
function support_index_iterator(
    backend::TranscriptionBackend,
    group_int_idxs::Vector{Int}
    )
    raw_supps = parameter_supports(backend)
    lens = map(i -> length(i), raw_supps)
    # prepare the indices of each support combo
    # note that the actual supports are from 1:length-1 and the placeholders are at the ends
    return CartesianIndices(ntuple(i -> i in group_int_idxs ? (1:lens[i]-1) : (lens[i]:lens[i]),
                                   length(raw_supps)))
end

"""
    index_to_support(backend::TranscriptionBackend, index::CartesianIndex)::Vector{Float64}

Given a particular support `index` generated via [`support_index_iterator`](@ref)
using `backend`, return the corresponding support from `TranscriptionData.supports`
using placeholder `NaN`s as appropriate for tuple elements that are unneeded.
"""
function index_to_support(
    backend::TranscriptionBackend,
    index::CartesianIndex
    )::Vector{Float64}
    raw_supps = parameter_supports(backend)
    return [j for i in eachindex(index.I) for j in raw_supps[i][index[i]]]
end

"""
    index_to_labels(backend::TranscriptionBackend, index::CartesianIndex)::Set{DataType}

Given a particular support `index` generated via [`support_index_iterator`](@ref)
using `backend`, return the corresponding support label set from `TranscriptionData.support_labels`.
"""
function index_to_labels(
    backend::TranscriptionBackend,
    index::CartesianIndex
    )
    raw_labels = transcription_data(backend).support_labels
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

## Helper functions to consistently get parameter group integer indices 
# Fallback
function _getparameter_group_int_indices(ref)
    return InfiniteOpt.parameter_group_int_indices(ref)
end 

# Expressions
function _getparameter_group_int_indices(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}
    )
    return sort(InfiniteOpt.parameter_group_int_indices(expr))
end 

"""
    make_narray(
        backend::TranscriptionBackend,
        ref::Union{JuMP.AbstractJuMPScalar, InfiniteOpt.InfOptConstraintRef},
        info::Vector,
        label::Type{<:InfiniteOpt.AbstractSupportLabel}
        )::Array 

Take the results `info` associated with `ref` and rearrange them into an 
n-dimensional array where the axes correspond to the infinite parameter dependencies 
in accordance with their creation. Note that this works by querying the object 
numbers. Thus, independent infinite parameters will each get their own dimension 
(even if they are defined at the same time in an array) and each dependent infinite 
parameter group will have its own dimension. 
"""
function make_ndarray(backend::TranscriptionBackend, ref, info::Vector, label::DataType)
    # get the parameter group integer indices
    group_int_idxs = _getparameter_group_int_indices(ref)
    # return result if it is from a finite object
    if isempty(group_int_idxs)
        return info
    end
    # determine the dimensions of the new array
    raw_supps = parameter_supports(backend)
    dims = Tuple(length(raw_supps[i]) - 1 for i in eachindex(raw_supps) if i in group_int_idxs)
    # check that the lengths match (otherwise we'll have some sparse set)
    # TODO add capability to avoid this problem (make reduced array by looking at the supports)
    if length(info) != prod(dims)
        error("Unable to make `ndarray`. This is likely due to the object being " * 
              "over a portion of the infinite-domain (e.g., bounded constraints and " * 
              "certain semi-infinite variables.")
    end
    # make and populate the array
    narray = Array{_get_array_type(info)}(undef, dims)
    for (i, idx) in enumerate(eachindex(narray))
        narray[idx] = info[i]
    end
    # rearrange the array as needed to match the object number order
    sorted_array = issorted(group_int_idxs) ? narray : permutedims(narray, sortperm(group_int_idxs)) 
    # consider the label specified (this will enforce the intersection of labels)
    if _ignore_label(backend, label)
        return sorted_array
    else 
        labels = transcription_data(backend).support_labels[group_int_idxs]
        inds = map(sets -> findall(s -> any(l -> l <: label, s), sets), labels)
        return sorted_array[inds...]
    end
end
