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
"""
mutable struct TranscriptionData
    # Variable information
    infvar_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, JuMP.VariableRef}}
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, Array{JuMP.VariableRef}}
    infvar_supports::Dict{InfiniteOpt.GeneralVariableRef, Array{Tuple}}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}

    # Metadata
    valid_indices::Dict{Any, Array{Bool}}

    # Internal variables (created via internal measure expansions)
    semi_infinite_vars::Vector{InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}}
    semi_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Dict{Int, Float64}}, InfiniteOpt.GeneralVariableRef}
    last_point_index::Int
    point_lookup::Dict{Tuple{InfiniteOpt.GeneralVariableRef, Vector{Float64}}, InfiniteOpt.GeneralVariableRef}

    # Measure information
    measure_lookup::Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}
    measure_mappings::Dict{InfiniteOpt.GeneralVariableRef, Array{JuMP.AbstractJuMPScalar}}
    measure_supports::Dict{InfiniteOpt.GeneralVariableRef, Array{Tuple}}

    # Constraint information
    constr_mappings::Dict{InfiniteOpt.InfOptConstraintRef,
                          Array{JuMP.ConstraintRef}}
    constr_supports::Dict{InfiniteOpt.InfOptConstraintRef,
                          Array{Tuple}}

    # Collected Supports
    supports::Tuple
    support_labels::Tuple
    has_internal_supports::Bool

    # Default constructor
    function TranscriptionData()
        return new( 
            # variable info
            Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, JuMP.VariableRef}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Array{JuMP.VariableRef}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Array{Tuple}}(),
            Dict{InfiniteOpt.GeneralVariableRef, JuMP.VariableRef}(),
            # meta data
            Dict{Any, Array{Bool}}(),
            # internal variables
            Vector{InfiniteOpt.SemiInfiniteVariable{InfiniteOpt.GeneralVariableRef}}(),
            Dict{Tuple{InfiniteOpt.GeneralVariableRef, Dict{Int, Float64}}, InfiniteOpt.GeneralVariableRef}(),
            0,
            Dict{Tuple{InfiniteOpt.GeneralVariableRef, Vector{Float64}}, InfiniteOpt.GeneralVariableRef}(),
            # measure info
            Dict{InfiniteOpt.GeneralVariableRef, Dict{Vector{Float64}, Int}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Array{JuMP.AbstractJuMPScalar}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Array{Tuple}}(),
            # constraint info
            Dict{InfiniteOpt.InfOptConstraintRef, Array{JuMP.ConstraintRef}}(),
            Dict{InfiniteOpt.InfOptConstraintRef, Array{Tuple}}(),
            # support storage
            (), 
            (), 
            false,
            )
    end
end

# Extend Base.empty!
function Base.empty!(data::TranscriptionData)
    empty!(data.infvar_lookup)
    empty!(data.infvar_mappings)
    empty!(data.infvar_supports)
    empty!(data.finvar_mappings)
    empty!(data.valid_indices)
    empty!(data.semi_infinite_vars)
    empty!(data.semi_lookup)
    data.last_point_index = 0
    empty!(data.point_lookup)
    empty!(data.measure_lookup)
    empty!(data.measure_mappings)
    empty!(data.measure_supports)
    empty!(data.constr_mappings)
    empty!(data.constr_supports)
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

# Get the solver name from MOI
# Inspired by https://github.com/jump-dev/JuMP.jl/blob/ce946b7092c45bdac916c9b531a13a5b929d45f0/src/print.jl#L281-L291
function _try_solver_name(model)
    if mode(model) != JuMP.DIRECT &&
       MOI.Utilities.state(backend(model)) == MOI.Utilities.NO_OPTIMIZER
        return "none"
    end
    try
        return MOI.get(backend(model), MOI.SolverName())
    catch
        return "unknown"
    end
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
    println(io, "  Solver: ", _try_solver_name(backend.model))
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

## truncate a collection according to a label
# 0-Array
function _truncate_by_label(
    arr::Array{T, 0},
    labels::Tuple{},
    label,
    ::Nothing
    ) where {T}
    return arr
end

# Vector (no valid indices to worry about)
function _truncate_by_label(
    arr::Vector,
    labels::Tuple{Vector{Set{DataType}}},
    label,
    ::Nothing
    )
    inds = map(s -> any(l -> l <: label, s), labels[1])
    return all(inds) ? arr : arr[inds]
end

# Vector (has 1-D valid indices to enforce)
function _truncate_by_label(
    arr::Vector,
    labels::Tuple{Vector{Set{DataType}}},
    label,
    valid_idxs::Vector{Bool}
    )
    new_labels = (labels[1][valid_idxs], )
    return _truncate_by_label(arr, new_labels, label, nothing)
end

# Vector (has N-D valid indices to enforce)
function _truncate_by_label(
    arr::Vector,
    labels::NTuple{N, Vector{Set{DataType}}},
    label,
    valid_idxs::Array{Bool, N}
    ) where {N}
    label_idx_array = zeros(Bool, size(valid_idxs)...)
    label_idxs = (map(s -> any(l -> l <: label, s), sets) for sets in labels)
    label_idx_array[label_idxs...] .= true
    return arr[label_idx_array[valid_idxs]]
end

# Array
function _truncate_by_label(
    arr::Array{T, N},
    labels::NTuple{N, Vector{Set{DataType}}},
    label,
    ::Nothing
    ) where {T, N}
    return arr[(map(s -> any(l -> l <: label, s), sets) for sets in labels)...]
end

# High-level
function _truncate_by_label(arr, ref, label, group_idxs, backend)
    data = backend.data
    labels = Tuple(data.support_labels[i][1:end-1] for i in group_idxs)
    valid_idxs = get(data.valid_indices, ref, nothing)
    return _truncate_by_label(arr, labels, label, valid_idxs)
end

"""
    transcription_variable(
        vref::InfiniteOpt.GeneralVariableRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
         )

Return the transcribed variable reference(s) corresponding to `vref`. Errors
if no transcription variable is found. Also can query via the syntax:
```julia
transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
    )
```
If the infinite model contains a built `TranscriptionBackend`. By default, this
method returns only transcribed variables associated with public supports. All the 
variables can be returned by setting `label = All`. 

If `vref` is infinite, then `label` will be used to search the intersection of variable 
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_variable(infvar, trans_backend)
2-element Array{VariableRef,1}:
 infvar[1]
 infvar[2]

julia> transcription_variable(hdvar, trans_backend)
hdvar

julia> transcription_variable(infvar)
2-element Array{VariableRef,1}:
 infvar[1]
 infvar[2]

julia> transcription_variable(hdvar)
hdvar
```
"""
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_variable(
        vref,
        InfiniteOpt._index_type(vref),
        backend,
        label
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
# FinVarIndex & FiniteParameterIndex
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel}
    ) where {V <: Union{FinVarIndex, InfiniteOpt.FiniteParameterIndex}}
    var = get(transcription_data(backend).finvar_mappings, vref, nothing)
    if isnothing(var)
        error("Variable reference $vref not used in transcription backend.")
    end
    return var
end

# InfVarIndex & ParameterFunctionIndex
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    label::Type{<:InfiniteOpt.AbstractSupportLabel}
    ) where {V <: Union{InfVarIndex, InfiniteOpt.ParameterFunctionIndex}}
    vars = get(transcription_data(backend).infvar_mappings, vref, nothing)
    if isnothing(vars)
        error("Variable reference $vref not used in transcription backend.")
    end
    if _ignore_label(backend, label)
        return vars
    else 
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
        return _truncate_by_label(vars, vref, label, group_idxs, backend)
    end
end

# Fallback
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type,
    backend::TranscriptionBackend,
    label
    )
    error("`transcription_variable` not defined for variables with indices of " *
          "type $(index_type) and/or is not defined for labels of type $(label).")
end

# Dispatch for internal backends
function transcription_variable(
    vref::InfiniteOpt.GeneralVariableRef; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_variable(
        vref,
        JuMP.owner_model(vref).backend,
        label = label
        )
end

"""
    InfiniteOpt.transformation_variable(
        vref::InfiniteOpt.GeneralVariableRef,
        [backend::TranscriptionBackend];
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Proper extension of [`InfiniteOpt.transformation_variable`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_variable`](@ref).
"""
function InfiniteOpt.transformation_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_variable(vref, backend, label = label)
end

"""
    InfiniteOpt.variable_supports(
        vref::InfiniteOpt.DecisionVariableRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Return the support alias mapping associated with `vref` in the transcription backend.
Errors if `vref` does not have transcripted variables.
"""
function InfiniteOpt.variable_supports(
    dvref::Union{
        InfiniteOpt.InfiniteVariableRef,
        InfiniteOpt.SemiInfiniteVariableRef, 
        InfiniteOpt.DerivativeRef
        },
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    vref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(dvref), JuMP.index(dvref))
    if !haskey(transcription_data(backend).infvar_mappings, vref)
        error("Variable reference $vref not used in transcription backend.")
    end
    supps = transcription_data(backend).infvar_supports[vref]
    if _ignore_label(backend, label)
        return supps
    else
        group_idxs = InfiniteOpt.parameter_group_int_indices(dvref)
        return _truncate_by_label(supps, vref, label, group_idxs, backend)
    end
end

# ParameterFunctionRef 
function InfiniteOpt.variable_supports(
    fref::InfiniteOpt.ParameterFunctionRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_idxs = InfiniteOpt.parameter_group_int_indices(fref)
    support_indices = support_index_iterator(backend, group_idxs)
    dims = size(support_indices)[group_idxs]
    supps = Array{Tuple, length(dims)}(undef, dims...)
    param_supps = parameter_supports(backend)
    # iterate over the indices and compute the values
    for idx in support_indices
        val_idx = idx.I[group_idxs]
        @inbounds supps[val_idx...] = Tuple(param_supps[j][idx[j]] for j in group_idxs)
    end
    # return the values
    if _ignore_label(backend, label)
        return supps
    else
        return _truncate_by_label(supps, fref, label, group_idxs, backend)
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

# InfiniteIndex & ParameterFunctionIndex
function lookup_by_support(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    support::Vector
    ) where {V <: Union{InfVarIndex, InfiniteOpt.ParameterFunctionIndex}}
    if !haskey(transcription_data(backend).infvar_lookup, vref)
        error("Variable reference $vref not used in transcription backend.")
    end
    return get(_supp_error, transcription_data(backend).infvar_lookup[vref], support)
end

# FiniteIndex
function lookup_by_support(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    backend::TranscriptionBackend,
    support::Vector
    ) where {V <: Union{FinVarIndex, InfiniteOpt.FiniteParameterIndex}}
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
    label::Type{<:InfiniteOpt.AbstractSupportLabel}
    )
    exprs = get(transcription_data(backend).measure_mappings, mref, nothing)
    if isnothing(exprs)
        error("Measure reference $mref not used in transcription backend.")
    end
    if length(exprs) > 1 && _ignore_label(backend, label)
        return exprs
    elseif length(exprs) > 1
        group_idxs = InfiniteOpt.parameter_group_int_indices(mref)
        return _truncate_by_label(exprs, mref, label, group_idxs, backend)
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
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    mref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(dmref), JuMP.index(dmref))
    if !haskey(transcription_data(backend).measure_mappings, mref)
        error("Measure reference $mref not used in transcription backend.")
    end
    supps = transcription_data(backend).measure_supports[mref]
    if length(supps) > 1 && _ignore_label(backend, label)
        return supps
    elseif length(supps) > 1
        group_idxs = InfiniteOpt.parameter_group_int_indices(mref)
        return _truncate_by_label(supps, mref, label, group_idxs, backend)
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
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Return the transcribed expression(s) corresponding to `expr`. Errors
if `expr` cannot be transcribed. Also can query via the syntax:
```julia
transcription_expression(
    expr::JuMP.AbstractJuMPScalar;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
    )
```
If the infinite model contains a built transcription backend. By default, this
method returns only transcribed expressions associated with public supports. All the 
expressions can be returned by setting `label = All`.

If `expr` is infinite, then `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_expression(my_expr, backend)
x[1] - y

julia> transcription_expression(my_expr)
x[1] - y
```
"""
function transcription_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_idxs = InfiniteOpt.parameter_group_int_indices(expr)
    support_indices = support_index_iterator(backend, group_idxs)
    dims = size(support_indices)[group_idxs]
    exprs = Array{JuMP.AbstractJuMPScalar, length(dims)}(undef, dims...)
    # iterate over the indices and compute the values
    for idx in support_indices
        supp = index_to_support(backend, idx)
        expr_idx = idx.I[group_idxs]
        @inbounds exprs[expr_idx...] = transcription_expression(expr, backend, supp)
    end
    # return the values
    if !_ignore_label(backend, label)
        exprs = _truncate_by_label(exprs, nothing, label, group_idxs, backend)
    end
    return length(support_indices) > 1 ? exprs : first(exprs)
end

# Define for variables
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel)
    return transcription_variable(vref, backend, label = label)
end

# Dispatch for internal backends
function transcription_expression(
    expr::JuMP.AbstractJuMPScalar; 
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    model = JuMP.owner_model(expr)
    isnothing(model) && return zero(JuMP.AffExpr) + JuMP.constant(expr)
    return transcription_expression(
        expr,
        model.backend,
        label = label
        )
end

"""
    InfiniteOpt.transformation_expression(
        expr::JuMP.AbstractJuMPScalar,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Proper extension of [`InfiniteOpt.transformation_expression`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_expression`](@ref).
"""
function InfiniteOpt.transformation_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_expression(expr, backend, label = label)
end

"""
    InfiniteOpt.expression_supports(
        expr::JuMP.AbstractJuMPScalar,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Return the support alias mappings associated with `expr`. Errors if `expr` cannot
be transcribed.
"""
function InfiniteOpt.expression_supports(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    # get the parameter group integer indices of the expression and form the support iterator
    group_idxs = InfiniteOpt.parameter_group_int_indices(expr)
    support_indices = support_index_iterator(backend, group_idxs)
    dims = size(support_indices)[group_idxs]
    supps = Array{Tuple, length(dims)}(undef, dims...)
    param_supps = parameter_supports(backend)
    # iterate over the indices and compute the values
    for idx in support_indices
        expr_idx = idx.I[group_idxs]
        @inbounds supps[expr_idx...] = Tuple(param_supps[j][idx[j]] for j in group_idxs)
    end
    # return the values
    if !_ignore_label(backend, label)
        supps = _truncate_by_label(supps, nothing, label, group_idxs, backend)
    end
    return length(support_indices) > 1 ? supps : first(supps)
end

################################################################################
#                             CONSTRAINT QUERIES
################################################################################
"""
    transcription_constraint(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Return the transcribed constraint reference(s) corresponding to `cref`. Errors
if `cref` has not been transcribed. Also can query via the syntax:
```julia
transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef;
    [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
    )
```
If the infinite model contains a built transcription backend. By default, this
method returns only transcribed constraints associated with public supports. All the 
constraints can be returned by setting `label = All`.

If `cref` is infinite, then `label` will be used to search the intersection of the
supports that use the label. This is defers from the default behavior which 
considers the union.

**Example**
```julia-repl
julia> transcription_constraint(fin_con, backend)
fin_con : x[1] - y <= 3.0

julia> transcription_constraint(fin_con)
fin_con : x[1] - y <= 3.0
```
"""
function transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    constr = get(transcription_data(backend).constr_mappings, cref, nothing)
    if isnothing(constr)
        error("Constraint reference $cref not used in transcription backend.")
    end
    if length(constr) > 1 && _ignore_label(backend, label)
        return constr
    elseif length(constr) > 1
        group_idxs = InfiniteOpt.parameter_group_int_indices(cref)
        return _truncate_by_label(constr, cref, label, group_idxs, backend)
    else 
        return first(constr)
    end
end

# Dispatch for internal backends
function transcription_constraint(
    cref::InfiniteOpt.InfOptConstraintRef;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_constraint(
        cref,
        JuMP.owner_model(cref).backend,
        label = label
        )
end

"""
    InfiniteOpt.transformation_constraint(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel]
        )

Proper extension of [`InfiniteOpt.transformation_constraint`](@ref) for
`TranscriptionBackend`s. This simply dispatches to [`transcription_constraint`](@ref).
"""
function InfiniteOpt.transformation_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    return transcription_constraint(cref, backend, label = label)
end

"""
    InfiniteOpt.constraint_supports(
        cref::InfiniteOpt.InfOptConstraintRef,
        backend::TranscriptionBackend;
        [label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel])

Return the support alias mappings associated with `cref`. Errors if `cref` is
not transcribed.
"""
function InfiniteOpt.constraint_supports(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::TranscriptionBackend;
    label::Type{<:InfiniteOpt.AbstractSupportLabel} = InfiniteOpt.PublicLabel
    )
    supps = get(transcription_data(backend).constr_supports, cref, nothing)
    if isnothing(supps)
        error("Constraint reference $cref not used in transcription backend.")
    end
    if length(supps) > 1 && _ignore_label(backend, label)
        return supps
    elseif length(supps) > 1
        group_idxs = InfiniteOpt.parameter_group_int_indices(cref)
        return _truncate_by_label(supps, cref, label, group_idxs, backend)
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
    )
    raw_supps = parameter_supports(backend)
    return Float64[j for i in eachindex(index.I) for j in raw_supps[i][index[i]]]
end
