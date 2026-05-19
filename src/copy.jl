################################################################################
#                          INFINITE REFERENCE MAP DISPATCH
################################################################################
# `InfiniteReferenceMap` is defined in `datatypes.jl`. Dispatched
# `Base.getindex` methods for individual ref / expression types live in
# their owning files (`general_variables.jl`, `constraints.jl`,
# `expressions.jl`). Only the catch-all and array fallbacks live here.

# Generic fallback: pure-data values pass through untouched
Base.getindex(::InfiniteReferenceMap, x) = x

# Extend Base.getindex for arrays via broadcast
function Base.getindex(m::InfiniteReferenceMap, val::AbstractArray)
    return getindex.(m, val)
end

Base.broadcastable(m::InfiniteReferenceMap) = Ref(m)

################################################################################
#                              MODEL COPY HELPERS
################################################################################
# Internal helpers used by `JuMP.copy_model(::InfiniteModel)` to rebuild
# values that carry references into the source model. Each returns a
# fresh value with refs rewritten via `ref_map`.

# Rebuild a VectorTuple of refs
function _rewrite_param_refs(
    vt::Collections.VectorTuple{GeneralVariableRef},
    ref_map::InfiniteReferenceMap
    )
    vals = [ref_map[v] for v in vt.values]
    return Collections.VectorTuple{eltype(vals)}(
        vals, copy(vt.ranges),
        copy(vt.dimensions),
        copy(vt.num_columns)
        )
end

# Rebuild a ParameterFunction (carries refs)
function _rewrite_param_function(
    pf::ParameterFunction,
    ref_map::InfiniteReferenceMap
    )
    return ParameterFunction(
        pf.func,
        _rewrite_param_refs(pf.parameter_refs, ref_map),
        copy(pf.group_int_idxs)
        )
end

## Rebuild a JuMP.VariableInfo bound — ParameterFunction gets refs
## rewritten, anything else (Float64, etc.) passes through

# ParameterFunction
function _rewrite_bound(b::ParameterFunction, ref_map::InfiniteReferenceMap)
    return _rewrite_param_function(b, ref_map)
end

# Generic fallback
_rewrite_bound(b, ::InfiniteReferenceMap) = b

# Rebuild a JuMP.VariableInfo (bounds may be ParameterFunctions)
function _rewrite_info(info::JuMP.VariableInfo, ref_map::InfiniteReferenceMap)
    return JuMP.VariableInfo(
        info.has_lb,    _rewrite_bound(info.lower_bound, ref_map),
        info.has_ub,    _rewrite_bound(info.upper_bound, ref_map),
        info.has_fix,   _rewrite_bound(info.fixed_value, ref_map),
        info.has_start, _rewrite_bound(info.start, ref_map),
        info.binary, info.integer
        )
end

# Rebuild a (ref, Vector{Float64})-keyed lookup dict
function _rewrite_var_lookup(d::Dict, ref_map::InfiniteReferenceMap)
    return Dict((ref_map[k[1]], copy(k[2])) => v for (k, v) in d)
end

## Rebuild a measure's data field

# DiscreteMeasureData
function _rewrite_measure_data(
    d::DiscreteMeasureData,
    ref_map::InfiniteReferenceMap
    )
    return DiscreteMeasureData(
        ref_map[d.parameter_refs],
        copy(d.coefficients), copy(d.supports),
        d.label, d.weight_function,
        copy(d.lower_bounds), copy(d.upper_bounds),
        d.is_expect
        )
end

# FunctionalDiscreteMeasureData
function _rewrite_measure_data(
    d::FunctionalDiscreteMeasureData,
    ref_map::InfiniteReferenceMap
    )
    return FunctionalDiscreteMeasureData(
        ref_map[d.parameter_refs],
        d.coeff_function,
        d.min_num_supports, d.label,
        d.generative_supp_info,
        d.weight_function,
        copy(d.lower_bounds),
        copy(d.upper_bounds),
        d.is_expect
        )
end

# Generic fallback
_rewrite_measure_data(d, ::InfiniteReferenceMap) = d

## Rebuild a JuMP constraint

# ScalarConstraint
function _rewrite_constraint(
    c::JuMP.ScalarConstraint,
    ref_map::InfiniteReferenceMap
    )
    return JuMP.ScalarConstraint(ref_map[c.func], c.set)
end

# VectorConstraint
function _rewrite_constraint(
    c::JuMP.VectorConstraint,
    ref_map::InfiniteReferenceMap
    )
    return JuMP.VectorConstraint(ref_map[c.func], c.set, c.shape)
end

# Generic fallback
_rewrite_constraint(c, ::InfiniteReferenceMap) = c

## Decide whether an `obj_dict` entry should survive the copy. An
## entry is dropped if it (or, for a container, any element of it)
## points at a variable that was deleted in the source or a constraint
## that was not copied (deleted in the source, or removed by
## `filter_constraints`). `source_to_new` is populated only for
## constraints actually copied, so checking membership there covers
## both cases.

# Scalar variable ref: drop if invalid in source
function _copyable_obj_dict_entry(
    model::InfiniteModel,
    ::InfiniteReferenceMap,
    v::GeneralVariableRef
    )
    return JuMP.is_valid(model, v)
end

# Scalar constraint ref: drop if its source index was not copied
function _copyable_obj_dict_entry(
    ::InfiniteModel,
    ref_map::InfiniteReferenceMap,
    v::InfOptConstraintRef
    )
    return haskey(ref_map.source_to_new, JuMP.index(v))
end

# Container of refs: drop the whole name if any element is uncopyable
function _copyable_obj_dict_entry(
    model::InfiniteModel,
    ref_map::InfiniteReferenceMap,
    v::AbstractArray
    )
    return all(
        _copyable_obj_dict_entry(model, ref_map, el) for el in v
        )
end

# Generic fallback: pure-data / unhandled types pass through
_copyable_obj_dict_entry(
    ::InfiniteModel, ::InfiniteReferenceMap, ::Any
    ) = true

################################################################################
#                              COPY MODEL
################################################################################
"""
    JuMP.copy_model(
        model::InfiniteModel;
        filter_constraints::Union{Nothing, Function} = nothing
        )

Return a copy of `model` and an [`InfiniteReferenceMap`](@ref) that
can be used to obtain the variable and constraint references of the
new model corresponding to references from the original model.

The copied model carries an empty backend of the same type as the
source's (built via [`copy_empty_backend`](@ref)), preserving the
source's solver and configuration but none of its transformation
data. `ready_to_optimize` is `false`.

The implementation walks `model`'s registries, rebuilds each entry
inline with `ref_map[...]` for ref-bearing pieces, and inserts the
result into `new_model` via the standard `_add_data_object` dispatch
(which calls `MOIUC.add_item`). The new index returned by each
insertion is recorded in `ref_map.source_to_new`, an explicit
`ObjectIndex => ObjectIndex` translation table. `InfiniteReferenceMap`
resolves copied refs by looking up the new index in that table, so the
copy remains valid even when the source has `CleverDict` gaps from
prior `delete!` calls (the new model's indices are dense and may
differ from the source's).

Entries left in the source's object dictionary that point at deleted
refs are skipped on the copy, so dangling registrations do not leak
into `new_model.obj_dict`.

## Keyword Arguments

- `filter_constraints::Union{Nothing, Function} = nothing`: when not
  `nothing`, a predicate that receives each source `InfOptConstraintRef`
  and returns `true` to copy it or `false` to skip it. Skipped
  constraints are absent from `new_model.constraints` and from
  `ref_map.source_to_new`, so `ref_map[dropped]` raises a `KeyError`
  (the same convention `JuMP.copy_model(::GenericModel)` uses for
  filtered constraints). Variables are not affected — a constraint's
  variables remain in `new_model` whether or not the constraint
  itself was copied. If a constraint registered by name in
  `model.obj_dict` is filtered out, its name registration is dropped
  from `new_model.obj_dict`; for a constraint container registered
  under a single name (e.g. `@constraint(model, c[1:n], ...)`), the
  whole name is dropped if **any** of its elements is filtered out.

**Example**
```julia-repl
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 1], num_supports = 10);

julia> @variable(model, x, Infinite(t));

julia> @constraint(model, c, x <= 1);

julia> new_model, ref_map = copy_model(model);

julia> new_x = ref_map[x];

julia> new_c = ref_map[c];

julia> # Skip `c` on the copy:
       new_model2, ref_map2 = copy_model(
           model;
           filter_constraints = cref -> cref != c
           );

julia> haskey(new_model2.obj_dict, :c)
false
```
"""
function JuMP.copy_model(
    model::InfiniteModel;
    filter_constraints::Union{Nothing, Function} = nothing
    )
    new_model = InfiniteModel(copy_empty_backend(model.backend))
    ref_map = InfiniteReferenceMap(model, new_model)

    # Parameters: iterate in original creation order so that the
    # `param_group_indices` vector is rebuilt by `_add_data_object`'s
    # `push!` in the same order as the source model.
    for grp_idx in model.param_group_indices
        if grp_idx isa IndependentParameterIndex
            data = model.independent_params[grp_idx]
            new_data = copy(data)
            new_data.parameter = copy(data.parameter)
        else
            data = model.dependent_params[grp_idx]
            new_data = copy(data)
            new_data.parameters = copy(data.parameters)
        end
        ref_map.source_to_new[grp_idx] =
            _add_data_object(new_model, new_data)
    end
    # Finite parameters (FiniteParameter is just Float64; aliasing is fine)
    for (idx, data) in model.finite_params
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, copy(data))
    end

    # Parameter functions
    for (idx, data) in model.param_functions
        new_data = copy(data)
        new_data.func = _rewrite_param_function(data.func, ref_map)
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end

    # piecewise_vars: Dict{Idx, Set{Idx}}, all immutable; clone the Sets
    new_model.piecewise_vars =
        Dict(k => copy(v) for (k, v) in model.piecewise_vars)

    # Infinite variables
    for (idx, data) in model.infinite_vars
        new_data = copy(data)
        var = data.variable
        new_data.variable = InfiniteVariable(
            _rewrite_info(var.info, ref_map),
            _rewrite_param_refs(var.parameter_refs, ref_map),
            copy(var.group_int_idxs)
            )
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end

    # Semi-infinite variables
    for (idx, data) in model.semi_infinite_vars
        new_data = copy(data)
        var = data.variable
        new_data.variable = SemiInfiniteVariable(
            var.info, ref_map[var.infinite_variable_ref],
            copy(var.eval_support), copy(var.group_int_idxs)
            )
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end
    new_model.semi_lookup = _rewrite_var_lookup(model.semi_lookup, ref_map)

    # Point variables
    for (idx, data) in model.point_vars
        new_data = copy(data)
        var = data.variable
        new_data.variable = PointVariable(
            var.info, ref_map[var.infinite_variable_ref],
            copy(var.parameter_values)
            )
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end
    new_model.point_lookup = _rewrite_var_lookup(model.point_lookup, ref_map)

    # Finite variables (ScalarVariable holds only Float64 bounds; no refs)
    for (idx, data) in model.finite_vars
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, copy(data))
    end

    # Derivatives
    for (idx, data) in model.derivatives
        new_data = copy(data)
        var = data.variable
        new_data.variable = Derivative(
            _rewrite_info(var.info, ref_map),
            ref_map[var.variable_ref],
            ref_map[var.parameter_ref], var.order
            )
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end
    new_model.deriv_lookup = Dict(
        (ref_map[k[1]], ref_map[k[2]], k[3]) => v
        for (k, v) in model.deriv_lookup
        )

    # Measures
    for (idx, data) in model.measures
        new_data = copy(data)
        meas = data.measure
        new_data.measure = Measure(
            ref_map[meas.func],
            _rewrite_measure_data(meas.data, ref_map),
            copy(meas.group_int_idxs), meas.constant_func
            )
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end

    # Constraints — honors `filter_constraints` if provided. Skipped
    # constraints are never inserted and never recorded in
    # `source_to_new`, so `ref_map[dropped]` raises `KeyError`.
    for (idx, data) in model.constraints
        if filter_constraints !== nothing &&
           !filter_constraints(InfOptConstraintRef(model, idx))
            continue
        end
        new_data = copy(data)
        new_data.constraint = _rewrite_constraint(data.constraint, ref_map)
        ref_map.source_to_new[idx] =
            _add_data_object(new_model, new_data)
    end

    # Objective
    new_model.objective_sense = model.objective_sense
    new_model.objective_function = ref_map[model.objective_function]
    new_model.objective_has_measures = model.objective_has_measures

    # NLP user operators (pure data, no model refs)
    new_model.operators = copy(model.operators)
    new_model.op_lookup = copy(model.op_lookup)

    # obj_dict — values whose types lack a getindex method on
    # InfiniteReferenceMap pass through unchanged via the catch-all
    # fallback. Entries that point at deleted variables or at
    # constraints that were not copied (deleted in source or filtered
    # out) are skipped, so the copy never carries dangling references.
    # For container registrations, the whole name is dropped if any
    # element is uncopyable (see `_copyable_obj_dict_entry`).
    new_model.obj_dict = Dict(
        k => ref_map[v]
        for (k, v) in model.obj_dict
        if _copyable_obj_dict_entry(model, ref_map, v)
        )

    # Backend / ready_to_optimize already set fresh by InfiniteModel()
    new_model.optimize_hook = model.optimize_hook

    # Extension data — JuMP's documented protocol; errors if an extension
    # stores data without defining copy_extension_data
    for (key, data) in model.ext
        new_model.ext[key] = JuMP.copy_extension_data(data, new_model, model)
    end

    return new_model, ref_map
end
