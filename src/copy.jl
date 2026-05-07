################################################################################
#                              MODEL COPYING
################################################################################
"""
    InfiniteReferenceMap

Maps variable and constraint references from an original `InfiniteModel`
to their counterparts in a copied model. Returned by
[`JuMP.copy_model(::InfiniteModel)`](@ref).

`ref_map[x]` rewrites any reference into the original model so it points
at the new model. Supports `GeneralVariableRef`, `InfOptConstraintRef`,
affine / quadratic / nonlinear expressions, and any `AbstractArray` of
those. Other values pass through unchanged.
"""
struct InfiniteReferenceMap
    old_model::InfiniteModel
    new_model::InfiniteModel
end

# Generic fallback: pure-data values pass through untouched
Base.getindex(::InfiniteReferenceMap, x) = x

# Extend Base.getindex for GeneralVariableRef
function Base.getindex(m::InfiniteReferenceMap, vref::GeneralVariableRef)
    return GeneralVariableRef(m.new_model, _raw_index(vref),
                              _index_type(vref), _param_index(vref))
end

# Extend Base.getindex for InfOptConstraintRef
function Base.getindex(m::InfiniteReferenceMap, cref::InfOptConstraintRef)
    return InfOptConstraintRef(m.new_model, JuMP.index(cref))
end

# Extend Base.getindex for affine expressions
function Base.getindex(
    m::InfiniteReferenceMap,
    expr::JuMP.GenericAffExpr{T, GeneralVariableRef}
    ) where {T}
    new_expr = zero(typeof(expr))
    for (coef, var) in JuMP.linear_terms(expr)
        JuMP.add_to_expression!(new_expr, coef, m[var])
    end
    new_expr.constant = expr.constant
    return new_expr
end

# Extend Base.getindex for quadratic expressions
function Base.getindex(
    m::InfiniteReferenceMap,
    expr::JuMP.GenericQuadExpr{T, GeneralVariableRef}
    ) where {T}
    new_aff = m[expr.aff]
    new_terms = DataStructures.OrderedDict(
        JuMP.UnorderedPair(m[k.a], m[k.b]) => v for (k, v) in expr.terms
    )
    return JuMP.GenericQuadExpr(new_aff, new_terms)
end

# Extend Base.getindex for nonlinear expressions
function Base.getindex(
    m::InfiniteReferenceMap,
    expr::JuMP.GenericNonlinearExpr{V}
    ) where {V}
    return JuMP.GenericNonlinearExpr{V}(expr.head,
                                        Any[m[a] for a in expr.args])
end

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
    vt::Collections.VectorTuple,
    ref_map::InfiniteReferenceMap
    )
    vals = [ref_map[v] for v in vt.values]
    return Collections.VectorTuple{eltype(vals)}(vals, copy(vt.ranges),
                                                 copy(vt.dimensions),
                                                 copy(vt.num_columns))
end

# Rebuild a ParameterFunction (carries refs)
function _rewrite_param_function(
    pf::ParameterFunction,
    ref_map::InfiniteReferenceMap
    )
    return ParameterFunction(pf.func,
                             _rewrite_param_refs(pf.parameter_refs, ref_map),
                             copy(pf.group_int_idxs))
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
        info.binary, info.integer)
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
    return DiscreteMeasureData(ref_map[d.parameter_refs],
                               copy(d.coefficients), copy(d.supports),
                               d.label, d.weight_function,
                               copy(d.lower_bounds), copy(d.upper_bounds),
                               d.is_expect)
end

# FunctionalDiscreteMeasureData
function _rewrite_measure_data(
    d::FunctionalDiscreteMeasureData,
    ref_map::InfiniteReferenceMap
    )
    return FunctionalDiscreteMeasureData(ref_map[d.parameter_refs],
                                         d.coeff_function,
                                         d.min_num_supports, d.label,
                                         d.generative_supp_info,
                                         d.weight_function,
                                         copy(d.lower_bounds),
                                         copy(d.upper_bounds),
                                         d.is_expect)
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

################################################################################
#                              COPY MODEL
################################################################################
"""
    JuMP.copy_model(model::InfiniteModel)

Return a copy of `model` and an [`InfiniteReferenceMap`](@ref) that
can be used to obtain the variable and constraint references of the
new model corresponding to references from the original model.

The copied model has a fresh `TranscriptionBackend` and
`ready_to_optimize = false`. An optimizer must be set on the new
model before calling `optimize!`.

The implementation walks `model`'s registries and reconstructs each
entry inline, using `ref_map[...]` for ref-bearing pieces and
`Base.copy(...)` for shallow clones of the data wrappers. It does
not call `Base.deepcopy` on the model.

!!! note
    Unlike `JuMP.copy_model(::JuMP.GenericModel)`, this method does
    not currently accept a `filter_constraints` keyword argument —
    all constraints are copied.

**Example**
```julia-repl
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 1], num_supports = 10);

julia> @variable(model, x, Infinite(t));

julia> @constraint(model, c, x <= 1);

julia> new_model, ref_map = copy_model(model);

julia> new_x = ref_map[x];

julia> new_c = ref_map[c];
```
"""
function JuMP.copy_model(model::InfiniteModel)
    new_model = InfiniteModel()
    ref_map = InfiniteReferenceMap(model, new_model)

    # Independent parameters
    for (idx, data) in model.independent_params
        new_data = copy(data)
        new_data.parameter = copy(data.parameter)
        new_model.independent_params[idx] = new_data
    end
    # Dependent parameters
    for (idx, data) in model.dependent_params
        new_data = copy(data)
        new_data.parameters = copy(data.parameters)
        new_model.dependent_params[idx] = new_data
    end
    # Finite parameters (FiniteParameter is just Float64; aliasing is fine)
    for (idx, data) in model.finite_params
        new_model.finite_params[idx] = copy(data)
    end
    new_model.param_group_indices = copy(model.param_group_indices)

    # Parameter functions
    for (idx, data) in model.param_functions
        new_data = copy(data)
        new_data.func = _rewrite_param_function(data.func, ref_map)
        new_model.param_functions[idx] = new_data
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
            copy(var.group_int_idxs))
        new_model.infinite_vars[idx] = new_data
    end

    # Semi-infinite variables
    for (idx, data) in model.semi_infinite_vars
        new_data = copy(data)
        var = data.variable
        new_data.variable = SemiInfiniteVariable(
            var.info, ref_map[var.infinite_variable_ref],
            copy(var.eval_support), copy(var.group_int_idxs))
        new_model.semi_infinite_vars[idx] = new_data
    end
    new_model.semi_lookup = _rewrite_var_lookup(model.semi_lookup, ref_map)

    # Point variables
    for (idx, data) in model.point_vars
        new_data = copy(data)
        var = data.variable
        new_data.variable = PointVariable(
            var.info, ref_map[var.infinite_variable_ref],
            copy(var.parameter_values))
        new_model.point_vars[idx] = new_data
    end
    new_model.point_lookup = _rewrite_var_lookup(model.point_lookup, ref_map)

    # Finite variables (ScalarVariable holds only Float64 bounds; no refs)
    for (idx, data) in model.finite_vars
        new_model.finite_vars[idx] = copy(data)
    end

    # Derivatives
    for (idx, data) in model.derivatives
        new_data = copy(data)
        var = data.variable
        new_data.variable = Derivative(
            _rewrite_info(var.info, ref_map),
            ref_map[var.variable_ref],
            ref_map[var.parameter_ref], var.order)
        new_model.derivatives[idx] = new_data
    end
    new_model.deriv_lookup = Dict(
        (ref_map[k[1]], ref_map[k[2]], k[3]) => v
        for (k, v) in model.deriv_lookup)

    # Measures
    for (idx, data) in model.measures
        new_data = copy(data)
        meas = data.measure
        new_data.measure = Measure(
            ref_map[meas.func],
            _rewrite_measure_data(meas.data, ref_map),
            copy(meas.group_int_idxs), meas.constant_func)
        new_model.measures[idx] = new_data
    end

    # Constraints
    for (idx, data) in model.constraints
        new_data = copy(data)
        new_data.constraint = _rewrite_constraint(data.constraint, ref_map)
        new_model.constraints[idx] = new_data
    end

    # Objective
    new_model.objective_sense = model.objective_sense
    new_model.objective_function = ref_map[model.objective_function]
    new_model.objective_has_measures = model.objective_has_measures

    # NLP user operators (pure data, no model refs)
    new_model.operators = copy(model.operators)
    new_model.op_lookup = copy(model.op_lookup)

    # obj_dict — values whose types lack a getindex method on
    # InfiniteReferenceMap pass through unchanged via the catch-all fallback
    new_model.obj_dict = Dict(k => ref_map[v] for (k, v) in model.obj_dict)

    # Backend / ready_to_optimize already set fresh by InfiniteModel()
    new_model.optimize_hook = model.optimize_hook

    # Extension data — JuMP's documented protocol; errors if an extension
    # stores data without defining copy_extension_data
    for (key, data) in model.ext
        new_model.ext[key] = JuMP.copy_extension_data(data, new_model, model)
    end

    return new_model, ref_map
end
