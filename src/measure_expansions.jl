"""
    make_point_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                            ivref::InfiniteVariableRef,
                            support::Union{Tuple, VectorTuple{Float64}}
                            )::PointVariableRef

Make a point variable for infinite variable `ivref` at `support`, add it to
the `write_model`, and return the `PointVariableRef`. This is an internal method
for point variables produced by expanding measures via [`expand_measure`](@ref).
This is also useful for those writing extension optimizer models and wish to
expand measures without modifiying the `InfiniteModel`. In such cases, `write_model`
should be the optimizer model and [`add_measure_variable`](@ref add_measure_variable(::JuMP.Model, ::Any, ::Any))
should be extended appropriately for point variables. Errors if `write_model` is
an optimizer model and `add_measure_variable` is not properly extended.
"""
function make_point_variable_ref(write_model::InfiniteModel,
                                 ivref::InfiniteVariableRef,
                                 support::VectorTuple{Float64})::PointVariableRef
    var = PointVariable(_variable_info(ivref), ivref, support)
    return JuMP.add_variable(write_model, var)
end

"""
    add_measure_variable(model::JuMP.Model, var,
                         key::Val{:ext_key_name})::GeneralVariableRef

Add a measure variable `var` to the optimizer model `model` (with `key`) and
return the correct `InfiniteOpt` variable reference. This is an internal method
used by [`make_point_variable_ref`](@ref) and [`make_reduced_variable_ref`](@ref)
to make point variables and reduced variables when the `write_model` is an
optimizer model. This is useful for extensions that wish to expand measures, but
without changing the original `InfiniteModel`. Thus, this should be extended for
adding `PointVariable`s and `ReducedInfiniteVariable`s for such extensions.
Otherwise, an error is thrown for unextended variable and/or optimizer model types.
"""
function add_measure_variable(model::JuMP.Model, var, key)
    error("`add_measure_variable` not defined for variable of type `$(typeof(var))` " *
          "and an optimizer model with key `$key`.")
end

# Store/add the variable to the optimizer model via add_measure_variable
# This avoids changing the InfiniteModel unexpectedly
function make_point_variable_ref(write_model::JuMP.Model, # this should be an optimizer model
                                 ivref::InfiniteVariableRef,
                                 support::VectorTuple{Float64})::PointVariableRef
    var = PointVariable(_variable_info(ivref), ivref, support)
    opt_key = optimizer_model_key(write_model)
    return add_measure_variable(write_model, var, Val(opt_key))
end

# Tuple support input and dispatch to appropriate VectorTuple function
function make_point_variable_ref(write_model::Union{JuMP.Model, InfiniteModel},
                                 ivref::InfiniteVariableRef,
                                 support::Tuple)::PointVariableRef
    vt_support = VectorTuple{Float64}(support)
    prefs = raw_parameter_refs(ivref)
    for i in eachindex(prefs.indices)
        vt_support.indices[i] = prefs.indices[i]
    end
    return make_point_variable_ref(write_model, ivref, vt_support)
end


"""
    make_reduced_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                              ivref::InfiniteVariableRef,
                              indices::Vector{Int},
                              values::Vector{Float64}
                              )::ReducedInfiniteVariableRef

Make a reduced variable for infinite variable `ivref` at `support`, add it to
the `write_model`, and return the `ReducedInfiniteVariableRef`. This is an internal method
for reduced variables produced by expanding measures via [`expand_measure`](@ref).
This is also useful for those writing extension optimizer models and wish to
expand measures without modifiying the `InfiniteModel`. In such cases, `write_model`
should be the optimizer model and [`add_measure_variable`](@ref add_measure_variable(::JuMP.Model, ::Any, ::Any))
should be extended appropriately for reduced variables. Errors if `write_model`
is an optimizer model and `add_measure_variable` is not properly extended.
"""
function make_reduced_variable_ref(write_model::InfiniteModel,
                                   ivref::InfiniteVariableRef,
                                   indices::Vector{Int},
                                   values::Vector{Float64}
                                   )::ReducedInfiniteVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    index = write_model.next_reduced_index += 1
    write_model.reduced_info[index] = ReducedInfiniteInfo(ivref, eval_supps)
    if haskey(write_model.infinite_to_reduced, JuMP.index(ivref))
        push!(write_model.infinite_to_reduced[JuMP.index(ivref)], index)
    else
        write_model.infinite_to_reduced[JuMP.index(ivref)] = [index]
    end
    return ReducedInfiniteVariableRef(write_model, index)
end

# Add reduced infinite variables in the optimizer model without modifying the InfiniteModel
function make_reduced_variable_ref(write_model::JuMP.Model,
                                   ivref::InfiniteVariableRef,
                                   indices::Vector{Int},
                                   values::Vector{Float64}
                                   )::ReducedInfiniteVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    var = ReducedInfiniteInfo(ivref, eval_supps)
    key = optimizer_model_key(write_model)
    return add_measure_variable(write_model, var, Val(key))
end

"""
    delete_internal_reduced_variable(write_model::Union{InfiniteModel, JuMP.Model},
                                     rvref::ReducedInfiniteVariableRef)

Delete the variable associated with `rvref` from `write_model` if it is purely
an internal variable only used for measure expansion and is no longer needed.
For `write_model`s that are an optimizer model,
[`delete_reduced_variable`](@ref delete_reduced_variable(::JuMP.Model, ::Any,::Any))
will need to be extended for this this to work. Otherwise, a warning will be thrown.
Note that this is intended as an internal method to assist with extensions to
[`expand_measure`](@ref).
"""
function delete_internal_reduced_variable(write_model::InfiniteModel,
                                          rvref::ReducedInfiniteVariableRef)
    if !used_by_measure(rvref) && !used_by_constraint(rvref)
        JuMP.delete(write_model, rvref)
    end
    return
end

"""
    delete_reduced_variable(model::JuMP.Model, vref, key::Val{:ext_key_name})

Delete the reduced variable associated with `vref` from the optimizer model
`model` with associated extension key `:ext_key_name`. A warning is thrown if this
is not properly extended. This is intended as a helper function for
[`delete_internal_reduced_variable`](@ref) which is used by [`expand_measure`](@ref).
"""
function delete_reduced_variable(write_model::JuMP.Model, vref, key)
    @warn "'delete_reduced_variable' not extended for reduced variable type " *
          "`$(typeof(vref))` and optimizer model with key `$(typeof(key).parameters[1])`."
    return
end

# Delete reduced infinite variable from optimizer model if it was not made by the InfiniteModel
function delete_internal_reduced_variable(write_model::JuMP.Model,
                                          rvref::ReducedInfiniteVariableRef)
    if !JuMP.is_valid(JuMP.owner_model(rvref), rvref)
        key = optimizer_model_key(write_model)
        delete_reduced_variable(write_model, rvref, Val(key))
    end
    return
end

"""
    expand_measure(expr, data::AbstractMeasureData,
                   write_model::JuMP.AbstractModel)::JuMP.AbstractJuMPScalar

Return the finite reformulation of a measure containing a variable/parameter
expression `expr` with measure data `data`. Here `write_model` is the target
model where this expanded expression will be used. Thus, any variables that need
to be created will be added to `write_model`. The methods [`make_point_variable_ref`](@ref)
and [`make_reduced_variable_ref`](@ref) should be used as appropriate to create
these variables. Developers might also choose to use
[`delete_internal_reduced_variable`](@ref) in order to remove reduced variables
once they are no longer needed. Note this is intended as an internal function,
but will need to be extended for unsupported `expr` types and for user-defined
measure data types. Principally, this is leveraged to enable the user methods
[`expand`](@ref) and [`expand_all_measures!`](@ref).
"""
function expand_measure end

# InfiniteVariableRef (DiscreteMeasureData)
function expand_measure(ivref::InfiniteVariableRef,
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    var_prefs = parameter_list(ivref)
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat variable as constant if doesn't have measure parameter
    if !(pref in var_prefs)
        var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, typeof(ivref)}(0, ivref => var_coef)
    # make point variables if var_prefs = pref (it is the only dependence)
    elseif length(var_prefs) == 1
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_point_variable_ref(write_model, ivref, (supps[i],))
                    for i in eachindex(coeffs)))
    # make reduced variables if the variable contains other parameters
    else
        index = findfirst(isequal(pref), var_prefs)
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_reduced_variable_ref(write_model, ivref, [index], [supps[i]])
                    for i in eachindex(coeffs)))
    end
end

# InfiniteVariableRef (MultiDiscreteMeasureData)
function expand_measure(ivref::InfiniteVariableRef,
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    var_prefs = raw_parameter_refs(ivref)
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # figure out the parameter groups
    group = group_id(first(prefs))
    groups = group_id.(var_prefs[:, 1])
    group_index = findfirst(isequal(group), groups)
    # treat variable as constant if doesn't have measure parameter
    if group_index isa Nothing
        var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, typeof(ivref)}(0, ivref => var_coef)
    # make point variables if all var_prefs are contained in prefs
    elseif size(var_prefs, 1) == 1 && length(var_prefs) == length(prefs)
        # var_prefs = prefs
        if parameter_list(ivref) == prefs
            return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                        make_point_variable_ref(write_model, ivref, (supps[:, i],))
                        for i in eachindex(coeffs)))
        # assume that var_prefs are prefs in a different order (otherwise will error)
        else
            indices = [findfirst(isequal(pref), prefs) for pref in var_prefs]
            indices isa Vector{<:Int} || error("Invalid high-dimensional measure, parameters partially " *
                                               "overlap with the multi-dimensional parameters in $ivref.")
            new_supps = supps[indices, :]
            return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                        make_point_variable_ref(write_model, ivref, (new_supps[:, i],))
                        for i in eachindex(coeffs)))
        end
    # make reduced variables if the variable contains other parameters
    elseif length(var_prefs[group_index, :]) == length(prefs)
        indices = [findfirst(isequal(pref), var_prefs) for pref in prefs]
        indices isa Vector{<:Int} || error("Invalid high-dimensional measure, parameters partially " *
                                           "overlap with the multi-dimensional parameters in $ivref.")
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, supps[:, i])
                    for i in eachindex(coeffs)))
    else
        error("Invalid high-dimensional measure, the dimensions of the parameters " *
              "do not agree with the related multi-dimensional parameters in $ivref. " *
              "Consider using nested 1-dimensional measures if this behavior is wanted.")
    end
end

# Write point support given reduced info and the support
function _make_point_support(orig_prefs::VectorTuple{ParameterRef},
                             support_dict::Dict{Int, Float64},
                             index::Int, value::Float64)::VectorTuple{Float64}
    values = [i == index ? value : support_dict[i] for i in eachindex(orig_prefs)]
    return VectorTuple{Float64}(values, orig_prefs.ranges, orig_prefs.indices)
end

# ReducedInfiniteVariableRef (DiscreteMeasureData)
function expand_measure(rvref::ReducedInfiniteVariableRef,
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    ivref = infinite_variable_ref(rvref)
    orig_prefs = raw_parameter_refs(ivref)
    var_prefs = parameter_list(rvref)
    eval_supps = eval_supports(rvref)
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat variable as constant if doesn't have measure parameter
    if !(pref in var_prefs)
        var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        expr = JuMP.GenericAffExpr{Float64, typeof(rvref)}(0, rvref => var_coef)
    # make point variables if var_prefs = pref (it is the only dependence)
    elseif length(var_prefs) == 1
        index = findfirst(isequal(pref), orig_prefs)
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_point_variable_ref(write_model, ivref,
                    _make_point_support(orig_prefs, eval_supps, index, supps[i]))
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, rvref)
    # make reduced variables if the variable contains other parameters
    else
        index = findfirst(isequal(pref), orig_prefs)
        indices = [collect(keys(eval_supps)); index]
        vals = collect(values(eval_supps))
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, Float64[vals; supps[i]])
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, rvref)
    end
    return expr
end

# Write point support given reduced info and the support
function _make_point_support(orig_prefs::VectorTuple{ParameterRef},
                             support_dict::Dict{Int, Float64},
                             group_index::Int,
                             values::Vector{Float64})::VectorTuple{Float64}
    offset = orig_prefs.ranges[group_index].start - 1
    new_values = [haskey(support_dict, i) ? support_dict[i] : values[i - offset]
                  for i in eachindex(orig_prefs)]
    return VectorTuple{Float64}(new_values, orig_prefs.ranges, orig_prefs.indices)
end

# ReducedInfiniteVariableRef (MultiDiscreteMeasureData)
function expand_measure(rvref::ReducedInfiniteVariableRef,
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    ivref = infinite_variable_ref(rvref)
    orig_prefs = raw_parameter_refs(ivref)
    var_prefs = raw_parameter_refs(rvref)
    eval_supps = eval_supports(rvref)
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # figure out the parameter groups
    group = group_id(first(prefs))
    groups = group_id.(var_prefs[:, 1])
    group_index = findfirst(isequal(group), groups)
    # treat variable as constant if doesn't have measure parameter
    if group_index isa Nothing
        var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        expr = JuMP.GenericAffExpr{Float64, typeof(rvref)}(0, rvref => var_coef)
    # make point variables if var_prefs = pref (it is the only dependence)
    elseif size(var_prefs, 1) == 1 && length(var_prefs) == length(prefs)
        orig_group_index = findfirst(isequal(group), group_id.(orig_prefs[:, 1]))
        # var_prefs = prefs
        if parameter_list(var_prefs) == prefs
            expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                        make_point_variable_ref(write_model, ivref,
                        _make_point_support(orig_prefs, eval_supps, orig_group_index, supps[:, i]))
                        for i in eachindex(coeffs)))
        # assume that var_prefs are prefs in a different order (otherwise will error)
        else
            indices = [findfirst(isequal(pref), prefs) for pref in var_prefs]
            indices isa Vector{<:Int} || error("Invalid high-dimensional measure, parameters partially " *
                                               "overlap with the multi-dimensional parameters in $rvref.")
            new_supps = supps[indices, :]
            expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                        make_point_variable_ref(write_model, ivref,
                        _make_point_support(orig_prefs, eval_supps, orig_group_index, new_supps[:, i]))
                        for i in eachindex(coeffs)))
        end
        delete_internal_reduced_variable(write_model, rvref)
    # make reduced variables if the variable contains other parameters
    elseif length(var_prefs[group_index, :]) == length(prefs)
        new_indices = [findfirst(isequal(pref), orig_prefs) for pref in prefs]
        check = any(new_indices[i] in keys(eval_supps) for i in eachindex(new_indices))
        if !(new_indices isa Vector{<:Int}) || check
            error("Invalid high-dimensional measure, parameters partially " *
                  "overlap with the multi-dimensional parameters in $rvref.")
        end
        indices = [collect(keys(eval_supps)); new_indices]
        vals = collect(values(eval_supps))
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, Float64[vals; supps[:, i]])
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, rvref)
    else
        error("Invalid high-dimensional measure, the dimensions of the parameters " *
              "do not agree with the related multi-dimensional parameters in $rvref. " *
              "Consider using nested 1-dimensional measures if this behavior is wanted.")
    end
    return expr
end

# FiniteVariableRef (DiscreteMeasureData)
function expand_measure(vref::FiniteVariableRef,
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the variable as a constant and build the expression
    var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
    return JuMP.GenericAffExpr{Float64, typeof(vref)}(0, vref => var_coef)
end

# FiniteVariableRef (MultiDiscreteMeasureData)
function expand_measure(vref::FiniteVariableRef,
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the variable as a constant and build the expression
    var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
    return JuMP.GenericAffExpr{Float64, typeof(vref)}(0, vref => var_coef)
end

# ParameterRef (DiscreteMeasureData)
function expand_measure(pref::ParameterRef,
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.GenericAffExpr, Float64}
    # pull in the needed information
    meas_pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the parameter
    if meas_pref != pref
        par_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, ParameterRef}(0, pref => par_coef)
    # replace the parameter with its value
    else
        return sum(coeffs[i] * w(supps[i]) * supps[i] for i in eachindex(coeffs))
    end
end

# ParameterRef (MultiDiscreteMeasureData)
function expand_measure(pref::ParameterRef,
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.GenericAffExpr, Float64}
    # pull in the needed information
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # find the position of pref if it is in the data
    index = findfirst(isequal(pref), prefs)
    # treat the parameter
    if index isa Nothing
        par_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, ParameterRef}(0, pref => par_coef)
    # replace the parameter with its value
    else
        return sum(coeffs[i] * w(supps[:, i]) * supps[index, i] for i in eachindex(coeffs))
    end
end

# Generic[Aff/Quad]Expr-FiniteVariableRef (DiscreteMeasureData)
function expand_measure(expr::JuMP._GenericAffOrQuadExpr{C, <:FiniteVariableRef},
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.AbstractJuMPScalar where {C}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat expr as a constant since it is finite
    expr_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
    return JuMP.@expression(write_model, expr_coef * expr)
end

# Generic[Aff/Quad]Expr-FiniteVariableRef (MultiDiscreteMeasureData)
function expand_measure(expr::JuMP._GenericAffOrQuadExpr{C, <:FiniteVariableRef},
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::JuMP.AbstractJuMPScalar where {C}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat expr as a constant since it is finite
    expr_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
    return JuMP.@expression(write_model, expr_coef * expr)
end

# GenericAffExpr-GeneralVariableRef (DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # expand each variable independently and add all together
    constant_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
    return JuMP.@expression(write_model, sum(coef * expand_measure(var, data, write_model)
                for (var, coef) in expr.terms) + expr.constant * constant_coef)
end

# GenericAffExpr-GeneralVariableRef (MultiDiscreteMeasureData)
function expand_measure(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # expand each variable independently and add all together
    constant_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
    return JuMP.@expression(write_model, sum(coef * expand_measure(var, data, write_model)
                for (var, coef) in expr.terms) + expr.constant * constant_coef)
end

# GenericQuadExpr-GeneralVariableRef (DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericQuadExpr{C, <:GeneralVariableRef},
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # get needed info
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    name = measure_name(data)
    # expand the affine expression
    new_expr = expand_measure(expr.aff, data, write_model)
    # expand the quadratic terms
    for i in eachindex(coeffs)
        # make viable data objects so we can multiply the terms
        coef_data = DiscreteMeasureData(pref, [coeffs[i]], [supps[i]], name, w)
        simple_data = DiscreteMeasureData(pref, [1], [supps[i]], name,
                                          default_weight)
        new_expr = JuMP.@expression(write_model, new_expr + sum(coef *
                        expand_measure(pair.a, coef_data, write_model) *
                        expand_measure(pair.b, simple_data, write_model)
                        for (pair, coef) in expr.terms))
    end
    return new_expr
end

# GenericQuadExpr-GeneralVariableRef (MultiDiscreteMeasureData)
function expand_measure(expr::JuMP.GenericQuadExpr{C, <:GeneralVariableRef},
                        data::MultiDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # get needed info
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    name = measure_name(data)
    # expand the affine expression
    new_expr = expand_measure(expr.aff, data, write_model)
    # expand the quadratic terms
    for i in eachindex(coeffs)
        # make viable data objects so we can multiply the terms
        coef_data = MultiDiscreteMeasureData(prefs, [coeffs[i]], supps[:, i:i],
                                             name, w)
        simple_data = MultiDiscreteMeasureData(prefs, [1], supps[:, i:i], name,
                                               default_weight)
        new_expr = JuMP.@expression(write_model, new_expr + sum(coef *
                        expand_measure(pair.a, coef_data, write_model) *
                        expand_measure(pair.b, simple_data, write_model)
                        for (pair, coef) in expr.terms))
    end
    return new_expr
end

# MeasureRef
function expand_measure(mref::MeasureRef,
                        data::Union{DiscreteMeasureData,
                                    MultiDiscreteMeasureData},
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64}
    # determine function and data of the inner measure
    deeper_func = measure_function(mref)
    deeper_data = measure_data(mref)
    # expand the inner measure (note this is recursive for nested measures)
    new_func = expand_measure(deeper_func, deeper_data, write_model)
    # expand current level with the inner measure now expanded
    return expand_measure(new_func, data, write_model)
end

# Catch all method for undefined behavior
function expand_measure(expr, data::AbstractMeasureData,
                        write_model::JuMP.AbstractModel)
    expr_type = typeof(expr)
    data_type = typeof(data)
    error("Undefined behavior to expand expression of type `$expr_type` with " *
          "measure data `$data_type`. If this functionality is needed consider " *
          "extending `expand_measure`.")
end

"""
    expand(mref::MeasureRef)::JuMP.AbstractJuMPScalar

Return a JuMP scalar function containing the explicit expansion of the measure
`mref`. This expansion is done according to the measure data. Note that
variables are added to the model as necessary to accomodate the expansion (i.e.,
point variables and reduced infinite variables are made as needed). Errors if
expansion is undefined for the measure data and/or the measure expression. If
desired this can be used in combination with [`measure`](@ref) to expand measures
on the fly.

This is useful for extensions that employ a custom optimizer model since it
can be used evaluate measures before expressions are translated to the new model.
This method can also be extended to handle custom measure data types by extending
[`expand_measure`](@ref).

**Example**
```julia-repl
julia> tdata = DiscreteMeasureData(t, [0.5, 0.5], [0, 1])

julia> expr = expand(measure(g + z + T - h - 2, tdata))
0.5 g(0) + 0.5 g(1) + z + 0.5 T(0, x) + 0.5 T(1, x) - h(x) - 2
```
"""
function expand(mref::MeasureRef)::JuMP.AbstractJuMPScalar
    return expand_measure(measure_function(mref), measure_data(mref),
                          JuMP.owner_model(mref))
end


"""
    expand_measures(expr, write_model::JuMP.AbstractModel)::JuMP.AbstractJuMPScalar

Expand all `MeasureRef`s in `expr` in-place via [`expand_measure`](@ref) and
return the expanded expression. This is an internal method used by
[`expand_all_measures!`](@ref) and `TranscriptionOpt` but can be useful for
user-defined optimizer model extensions that add implement
[`add_measure_variable`](@ref) in combination with `expand_measure`. `write_model`
is the model that the measure variables are added to as described in
[`expand_measure`](@ref).
"""
function expand_measures end

# MeasureRef
function expand_measures(mref::MeasureRef,
                         write_model::JuMP.AbstractModel
                         )::JuMP.AbstractJuMPScalar
    return expand_measure(measure_function(mref), measure_data(mref),
                           write_model)
end

# GeneralVariableRef
function expand_measures(vref::GeneralVariableRef,
                         write_model::JuMP.AbstractModel
                         )::GeneralVariableRef
    return vref
end

# Generic[Aff/Quad]Expr (FiniteVariableRef)
function expand_measures(expr::JuMP._GenericAffOrQuadExpr{C, <:FiniteVariableRef},
                         write_model::JuMP.AbstractModel
                         )::JuMP.AbstractJuMPScalar where {C}
    return expr
end

# GenericAffExpr (GeneralVariableRef)
function expand_measures(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                         write_model::JuMP.AbstractModel
                         )::JuMP.AbstractJuMPScalar where {C}
    return JuMP.@expression(write_model, sum(coef * expand_measures(var, write_model)
                            for (var, coef) in expr.terms) + expr.constant)
end

# GenericQuadExpr (GeneralVariableRef)
function expand_measures(expr::JuMP.GenericQuadExpr{C, <:GeneralVariableRef},
                         write_model::JuMP.AbstractModel
                         )::JuMP.AbstractJuMPScalar where {C}
    return JuMP.@expression(write_model, sum(coef * expand_measures(pair.a, write_model) *
                expand_measures(pair.b, write_model) for (pair, coef) in expr.terms) +
                expand_measures(expr.aff, write_model))
end

# Fallback
function expand_measures(expr, write_model::JuMP.AbstractModel)
    error("`expand_measures` not defined for expressions of type `$(typeof(expr))`.")
end

"""
    expand_all_measures!(model::InfiniteModel)

Expand all of the measures used in the objective and/or constraints of `model`.
The objective and constraints are updated accordingly. Note that
variables are added to the model as necessary to accomodate the expansion (i.e.,
point variables and reduced infinite variables are made as needed). Errors if
expansion is undefined for the measure data and/or the measure expression. Also
errors if the expanded objective function is not finite.

This is useful for extensions that employ a custom optimizer model since it
can be used evaluate measures before `model` is translated into the new model.
This method can also be extended to handle custom measure data types by extending
[`expand_measure`](@ref). Note that this method leverages `expand_measure` via
[`expand_measures`](@ref).

**Example**
```julia-repl
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 T(t, xi) >= 0.0
 z >= 0.0
 g(t) + z >= 42.0
 measure(T(t, xi)) >= 0.0, for all xi in [-1, 1]
 t in [0, 6]
 xi in Normal(μ=0.0, σ=1.0)

julia> expand_all_measures!(model)

julia> print(model)
Min 3 g(6) + z
Subject to
 T(t, xi) >= 0.0
 z >= 0.0
 g(t) + z >= 42.0
 0.5 T(0, xi) + 0.5 T(6, xi) >= 0.0, for all xi in [-1, 1]
 t in [0, 6]
 xi in Normal(μ=0.0, σ=1.0)
```
"""
function expand_all_measures!(model::InfiniteModel)
    # expand the objective if it contains measures
    if JuMP.objective_function_type(model) <: MeasureExpr
        new_obj = expand_measures(JuMP.objective_function(model), model)
        isa(new_obj, InfiniteExpr) && error("Objective is not finite, ensure " *
                                            "all infinite variables/parameters " *
                                            "in measures are evaluated " *
                                            "completely.")
        JuMP.set_objective_function(model, new_obj)
    end
    # expand all of the constraints that contain measures
    for cindex in keys(model.constr_to_meas)
        # expand the expression
        new_func = expand_measures(model.constrs[cindex].func, model)
        # get the necessary info
        cref = _make_constraint_ref(model, cindex)
        name = JuMP.name(cref)
        set = model.constrs[cindex].set
        curr_index = model.next_constr_index
        # delete the old cosntraint and replace it with the expanded version
        model.next_constr_index = cindex - 1
        if isa(model.constrs[cindex], BoundedScalarConstraint) && isa(new_func,
                                                                   InfiniteExpr)
            orig_bounds = model.constrs[cindex].orig_bounds
            JuMP.delete(model, cref)
            JuMP.add_constraint(model, JuMP.build_constraint(error, new_func,
                                set; parameter_bounds = orig_bounds), name)
        else
            JuMP.delete(model, cref)
            JuMP.add_constraint(model, JuMP.build_constraint(error, new_func,
                                set), name)
        end
        model.next_constr_index = curr_index
    end
    return
end
