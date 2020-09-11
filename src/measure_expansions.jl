################################################################################
#                       MEASURE VARIABLE CREATION METHODS
################################################################################
"""
    make_point_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                            ivref::GeneralVariableRef,
                            support::Vector{Float64}
                            )::GenealVariableRef

Make a point variable for infinite variable `ivref` at `support`, add it to
the `write_model`, and return the `GeneralVariableRef`. This is an internal method
for point variables produced by expanding measures via [`expand_measure`](@ref).
This is also useful for those writing extension optimizer models and wish to
expand measures without modifiying the `InfiniteModel`. In such cases, `write_model`
should be the optimizer model and [`add_measure_variable`](@ref add_measure_variable(::JuMP.Model, ::Any, ::Any))
should be extended appropriately for point variables. Errors if `write_model` is
an optimizer model and `add_measure_variable` is not properly extended.
"""
function make_point_variable_ref(write_model::InfiniteModel,
                                 ivref::GeneralVariableRef,
                                 support::Vector{Float64})::GeneralVariableRef
    prefs = parameter_list(ivref)
    for i in eachindex(support)
        support[i] = round(support[i], sigdigits = significant_digits(prefs[i]))
    end
    base_info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false,
                                  NaN, false, false)
    new_info = _update_point_info(base_info, dispatch_variable_ref(ivref), support)
    var = PointVariable(_make_float_info(new_info), ivref, support)
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
adding `PointVariable`s and `ReducedVariable`s for such extensions.
Otherwise, an error is thrown for unextended variable and/or optimizer model types.
Note if this is extended, than [`internal_reduced_variable`](@ref) should also
be extended in order to direct reduced variables references to the underlying
[`ReducedVariable`](@ref).
"""
function add_measure_variable(model::JuMP.Model, var, key)
    error("`add_measure_variable` not defined for variable of type `$(typeof(var))` " *
          "and an optimizer model with key `$(typeof(key).parameters[1])`.")
end

# Store/add the variable to the optimizer model via add_measure_variable
# This avoids changing the InfiniteModel unexpectedly
function make_point_variable_ref(write_model::JuMP.Model, # this should be an optimizer model
                                 ivref::GeneralVariableRef,
                                 support::Vector{Float64})::GeneralVariableRef
    prefs = parameter_list(ivref)
    for i in eachindex(support)
        support[i] = round(support[i], sigdigits = significant_digits(prefs[i]))
    end
    base_info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false,
                                  NaN, false, false)
    new_info = _update_point_info(base_info, dispatch_variable_ref(ivref), support)
    var = PointVariable(_make_float_info(new_info), ivref, support)
    opt_key = optimizer_model_key(write_model)
    return add_measure_variable(write_model, var, Val(opt_key))
end

"""
    make_reduced_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                              ivref::GeneralVariableRef,
                              indices::Vector{Int},
                              values::Vector{Float64}
                              )::GeneralVariableRef

Make a reduced variable for infinite variable `ivref` at `support`, add it to
the `write_model`, and return the `GeneralVariableRef`. This is an internal method
for reduced variables produced by expanding measures via [`expand_measure`](@ref).
This is also useful for those writing extension optimizer models and wish to
expand measures without modifiying the `InfiniteModel`. In such cases, `write_model`
should be the optimizer model and [`add_measure_variable`](@ref add_measure_variable(::JuMP.Model, ::Any, ::Any))
should be extended appropriately for reduced variables. Errors if `write_model`
is an optimizer model and `add_measure_variable` is not properly extended.
Note this is only intended for optimizer models that are currently stored in
`InfiniteModel.optimizer_model`.
"""
function make_reduced_variable_ref(write_model::InfiniteModel,
                                   ivref::GeneralVariableRef,
                                   indices::Vector{Int},
                                   values::Vector{Float64}
                                   )::GeneralVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    var = JuMP.build_variable(error, ivref, eval_supps, check = false)
    return JuMP.add_variable(write_model, var)
end

# Add reduced infinite variables in the optimizer model without modifying the InfiniteModel
function make_reduced_variable_ref(write_model::JuMP.Model,
                                   ivref::GeneralVariableRef,
                                   indices::Vector{Int},
                                   values::Vector{Float64}
                                   )::GeneralVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    var = JuMP.build_variable(error, ivref, eval_supps, check = false)
    key = optimizer_model_key(write_model)
    return add_measure_variable(write_model, var, Val(key))
end

"""
    delete_internal_reduced_variable(write_model::Union{InfiniteModel, JuMP.Model},
                                     rvref::ReducedVariableRef)::Nothing

Delete the variable associated with `rvref` from `write_model` if it is purely
an internal variable only used for measure expansion and is no longer needed.
For `write_model`s that are an optimizer model,
[`delete_reduced_variable`](@ref delete_reduced_variable(::JuMP.Model, ::Any,::Any))
will need to be extended for this this to work. Otherwise, a warning will be thrown.
Note that this is intended as an internal method to assist with extensions to
[`expand_measure`](@ref).
"""
function delete_internal_reduced_variable(write_model::InfiniteModel,
                                          rvref::ReducedVariableRef)::Nothing
    if !used_by_measure(rvref) && !used_by_constraint(rvref)
        JuMP.delete(write_model, rvref)
    end
    return
end

"""
    delete_reduced_variable(model::JuMP.Model, vref, key::Val{:ext_key_name})::Nothing

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
                                          rvref::ReducedVariableRef)::Nothing
    if !JuMP.is_valid(JuMP.owner_model(rvref), rvref)
        key = optimizer_model_key(write_model)
        delete_reduced_variable(write_model, rvref, Val(key))
    end
    return
end

################################################################################
#                          EXPAND_MEASURE DEFINITIONS
################################################################################
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

# GeneralVariableRef
function expand_measure(vref::GeneralVariableRef,
                        data::DiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64}
    return expand_measure(vref, _index_type(vref), data, write_model)
end

# InfiniteVariableRef (1D DiscreteMeasureData)
function expand_measure(ivref::GeneralVariableRef,
                        index_type::Type{InfiniteVariableIndex},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
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
        return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, ivref => var_coef)
    # make point variables if var_prefs = pref (it is the only dependence)
    elseif length(var_prefs) == 1
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_point_variable_ref(write_model, ivref, [supps[i]])
                    for i in eachindex(coeffs)))
    # make reduced variables if the variable contains other parameters
    else
        index = [findfirst(isequal(pref), var_prefs)]
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_reduced_variable_ref(write_model, ivref, index, [supps[i]])
                    for i in eachindex(coeffs)))
    end
end

# InfiniteVariableRef (Multi DiscreteMeasureData)
function expand_measure(ivref::GeneralVariableRef,
                        index_type::Type{InfiniteVariableIndex},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    var_prefs = parameter_list(ivref)
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # var_prefs == prefs so let's make a point variable
    if var_prefs == prefs
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_point_variable_ref(write_model, ivref, supps[:, i])
                    for i in eachindex(coeffs)))
    # treat variable as constant if doesn't have measure parameter
    elseif !any(pref in var_prefs for pref in prefs)
        var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, ivref => var_coef)
    # make point variables if all var_prefs are contained in prefs
    elseif all(pref in prefs for pref in var_prefs)
        indices = [findfirst(isequal(pref), prefs) for pref in var_prefs]
        new_supps = supps[indices, :]
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_point_variable_ref(write_model, ivref, new_supps[:, i])
                    for i in eachindex(coeffs)))
    # make reduced variables if the variable contains other parameters
    else
        # get indices of each pref to map properly
        indices = [findfirst(isequal(pref), var_prefs) for pref in prefs]
        # check that if any of the indices are empty and truncate as needed
        empty = map(i -> i === nothing, indices)
        if any(empty)
            indices = convert(Vector{Int}, deleteat!(indices, empty))
            supps = supps[.!empty, :]
        end
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, supps[:, i])
                    for i in eachindex(coeffs)))
    end
end

# Write point support given reduced info and the support
function _make_point_support(orig_prefs::Vector{GeneralVariableRef},
                             support_dict::Dict{Int, Float64},
                             index::Int, value::Float64)::Vector{Float64}
    return [i == index ? value : support_dict[i] for i in eachindex(orig_prefs)]
end

# ReducedVariableRef (1D DiscreteMeasureData)
function expand_measure(rvref::GeneralVariableRef,
                        index_type::Type{ReducedVariableIndex},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    drvref = dispatch_variable_ref(rvref)
    ivref = infinite_variable_ref(drvref)
    orig_prefs = parameter_list(ivref)
    var_prefs = parameter_list(drvref)
    eval_supps = eval_supports(drvref)
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat variable as constant if doesn't have measure parameter
    if !(pref in var_prefs)
        var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        expr = JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, rvref => var_coef)
    # make point variables if var_prefs = pref (it is the only dependence)
    elseif length(var_prefs) == 1
        index = findfirst(isequal(pref), orig_prefs)
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_point_variable_ref(write_model, ivref,
                    _make_point_support(orig_prefs, eval_supps, index, supps[i]))
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, drvref) # TODO not sure this is helpful
    # make reduced variables if the variable contains other parameters
    else
        index = findfirst(isequal(pref), orig_prefs)
        collected_indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], collected_indices) # a support will be appended on the fly
        indices = push!(collected_indices, index)
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, vcat(vals, supps[i]))
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, drvref) # TODO not sure this is helpful
    end
    return expr
end

# Write point support given reduced info and the support
function _make_point_support(orig_prefs::Vector{GeneralVariableRef},
                             support_dict::Dict{Int, Float64},
                             new_indices::Vector{Int},
                             values::Vector{Float64})::Vector{Float64}
    # these might overlap with the old dict, so we favor the old dict
    new_dict = Dict(new_indices[i] => values[i] for i in eachindex(values))
    return [haskey(support_dict, i) ? support_dict[i] : new_dict[i]
            for i in eachindex(orig_prefs)]
end

# ReducedVariableRef (Multi DiscreteMeasureData)
function expand_measure(rvref::GeneralVariableRef,
                        index_type::Type{ReducedVariableIndex},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr
    # pull in the needed information
    drvref = dispatch_variable_ref(rvref)
    ivref = infinite_variable_ref(drvref)
    orig_prefs = parameter_list(ivref)
    var_prefs = parameter_list(drvref)
    eval_supps = eval_supports(drvref)
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat variable as constant if doesn't have measure parameter
    if !any(pref in var_prefs for pref in prefs)
        var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        expr = JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, rvref => var_coef)
    # make point variables if prefs includes all var_prefs
    elseif all(pref in prefs for pref in var_prefs)
        # get the indices of measure prefs to reorder/truncate as needed
        supp_indices = [findfirst(isequal(pref), prefs) for pref in var_prefs]
        # reorder/truncate if necesary
        if supp_indices != 1:size(supps, 1)
            supps = supps[supp_indices, :]
        end
        # get the parameter indices of the variable parameters to be reduced
        indices = [findfirst(isequal(pref), orig_prefs) for pref in var_prefs]
        # make the expression
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_point_variable_ref(write_model, ivref,
                    _make_point_support(orig_prefs, eval_supps, indices, supps[:, i]))
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, drvref) # TODO not sure this is helpful
    # make reduced variables if the variable contains other parameters
    else
        # get the indices of prefs in terms of the ivref
        new_indices = [findfirst(isequal(pref), orig_prefs) for pref in prefs]
        # check that if any of the indices are empty or already reduced and truncate as needed
        bad_index = map(i -> i === nothing || i in keys(eval_supps), new_indices)
        if any(bad_index)
            deleteat!(new_indices, bad_index)
            supps = supps[.!bad_index, :]
        end
        # prepare the indices and values for reduced variable construction
        collected_indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], collected_indices) # a support will be appended on the fly
        indices = append!(collected_indices, new_indices)
        # make the expression
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_reduced_variable_ref(write_model, ivref, indices, vcat(vals, supps[:, i]))
                    for i in eachindex(coeffs)))
        delete_internal_reduced_variable(write_model, drvref) # TODO not sure this is helpful
    end
    return expr
end

# FiniteVariableRef (1D DiscreteMeasureData)
function expand_measure(vref::GeneralVariableRef,
                        index_type::Type{V},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr where {V <: FiniteVariableIndex}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the variable as a constant and build the expression
    var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
    return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, vref => var_coef)
end

# FiniteVariableRef (Multi DiscreteMeasureData)
function expand_measure(vref::GeneralVariableRef,
                        index_type::Type{V},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr where {V <: FiniteVariableIndex}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the variable as a constant and build the expression
    var_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
    return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, vref => var_coef)
end

# InfiniteParameterRef (1D DiscreteMeasureData)
function expand_measure(pref::GeneralVariableRef,
                        index_type::Type{P},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.GenericAffExpr, Float64} where {P <: InfiniteParameterIndex}
    # pull in the needed information
    meas_pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the parameter
    if meas_pref != pref
        par_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, pref => par_coef)
    # replace the parameter with its value
    else
        return sum(coeffs[i] * w(supps[i]) * supps[i] for i in eachindex(coeffs))
    end
end

# InfiniteParameterRef (Multi DiscreteMeasureData)
function expand_measure(pref::GeneralVariableRef,
                        index_type::Type{P},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.GenericAffExpr, Float64} where {P <: InfiniteParameterIndex}
    # pull in the needed information
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # find the position of pref if it is in the data
    index = findfirst(isequal(pref), prefs)
    # treat the parameter as a constant
    if index === nothing
        par_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, pref => par_coef)
    # replace the parameter with its value
    else
        return sum(coeffs[i] * w(supps[:, i]) * supps[index, i] for i in eachindex(coeffs))
    end
end

# GenericAffExpr (1D DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericAffExpr{C, GeneralVariableRef},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
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

# GenericAffExpr (Multi DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericAffExpr{C, GeneralVariableRef},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
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

# GenericQuadExpr (1D DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericQuadExpr{C, GeneralVariableRef},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # get needed info
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lb = JuMP.lower_bound(data)
    ub = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # expand the affine expression
    new_expr = expand_measure(expr.aff, data, write_model)
    # make viable data objects so we can multiply the terms
    coef_data = DiscreteMeasureData(pref, ones(1), ones(1), label, w, lb, ub,
                                    is_expect)
    simple_data = DiscreteMeasureData(pref, ones(1), ones(1), label,
                                      default_weight, lb, ub, is_expect)
    # expand the quadratic terms
    for i in eachindex(coeffs)
        # update the temp data
        coefficients(coef_data)[1] = coeffs[i]
        supports(coef_data)[1] = supps[i]
        supports(simple_data)[1] = supps[i]
        new_expr = JuMP.@expression(write_model, new_expr + sum(coef *
                        expand_measure(pair.a, coef_data, write_model) *
                        expand_measure(pair.b, simple_data, write_model)
                        for (pair, coef) in expr.terms))
    end
    return new_expr
end

# GenericQuadExpr(Multi DiscreteMeasureData)
function expand_measure(expr::JuMP.GenericQuadExpr{C, GeneralVariableRef},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64} where {C}
    # get needed info
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # expand the affine expression
    new_expr = expand_measure(expr.aff, data, write_model)
    # make viable data objects so we can multiply the terms
    coef_data = DiscreteMeasureData(prefs, ones(1), ones(length(prefs), 1),
                                    label, w, lbs, ubs, is_expect)
    simple_data = DiscreteMeasureData(prefs, ones(1), ones(length(prefs), 1),
                                      label, default_weight, lbs, ubs, is_expect)
    # expand the quadratic terms
    for i in eachindex(coeffs)
        # update the temp data
        coefficients(coef_data)[1] = coeffs[i]
        supports(coef_data)[:, 1] = @view(supps[:, i])
        supports(simple_data)[:, 1] = @view(supps[:, i])
        new_expr = JuMP.@expression(write_model, new_expr + sum(coef *
                        expand_measure(pair.a, coef_data, write_model) *
                        expand_measure(pair.b, simple_data, write_model)
                        for (pair, coef) in expr.terms))
    end
    return new_expr
end

# MeasureRef
function expand_measure(mref::GeneralVariableRef,
                        index_type::Type{MeasureIndex},
                        data::DiscreteMeasureData,
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

# FunctionalDiscreteMeasureData
function expand_measure(expr,
                        data::FunctionalDiscreteMeasureData,
                        write_model::JuMP.AbstractModel
                        )::Union{JuMP.AbstractJuMPScalar, Float64}
    # get the info
    prefs = parameter_refs(data)
    supps = supports(data)
    coef_func = coefficient_function(data)
    coeffs = coef_func(supps)
    label = support_label(data)
    w = weight_function(data)
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # prepare the DiscreteMeasureData and dispatch
    new_data = DiscreteMeasureData(prefs, coeffs, supps, label, w, lbs, ubs,
                                   is_expect)
    return expand_measure(expr, new_data, write_model)
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

################################################################################
#                          ANALYTIC EXPANSION METHODS
################################################################################
"""
    analytic_expansion(expr, data::AbstractMeasureData,
                       write_model::JuMP.AbstractModel)::JuMP.AbstractJuMPScalar

Analytically, evaluate measure in the simple case where the measure expression
`expr` doesn't depend on `data` and thus `expr` can be treated as a constant in
conjunction with an analytic result of the `data`. This is intended as an internal
method that is used by [`expand`](@ref) and [`expand_measures`](@ref). For
unrecognized `data` types, `expand_measure` is called instead. User defined
measure data type may choose to extend this method if desired. This is triggered
when `is_analytic(mref) = true`.
"""
function analytic_expansion end

# 1D DiscreteMeasureData/FunctionalDiscreteMeasureData
function analytic_expansion(expr::JuMP.AbstractJuMPScalar,
                            data::Union{DiscreteMeasureData{GeneralVariableRef, 1},
                                        FunctionalDiscreteMeasureData{GeneralVariableRef}},
                            write_model::JuMP.AbstractModel # needed for fallback
                            )::JuMP.AbstractJuMPScalar
    # get the bounds and expect
    lb = JuMP.lower_bound(data)
    ub = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    w = weight_function(data)
    # build the expression
    if w == default_weight && !isnan(lb) && !isnan(ub)
        return (ub - lb) * expr
    elseif w == default_weight && is_expect
        return expr # this is an expectation and should be equal to 1 * expr
    else
        supps = supports(data)
        coeffs = coefficients(data)
        constant_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
        return constant_coef * expr
    end
end

# Multi DiscreteMeasureData/FunctionalDiscreteMeasureData
function analytic_expansion(expr::JuMP.AbstractJuMPScalar,
                            data::Union{DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                                        FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}},
                            write_model::JuMP.AbstractModel # needed for fallback
                            )::JuMP.AbstractJuMPScalar
    # get the bounds and expect
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    w = weight_function(data)
    # build the expression
    if w == default_weight && !isnan(first(lbs)) && !isnan(first(ubs))
        return prod(ubs .- lbs) * expr
    elseif w == default_weight && is_expect
        return expr # this is an expectation and should be equal to 1 * expr
    else
        supps = supports(data)
        coeffs = coefficients(data)
        constant_coef = sum(coeffs[i] * w(supps[:, i]) for i in eachindex(coeffs))
        return constant_coef * expr
    end
end

# Fallback
function analytic_expansion(expr, data::AbstractMeasureData,
                            write_model::JuMP.AbstractModel
                            )::JuMP.AbstractJuMPScalar
    return expand_measure(expr, data, write_model)
end

################################################################################
#                              EXPANSION METHODS
################################################################################
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
[`expand_measure`](@ref). Optionally, [`analytic_expansion`](@ref) can also
be extended which is triggered by [`is_analytic`](@ref) for such types if
analytic expansion is possible in certain cases.

**Example**
```julia-repl
julia> tdata = DiscreteMeasureData(t, [0.5, 0.5], [0, 1])

julia> expr = expand(measure(g + z + T - h - 2, tdata))
0.5 g(0) + 0.5 g(1) + z + 0.5 T(0, x) + 0.5 T(1, x) - h(x) - 2
```
"""
function expand(mref::MeasureRef)::Union{JuMP.AbstractJuMPScalar, Float64}
    if is_analytic(mref)
        return analytic_expansion(measure_function(mref), measure_data(mref),
                                  JuMP.owner_model(mref))
    else
        return expand_measure(measure_function(mref), measure_data(mref),
                          JuMP.owner_model(mref))
    end
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
function expand_measures(mref::GeneralVariableRef,
                         index_type::Type{MeasureIndex},
                         write_model::JuMP.AbstractModel
                         )::Union{JuMP.AbstractJuMPScalar, Float64}
    if is_analytic(mref)
        return analytic_expansion(measure_function(mref), measure_data(mref),
                                  write_model)
    else
        return expand_measure(measure_function(mref), measure_data(mref),
                              write_model)
    end
end

# NonMeasureRef
function expand_measures(vref::GeneralVariableRef,
                         index_type::Type{V},
                         write_model::JuMP.AbstractModel
                         )::GeneralVariableRef where {V <: AbstractInfOptIndex}
    return vref
end

# GeneralVariableRef
function expand_measures(vref::GeneralVariableRef,
                         write_model::JuMP.AbstractModel
                         )::Union{JuMP.AbstractJuMPScalar, Float64}
    return expand_measures(vref, _index_type(vref), write_model)
end

# GenericAffExpr
function expand_measures(expr::JuMP.GenericAffExpr{C, GeneralVariableRef},
                         write_model::JuMP.AbstractModel
                         )::JuMP.AbstractJuMPScalar where {C}
    return JuMP.@expression(write_model, sum(coef * expand_measures(var, write_model)
                            for (var, coef) in expr.terms) + expr.constant)
end

# GenericQuadExpr
function expand_measures(expr::JuMP.GenericQuadExpr{C, GeneralVariableRef},
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
    expand_all_measures!(model::InfiniteModel)::Nothing

Expand all of the measures used in the objective and/or constraints of `model`.
The objective and constraints are updated accordingly. Note that
variables are added to the model as necessary to accomodate the expansion (i.e.,
point variables and reduced infinite variables are made as needed). Errors if
expansion is undefined for the measure data and/or the measure expression.

This is useful for extensions that employ a custom optimizer model since it
can be used evaluate measures before `model` is translated into the new model.
This method can also be extended to handle custom measure data types by extending
[`expand_measure`](@ref). Note that this method leverages `expand_measure` via
[`expand_measures`](@ref). Optionally, [`analytic_expansion`](@ref) can also
be extended which is triggered by [`is_analytic`](@ref) for such types if
analytic expansion is possible in certain cases.

**Example**
```julia-repl
julia> print(model)
Min integral{t ∈ [0, 6]}[g(t)*t] + z
Subject to
 T(t, x) ≥ 0.0, ∀ t ∈ [0, 6], xi ∈ [-1, 1]
 z ≥ 0.0
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]
 integral{t ∈ [0, 6]}[T(t, x)] ≥ 0.0, ∀ x ∈ [-1, 1]

julia> expand_all_measures!(model)

julia> print(model)
Min 3 g(6) + z
Subject to
 T(t, x) ≥ 0.0, ∀ t ∈ [0, 6], xi ∈ [-1, 1]
 z ≥ 0.0
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]
 0.5 T(0, x) + 0.5 T(6, xi) ≥ 0.0, ∀ x ∈ [-1, 1]
```
"""
function expand_all_measures!(model::InfiniteModel)::Nothing
    # expand the objective if it contains measures
    if objective_has_measures(model)
        new_obj = expand_measures(JuMP.objective_function(model), model)
        JuMP.set_objective_function(model, new_obj)
    end
    # expand all of the constraints that contain measures
    for (cindex, object) in model.constraints
        if !isempty(object.measure_indices)
            old_constr = object.constraint
            # clear the old dependencies
            old_vrefs = _all_function_variables(JuMP.jump_function(old_constr))
            for vref in old_vrefs
                filter!(e -> e != cindex, _constraint_dependencies(vref))
            end
            # expand the expression
            new_func = expand_measures(JuMP.jump_function(old_constr), model)
            offset = JuMP.constant(new_func)
            JuMP.add_to_expression!(new_func, -offset)
            new_set = MOIU.shift_constant(JuMP.moi_set(old_constr), -offset)
            vrefs = _all_function_variables(new_func)
            # make the new constraint object
            if old_constr isa BoundedScalarConstraint
                orig_bounds = original_parameter_bounds(old_constr)
                new_constr = BoundedScalarConstraint(new_func, new_set,
                                                     copy(orig_bounds),
                                                     orig_bounds)
            else
                new_constr = JuMP.ScalarConstraint(new_func, new_set)
            end
            # update the bounds if there are bounded hold variables in the model
            if model.has_hold_bounds
                new_constr = _check_and_update_bounds(model, new_constr, vrefs)
            end
            # update the constraint data
            cref = _temp_constraint_ref(model, cindex)
            _set_core_constraint_object(cref, new_constr)
            empty!(object.measure_indices)
            _update_var_constr_mapping(vrefs, cref)
        end
    end
    return
end
