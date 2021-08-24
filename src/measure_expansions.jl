################################################################################
#                       MEASURE VARIABLE CREATION METHODS
################################################################################
"""
    make_point_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                            ivref::GeneralVariableRef,
                            support::Vector{Float64}
                            )::GeneralVariableRef

Make a point variable for infinite variable/derivative `ivref` at `support`, add it to
the `write_model`, and return the `GeneralVariableRef`. This is an internal method
for point variables produced by expanding measures via [`expand_measure`](@ref).
This is also useful for those writing extension optimizer models and wish to
expand measures without modifiying the `InfiniteModel`. In such cases, `write_model`
should be the optimizer model and 
[`add_point_variable`](@ref add_point_variable(::JuMP.Model, ::Any, ::Any, ::Any)) 
should be extended appropriately for point variables. Errors if `write_model` is
an optimizer model and `add_point_variable` is not properly extended. 

Note this is also accomodates infinite parameter functions, in which case the 
infinite parameter function is called with the support as input. 
"""
function make_point_variable_ref(
    write_model::InfiniteModel,
    ivref::GeneralVariableRef,
    support::Vector{Float64}
    )
    return make_point_variable_ref(write_model, ivref, support, _index_type(ivref))
end

# Infinite variable index
function make_point_variable_ref(
    write_model::InfiniteModel, 
    ivref, 
    support, 
    ::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}}
    )::GeneralVariableRef
    prefs = parameter_list(ivref)
    for i in eachindex(support)
        support[i] = round(support[i], sigdigits = significant_digits(prefs[i]))
    end
    pindex = get(write_model.point_lookup, (ivref, support), nothing)
    if pindex === nothing
        base_info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false,
                                      NaN, false, false)
        info = _update_point_info(base_info, dispatch_variable_ref(ivref), support)
        var = PointVariable(_make_float_info(info), ivref, support)
        return JuMP.add_variable(write_model, var; add_support = false)
    else
        return _make_variable_ref(write_model, pindex)
    end
end

# Infinite parameter function index (this works for both optimizer and infinite models)
function make_point_variable_ref(
    write_model, 
    fref, 
    support, 
    ::Type{ParameterFunctionIndex}
    )::Float64
    prefs = raw_parameter_refs(fref)
    supp_tuple = Tuple(support, prefs)
    func = raw_function(fref)
    return func(supp_tuple...)
end

"""
    add_point_variable(model::JuMP.Model, ivref::GeneralVariableRef, 
                       support::Vector{Float64}, key::Val{:ext_key_name}
                       )::GeneralVariableRef

Add a point variable (defined by restricting `ivref` to `support`) to the 
optimizer model `model` (with `key`) and return the correct `InfiniteOpt` 
variable reference. This is an internal method used by 
[`make_point_variable_ref`](@ref) to make point variables when the `write_model` 
is an optimizer model. This is useful for extensions that wish to expand 
measures, but without changing the original `InfiniteModel`. An error is thrown 
for unextended optimizer model types.
"""
function add_point_variable(model::JuMP.Model, ivref, supp, key)
    error("`add_point_variable` not defined for an optimizer model with key ",
          "`$(typeof(key).parameters[1])`.")
end

# Store/add the variable to the optimizer model via add_measure_variable
# This avoids changing the InfiniteModel unexpectedly
function make_point_variable_ref(
    write_model::JuMP.Model, # this should be an optimizer model
    ivref::GeneralVariableRef,
    support::Vector{Float64}
    )
    return make_point_variable_ref(write_model, ivref, support, _index_type(ivref))
end

# Infinite variable index
function make_point_variable_ref(
    write_model::JuMP.Model, 
    ivref, 
    support, 
    ::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}}
    )::GeneralVariableRef
    prefs = parameter_list(ivref)
    for i in eachindex(support)
        support[i] = round(support[i], sigdigits = significant_digits(prefs[i]))
    end 
    opt_key = optimizer_model_key(write_model)
    return add_point_variable(write_model, ivref, support, Val(opt_key))
end

"""
    make_semi_infinite_variable_ref(write_model::Union{InfiniteModel, JuMP.Model},
                                    ivref::GeneralVariableRef,
                                    indices::Vector{Int},
                                    values::Vector{Float64}
                                    )::GeneralVariableRef

Make a semi-infinite variable for infinite variable/derivative/parameter 
function `ivref` at `support`, add it to the `write_model`, and return the 
`GeneralVariableRef`. This is an internal method for semi-infinite variables 
produced by expanding measures via [`expand_measure`](@ref). This is also useful 
for those writing extension optimizer models and wish to expand measures without 
modifiying the `InfiniteModel`. In such cases, `write_model` should be the 
optimizer model and 
[`add_semi_infinite_variable`](@ref add_semi_infinite_variable(::JuMP.Model, ::Any, ::Any)) 
should be extended appropriately for semi-infinite variables. Errors if 
`write_model` is an optimizer model and `add_semi_infinite_variable` is not 
properly extended. Note this is only intended for optimizer models that are 
currently stored in `InfiniteModel.optimizer_model`.
"""
function make_semi_infinite_variable_ref(
    write_model::InfiniteModel,
    ivref::GeneralVariableRef,
    indices::Vector{Int},
    values::Vector{Float64}
    )::GeneralVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    existing_index = get(write_model.semi_lookup, (ivref, eval_supps), nothing)
    if existing_index === nothing
        var = JuMP.build_variable(error, ivref, eval_supps, check = false)
        return JuMP.add_variable(write_model, var, add_support = false)
    else 
        return _make_variable_ref(write_model, existing_index)
    end
end

"""
    add_semi_infinite_variable(model::JuMP.Model, var::SemiInfiniteVariable, 
                               key::Val{:ext_key_name})::GeneralVariableRef

Add a semi-infinite variable `var` to the optimizer model `model` (with `key`) 
and return the correct `InfiniteOpt` variable reference. This is an internal 
method used by [`make_semi_infinite_variable_ref`](@ref) to make semi-infinite 
variables when the `write_model` is an optimizer model. This is useful for 
extensions that wish to expand measures, but without changing the original 
`InfiniteModel`. An error is thrown for optimizer model types. Note if this is 
extended, than [`internal_semi_infinite_variable`](@ref) should also be extended 
in order to direct semi-infinite variables references to the underlying 
[`SemiInfiniteVariable`](@ref).
"""
function add_semi_infinite_variable(model::JuMP.Model, var, key)
    error("`add_semi_infinite_variable` not defined for an optimizer model ",
          "with key `$(typeof(key).parameters[1])`.")
end


# Add semi-infinite infinite variables in the optimizer model without modifying the InfiniteModel
function make_semi_infinite_variable_ref(
    write_model::JuMP.Model,
    ivref::GeneralVariableRef,
    indices::Vector{Int},
    values::Vector{Float64}
    )::GeneralVariableRef
    eval_supps = Dict(indices[i] => values[i] for i in eachindex(indices))
    var = JuMP.build_variable(error, ivref, eval_supps, check = false)
    key = optimizer_model_key(write_model)
    return add_semi_infinite_variable(write_model, var, Val(key))
end

# Map a variable to a new one given a new one (TODO make more efficient)
function _map_variable(vref, data, supp::Float64, write_model)
    data.supports[1] = supp
    return expand_measure(vref, data, write_model)
end

function _map_variable(vref, data, supp, write_model)
    data.supports[:, 1] = supp
    return expand_measure(vref, data, write_model)
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
and [`make_semi_infinite_variable_ref`](@ref) should be used as appropriate to create
these variables. Note this is intended as an internal function,
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

# InfiniteVariableRef/DerivativeRef/ParameterFunctionRef (1D DiscreteMeasureData)
function expand_measure(ivref::GeneralVariableRef,
                        index_type::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}, Type{ParameterFunctionIndex}},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel)
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
    # make semi-infinite variables if the variable contains other parameters
    else
        index = [findfirst(isequal(pref), var_prefs)]
        return JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_semi_infinite_variable_ref(write_model, ivref, index, [supps[i]])
                    for i in eachindex(coeffs)))
    end
end

# InfiniteVariableRef/DerivativeRef/ParameterFunctionRef (Multi DiscreteMeasureData)
function expand_measure(ivref::GeneralVariableRef,
                        index_type::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}, Type{ParameterFunctionIndex}},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel)
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
    # make semi-infinite variables if the variable contains other parameters
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
                    make_semi_infinite_variable_ref(write_model, ivref, indices, supps[:, i])
                    for i in eachindex(coeffs)))
    end
end

# Write point support given semi_infinite info and the support
function _make_point_support(
    orig_prefs::Vector{GeneralVariableRef},
    support_dict::Dict{Int, Float64},
    idx::Int, 
    value::Float64
    )::Vector{Float64}
    return [i == idx ? value : support_dict[i] for i in eachindex(orig_prefs)]
end

# SemiInfiniteVariableRef (1D DiscreteMeasureData)
function expand_measure(rvref::GeneralVariableRef,
                        index_type::Type{SemiInfiniteVariableIndex},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel)
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
    # make semi-infinite variables if the variable contains other parameters
    else
        index = findfirst(isequal(pref), orig_prefs)
        collected_indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], collected_indices) # a support will be appended on the fly
        indices = push!(collected_indices, index)
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[i]) *
                    make_semi_infinite_variable_ref(write_model, ivref, indices, vcat(vals, supps[i]))
                    for i in eachindex(coeffs)))
    end
    return expr
end

# Write point support given semi_infinite info and the support
function _make_point_support(orig_prefs::Vector{GeneralVariableRef},
                             support_dict::Dict{Int, Float64},
                             new_indices::Vector{Int},
                             values::Vector{Float64})::Vector{Float64}
    # these might overlap with the old dict, so we favor the old dict
    new_dict = Dict(new_indices[i] => values[i] for i in eachindex(values))
    return [haskey(support_dict, i) ? support_dict[i] : new_dict[i]
            for i in eachindex(orig_prefs)]
end

# SemiInfiniteVariableRef (Multi DiscreteMeasureData)
function expand_measure(rvref::GeneralVariableRef,
                        index_type::Type{SemiInfiniteVariableIndex},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel)
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
    # make semi-infinite variables if the variable contains other parameters
    else
        # get the indices of prefs in terms of the ivref
        new_indices = [findfirst(isequal(pref), orig_prefs) for pref in prefs]
        # check that if any of the indices are empty or already reduced and truncate as needed
        bad_index = map(i -> i === nothing || i in keys(eval_supps), new_indices)
        if any(bad_index)
            deleteat!(new_indices, bad_index)
            supps = supps[.!bad_index, :]
        end
        # prepare the indices and values for semi-infinite variable construction
        collected_indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], collected_indices) # a support will be appended on the fly
        indices = append!(collected_indices, new_indices)
        # make the expression
        expr = JuMP.@expression(write_model, sum(coeffs[i] * w(supps[:, i]) *
                    make_semi_infinite_variable_ref(write_model, ivref, indices, vcat(vals, supps[:, i]))
                    for i in eachindex(coeffs)))
    end
    return expr
end

# FiniteRef (1D DiscreteMeasureData)
function expand_measure(vref::GeneralVariableRef,
                        index_type::Type{V},
                        data::DiscreteMeasureData{GeneralVariableRef, 1},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr where {V <: FiniteIndex}
    # pull in the needed information
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    # treat the variable as a constant and build the expression
    var_coef = sum(coeffs[i] * w(supps[i]) for i in eachindex(coeffs))
    return JuMP.GenericAffExpr{Float64, GeneralVariableRef}(0, vref => var_coef)
end

# FiniteRef (Multi DiscreteMeasureData)
function expand_measure(vref::GeneralVariableRef,
                        index_type::Type{V},
                        data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                        write_model::JuMP.AbstractModel
                        )::JuMP.GenericAffExpr where {V <: FiniteIndex}
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
function expand_measure(
    expr::JuMP.GenericQuadExpr,
    data::DiscreteMeasureData{GeneralVariableRef, 1},
    write_model::JuMP.AbstractModel
    )
    # get needed info
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lb = JuMP.lower_bound(data)
    ub = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # make the expression
    simple_data = DiscreteMeasureData(pref, ones(1), ones(1), label,
                                      default_weight, lb, ub, is_expect)
    return _MA.@rewrite(sum(sum(coeffs[i] * w(supps[i]) * c *
                        _map_variable(p.a, simple_data, supps[i], write_model) * 
                        _map_variable(p.b, simple_data, supps[i], write_model) 
                        for (p, c) in expr.terms) for i in eachindex(coeffs)) + 
                        expand_measure(expr.aff, data, write_model))
end

# GenericQuadExpr(Multi DiscreteMeasureData)
function expand_measure(
    expr::JuMP.GenericQuadExpr,
    data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
    write_model::JuMP.AbstractModel
    )
    # get needed info
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # make the expression
    simple_data = DiscreteMeasureData(prefs, ones(1), ones(length(prefs), 1),
                                      label, default_weight, lbs, ubs, is_expect)
    return _MA.@rewrite(sum(sum(coeffs[i] * w(@view(supps[:, i])) * c *
                        _map_variable(p.a, simple_data, @view(supps[:, i]), write_model) * 
                        _map_variable(p.b, simple_data, @view(supps[:, i]), write_model) 
                        for (p, c) in expr.terms) for i in eachindex(coeffs)) + 
                        expand_measure(expr.aff, data, write_model))
end

# NLPExpr (1D DiscreteMeasureData)
function expand_measure(
    expr::NLPExpr,
    data::DiscreteMeasureData{GeneralVariableRef, 1},
    write_model::JuMP.AbstractModel
    )
    # get needed info
    pref = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lb = JuMP.lower_bound(data)
    ub = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # make the expression
    simple_data = DiscreteMeasureData(pref, ones(1), ones(1), label,
                                      default_weight, lb, ub, is_expect)
    return _MA.@rewrite(nlsum(coeffs[i] * w(supps[i]) * 
            map_expression(v -> _map_variable(v, simple_data, supps[i], write_model), expr) 
            for i in eachindex(supps)))
end

# NLPExpr (Multi DiscreteMeasureData)
function expand_measure(
    expr::NLPExpr,
    data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
    write_model::JuMP.AbstractModel
    )
    # get needed info
    prefs = parameter_refs(data)
    supps = supports(data)
    coeffs = coefficients(data)
    w = weight_function(data)
    label = support_label(data)
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    is_expect = _is_expect(data)
    # make the expression
    simple_data = DiscreteMeasureData(prefs, ones(1), ones(length(prefs), 1),
                                      label, default_weight, lbs, ubs, is_expect)
    return _MA.@rewrite(nlsum(coeffs[i] * w(@view(supps[:, i])) * 
            map_expression(v -> _map_variable(v, simple_data, @view(supps[:, i]), write_model), expr) 
            for i in eachindex(coeffs)))
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

# Call add_generative_supports if needed 
function _prep_generative_supps(prefs, info_type::Type{NoGenerativeSupports})::Nothing 
    return
end
function _prep_generative_supps(pref, info_type)::Nothing 
    add_generative_supports(pref)
    return
end

# FunctionalDiscreteMeasureData
function expand_measure(
    expr,
    data::FunctionalDiscreteMeasureData{P, B, I},
    write_model::JuMP.AbstractModel
    ) where {P, B, I}
    # get the info
    prefs = parameter_refs(data)
    _prep_generative_supps(prefs, I)
    supps = supports(data)
    coeffs = coefficients(data)
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
function expand_measure(
    expr, 
    data::AbstractMeasureData,
    ::JuMP.AbstractModel
    )
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
function analytic_expansion(
    expr::JuMP.AbstractJuMPScalar,
    data::Union{DiscreteMeasureData{GeneralVariableRef, 1},
                FunctionalDiscreteMeasureData{GeneralVariableRef}},
    write_model::JuMP.AbstractModel # needed for fallback
    )
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
function analytic_expansion(
    expr::JuMP.AbstractJuMPScalar,
    data::Union{DiscreteMeasureData{Vector{GeneralVariableRef}, 2},
                FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}},
    write_model::JuMP.AbstractModel # needed for fallback
    )
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
function analytic_expansion(
    expr, 
    data::AbstractMeasureData,
    write_model::JuMP.AbstractModel
    )
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
point variables and semi-infinite variables are made as needed). Errors if
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
function expand(mref::MeasureRef)
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
[`add_point_variable`](@ref)/[`add_semi_infinite_variable`](@ref) in combination 
with `expand_measure`. `write_model` is the model that the measure variables are 
added to as described in [`expand_measure`](@ref).
"""
function expand_measures end

# MeasureRef
function expand_measures(
    mref::GeneralVariableRef,
    ::Type{MeasureIndex},
    write_model::JuMP.AbstractModel
    )
    if is_analytic(mref)
        return analytic_expansion(measure_function(mref), measure_data(mref),
                                  write_model)
    else
        return expand_measure(measure_function(mref), measure_data(mref),
                              write_model)
    end
end

# NonMeasureRef
function expand_measures(
    vref::GeneralVariableRef,
    ::Type{V},
    ::JuMP.AbstractModel
    ) where {V <: AbstractInfOptIndex}
    return vref
end

# GeneralVariableRef
function expand_measures(
    vref::GeneralVariableRef,
    write_model::JuMP.AbstractModel
    )
    return expand_measures(vref, _index_type(vref), write_model)
end

# Expressions
function expand_measures(
    expr::AbstractInfOptExpr,
    write_model::JuMP.AbstractModel
    )
    return map_expression(v -> expand_measures(v, write_model), expr)
end

# AbstractArray of expressions 
function expand_measures(arr::AbstractArray,
                         write_model::JuMP.AbstractModel
                         )
    return map(e -> expand_measures(e, write_model), arr)
end

# Fallback
function expand_measures(expr, ::JuMP.AbstractModel)
    error("`expand_measures` not defined for expressions of type ",
          "`$(typeof(expr))`.")
end

"""
    expand_all_measures!(model::InfiniteModel)::Nothing

Expand all of the measures used in the objective and/or constraints of `model`.
The objective and constraints are updated accordingly. Note that
variables are added to the model as necessary to accomodate the expansion (i.e.,
point variables and semi-infinite variables are made as needed). Errors if
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
function expand_all_measures!(model::InfiniteModel)
    # expand the objective if it contains measures
    if objective_has_measures(model)
        new_obj = expand_measures(JuMP.objective_function(model), model)
        JuMP.set_objective_function(model, new_obj)
    end
    # expand all of the constraints that contain measures
    for (cindex, object) in model.constraints
        if !isempty(object.measure_indices)
            old_constr = object.constraint
            set = JuMP.moi_set(old_constr)
            # clear the old dependencies
            old_vrefs = _all_function_variables(JuMP.jump_function(old_constr))
            for vref in old_vrefs
                filter!(e -> e != cindex, _constraint_dependencies(vref))
            end
            # expand the expression
            new_func = expand_measures(JuMP.jump_function(old_constr), model)
            vrefs = _all_function_variables(new_func)
            # make the new constraint object
            new_constr = JuMP.build_constraint(error, new_func, set)
            # update the constraint data
            cref = _make_constraint_ref(model, cindex)
            _set_core_constraint_object(cref, new_constr)
            empty!(object.measure_indices)
            _update_var_constr_mapping(vrefs, cref)
        end
    end
    return
end
