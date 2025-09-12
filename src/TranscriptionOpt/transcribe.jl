################################################################################
#                           SUPPORT ITERATION METHODS
################################################################################
"""
    set_parameter_supports(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Collect the infinite parameter supports stored in their respective dictionaries
form `model` and process them into a tuple of vectors where each vector
contains the collected supports of a particular infinite parameter. These support
collections are ordered in accordance with the definition order of the
parameters (i.e., their group integer indices). A support collection assocciated with
an independent will be a `Vector{Float64}` and a support collection associated
with a group of dependent parameters will be a `Vector{Vector{Float64}}`. Note
that each collection vector will include an extra final placeholder element
comprised of `NaN`s for convenience in generating support indices via
[`support_index_iterator`](@ref). This also gathers the associated support labels. 

Before this is all done, `InfiniteOpt.add_generative_supports` is invoked as needed.
"""
function set_parameter_supports(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    # gather the basic information
    prefs = InfiniteOpt.parameter_refs(model)
    data = transcription_data(backend)
    # check and add supports to prefs as needed
    for group in prefs 
        InfiniteOpt.add_generative_supports(first(group))
        if InfiniteOpt.has_internal_supports(first(group))
            data.has_internal_supports = true
        end
    end
    # build and add the support/label tuples
    supps = Tuple(begin 
        supp_dict = InfiniteOpt.core_object(first(group)).supports
        supp_list = collect(keys(supp_dict))
        # append a placeholder NaN support at the end to be used for efficient combinatorics
        push!(supp_list, map(i -> NaN, first(supp_list)))
    end for group in prefs)
    labels = Tuple(begin
        supp_dict = InfiniteOpt.core_object(first(group)).supports
        default = Set{DataType}()
        map(k -> get(supp_dict, k, default), supps[i])
    end for (i, group) in enumerate(prefs))
    data.supports = supps
    data.support_labels = labels
    return
end

################################################################################
#                        VARIABLE INITIALIZATION METHODS
################################################################################
"""
    transcribe_finite_parameters!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Create a transcription variable (i.e., a JuMP Parameter) for each `FiniteParameter`
stored in `model` and add it to `backend`. The variable mapping is
also stored in `TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref).
"""

function transcribe_finite_parameters!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.FiniteParameter)
        fpref = InfiniteOpt.GeneralVariableRef(model, idx)

        # Prepare arguments for building JuMP parameter
        name = object.name
        pValue = object.parameter.value
        
        # Build the JuMP variable and add it to the backend
        pref = JuMP.@variable(backend.model, base_name = name, set = MOI.Parameter(pValue))
        transcription_data(backend).finvar_mappings[fpref] = pref
    end
    return
end
"""
    transcribe_finite_variables!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Create a transcription variable (i.e., a JuMP variable) for each `FiniteVariable`
stored in `model` and add it to `backend`. The variable mapping is
also stored in `TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref).
"""
function transcribe_finite_variables!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.FiniteVariable)
        hvref = InfiniteOpt.GeneralVariableRef(model, idx)
        v = JuMP.ScalarVariable(object.variable.info)
        vref = JuMP.add_variable(backend.model, v, object.name)
        transcription_data(backend).finvar_mappings[hvref] = vref
    end
    return
end

# Return a proper scalar variable info object given on with a start value function
function _format_infinite_info(
    var::InfiniteOpt.InfiniteVariable,
    support::Vector{Float64}
    )
    # generate the start value
    info = var.info
    if var.is_vector_start
        start = info.start(support)
    else
        start = info.start(Tuple(support, var.parameter_refs)...)
    end
    # make the info and return
    return JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                             info.upper_bound, info.has_fix, info.fixed_value,
                             info.has_start, start, info.binary, info.integer)
end

# Set the cutoff for number of infinite parameters to be included in a variable name
const _MaxNumParamsForPrinting = 4

# Make variable name with infinite parameters values directly if possible
function _make_var_name(base_name, param_nums, tuple_supp, var_idx)
    if length(param_nums) <= _MaxNumParamsForPrinting
        return string(base_name, "(", join(tuple_supp, ", "), ")")
    else
        return string(base_name, "[", join(var_idx, ", "), "]")
    end
end

"""
    transcribe_infinite_variables!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Create transcription variables (i.e., JuMP variables) for each `InfiniteVariable`
stored in `model` and add them to `backend`. The variable mappings are
also stored in `TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_infinite_variables!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.InfiniteVariable)
        # get the basic variable information
        var = object.variable
        base_name = object.name
        param_nums = var.parameter_nums
        group_idxs = var.group_int_idxs
        prefs = var.parameter_refs
        # prepare for iterating over its supports
        supp_indices = support_index_iterator(backend, group_idxs)
        dims = size(supp_indices)[group_idxs]
        vrefs = Array{JuMP.VariableRef, length(dims)}(undef, dims...)
        supp_type = typeof(Tuple(ones(length(prefs)), prefs))
        supps = Array{supp_type, length(dims)}(undef, dims...)
        lookup_dict = sizehint!(Dict{Vector{Float64}, JuMP.VariableRef}(), length(vrefs))
        # create a variable for each support
        for i in supp_indices
            supp = index_to_support(backend, i)[param_nums]
            info = _format_infinite_info(var, supp)
            var_idx = i.I[group_idxs]
            tuple_supp = Tuple(supp, prefs)
            v_name = _make_var_name(base_name, param_nums, tuple_supp, var_idx)
            v = JuMP.ScalarVariable(info)
            jump_vref = JuMP.add_variable(backend.model, v, v_name)
            @inbounds vrefs[var_idx...] = jump_vref
            lookup_dict[supp] = jump_vref
            @inbounds supps[var_idx...] = tuple_supp
        end
        # save the transcription information
        ivref = InfiniteOpt.GeneralVariableRef(model, idx)
        data = transcription_data(backend)
        data.infvar_lookup[ivref] = lookup_dict
        data.infvar_mappings[ivref] = vrefs
        data.infvar_supports[ivref] = supps
    end
    return
end

# Return a proper scalar variable info object given on with a start value function
function _format_derivative_info(d::InfiniteOpt.Derivative, support::Vector{Float64})
    # generate the start value
    info = d.info
    if d.is_vector_start
        start = info.start(support)
    else
        prefs = InfiniteOpt.raw_parameter_refs(d.variable_ref)
        start = info.start(Tuple(support, prefs)...)
    end
    # make the info and return
    return JuMP.VariableInfo(info.has_lb, info.lower_bound, info.has_ub,
                             info.upper_bound, info.has_fix, info.fixed_value,
                             info.has_start, start, info.binary, info.integer)
end

function _transcribe_derivative_variable(dref, d, backend)
    base_name = InfiniteOpt.variable_string(MIME("text/plain"), dispatch_variable_ref(dref))
    param_nums = InfiniteOpt._parameter_numbers(d.variable_ref)
    group_idxs = InfiniteOpt.parameter_group_int_indices(d.variable_ref)
    prefs = InfiniteOpt.raw_parameter_refs(dref)
    # prepare for iterating over its supports
    supp_indices = support_index_iterator(backend, group_idxs)
    dims = size(supp_indices)[group_idxs]
    vrefs = Array{JuMP.VariableRef, length(dims)}(undef, dims...)
    supp_type = typeof(Tuple(ones(length(prefs)), prefs))
    supps = Array{supp_type, length(dims)}(undef, dims...)
    lookup_dict = sizehint!(Dict{Vector{Float64}, JuMP.VariableRef}(), length(vrefs))
    # create a variable for each support
    for i in supp_indices
        supp = index_to_support(backend, i)[param_nums]
        info = _format_derivative_info(d, supp)
        var_idx = i.I[group_idxs]
        tuple_supp = Tuple(supp, prefs)
        d_name = _make_var_name(base_name, param_nums, tuple_supp, var_idx)
        d_var = JuMP.ScalarVariable(info)
        jump_vref = JuMP.add_variable(backend.model, d_var, d_name)
        @inbounds vrefs[var_idx...] = jump_vref
        lookup_dict[supp] = jump_vref
        @inbounds supps[var_idx...] = tuple_supp
    end
    # save the transcription information
    data = transcription_data(backend)
    data.infvar_lookup[dref] = lookup_dict
    data.infvar_mappings[dref] = vrefs
    data.infvar_supports[dref] = supps
    return
end

"""
    transcribe_derivative_variables!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Create transcription variables (i.e., JuMP variables) for each `Derivative`
stored in `model` and add them to `backend`. The variable mappings are
also stored in `TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. The 
futher derivative evaluation constraints are added when 
`transcribe_derivative_evaluations!` is invoked. Note that 
`TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_derivative_variables!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.Derivative)
        # get the basic derivative information
        dref = InfiniteOpt.GeneralVariableRef(model, idx)
        d = object.variable
        method = InfiniteOpt.derivative_method(dref)
        # if needed process lower order derivatives
        if !InfiniteOpt.allows_high_order_derivatives(method) && d.order > 1
            for o in d.order-1:-1:1
                if !haskey(model.deriv_lookup, (d.variable_ref, d.parameter_ref, o))
                    info = JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false, 
                                             s -> NaN, false, false)
                    new_d = InfiniteOpt.Derivative(info, true, d.variable_ref, d.parameter_ref, o)
                    new_dref = InfiniteOpt.add_derivative(model, new_d)
                    _transcribe_derivative_variable(new_dref, d, backend)
                end
            end
        end
        # process the derivative
        _transcribe_derivative_variable(dref, d, backend)
    end
    return
end

"""
    transcribe_semi_infinite_variables!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Map each `SemiInfiniteVariable` in `model` to transcription variables stored in
`backend`. The variable mappings are also stored in
`TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that [`transcribe_infinite_variables!`](@ref)
must be called first. Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_semi_infinite_variables!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.SemiInfiniteVariable)
        # get the basic variable information
        var = object.variable
        rvref = InfiniteOpt.GeneralVariableRef(model, idx)
        # setup the mappings
        ivref = var.infinite_variable_ref
        if InfiniteOpt._index_type(ivref) != InfiniteOpt.ParameterFunctionIndex
            param_nums = var.parameter_nums
            ivref_param_nums = InfiniteOpt._parameter_numbers(ivref)
            eval_supps = var.eval_supports
            group_idxs = var.group_int_idxs
            prefs = InfiniteOpt.raw_parameter_refs(var)
            # prepare for iterating over its supports
            supp_indices = support_index_iterator(backend, group_idxs)
            dims = size(supp_indices)[group_idxs]
            vrefs = Array{JuMP.VariableRef, length(dims)}(undef, dims...)
            supp_type = typeof(Tuple(ones(length(prefs)), prefs))
            supps = Array{supp_type, length(dims)}(undef, dims...)
            lookup_dict = sizehint!(Dict{Vector{Float64}, JuMP.VariableRef}(), length(vrefs))
            valid_idxs = ones(Bool, dims...)
            # map a variable for each support
            for i in supp_indices
                raw_supp = index_to_support(backend, i)
                var_idx = i.I[group_idxs]
                # ensure this support is valid with the reduced restriction
                if any(!isnan(raw_supp[ivref_param_nums[k]]) && raw_supp[ivref_param_nums[k]] != v for (k, v) in eval_supps)
                    valid_idxs[var_idx...] = false
                    continue
                end
                # map to the current transcription variable
                supp = raw_supp[param_nums]
                ivref_supp = [haskey(eval_supps, j) ? eval_supps[j] : raw_supp[k] 
                            for (j, k) in enumerate(ivref_param_nums)]
                jump_vref = lookup_by_support(ivref, backend, ivref_supp)
                @inbounds vrefs[var_idx...] = jump_vref
                lookup_dict[supp] = jump_vref
                @inbounds supps[var_idx...] = Tuple(supp, prefs)
            end
            # truncate vrefs if any supports were skipped because of dependent parameter supps and save
            data = transcription_data(backend)
            if !all(valid_idxs)
                data.infvar_mappings[rvref] = vrefs[valid_idxs]
                data.infvar_supports[rvref] = supps[valid_idxs]
                data.valid_indices[rvref] = valid_idxs
            else
                data.infvar_mappings[rvref] = vrefs
                data.infvar_supports[rvref] = supps
            end
            data.infvar_lookup[rvref] = lookup_dict
        end
    end
    return
end

# Override the info of the jump variable with the point variable's if it is different
function _update_point_info(
    gvref::InfiniteOpt.GeneralVariableRef,
    vref::JuMP.VariableRef
    )
    pvref = InfiniteOpt.dispatch_variable_ref(gvref)
    if JuMP.has_lower_bound(pvref)
        if JuMP.is_fixed(vref)
            JuMP.unfix(vref)
        end
        JuMP.set_lower_bound(vref, JuMP.lower_bound(pvref))
    end
    if JuMP.has_upper_bound(pvref)
        if JuMP.is_fixed(vref)
            JuMP.unfix(vref)
        end
        JuMP.set_upper_bound(vref, JuMP.upper_bound(pvref))
    end
    if JuMP.is_fixed(pvref)
        JuMP.fix(vref, JuMP.fix_value(pvref), force = true)
    end
    if JuMP.is_binary(pvref) && !JuMP.is_binary(vref)
        if JuMP.is_integer(vref)
            JuMP.unset_integer(vref)
        end
        JuMP.set_binary(vref)
    elseif JuMP.is_integer(pvref) && !JuMP.is_integer(vref)
        if JuMP.is_binary(vref)
            JuMP.unset_binary(vref)
        end
        JuMP.set_integer(vref)
    end
    if JuMP.start_value(pvref) != JuMP.start_value(vref)
        JuMP.set_start_value(vref, JuMP.start_value(pvref))
    end
    return
end

"""
    transcribe_point_variables!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Map each `PointVariable` in `model` to a transcription variable stored in
`backend`. The variable mapping is also stored in
`TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref). Note that
[`transcribe_infinite_variables!`](@ref) must be called first and that the
info constraints associated with the transcription variable will be updated in
accordance with the point variable.
"""
function transcribe_point_variables!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.PointVariable)
        # get the basic variable information
        var = object.variable
        ivref = var.infinite_variable_ref
        supp = var.parameter_values
        # find the corresponding variable record the mapping
        vref = lookup_by_support(ivref, backend, supp)
        pvref = InfiniteOpt.GeneralVariableRef(model, idx)
        transcription_data(backend).finvar_mappings[pvref] = vref
        # update the info constraints as needed
        _update_point_info(pvref, vref)
    end
    return
end

################################################################################
#                       TRANSCRIPTION EXPRESSION METHODS
################################################################################
"""
    transcription_expression(
        expr,
        backend::TranscriptionBackend,
        support::Vector{Float64}
        )

Given the `expr` from an `InfiniteModel`, form its transcripted version in
accordance with the variable mappings available in `backend` defined at
`support`. This should only be used once all variables and measures have been
transcribed (e.g., via [`transcribe_finite_variables!`](@ref)).
"""
function transcription_expression(expr, backend::TranscriptionBackend, support)
    error("Unsupported expression type `$(typeof(expr))` for automated " *
          "transcription.")
end

# GeneralVariableRef
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::TranscriptionBackend,
    support::Vector{Float64}
    )
    idx_type = InfiniteOpt._index_type(vref)
    return transcription_expression(vref, idx_type, backend, support)
end

# Infinite variables, infinite parameter functions, and measures
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    backend::TranscriptionBackend,
    support::Vector{Float64}
    ) where {V <: Union{InfVarIndex, InfiniteOpt.ParameterFunctionIndex, InfiniteOpt.MeasureIndex}}
    param_nums = InfiniteOpt._parameter_numbers(vref)
    return lookup_by_support(vref, index_type, backend, support[param_nums])
end

# Semi-Infinite variables
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.SemiInfiniteVariableIndex},
    backend::TranscriptionBackend,
    support::Vector{Float64}
    )
    ivref = InfiniteOpt.infinite_variable_ref(vref)
    if InfiniteOpt._index_type(ivref) == InfiniteOpt.ParameterFunctionIndex 
        prefs = InfiniteOpt.raw_parameter_refs(ivref)
        param_nums = InfiniteOpt._parameter_numbers(ivref)
        func = InfiniteOpt.raw_function(ivref)
        return func(Tuple(support[param_nums], prefs)...)
    else 
        param_nums = InfiniteOpt._parameter_numbers(vref)
        return lookup_by_support(vref, index_type, backend, support[param_nums])
    end
end

# Point variables, finite variables and finite parameters
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    backend::TranscriptionBackend,
    support::Vector{Float64}
    ) where {V <: Union{FinVarIndex, InfiniteOpt.FiniteParameterIndex}}
    return lookup_by_support(vref, index_type, backend, support)
end

# Infinite parameters
function transcription_expression(
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    backend::TranscriptionBackend,
    support::Vector{Float64}
    ) where {V <: InfiniteOpt.InfiniteParameterIndex}
    param_num = InfiniteOpt._parameter_number(vref)
    return support[param_num]
end


# AffExpr and QuadExpr and NonlinearExpr
function transcription_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::TranscriptionBackend,
    support::Vector{Float64}
    )
    return InfiniteOpt.map_expression(
        v -> transcription_expression(v, backend, support), 
        expr)
end

# Real Number 
function transcription_expression(
    num::Real,
    backend::TranscriptionBackend,
    support::Vector{Float64}
    )
    return zero(JuMP.AffExpr) + num
end

################################################################################
#                         MEASURE TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_measures!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

For each `Measure` in `model` expand it via `InfiniteOpt.expand_measure` or
`analytic_expansion` as appropriate and transcribe the expanded expression via
[`transcription_expression`](@ref). Then store the measure to transcripted
expression mappings in `TranscriptionData.measure_mappings` and
`TranscriptionData.measure_lookup` to enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.measure_support_labels` is also populated.
"""
function transcribe_measures!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.Measure)
        # get the basic information
        meas = object.measure
        group_idxs = meas.group_int_idxs
        mref = InfiniteOpt.GeneralVariableRef(model, idx)
        prefs = InfiniteOpt.raw_parameter_refs(mref)
        # expand the measure
        if meas.constant_func
            new_expr = InfiniteOpt.analytic_expansion(meas.func, meas.data, backend)
        else
            new_expr = InfiniteOpt.expand_measure(meas.func, meas.data, backend)
        end
        # prepare to transcribe over the supports
        supp_indices = support_index_iterator(backend, group_idxs)
        dims = size(supp_indices)[group_idxs]
        exprs = Array{JuMP.AbstractJuMPScalar, length(dims)}(undef, dims...)
        supp_type = typeof(Tuple(ones(length(prefs)), prefs))
        supps = Array{supp_type, length(dims)}(undef, dims...)
        lookup_dict = sizehint!(Dict{Vector{Float64}, Int}(), length(exprs))
        # map a variable for each support
        for (lin_idx, i) in enumerate(supp_indices)
            raw_supp = index_to_support(backend, i)
            expr_idx = i.I[group_idxs]
            @inbounds exprs[expr_idx...] = transcription_expression(new_expr, backend, raw_supp)
            supp = raw_supp[meas.parameter_nums]
            lookup_dict[supp] = lin_idx
            @inbounds supps[expr_idx...] = Tuple(supp, prefs)
        end
        # save the transcription information
        data = transcription_data(backend)
        data.measure_lookup[mref] = lookup_dict
        data.measure_mappings[mref] = exprs
        data.measure_supports[mref] = supps
    end
    return
end

################################################################################
#                        OBJECTIVE TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_objective!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Form the transcripted version of the objective stored in `model` and add it
to `backend`. Note that all the variables and measures in `model` must
by transcripted first (e.g., via [`transcribe_infinite_variables!`](@ref)).
"""
function transcribe_objective!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    expr = JuMP.objective_function(model)
    sense = JuMP.objective_sense(model)
    trans_expr = transcription_expression(expr, backend, Float64[])
    JuMP.set_objective(backend.model, sense, trans_expr)
    return
end

################################################################################
#                        CONSTRAINT TRANSCRIPTION METHODS
################################################################################
## Given a variable and its set from an info constraint, get the transcribed version
# Lower bound constraint
function _get_info_constr_from_var(
    backend::TranscriptionBackend,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.GreaterThan,
    support::Vector{Float64}
    )
    trans_vref = transcription_expression(vref, backend, support)
    return JuMP.has_lower_bound(trans_vref) ? JuMP.LowerBoundRef(trans_vref) : nothing
end

# Upper bound constraint
function _get_info_constr_from_var(
    backend::TranscriptionBackend,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.LessThan,
    support::Vector{Float64}
    )
    trans_vref = transcription_expression(vref, backend, support)
    return JuMP.has_upper_bound(trans_vref) ? JuMP.UpperBoundRef(trans_vref) : nothing
end

# Fix constraint
function _get_info_constr_from_var(
    backend::TranscriptionBackend,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.EqualTo,
    support::Vector{Float64}
    )
    trans_vref = transcription_expression(vref, backend, support)
    return JuMP.is_fixed(trans_vref) ? JuMP.FixRef(trans_vref) : nothing
end

# Binary constraint
function _get_info_constr_from_var(
    backend::TranscriptionBackend,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.ZeroOne,
    support::Vector{Float64}
    )
    trans_vref = transcription_expression(vref, backend, support)
    return JuMP.is_binary(trans_vref) ? JuMP.BinaryRef(trans_vref) : nothing
end

# Integer constraint
function _get_info_constr_from_var(
    backend::TranscriptionBackend,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.Integer,
    support::Vector{Float64}
    )
    trans_vref = transcription_expression(vref, backend, support)
    return JuMP.is_integer(trans_vref) ? JuMP.IntegerRef(trans_vref) : nothing
end

# Determine if a given raw support satisfies constraint domain restrictions
function _support_in_restrictions(
    support::Vector{Float64},
    indices::Vector{Int},
    domains::Vector{InfiniteOpt.IntervalDomain}
    )
    for i in eachindex(indices)
        s = support[indices[i]]
        if !isnan(s) && (s < JuMP.lower_bound(domains[i]) || 
            s > JuMP.upper_bound(domains[i]))
            return false
        end
    end
    return true
end

## Process constraint objects at a particular support value and make the new 
## transcribed version
# JuMP.ScalarConstraint
function _process_constraint(
    backend::TranscriptionBackend, 
    constr::JuMP.ScalarConstraint, 
    func::JuMP.AbstractJuMPScalar, 
    set::MOI.AbstractScalarSet, 
    raw_supp::Vector{Float64}, 
    name::String
    )
    new_func = transcription_expression(func, backend, raw_supp)
    trans_constr = JuMP.build_constraint(error, new_func, set)
    return JuMP.add_constraint(backend.model, trans_constr, name)
end

# JuMP.VectorConstraint
function _process_constraint(
    backend::TranscriptionBackend, 
    constr::JuMP.VectorConstraint, 
    func::Vector{<:JuMP.AbstractJuMPScalar}, 
    set::MOI.AbstractVectorSet, 
    raw_supp::Vector{Float64}, 
    name::String
    )
    new_func = map(f -> transcription_expression(f, backend, raw_supp), func)
    shape = JuMP.shape(constr)
    shaped_func = JuMP.reshape_vector(new_func, shape)
    shaped_set = JuMP.reshape_set(set, shape)
    trans_constr = JuMP.build_constraint(error, shaped_func, shaped_set)
    return JuMP.add_constraint(backend.model, trans_constr, name)
end

# Fallback 
function _process_constraint(
    backend::TranscriptionBackend, 
    constr, 
    func, 
    set, 
    raw_supp, 
    name
    )
    error("Transcribing constraints of type `$(typeof(constr))` is not ", 
          "currently supported.")
end

"""
    transcribe_constraints!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel
        )::Nothing

For each constraint in `model` form its transcripted version(s) and add them
to `backend`. The mappings are stored in `TranscriptionData.constr_mappings`
and the associated supports are stored in `TranscriptionData.constr_supports`
to enable [`transcription_constraint`](@ref) and `InfiniteOpt.constraint_supports`.
Note that variable info constraints are simply mapped to the existing info
constraints already generated along with the transcription variables. Note that
the variables and measures must all first be transcripted (e.g., via
[`transcribe_measures!`](@ref)). Note that 
`TranscriptionData.constr_support_labels` is also populated.
"""
function transcribe_constraints!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel
    )
    param_supps = parameter_supports(backend)
    for (idx, object) in model.constraints
        # get the basic information
        constr = object.constraint
        func = JuMP.jump_function(constr)
        set = JuMP.moi_set(constr)
        group_idxs = object.group_int_idxs
        cref = InfiniteOpt.InfOptConstraintRef(model, idx)
        # prepare the iteration helpers
        supp_indices = support_index_iterator(backend, group_idxs)
        dims = size(supp_indices)[group_idxs]
        crefs = Array{JuMP.ConstraintRef, length(dims)}(undef, dims...)
        supps = Array{Tuple, length(dims)}(undef, dims...)
        valid_idxs = ones(Bool, dims...)
        # iterate over the support indices for the info constraints
        if object.is_info_constraint
            for i in supp_indices
                raw_supp = index_to_support(backend, i)
                info_ref = _get_info_constr_from_var(backend, func, set,
                                                     raw_supp)
                # not all supports may be defined if overwritten by a point variable
                con_idx = i.I[group_idxs]
                if !isnothing(info_ref)
                    @inbounds crefs[con_idx...] = info_ref
                    @inbounds supps[con_idx...] = Tuple(param_supps[j][i[j]]
                                                        for j in group_idxs)
                else
                    valid_idxs[con_idx...] = false
                end
            end
        # iterate over the supports for regular constraints
        else
            # get basic setup information
            restrictions = InfiniteOpt.domain_restrictions(cref)
            prefs = collect(keys(restrictions))
            restrict_indices = map(p -> InfiniteOpt._parameter_number(p), prefs)
            restrict_domains = map(p -> restrictions[p], prefs)
            name = object.name
            for i in supp_indices
                raw_supp = index_to_support(backend, i)
                # ensure the support satisfies parameter bounds and then add it
                con_idx = i.I[group_idxs]
                if _support_in_restrictions(raw_supp, restrict_indices, 
                                            restrict_domains)
                    new_name = if isempty(name)
                        ""
                    elseif isempty(group_idxs)
                        name
                    else
                        string(name, "[", join(con_idx, ", "), "]")
                    end
                    new_cref = _process_constraint(backend, constr, func, 
                                                   set, raw_supp, new_name) 
                    @inbounds crefs[con_idx...] = new_cref
                    @inbounds supps[con_idx...] = Tuple(param_supps[j][i[j]]
                                                        for j in group_idxs)
                else
                    valid_idxs[con_idx...] = false
                end
            end
        end
        # truncate the arrays in case not all the supports satisfied the bounds
        # and save
        data = transcription_data(backend)
        if !all(valid_idxs)
            data.constr_mappings[cref] = crefs[valid_idxs]
            data.constr_supports[cref] = supps[valid_idxs]
            data.valid_indices[cref] = valid_idxs
        else
            data.constr_mappings[cref] = crefs
            data.constr_supports[cref] = supps
        end
    end
    return
end

################################################################################
#                    DERIVATIVE CONSTRAINT TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_derivative_evaluations!(
        backend::TranscriptionBackend, 
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Generate the auxiliary derivative evaluation equations and transcribe them 
appropriately for all the derivatives in `model`. These are in turn added to 
`backend`. Note that no mapping information is recorded since the InfiniteModel 
won't have any constraints that correspond to these equations. Also note that
the variables and measures must all first be transcripted (e.g., via
[`transcribe_derivative_variables!`](@ref)).
"""
function transcribe_derivative_evaluations!(
    backend::TranscriptionBackend, 
    model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(model, InfiniteOpt.Derivative)
        # get the basic variable information
        dref = InfiniteOpt.GeneralVariableRef(model, idx)
        pref = dispatch_variable_ref(object.variable.parameter_ref)
        method = InfiniteOpt.derivative_method(pref)
        order = object.variable.order
        # generate the transcription constraints as needed
        if !InfiniteOpt.has_derivative_constraints(dref)
            # generate the evaluation expressions
            vref = object.variable.variable_ref
            if !InfiniteOpt.allows_high_order_derivatives(method) && order > 1
                d_idx = model.deriv_lookup[vref, object.variable.parameter_ref, order - 1]
                vref = InfiniteOpt.GeneralVariableRef(model, d_idx)
            end
            exprs = InfiniteOpt.evaluate_derivative(dref, vref, method, backend)
            # prepare the iteration helpers
            param_group_int_idx = InfiniteOpt.parameter_group_int_index(pref)
            group_int_idxs = filter(!isequal(param_group_int_idx), InfiniteOpt.parameter_group_int_indices(dref))
            supp_indices = support_index_iterator(backend, group_int_idxs)
            # transcribe the constraints 
            set = MOI.EqualTo(0.0)
            for i in supp_indices 
                raw_supp = index_to_support(backend, i)
                for expr in exprs
                    new_expr = transcription_expression(expr, backend, raw_supp)
                    trans_constr = JuMP.build_constraint(error, new_expr, set)
                    JuMP.add_constraint(backend.model, trans_constr) # TODO maybe add name?
                end 
            end
        end
    end
    return
end

"""
    transcribe_variable_collocation_restictions!(
        backend::TranscriptionBackend, 
        model::InfiniteOpt.InfiniteModel
        )::Nothing

Add constraints to `backend` that make infinite variables constant over collocation 
points following the calls made to [`InfiniteOpt.constant_over_collocation`](@ref). Note that 
[`set_parameter_supports`](@ref) and [`transcribe_infinite_variables!`](@ref) must be called first.
"""
function transcribe_variable_collocation_restictions!(
    backend::TranscriptionBackend, 
    model::InfiniteOpt.InfiniteModel
    )
    data = transcription_data(backend)
    set = MOI.EqualTo(0.0)
    for (pidx, vidxs) in model.piecewise_vars
        pref = InfiniteOpt.GeneralVariableRef(model, pidx)
        if !InfiniteOpt.has_generative_supports(pref)
            continue
        end
        group_int_idx = InfiniteOpt.parameter_group_int_index(pref)
        supps = reverse!(data.supports[group_int_idx][1:end-1])
        labels = reverse!(data.support_labels[group_int_idx][1:end-1])
        @assert any(l -> l <: InfiniteOpt.PublicLabel, first(labels))
        v_manip = GeneralVariableRef(model, -1, IndependentParameterIndex) # placeholder
        for vidx in vidxs
            vref = InfiniteOpt.GeneralVariableRef(model, vidx)
            group_int_idxs = filter(!isequal(group_int_idx), InfiniteOpt.parameter_group_int_indices(vref))
            supp_indices = support_index_iterator(backend, group_int_idxs)
            for (s, ls) in zip(supps, labels)
                if any(l -> l <: InfiniteOpt.PublicLabel, ls)
                    v_manip = InfiniteOpt.make_reduced_expr(vref, pref, s, backend)
                else
                    inf_expr = v_manip - InfiniteOpt.make_reduced_expr(vref, pref, s, backend)
                    for i in supp_indices 
                        raw_supp = index_to_support(backend, i)
                        new_expr = transcription_expression(inf_expr, backend, raw_supp)
                        trans_constr = JuMP.build_constraint(error, new_expr, set)
                        JuMP.add_constraint(backend.model, trans_constr)
                    end
                end
            end
        end
    end
    return
end

################################################################################
#                      INFINITEMODEL TRANSCRIPTION METHODS
################################################################################
"""
    build_transcription_backend!(
        backend::TranscriptionBackend,
        model::InfiniteOpt.InfiniteModel;
        [check_support_dims::Bool = true]
        )::Nothing

Given an empty `backend` build it using the information stored in `model`.
This is intended for a `TranscriptionModel` that serves as a internal transformation backend
of `model`. This detail is important to correctly enable internally generated
semi-infinite variables during the transcription process such that `model` is not
modified. Note that this will add supports to `model` via
[`InfiniteOpt.fill_in_supports!`](@ref) for infinite parameters that contain
no supports. Also a warning is thrown when the transcription backend contains
more than 15,000 support points to alert users when they may naively have
a few independent supports whose product quickly yields a very large grid.
For example having 3 independent parameters with 100 supports each would result
in 1,000,000 supports if all three are together in at least 1 constraint. This 
behavior can be overcome using dependent parameters. The warning can be turned off 
via `check_support_dims = false`.
"""
function build_transcription_backend!(
    backend::TranscriptionBackend,
    model::InfiniteOpt.InfiniteModel;
    check_support_dims::Bool = true
    )
    # ensure there are supports to add and add them to the trans model
    InfiniteOpt.fill_in_supports!(model, modify = false)
    set_parameter_supports(backend, model)
    # check that there isn't a crazy amount of supports from taking the product
    supps = parameter_supports(backend)
    num_supps = isempty(supps) ? 0 : prod(length, supps)
    if check_support_dims && length(supps) > 1 && num_supps > 15000 # NOTE this is an arbitrary cutoff
        @warn("Due to necessarily considering the combinatorics of independent " *
              "parameter supports, the model will be transcripted over up to $(num_supps) " *
              "supports (if any constraint uses all the infinite parameters) "  *
              "and thus naive solution of the discretized problem may be slow. " * 
              "This warning can be turned off via `check_support_dims = false`.")
    end
    # add nonlinear operators as needed 
    InfiniteOpt.add_operators_to_jump(backend.model, model)
    # define the variables
    transcribe_finite_parameters!(backend, model)
    transcribe_finite_variables!(backend, model)
    transcribe_infinite_variables!(backend, model)
    transcribe_derivative_variables!(backend, model)
    transcribe_semi_infinite_variables!(backend, model)
    transcribe_point_variables!(backend, model)
    transcribe_measures!(backend, model)
    # define the objective
    transcribe_objective!(backend, model)
    # define the constraints
    transcribe_constraints!(backend, model)
    # define the derivative evaluation constraints
    transcribe_derivative_evaluations!(backend, model)
    # define constraints for variables that are constant over collocation points
    transcribe_variable_collocation_restictions!(backend, model)
    return
end

"""
    InfiniteOpt.build_transformation_backend!(
        model::InfiniteOpt.InfiniteModel,
        backend::TranscriptionBackend;
        check_support_dims::Bool = true
        )::Nothing

Build `backend` and set it as the transformation backend to `model`.
Ths clears out the existing `backend` and rebuilds it. Optionally, 
the dimension check to through a warning if there is potentially 
a very large number of supports can be turned off via 
`check_support_dims = false`.
"""
function InfiniteOpt.build_transformation_backend!(
    model::InfiniteOpt.InfiniteModel,
    backend::TranscriptionBackend;
    check_support_dims::Bool = true,
    extra_kwargs...
    )
    # throw error for extra keywords
    for (kw, _) in extra_kwargs
        error("Unrecognized keyword argument `$kw` for building transcription backends.")
    end
    # clear the the backend model contents
    empty!(backend)
    backend.model.operator_counter = 0
    # build the transcription backend based on model
    build_transcription_backend!(
        backend,
        model, 
        check_support_dims = check_support_dims
        )
    return
end

# TODO make TranscriptionModel(::InfiniteModel) --> have to deal with accessing semi_infinite vars
