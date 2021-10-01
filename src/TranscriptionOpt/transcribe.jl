################################################################################
#                           SUPPORT ITERATION METHODS
################################################################################
## Form a placeholder parameter reference given the object index
# IndependentParameterIndex
function _temp_parameter_ref(
    model::InfiniteOpt.InfiniteModel,
    idx::InfiniteOpt.IndependentParameterIndex
    )::InfiniteOpt.IndependentParameterRef
    return InfiniteOpt.IndependentParameterRef(model, idx)
end

# DependentParametersIndex
function _temp_parameter_ref(
    model::InfiniteOpt.InfiniteModel,
    idx::InfiniteOpt.DependentParametersIndex
    )::InfiniteOpt.DependentParameterRef
    idx = InfiniteOpt.DependentParameterIndex(idx, 1)
    return InfiniteOpt.DependentParameterRef(model, idx)
end

# Return the collected supports of an infinite parameter
function _collected_supports(
    pref::Union{InfiniteOpt.IndependentParameterRef, InfiniteOpt.DependentParameterRef}
    )::Vector
    supp_dict = InfiniteOpt._parameter_supports(pref)
    supp_list = collect(keys(supp_dict))
    # append a placeholder NaN support at the end to be used for efficient combinatorics
    return push!(supp_list, map(i -> NaN, first(supp_list)))
end

# Return the collected support labels of an infinite parameter
function _collected_support_labels(
    pref::Union{InfiniteOpt.IndependentParameterRef, InfiniteOpt.DependentParameterRef},
    supports::Vector
    )::Vector{Set{DataType}}
    supp_dict = InfiniteOpt._parameter_supports(pref)
    default = Set{DataType}()
    return map(k -> get(supp_dict, k, default), supports)
end

"""
    set_parameter_supports(trans_model::JuMP.Model,
                           inf_model::InfiniteOpt.InfiniteModel)::Nothing

Collect the infinite parameter supports stored in their respective dictionaries
form `inf_model` and process them into a tuple of vectors where each vector
contains the collected supports of a particular infinite parameter. These support
collections are ordered in accordance with the definition order of the
parameters (i.e., their object numbers). A support collection assocciated with
an independent will be a `Vector{Float64}` and a support collection associated
with a group of dependent parameters will be a `Vector{Vector{Float64}}`. Note
that each collection vector will include an extra final placeholder element
comprised of `NaN`s for convenience in generating support indices via
[`support_index_iterator`](@ref). This also gathers the associated support labels. 

Before this is all done, `InfiniteOpt.add_generative_supports` is invoked as needed.
"""
function set_parameter_supports(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    # gather the basic information
    param_indices = InfiniteOpt._param_object_indices(inf_model)
    prefs = map(idx -> _temp_parameter_ref(inf_model, idx), param_indices)
    data = transcription_data(trans_model)
    # check and add supports to prefs as needed
    for pref in prefs 
        InfiniteOpt.add_generative_supports(pref)
        if InfiniteOpt.has_internal_supports(pref)
            data.has_internal_supports = true
        end 
    end
    # build and add the support/label tuples
    supps = Tuple(_collected_supports(pref) for pref in prefs)
    labels = Tuple(_collected_support_labels(pref, supps[i]) 
                   for (i, pref) in enumerate(prefs))
    data.supports = supps
    data.support_labels = labels
    return
end

################################################################################
#                        VARIABLE INITIALIZATION METHODS
################################################################################
"""
    transcribe_finite_variables!(trans_model::JuMP.Model,
                               inf_model::InfiniteOpt.InfiniteModel)::Nothing

Create a transcription variable (i.e., a JuMP variable) for each `FiniteVariable`
stored in `inf_model` and add it to `trans_model`. The variable mapping is
also stored in `TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref).
"""
function transcribe_finite_variables!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.FiniteVariable)
        hvref = InfiniteOpt._make_variable_ref(inf_model, idx)
        vref = JuMP.add_variable(trans_model,
                                 JuMP.ScalarVariable(object.variable.info),
                                 object.name)
        transcription_data(trans_model).finvar_mappings[hvref] = vref
    end
    return
end

# Return a proper scalar variable info object given on with a start value function
function _format_infinite_info(
    var::InfiniteOpt.InfiniteVariable,
    support::Vector{Float64}
    )::JuMP.VariableInfo
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

"""
    transcribe_infinite_variables!(trans_model::JuMP.Model,
                                   inf_model::InfiniteOpt.InfiniteModel)::Nothing

Create transcription variables (i.e., JuMP variables) for each `InfiniteVariable`
stored in `inf_model` and add them to `trans_model`. The variable mappings are
also stored in `TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_infinite_variables!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.InfiniteVariable)
        # get the basic variable information
        var = object.variable
        base_name = object.name
        param_nums = var.parameter_nums
        # prepare for iterating over its supports
        supp_indices = support_index_iterator(trans_model, var.object_nums)
        vrefs = Vector{JuMP.VariableRef}(undef, length(supp_indices))
        labels = Vector{Set{DataType}}(undef, length(supp_indices))
        lookup_dict = Dict{Vector{Float64}, Int}()
        # create a variable for each support
        for (counter, i) in enumerate(supp_indices)
            raw_supp = index_to_support(trans_model, i)
            supp = raw_supp[param_nums]
            info = _format_infinite_info(var, supp)
            var_name = string(base_name, "(support: ", counter, ")")
            @inbounds vrefs[counter] = JuMP.add_variable(trans_model,
                                                         JuMP.ScalarVariable(info),
                                                         var_name)
            lookup_dict[supp] = counter
            @inbounds labels[counter] = index_to_labels(trans_model, i)
        end
        # save the transcription information
        ivref = InfiniteOpt._make_variable_ref(inf_model, idx)
        data = transcription_data(trans_model)
        data.infvar_lookup[ivref] = lookup_dict
        data.infvar_mappings[ivref] = vrefs
        data.infvar_support_labels[ivref] = labels
    end
    return
end

# Return a proper scalar variable info object given on with a start value function
function _format_derivative_info(d::InfiniteOpt.Derivative,
    support::Vector{Float64}
    )::JuMP.VariableInfo
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

"""
    transcribe_derivative_variables!(trans_model::JuMP.Model,
                                     inf_model::InfiniteOpt.InfiniteModel)::Nothing

Create transcription variables (i.e., JuMP variables) for each `Derivative`
stored in `inf_model` and add them to `trans_model`. The variable mappings are
also stored in `TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. The 
futher derivative evaluation constraints are added when 
`transcribe_derivative_evaluations!` is invoked. Note that 
`TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_derivative_variables!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.Derivative)
        # get the basic variable information
        dref = InfiniteOpt._make_variable_ref(inf_model, idx)
        d = object.variable
        base_name = InfiniteOpt.variable_string(JuMP.REPLMode, dispatch_variable_ref(dref))
        param_nums = InfiniteOpt._parameter_numbers(d.variable_ref)
        obj_nums = InfiniteOpt._object_numbers(d.variable_ref)
        # prepare for iterating over its supports
        supp_indices = support_index_iterator(trans_model, obj_nums)
        vrefs = Vector{JuMP.VariableRef}(undef, length(supp_indices))
        labels = Vector{Set{DataType}}(undef, length(supp_indices))
        lookup_dict = Dict{Vector{Float64}, Int}()
        # create a variable for each support
        for (counter, i) in enumerate(supp_indices)
            raw_supp = index_to_support(trans_model, i)
            supp = raw_supp[param_nums]
            info = _format_derivative_info(d, supp)
            deriv_name = string(base_name, "(support: ", counter, ")")
            @inbounds vrefs[counter] = JuMP.add_variable(trans_model,
                                                         JuMP.ScalarVariable(info),
                                                         deriv_name)
            lookup_dict[supp] = counter
            @inbounds labels[counter] = index_to_labels(trans_model, i)
        end
        # save the transcription information
        data = transcription_data(trans_model)
        data.infvar_lookup[dref] = lookup_dict
        data.infvar_mappings[dref] = vrefs
        data.infvar_support_labels[dref] = labels
    end
    return
end

# Setup the mapping for a given semi_infinite variable
function _set_semi_infinite_variable_mapping(
    trans_model::JuMP.Model,
    var::InfiniteOpt.SemiInfiniteVariable,
    rvref::InfiniteOpt.GeneralVariableRef,
    index_type
    )::Nothing
    param_nums = var.parameter_nums
    ivref = var.infinite_variable_ref
    ivref_param_nums = InfiniteOpt._parameter_numbers(ivref)
    eval_supps = var.eval_supports
    # prepare for iterating over its supports
    supp_indices = support_index_iterator(trans_model, var.object_nums)
    vrefs = Vector{JuMP.VariableRef}(undef, length(supp_indices))
    labels = Vector{Set{DataType}}(undef, length(supp_indices))
    lookup_dict = Dict{Vector{Float64}, Int}()
    counter = 1
    # map a variable for each support
    for i in supp_indices
        raw_supp = index_to_support(trans_model, i)
        # ensure this support is valid with the reduced restriction
        if any(!isnan(raw_supp[ivref_param_nums[k]]) && raw_supp[ivref_param_nums[k]] != v 
               for (k, v) in eval_supps)
            continue
        end
        # map to the current transcription variable
        supp = raw_supp[param_nums]
        ivref_supp = [haskey(eval_supps, j) ? eval_supps[j] : raw_supp[k] 
                      for (j, k) in enumerate(ivref_param_nums)]
        @inbounds vrefs[counter] = lookup_by_support(trans_model, ivref, ivref_supp)
        lookup_dict[supp] = counter
        @inbounds labels[counter] = index_to_labels(trans_model, i)
        counter += 1
    end
    # truncate vrefs if any supports were skipped because of dependent parameter supps
    deleteat!(vrefs, counter:length(vrefs))
    deleteat!(labels, counter:length(vrefs))
    # save the transcription information
    data = transcription_data(trans_model)
    data.infvar_lookup[rvref] = lookup_dict
    data.infvar_mappings[rvref] = vrefs
    data.infvar_support_labels[rvref] = labels
    return
end

# Empty mapping dispatch for infinite parameter functions
function _set_semi_infinite_variable_mapping(
    trans_model::JuMP.Model,
    var::InfiniteOpt.SemiInfiniteVariable,
    rvref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.ParameterFunctionIndex}
    )::Nothing
    return
end

"""
    transcribe_semi_infinite_variables!(trans_model::JuMP.Model,
                                  inf_model::InfiniteOpt.InfiniteModel)::Nothing

Map each `SemiInfiniteVariable` in `inf_model` to transcription variables stored in
`trans_model`. The variable mappings are also stored in
`TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that [`transcribe_infinite_variables!`](@ref)
must be called first. Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.infvar_support_labels` is also populated.
"""
function transcribe_semi_infinite_variables!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.SemiInfiniteVariable)
        # get the basic variable information
        var = object.variable
        rvref = InfiniteOpt._make_variable_ref(inf_model, idx)
        # setup the mappings
        ivref = InfiniteOpt.infinite_variable_ref(rvref)
        _set_semi_infinite_variable_mapping(trans_model, var, rvref, 
                                            InfiniteOpt._index_type(ivref))
    end
    return
end

# Override the info of the jump variable with the point variable's if it is different
function _update_point_info(
    gvref::InfiniteOpt.GeneralVariableRef,
    vref::JuMP.VariableRef
    )::Nothing
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
    transcribe_point_variables!(trans_model::JuMP.Model,
                               inf_model::InfiniteOpt.InfiniteModel)::Nothing

Map each `PointVariable` in `inf_model` to a transcription variable stored in
`trans_model`. The variable mapping is also stored in
`TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref). Note that
[`transcribe_infinite_variables!`](@ref) must be called first and that the
info constraints associated with the transcription variable will be updated in
accordance with the point variable.
"""
function transcribe_point_variables!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.PointVariable)
        # get the basic variable information
        var = object.variable
        ivref = var.infinite_variable_ref
        supp = var.parameter_values
        # find the corresponding variable record the mapping
        vref = lookup_by_support(trans_model, ivref, supp)
        pvref = InfiniteOpt._make_variable_ref(inf_model, idx)
        transcription_data(trans_model).finvar_mappings[pvref] = vref
        # update the info constraints as needed
        _update_point_info(pvref, vref)
    end
    return
end

################################################################################
#                       TRANSCRIPTION EXPRESSION METHODS
################################################################################
"""
    transcription_expression(trans_model::JuMP.Model, expr, support::Vector{Float64})

Given the `expr` from an `InfiniteModel`, form its transcripted version in
accordance with the variable mappings available in `trans_model` defined at
`support`. This should only be used once all variables and measures have been
transcribed (e.g., via [`transcribe_finite_variables!`](@ref)).
"""
function transcription_expression(trans_model::JuMP.Model, expr, support)
    error("Unsupported expression type `$(typeof(expr))` for automated " *
          "transcription.")
end

# GeneralVariableRef
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    support::Vector{Float64}
    )
    return transcription_expression(trans_model, vref,
                                    InfiniteOpt._index_type(vref), support)
end

# Infinite variables, infinite parameter functions, and measures
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    ) where {V <: Union{InfVarIndex, InfiniteOpt.ParameterFunctionIndex, InfiniteOpt.MeasureIndex}}
    param_nums = InfiniteOpt._parameter_numbers(vref)
    return lookup_by_support(trans_model, vref, index_type, support[param_nums])
end

# Semi-Infinite variables
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.SemiInfiniteVariableIndex},
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
        return lookup_by_support(trans_model, vref, index_type, support[param_nums])
    end
end

# Point variables and finite variables
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    ) where {V <: FinVarIndex}
    return lookup_by_support(trans_model, vref, index_type, support)
end

# Infinite parameters
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    ) where {V <: InfiniteOpt.InfiniteParameterIndex}
    param_num = InfiniteOpt._parameter_number(vref)
    return support[param_num]
end

# Finite parameters
function transcription_expression(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.FiniteParameterIndex},
    support::Vector{Float64}
    )
    return InfiniteOpt.parameter_value(vref)
end

# AffExpr and QuadExpr
function transcription_expression(
    trans_model::JuMP.Model,
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},
    support::Vector{Float64}
    )
    return InfiniteOpt.map_expression(
        v -> transcription_expression(trans_model, v, support), 
        expr)
end

# NLPExpr
function transcription_expression(
    trans_model::JuMP.Model,
    nlp::InfiniteOpt.NLPExpr,
    support::Vector{Float64}
    )
    ast = InfiniteOpt.map_nlp_to_ast(
        v -> transcription_expression(trans_model, v, support), 
        nlp)
    return JuMP.add_NL_expression(trans_model, ast)
end

# Real Number 
function transcription_expression(
    trans_model::JuMP.Model,
    num::Real,
    support::Vector{Float64}
    )
    return zero(JuMP.AffExpr) + num
end

################################################################################
#                         MEASURE TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_measures!(trans_model::JuMP.Model,
                         inf_model::InfiniteOpt.InfiniteModel)::Nothing

For each `Measure` in `inf_model` expand it via `InfiniteOpt.expand_measure` or
`analytic_expansion` as appropriate and transcribe the expanded expression via
[`transcription_expression`](@ref). Then store the measure to transcripted
expression mappings in `TranscriptionData.measure_mappings` and
`TranscriptionData.measure_lookup` to enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`. 
Note that `TranscriptionData.measure_support_labels` is also populated.
"""
function transcribe_measures!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.Measure)
        # get the basic information
        meas = object.measure
        # expand the measure
        if meas.constant_func
            new_expr = InfiniteOpt.analytic_expansion(meas.func, meas.data, trans_model)
        else
            new_expr = InfiniteOpt.expand_measure(meas.func, meas.data, trans_model)
        end
        # prepare to transcribe over the supports
        supp_indices = support_index_iterator(trans_model, meas.object_nums)
        exprs = Vector{Any}(undef, length(supp_indices))
        labels = Vector{Set{DataType}}(undef, length(supp_indices))
        lookup_dict = Dict{Vector{Float64}, Int}()
        # map a variable for each support
        for (counter, i) in enumerate(supp_indices)
            raw_supp = index_to_support(trans_model, i)
            @inbounds exprs[counter] = transcription_expression(trans_model,
                                                                new_expr, raw_supp)
            supp = raw_supp[meas.parameter_nums]
            lookup_dict[supp] = counter
            @inbounds labels[counter] = index_to_labels(trans_model, i)
        end
        # save the transcription information
        mref = InfiniteOpt._make_variable_ref(inf_model, idx)
        data = transcription_data(trans_model)
        data.measure_lookup[mref] = lookup_dict
        data.measure_mappings[mref] = exprs
        data.measure_support_labels[mref] = labels
    end
    return
end

################################################################################
#                        OBJECTIVE TRANSCRIPTION METHODS
################################################################################
## Dispatch functions for setting the objective
# Normal Expr 
function _set_objective(trans_model, sense, expr)
    return JuMP.set_objective(trans_model, sense, expr)
end

# NonlinearExpression 
function _set_objective(trans_model, sense, expr::JuMP.NonlinearExpression)
    return JuMP.set_NL_objective(trans_model, sense, expr)
end

"""
    transcribe_objective!(trans_model::JuMP.Model,
                          inf_model::InfiniteOpt.InfiniteModel)::Nothing

Form the transcripted version of the objective stored in `inf_model` and add it
to `trans_model`. Note that all the variables and measures in `inf_model` must
by transcripted first (e.g., via [`transcribe_infinite_variables!`](@ref)).
"""
function transcribe_objective!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    expr = JuMP.objective_function(inf_model)
    sense = JuMP.objective_sense(inf_model)
    trans_expr = transcription_expression(trans_model, expr, Float64[])
    _set_objective(trans_model, sense, trans_expr)
    return
end

################################################################################
#                        CONSTRAINT TRANSCRIPTION METHODS
################################################################################
## Given a variable and its set from an info constraint, get the transcribed version
# Lower bound constraint
function _get_info_constr_from_var(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.GreaterThan,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.has_lower_bound(trans_vref) ? JuMP.LowerBoundRef(trans_vref) : nothing
end

# Upper bound constraint
function _get_info_constr_from_var(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.LessThan,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.has_upper_bound(trans_vref) ? JuMP.UpperBoundRef(trans_vref) : nothing
end

# Fix constraint
function _get_info_constr_from_var(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.EqualTo,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_fixed(trans_vref) ? JuMP.FixRef(trans_vref) : nothing
end

# Binary constraint
function _get_info_constr_from_var(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.ZeroOne,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_binary(trans_vref) ? JuMP.BinaryRef(trans_vref) : nothing
end

# Integer constraint
function _get_info_constr_from_var(
    trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.Integer,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_integer(trans_vref) ? JuMP.IntegerRef(trans_vref) : nothing
end

# Determine if a given raw support satisfies constraint domain restrictions
function _support_in_restrictions(
    support::Vector{Float64},
    indices::Vector{Int},
    domains::Vector{InfiniteOpt.IntervalDomain}
    )::Bool
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
    trans_model::JuMP.Model, 
    constr::JuMP.ScalarConstraint, 
    func::JuMP.AbstractJuMPScalar, 
    set::MOI.AbstractScalarSet, 
    raw_supp::Vector{Float64}, 
    name::String
    )
    new_func = transcription_expression(trans_model, func, raw_supp)
    trans_constr = JuMP.build_constraint(error, new_func, set)
    return JuMP.add_constraint(trans_model, trans_constr, name)
end

# MOI.LessThan expr 
function _make_constr_ast(ref, set::MOI.LessThan)
    return :($ref <= $(set.upper))
end

# MOI.LessGreat expr 
function _make_constr_ast(ref, set::MOI.GreaterThan)
    return :($ref >= $(set.lower))
end

# MOI.EqualTo expr 
function _make_constr_ast(ref, set::MOI.EqualTo)
    return :($ref == $(set.value))
end

# MOI.Interval expr 
function _make_constr_ast(ref, set::MOI.Interval)
    return :($(set.lower) <= $ref <= $(set.upper))
end

# MOI.Set fallback 
function _make_constr_ast(ref, set)
    error("TranscriptionOpt does not support constraint sets of type " * 
          "`$(typeof(set))` for general nonlinear constraints because this " *
          "is not yet supported by JuMP.")
end

# JuMP.ScalarConstraint with NLPExpr
function _process_constraint(
    trans_model::JuMP.Model, 
    constr::JuMP.ScalarConstraint, 
    func::InfiniteOpt.NLPExpr, 
    set::MOI.AbstractScalarSet, 
    raw_supp::Vector{Float64}, 
    name::String
    )
    nlp_ref = transcription_expression(trans_model, func, raw_supp)
    return JuMP.add_NL_constraint(trans_model, _make_constr_ast(nlp_ref, set))
end

# JuMP.VectorConstraint
function _process_constraint(
    trans_model::JuMP.Model, 
    constr::JuMP.VectorConstraint, 
    func::Vector{<:JuMP.AbstractJuMPScalar}, 
    set::MOI.AbstractVectorSet, 
    raw_supp::Vector{Float64}, 
    name::String
    )
    new_func = map(f -> transcription_expression(trans_model, f, raw_supp), func)
    if any(f -> f isa JuMP.NonlinearExpression, new_func)
        error("TranscriptionOpt does not support vector constraints of general " * 
              "nonlinear expressions because this is not yet supported by JuMP.")
    end
    shape = JuMP.shape(constr)
    shaped_func = JuMP.reshape_vector(new_func, shape)
    shaped_set = JuMP.reshape_set(set, shape)
    trans_constr = JuMP.build_constraint(error, shaped_func, shaped_set)
    return JuMP.add_constraint(trans_model, trans_constr, name)
end

# Fallback 
function _process_constraint(
    trans_model::JuMP.Model, 
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
    transcribe_constraints!(trans_model::JuMP.Model,
                            inf_model::InfiniteOpt.InfiniteModel)::Nothing

For each constraint in `inf_model` form its transcripted version(s) and add them
to `trans_model`. The mappings are stored in `TranscriptionData.constr_mappings`
and the associated supports are stored in `TranscriptionData.constr_supports`
to enable [`transcription_constraint`](@ref) and `InfiniteOpt.constraint_supports`.
Note that variable info constraints are simply mapped to the existing info
constraints already generated along with the transcription variables. Note that
the variables and measures must all first be transcripted (e.g., via
[`transcribe_measures!`](@ref)). Note that 
`TranscriptionData.constr_support_labels` is also populated.
"""
function transcribe_constraints!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    param_supps = parameter_supports(trans_model)
    for (idx, object) in inf_model.constraints
        # get the basic information
        constr = object.constraint
        func = JuMP.jump_function(constr)
        set = JuMP.moi_set(constr)
        obj_nums = object.object_nums
        cref = InfiniteOpt._make_constraint_ref(inf_model, idx)
        # prepare the iteration helpers
        supp_indices = support_index_iterator(trans_model, obj_nums)
        crefs = Vector{JuMP.ConstraintRef}(undef, length(supp_indices))
        supps = Vector{Tuple}(undef, length(supp_indices))
        labels = Vector{Set{DataType}}(undef, length(supp_indices))
        counter = 1
        # iterate over the support indices for the info constraints
        if object.is_info_constraint
            for i in supp_indices
                raw_supp = index_to_support(trans_model, i)
                info_ref = _get_info_constr_from_var(trans_model, func, set,
                                                     raw_supp)
                # not all supports may be defined if overwritten by a point variable
                if info_ref !== nothing
                    @inbounds crefs[counter] = info_ref
                    @inbounds supps[counter] = Tuple(param_supps[j][i[j]]
                                                     for j in obj_nums)
                    @inbounds labels[counter] = index_to_labels(trans_model, i)
                    counter += 1
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
                raw_supp = index_to_support(trans_model, i)
                # ensure the support satisfies parameter bounds and then add it
                if _support_in_restrictions(raw_supp, restrict_indices, 
                                            restrict_domains)
                    new_name = isempty(name) ? "" : string(name, "(support: ", counter, ")")
                    new_cref = _process_constraint(trans_model, constr, func, 
                                                   set, raw_supp, new_name) 
                    @inbounds crefs[counter] = new_cref
                    @inbounds supps[counter] = Tuple(param_supps[j][i[j]]
                                                     for j in obj_nums)
                    @inbounds labels[counter] = index_to_labels(trans_model, i)
                    counter += 1
                end
            end
        end
        # truncate the arrays in case not all the supports satisfied the bounds
        deleteat!(crefs, counter:length(crefs))
        deleteat!(supps, counter:length(supps))
        deleteat!(labels, counter:length(supps))
        # add the constraint mappings to the trans model
        data = transcription_data(trans_model)
        data.constr_mappings[cref] = crefs
        data.constr_supports[cref] = supps
        data.constr_support_labels[cref] = labels
    end
    return
end

################################################################################
#                      DERIVATIVE CONSTRAINT TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_derivative_evaluations!(trans_model::JuMP.Model, 
                                       inf_model::InfiniteOpt.InfiniteModel)::Nothing

Generate the auxiliary derivative evaluation equations and transcribe them 
appropriately for all the derivatives in `inf_model`. These are in turn added to 
`trans_model`. Note that no mapping information is recorded since the InfiniteModel 
won't have any constraints that correspond to these equations. Also Note that
the variables and measures must all first be transcripted (e.g., via
[`transcribe_derivative_variables!`](@ref)).
"""
function transcribe_derivative_evaluations!(
    trans_model::JuMP.Model, 
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (idx, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.Derivative)
        # get the basic variable information
        dref = InfiniteOpt._make_variable_ref(inf_model, idx)
        pref = dispatch_variable_ref(object.variable.parameter_ref)
        method = InfiniteOpt.derivative_method(pref)
        # generate the transcription constraints as needed
        if !InfiniteOpt.has_derivative_constraints(dref)
            # generate the evaluation expressions
            exprs = InfiniteOpt.evaluate_derivative(dref, method, trans_model)
            # prepare the iteration helpers
            param_obj_num = InfiniteOpt._object_number(pref)
            obj_nums = filter(!isequal(param_obj_num), InfiniteOpt._object_numbers(dref))
            supp_indices = support_index_iterator(trans_model, obj_nums)
            # transcribe the constraints 
            set = MOI.EqualTo(0.0)
            for i in supp_indices 
                raw_supp = index_to_support(trans_model, i)
                for expr in exprs
                    new_expr = transcription_expression(trans_model, expr, raw_supp)
                    trans_constr = JuMP.build_constraint(error, new_expr, set)
                    JuMP.add_constraint(trans_model, trans_constr) # TODO maybe add name?
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
    build_transcription_model!(trans_model::JuMP.Model,
                               inf_model::InfiniteOpt.InfiniteModel;
                               [check_support_dims::Bool = true])::Nothing

Given an empty `trans_model` build it using the information stored in `inf_model`.
This is intended for a `TranscriptionModel` that serves as a internal optimizer model
of `inf_model`. This detail is important to correctly enable internally generated
semi-infinite variables during the transcription process such that `inf_model` is not
modified. Note that this will add supports to `inf_model` via
[`InfiniteOpt.fill_in_supports!`](@ref) for infinite parameters that contain
no supports. Also a warning is thrown when the transcription model contains
more than 15,000 support points to alert users when they may naively have
a few independent supports whose product quickly yields a very large grid.
For example having 3 independent parameters with 100 supports each would result
in 1,000,000 supports if all three are together in at least 1 constraint. This 
behavior can be overcome using dependent parameters. The warning can be turned off 
via `check_support_dims = false`.
"""
function build_transcription_model!(
    trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel; check_support_dims::Bool = true
    )::Nothing
    # ensure there are supports to add and add them to the trans model
    InfiniteOpt.fill_in_supports!(inf_model, modify = false)
    set_parameter_supports(trans_model, inf_model)
    # check that there isn't a crazy amount of supports from taking the product
    supps = parameter_supports(trans_model)
    num_supps = isempty(supps) ? 0 : prod(length, supps)
    if check_support_dims && length(supps) > 1 && num_supps > 15000 # NOTE this is an arbitrary cutoff
        @warn("Due to necessarily considering the combinatorics of independent " *
              "parameter supports, the model will be transcripted over up to $(num_supps) " *
              "supports (if any constraint uses all the infinite parameters) "  *
              "and thus naive solution of the discretized problem may be slow. " * 
              "This warning can be turned off via `check_support_dims = false`.")
    end
    # define the variables
    transcribe_finite_variables!(trans_model, inf_model)
    transcribe_infinite_variables!(trans_model, inf_model)
    transcribe_derivative_variables!(trans_model, inf_model)
    transcribe_semi_infinite_variables!(trans_model, inf_model)
    transcribe_point_variables!(trans_model, inf_model)
    transcribe_measures!(trans_model, inf_model)
    # define the objective
    transcribe_objective!(trans_model, inf_model)
    # define the constraints
    transcribe_constraints!(trans_model, inf_model)
    # define the derivative evaluation constraints
    transcribe_derivative_evaluations!(trans_model, inf_model)
    return
end

# TODO make TranscriptionModel(::InfiniteModel) --> have to deal with accessing semi_infinite vars
