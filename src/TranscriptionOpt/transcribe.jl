################################################################################
#                           SUPPORT ITERATION METHODS
################################################################################
## Form a placeholder parameter reference given the object index
# IndependentParameterIndex
function _temp_parameter_ref(model::InfiniteOpt.InfiniteModel,
    index::InfiniteOpt.IndependentParameterIndex
    )::InfiniteOpt.IndependentParameterRef
    return InfiniteOpt.IndependentParameterRef(model, index)
end

# DependentParametersIndex
function _temp_parameter_ref(model::InfiniteOpt.InfiniteModel,
    index::InfiniteOpt.DependentParametersIndex
    )::InfiniteOpt.DependentParameterRef
    idx = InfiniteOpt.DependentParameterIndex(index, 1)
    return InfiniteOpt.DependentParameterRef(model, idx)
end

# Return the collected supports of an infinite parameter
function _collected_supports(model::InfiniteOpt.InfiniteModel,
    index::InfiniteOpt.InfiniteParameterIndex
    )::Vector
    pref = _temp_parameter_ref(model, index)
    supps = InfiniteOpt._parameter_supports(pref)
    supp_list = collect(keys(supps))
    # append a placeholder NaN support at the end to be used for efficient combinatorics
    return push!(supp_list, map(i -> NaN, first(supp_list)))
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
[`support_index_iterator`](@ref).
"""
function set_parameter_supports(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    param_indices = InfiniteOpt._param_object_indices(inf_model)
    supps = Tuple(_collected_supports(inf_model, i) for i in param_indices)
    transcription_data(trans_model).supports = supps
    return
end

################################################################################
#                        VARIABLE INITIALIZATION METHODS
################################################################################
"""
    transcribe_hold_variables!(trans_model::JuMP.Model,
                               inf_model::InfiniteOpt.InfiniteModel)::Nothing

Create a transcription variable (i.e., a JuMP variable) for each `HoldVariable`
stored in `inf_model` and add it to `trans_model`. The variable mapping is
also stored in `TranscriptionData.finvar_mappings` which enables
[`transcription_variable`](@ref) and [`lookup_by_support`](@ref).
"""
function transcribe_hold_variables!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (index, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.HoldVariable)
        hvref = InfiniteOpt._make_variable_ref(inf_model, index)
        vref = JuMP.add_variable(trans_model,
                                 JuMP.ScalarVariable(object.variable.info),
                                 object.name)
        transcription_data(trans_model).finvar_mappings[hvref] = vref
    end
    return
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
"""
function transcribe_infinite_variables!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (index, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.InfiniteVariable)
        # get the basic variable information
        var = object.variable
        base_name = InfiniteOpt._root_name(object.name)
        param_nums = var.parameter_nums
        # prepare for iterating over its supports
        supp_indices = support_index_iterator(trans_model, var.object_nums)
        vrefs = Vector{JuMP.VariableRef}(undef, length(supp_indices))
        lookup_dict = Dict{Vector{Float64}, Int}()
        counter = 1
        # create a variable for each support
        for i in supp_indices
            raw_supp = index_to_support(trans_model, i)
            supp = raw_supp[param_nums]
            var_name = string(base_name, "(support: ", counter, ")")
            @inbounds vrefs[counter] = JuMP.add_variable(trans_model,
                                               JuMP.ScalarVariable(var.info),
                                               var_name)
            lookup_dict[supp] = counter
            counter += 1
        end
        # save the transcription information
        ivref = InfiniteOpt._make_variable_ref(inf_model, index)
        transcription_data(trans_model).infvar_lookup[ivref] = lookup_dict
        transcription_data(trans_model).infvar_mappings[ivref] = vrefs
    end
    return
end

# Setup the mapping for a given reduced variable
function _set_reduced_variable_mapping(trans_model::JuMP.Model,
    var::InfiniteOpt.ReducedVariable,
    rvref::InfiniteOpt.GeneralVariableRef
    )::Nothing
    param_nums = var.parameter_nums
    ivref = var.infinite_variable_ref
    ivref_param_nums = InfiniteOpt._parameter_numbers(ivref)
    # prepare for iterating over its supports
    supp_indices = support_index_iterator(trans_model, var.object_nums)
    vrefs = Vector{JuMP.VariableRef}(undef, length(supp_indices))
    lookup_dict = Dict{Vector{Float64}, Int}()
    counter = 1
    # map a variable for each support
    for i in supp_indices
        raw_supp = index_to_support(trans_model, i)
        supp = raw_supp[param_nums]
        ivref_supp = raw_supp[ivref_param_nums]
        @inbounds vrefs[counter] = lookup_by_support(trans_model, ivref, ivref_supp)
        lookup_dict[supp] = counter
        counter += 1
    end
    # save the transcription information
    transcription_data(trans_model).infvar_lookup[rvref] = lookup_dict
    transcription_data(trans_model).infvar_mappings[rvref] = vrefs
    return
end

"""
    transcribe_reduced_variables!(trans_model::JuMP.Model,
                                  inf_model::InfiniteOpt.InfiniteModel)::Nothing

Map each `ReducedVariable` in `inf_model` to transcription variables stored in
`trans_model`. The variable mappings are also stored in
`TranscriptionData.infvar_mappings` in accordance with
`TranscriptionData.infvar_lookup` which enable [`transcription_variable`](@ref)
and [`lookup_by_support`](@ref). Note that [`transcribe_infinite_variables!`](@ref)
must be called first. Note that the supports will not be generated
until `InfiniteOpt.variable_supports` is invoked via `InfiniteOpt.supports`.
"""
function transcribe_reduced_variables!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (index, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.ReducedVariable)
        # get the basic variable information
        var = object.variable
        rvref = InfiniteOpt._make_variable_ref(inf_model, index)
        # setup the mappings
        _set_reduced_variable_mapping(trans_model, var, rvref)
    end
    return
end

# Override the info of the jump variable with the point variable's if it is different
function _update_point_info(gvref::InfiniteOpt.GeneralVariableRef,
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
function transcribe_point_variables!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (index, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.PointVariable)
        # get the basic variable information
        var = object.variable
        ivref = var.infinite_variable_ref
        param_nums = InfiniteOpt._parameter_numbers(ivref)
        supp = var.parameter_values
        # find the corresponding variable record the mapping
        vref = lookup_by_support(trans_model, ivref, supp)
        pvref = InfiniteOpt._make_variable_ref(inf_model, index)
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
transcribed (e.g., via [`transcribe_hold_variables!`](@ref)).
"""
function transcription_expression(trans_model::JuMP.Model, expr, support)
    error("Unsupported expression type `$(typeof(expr))` for automated " *
          "transcription.")
end

# GeneralVariableRef
function transcription_expression(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    support::Vector{Float64}
    )
    return transcription_expression(trans_model, vref,
                                    InfiniteOpt._index_type(vref), support)
end

# Infinite variables, reduced variables, and measures
function transcription_expression(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    )::JuMP.AbstractJuMPScalar where {V <: Union{InfVarIndex, InfiniteOpt.MeasureIndex}}
    param_nums = InfiniteOpt._parameter_numbers(vref)
    return lookup_by_support(trans_model, vref, index_type, support[param_nums])
end

# Point variables and hold variables
function transcription_expression(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    )::JuMP.VariableRef where {V <: FinVarIndex}
    return lookup_by_support(trans_model, vref, index_type, support)
end

# Infinite parameters
function transcription_expression(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{V},
    support::Vector{Float64}
    )::Float64 where {V <: InfiniteOpt.InfiniteParameterIndex}
    param_num = InfiniteOpt._parameter_num(vref)
    return support[param_num]
end

# Finite parameters
function transcription_expression(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    index_type::Type{InfiniteOpt.FiniteParameterIndex},
    support::Vector{Float64}
    )::Float64
    return JuMP.value(vref)
end

# AffExpr
function transcription_expression(trans_model::JuMP.Model,
    aff::JuMP.GenericAffExpr,
    support::Vector{Float64}
    )
    return JuMP.@expression(trans_model,
                 sum(coef * transcription_expression(trans_model, vref, support)
                     for (vref, coef) in aff.terms) + aff.constant)
end

# QuadExpr
function transcription_expression(trans_model::JuMP.Model,
    quad::JuMP.GenericQuadExpr,
    support::Vector{Float64}
    )
    return JuMP.@expression(trans_model,
                 sum(coef * transcription_expression(trans_model, pair.a, support) *
                     transcription_expression(trans_model, pair.b, support)
                     for (pair, coef) in quad.terms) +
                 transcription_expression(trans_model, quad.aff, support))
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
"""
function transcribe_measures!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    for (index, object) in InfiniteOpt._data_dictionary(inf_model, InfiniteOpt.Measure)
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
        exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(supp_indices))
        lookup_dict = Dict{Vector{Float64}, Int}()
        counter = 1
        # map a variable for each support
        for i in supp_indices
            raw_supp = index_to_support(trans_model, i)
            @inbounds exprs[counter] = transcription_expression(trans_model,
                                                                new_expr, raw_supp)
            supp = raw_supp[meas.param_nums]
            lookup_dict[supp] = counter
            counter += 1
        end
        # save the transcription information
        mref = InfiniteOpt._make_variable_ref(inf_model, index)
        transcription_data(trans_model).measure_lookup[mref] = lookup_dict
        transcription_data(trans_model).measure_mappings[mref] = exprs
    end
    return
end

################################################################################
#                        OBJECTIVE TRANSCRIPTION METHODS
################################################################################
"""
    transcribe_objective!(trans_model::JuMP.Model,
                          inf_model::InfiniteOpt.InfiniteModel)::Nothing

Form the transcripted version of the objective stored in `inf_model` and add it
to `trans_model`. Note that all the variables and measures in `inf_model` must
by transcripted first (e.g., via [`transcribe_infinite_variables!`](@ref)).
"""
function transcribe_objective!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    expr = JuMP.objective_function(inf_model)
    sense = JuMP.objective_sense(inf_model)
    trans_expr = transcription_expression(trans_model, expr, Float64[])
    JuMP.set_objective(trans_model, sense, trans_expr)
    return
end

################################################################################
#                        CONSTRAINT TRANSCRIPTION METHODS
################################################################################
## Given a variable and its set from an info constraint, get the transcribed version
# Lower bound constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.GreaterThan,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.has_lower_bound(trans_vref) ? JuMP.LowerBoundRef(trans_vref) : nothing
end

# Upper bound constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.LessThan,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.has_upper_bound(trans_vref) ? JuMP.UpperBoundRef(trans_vref) : nothing
end

# Fix constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.EqualTo,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_fixed(trans_vref) ? JuMP.FixRef(trans_vref) : nothing
end

# Binary constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.ZeroOne,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_binary(trans_vref) ? JuMP.BinaryRef(trans_vref) : nothing
end

# Integer constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.Integer,
    support::Vector{Float64}
    )::Union{JuMP.ConstraintRef, Nothing}
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.is_integer(trans_vref) ? JuMP.IntegerRef(trans_vref) : nothing
end

## Extract the parameter bounds from a constr if it has them and collect them
## to return them in a convenient form for iteration
# ScalarConstraint
function _get_constr_bounds(constr::JuMP.ScalarConstraint
    )::Tuple{Vector{Int}, Vector{InfiniteOpt.IntervalSet}}
    return Int[], InfiniteOpt.IntervalSet[]
end

# BoundedScalarConstraint
function _get_constr_bounds(constr::InfiniteOpt.BoundedScalarConstraint
    )::Tuple{Vector{Int}, Vector{InfiniteOpt.IntervalSet}}
    bounds = InfiniteOpt.parameter_bounds(constr)
    prefs = collect(keys(bounds))
    sets = map(p -> bounds[p], prefs)
    indices = map(p -> InfiniteOpt._parameter_number(p), prefs)
    return indices, sets
end

## Determine if a given raw support satisfies constraint parameter bounds
function _support_in_bounds(support::Vector{Float64},
    indices::Vector{Int},
    sets::Vector{InfiniteOpt.IntervalSet}
    )::Bool
    for i in eachindex(indices)
        s = support[i]
        if !isnan(s) && (s < JuMP.lower_bound(sets[i]) || s > JuMP.upper_bound(sets[i]))
            return false
        end
    end
    return true
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
[`transcribe_measures!`](@ref)).
"""
function transcribe_constraints!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    param_supps = parameter_supports(trans_model)
    for (index, object) in inf_model.constraints
        # get the basic information
        constr = object.constraint
        func = JuMP.jump_function(constr)
        set = JuMP.moi_set(constr)
        obj_nums = constr.object_nums
        # prepare the iteration helpers
        supp_indices = support_index_iterator(trans_model, obj_nums)
        crefs = Vector{JuMP.ConstraintRef}(undef, length(supp_indices))
        supps = Vector{Tuple}(undef, length(supp_indices))
        counter = 1
        # iterate over the support indices for the info constraints
        if object.is_info_constraint
            for i in supp_indices
                raw_supp = index_to_support(trans_model, i)
                info_ref = _get_info_constr_from_var(trans_model, func, set,
                                                     raw_supp)
                # not all supports may be defined if overwritten by a point variable
                if !isnothing(info_ref)
                    @inbounds crefs[counter] = info_ref
                    @inbounds supps[counter] = Tuple(param_supps[j][i[j]]
                                                     for j in obj_nums)
                    counter += 1
                end
            end
        # iterate over the supports for regular constraints
        else
            # get basic setup information
            bound_info = _get_constr_bounds(constr)
            name = object.name
            for i in supp_indices
                raw_supp = index_to_support(trans_model, i)
                # ensure the support satisfies parameter bounds and then add it
                if _support_in_bounds(raw_supp, bound_info[1], bound_info[2])
                    new_func = transcription_expression(trans_model, func, raw_supp)
                    trans_constr = JuMP.build_constraint(error, new_func, set)
                    new_name = string(name, "(support: ", counter, ")")
                    @inbounds crefs[counter] = JuMP.add_constraint(trans_model,
                                                         trans_constr, new_name)
                    @inbounds supps[counter] = Tuple(param_supps[j][i[j]]
                                                     for j in obj_nums)
                    counter += 1
                end
            end
        end
        # truncate the arrays in case not all the supports satisfied the bounds
        deleteat!(crefs, counter:length(crefs))
        deleteat!(supps, counter:length(supps))
        # add the constraint mappings to the trans model
        cref = InfiniteOpt._make_constraint_ref(inf_model, index)
        transcription_data(trans_model).constr_mappings[cref] = crefs
        transcription_data(trans_model).constr_supports[cref] = supps
    end
    return
end

################################################################################
#                      INFINITEMODEL TRANSCRIPTION METHODS
################################################################################
"""
    build_transcription_model!(trans_model::JuMP.Model,
                               inf_model::InfiniteOpt.InfiniteModel)::Nothing

Given an empty `trans_model` build it using the information stored in `inf_model`.
This is intended for a `TranscriptionModel` that serves as a internal optimizer model
of `inf_model`. This detail is important to correctly enable internally generated
reduced variables during the transcription process such that `inf_model` is not
modified. Note that this will add supports to `inf_model` via
[`InfiniteOpt.fill_in_supports!`](@ref) for infinite parameters that contain
no supports. Also a warning is thrown when the transcription model contains
more than 15,000 support points to alert users when they may naively have
a few independent supports whose product wuickly yields a very large grid.
For example having 3 independent parameters with 100 supports each would result
in 1,000,000 supports. This behavior can be overcome using dependent parameters.
"""
function build_transcription_model!(trans_model::JuMP.Model,
    inf_model::InfiniteOpt.InfiniteModel
    )::Nothing
    # ensure there are supports to add and add them to the trans model
    InfiniteOpt.fill_in_supports!(inf_model, modify = false)
    set_parameter_supports(trans_model, inf_model)
    # check that there isn't a crazy amount of supports from taking the product
    supps = parameter_supports(trans_model)
    num_supps = prod(length, supps)
    if length(supps) > 1 && num_supps > 15000 # NOTE this is an arbitrary cutoff
        @warn("Due to necessarily considering the combinatorics of independent " *
              "parameter supports, the model will be transcripted over $(num_supps) " *
              "supports and naive solution of the discretized problem may be slow.")
    end
    # define the variables
    transcribe_hold_variables!(trans_model, inf_model)
    transcribe_infinite_variables!(trans_model, inf_model)
    transcribe_reduced_variables!(trans_model, inf_model)
    transcribe_point_variables!(trans_model, inf_model)
    transcribe_measures!(trans_model, inf_model)
    # define the objective
    transcribe_objective!(trans_model, inf_model)
    # define the constraints
    transcribe_constraints!(trans_model, inf_model)
    return
end

# TODO make TranscriptionModel(::InfiniteModel) --> have to deal with accessing reduced vars
