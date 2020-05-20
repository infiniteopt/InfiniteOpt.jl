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

# Build the parameter supports
"""
""" # TODO add docstring
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
# Make jump variables for hold variables
"""
""" # TODO add docstring
function transcribe_hold_variables(trans_model::JuMP.Model,
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

# Make jump variables for infinite variables
"""
""" # TODO add docstring
function transcribe_infinite_variables(trans_model::JuMP.Model,
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
    var::InfiniteOpt.ReducedVariable
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

# Map reduced variables to the correct infinite variable transcription variables
"""
""" # TODO add docstring
function transcribe_reduced_variables(trans_model::JuMP.Model,
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

# Override the info of the jump variable with the point variable's if any is provided
function _update_point_info(pvref::InfiniteOpt.PointVariableRef,
    vref::JuMP.VariableRef
    )::Nothing
    if JuMP.has_lower_bound(pvref) && !JuMP.has_lower_bound(vref)
        if JuMP.is_fixed(vref)
            JuMP.unfix(vref)
        end
        JuMP.set_lower_bound(vref, JuMP.lower_bound(pvref))
    end
    if JuMP.has_upper_bound(pvref) && !JuMP.has_upper_bound(vref)
        if JuMP.is_fixed(vref)
            JuMP.unfix(vref)
        end
        JuMP.set_upper_bound(vref, JuMP.upper_bound(pvref))
    end
    if JuMP.is_fixed(pvref) && !JuMP.is_fixed(vref)
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

# Map point variables to the correct infinite variable transcription variable
"""
""" # TODO add docstring
function transcribe_point_variables(trans_model::JuMP.Model,
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
""" # TODO add docstring
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

""" # TODO add doc string
function transcribe_measures(trans_model::JuMP.Model,
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
""" # TODO add doc string
function transcribe_objective(trans_model::JuMP.Model,
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
    )::JuMP.ConstraintRef
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.LowerBoundRef(trans_vref)
end

# Upper bound constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.LessThan,
    support::Vector{Float64}
    )::JuMP.ConstraintRef
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.UpperBoundRef(trans_vref)
end

# Fix constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.EqualTo,
    support::Vector{Float64}
    )::JuMP.ConstraintRef
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.FixRef(trans_vref)
end

# Binary constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.ZeroOne,
    support::Vector{Float64}
    )::JuMP.ConstraintRef
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.BinaryRef(trans_vref)
end

# Integer constraint
function _get_info_constr_from_var(trans_model::JuMP.Model,
    vref::InfiniteOpt.GeneralVariableRef,
    set::MOI.Integer,
    support::Vector{Float64}
    )::JuMP.ConstraintRef
    trans_vref = transcription_expression(trans_model, vref, support)
    return JuMP.IntegerRef(trans_vref)
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

""" # TODO add doc string
function transcribe_constraints(trans_model::JuMP.Model,
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
                crefs[counter] = _get_info_constr_from_var(trans_model, func,
                                                           set, raw_supp)
                supps[counter] = Tuple(param_supps[j][i[j]] for j in obj_nums)
                counter += 1
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
                    crefs[counter] = JuMP.add_constraint(trans_model, trans_constr,
                                                         new_name)
                    supps[counter] = Tuple(param_supps[j][i[j]] for j in obj_nums)
                    counter += 1
                end
            end
            # truncate the crefs array in case not all the supports satisfied the bounds
            deleteat!(crefs, counter:length(crefs))
            deleteat!(supps, counter:length(supps))
        end
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
""" # TODO add doc string
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
    transcribe_hold_variables(trans_model, inf_model)
    transcribe_infinite_variables(trans_model, inf_model)
    transcribe_reduced_variables(trans_model, inf_model)
    transcribe_point_variables(trans_model, inf_model)
    transcribe_measures(trans_model, inf_model)
    # define the objective
    transcribe_objective(trans_model, inf_model)
    # define the constraints
    transcribe_constraints(trans_model, inf_model)
    return
end

# TODO make TranscriptionModel(::InfiniteModel) --> have to deal with accessing reduced vars
