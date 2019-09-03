# Make jump variables and a dict mapping global vars to jump vars
function _initialize_global_variables(trans_model::JuMP.Model,
                                      inf_model::InfiniteOpt.InfiniteModel)
    # search inf_model for global vars and make a jump var for one that is used
    for index in sort(collect(keys(inf_model.vars)))
        var = inf_model.vars[index]
        if isa(var, InfiniteOpt.GlobalVariable)
            gvref = InfiniteOpt.GlobalVariableRef(inf_model, index)
            if is_used(gvref)
                vref = JuMP.add_variable(trans_model,
                                         JuMP.ScalarVariable(var.info),
                                         JuMP.name(gvref))
                transcription_data(trans_model).global_to_var[gvref] = vref
            end
        end
    end
    return
end

# Return a vector of arrays containing the supports
function _list_supports(prefs::Tuple)::Vector
    return [InfiniteOpt.supports(pref) for pref in prefs]
end

# Make an index mapping for parameter support combinations for an infinite variable
function _make_supports(prefs::Tuple)::Vector
    support_list = _list_supports(prefs)
    combos = Iterators.product(support_list...)
    return [combo for combo in Iterators.take(combos, length(combos))]
end

# Make jump variables and a dict mapping infinite/point vars to jump vars
function _initialize_infinite_variables(trans_model::JuMP.Model,
                                        inf_model::InfiniteOpt.InfiniteModel)
    # search inf_model for infinite vars and make a jump var for all of its supports
    for index in sort(collect(keys(inf_model.vars)))
        var = inf_model.vars[index]
        if isa(var, InfiniteOpt.InfiniteVariable)
            ivref = InfiniteOpt.InfiniteVariableRef(inf_model, index)
            if InfiniteOpt.is_used(ivref)
                prefs = InfiniteOpt.parameter_refs(ivref)
                supports = _make_supports(prefs)
                transcription_data(trans_model).infvar_to_supports[ivref] = supports
                vrefs = Vector{JuMP.VariableRef}(undef, length(supports))
                name = InfiniteOpt._root_name(ivref)
                for i in eachindex(vrefs)
                    # TODO Perhaps add different naming options...
                    vrefs[i] = JuMP.add_variable(trans_model,
                                                 JuMP.ScalarVariable(var.info),
                                                 string(name, "(support: ", i, ")"))
                end
                transcription_data(trans_model).infinite_to_vars[ivref] = vrefs
            end
        end
    end
    return
end

# Map the point variable reference to a transcribed variable based off of the support
function _update_point_mapping(trans_model::JuMP.Model,
                               pvref::InfiniteOpt.PointVariableRef,
                               ivref::InfiniteOpt.InfiniteVariableRef,
                               support::Tuple)
    supps = InfiniteOpt.supports(trans_model, ivref)
    # search for the support and update mapping
    for i in eachindex(supps)
       if all(isapprox.(support, supps[i]))
           vrefs = transcription_variable(trans_model, ivref)
           transcription_data(trans_model).point_to_var[pvref] = vrefs[i]
           return
       end
    end
    # fallback that should not be needed
    error("Couldn't find variable to map $pvref to.")
    return
end

# Override the info of the jump variable with the point variable's if any is provided
function _update_point_info(trans_model::JuMP.Model,
                            pvref::InfiniteOpt.PointVariableRef)
    vref = transcription_variable(trans_model, pvref)
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

# Map point variables to the correct transcribed infinite variable
function _map_point_variables(trans_model::JuMP.Model,
                              inf_model::InfiniteOpt.InfiniteModel)
    # search inf_model for point vars and map them to jump vars if they are used
    for (index, var) in inf_model.vars
        if isa(var, InfiniteOpt.PointVariable)
            pvref = InfiniteOpt.PointVariableRef(inf_model, index)
            if InfiniteOpt.is_used(pvref)
                ivref = InfiniteOpt.infinite_variable_ref(pvref)
                support = InfiniteOpt.parameter_values(pvref)
                _update_point_mapping(trans_model, pvref, ivref, support)
                _update_point_info(trans_model, pvref)
            end
        end
    end
    return
end

# Return the tuple index of a parameter based off of group_id
function _parameter_tuple_index(pref::InfiniteOpt.ParameterRef,
                                prefs::Tuple)::Int
    group = InfiniteOpt.group_id(pref)
    groups = InfiniteOpt._group.(prefs)
    if !(group in groups)
        return -1
    else
        return findfirst(isequal(group), groups)
    end
end

# Return the support value corresponding to a parameter reference
# return NaN is the parameter is not contained in prefs
function _parameter_value(pref::InfiniteOpt.ParameterRef, support::Tuple,
                          prefs::Tuple)::Float64
    pref_index = _parameter_tuple_index(pref, prefs)
    if pref_index == -1
        return NaN
    elseif isa(prefs[pref_index], InfiniteOpt.ParameterRef)
        return support[pref_index]
    else
        for (k, v) in prefs[pref_index].data
            if v == pref
                return support[pref_index].data[k]
            end
        end
    end
    return NaN
end

## Helper function for mapping InfiniteOpt variables to jump variables
# FiniteVariableRef
function _map_to_variable(fvref::InfiniteOpt.FiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)::JuMP.VariableRef
    return transcription_variable(trans_model, fvref)
end

# InfiniteVariableRef
function _map_to_variable(ivref::InfiniteOpt.InfiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)::JuMP.VariableRef
    # reduce support to only include the relavent parameter id groups
    ivref_groups = InfiniteOpt._group.(InfiniteOpt.parameter_refs(ivref))
    support_groups = InfiniteOpt._group.(prefs)
    reduced_support = Tuple(support[findfirst(isequal(group), support_groups)]
                            for group in ivref_groups)
    # find the jump variable associated with the support
    supps = InfiniteOpt.supports(trans_model, ivref)
    for i in eachindex(supps)
        if all(isapprox.(reduced_support, supps[i]))
            return transcription_variable(trans_model, ivref)[i]
        end
    end
    # this shouldn't be needed, but provided as a backup
    error("Couldn't find JuMP variable corresponding to $ivref.")
end

# ReducedInfiniteVariableRef
function _map_to_variable(rvref::InfiniteOpt.ReducedInfiniteVariableRef,
                          support::Tuple, prefs::Tuple,
                          trans_model::JuMP.Model)::JuMP.VariableRef
    # parse the reduced parameters and modify the support to include them
    ivref = InfiniteOpt.infinite_variable_ref(rvref)
    orig_groups = InfiniteOpt._group.(InfiniteOpt.parameter_refs(ivref))
    support_groups = InfiniteOpt._group.(prefs)
    support = Any[i for i in support]
    prefs = Any[pref for pref in prefs]
    for (index, value) in InfiniteOpt.eval_supports(rvref)
        if orig_groups[index] in support_groups
            support[findfirst(isequal(orig_groups[index]), support_groups)] = value
        else
            push!(support, value)
            push!(prefs, InfiniteOpt.parameter_refs(ivref)[index])
        end
    end
    prefs = Tuple(pref for pref in prefs)
    support = Tuple(i for i in support)
    # call the infinite variable function with the updated support and prefs
    return _map_to_variable(ivref, support, prefs, trans_model)
end

# ParameterRef
function _map_to_variable(pref::InfiniteOpt.ParameterRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)::Float64
    # find pref in prefs and return associated support value
    value = _parameter_value(pref, support, prefs)
    # this error is a backup that should not be needed
    value === NaN && error("Couldn't find support corresponding to $pref.")
    return value
end

# Return a Bool whether the support satisfies the bounds or not
function _support_in_bounds(support::Tuple, prefs::Tuple, bounds::Dict)::Bool
    for (pref, set) in bounds
        value = _parameter_value(pref, support, prefs)
        if value === NaN
            continue
        end
        if value < set.lower_bound || value > set.upper_bound
            return false
        end
    end
    return true
end

# Return a bool array indicating which supports are in bounds
function _supports_in_bounds(supports::Vector, prefs::Tuple, bounds::Dict)::Vector
    return [_support_in_bounds(supports[i], prefs, bounds)
            for i in eachindex(supports)]
end

## Convert jump scalar expressions with InfiniteOpt variables into transcribed relations
# PointVariableRef and GlobalVariableRef --> return scalar jump object
function _make_transcription_function(vref::InfiniteOpt.FiniteVariableRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple
    return transcription_variable(trans_model, vref), ()
end

# InfiniteVariableRef --> return tuple of expressions, prefs, and support mapping
function _make_transcription_function(vref::InfiniteOpt.InfiniteVariableRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple
    # Update the supports to include only those in the bounds and return vars
    if length(bounds) != 0
        supports = InfiniteOpt.supports(trans_model, vref)
        prefs = InfiniteOpt.parameter_refs(vref)
        old_support_indices = _supports_in_bounds(supports, prefs, bounds)
        new_supports = supports[old_support_indices]
        return transcription_variable(trans_model, vref)[old_support_indices],
               InfiniteOpt.parameter_refs(vref), new_supports
    # easy case just return the info
    else
        return transcription_variable(trans_model, vref),
               InfiniteOpt.parameter_refs(vref),
               InfiniteOpt.supports(trans_model, vref)
    end
end

# ParameterRef --> return tuple of numbers, pref, and support mapping
function _make_transcription_function(pref::InfiniteOpt.ParameterRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple
    supports = InfiniteOpt.supports(pref)
    # truncate the supports if out of bounds
    if pref in keys(bounds)
        old_support_indices = bounds[pref].lower_bound .<= supports .<= bounds[pref].upper_bound
        new_supports = supports[old_support_indices]
        return new_supports, (pref, ), new_supports
    # easy case just return the supports
    else
        return supports, (pref, ), supports
    end
end

# MeasureRef --> return depends if finite or infinite
function _make_transcription_function(mref::InfiniteOpt.MeasureRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple
    # expand the measure
    func = InfiniteOpt.measure_function(mref)
    data = InfiniteOpt.measure_data(mref)
    new_func = InfiniteOpt._possible_convert(InfiniteOpt.FiniteVariableRef,
                                         InfiniteOpt._expand_measure(func, data,
                                            trans_model, _update_point_mapping))
    # will either call finite variable function or general variable function
    return _make_transcription_function(new_func, trans_model, bounds)
end

# GenericAffExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericAffExpr{C,
                                               <:InfiniteOpt.FiniteVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    # replace finite vars with jump vars
    pairs = [transcription_variable(trans_model, var) => coef
             for (var, coef) in expr.terms]
    return JuMP.GenericAffExpr(expr.constant,
                         JuMP._new_ordered_dict(JuMP.VariableRef, C, pairs)), ()
end

# GenericQuadExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C,
                                               <:InfiniteOpt.FiniteVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    # replace finite vars with jump vars
    pairs = Vector{Pair{JuMP.UnorderedPair{JuMP.VariableRef}, C}}(undef,
                                                             length(expr.terms))
    counter = 1
    for k in keys(expr.terms)
        a = transcription_variable(trans_model, k.a)
        b = transcription_variable(trans_model, k.b)
        pairs[counter] = JuMP.UnorderedPair(a, b) => expr.terms[k]
        counter += 1
    end
    aff = _make_transcription_function(expr.aff, trans_model)[1]
    return JuMP.GenericQuadExpr(aff,
     JuMP._new_ordered_dict(JuMP.UnorderedPair{JuMP.VariableRef}, C, pairs)), ()
end

# GenericAffExpr of GeneralVariableRefs --> return tuple of numbers, pref,
# and support mapping
function _make_transcription_function(expr::JuMP.GenericAffExpr{C,
                                              <:InfiniteOpt.GeneralVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 1
        var = first(keys(expr.terms))
        results = _make_transcription_function(var, trans_model, bounds)
        # results[2] contains prefs if is infinite or it is empty otherwise
        if length(results[2]) != 0
            exprs = expr.terms[var] * results[1] .+ expr.constant
            return exprs, results[2], results[3]
        else
            return expr.terms[var] * results[1] + expr.constant, ()
        end
    end
    # check to see if there are measures and expand them
    # TODO rework paradigm so expressions don't need to be checked for measures
    has_measure = false
    for var in keys(expr.terms)
        if var isa InfiniteOpt.MeasureRef
            has_measure = true
            break
        end
    end
    if has_measure
        expr = InfiniteOpt._expand_measures(expr, trans_model,
                                            _update_point_mapping)
    end
    # dispatch to quadratic method if the measures contained quadratic terms
    if isa(expr, JuMP.GenericQuadExpr)
        return _make_transcription_function(expr, trans_model, bounds)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = InfiniteOpt._all_parameter_refs(expr)
    supports = _make_supports(prefs)
    if length(bounds) != 0
        old_support_indices = _supports_in_bounds(supports, prefs, bounds)
        supports = supports[old_support_indices]
    end
    exprs = [zero(JuMP.GenericAffExpr{C, JuMP.VariableRef})
             for i in eachindex(supports)]
    # make an expression for each support
    for i in eachindex(exprs)
        exprs[i].constant = expr.constant
        for (var, coef) in expr.terms
            # replace each variable with appropriate jump var
            new_var = _map_to_variable(var, supports[i], prefs,
                                       trans_model)
            if isa(new_var, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, new_var)
            else
                JuMP.add_to_expression!(exprs[i], coef * new_var)
            end
        end
    end
    return exprs, prefs, supports
end

# GenericQuadExpr of GeneralVariableRefs
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C,
                                              <:InfiniteOpt.GeneralVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 0 && length(expr.aff.terms) == 1
      var = first(keys(expr.aff.terms))
      results = _make_transcription_function(var, trans_model, bounds)
      # results[2] contains prefs if is infinite or it is empty otherwise
      if length(results[2]) != 0
          exprs = expr.aff.terms[var] * results[1] .+ expr.aff.constant
          return exprs, results[2], results[3]
      else
          return expr.aff.terms[var] * results[1] + expr.aff.constant, ()
      end
    end
    # check to see if there are measures and expand them
    # TODO rework paradigm so expressions don't need to be checked for measures
    all_variables = InfiniteOpt._all_function_variables(expr)
    has_measure = false
    for var in all_variables
        if var isa InfiniteOpt.MeasureRef
            has_measure = true
            break
        end
    end
    if has_measure
        # TODO expand measures in place by deleting entry and appending to the end
        # TODO thus the measure search can be done by expand measures
        expr = InfiniteOpt._expand_measures(expr, trans_model,
                                            _update_point_mapping)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = InfiniteOpt._all_parameter_refs(expr)
    supports = _make_supports(prefs)
    if length(bounds) != 0
        old_support_indices = _supports_in_bounds(supports, prefs, bounds)
        supports = supports[old_support_indices]
    end
    exprs = [zero(JuMP.GenericQuadExpr{C, JuMP.VariableRef})
             for i in eachindex(supports)]
    # make an expression for each support
    for i in eachindex(exprs)
        exprs[i].aff.constant = expr.aff.constant
        for (var, coef) in expr.aff.terms
            # replace each variable with appropriate jump var
            new_var = _map_to_variable(var, supports[i], prefs,
                                       trans_model)
            if isa(new_var, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, new_var)
            else
                JuMP.add_to_expression!(exprs[i], coef * new_var)
            end
        end
        for (pair, coef) in expr.terms
            # replace each variable with appropriate jump var
            var_a = _map_to_variable(pair.a, supports[i], prefs,
                                     trans_model)
            var_b = _map_to_variable(pair.b, supports[i], prefs,
                                     trans_model)
            if isa(var_a, JuMP.VariableRef) && isa(var_b, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, var_a, var_b)
            elseif isa(var_b, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef * var_a, var_b)
            elseif isa(var_a, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef * var_b, var_a)
            else
                JuMP.add_to_expression!(exprs[i], coef * var_a * var_b)
            end
        end
    end
    return exprs, prefs, supports
end

# GenericAffExpr and GenericQuadExpr of MeasureFiniteVariableRefs --> return
# depends on whether finite
function _make_transcription_function(expr::Union{JuMP.GenericAffExpr{C,
                                        <:InfiniteOpt.MeasureFiniteVariableRef},
                                        JuMP.GenericQuadExpr{C,
                                        <:InfiniteOpt.MeasureFiniteVariableRef}},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    expr = InfiniteOpt._possible_convert(FiniteVariableRef,
                                         InfiniteOpt._expand_measures(expr,
                                         trans_model, _update_point_mapping))
    return _make_transcription_function(expr, trans_model, bounds)
end

# Empty jump variable expr (for constraints of form number <= number)
function _make_transcription_function(expr::JuMP.GenericAffExpr{C,
                                                              JuMP.VariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple where {C}
    return expr, ()
end

# Fall back function for other jump objects
function _make_transcription_function(expr::JuMP.AbstractJuMPScalar,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())::Tuple
    type = typeof(expr)
    error("Unsupported transcription of expression of type $type.")
end

## Construct the objective and error is contains non finite variables
function _set_objective(trans_model::JuMP.Model,
                        inf_model::InfiniteOpt.InfiniteModel)
    trans_obj, = _make_transcription_function(JuMP.objective_function(inf_model),
                                              trans_model)
    isa(trans_obj, Vector) && error("Objective is not finite, ensure all " *
                                    "infinite variables/parameters in measures " *
                                    "are evaluated completely.")
    JuMP.set_objective(trans_model, JuMP.objective_sense(inf_model), trans_obj)
    return
end

## Define helper functions for setting constraint mappings
# InfiniteConstraintRef
function _set_mapping(icref::InfiniteOpt.InfiniteConstraintRef,
                      crefs::Vector{<:JuMP.ConstraintRef})
    transcription_data(JuMP.owner_model(crefs[1])).infinite_to_constrs[icref] = crefs
    return
end

# MeasureConstraintRef (infinite)
function _set_mapping(mcref::InfiniteOpt.MeasureConstraintRef,
                      crefs::Vector{<:JuMP.ConstraintRef})
    transcription_data(JuMP.owner_model(crefs[1])).measure_to_constrs[mcref] = crefs
    return
end

# MeasureConstraintRef (finite)
function _set_mapping(mcref::InfiniteOpt.MeasureConstraintRef,
                      cref::JuMP.ConstraintRef)
    transcription_data(JuMP.owner_model(cref)).measure_to_constrs[mcref] = [cref]
    return
end

# FiniteConstraintRef
function _set_mapping(fcref::InfiniteOpt.FiniteConstraintRef,
                      cref::JuMP.ConstraintRef)
    transcription_data(JuMP.owner_model(cref)).finite_to_constr[fcref] = cref
    return
end

## Define helper functions for setting constraint supports
# InfiniteConstraintRef
function _set_supports(trans_model::JuMP.Model,
                       icref::InfiniteOpt.InfiniteConstraintRef,
                       supports::Vector)
    transcription_data(trans_model).infconstr_to_supports[icref] = supports
    return
end

# MeasureConstraintRef
function _set_supports(trans_model::JuMP.Model,
                       mcref::InfiniteOpt.MeasureConstraintRef,
                       supports::Vector)
    transcription_data(trans_model).measconstr_to_supports[mcref] = supports
    return
end

## Define helper functions for setting constraint parameter reference tuples
# InfiniteConstraintRef
function _set_parameter_refs(trans_model::JuMP.Model,
                             icref::InfiniteOpt.InfiniteConstraintRef,
                             prefs::Tuple)
    transcription_data(trans_model).infconstr_to_params[icref] = prefs
    return
end

# MeasureConstraintRef
function _set_parameter_refs(trans_model::JuMP.Model,
                             mcref::InfiniteOpt.MeasureConstraintRef,
                             prefs::Tuple)
    transcription_data(trans_model).measconstr_to_params[mcref] = prefs
    return
end

# Extract the root name of a constraint
function InfiniteOpt._root_name(cref::InfiniteOpt.GeneralConstraintRef)
    name = JuMP.name(cref)
    if length(name) == 0
        return "noname"
    else
        first_bracket = findfirst(isequal('['), name)
        if first_bracket == nothing
            return name
        else
            return name[1:first_bracket-1]
        end
    end
end

## leverage _make_transcription_function to transcribe the constraints
function _set_constraints(trans_model::JuMP.Model, inf_model::InfiniteOpt.InfiniteModel)
    # transform and add constraints that haven't already been added through add_variable
    ordered_indices = sort(collect(keys(inf_model.constrs)))
    for index in ordered_indices
        constr = inf_model.constrs[index]
        if !inf_model.constr_in_var_info[index]
            # extract the reference and transcribe the jump object function
            icref = InfiniteOpt._make_constraint_ref(inf_model, index)
            if isa(constr, BoundedScalarConstraint)
                results = _make_transcription_function(constr.func, trans_model,
                                                       constr.bounds)
            else
                results = _make_transcription_function(constr.func, trans_model)
            end
            # extract the name and create a transcribed constraint
            root_name = InfiniteOpt._root_name(icref)
            # results[2] contains prefs if is infinite or it is empty otherwise
            if length(results[2]) != 0
                crefs = Vector{JuMP.ConstraintRef}(undef, length(results[1]))
                for i in eachindex(crefs)
                    con = JuMP.build_constraint(error, results[1][i], constr.set)
                    # TODO Perhaps improve naming
                    name = string(root_name, "(Support: ", i, ")")
                    crefs[i] = JuMP.add_constraint(trans_model, con, name)
                end
                # update the mappings
                _set_mapping(icref, crefs)
                _set_parameter_refs(trans_model, icref, results[2])
                _set_supports(trans_model, icref, results[3])
            else
                # We have a finite constraint, just add normally and update mapping.
                con = JuMP.build_constraint(error, results[1], constr.set)
                cref = JuMP.add_constraint(trans_model, con, root_name)
                _set_mapping(icref, cref)
            end
        end
    end
    return
end

## Helper functions for mapping InfiniteOpt var info constrs to the transcribed ones
# FiniteVariableRef
function _map_info_constraints(fvref::InfiniteOpt.FiniteVariableRef,
                               trans_model::JuMP.Model)
    vref = transcription_variable(trans_model, fvref)
    # Check if variables have info constraints and map if they do
    if JuMP.has_lower_bound(fvref)
        _set_mapping(JuMP.LowerBoundRef(fvref), JuMP.LowerBoundRef(vref))
    end
    if JuMP.has_upper_bound(fvref)
        _set_mapping(JuMP.UpperBoundRef(fvref), JuMP.UpperBoundRef(vref))
    end
    if JuMP.is_fixed(fvref)
        _set_mapping(JuMP.FixRef(fvref), JuMP.FixRef(vref))
    end
    if JuMP.is_integer(fvref)
        _set_mapping(JuMP.IntegerRef(fvref), JuMP.IntegerRef(vref))
    end
    if JuMP.is_binary(fvref)
        _set_mapping(JuMP.BinaryRef(fvref), JuMP.BinaryRef(vref))
    end
    return
end

# InfiniteVariableRef
function _map_info_constraints(ivref::InfiniteOpt.InfiniteVariableRef,
                               trans_model::JuMP.Model)
    vrefs = transcription_variable(trans_model, ivref)
    # TODO Prealocate the size
    # Check if both variables have a constraint and map if they do
    if JuMP.has_lower_bound(ivref)
        crefs = JuMP.ConstraintRef[]
        supports = Tuple[]
        for i in eachindex(vrefs)
            if JuMP.has_lower_bound(vrefs[i])
                push!(crefs, JuMP.LowerBoundRef(vrefs[i]))
                push!(supports, InfiniteOpt.supports(trans_model, ivref)[i])
            end
        end
        if length(crefs) != 0
            _set_mapping(JuMP.LowerBoundRef(ivref), crefs)
            _set_parameter_refs(trans_model, JuMP.LowerBoundRef(ivref),
                                InfiniteOpt.parameter_refs(ivref))
            _set_supports(trans_model, JuMP.LowerBoundRef(ivref), supports)
        end
    end
    if JuMP.has_upper_bound(ivref)
        crefs = JuMP.ConstraintRef[]
        supports = Tuple[]
        for i in eachindex(vrefs)
            if JuMP.has_upper_bound(vrefs[i])
                push!(crefs, JuMP.UpperBoundRef(vrefs[i]))
                push!(supports, InfiniteOpt.supports(trans_model, ivref)[i])
            end
        end
        if length(crefs) != 0
            _set_mapping(JuMP.UpperBoundRef(ivref), crefs)
            _set_parameter_refs(trans_model, JuMP.UpperBoundRef(ivref),
                                InfiniteOpt.parameter_refs(ivref))
            _set_supports(trans_model, JuMP.UpperBoundRef(ivref), supports)
        end
    end
    if JuMP.is_fixed(ivref)
        crefs = JuMP.ConstraintRef[]
        supports = Tuple[]
        for i in eachindex(vrefs)
            if JuMP.is_fixed(vrefs[i])
                push!(crefs, JuMP.FixRef(vrefs[i]))
                push!(supports, InfiniteOpt.supports(trans_model, ivref)[i])
            end
        end
        if length(crefs) != 0
            _set_mapping(JuMP.FixRef(ivref), crefs)
            _set_parameter_refs(trans_model, JuMP.FixRef(ivref),
                                InfiniteOpt.parameter_refs(ivref))
            _set_supports(trans_model, JuMP.FixRef(ivref), supports)
        end
    end
    if JuMP.is_integer(ivref)
        crefs = JuMP.ConstraintRef[]
        supports = Tuple[]
        for i in eachindex(vrefs)
            if JuMP.is_integer(vrefs[i])
                push!(crefs, JuMP.IntegerRef(vrefs[i]))
                push!(supports, InfiniteOpt.supports(trans_model, ivref)[i])
            end
        end
        if length(crefs) != 0
            _set_mapping(JuMP.IntegerRef(ivref), crefs)
            _set_parameter_refs(trans_model, JuMP.IntegerRef(ivref),
                                InfiniteOpt.parameter_refs(ivref))
            _set_supports(trans_model, JuMP.IntegerRef(ivref), supports)
        end
    end
    if JuMP.is_binary(ivref)
        crefs = JuMP.ConstraintRef[]
        supports = Tuple[]
        for i in eachindex(vrefs)
            if JuMP.is_binary(vrefs[i])
                push!(crefs, JuMP.BinaryRef(vrefs[i]))
                push!(supports, InfiniteOpt.supports(trans_model, ivref)[i])
            end
        end
        if length(crefs) != 0
            _set_mapping(JuMP.BinaryRef(ivref), crefs)
            _set_parameter_refs(trans_model, JuMP.BinaryRef(ivref),
                                InfiniteOpt.parameter_refs(ivref))
            _set_supports(trans_model, JuMP.BinaryRef(ivref), supports)
        end
    end
    return
end

## Map the variable info constraints between the two models
function _map_variable_info_constraints(trans_model::JuMP.Model,
                                        inf_model::InfiniteOpt.InfiniteModel)
    for index in keys(inf_model.vars)
        ivref = InfiniteOpt._make_variable_ref(inf_model, index)
        if InfiniteOpt.is_used(ivref)
            _map_info_constraints(ivref, trans_model)
        end
    end
    return
end

"""
    TranscriptionModel(model::InfiniteModel, args...)

Return a `TranscriptionModel` of `model`. This transcribes all of the variables,
constraints, and objective.

**Example**
```julia
julia> TranscriptionModel(model)
A JuMP Model
Feasibility problem with:
Variables: 130
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 25 constraints
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.GreaterThan{Float64}`: 100 constraint
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.LessThan{Float64}`: 84 constraints
`VariableRef`-in-`MathOptInterface.EqualTo{Float64}`: 40 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 40 constraints
`VariableRef`-in-`MathOptInterface.LessThan{Float64}`: 25 constraints
`VariableRef`-in-`MathOptInterface.Integer`: 40 constraints
`VariableRef`-in-`MathOptInterface.ZeroOne`: 40 constraints
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
function TranscriptionModel(inf_model::InfiniteOpt.InfiniteModel,
                            optimizer_factory::Union{JuMP.OptimizerFactory,
                                                     Nothing} = nothing;
                            kwargs...)::JuMP.Model
    if isa(optimizer_factory, Nothing)
        trans_model = TranscriptionModel(; kwargs...)
    else
        trans_model = TranscriptionModel(optimizer_factory; kwargs...)
    end
    _initialize_global_variables(trans_model, inf_model)
    _initialize_infinite_variables(trans_model, inf_model)
    _map_point_variables(trans_model, inf_model)
    if JuMP.objective_sense(inf_model) != MOI.FEASIBILITY_SENSE
        _set_objective(trans_model, inf_model)
    end
    _map_variable_info_constraints(trans_model, inf_model)
    # TODO optimize performance --> bottlenecked with _make_transcription_function
    _set_constraints(trans_model, inf_model)
    return trans_model
end
