# Make jump variables and a dict mapping global vars to jump vars
function _initialize_global_variables(trans_model::JuMP.Model,
                                      inf_model::InfiniteOpt.InfiniteModel)
    # search inf_model for global vars and make a jump var for one that is used
    ordered_indices = sort(collect(keys(inf_model.vars)))
    for index in ordered_indices
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
function _list_supports(prefs::Tuple)
    support_list = Vector{Array}(undef, length(prefs))
    for i = 1:length(prefs)
        support_list[i] = InfiniteOpt.supports(prefs[i])
    end
    return support_list
end

# Make an index mapping for parameter support combinations for an infinite variable
function _make_support_indices(prefs::Tuple)
    support_indices = Dict{Int, Tuple}()
    support_list = _list_supports(prefs)
    index = 1
    for combo in Iterators.product(support_list...)
        support_indices[index] = combo
        index += 1
    end
    return support_indices
end

# Make jump variables and a dict mapping infinite/point vars to jump vars
function _initialize_infinite_variables(trans_model::JuMP.Model,
                                        inf_model::InfiniteOpt.InfiniteModel)
    # search inf_model for infinite vars and make a jump var for all of its supports
    ordered_indices = sort(collect(keys(inf_model.vars)))
    for index in ordered_indices
        var = inf_model.vars[index]
        if isa(var, InfiniteOpt.InfiniteVariable)
            ivref = InfiniteOpt.InfiniteVariableRef(inf_model, index)
            if InfiniteOpt.is_used(ivref)
                prefs = InfiniteOpt.parameter_refs(ivref)
                transcription_data(trans_model).infvar_to_supports[ivref] = _make_support_indices(prefs)
                vrefs = Vector{JuMP.VariableRef}(undef, length(keys(supports(trans_model, ivref))))
                for i = 1:length(vrefs)
                    # TODO Perhaps add different naming options...
                    name = InfiniteOpt._root_name(ivref)
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
function _update_point_mapping(trans_model::JuMP.Model, pvref::InfiniteOpt.PointVariableRef,
                               ivref::InfiniteOpt.InfiniteVariableRef, support::Tuple)
    for (k, v) in InfiniteOpt.supports(trans_model, ivref)
       # if isequal(v, support)
       if v == support
           vrefs = transcription_variable(trans_model, ivref)
           transcription_data(trans_model).point_to_var[pvref] = vrefs[k]
           return
       end
    end
    error("Couldn't find variable to map $pvref to.")
    return
end

# Override the info of the jump variable with the point variable's if any is provided
function _update_point_info(trans_model::JuMP.Model, pvref::InfiniteOpt.PointVariableRef)
    vref = transcription_variable(trans_model, pvref)
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
    if JuMP.is_binary(pvref)
        if JuMP.is_integer(vref)
            JuMP.unset_integer(vref)
        end
        JuMP.set_binary(vref)
    elseif JuMP.is_integer(pvref)
        if JuMP.is_binary(vref)
            JuMP.unset_binary(vref)
        end
        JuMP.set_integer(vref)
    end
    if !(JuMP.start_value(pvref) === NaN)
        JuMP.set_start_value(vref, JuMP.start_value(pvref))
    end
end

# Map point variables to the correct transcribed infinite variable
function _map_point_variables(trans_model::JuMP.Model, inf_model::InfiniteOpt.InfiniteModel)
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

## Helper functions for expanding the measure references in expressions
# GenericAffExpr
function _expand_measures(expr::JuMP.GenericAffExpr{C, <:InfiniteOpt.GeneralVariableRef},
                          trans_model::JuMP.Model) where {C}
    # use a QuadExpr in case measures contain quadratic espressions
    quad = zero(JuMP.GenericQuadExpr{C, InfiniteOpt.GeneralVariableRef})
    quad.aff.constant = expr.constant
    # add the variables to the expr, converting measures into expanded exprs
    for (var, coef) in expr.terms
        if isa(var, InfiniteOpt.MeasureRef)
            func = InfiniteOpt.measure_function(var)
            data = InfiniteOpt.measure_data(var)
            new_func = InfiniteOpt._expand_measure(func, data, trans_model)
            JuMP.add_to_expression!(quad, coef, new_func)
        else
            JuMP.add_to_expression!(quad, coef, var)
        end
    end
    # return a AffExpr if there are no quadratic terms
    if length(quad.terms) == 0
        return quad.aff
    else
        return quad
    end
end

## Helper functions for expanding the measure references in expressions
# GenericQuadExpr
function _expand_measures(expr::JuMP.GenericQuadExpr{C, <:InfiniteOpt.GeneralVariableRef},
                          trans_model::JuMP.Model) where {C}
    quad = zero(JuMP.GenericQuadExpr{C, InfiniteOpt.GeneralVariableRef})
    quad.aff = _expand_measures(expr.aff, trans_model)
    # add the quadratic terms to the expr, converting measures into expanded exprs
    # note that this will error if the expanded forms are not linear or quadratic
    for (pair, coef) in expr.terms
        var_a = pair.a
        var_b = pair.b
        if isa(var_a, InfiniteOpt.MeasureRef)
            func = InfiniteOpt.measure_function(var_a)
            data = InfiniteOpt.measure_data(var_a)
            var_a = InfiniteOpt._expand_measure(func, data, trans_model)
        end
        if isa(var_b, InfiniteOpt.MeasureRef)
            func = InfiniteOpt.measure_function(var_b)
            data = InfiniteOpt.measure_data(var_b)
            var_b = InfiniteOpt._expand_measure(func, data, trans_model)
        end
        JuMP.add_to_expression!(quad, coef * var_a * var_b)
    end
    return quad
end

# Return the support value corresponding to a parameter reference
function _parameter_value(pref::InfiniteOpt.ParameterRef, support::Tuple, prefs::Tuple)
    group = InfiniteOpt.group_id(pref)
    pref_index = findfirst(isequal(group), InfiniteOpt._groups(prefs))
    if isa(prefs[pref_index], InfiniteOpt.ParameterRef)
        return support[pref_index]
    else
        for (k, v) in prefs[pref_index].data
            if v == pref
                return support[pref_index].data[k]
            end
        end
    end
    return
end

## Helper function for mapping infinite variables to jump variables
# FiniteVariableRef
function _map_to_variable(fvref::InfiniteOpt.FiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    return transcription_variable(trans_model, fvref)
end

# InfiniteVariableRef
function _map_to_variable(ivref::InfiniteOpt.InfiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    # reduce support to only include the relavent parameter id groups
    ivref_groups = InfiniteOpt._groups(InfiniteOpt.parameter_refs(ivref))
    support_groups = InfiniteOpt._groups(prefs)
    reduced_support = Tuple(support[findfirst(isequal(group), support_groups)] for group in ivref_groups)
    # find the jump variable associated with the support
    for (index, value) in InfiniteOpt.supports(trans_model, ivref)
        # if isequal(value, reduced_support)
        if value == reduced_support
            return transcription_variable(trans_model, ivref)[index]
        end
    end
    # this error can likely be eliminated
    error("Couldn't find JuMP variable corresponding to $ivref.")
end

# _ReducedInfiniteRef
function _map_to_variable(rvref::InfiniteOpt._ReducedInfiniteRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    # parse the reduced parameters and modify the support to include them
    ivref = rvref.original
    orig_groups = InfiniteOpt._groups(InfiniteOpt.parameter_refs(ivref))
    support_groups = InfiniteOpt._groups(prefs)
    support = [i for i in support]
    prefs = [pref for pref in prefs]
    for (index, value) in rvref.supports
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
                          prefs::Tuple, trans_model::JuMP.Model)
    # find pref in prefs and return associated support value
    value = _parameter_value(pref, support, prefs)
    # this error can likely be eliminated
    isa(value, Nothing) && error("Couldn't find support corresponding to $pref.")
    return value
end

# Return a truncated support dict based on the specified parameter bounds
function _truncate_supports(supports::Dict, prefs::Tuple, bounds::Dict)
    old_support_indices = Dict{Int, Tuple}()
    new_support_indices = Dict{Int, Tuple}()
    eliminate_index = false
    counter = 1
    for (index, support) in supports
        for (pref, set) in bounds
            value = _parameter_value(pref, support, prefs)
            if isa(value, Nothing)
                continue
            end
            if value < set.lower_bound || value > set.upper_bound
                eliminate_index = true
                break
            end
        end
        if !eliminate_index
            old_support_indices[index] = support
            new_support_indices[counter] = support
            counter += 1
        end
        eliminate_index = false
    end
    return new_support_indices, old_support_indices
end

## Convert jump scalar expressions with InfiniteOpt variables into transcribed relations
# PointVariableRef and GlobalVariableRef --> return scalar jump object
function _make_transcription_function(vref::InfiniteOpt.FiniteVariableRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())
    return transcription_variable(trans_model, vref)
end

# InfiniteVariableRef --> return vector of expressions, prefs, and support mapping
function _make_transcription_function(vref::InfiniteOpt.InfiniteVariableRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())
    if length(bounds) != 0
        new_supports, old_supports = _truncate_supports(InfiniteOpt.supports(trans_model, vref),
                                                        InfiniteOpt.parameter_refs(vref),
                                                        bounds)
        ordered_indices = sort(collect(keys(old_supports)))
        return transcription_variable(trans_model, vref)[ordered_indices],
               InfiniteOpt.parameter_refs(vref), new_supports
    else
        return transcription_variable(trans_model, vref), InfiniteOpt.parameter_refs(vref),
               InfiniteOpt.supports(trans_model, vref)
    end
end

# ParameterRef --> return vector of numbers, pref, and support mapping
function _make_transcription_function(pref::InfiniteOpt.ParameterRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())
    support_indices = Dict(i => InfiniteOpt.supports(pref)[i] for i = 1:length(InfiniteOpt.supports(pref)))
    if length(bounds) != 0
        new_supports, old_supports = _truncate_supports(support_indices,
                                                        (pref, ), bounds)
        ordered_indices = sort(collect(keys(old_supports)))
        return InfiniteOpt.supports(pref)[ordered_indices], (pref, ), new_supports
    else
        return InfiniteOpt.supports(pref), (pref, ), support_indices
    end
end

# MeasureRef --> return depends if finite or infinite
function _make_transcription_function(mref::InfiniteOpt.MeasureRef,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())
    func = InfiniteOpt.measure_function(mref)
    data = InfiniteOpt.measure_data(mref)
    new_func = InfiniteOpt._possible_convert(InfiniteOpt.FiniteVariableRef,
                                 InfiniteOpt._expand_measure(func, data, trans_model))
    # will either call finite variable function or general variable function
    return _make_transcription_function(new_func, trans_model, bounds)
end

# GenericAffExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:InfiniteOpt.FiniteVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    # replace finite vars with jump vars
    pairs = [transcription_variable(trans_model, var) => coef for (var, coef) in expr.terms]
    return JuMP.GenericAffExpr(expr.constant,
                               JuMP._new_ordered_dict(JuMP.VariableRef, C, pairs))
end

# GenericQuadExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C, <:InfiniteOpt.FiniteVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    # replace finite vars with jump vars
    pairs = Vector{Pair{JuMP.UnorderedPair{JuMP.VariableRef}, C}}(undef, length(expr.terms))
    counter = 1
    for k in keys(expr.terms)
        a = transcription_variable(trans_model, k.a)
        b = transcription_variable(trans_model, k.b)
        pairs[counter] = JuMP.UnorderedPair(a, b) => expr.terms[k]
        counter += 1
    end
    aff = _make_transcription_function(expr.aff, trans_model)
    return JuMP.GenericQuadExpr(aff, JuMP._new_ordered_dict(JuMP.UnorderedPair{JuMP.VariableRef}, C, pairs))
end

# GenericAffExpr of GeneralVariableRefs --> return vector of numbers, pref, and support mapping
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:InfiniteOpt.GeneralVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 1
        var = collect(keys(expr.terms))[1]
        results = _make_transcription_function(var, trans_model, bounds)
        # results is a tuple if the var is infinite or it is an expr otherwise
        if isa(results, Tuple)
            exprs = expr.terms[var] * results[1] + ones(length(results[1])) * expr.constant
            return exprs, results[2], results[3]
        else
            return expr.terms[var] * results + expr.constant
        end
    end
    # check to see if there are measures and expand them
    is_measure = [var isa InfiniteOpt.MeasureRef for var in InfiniteOpt._all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    # dispatch to quadratic method if the measures contained quadratic terms
    if isa(expr, JuMP.GenericQuadExpr)
        return _make_transcription_function(expr, trans_model, bounds)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = InfiniteOpt._all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
    if length(bounds) != 0
        support_indices =  _truncate_supports(support_indices, prefs, bounds)[1]
    end
    exprs = [zero(JuMP.GenericAffExpr{C, JuMP.VariableRef}) for i = 1:length(support_indices)]
    # make an expression for each support
    for i = 1:length(exprs)
        exprs[i].constant = expr.constant
        for (var, coef) in expr.terms
            # replace each variable with appropriate jump var
            new_var = _map_to_variable(var, support_indices[i], prefs,
                                       trans_model)
            if isa(new_var, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, new_var)
            else
                JuMP.add_to_expression!(exprs[i], coef * new_var)
            end
        end
    end
    return exprs, prefs, support_indices
end

# GenericQuadExpr of GeneralVariableRefs
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C, <:InfiniteOpt.GeneralVariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 0 && length(expr.aff.terms) == 1
      var = collect(keys(expr.aff.terms))[1]
      results = _make_transcription_function(var, trans_model, bounds)
      # results is a tuple if the var is infinite or it is an expr otherwise
      if isa(results, Tuple)
          exprs = expr.aff.terms[var] * results[1] + ones(length(results[1])) * expr.aff.constant
          return exprs, results[2], results[3]
      else
          return expr.aff.terms[var] * results + expr.aff.constant
      end
    end
    # check to see if there are measures and expand them
    is_measure = [var isa MeasureRef for var in InfiniteOpt._all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = InfiniteOpt._all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
    if length(bounds) != 0
        support_indices =  _truncate_supports(support_indices, prefs, bounds)[1]
    end
    exprs = [zero(JuMP.GenericQuadExpr{C, JuMP.VariableRef}) for i = 1:length(support_indices)]
    # make an expression for each support
    for i = 1:length(exprs)
        exprs[i].aff.constant = expr.aff.constant
        for (var, coef) in expr.aff.terms
            # replace each variable with appropriate jump var
            new_var = _map_to_variable(var, support_indices[i], prefs,
                                       trans_model)
            if isa(new_var, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, new_var)
            else
                JuMP.add_to_expression!(exprs[i], coef * new_var)
            end
        end
        for (pair, coef) in expr.terms
            # replace each variable with appropriate jump var
            var_a = _map_to_variable(pair.a, support_indices[i], prefs,
                                     trans_model)
            var_b = _map_to_variable(pair.b, support_indices[i], prefs,
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
    return exprs, prefs, support_indices
end

# GenericAffExpr and GenericQuadExpr of MeasureFiniteVariableRefs --> returns depends on whether finite
function _make_transcription_function(expr::Union{JuMP.GenericAffExpr{C, <:InfiniteOpt.MeasureFiniteVariableRef},
                                      JuMP.GenericQuadExpr{C, <:InfiniteOpt.MeasureFiniteVariableRef}},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    expr = InfiniteOpt._possible_convert(FiniteVariableRef, _expand_measures(expr, trans_model))
    return _make_transcription_function(expr, trans_model, bounds)
end

# Empty jump variable expr (for constraints of form number <= number)
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, JuMP.VariableRef},
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict()) where {C}
    return expr
end

# Fall back function for other jump objects
function _make_transcription_function(expr::JuMP.AbstractJuMPScalar,
                                      trans_model::JuMP.Model,
                                      bounds::Dict = Dict())
    type = typeof(expr)
    error("Unsupported transcription of expression of type $type.")
    return
end

## Construct the objective and error is contains non finite variables
function _set_objective(trans_model::JuMP.Model, inf_model::InfiniteOpt.InfiniteModel)
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
function _set_mapping(icref::InfiniteOpt.InfiniteConstraintRef, crefs::Vector{JuMP.ConstraintRef})
    transcription_data(JuMP.owner_model(crefs[1])).infinite_to_constr[icref] = crefs
    return
end

# MeasureConstraintRef (infinite)
function _set_mapping(mcref::InfiniteOpt.MeasureConstraintRef, crefs::Vector{JuMP.ConstraintRef})
    transcription_data(JuMP.owner_model(crefs[1])).measure_to_constr[mcref] = crefs
    return
end

# MeasureConstraintRef (finite)
function _set_mapping(mcref::InfiniteOpt.MeasureConstraintRef, cref::JuMP.ConstraintRef)
    transcription_data(JuMP.owner_model(cref)).measure_to_constr[mcref] = [cref]
    return
end

# FiniteConstraintRef
function _set_mapping(fcref::InfiniteOpt.FiniteConstraintRef, cref::JuMP.ConstraintRef)
    transcription_data(JuMP.owner_model(cref)).finite_to_constr[fcref] = cref
    return
end

## Define helper functions for setting constraint supports
# InfiniteConstraintRef
function _set_supports(trans_model::JuMP.Model, icref::InfiniteOpt.InfiniteConstraintRef,
                       supports::Dict{})
    transcription_data(trans_model).infconstr_to_supports[icref] = supports
    return
end

# MeasureConstraintRef
function _set_supports(trans_model::JuMP.Model, mcref::InfiniteOpt.MeasureConstraintRef,
                       supports::Dict{})
    transcription_data(trans_model).measconstr_to_supports[mcref] = supports
    return
end

## Define helper functions for setting constraint parameter reference tuples
# InfiniteConstraintRef
function _set_parameter_refs(trans_model::JuMP.Model, icref::InfiniteOpt.InfiniteConstraintRef,
                             prefs::Tuple)
    transcription_data(trans_model).infconstr_to_params[icref] = prefs
    return
end

# MeasureConstraintRef
function _set_parameter_refs(trans_model::JuMP.Model, mcref::InfiniteOpt.MeasureConstraintRef,
                             prefs::Tuple)
    transcription_data(trans_model).measconstr_to_params[icref] = prefs
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
            # Tuple results correspond to transcribed infinite constraints
            if isa(results, Tuple)
                crefs = Vector{JuMP.ConstraintRef}(undef, length(results[1]))
                for i = 1:length(crefs)
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
                con = JuMP.build_constraint(error, results, constr.set)
                cref = JuMP.add_constraint(trans_model, con, root_name)
                _set_mapping(icref, cref)
            end
        end
    end
    return
end

## Helper functions for mapping InfiniteOpt var info constrs to the transcribed ones
# FiniteVariableRef
function _map_info_constraints(fvref::InfiniteOpt.FiniteVariableRef, trans_model::JuMP.Model)
    vref = transcription_variable(trans_model, fvref)
    # Check if both variables have a constraint and map if they do
    if JuMP.has_lower_bound(fvref) && JuMP.has_lower_bound(vref)
        _set_mapping(JuMP.LowerBoundRef(fvref), JuMP.LowerBoundRef(vref))
    end
    if JuMP.has_upper_bound(fvref) && JuMP.has_upper_bound(vref)
        _set_mapping(JuMP.UpperBoundRef(fvref), JuMP.UpperBoundRef(vref))
    end
    if JuMP.is_fixed(fvref) && JuMP.is_fixed(vref)
        _set_mapping(JuMP.FixRef(fvref), JuMP.FixRef(vref))
    end
    if JuMP.is_integer(fvref) && JuMP.is_integer(vref)
        _set_mapping(JuMP.IntegerRef(fvref), JuMP.IntegerRef(vref))
    end
    if JuMP.is_binary(fvref) && JuMP.is_binary(vref)
        _set_mapping(JuMP.BinaryRef(fvref), JuMP.BinaryRef(vref))
    end
    return
end

# InfiniteVariableRef
function _map_info_constraints(ivref::InfiniteOpt.InfiniteVariableRef, trans_model::JuMP.Model)
    vrefs = transcription_variable(trans_model, ivref)
    # Check if both variables have a constraint and map if they do
    if JuMP.has_lower_bound(ivref)
        crefs = JuMP.ConstraintRef[]
        support_indices = Dict{Int, Tuple}()
        counter = 1
        for i = 1:length(vrefs)
            if JuMP.has_lower_bound(vrefs[i])
                push!(crefs, JuMP.LowerBoundRef(vrefs[i]))
                support_indices[counter] = InfiniteOpt.supports(trans_model, ivref)[i]
                counter += 1
            end
        end
        _set_mapping(JuMP.LowerBoundRef(ivref), crefs)
        _set_parameter_refs(trans_model, JuMP.LowerBoundRef(ivref),
                            InfiniteOpt.parameter_refs(ivref))
        _set_supports(trans_model, JuMP.LowerBoundRef(ivref), support_indices)
    end
    if JuMP.has_upper_bound(ivref)
        crefs = JuMP.ConstraintRef[]
        support_indices = Dict{Int, Tuple}()
        counter = 1
        for i = 1:length(vrefs)
            if JuMP.has_upper_bound(vrefs[i])
                push!(crefs, JuMP.UpperBoundRef(vrefs[i]))
                support_indices[counter] = InfiniteOpt.supports(trans_model, ivref)[i]
                counter += 1
            end
        end
        _set_mapping(JuMP.UpperBoundRef(ivref), crefs)
        _set_parameter_refs(trans_model, JuMP.UpperBoundRef(ivref),
                            InfiniteOpt.parameter_refs(ivref))
        _set_supports(trans_model, JuMP.UpperBoundRef(ivref), support_indices)
    end
    if JuMP.is_fixed(ivref)
        crefs = JuMP.ConstraintRef[]
        support_indices = Dict{Int, Tuple}()
        counter = 1
        for i = 1:length(vrefs)
            if JuMP.is_fixed(vrefs[i])
                push!(crefs, JuMP.FixRef(vrefs[i]))
                support_indices[counter] = InfiniteOpt.supports(trans_model, ivref)[i]
                counter += 1
            end
        end
        _set_mapping(JuMP.FixRef(ivref), crefs)
        _set_parameter_refs(trans_model, JuMP.FixRef(ivref),
                            InfiniteOpt.parameter_refs(ivref))
        _set_supports(trans_model, JuMP.FixRef(ivref), support_indices)
    end
    if JuMP.is_integer(ivref)
        crefs = JuMP.ConstraintRef[]
        support_indices = Dict{Int, Tuple}()
        counter = 1
        for i = 1:length(vrefs)
            if JuMP.is_integer(vrefs[i])
                push!(crefs, JuMP.IntegerRef(vrefs[i]))
                support_indices[counter] = InfiniteOpt.supports(trans_model, ivref)[i]
                counter += 1
            end
        end
        _set_mapping(JuMP.IntegerRef(ivref), crefs)
        _set_parameter_refs(trans_model, JuMP.IntegerRef(ivref),
                            InfiniteOpt.parameter_refs(ivref))
        _set_supports(trans_model, JuMP.IntegerRef(ivref), support_indices)
    end
    if JuMP.is_binary(ivref)
        crefs = JuMP.ConstraintRef[]
        support_indices = Dict{Int, Tuple}()
        counter = 1
        for i = 1:length(vrefs)
            if JuMP.is_binary(vrefs[i])
                push!(crefs, JuMP.BinaryRef(vrefs[i]))
                support_indices[counter] = InfiniteOpt.supports(trans_model, ivref)[i]
                counter += 1
            end
        end
        _set_mapping(JuMP.BinaryRef(ivref), crefs)
        _set_parameter_refs(trans_model, JuMP.BinaryRef(ivref),
                            InfiniteOpt.parameter_refs(ivref))
        _set_supports(trans_model, JuMP.BinaryRef(ivref), support_indices)
    end
    return
end

## Map the variable info constraints between the two models
function _map_variable_info_constraints(trans_model::JuMP.Model,
                                        inf_model::InfiniteOpt.InfiniteModel)
    for (index, var) in inf_model.vars
        ivref = InfiniteOpt._make_variable_ref(inf_model, index)
        if InfiniteOpt.is_used(ivref)
            _map_info_constraints(ivref, trans_model)
        end
    end
    return
end

"""
    TranscriptionModel(model::InfiniteModel, args...)
Return a transcribed version of the infinite model.
"""
function TranscriptionModel(inf_model::InfiniteOpt.InfiniteModel,
                            optimizer_factory::Union{JuMP.OptimizerFactory, Nothing} = nothing;
                            kwargs...)
    if isa(optimizer_factory, Nothing)
        trans_model = TranscriptionModel(; kwargs...)
    else
        trans_model = TranscriptionModel(optimizer_factory; kwargs...)
    end

    @time _initialize_global_variables(trans_model, inf_model)
    @time _initialize_infinite_variables(trans_model, inf_model)
    @time _map_point_variables(trans_model, inf_model)
    if JuMP.objective_sense(inf_model) != MOI.FEASIBILITY_SENSE
        @time _set_objective(trans_model, inf_model)
    end
    @time _map_variable_info_constraints(trans_model, inf_model)
    # TODO optimize performance --> bottlenecked with _make_transcription_function
    @time _set_constraints(trans_model, inf_model)
    return trans_model
end
