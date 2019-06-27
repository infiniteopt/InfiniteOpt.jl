"""
    TranscriptionModel(args...)
Return a JuMP `Model` with `TranscriptionData` included in the extension
data field.
"""
function TranscriptionModel(args...)
    model = JuMP.Model(args...)
    model.ext[:TransData] = TranscriptionData()
    return model
end

"""
    is_transcription_model(model::JuMP.Model)
Return true if JuMP `Model` is a `TranscriptionModel` or false otherwise.
"""
function is_transcription_model(model::JuMP.Model)
    return haskey(model.ext, :TransData)
end

"""
    transcription_data(model::JuMP.Model)
Return the `TranscriptionData` from a `TranscriptionModel`. Error if it is not a
a `TranscriptionModel`.
"""
function transcription_data(model::JuMP.Model)
    !is_transcription_model(model) && error("Model is not a transcription model.")
    return model.ext[:TransData]
end

"""
    transcription_variable(model::JuMP.Model, vref::GlobalVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::GlobalVariableRef)
    !haskey(transcription_data(model).global_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).global_to_var[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::InfiniteVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::InfiniteVariableRef)
    !haskey(transcription_data(model).infinite_to_vars, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infinite_to_vars[vref]
end

"""
    transcription_variable(model::JuMP.Model, vref::PointVariableRef)
Return the transcribed variable reference corresponding to `vref`.
"""
function transcription_variable(model::JuMP.Model, vref::PointVariableRef)
    !haskey(transcription_data(model).point_to_var, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).point_to_var[vref]
end

"""
    supports(model::JuMP.Model, vref::InfiniteVariableRef)
Return the supports associated `vref` in the transcribed model.
"""
function supports(model::JuMP.Model, vref::InfiniteVariableRef)
    !haskey(transcription_data(model).infinite_to_supports, vref) && error("Variable reference $vref not used in transcription model.")
    return transcription_data(model).infinite_to_supports[vref]
end

# Make jump variables and a dict mapping global vars to jump vars
function _initialize_global_variables(trans_model::JuMP.Model,
                                      inf_model::InfiniteModel)
    # search inf_model for global vars and make a jump var for one that is used
    for (index, var) in inf_model.vars
        if isa(var, GlobalVariable)
            gvref = GlobalVariableRef(inf_model, index)
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
        support_list[i] = supports(prefs[i])
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
                                        inf_model::InfiniteModel)
    # search inf_model for infinite vars and make a jump var for all of its supports
    for (index, var) in inf_model.vars
        if isa(var, InfiniteVariable)
            ivref = InfiniteVariableRef(inf_model, index)
            if is_used(ivref)
                prefs = parameter_refs(ivref)
                transcription_data(trans_model).infinite_to_supports[ivref] = _make_support_indices(prefs)
                vrefs = Vector{JuMP.VariableRef}(undef, length(keys(supports(trans_model, ivref))))
                for i = 1:length(vrefs)
                    # TODO Perhaps add different naming options...
                    name = _root_name(ivref)
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
function _update_point_mapping(trans_model::JuMP.Model, pvref::PointVariableRef,
                               ivref::InfiniteVariableRef, support::Tuple)
    for (k, v) in supports(trans_model, ivref)
       if v == support
           vrefs = transcription_variable(trans_model, ivref)
           transcription_data(trans_model).point_to_var[pvref] = vrefs[k]
           break
       end
    end
    return
end

# Override the info of the jump variable with the point variable's if any is provided
function _update_point_info(trans_model::JuMP.Model, pvref::PointVariableRef)
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
function _map_point_variables(trans_model::JuMP.Model, inf_model::InfiniteModel)
    # search inf_model for point vars and map them to jump vars if they are used
    for (index, var) in inf_model.vars
        if isa(var, PointVariable)
            pvref = PointVariableRef(inf_model, index)
            if is_used(pvref)
                ivref = infinite_variable_ref(pvref)
                support = parameter_values(pvref)
                _update_point_mapping(trans_model, pvref, ivref, support)
                _update_point_info(trans_model, pvref)
            end
        end
    end
    return
end

## Helper functions for expanding the measure references in expressions
# GenericAffExpr
function _expand_measures(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                          trans_model::JuMP.Model) where {C}
    # use a QuadExpr in case measures contain quadratic espressions
    quad = zero(JuMP.GenericQuadExpr{C,GeneralVariableRef})
    quad.aff.constant = expr.constant
    # add the variables to the expr, converting measures into expanded exprs
    for (var, coef) in expr.terms
        if isa(var, MeasureRef)
            func = measure_function(var)
            data = measure_data(var)
            new_func = _expand_measure(func, data, trans_model)
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
function _expand_measures(expr::JuMP.GenericQuadExpr{C, <:GeneralVariableRef},
                          trans_model::JuMP.Model) where {C}
    quad = zero(JuMP.GenericQuadExpr{C, GeneralVariableRef})
    quad.aff = _expand_measures(expr.aff, trans_model)
    # add the quadratic terms to the expr, converting measures into expanded exprs
    # note that this will error if the expanded forms are not linear or quadratic
    for (pair, coef) in expr.terms
        var_a = pair.a
        var_b = pair.b
        if isa(var_a, MeasureRef)
            func = measure_function(var_a)
            data = measure_data(var_a)
            var_a = _expand_measure(func, data, trans_model)
        end
        if isa(var_b, MeasureRef)
            func = measure_function(var_b)
            data = measure_data(var_b)
            var_b = _expand_measure(func, data, trans_model)
        end
        JuMP.add_to_expression!(quad, coef * var_a * var_b)
    end
    return quad
end

## Helper function for mapping infinite variables to jump variables
# FiniteVariableRef
function _map_to_variable(fvref::FiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    return transcription_variable(trans_model, fvref)
end

# InfiniteVariableRef
function _map_to_variable(ivref::InfiniteVariableRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    # reduce support to only include the relavent parameter id groups
    ivref_groups = _groups(parameter_refs(ivref))
    support_groups = _groups(prefs)
    reduced_support = Tuple(support[findfirst(isequal(group), support_groups)] for group in ivref_groups)
    # find the jump variable associated with the support
    for (index, value) in supports(trans_model, ivref)
        if value == reduced_support
            return transcription_variable(trans_model, ivref)[index]
        end
    end
    # this error can likely be eliminated
    error("Couldn't find JuMP variable corresponding to $ivref.")
end

# _ReducedInfiniteRef
function _map_to_variable(rvref::_ReducedInfiniteRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    # parse the reduced parameters and modify the support to include them
    ivref = rvref.original
    orig_groups = _groups(parameter_refs(ivref))
    support_groups = _groups(prefs)
    support = [i for i in support]
    prefs = [pref for pref in prefs]
    for (index, value) in rvref.supports
        if orig_groups[index] in support_groups
            support[findfirst(isequal(orig_groups[index]), support_groups)] = value
        else
            push!(support, value)
            push!(prefs, parameter_refs(ivref)[index])
        end
    end
    prefs = Tuple(pref for pref in prefs)
    support = Tuple(i for i in support)
    # call the infinite variable function with the updated support and prefs
    return _map_to_variable(ivref, support, prefs, trans_model)
end

# ParameterRef
function _map_to_variable(pref::ParameterRef, support::Tuple,
                          prefs::Tuple, trans_model::JuMP.Model)
    # find pref in prefs and return associated support value
    group = group_id(pref)
    pref_index = findfirst(isequal(group), _groups(prefs))
    if isa(prefs[pref_index], ParameterRef)
        return support[pref_index]
    else
        for (k, v) in support[pref_index].data
            if v == pref
                return support[pref_index].data[k]
            end
        end
    end
    # this error can likely be eliminated
    error("Couldn't find support corresponding to $pref.")
end

## Convert jump scalar expressions with InfOpt variables into transcribed relations
# PointVariableRef and GlobalVariableRef --> return scalar jump object
function _make_transcription_function(vref::FiniteVariableRef,
                                      trans_model::JuMP.Model)
    return transcription_variable(trans_model, vref)
end

# InfiniteVariableRef --> return vector of expressions, prefs, and support mapping
function _make_transcription_function(vref::InfiniteVariableRef,
                                      trans_model::JuMP.Model)
    return transcription_variable(trans_model, vref), parameter_refs(vref),
           supports(trans_model, vref)
end

# ParameterRef --> return vector of numbers, pref, and support mapping
function _make_transcription_function(pref::ParameterRef,
                                      trans_model::JuMP.Model)
    return supports(pref), (pref, ),
           Dict(i => supports(pref)[i] for i = 1:length(supports(pref)))
end

# MeasureRef --> return depends if finite or infinite
function _make_transcription_function(mref::MeasureRef,
                                      trans_model::JuMP.Model)
    func = measure_function(mref)
    data = measure_data(mref)
    new_func = _possible_convert(FiniteVariableRef,
                                 _expand_measure(func, data, trans_model))
    # will either call finite variable function or general variable function
    return _make_transcription_function(new_func, trans_model)
end

# GenericAffExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:FiniteVariableRef},
                                      trans_model::JuMP.Model) where {C}
    # replace finite vars with jump vars
    pairs = [transcription_variable(trans_model, var) => coef for (var, coef) in expr.terms]
    return JuMP.GenericAffExpr(expr.constant,
                               JuMP._new_ordered_dict(JuMP.VariableRef, C, pairs))
end

# GenericQuadExpr of FiniteVariableRefs --> return scalar jump object
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C, <:FiniteVariableRef},
                                      trans_model::JuMP.Model) where {C}
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
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                                      trans_model::JuMP.Model) where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 1
        var = collect(keys(expr.terms))[1]
        results = _make_transcription_function(var, trans_model)
        # results is a tuple if the var is infinite or it is an expr otherwise
        if isa(results, Tuple)
            exprs = expr.terms[var] * results[1] + ones(length(results[1])) * expr.constant
            return exprs, results[2], results[3]
        else
            return expr.terms[var] * results + expr.constant
        end
    end
    # check to see if there are measures and expand them
    is_measure = [var isa MeasureRef for var in _all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    # dispatch to quadratic method if the measures contained quadratic terms
    if isa(expr, JuMP.GenericQuadExpr)
        return _make_transcription_function(expr, trans_model)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = _all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
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
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C, <:GeneralVariableRef},
                                      trans_model::JuMP.Model) where {C}
    # check if there is only 1 var to dispatch to that transcription method
    if length(expr.terms) == 0 && length(expr.aff.terms) == 1
      var = collect(keys(expr.aff.terms))[1]
      results = _make_transcription_function(var, trans_model)
      # results is a tuple if the var is infinite or it is an expr otherwise
      if isa(results, Tuple)
          exprs = expr.aff.terms[var] * results[1] + ones(length(results[1])) * expr.aff.constant
          return exprs, results[2], results[3]
      else
          return expr.aff.terms[var] * results + expr.aff.constant
      end
    end
    # check to see if there are measures and expand them
    is_measure = [var isa MeasureRef for var in _all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    # determine the common set of prefs and make all of the support combos
    prefs = _all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
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
function _make_transcription_function(expr::Union{JuMP.GenericAffExpr{C, <:MeasureFiniteVariableRef},
                                      JuMP.GenericQuadExpr{C, <:MeasureFiniteVariableRef}},
                                      trans_model::JuMP.Model) where {C}
    expr = _possible_convert(FiniteVariableRef, _expand_measures(expr, trans_model))
    return _make_transcription_function(expr, trans_model)
end

# Empty jump variable expr (for constraints of form number <= number)
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, VariableRef},
                                      trans_model::JuMP.Model) where {C}
    return expr
end

# Fall back function for other jump objects
function _make_transcription_function(expr::JuMP.AbstractJuMPScalar,
                                      trans_model::JuMP.Model)
    type = typeof(expr)
    error("Unsupported transcription of expression of type $type.")
end

## Construct the objective and error is contains non finite variables
function _set_objective(trans_model::JuMP.Model, inf_model::InfiniteModel)
    trans_obj, = _make_transcription_function(JuMP.objective_function(inf_model), trans_model)
    isa(trans_obj, Vector) && error("Objective is not finite, ensure all " *
                                    "infinite variables/parameters in measures " *
                                    "are evaluated completely.")
    JuMP.set_objective(trans_model, JuMP.objective_sense(inf_model), trans_obj)
    return
end

# TODO Make constraint intializer

"""
    generate_transcribed_model(model::InfiniteModel)
Return a transcribed version of the model.
"""
function generate_transcription_model(inf_model::InfiniteModel)
    trans_model = TranscriptionModel()
    _initialize_global_variables(trans_model, inf_model)
    _initialize_infinite_variables(trans_model, inf_model)
    _map_point_variables(trans_model, inf_model)
    if JuMP.objective_sense(inf_model) != MOI.FEASIBILITY_SENSE
        _set_objective(trans_model, inf_model)
    end
    # TODO Add constraints if there are constraints
    return trans_model
end
