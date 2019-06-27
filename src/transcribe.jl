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
    for (index, var) in inf_model.vars
        if isa(var, InfiniteVariable)
            ivref = InfiniteVariableRef(inf_model, index)
            if is_used(ivref)
                prefs = parameter_refs(ivref)
                transcription_data(trans_model).infinite_to_supports[ivref] = _make_support_indices(prefs)
                vrefs = Vector{JuMP.VariableRef}(undef, length(keys(supports(trans_model, ivref))))
                for i = 1:length(vrefs)
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

# Map point variables to the correct transcribed infinite variable
function _map_point_variables(trans_model::JuMP.Model, inf_model::InfiniteModel)
    for (index, var) in inf_model.vars
        if isa(var, PointVariable)
            pvref = PointVariableRef(inf_model, index)
            if is_used(pvref)
                ivref = infinite_variable_ref(pvref)
                support = parameter_values(pvref)
                _update_point_mapping(trans_model, pvref, ivref, support)
            end
        end
    end
    return
end

## Helper functions for expanding the measure references in expressions
# GenericAffExpr
function _expand_measures(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                          trans_model::JuMP.Model) where {C}
    quad = zero(JuMP.GenericQuadExpr{C,GeneralVariableRef})
    quad.aff.constant = expr.constant
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
    error("Couldn't find support corresponding to $pref.")
end

## Convert jump scalar expressions with InfOpt variables into transcribed relations
# InfiniteVariableRef, PointVariableRef, GlobalVariableRef
function _make_transcription_function(vref::InfOptVariableRef,
                                      trans_model::JuMP.Model)
    return transcription_variable(trans_model, vref), parameter_refs(ivref),
           supports(trans_model, ivref)
end

# ParameterRef
function _make_transcription_function(pref::ParameterRef,
                                      trans_model::JuMP.Model)
    return supports(pref), (pref, ),
           Dict(i => supports(pref)[i] for i = 1:length(supports(pref)))
end

# MeasureRef
function _make_transcription_function(mref::MeasureRef,
                                      trans_model::JuMP.Model)
    func = measure_function(mref)
    data = measure_data(mref)
    new_func = _possible_convert(FiniteVariableRef,
                                 _expand_measure(func, data, trans_model))
    return _make_transcription_function(new_func, trans_model)
end

# GenericAffExpr of FiniteVariableRefs
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:FiniteVariableRef},
                                      trans_model::JuMP.Model) where {C}
    pairs = [transcription_variable(trans_model, var) => coef for (var, coef) in expr.terms]
    return JuMP.GenericAffExpr(expr.constant,
                               JuMP._new_ordered_dict(JuMP.VariableRef, C, pairs))
end

# GenericQuadExpr of FiniteVariableRefs
function _make_transcription_function(expr::JuMP.GenericQuadExpr{C, <:FiniteVariableRef},
                                      trans_model::JuMP.Model) where {C}
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

# GenericAffExpr of GeneralVariableRefs
function _make_transcription_function(expr::JuMP.GenericAffExpr{C, <:GeneralVariableRef},
                                      trans_model::JuMP.Model) where {C}
    is_measure = [var isa MeasureRef for var in _all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    if isa(expr, JuMP.GenericQuadExpr)
        return _make_transcription_function(expr, trans_model)
    end
    prefs = _all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
    exprs = [zero(JuMP.GenericAffExpr{C, JuMP.VariableRef}) for i = 1:length(support_indices)]
    for i = 1:length(exprs)
        exprs[i].constant = expr.constant
        for (var, coef) in expr.terms
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
    is_measure = [var isa MeasureRef for var in _all_function_variables(expr)]
    if any(is_measure)
        expr = _expand_measures(expr, trans_model)
    end
    prefs = _all_parameter_refs(expr)
    support_indices = _make_support_indices(prefs)
    exprs = [zero(JuMP.GenericQuadExpr{C, JuMP.VariableRef}) for i = 1:length(support_indices)]
    for i = 1:length(exprs)
        exprs[i].aff.constant = expr.aff.constant
        for (var, coef) in expr.aff.terms
            new_var = _map_to_variable(var, support_indices[i], prefs,
                                       trans_model)
            if isa(new_var, JuMP.VariableRef)
                JuMP.add_to_expression!(exprs[i], coef, new_var)
            else
                JuMP.add_to_expression!(exprs[i], coef * new_var)
            end
        end
        for (pair, coef) in expr.terms
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

# GenericAffExpr and GenericQuadExpr of MeasureFiniteVariableRefs
function _make_transcription_function(expr::Union{JuMP.GenericAffExpr{C, MeasureFiniteVariableRef},
                                      JuMP.GenericQuadExpr{C, MeasureFiniteVariableRef}},
                                      trans_model::JuMP.Model) where {C}
    expr = _possible_convert(FiniteVariableRef, _expand_measures(expr, trans_model))
    return _make_transcription_function(expr, trans_model)
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
