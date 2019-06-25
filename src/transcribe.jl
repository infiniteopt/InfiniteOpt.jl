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
function _make_support_indices(ivref::InfiniteVariableRef)
    support_indices = Dict{Int, Tuple}()
    support_list = _list_supports(parameter_refs(ivref))
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
                transcription_data(trans_model).infinite_to_supports[ivref] = _make_support_indices(ivref)
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

# Convert jump scalar expressions with InfOpt variables into transcribed relations
function _make_transcription_function(vref::InfOptVariableRef,
                                      trans_model::JuMP.Model)
    return transcription_variable(trans_model, vref)
end

function _make_transcription_function(pref::ParameterRef,
                                      trans_model::JuMP.Model)
    return supports(pref)
end

# TODO Implement for jump expressions
function _make_transcription_function(expr::JuMP.GenericAffExpr{V, <:GeneralVariableRef},
                                      trans_model::JuMP.Model) where {V}
    constant = expr.constant
    pairs = [transcription_variable(trans_model, k) => expr.terms[k] for k in keys(expr.terms)]
    return JuMP.GenericAffExpr(expr.constant,
                               JuMP._new_ordered_dict(JuMP.VariableRef, V, pairs))
end

function _make_transcription_function(expr::JuMP.GenericQuadExpr{V, <:FiniteVariableRef},
                                      trans_model::JuMP.Model) where {V}
    pairs = Vector{Pair{JuMP.UnorderedPair{JuMP.VariableRef}, V}}(undef, length(expr.terms))
    counter = 1
    for k in keys(expr.terms)
        a = transcription_variable(trans_model, k.a)
        b = transcription_variable(trans_model, k.b)
        pairs[counter] = JuMP.UnorderedPair(a, b) => expr.terms[k]
        counter += 1
    end
    aff = _make_transcription_function(expr.aff, trans_model)
    return JuMP.GenericQuadExpr(aff, JuMP._new_ordered_dict(JuMP.UnorderedPair{JuMP.VariableRef}, V, pairs))
end

function _make_transcription_function(expr::JuMP.AbstractJuMPScalar,
                                      trans_model::JuMP.Model)
    type = typeof(expr)
    error("Unsupported transcription of expression of type $type.")
end

function _make_transcription_function(mref::MeasureRef,
                                      trans_model::JuMP.Model)
    func = measure_function(mref)
    data = measure_data(mref)
    new_func =  _expand_measure_function(func, data, trans_model)
    # TODO transcribe properly
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
    # TODO Add objective if there is one
    # TODO Add constraints if there are constraints
    return trans_model
end
