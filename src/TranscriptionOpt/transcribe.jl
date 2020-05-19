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
    return collect(keys(supps))
end

# Build the parameter supports
function set_parameter_supports(model::JuMP.Model,
                                inf_model::InfiniteOpt.InfiniteModel)::Nothing
    param_indices = InfiniteOpt._param_object_indices(inf_model)
    supps = Tuple(_collected_supports(inf_model, i) for i in param_indices)
    transcription_data(model).supports = supps
    return
end


################################################################################
#                        VARIABLE INITIALIZATION METHODS
################################################################################


################################################################################
#                       TRANSCRIPTION EXPRESSION METHODS
################################################################################


################################################################################
#                         MEASURE TRANSCRIPTION METHODS
################################################################################


################################################################################
#                        CONSTRAINT TRANSCRIPTION METHODS
################################################################################


################################################################################
#                      INFINITEMODEL TRANSCRIPTION METHODS
################################################################################
