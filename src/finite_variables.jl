################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::FiniteVariableIndex
    )::FiniteVariableRef
    return FiniteVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:JuMP.ScalarVariable}
    )::FiniteVariableIndex
    return MOIUC.add_item(model.finite_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel, 
    ::Type{FiniteVariable}
    )::MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}
    return model.finite_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(
    vref::FiniteVariableRef
    )::MOIUC.CleverDict{FiniteVariableIndex, VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}}
    return JuMP.owner_model(vref).finite_vars
end

# Extend _data_object
function _data_object(
    vref::FiniteVariableRef
    )::VariableData{JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}
  object = get(_data_dictionary(vref), JuMP.index(vref), nothing)
  object === nothing && error("Invalid finite variable reference, cannot find " *
                        "corresponding variable in the model. This is likely " *
                        "caused by using the reference of a deleted variable.")
  return object
end

# Extend _core_variable_object
function _core_variable_object(
    vref::FiniteVariableRef
    )::JuMP.ScalarVariable{Float64, Float64, Float64, Float64}
    return _data_object(vref).variable
end

################################################################################
#                             DEFINTION METHODS
################################################################################
# Define _check_and_make_variable_ref (used by JuMP.add_variable)
function _check_and_make_variable_ref(
    model::InfiniteModel,
    v::JuMP.ScalarVariable{Float64, Float64, Float64, Float64},
    name::String
    )::FiniteVariableRef
    data_object = VariableData(v, name)
    vindex = _add_data_object(model, data_object)
    vref = FiniteVariableRef(model, vindex)
    return vref
end

# Define for scalar variables without strictly Float64 info
function _check_and_make_variable_ref(
    model::InfiniteModel,
    v::JuMP.ScalarVariable,
    name::String
    )::FiniteVariableRef
    new_v = JuMP.ScalarVariable(_make_float_info(v.info))
    return _check_and_make_variable_ref(model, new_v, name)
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for finite variables
function _update_variable_info(
    vref::FiniteVariableRef,
    info::JuMP.VariableInfo
    )::Nothing
    _set_core_variable_object(vref, JuMP.ScalarVariable(info))
    return
end

################################################################################
#                                 DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::FiniteVariableRef)::Nothing
    # remove variable info constraints associated with vref
    _delete_info_constraints(vref)
    return
end
