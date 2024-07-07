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
    if isnothing(object) 
        error("Invalid finite variable reference, cannot find ",
              "corresponding variable in the model. This is likely ",
              "caused by using the reference of a deleted variable.")
    end
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
## Process the variable to use the correct info 
# All Floats 
function _process_scalar_var(
    v::V
    )::V where {V <: JuMP.ScalarVariable{Float64, Float64, Float64, Float64}}
    return v
end

# Other
function _process_scalar_var(
    v::JuMP.ScalarVariable
    )::JuMP.ScalarVariable{Float64, Float64, Float64, Float64}
    return JuMP.ScalarVariable(_make_float_info(v.info))
end

"""
    JuMP.add_variable(model::InfiniteModel, var::JuMP.ScalarVariable,
                      [name::String = ""])::GeneralVariableRef

Extend the `JuMP.add_variable` function to accomodate finite (scalar) variable 
types. Adds a variable to an infinite model `model` and returns a 
[`GeneralVariableRef`](@ref). Primarily intended to be an internal function of 
the constructor macro `@variable`. However, it can be used in combination with
`JuMP.build_variable` to add finite variables to an infinite model object.

**Examples**
```julia-repl
julia> f_var = build_variable(error, info);

julia> fvref = add_variable(m, f_var, "var_name")
var_name
```
"""
function JuMP.add_variable(
    model::InfiniteModel, 
    var::JuMP.ScalarVariable,
    name::String = ""
    )::GeneralVariableRef
    new_var = _process_scalar_var(var)
    data_object = VariableData(new_var, name)
    vindex = _add_data_object(model, data_object)
    vref = FiniteVariableRef(model, vindex)
    gvref = GeneralVariableRef(model, vindex)
    _set_info_constraints(new_var.info, gvref, vref)
    model.name_to_var = nothing
    return gvref
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
