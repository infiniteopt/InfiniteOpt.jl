# Make jump variables and return a dict mapping global vars to jump vars
function _initialize_global_variables(trans_model::JuMP.Model,
                                      inf_model::InfiniteModel)
    global_to_var = Dict{GlobalVariableRef, JuMP.VariableRef}()
    for (index, var) in inf_model.vars
        if isa(var, GlobalVariable)
            gvref = GlobalVariableRef(inf_model, index)
            if is_used(gvref)
                vref = JuMP.add_variable(trans_model, JuMP.ScalarVariable(var.info),
                                         JuMP.name(gvref))
                global_to_var[gvref] = vref
            end
        end
    end
    return global_to_var
end

# Return a vector of arrays containing the supports
function _list_supports(prefs::Tuple)
    support_list = Vector{Array}(undef, length(prefs))
    for i = 1:length(prefs)
        support_list[i] = supports(prefs[i])
    end
end

# Make an index mapping for parameter support combinations for an infinite variable
function _make_support_indices(vref::InfiniteVariableRef)
    support_indices = Dict{Int, Tuple}()
    prefs = get_parameter_refs(vref)

end

# Make jump variables and return a dict mapping infinite/point vars to jump vars
function _initialize_infinite_variables(trans_model::JuMP.Model,
                                        inf_model::InfiniteModel)
    infinite_to_var = Dict{GlobalVariableRef, Vector{JuMP.VariableRef}}()
    for (index, var) in inf_model.vars
        if isa(var, InfiniteVariable)

        end
    end
    return
end

"""
    generate_transcribed_model(model::InfiniteModel)
Return a transcribed version of the model.
"""
function generate_transcribed_model(inf_model::InfiniteModel)
    trans_model = JuMP.Model()
    global_to_var = _initialize_global_variables(trans_model, inf_model)
    return trans_model
end
