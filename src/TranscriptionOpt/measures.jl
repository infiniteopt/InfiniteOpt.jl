"""
    InfiniteOpt.add_point_variable(model::JuMP.Model,
                                   var::InfiniteOpt.PointVariable,
                                   key::Val{:TransData}
                                   )::InfiniteOpt.GeneralVariableRef

Make a `PointVariableRef` and map it to the appropriate transcription variable
and return the `GeneralVariableRef`. This is an extension of
[`add_point_variable`](@ref InfiniteOpt.add_point_variable(::JuMP.Model,::Any,::Any, ::Any))
for `TranscriptionOpt`.
"""
function InfiniteOpt.add_point_variable(
    model::JuMP.Model,
    ivref::InfiniteOpt.GeneralVariableRef,
    support::Vector{Float64},
    ::Val{:TransData}
    )::InfiniteOpt.GeneralVariableRef
    # check if an internal variable was already created
    data = transcription_data(model)
    internal_vref = get(data.point_lookup, (ivref, support), nothing)
    if !isnothing(internal_vref)
        return internal_vref
    end
    # check if whether the infinite model already has one, else create one
    inf_model = JuMP.owner_model(ivref)
    inf_model_index = get(inf_model.point_lookup, (ivref, support), nothing)
    if !isnothing(inf_model_index)
        return InfiniteOpt._make_variable_ref(inf_model, inf_model_index)
    else
        # make negative index to not conflict with the InfiniteModel
        raw_index = data.last_point_index -= 1
        # make the reference and map it to a transcription variable
        pvref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(ivref), raw_index,
                                               InfiniteOpt.PointVariableIndex)
        trans_var = lookup_by_support(model, ivref, support)
        data.finvar_mappings[pvref] = trans_var
        data.point_lookup[(ivref, support)] = pvref
        return pvref
    end
end

"""
    InfiniteOpt.add_semi_infinite_variable(model::JuMP.Model,
                                     var::InfiniteOpt.SemiInfiniteVariable,
                                     key::Val{:TransData}
                                     )::InfiniteOpt.GeneralVariableRef

Make a `SemiInfiniteVariableRef` and add `var` to the transcription data 
and return the `GeneralVariableRef`. This is an extension of 
[`add_semi_infinite_variable`](@ref InfiniteOpt.add_semi_infinite_variable(::JuMP.Model,::Any,::Any)) 
for `TranscriptionOpt`. Note that `internal_semi_infinite_variable` is also 
extended to be able to access the `var`.
"""
function InfiniteOpt.add_semi_infinite_variable(
    model::JuMP.Model,
    var::InfiniteOpt.SemiInfiniteVariable,
    ::Val{:TransData}
    )::InfiniteOpt.GeneralVariableRef
    # check if an internal variable was already created
    ivref = var.infinite_variable_ref
    eval_supps = var.eval_supports
    data = transcription_data(model)
    internal_vref = get(data.semi_lookup, (ivref, eval_supps), nothing)
    if !isnothing(internal_vref)
        return internal_vref
    end
    # check if whether the infinite model already has one, else create one
    inf_model = JuMP.owner_model(ivref)
    inf_model_index = get(inf_model.semi_lookup, (ivref, eval_supps), nothing)
    if !isnothing(inf_model_index)
        return InfiniteOpt._make_variable_ref(inf_model, inf_model_index)
    else
        # make negative index to not conflict with the InfiniteModel
        semi_infinite_vars = transcription_data(model).semi_infinite_vars
        raw_index = -1 * (length(semi_infinite_vars) + 1)
        # make the reference and map it to a transcription variable
        rvref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(ivref), raw_index,
                                               InfiniteOpt.SemiInfiniteVariableIndex)
        push!(semi_infinite_vars, var)
        _set_semi_infinite_variable_mapping(model, var, rvref, InfiniteOpt._index_type(ivref))
        data.semi_lookup[(ivref, eval_supps)] = rvref
        return rvref
    end
end
