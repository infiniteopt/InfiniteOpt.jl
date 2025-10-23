"""
    InfiniteOpt.add_point_variable(
        backend::TranscriptionBackend,
        var::InfiniteOpt.PointVariable,
        support::Vector{Float64}
        )::InfiniteOpt.GeneralVariableRef

Make a `PointVariableRef` and map it to the appropriate transcription variable
and return the `GeneralVariableRef`. This is an extension of
[`add_point_variable`](@ref InfiniteOpt.add_point_variable(::InfiniteOpt.AbstractTransformationBackend,::Any,::Any))
for `TranscriptionOpt`.
"""
function InfiniteOpt.add_point_variable(
    backend::TranscriptionBackend,
    ivref::InfiniteOpt.GeneralVariableRef,
    support::Vector{Float64},
    )
    # check if an internal variable was already created
    data = transcription_data(backend)
    internal_vref = get(data.point_lookup, (ivref, support), nothing)
    if !isnothing(internal_vref)
        return internal_vref
    end
    # check if whether the infinite model already has one, else create one
    inf_model = JuMP.owner_model(ivref)
    inf_model_index = get(inf_model.point_lookup, (ivref, support), nothing)
    if !isnothing(inf_model_index)
        return InfiniteOpt.GeneralVariableRef(inf_model, inf_model_index)
    else
        # make negative index to not conflict with the InfiniteModel
        raw_index = data.last_point_index -= 1
        # make the reference and map it to a transcription variable
        pvref = InfiniteOpt.GeneralVariableRef(inf_model, raw_index, InfiniteOpt.PointVariableIndex)
        trans_var = lookup_by_support(ivref, backend, support)
        if trans_var isa Float64
            data.point_pfunc_mappings[pvref] = trans_var
        else
            data.finvar_mappings[pvref] = trans_var
        end
        data.point_lookup[(ivref, support)] = pvref
        return pvref
    end
end

"""
    InfiniteOpt.add_semi_infinite_variable(
        backend::TranscriptionBackend,
        var::InfiniteOpt.SemiInfiniteVariable
        )::InfiniteOpt.GeneralVariableRef

Make a `SemiInfiniteVariableRef` and add `var` to the transcription data 
and return the `GeneralVariableRef`. This is an extension of 
[`add_semi_infinite_variable`](@ref InfiniteOpt.add_semi_infinite_variable(::InfiniteOpt.AbstractTransformationBackend,::Any)) 
for `TranscriptionOpt`. Note that `internal_semi_infinite_variable` is also 
extended to be able to access the `var`.
"""
function InfiniteOpt.add_semi_infinite_variable(
    backend::TranscriptionBackend,
    var::InfiniteOpt.SemiInfiniteVariable
    )
    # check if an internal variable was already created
    ivref = var.infinite_variable_ref
    eval_supp = var.eval_support
    data = transcription_data(backend)
    internal_vref = get(data.semi_lookup, (ivref, eval_supp), nothing)
    if !isnothing(internal_vref)
        return internal_vref
    end
    # check if whether the infinite model already has one, else create one
    inf_model = JuMP.owner_model(ivref)
    inf_model_index = get(inf_model.semi_lookup, (ivref, eval_supp), nothing)
    if !isnothing(inf_model_index)
        return InfiniteOpt.GeneralVariableRef(inf_model, inf_model_index)
    else
        # make negative index to not conflict with the InfiniteModel
        semi_infinite_vars = data.semi_infinite_vars
        raw_index = -1 * (length(semi_infinite_vars) + 1)
        # make the reference and map it to a transcription variable
        rvref = InfiniteOpt.GeneralVariableRef(
            inf_model,
            raw_index,
            InfiniteOpt.SemiInfiniteVariableIndex
        )
        push!(semi_infinite_vars, var)
        ivref_param_nums = InfiniteOpt._parameter_numbers(ivref)
        param_nums = var.parameter_nums
        supp_indices = support_index_iterator(backend, var.group_int_idxs)
        if ivref.index_type == InfiniteOpt.ParameterFunctionIndex && 
           !data.update_parameter_functions
            val_type = Float64
        else
            val_type = JuMP.VariableRef
        end
        lookup_dict = Dict{Vector{Float64}, val_type}()
        sizehint!(lookup_dict, length(supp_indices))
        for i in supp_indices
            raw_supp = index_to_support(backend, i)
            ivref_supp = [isnan(s) ? raw_supp[ivref_param_nums[j]] : s 
                          for (j, s) in enumerate(eval_supp)]
            supp = raw_supp[param_nums]
            lookup_dict[supp] = lookup_by_support(ivref, backend, ivref_supp)
        end
        if val_type == JuMP.VariableRef
            data.infvar_lookup[rvref] = lookup_dict
        else
            data.pfunc_lookup[rvref] = lookup_dict
        end
        data.semi_lookup[(ivref, eval_supp)] = rvref
        return rvref
    end
end
