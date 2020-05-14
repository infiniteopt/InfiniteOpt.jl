"""
    InfiniteOpt.add_measure_variable(model::JuMP.Model,
                                     var::InfiniteOpt.PointVariable,
                                     key::Val{:TransData}
                                     )::InfiniteOpt.GeneralVariableRef

Make a `PointVariableRef` and map it to the appropriate transcription variable
and return the `PointVariableRef`. This is an extension of
[`add_measure_variable`](@ref InfiniteOpt.add_measure_variable(::JuMP.Model,::Any,::Any))
for `TranscriptionOpt`.
"""
function InfiniteOpt.add_measure_variable(model::JuMP.Model,
                                          var::InfiniteOpt.PointVariable,
                                          key::Val{:TransData}
                                          )::InfiniteOpt.GeneralVariableRef
    # make negative index to not conflict with the InfiniteModel
    raw_index = transcription_data(model).last_point_index -= 1
    # make the reference and map it to a transcription variable
    ivref = var.infinite_variable_ref
    pvref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(ivref), raw_index,
                                           InfiniteOpt.PointVariableIndex)
    trans_var = lookup_by_support(model, ivref, var.parameter_values)
    transcription_data(model).finvar_mappings[pvref] = trans_var
    return pvref
end

"""
    InfiniteOpt.add_measure_variable(model::JuMP.Model,
                                     var::InfiniteOpt.ReducedVariable,
                                     key::Val{:TransData}
                                     )::InfiniteOpt.GeneralVariableRef

Make a `ReducedVariableRef` and add `var` to the transcription data
and return the `ReducedVariableRef`. This is an extension of
[`add_measure_variable`](@ref InfiniteOpt.add_measure_variable(::JuMP.Model,::Any,::Any))
for `TranscriptionOpt`. Note that `internal_reduced_variable` is also extended
to be able to access the `var`.
"""
function InfiniteOpt.add_measure_variable(model::JuMP.Model,
                                          var::InfiniteOpt.ReducedVariable,
                                          key::Val{:TransData}
                                          )::InfiniteOpt.GeneralVariableRef
    # make negative index to not conflict with the InfiniteModel
    reduced_vars = transcription_data(model).reduced_vars
    raw_index = -1 * (length(reduced_vars) + 1)
    # make the reference and map it to a transcription variable
    rvref = InfiniteOpt.GeneralVariableRef(JuMP.owner_model(ivref), raw_index,
                                           InfiniteOpt.ReducedVariableIndex)
    push!(reduced_vars, var)
    # TODO map the reduced variable (will probably make function to do this)
    return rvref
end

"""
    InfiniteOpt.delete_reduced_variable(model::JuMP.Model,
                                        vref::InfiniteOpt.ReducedVariableRef,
                                        key::Val{:TransData})::Nothing

This is an extension of
[`delete_reduced_variable`](@ref InfiniteOpt.delete_reduced_variable(::JuMP.Model, ::Any, ::Any))
for use in `TranscriptionOpt`. Here we do not delete reduced variables once they
have been used since there is no performance gain for this paradigm and the
memory saving is small. Note this may change in the future.
"""
function InfiniteOpt.delete_reduced_variable(model::JuMP.Model,
                                             rvref::InfiniteOpt.ReducedVariableRef,
                                             key::Val{:TransData})::Nothing
    return
end
