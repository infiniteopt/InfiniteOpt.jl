"""
    DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                        Vector{<:JuMP.AbstractVariableRef}},
                        N, B <: Union{Float64, Vector{Float64}},
                        F <: Function
                        } <: AbstractMeasureData

A DataType for immutable measure abstraction data where the
abstraction is of the form:
``measure = \\int_{\\tau \\in T} f(\\tau) w(\\tau) d\\tau \\approx \\sum_{i = 1}^N \\alpha_i f(\\tau_i) w(\\tau_i)``.
The supports and coefficients are immutable (i.e., they will not change
even if supports are changed for the underlying infinite parameter.) This
type can be used for both 1-dimensional and multi-dimensional measures.

**Fields**
- `parameter_refs::P`: The infinite parameter(s) over which the integration occurs.
                       These can be comprised of multiple independent parameters,
                       but dependent parameters cannot be mixed with other types.
- `coefficients::Vector{Float64}`: Coefficients ``\\alpha_i`` for the above
                                   measure abstraction.
- `supports::Array{Float64, N}`: Supports points ``\\tau_i``. This is a `Vector`
                                 if only one parameter is given, otherwise it is
                                 a `Matrix` where the supports are stored column-wise.
- `label::DataType`: Label for the support points ``\\tau_i`` when stored in the
                   infinite parameter(s), stemming from [`AbstractSupportLabel`](@ref).
- `weight_function::F`: Weighting function ``w`` must map an individual
                               support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bound in accordance with ``T``, this denotes the
                    intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bound.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct DiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                           Vector{<:JuMP.AbstractVariableRef}},
                           N, B <: Union{Float64, Vector{Float64}},
                           F <: Function
                           } <: AbstractMeasureData
    parameter_refs::P
    coefficients::Vector{Float64}
    supports::Array{Float64, N} # supports are stored column-wise
    label::DataType # label that will used when the supports are added to the model
    weight_function::F # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function DiscreteMeasureData(
        param_ref::V, coeffs::Vector{<:Real},
        supps::Vector{<:Real}, 
        label::DataType,
        weight_func::F,
        lower_bound::Real,
        upper_bound::Real,
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F <: Function}
        return new{V, 1, Float64, F}(param_ref, coeffs, supps, label, weight_func,
                                     lower_bound, upper_bound, expect)
    end
    # multi constructor
    function DiscreteMeasureData(
        param_refs::Vector{V}, 
        coeffs::Vector{<:Real},
        supps::Matrix{<:Real}, 
        label::DataType,
        weight_func::F,
        lower_bound::Vector{<:Real},
        upper_bound::Vector{<:Real},
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F <: Function}
        return new{Vector{V}, 2, Vector{Float64}, F}(param_refs, coeffs, supps,
                                                     label, weight_func, lower_bound,
                                                     upper_bound, expect)
    end
end

"""
    FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                  Vector{<:JuMP.AbstractVariableRef}},
                                  B <: Union{Float64, Vector{Float64}},
                                  I <: AbstractGenerativeInfo,
                                  F1 <: Function,
                                  F2 <: Function
                                  } <: AbstractMeasureData

A DataType for mutable measure abstraction data where the
abstraction is of the form:
``measure = \\int_{\\tau \\in T} f(\\tau) w(\\tau) d\\tau \\approx \\sum_{i = 1}^N \\alpha_i f(\\tau_i) w(\\tau_i)``.
This abstraction is equivalent to that of [`DiscreteMeasureData`](@ref), but
the difference is that the supports are not fully known at the time of measure
creation. Thus, functions are stored that will be used to generate the
concrete support points ``\\tau_i`` and their coefficients ``\\alpha_i`` when
the measure is evaluated (expanded). These supports are identified/generated
in accordance with the `label` with a gaurantee that at least `num_supports` are
generated. For example, if `label = MCSample` and `num_supports = 100` then
the measure will use all of the supports stored in the `parameter_refs` with the
label `MCSample` and will ensure there are at least 100 are generated. This
type can be used for both 1-dimensional and multi-dimensional measures.

For 1-dimensional measures over independent infinite parameters, the 
`generative_supp_info` specifies the info needed to make generative supports based 
on those with that exist with `label`. Note that only 1 kind of generative 
supports are allowed for each infinite parameter.

**Fields**
- `parameter_refs::P`: The infinite parameter(s) over which the integration occurs.
                     These can be comprised of multiple independent parameters,
                     but dependent parameters cannot be mixed with other types.
- `coeff_function::F1`: Coefficient generation function making ``\\alpha_i``
                              for the above measure abstraction. It should take
                              all the supports as input (formatted as an Array)
                              and return the corresponding vector of coefficients.
- `min_num_supports::Int`: Specifies the minimum number of supports ``\\tau_i``
                       desired in association with `parameter_refs` and `label`.
- `label::DataType`: Label for the support points ``\\tau_i`` which are/will be
                   stored in the infinite parameter(s), stemming from [`AbstractSupportLabel`](@ref).
- `generative_supp_info::I`: Information needed to generate supports based on other 
   existing ones.
- `weight_function::F2`: Weighting function ``w`` must map an individual
                              support value to a `Real` scalar value.
- `lower_bounds::B`: Lower bounds in accordance with ``T``, this denotes the
                  intended interval of the measure and should be `NaN` if ignored
- `upper_bounds::B`: Same as above but the upper bounds.
- `is_expect::Bool`: Is this data associated with an expectation call?
"""
struct FunctionalDiscreteMeasureData{P <: Union{JuMP.AbstractVariableRef,
                                     Vector{<:JuMP.AbstractVariableRef}},
                                     B <: Union{Float64, Vector{Float64}},
                                     I <: AbstractGenerativeInfo,
                                     F1 <: Function,
                                     F2 <: Function
                                     } <: AbstractMeasureData
    parameter_refs::P
    coeff_function::F1 # supports (excluding generative)--> coefficient vector (includes generative)
    min_num_supports::Int # minimum number of supports
    label::DataType # support label of included supports
    generative_supp_info::I
    weight_function::F2 # single support --> weight value
    lower_bounds::B
    upper_bounds::B
    is_expect::Bool
    # scalar constructor
    function FunctionalDiscreteMeasureData(
        param_ref::V, 
        coeff_func::F1,
        num_supps::Int, 
        label::DataType,
        gen_info::I,
        weight_func::F2,
        lower_bound::Real,
        upper_bound::Real,
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, I <: AbstractGenerativeInfo,
                 F1 <: Function, F2 <: Function}
        return new{V, Float64, I, F1, F2}(param_ref, coeff_func, num_supps, label, 
                                          gen_info, weight_func, lower_bound, 
                                          upper_bound, expect)
    end
    # multi constructor
    function FunctionalDiscreteMeasureData(
        param_refs::Vector{V},
        coeff_func::F1,
        num_supps::Int, 
        label::DataType,
        weight_func::F2,
        lower_bound::Vector{<:Real},
        upper_bound::Vector{<:Real},
        expect::Bool
        ) where {V <: JuMP.AbstractVariableRef, F1 <: Function, F2 <: Function}
        return new{Vector{V}, Vector{Float64}, NoGenerativeSupports, F1, F2}(
            param_refs, coeff_func, num_supps, label, NoGenerativeSupports(), 
            weight_func, lower_bound, upper_bound, expect)
    end
end

# Convenient Dispatch constructor 
function FunctionalDiscreteMeasureData(
    param_refs::Vector{V},
    coeff_func::Function,
    num_supps::Int, 
    label::DataType,
    info::NoGenerativeSupports,
    weight_func::Function,
    lower_bound::Vector{<:Real},
    upper_bound::Vector{<:Real},
    expect::Bool
    ) where {V <: JuMP.AbstractVariableRef}
    return FunctionalDiscreteMeasureData(param_refs, coeff_func, num_supps,
                                         label, weight_func, lower_bound,
                                         upper_bound, expect)
end

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
