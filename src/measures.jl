struct _DumbyModel <: JuMP.AbstractModel end

################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(model::InfiniteModel,
                               index::MeasureIndex)::MeasureRef
    return MeasureRef(model, index)
end

# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::MeasureData)::MeasureIndex
    return MOIUC.add_item(model.measures, object)
end

# Extend _data_dictionary (reference based)
function _data_dictionary(mref::MeasureRef
    )::MOIUC.CleverDict{MeasureIndex, MeasureData}
    return JuMP.owner_model(mref).measures
end

# Extend _data_object
function _data_object(mref::MeasureRef)::MeasureData
    return _data_dictionary(mref)[JuMP.index(mref)]
end

# Extend _core_variable_object
function _core_variable_object(mref::MeasureRef)::Measure
    return _data_object(mref).measure
end

# Extend _object_numbers
function _object_numbers(mref::MeasureRef)::Vector{Int}
    return _core_variable_object(mref).object_nums
end

# Extend _parameter_numbers
function _parameter_numbers(mref::MeasureRef)::Vector{Int}
    return _core_variable_object(mref).parameter_nums
end

# Extend _set_core_variable_object
function _set_core_variable_object(mref::MeasureRef, object::Measure)::Nothing
    _data_object(mref).measure = object
    return
end

# Extend _measure_dependencies
function _measure_dependencies(mref::MeasureRef)::Vector{MeasureIndex}
    return _data_object(mref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(mref::MeasureRef)::Vector{ConstraintIndex}
    return _data_object(mref).constraint_indices
end

################################################################################
#                              MEASURE DATA METHODS
################################################################################
# Extend Base.copy for DiscreteMeasureData
function Base.copy(data::DiscreteMeasureData)::DiscreteMeasureData
    return DiscreteMeasureData(copy(data.parameter_refs), copy(data.coefficients),
                               copy(data.supports), data.label,
                               data.weight_function)
end

# Extend Base.copy for FunctionalDiscreteMeasureData
function Base.copy(data::FunctionalDiscreteMeasureData)::FunctionalDiscreteMeasureData
    return FunctionalDiscreteMeasureData(copy(data.parameter_refs),
                                         data.coeff_function,
                                         data.num_supports, data.label,
                                         data.weight_function)
end

"""
    default_weight(t) = 1

Default weight function for [`DiscreteMeasureData`](@ref) and
[`FunctionalDiscreteMeasureData`](@ref). Returns 1 regardless of the input value.
"""
default_weight(t) = 1

## NOTE THE BELOW CONSTRUCTOR IS INTENDED FOR USERS ONLY, NOT INTERNAL USE
"""
    DiscreteMeasureData(pref::GeneralVariableRef,
                        coefficients::Vector{<:Real},
                        supports::Vector{<:Real}; label::Symbol = gensym(),
                        weight_function::Function = [`default_weight`](@ref)
                        )::DiscreteMeasureData

Returns a 1-dimensional `DiscreteMeasureData` object that can be utilized
to define measures using [`measure`](@ref). This accepts input for a scalar (single)
infinite parameter. A description of the other arguments is provided in the
documentation for [`DiscreteMeasureData`](@ref). Errors if supports are out
bounds or an unequal number of supports and coefficients are given. Note that by
default a unique `label` is generated via `gensym` to ensure the supports can
be located in the infinite parameter support storage. Advanced implementations,
may choose a different behavior but should do so with caution.

**Example**
```julia-repl
julia> data = DiscreteMeasureData(pref, [0.5, 0.5], [1, 2])
DiscreteMeasureData{GeneralVariableRef,1}(pref, [0.5, 0.5], [1.0, 2.0], Symbol("##373"), default_weight)
```
"""
function DiscreteMeasureData(pref::GeneralVariableRef,
                             coefficients::Vector{<:Real},
                             supports::Vector{<:Real};
                             label::Symbol = gensym(), # makes a unique symbol
                             weight_function::Function = default_weight,
                             lower_bound::Float64 = NaN,
                             upper_bound::Float64 = NaN,
                             expect::Bool = false
                             )::DiscreteMeasureData{GeneralVariableRef, 1}
    if !(_index_type(pref) <: IndependentParameterIndex)
        error("$pref is not an independent parameter.")
    elseif length(coefficients) != length(supports)
        error("The amount of coefficients must match the amount of " *
              "support points.")
    end
    set = infinite_set(pref)
    if !supports_in_set(supports, set)
        error("Support points violate parameter domain.")
    end
    return DiscreteMeasureData(pref, coefficients, supports, label,
                               weight_function, lower_bound, upper_bound, expect)
end

# NOTE This will be useful for most user measure methods
# Ensure that an array of prefs is valid for use in multi-dimensional data
function _check_multidim_params(prefs::AbstractArray{GeneralVariableRef})::Nothing
    only_dep = all(_index_type(pref) == DependentParameterIndex for pref in prefs)
    only_indep = !only_dep && all(_index_type(pref) == IndependentParameterIndex
                                  for pref in prefs)
    if !only_dep && !only_indep
        error("Cannot specify a mixture of infinite parameter types.")
    elseif only_dep && any(_raw_index(pref) != _raw_index(first(prefs)) for pref in prefs)
        error("Cannot specify multiple dependent parameter groups for one measure.")
    elseif only_dep && length(prefs) != _num_parameters(first(dispatch_variable_ref.(prefs)))
        error("Cannot specify a subset of dependent parameters, consider using " *
              "nested one-dimensional measures instead.")
    end
    return
end

## Define function to intelligently format prefs and supports properly
# Vector
_convert_param_refs(prefs::Vector{GeneralVariableRef})::Vector{GeneralVariableRef} = prefs
function _convert_param_refs_and_supports(
    prefs::Vector{GeneralVariableRef},
    supports::Vector{<:Vector{T}}
    )::Tuple{Vector{GeneralVariableRef}, Matrix{T}} where {T <: Real}
    return prefs, reduce(hcat, supports)
end

# Array
function _convert_param_refs(
    prefs::Union{Array{GeneralVariableRef},
                 JuMPC.DenseAxisArray{GeneralVariableRef}}
    )::Vector{GeneralVariableRef}
    return reduce(vcat, prefs)
end

function _convert_param_refs_and_supports(
    prefs::Array{GeneralVariableRef},
    supports::Vector{<:Array{T}}
    )::Tuple{Vector{GeneralVariableRef}, Matrix{T}} where {T <: Real}
    return reduce(vcat, prefs), [supp[i] for i in eachindex(supports[1]), supp in supports]
end

# DenseAxisArray
function _convert_param_refs_and_supports(
    prefs::JuMPC.DenseAxisArray{GeneralVariableRef},
    supports::Vector{<:JuMPC.DenseAxisArray{T}}
    )::Tuple{Vector{GeneralVariableRef}, Matrix{T}} where {T <: Real}
    return reduce(vcat, prefs), [supp[i] for i in eachindex(first(supports)), supp in supports]
end

# SparseAxisArray
function _convert_param_refs(
    prefs::JuMPC.SparseAxisArray{GeneralVariableRef}
    )::Vector{GeneralVariableRef}
    indices = Collections._get_indices(prefs)
    return Collections._make_ordered(prefs, indices)
end

function _convert_param_refs_and_supports(
    prefs::JuMPC.SparseAxisArray{GeneralVariableRef},
    supports::Vector{<:JuMPC.SparseAxisArray{T}}
    )::Tuple{Vector{GeneralVariableRef}, Matrix{T}} where {T <: Real}
    indices = Collections._get_indices(prefs)
    ordered_prefs = Collections._make_ordered(prefs, indices)
    supps = [Collections._make_ordered(supp, indices) for supp in supports]
    return _convert_param_refs_and_supports(ordered_prefs, supps)
end

# NOTE THIS IS INTENDED ONLY FOR USER USE ONLY (IT IS NOT EFFICIENT FOR INTERNAL USE)
# NOTE INTERNAL METHODS SHOULD JUST GENERATE MULTI-DIM DATA DIRECTLY
"""
    DiscreteMeasureData(prefs::AbstractArray{GeneralVariableRef},
                        coefficients::Vector{<:Real},
                        supports::Vector{<:AbstractArray{<:Real}};
                        label::Symbol = gensym(),
                        weight_function::Function = [`default_weight`](@ref)
                        )::DiscreteMeasureData

Returns a `MultiDiscreteMeasureData` object that can be utilized to
define measures using [`measure`](@ref). This accepts input for an array (multi)
parameter. The inner arrays in the supports vector need to match the formatting
of the array used for `parameter_refs`. Note that `name` is used for printing
purposes and a description of the other arguments is provided in the
documentation for [`MultiDiscreteMeasureData`](@ref). Errors if supports are out
bounds, an unequal number of supports and coefficients are given, the array
formats do not match, or the parameters have different group IDs.

**Example**
```julia-repl
julia> data = DiscreteMeasureData(prefs, [0.5, 0.5], [[1, 1], [2, 2]]);

julia> typeof(data)
DiscreteMeasureDataDiscreteMeasureData{Vector{GeneralVariableRef},2}
```
"""
function DiscreteMeasureData(prefs::AbstractArray{GeneralVariableRef},
                             coefficients::Vector{<:Real},
                             supports::Vector{<:AbstractArray{<:Real}};
                             label::Symbol = gensym(),
                             weight_function::Function = default_weight,
                             lower_bound::Vector{Float64} = NaN * ones(length(prefs)),
                             upper_bound::Vector{Float64} = NaN * ones(length(prefs)),
                             expect::Bool = false
                             )::DiscreteMeasureData{Vector{GeneralVariableRef}, 2}
    _check_multidim_params(prefs) # ensures that prefs are valid
    if _keys(prefs) != _keys(first(supports))
        error("Parameter references and supports must use same container type.")
    end
    if length(coefficients) != size(supports, 2)
        error("The amount of coefficients must match the amount of " *
              "support points.")
    end
    vector_prefs, supps = _convert_param_refs_and_supports(prefs, supports)
    for i in eachindex(vector_prefs)
        set = infinite_set(vector_prefs[i])
        if !supports_in_set(supps[i, :], set)
            error("Support points violate parameter domain.")
        end
    end
    return DiscreteMeasureData(vector_prefs, coefficients, supps, label,
                               weight_function, lower_bound, upper_bound, expect)
end

# (done) TODO Make user constructor(s) for FunctionalDiscreteMeasureData
function FunctionalDiscreteMeasureData(
    pref::GeneralVariableRef,
    coeff_func::Function,
    num_supports::Int,
    label::Symbol;
    weight_function::Function = default_weight,
    lower_bound::Float64 = NaN,
    upper_bound::Float64 = NaN,
    expect::Bool = false
    )::FunctionalDiscreteMeasureData{GeneralVariableRef}
    # NOTE Might need to add some check for valid label input here or later...
    if !(_index_type(pref) <: IndependentParameterIndex)
        error("$pref is not an independent parameter.")
    end
    if num_supports <= 0
        error("Number of supports must be positive")
    end
    return FunctionalDiscreteMeasureData(pref, coeff_func, num_supports,
                                         label, weight_function,
                                         lower_bound, upper_bound, expect)
end

function FunctionalDiscreteMeasureData(
    prefs::AbstractArray{GeneralVariableRef},
    coeff_func::Function,
    num_supports::Int,
    label::Symbol;
    weight_function::Function = default_weight,
    lower_bound::Vector{Float64} = NaN * ones(length(prefs)),
    upper_bound::Vector{Float64} = NaN * ones(length(prefs)),
    expect::Bool = false
    )::FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}
    # NOTE Might need to add some check for valid label input here or later...
    _check_multidim_params(prefs)
    vector_prefs = _convert_param_refs(prefs)
    if num_supports <= 0
        error("Number of supports must be positive")
    end
    return FunctionalDiscreteMeasureData(vector_prefs, coeff_func, num_supports,
                                         label, weight_function, lower_bound,
                                         upper_bound, expect)
end

"""
    parameter_refs(data::AbstractMeasureData)::Union{GeneralVariableRef,
                                                     AbstractArray{GeneralVariableRef}}

Return the infinite parameter reference(s) in `data`. This is intended as an
internal function to be used with measure addition. User-defined measure data types
will need to extend this function otherwise an error is thrown.
"""
function parameter_refs(data::AbstractMeasureData)
    error("Function `parameter_refs` not extended for measure data of type $(typeof(data)).")
end

# DiscreteMeasureData
function parameter_refs(data::DiscreteMeasureData{T, N})::T where {T, N}
    return data.parameter_refs
end

# FunctionalDiscreteMeasureData
function parameter_refs(data::FunctionalDiscreteMeasureData{T})::T where {T}
    return data.parameter_refs
end

# NOTE Might need to change fallback behavior depending on how we implement measures
"""
    support_label(data::AbstractMeasureData)::Symbol

Return the label stored in `data` associated with its supports.
This is intended as en internal method for measure creation and ensures any
new supports are added to parameters with such a label.
User-defined measure data types should extend this function, otherwise an error
is thrown.
"""
function support_label(data::AbstractMeasureData)
    error("Function `support_label` not defined for measure data of type " *
          "$(typeof(data)).")
end

# DiscreteMeasureData and FunctionalDiscreteMeasureData
function support_label(data::Union{FunctionalDiscreteMeasureData,
               DiscreteMeasureData})::Symbol
    return data.label
end

# NOTE Might need to change fallback behavior depending on how we implement measures
"""
    supports(data::AbstractMeasureData)::Array{Float64}

Return the supports associated with `data` and its infinite parameters.
This is intended as en internal method for measure creation and ensures any
new supports are added to parameters. User-defined measure data types should
extend this function if appropriate, otherwise an empty vector is returned.
"""
function supports(data::AbstractMeasureData)::Vector{Float64}
    return Float64[]
end

# DiscreteMeasureData
function supports(data::DiscreteMeasureData{T, N})::Array{Float64, N} where {T, N}
    return data.supports
end

# FunctionalDiscreteMeasureData
function supports(data::FunctionalDiscreteMeasureData)::Array{Float64}
    return supports(parameter_refs(data), label = support_label(data))
end

# NOTE Might need to change fallback behavior depending on how we implement measures
"""
    num_supports(data::AbstractMeasureData)::Int

Return the number supports associated with `data` and its infinite parameters.
This is intended as en internal method for measure creation. User-defined
measure data types should extend this function if appropriate, otherwise
0 is returned.
"""
function num_supports(data::AbstractMeasureData)::Int
    return 0
end

# DiscreteMeasureData
function num_supports(data::DiscreteMeasureData)::Int
    return size(supports(data))[end]
end

# FunctionalDiscreteMeasureData
function num_supports(data::FunctionalDiscreteMeasureData)::Int
    return data.num_supports
end

"""
    coefficient_function(data::AbstractMeasureData)::Function

Return the coefficient function stored in `data` associated with its
expansion abstraction is there is such a function. This is intended as an
internal method for measure creation. User-defined measure
data types should extend this function if appropriate, otherwise an error is
thrown for unsupported types.
"""
function coefficient_function(data::AbstractMeasureData)::Function
    error("Function `coefficient_function` not defined for measure data of type " *
          "$(typeof(data)).")
end

# FunctionalDiscreteMeasureData
function coefficient_function(data::FunctionalDiscreteMeasureData)::Function
    return data.coeff_function
end

# NOTE Might need to change fallback behavior depending on how we implement measures
"""
    coefficients(data::AbstractMeasureData)::Vector{<:Real}

Return the coefficients associated with `data` associated with its expansion abstraction.
This is intended as en internal method for measure creation. User-defined measure
data types should extend this function if appropriate, otherwise an empty vector
is returned.
"""
function coefficients(data::AbstractMeasureData)::Vector{Float64}
    return Float64[]
end

# DiscreteMeasureData
function coefficients(data::DiscreteMeasureData)::Vector{Float64}
    return data.coefficients
end

# FunctionalDiscreteMeasureData (1-Dim)
function coefficients(data::FunctionalDiscreteMeasureData{GeneralVariableRef}
                      )::Vector{Float64}
    return [data.coeff_function(supp) for supp in supports(data)]
end

# FunctionalDiscreteMeasureData (n-Dim)
function coefficients(
    data::FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}
    )::Vector{Float64}
    supps = supports(data)
    coeffs = coefficient_function(data)
    return [coeffs(supps[:, i]) for i in 1:size(supps, 2)]
end

# NOTE Might need to change fallback behavior depending on how we implement measures
"""
    weight_function(data::AbstractMeasureData)::Function

Return the weight function stored in `data` associated with its expansion abstraction.
This is intended as en internal method for measure creation. User-defined measure
data types should extend this function if appropriate, otherwise an error is thrown.
"""
function weight_function(data::AbstractMeasureData)::Function
    error("Function `weight_function` not defined for measure data of type " *
          "$(typeof(data)).")
end

# DiscreteMeasureData and FunctionalDiscreteMeasureData
function weight_function(data::Union{DiscreteMeasureData,
                                     FunctionalDiscreteMeasureData})::Function
    return data.weight_function
end

################################################################################
#                            MEASURE CONSTRUCTION METHODS
################################################################################
"""
    measure_data_in_hold_bounds(data::AbstractMeasureData,
                                bounds::ParameterBounds)::Bool

Return a `Bool` whether the domain of `data` is valid in accordance with
`bounds`. This is intended as an internal method and is used to check hold
variables used in measures. User-defined measure data types will need to
extend this function to enable this error checking, otherwise it is skipped and
a warning is given.
"""
function measure_data_in_hold_bounds(data::AbstractMeasureData,
                                     bounds::ParameterBounds)::Bool
    @warn "Unable to check if hold variables bounds are valid in measure " *
           "with measure data type `$(typeof(data))`. This can be resolved by " *
           "extending `measure_data_in_hold_bounds`."
    return true
end

# Scalar DiscreteMeasureData and FunctionalDiscreteMeasureData
function measure_data_in_hold_bounds(
    data::Union{DiscreteMeasureData{P, 1}, FunctionalDiscreteMeasureData{P}},
    bounds::ParameterBounds{GeneralVariableRef}
    )::Bool where {P <: GeneralVariableRef}
    pref = parameter_refs(data)
    supps = supports(data)
    if haskey(bounds, pref) && length(supps) != 0
        return supports_in_set(supps, bounds[pref])
    end
    return true
end

# Multi-dimensional DiscreteMeasureData and FunctionalDiscreteMeasureData
function measure_data_in_hold_bounds(
    data::Union{DiscreteMeasureData{P, 2}, FunctionalDiscreteMeasureData{P}},
    bounds::ParameterBounds{GeneralVariableRef}
    )::Bool where {P <: Vector{GeneralVariableRef}}
    prefs = parameter_refs(data)
    supps = supports(data)
    if length(supps) != 0
        for i in eachindex(prefs)
            if haskey(bounds, prefs[i])
                if !supports_in_set(supps[i, :], bounds[prefs[i]])
                    return false
                end
            end
        end
    end
    return true
end

# Check that variables don't violate the parameter bounds
function _check_var_bounds(vref::GeneralVariableRef, data::AbstractMeasureData)
    if _index_type(vref) == HoldVariableIndex
        bounds = parameter_bounds(vref)
        if !measure_data_in_hold_bounds(data, bounds)
            error("Measure bounds violate hold variable bounds.")
        end
    elseif _index_type(vref) == MeasureIndex
        vrefs = _all_function_variables(measure_function(vref))
        for vref in vrefs
            _check_var_bounds(vref, data)
        end
    end
    return
end

# (done) TODO finish build_measure
"""
    build_measure(expr::JuMP.AbstractJuMPScalar,
                  data::AbstractMeasureData)::Measure

Build and return a [`Measure`](@ref) given the expression to be measured `expr`
using measure data `data`. This principally serves as an internal method for
measure definition. Errors if the supports associated with `data` violate
an hold variable parameter bounds of hold variables that are included in the
measure.
"""
function build_measure(model::InfiniteModel, expr::T, data::D;
    )::Measure{T, D} where {T <: JuMP.AbstractJuMPScalar, D <: AbstractMeasureData}
    vrefs = _all_function_variables(expr)
    if model.has_hold_bounds
        for vref in vrefs
            _check_var_bounds(vref, data)
        end
    end
    expr_obj_nums = _object_numbers(expr)
    expr_param_nums = _parameter_numbers(expr)
    prefs = parameter_refs(data)
    data_obj_num = _object_number(first(prefs))
    data_param_nums = [_parameter_number(pref) for pref in prefs]
    obj_nums = setdiff(expr_obj_nums, data_obj_num)
    param_nums = setdiff(expr_param_nums, data_param_nums)
    # check if analytic method should be applied
    constant_func = false
    lb_unique = unique(data.lower_bound)
    ub_unique = unique(data.upper_bound)
    if length(param_nums) == length(expr_param_nums) &&
        ((isequal(lb_unique, [NaN]) && isequal(ub_unique, [NaN])) || expect)
        constant_func = true
    end
    return Measure(expr, data, obj_nums, param_nums, constant_func)
end

################################################################################
#                               DEFINITION METHODS
################################################################################
function _add_supports_to_multiple_parameters(
    prefs::Vector{<:DependentParameterRef},
    supps::Array{Float64, 2},
    label::Symbol
    )::Nothing
    add_supports(prefs, supps, label = label, check = false)
    return
end

function _add_supports_to_multiple_parameters(
    prefs::Vector{<:IndependentParameterRef},
    supps::Array{Float64, 2},
    label::Symbol
    )::Nothing
    for i in eachindex(prefs)
        add_supports(prefs[i], supps[i, :], label = label, check = false)
    end
    return
end
"""
    add_supports_to_parameters(data::AbstractMeasureData,
                               [constant_func::Bool = false])::Nothing

Add supports as appropriate with `data` to the underlying infinite parameters.
This is an internal method with by [`add_measure`](@ref) and should be defined
for user-defined measure data types. `constant_func` indicates if the measure
expression depends on `data` or not.
"""
function add_supports_to_parameters(data::AbstractMeasureData)
    # TODO should it be an error or a warning?
    error("`add_supports_to_parameters` not defined for measures with " *
          "measure data type $(typeof(data)).")
end

## Internal functions for adding measure data supports to the parameter supports
# scalar DiscreteMeasureData
function add_supports_to_parameters(
    data::DiscreteMeasureData{GeneralVariableRef, 1}
    )::Nothing
    pref = parameter_refs(data)
    supps = supports(data)
    label = support_label(data)
    add_supports(pref, supps, label = label, check = false)
    return
end

# multi-dimensional DiscreteMeasureData
function add_supports_to_parameters(
    data::DiscreteMeasureData{Vector{GeneralVariableRef}, 2}
    )::Nothing
    prefs = map(p->dispatch_variable_ref(p), parameter_refs(data))
    supps = supports(data)
    label = support_label(data)
    _add_supports_to_multiple_parameters(prefs, supps, label)
    return
end

# scalar FunationalDiscreteMeasureData
function add_supports_to_parameters(
    data::FunctionalDiscreteMeasureData{GeneralVariableRef}
    )::Nothing
    pref = dispatch_variable_ref(parameter_refs(data))
    num_supps = num_supports(data)
    label = support_label(data)
    curr_num_supps = num_supports(pref, label = label)
    if curr_num_supps < num_supps
        # Assuming pref could not be DependentParameterRef
        generate_and_add_supports!(pref, infinite_set(pref), label,
                                   num_supports = num_supps - curr_num_supps)
    end
    return
end

# multi-dimensional FunationalDiscreteMeasureData
function add_supports_to_parameters(
    data::FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}
    )::Nothing
    prefs = dispatch_variable_ref.(parameter_refs(data))
    num_supps = num_supports(data)
    label = support_label(data)
    curr_num_supps = num_supports(prefs, label = label)
    if curr_num_supps < num_supps
        generate_and_add_supports!(prefs, infinite_set(prefs), label,
                                   num_supports = num_supps - curr_num_supps)
    end
    return
end

"""
    add_measure(model::InfiniteModel, meas::Measure,
                name::String = "measure")::GeneralVariableRef

Add a measure to `model` and return the corresponding measure reference. This
operates in a manner similar to `JuMP.add_variable`. Note this intended
as an internal method.
"""
function add_measure(model::InfiniteModel, meas::Measure,
                     name::String = "measure")::GeneralVariableRef
    # get the expression variables and check validity
    vrefs = _all_function_variables(meas.func)
    for vref in vrefs
        JuMP.check_belongs_to_model(vref, model)
    end
    # get the measure data info and check validity
    data = meas.data
    prefs = parameter_refs(data)
    for pref in prefs
        JuMP.check_belongs_to_model(pref, model)
    end
    # add supports to the model as needed
    if !meas.constant_func
        add_supports_to_parameters(data)
    end
    # add the measure to the model
    object = MeasureData(meas, name)
    mindex = _add_data_object(model, object)
    mref = _make_variable_ref(model, mindex)
    # update mappings
    for vref in union!(vrefs, prefs)
        push!(_measure_dependencies(vref), mindex)
    end
    return mref
end

"""
    measure_function(mref::MeasureRef)::JuMP.AbstractJuMPScalar

Return the function associated with `mref`.
"""
function measure_function(mref::MeasureRef)::JuMP.AbstractJuMPScalar
    return _core_variable_object(mref).func
end

"""
    measure_data(mref::MeasureRef)::AbstractMeasureData

Return the measure data associated with `mref`.
"""
function measure_data(mref::MeasureRef)::AbstractMeasureData
    return _core_variable_object(mref).data
end

"""
    is_analytic(mref::MeasureRef)::AbstractMeasureData

Return if `mref` is evaluated analytically.
"""
function is_analytic(mref::MeasureRef)::Bool
    return _core_variable_object(mref).constant_func
end

# NOTE This will be the preferred user method for small expressions
"""
    measure(expr::JuMP.AbstractJuMPScalar,
            data::AbstractMeasureData;
            name::String = "measure")::GeneralVariableRef

Return a measure reference that evaluates `expr` using according to `data`.
The measure data `data` determines how the measure is to be evaluated.
Typically, the [`DiscreteMeasureData`](@ref)
constructor can be used to for `data`. The variable expression `expr` can contain
`InfiniteOpt` variables, infinite parameters, other measure references (meaning
measures can be nested), and constants. Errors if `expr` does not contain
infinite variables, infinite parameters, or measure references. Also errors if
the measure parameter specified in `data` is not in `expr` and is not in any
the nested measure references. Typically, this is called inside of
[`JuMP.@expression`](@ref), [`JuMP.@objective`](@ref), and
[`JuMP.@constraint`](@ref) in a manner similar to `sum`. Note measures are not
explicitly evaluated until [`build_optimizer_model!`](@ref) is called or unless
they are expanded via [`expand`](@ref) or [`expand_all_measures!`](@ref).

**Example**
```julia-repl
julia> tdata = DiscreteMeasureData(t, [0.5, 0.5], [1, 2]);

julia> xdata = DiscreteMeasureData(xs, [0.5, 0.5], [[-1, -1], [1, 1]]);

julia> constr_RHS = @expression(model, measure(g - s + 2, tdata) + s^2)
measure(g(t) - s + 2) + s²

julia> @objective(model, Min, measure(g - 1  + measure(T, xdata), tdata))
measure(g(t) - 1 + name2(T(t, x)))
```
"""
function measure(expr::JuMP.AbstractJuMPScalar,
                 data::AbstractMeasureData;
                 name::String = "measure")::GeneralVariableRef
    model = _model_from_expr(expr)
    if model === nothing
        error("Expression contains no variables or parameters.")
    end
    # TODO make meaningful error function
    meas = build_measure(model, expr, data)
    return add_measure(model, meas, name)
end

# NOTE This will be the preferred user method for medium to large expressions
# (done) TODO Make @measure that will call measure but will build expr via @expression
#      Maybe this will also pass an error function to measure ofr build_measure
macro measure(expr, args...)
    _error(str...) = JuMP._macro_error(:measure, (expr, args...), str...)
    extra, kw_args, requestedcontainer = JuMPC._extract_kw_args(args)
    if length(extra) != 1
        _error("Incorrect number of arguments. Must be of form " *
               "@measure(expr, data).")
    elseif length(kw_args) > 0
        _error("No keyword arguments are accepted. Must be of form " *
               "@measure(expr, data).")
    end
    data = first(extra)
    expression = :( JuMP.@expression(InfiniteOpt._DumbyModel(), $expr) )
    mref = :( measure($expression, $data) )
    return esc(mref)
end

################################################################################
#                               NAMING METHODS
################################################################################
"""
    JuMP.name(mref::MeasureRef)::String

Extend `JuMP.name` to return the name associated with a measure
reference.
"""
function JuMP.name(mref::MeasureRef)::String
    return _data_object(mref).name
end

"""
    JuMP.set_name(mref::MeasureRef, name::String)::Nothing

Extend `JuMP.set_name` to specify the name of a measure reference.
"""
function JuMP.set_name(mref::MeasureRef, name::String)::Nothing
    _data_object(mref).name = name
    return
end

################################################################################
#                              DEPENDENCY METHODS
################################################################################
"""
    used_by_measure(mref::MeasureRef)::Bool

Return a `Bool` indicating if `mref` is used by a measure.

**Example**
```julia-repl
julia> used_by_measure(mref)
true
```
"""
function used_by_measure(mref::MeasureRef)::Bool
    return !isempty(_measure_dependencies(mref))
end

"""
    used_by_constraint(mref::MeasureRef)::Bool

Return a `Bool` indicating if `mref` is used by a constraint.

**Example**
```julia-repl
julia> used_by_constraint(mref)
false
```
"""
function used_by_constraint(mref::MeasureRef)::Bool
    return !isempty(_constraint_dependencies(mref))
end

"""
    used_by_objective(vmref::MeasureRef)::Bool

Return a `Bool` indicating if `mref` is used by the objective.

**Example**
```julia-repl
julia> used_by_objective(mref)
true
```
"""
function used_by_objective(mref::MeasureRef)::Bool
    return _data_object(mref).in_objective
end

"""
    is_used(mref::MeasureRef)::Bool

Return a `Bool` indicating if `mref` is used in the model.

**Example**
```julia-repl
julia> is_used(mref)
true
```
"""
function is_used(mref::MeasureRef)::Bool
    return used_by_measure(mref) || used_by_constraint(mref) ||
           used_by_objective(mref)
end

################################################################################
#                           HIGH LEVEL USER DEFINITION
################################################################################
# Define symbol inputs for general measure method types
const sampling = :sampling
const quadrature = :quadrature
const IntegralDefaults = Dict(:eval_method => sampling,
                              :num_supports => DefaultNumSupports,
                              :weight_func => default_weight,
                              :name => "integral",
                              :use_existing_supports => false)

"""
    integral_defaults()

Get the default keyword argument values for defining integrals in `model`.

```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> integral_defaults(model)
Dict{Symbol,Any} with 6 entries:
  :num_supports          => 10
  :eval_method           => nothing
  :name                  => "integral"
  :weight_func           => default_weight
  :use_existing_supports => false
```
"""
integral_defaults() = IntegralDefaults

"""
    set_integral_defaults(; kwargs...)

Set the default keyword argument settings for integrals of the specified model.
The keyword arguments of this function will be recorded in the default keyword
argument values of the model. If the keyword argument has been defined in the
model default, it will be overwritten with the new keyword argument value.
Otherwise, the default will record the new keyword argument and its value for
measures. The default values will be used by integrals constructed from
[`integral`](@ref) calls.

**Example**
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> integral_defaults(model)
Dict{Symbol,Any} with 6 entries:
  :num_supports          => 10
  :eval_method           => nothing
  :name                  => "integral"
  :weight_func           => default_weight
  :use_existing_supports => false

julia> set_integral_default(m, num_supports = 5, eval_method = quadrature, new_kwarg = true)

julia> integral_defaults(m)
Dict{Symbol,Any} with 6 entries:
  :new_kwarg             => true
  :num_supports          => 5
  :eval_method           => :quadrature
  :name                  => "integral"
  :weight_func           => default_weight
  :use_existing_supports => false
```
"""
function set_integral_defaults(; kwargs...)
    merge!(IntegralDefaults, kwargs)
    return
end

# TODO update this method to work
# NOTE this probably needs to be significantly redone and should not convert things
#      SparseAxisArrays as done previously, instead parameters should be expressed
#      a vector using the methods shown in the DiscreteMeasureData constructor
#      and supports should be stored as a vector or matrix since that is what
#      the measure data type use.
# NOTE Multi-dimensional params can be checked with _check_multidim_params
"""
    integral(expr::JuMP.AbstractJuMPScalar,
             [params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}},
             lb::Union{Number, AbstractArray{<:Number}, Nothing} = nothing,
             ub::Union{Number, AbstractArray{<:Number}, Nothing} = nothing];
             [eval_method::Symbol = sampling,
             num_supports::Int = 10,
             weight_func::Function = [`default_weight`](@ref),
             name = "integral",
             use_existing_supports::Bool = false,
             kwargs...])::GeneralVariableRef

Returns a measure reference that evaluates the integral of `expr` with respect
to infinite parameter(s) from `lb` to `ub`. This thus considers integrals of the
form: ``\\int_{p \\in P} expr(p) w(p) dp`` where ``p`` is an infinite parameter
(scalar or vector) and ``w`` is the weight function is 1 by default. This
function provides a high-level interface that ultimately constructs a
[`DiscreteMeasureData`](@ref) via `eval_method` that is used to call
[`measure`](@ref).

The arugments are as follows:
- `params`: the integral infinite parameter(s) in ``dp``
- `lb` & `ub`: integral upper and lower bounds (defaults to entire parameter domain)
Note that `params` is required if `expr` contains multiple parameters.

The keyword arguments are as follows:
- `eval_method`: method tha generates the supports and coefficients
    - `sampling` --> dispatch to appropriate sampling method (default)
    - `quadrature` --> dispatch to appropriate quadrature method
    - `mc_sampling` --> Monte Carlo sampling
    - `gauss_hermite` --> Gaussian Hermite quadrature (infinite interval)
    - `gauss_legendre` --> Gaussian Legendre quadrature (finite interval)
    - `gauss_laguerre` --> Gaussian Laguerre quadrature (semi-infinite interval)
    - `trapezoid` --> trapezoidal quadrature (finite interval)
- `num_supports`: The number of supports to be generated
- `weight_func`: ``w(p)`` above with parameter value inputs and scalar output
- `name`: the name used in printing
- `use_existing_supports`: Use all supports currently stored in `params`
Note that `mc_sampling` samples from the parameter interval for [`IntervalSet`](@ref)s,
and from the underlying distribution for [`DistributionSet`](@ref)s. Also,
`use_existing_supports` is useful for subsequent integral calls when using
`mc_sampling` such that all the measures use the same supports.

See [`set_integral_defaults`](@ref) to update the default keyword argument values
for all integral calls.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, x in [0., 1.])
x

julia> @infinite_variable(model, f(x))
f(x)

julia> int = integral(f, num_supports = 5)
(f(x))

julia> expand(int)
0.2 f(0.8236475079774124) + 0.2 f(0.9103565379264364) + 0.2 f(0.16456579813368521) + 0.2 f(0.17732884646626457) + 0.2 f(0.278880109331201)
```
"""
function integral(expr::JuMP.AbstractJuMPScalar,
                  params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}},
                  lb::Union{Number, AbstractArray{<:Number}, Nothing} = nothing,
                  ub::Union{Number, AbstractArray{<:Number}, Nothing} = nothing;
                  kwargs...)::GeneralVariableRef
    # NOTE: params is required
    # NOTE: no more SparseAxisArray

    # count number of parameters check array formatting
    num_params = length(params)
    if num_params == 0
        # error if empty GeneralVariableRef array is provided
        error("No infinite parameter is provided.")
    elseif num_params == 1
        params = first(params)
    else
        # check that all parameters are from the same group
        _allequal(group_id.(params)) || error("Multiple groups of parameters " *
                                              "are specified.")
        # Use SparseAxisArray for params, lb, ub if multiple parameters
        # TODO change this so it doesn't have to needlessly convert to a SparseAxisArray
        # and do so via typed functions
        params = convert(JuMPC.SparseAxisArray, params)
        if isa(lb, AbstractArray)
            lb = convert(JuMPC.SparseAxisArray, lb)
        end
        if isa(ub, AbstractArray)
            ub = convert(JuMPC.SparseAxisArray, ub)
        end
    end

    # collect model of the measure
    model = JuMP.owner_model(first(params))

    # collect keyword arguments
    kwargs = merge(model.integral_defaults, kwargs)
    eval_method = kwargs[:eval_method]
    num_supports = kwargs[:num_supports]
    name = kwargs[:name]
    weight_func = kwargs[:weight_func]
    use_existing_supports = kwargs[:use_existing_supports]

    # parse the sampling key if given
    if eval_method == sampling
        eval_method = mc_sampling
        kwargs[:eval_method] = mc_sampling
    end

    # delete unneeded keyword arguments
    delete!(kwargs, :use_existing_supports)
    delete!(kwargs, :num_supports)

    # make bounds arrays if needed
    if num_params > 1
        if isa(lb, Number)
            lb = JuMPC.SparseAxisArray(Dict(k => lb for k in keys(params)))
        end
        if isa(ub, Number)
            ub = JuMPC.SparseAxisArray(Dict(k => ub for k in keys(params)))
        end
    end

    # check lower bounds TODO replace with typed functions
    bounds_valid = true
    if isa(lb, Number)
        bounds_valid = supports_in_set(lb, _parameter_set(first(params)))
    elseif isa(lb, JuMPC.SparseAxisArray)
        bounds_valid = all(supports_in_set.(lb, _parameter_set.(params)))
    end
    bounds_valid || error("Lower bound(s) violate(s) the infinite set domain.")

    # check upper bounds TODO replace with typed functions
    if isa(ub, Number)
        bounds_valid = supports_in_set(ub, _parameter_set(first(params)))
    elseif isa(ub, JuMPC.SparseAxisArray)
        bounds_valid = all(supports_in_set.(ub, _parameter_set.(params)))
    end
    bounds_valid || error("Upper bound(s) violate(s) the infinite set domain.")

    # fill in bounds if needed
    set = _parameter_set(first(params))
    if JuMP.has_lower_bound(set)
        # Fill in lower bounds and upper bounds if not given
        if isa(lb, Nothing) || length(lb) == 0
            lb = JuMP.lower_bound.(params)
        end
        if isa(ub, Nothing) || length(ub) == 0
            ub = JuMP.upper_bound.(params)
        end

        # Check the dimension of lb and ub matches number of parameters
        if length(lb) != num_params || length(ub) != num_params
            error("Number of parameters do not match number of lower bounds or " *
                  "upper bounds.")
        end

        if num_params == 1
            if isa(lb, AbstractArray)
                lb = lb[findfirst(x->isa(x, Number), lb)]
            end
            if isa(ub, AbstractArray)
                ub = ub[findfirst(x->isa(x, Number), ub)]
            end
        end

        # Check the input lower bounds and upper bounds are reasonable
        for i in eachindex(lb)
            if lb[i] >= ub[i]
                error("Lower bound is not less than upper bound for some " *
                      "parameter. Please check the input lower bounds " *
                      "and upper bounds.")
            end
        end
    end

    # construct data and return measure if we use existing supports
    if use_existing_supports
        if eval_method == quadrature || eval_method == gauss_legendre ||
           eval_method == gauss_hermite || eval_method == gauss_laguerre
            @warn("Quadrature method will not be used because " *
                  "use_existing_supports is set as true.")
        end
        support = supports(params)
        if !isa(lb, Nothing) && !isa(ub, Nothing)
            support = [i for i in support if all(i .>= lb) && all(i .<= ub)]
        end

        # TODO: generate reasonable coefficients for given supports
        # TODO replace with typed functions
        len = length(support)
        if isa(set, DistributionSet)
            data = DiscreteMeasureData(params, ones(len) ./ len, support,
                                       name = name, weight_function = weight_func)
        elseif isa(set, IntervalSet) && eval_method == mc_sampling
            coeffs = ones(len) / len * prod(ub .- lb)
            data = DiscreteMeasureData(params, coeffs, support, name = name,
                                       weight_function = weight_func)
        else
            data = DiscreteMeasureData(params, ones(len), support, name = name,
                                       weight_function = weight_func)
        end
        return measure(expr, data)
    end

    if eval_method == quadrature || eval_method == gauss_legendre ||
       eval_method == gauss_hermite || eval_method == gauss_laguerre
        if num_params > 1
            error("Quadrature method is not supported for multivariate measures.")
        end
        inf_bound_num = (lb == -Inf) + (ub == Inf)
        if inf_bound_num == 0
            kwargs[:eval_method] = gauss_legendre
        elseif inf_bound_num == 1
            kwargs[:eval_method] = gauss_laguerre
        else
            kwargs[:eval_method] = gauss_hermite
        end
    end

    # construct DiscreteMeasureData as data
    data = generate_measure_data(params, num_supports, lb, ub; kwargs...)

    # call measure function to construct the measure
    return measure(expr, data)
end

# TODO update this method
"""
    expect(expr::JuMP.AbstractJuMPScalar,
           [params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef},
                          Nothing} = nothing];
           [num_supports::Int = 10,
            use_existing_supports::Bool = false])::GeneralVariableRef

Creates a measure that represents the expected value of an expression in
a random parameter involved in the expression. Return the [`MeasureRef`](@ref)
of the created measure.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP, Distributions, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, x in Normal(0., 1.))
x

julia> @infinite_variable(model, f(x))
f(x)

julia> meas = expect(f, num_supports = 2)
expect(f(x))

julia> expand(meas)
0.5 f(0.6791074260357777) + 0.5 f(0.8284134829000359)
```
"""
function expect(expr::JuMP.AbstractJuMPScalar,
                params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef},
                              Nothing} = nothing;
                num_supports::Int = 10,
                use_existing_supports::Bool = false)::GeneralVariableRef
    # expectation measure
    return integral(expr, params, num_supports = num_supports,
                    name = "expect", use_existing_supports = use_existing_supports)
end

# TODO update this method
"""
    support_sum(expr::JuMP.AbstractJuMPScalar,
                [params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef},
                               Nothing} = nothing])::GeneralVariableRef

Creates a measure that represents the sum of the expression over a parameter
using its existing supports. Return the [`MeasureRef`](@ref) of the created
measure.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, x in [0, 1], supports = [0.3, 0.7])
x

julia> @infinite_variable(model, f(x))
f(x)

julia> meas = support_sum(f)
sum(f(x))

julia> expand(meas)
f(0.3) + f(0.7)
```
"""
function support_sum(expr::JuMP.AbstractJuMPScalar,
                     params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef},
                                   Nothing} = nothing)::GeneralVariableRef
    # sum measure
    return integral(expr, params, use_existing_supports = true, name = "sum",
                    eval_method = :not_default)
end

# TODO make macros for integral, expect, and support_sum like @measure

################################################################################
#                                   DELETION
################################################################################
# TODO update this to work
# TODO refer to variable and constraint deletion
"""
    JuMP.delete(model::InfiniteModel, mref::MeasureRef)::Nothing

Extend [`JuMP.delete`](@ref) to delete measures. Errors if measure is invalid,
meaning it does not belong to the model or it has already been deleted.

**Example**
```julia-repl
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z >= 0.0
 measure(g(t)) == 0
 g(t) + z >= 42.0
 g(0.5) == 0
 t in [0, 6]

julia> delete(model, meas)

julia> print(model)
Min z
Subject to
 z >= 0.0
 0 == 0
 g(t) + z >= 42.0
 g(0.5) == 0
 t in [0, 6]
```
"""
function JuMP.delete(model::InfiniteModel, mref::MeasureRef)::Nothing
    @assert JuMP.is_valid(model, mref) "Invalid measure reference."
    # Reset the transcription status
    if is_used(mref)
        set_optimizer_model_ready(model, false)
    end
    # Remove from dependent measures if there are any
    if used_by_measure(mref)
        for mindex in model.meas_to_meas[JuMP.index(mref)]
            if isa(model.measures[mindex].func, MeasureRef)
                data = model.measures[mindex].data
                model.measures[mindex] = Measure(zero(JuMP.AffExpr), data)
            else
                _remove_variable(model.measures[mindex].func, mref)
            end
            JuMP.set_name(MeasureRef(model, mindex),
                          _make_meas_name(model.measures[mindex]))
        end
        delete!(model.meas_to_meas, JuMP.index(mref))
    end
    # Remove from dependent constraints if there are any
    if used_by_constraint(mref)
        for cindex in model.meas_to_constrs[JuMP.index(mref)]
            if isa(model.constrs[cindex].func, MeasureRef)
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr),
                                                      model.constrs[cindex].set)
            else
                _remove_variable(model.constrs[cindex].func, mref)
            end
        end
        delete!(model.meas_to_constrs, JuMP.index(mref))
    end
    # Remove from objective if used there
    if used_by_objective(mref)
        if isa(model.objective_function, MeasureRef)
            model.objective_function = zero(JuMP.AffExpr)
        else
            _remove_variable(model.objective_function, mref)
        end
    end
    # Update that the parameters from the data are no longer dependent
    prefs = parameter_refs(measure_data(mref))
    for pref in prefs
        filter!(e -> e != JuMP.index(mref), _measure_dependencies(pref))
    end
    # Update that the variables used by it are no longer used by it
    vrefs = _all_function_variables(measure_function(mref))
    for vref in vrefs
        filter!(e -> e != JuMP.index(mref), _measure_dependencies(vref))
    end
    # delete remaining measure information
    _delete_data_object(mref)
    return
end
