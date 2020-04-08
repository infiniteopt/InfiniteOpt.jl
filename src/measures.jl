# Define symbol inputs for general measure method types
const sampling = :sampling
const quadrature = :quadrature

# Extend Base.copy for new variable types
Base.copy(v::MeasureRef, new_model::InfiniteModel) = MeasureRef(new_model, v.index)

# Extend Base.copy for DiscreteMeasureData
function Base.copy(data::DiscreteMeasureData)::DiscreteMeasureData
    return DiscreteMeasureData(copy(data.parameter_ref), copy(data.coefficients),
                               copy(data.supports), data.name,
                               data.weight_function)
end

# Extend Base.copy for MultiDiscreteMeasureData
function Base.copy(data::MultiDiscreteMeasureData)::MultiDiscreteMeasureData
    return MultiDiscreteMeasureData(copy(data.parameter_refs), copy(data.coefficients),
                                    copy(data.supports), data.name,
                                    data.weight_function)
end

"""
    JuMP.name(mref::MeasureRef)::String

Extend [`JuMP.name`](@ref) to return the name associated with a measure
reference.
"""
function JuMP.name(mref::MeasureRef)::String
    return JuMP.owner_model(mref).meas_to_name[mref.index]
end

"""
    JuMP.set_name(mref::MeasureRef, name::String)

Extend [`JuMP.set_name`](@ref) to specify the name of a measure reference.
"""
function JuMP.set_name(mref::MeasureRef, name::String)
    JuMP.owner_model(mref).meas_to_name[JuMP.index(mref)] = name
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, mref::MeasureRef)::Bool

Extend [`JuMP.is_valid`](@ref) to return `Bool` whether `mref` is valid.
"""
function JuMP.is_valid(model::InfiniteModel, mref::MeasureRef)::Bool
    return (model === JuMP.owner_model(mref) && JuMP.index(mref) in keys(model.measures))
end

"""
    default_weight(t)

Default weight function for [`DiscreteMeasureData`](@ref) and
[`MultiDiscreteMeasureData`](@ref). Returns 1 regardless of the input value.
"""
default_weight(t) = 1

"""
    DiscreteMeasureData(parameter_ref::ParameterRef,
                        coefficients::Vector{<:Real},
                        supports::Vector{<:Real}; name::String = "measure",
                        weight_function::Function = [`default_weight`](@ref))::DiscreteMeasureData

Returns a `DiscreteMeasureData` object that can be utilized to define
measures using [`measure`](@ref). This accepts input for a scalar (single)
parameter. Note that `name` is used for printing purposes and a description of
the other arguments is provided in the documentation for
[`DiscreteMeasureData`](@ref). Errors if supports are out bounds or an unequal
number of supports and coefficients are given.

**Example**
```julia-repl
julia> data = DiscreteMeasureData(pref, [0.5, 0.5], [1, 2], name = "example")
DiscreteMeasureData(pref, [0.5, 0.5], [1.0, 2.0], "example", default_weight)
```
"""
function DiscreteMeasureData(parameter_ref::ParameterRef,
                             coefficients::Vector{<:Real},
                             supports::Vector{<:Real};
                             name::String = "measure",
                             weight_function::Function = default_weight
                             )::DiscreteMeasureData
    is_finite_parameter(parameter_ref) && error("Measure parameter cannot be " *
                                                "finite.")
    if length(coefficients) != length(supports)
        error("The amount of coefficients must match the amount of " *
              "support points.")
    end
    set = infinite_set(parameter_ref)
    if !supports_in_set(supports, set)
        error("Support points violate parameter domain.")
    end
    return DiscreteMeasureData(parameter_ref, coefficients, supports, name,
                               weight_function)
end

## Define function to intelligently format prefs and supports properly
# Vector
function _convert_param_refs_and_supports(parameter_refs::Vector{ParameterRef},
                                          supports::Vector{<:Vector{<:Real}})::Tuple
    return parameter_refs, reduce(hcat, supports)
end

# Array
function _convert_param_refs_and_supports(parameter_refs::Array{ParameterRef},
                                          supports::Vector{<:Array{<:Real}})::Tuple
    return [parameter_refs...], [supp[i] for i in eachindex(supports[1]), supp in supports]
end

# DenseAxisArray
function _convert_param_refs_and_supports(parameter_refs::JuMPC.DenseAxisArray{ParameterRef},
                                          supports::Vector{<:JuMPC.DenseAxisArray{<:Real}})::Tuple
    return [parameter_refs...], [supp[i] for i in eachindex(first(supports)), supp in supports]
end

# SparseAxisArray
function _convert_param_refs_and_supports(parameter_refs::JuMPC.SparseAxisArray{ParameterRef},
                                          supports::Vector{<:JuMPC.SparseAxisArray{<:Real}})::Tuple
    indices = Collections._get_indices(parameter_refs)
    parameter_refs = Collections._make_ordered(parameter_refs, indices)
    supports = [Collections._make_ordered(supp, indices) for supp in supports]
    return _convert_param_refs_and_supports(parameter_refs, supports)
end

"""
    DiscreteMeasureData(parameter_refs::AbstractArray{<:ParameterRef},
                        coefficients::Vector{<:Real},
                        supports::Vector{<:AbstractArray{<:Real}};
                        name::String = "measure",
                        weight_function::Function = [`default_weight`](@ref)
                        )::MultiDiscreteMeasureData

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
julia> data = DiscreteMeasureData(prefs, [0.5, 0.5], [[1, 1], [2, 2]], name = "example");

julia> typeof(data)
MultiDiscreteMeasureData
```
"""
function DiscreteMeasureData(parameter_refs::AbstractArray{<:ParameterRef},
                             coefficients::Vector{<:Real},
                             supports::Vector{<:AbstractArray{<:Real}};
                             name::String = "measure",
                             weight_function::Function = default_weight
                             )::MultiDiscreteMeasureData
    any(is_finite_parameter.(parameter_refs)) && error("Measure parameter cannot be finite.")
    keys(parameter_refs) == keys(first(supports)) || error("Parameter references and supports must " *
                                                           "use same container type.")
    parameter_refs, supports = _convert_param_refs_and_supports(parameter_refs, supports)
    if length(coefficients) != size(supports, 2)
        error("The amount of coefficients must match the amount of " *
              "support points.")
    elseif !_allequal(group_id.(parameter_refs))
        error("Parameters must all have same group ID.")
    # elseif length(supports) == 0
    #     error("Must specify at least one support.")
    # elseif size(supports, 1) != length(parameter_refs)
    #     error("The dimensions of the support points and parameters " *
    #           "do not match.")
    end
    for i in eachindex(parameter_refs)
        set = infinite_set(parameter_refs[i])
        if !supports_in_set(supports[i, :], set)
            error("Support points violate parameter domain.")
        end
    end
    return MultiDiscreteMeasureData(parameter_refs, coefficients, supports, name,
                                    weight_function)
end

"""
    measure_name(data::AbstractMeasureData)::String

Return the measure name string stored in `data`. This is intended as an internal
function to be used with measure addition. User-defined measure data types will
need to extend this function otherwise the measure names default to `"measure"`.
"""
function measure_name(data::AbstractMeasureData)::String
    return "measure"
end

# DiscreteMeasureData and MultiDiscreteMeasureData
function measure_name(data::Union{DiscreteMeasureData, MultiDiscreteMeasureData})::String
    return data.name
end

"""
    parameter_refs(data::AbstractMeasureData)::Union{ParameterRef,
                                                     AbstractArray{<:ParameterRef}}

Return the infinite parameter reference(s) in `data`. This is intended as an
internal function to be used with measure addition. User-defined measure data types
will need to extend this function otherwise an error is thrown.
"""
function parameter_refs(data::AbstractMeasureData)
    error("Function `parameter_refs` not extended for measure data of type $(typeof(data)).")
end

# DiscreteMeasureData
function parameter_refs(data::DiscreteMeasureData)::ParameterRef
    return data.parameter_ref
end

# MultiDiscreteMeasureData
function parameter_refs(data::MultiDiscreteMeasureData)::Vector{ParameterRef}
    return data.parameter_refs
end

"""
    supports(data::AbstractMeasureData)::Vector

Return the supports stored in `data` associated with its infinite parameters.
This is intended as en internal method for measure creation and ensures any
new supports are added to parameters. User-defined measure data types should
extend this function if appropriate, otherwise an empty vector is returned.
"""
function supports(data::AbstractMeasureData)::Vector{Float64}
    return Float64[]
end

# DiscreteMeasureData
function supports(data::DiscreteMeasureData)::Vector{Float64}
    return data.supports
end

# MultiDiscreteMeasureData
function supports(data::MultiDiscreteMeasureData)::Array{Float64, 2}
    return data.supports
end

"""
    coefficients(data::AbstractMeasureData)::Vector

Return the coefficients stored in `data` associated with its expansion abstraction.
This is intended as en internal method for measure creation. User-defined measure
data types should extend this function if appropriate, otherwise an empty vector
is returned.
"""
function coefficients(data::AbstractMeasureData)::Vector{Float64}
    return Float64[]
end

# DiscreteMeasureData and MultiDiscreteMeasureData
function coefficients(data::Union{DiscreteMeasureData,
                                  MultiDiscreteMeasureData})::Vector{Float64}
    return data.coefficients
end

"""
    weight_function(data::AbstractMeasureData)::Function

Return the weight function stored in `data` associated with its expansion abstraction.
This is intended as en internal method for measure creation. User-defined measure
data types should extend this function if appropriate, otherwise the default
function [`default_weight`](@ref) is returned.
"""
function weight_function(data::AbstractMeasureData)::Function
    return default_weight
end

# DiscreteMeasureData and MultiDiscreteMeasureData
function weight_function(data::Union{DiscreteMeasureData,
                                     MultiDiscreteMeasureData})::Function
    return data.weight_function
end

# Parse the string for displaying a measure
function _make_meas_name(meas::Measure)::String
    return string(measure_name(meas.data), "(", JuMP.function_string(JuMP.REPLMode,
                                                            meas.func), ")")
end

# Used to update the model.var_to_meas and model.param_tomes fields
# this is needed to update measures if variables are deleted
function _update_var_meas_mapping(vrefs::Vector{<:GeneralVariableRef},
                                  mindex::Int)
    for vref in vrefs
        model = JuMP.owner_model(vref)
        if isa(vref, InfOptVariableRef)
            if haskey(model.var_to_meas, JuMP.index(vref))
                push!(model.var_to_meas[JuMP.index(vref)], mindex)
            else
                model.var_to_meas[JuMP.index(vref)] = [mindex]
            end
        elseif isa(vref, ParameterRef)
            if haskey(model.param_to_meas, JuMP.index(vref))
                push!(model.param_to_meas[JuMP.index(vref)], mindex)
            else
                model.param_to_meas[JuMP.index(vref)] = [mindex]
            end
        elseif isa(vref, MeasureRef)
            if haskey(model.meas_to_meas, JuMP.index(vref))
                push!(model.meas_to_meas[JuMP.index(vref)], mindex)
            else
                model.meas_to_meas[JuMP.index(vref)] = [mindex]
            end
        elseif isa(vref, ReducedInfiniteVariableRef)
            if haskey(model.reduced_to_meas, JuMP.index(vref))
                push!(model.reduced_to_meas[JuMP.index(vref)], mindex)
            else
                model.reduced_to_meas[JuMP.index(vref)] = [mindex]
            end
        end
    end
    return
end

## Used to add the measure index to param_to_meas for parameters that are used
## in the evaluation data
# DiscreteMeasureData
function _update_param_data_mapping(model::InfiniteModel,
                                    pref::ParameterRef,
                                    mindex::Int)
    if haskey(model.param_to_meas, JuMP.index(pref))
        if !(mindex in model.param_to_meas[JuMP.index(pref)])
            push!(model.param_to_meas[JuMP.index(pref)], mindex)
        end
    else
        model.param_to_meas[JuMP.index(pref)] = [mindex]
    end
    return
end

# MultiDiscreteMeasureData
function _update_param_data_mapping(model::InfiniteModel,
                                    prefs::AbstractArray{<:ParameterRef},
                                    mindex::Int)
    for pref in prefs
        if haskey(model.param_to_meas, JuMP.index(pref))
            if !(mindex in model.param_to_meas[JuMP.index(pref)])
                push!(model.param_to_meas[JuMP.index(pref)], mindex)
            end
        else
            model.param_to_meas[JuMP.index(pref)] = [mindex]
        end
    end
    return
end

## Internal functions for adding measure data supports to the parameter supports
# scalar pref
function _add_supports_to_parameters(pref::ParameterRef,
                                     supports::Vector{<:Number})
    add_supports(pref, supports)
    return
end

# array pref (from MultiDiscreteMeasureData)
function _add_supports_to_parameters(prefs::Vector{ParameterRef},
                                     supports::Array{<:Number})
    for i = eachindex(prefs)
        add_supports(prefs[i], supports[i, :])
    end
    return
end

# array fallback for extensions
function _add_supports_to_parameters(prefs::AbstractArray{<:ParameterRef},
                                     supports::Array{<:AbstractArray{<:Number}})
    for i = eachindex(supports)
        for key in keys(prefs)
            add_supports(prefs[key], supports[i][key])
        end
    end
    return
end

"""
    add_measure(model::InfiniteModel, meas::Measure)::MeasureRef

Add a measure to `model` and return the corresponding measure reference. This
operates in a manner similar to [`JuMP.add_variable`](@ref). Note this intended
as an internal method.
"""
function add_measure(model::InfiniteModel, meas::Measure)::MeasureRef
    model.next_meas_index += 1
    index = model.next_meas_index
    JuMP.check_belongs_to_model(meas.func, model)
    prefs = parameter_refs(meas.data)
    supps = supports(meas.data)
    all(JuMP.is_valid.(model, prefs)) || error("Invalid parameter dependence in " *
                                               "measure data")
    if length(supps) != 0
        _add_supports_to_parameters(prefs, supps)
    end
    vrefs = _all_function_variables(meas.func)
    _update_var_meas_mapping(vrefs, index)
    _update_param_data_mapping(model, prefs, index)
    mref = MeasureRef(model, model.next_meas_index)
    model.measures[mref.index] = meas
    JuMP.set_name(mref, _make_meas_name(meas))
    model.meas_in_objective[index] = false
    return mref
end

"""
    measure_function(mref::MeasureRef)::JuMP.AbstractJuMPScalar

Return the function associated with `mref`.
"""
function measure_function(mref::MeasureRef)::JuMP.AbstractJuMPScalar
    return JuMP.owner_model(mref).measures[JuMP.index(mref)].func
end

"""
    measure_data(mref::MeasureRef)::AbstractMeasureData

Return the measure data associated with `mref`.
"""
function measure_data(mref::MeasureRef)::AbstractMeasureData
    return JuMP.owner_model(mref).measures[JuMP.index(mref)].data
end

# TODO this is not used, but should be used when adding analytic measure checks
# Check a measure function for a particular parameter and return Bool
function _has_parameter(vrefs::Vector{<:GeneralVariableRef},
                        pref::ParameterRef)::Bool
    if _has_variable(vrefs, pref)
        return true
    end
    model = JuMP.owner_model(pref)
    relavent_ivindices = model.param_to_vars[JuMP.index(pref)]
    relavent_ivrefs = [InfiniteVariableRef(model, vindex) for vindex in relavent_ivindices]
    for ivref in relavent_ivrefs
        if _has_variable(vrefs, ivref)
            return true
        elseif used_by_reduced_variable(ivref)
            for index in model.infinite_to_reduced[JuMP.index(ivref)]
                if _has_variable(vrefs, ReducedInfiniteVariableRef(model, index))
                    return true
                end
            end
        end
    end
    return false
end

## Check if expr contains a parameter directly or via an infinite variable
# scalar pref
function _check_has_parameter(vrefs::Vector{<:GeneralVariableRef},
                              pref::ParameterRef)
    if !_has_parameter(vrefs, pref)
        error("Measure expression is not parameterized by the parameter " *
              "specified in the measure data.")
    end
    return
end

# array pref
function _check_has_parameter(vrefs::Vector{<:GeneralVariableRef},
                              prefs::AbstractArray{<:ParameterRef})
    for pref in prefs
        if !_has_parameter(vrefs, pref)
            error("Measure expression is not parameterized by the parameter " *
                  "specified in the measure data.")
        end
    end
    return
end

# Parse the model pertaining to an expression
function _model_from_expr(vrefs::Vector{<:GeneralVariableRef})
    if length(vrefs) > 0
        return JuMP.owner_model(vrefs[1])
    else
        return
    end
end

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

# DiscreteMeasureData
function measure_data_in_hold_bounds(data::DiscreteMeasureData,
                                     bounds::ParameterBounds)::Bool
    pref = parameter_refs(data)
    supps = supports(data)
    if haskey(bounds.intervals, pref)
        return supports_in_set(supps, bounds.intervals[pref])
    end
    return true
end

# MultiDiscreteMeasureData
function measure_data_in_hold_bounds(data::MultiDiscreteMeasureData,
                                     bounds::ParameterBounds)::Bool
    prefs = parameter_refs(data)
    supps = supports(data)
    for i in eachindex(prefs)
        if haskey(bounds.intervals, prefs[i])
            if !supports_in_set(supps[i, :], bounds.intervals[prefs[i]])
                return false
            end
        end
    end
    return true
end

## Check that variables don't violate the parameter bounds
# GeneralVariableRef
function _check_var_bounds(vref::GeneralVariableRef, data::AbstractMeasureData)
    return
end

# HoldVariableRef (single parameter)
function _check_var_bounds(vref::HoldVariableRef, data::AbstractMeasureData)
    bounds = parameter_bounds(vref)
    if !measure_data_in_hold_bounds(data, bounds)
        error("Measure bounds violate hold variable bounds.")
    end
    return
end

# MeasureRef
function _check_var_bounds(mref::MeasureRef, data::AbstractMeasureData)
    vrefs = _all_function_variables(measure_function(mref))
    for vref in vrefs
        _check_var_bounds(vref, data)
    end
    return
end

"""
    measure(expr::JuMP.AbstractJuMPScalar, data::AbstractMeasureData)::MeasureRef

Return a measure reference that evaluates `expr` using according to `data`. This
is the preferred method for implementing measures which follow the form:
``\\int_{p \\in P} expr(p) w(p) dp`` where ``p`` is an infinite parameter (scalar
or vector) and ``w`` is the weight function. The measure data `data` determines
how the measure is to be evaluated. Typically, the [`DiscreteMeasureData`](@ref)
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
julia> tdata = DiscreteMeasureData(t, [0.5, 0.5], [1, 2], name = "name1");

julia> xdata = DiscreteMeasureData(xs, [0.5, 0.5], [[-1, -1], [1, 1]],
                                   name = "name2");

julia> constr_RHS = @expression(model, measure(g - s + 2, tdata) + s^2)
name1(g(t) - s + 2) + s²

julia> @objective(model, Min, measure(g - 1  + measure(T, xdata), tdata))
name1(g(t) - 1 + name2(T(t, x)))
```
"""
function measure(expr::JuMP.AbstractJuMPScalar,
                 data::AbstractMeasureData)::MeasureRef
    vrefs = _all_function_variables(expr)
    model = _model_from_expr(vrefs)
    if model == nothing
        error("Expression contains no variables or parameters.")
    end
    pref = parameter_refs(data)
    if model.has_hold_bounds
        for vref in vrefs
            _check_var_bounds(vref, data)
        end
    end
    meas = Measure(expr, data)
    return add_measure(model, meas)
end

"""
    integral(expr::JuMP.AbstractJuMPScalar,
             [params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                            Nothing} = nothing,
             lb::Union{Number, AbstractArray{<:Number}, Nothing} = nothing,
             ub::Union{Number, AbstractArray{<:Number}, Nothing} = nothing];
             [eval_method::Symbol = sampling,
             num_supports::Int = 10,
             weight_func::Function = [`default_weight`](@ref),
             name = "integral",
             use_existing_supports::Bool = false,
             kwargs...])::MeasureRef

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
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(seed = true))
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
                  params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                                                            Nothing} = nothing,
                  lb::Union{Number, AbstractArray{<:Number}, Nothing} = nothing,
                  ub::Union{Number, AbstractArray{<:Number}, Nothing} = nothing;
                  kwargs...)::MeasureRef

    # collect parameters from expression if they are not provided
    if isa(params, Nothing)
        if isa(expr, MeasureRef)
            error("Nested call of intregral must specify parameters.")
        end
        params = _all_parameter_refs(expr)
        if length(params) == 0
            error("No infinite parameters in the expression.")
        elseif length(params) == 1
            params = params[1]
        else
            error("Multiple groups of parameters are in the expression. Need to " *
                  "specify one group of parameters only.")
        end
    end

    # count number of parameters check array formatting
    num_params = length(params)
    if num_params == 0
        # error if empty ParameterRef array is provided
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

"""
    integral_defaults(model::InfiniteModel)

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
function integral_defaults(model::InfiniteModel)
    return model.integral_defaults
end

"""
    set_integral_defaults(model::InfiniteModel; kwargs...)

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
function set_integral_defaults(model::InfiniteModel; kwargs...)
    merge!(model.integral_defaults, kwargs)
    return
end

"""
    expect(expr::JuMP.AbstractJuMPScalar,
           [params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                          Nothing} = nothing];
           [num_supports::Int = 10,
            use_existing_supports::Bool = false])::MeasureRef

Creates a measure that represents the expected value of an expression in
a random parameter involved in the expression. Return the [`MeasureRef`](@ref)
of the created measure.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP, Distributions; model = InfiniteModel(seed = true))
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
                params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                              Nothing} = nothing;
                num_supports::Int = 10,
                use_existing_supports::Bool = false)::MeasureRef
    # expectation measure
    return integral(expr, params, num_supports = num_supports,
                    name = "expect", use_existing_supports = use_existing_supports)
end

"""
    support_sum(expr::JuMP.AbstractJuMPScalar,
                [params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                               Nothing} = nothing])::MeasureRef

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
                     params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                                   Nothing} = nothing)::MeasureRef
    # sum measure
    return integral(expr, params, use_existing_supports = true, name = "sum",
                    eval_method = :not_default)
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
    return haskey(JuMP.owner_model(mref).meas_to_constrs, JuMP.index(mref))
end

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
    return haskey(JuMP.owner_model(mref).meas_to_meas, JuMP.index(mref))
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
    return JuMP.owner_model(mref).meas_in_objective[JuMP.index(mref)]
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
    return used_by_measure(mref) || used_by_constraint(mref) || used_by_objective(mref)
end

"""
    JuMP.delete(model::InfiniteModel, mref::MeasureRef)

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
function JuMP.delete(model::InfiniteModel, mref::MeasureRef)
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
    if prefs isa ParameterRef
        prefs = [prefs]
    end
    for pref in prefs
        filter!(e -> e != JuMP.index(mref),
                model.param_to_meas[JuMP.index(pref)])
        if length(model.param_to_meas[JuMP.index(pref)]) == 0
            delete!(model.param_to_meas, JuMP.index(pref))
        end
    end
    # Update that the variables used by it are no longer used by it
    vrefs = _all_function_variables(measure_function(mref))
    for vref in vrefs
        if isa(vref, InfOptVariableRef)
            filter!(e -> e != JuMP.index(mref),
                    model.var_to_meas[JuMP.index(vref)])
            if length(model.var_to_meas[JuMP.index(vref)]) == 0
                delete!(model.var_to_meas, JuMP.index(vref))
            end
        elseif isa(vref, ParameterRef)
            filter!(e -> e != JuMP.index(mref),
                    model.param_to_meas[JuMP.index(vref)])
            if length(model.param_to_meas[JuMP.index(vref)]) == 0
                delete!(model.param_to_meas, JuMP.index(vref))
            end
        elseif isa(vref, MeasureRef)
            filter!(e -> e != JuMP.index(mref),
                    model.meas_to_meas[JuMP.index(vref)])
            if length(model.meas_to_meas[JuMP.index(vref)]) == 0
                delete!(model.meas_to_meas, JuMP.index(vref))
            end
        elseif isa(vref, ReducedInfiniteVariableRef)
            filter!(e -> e != JuMP.index(mref),
                    model.reduced_to_meas[JuMP.index(vref)])
            if length(model.reduced_to_meas[JuMP.index(vref)]) == 0
                delete!(model.reduced_to_meas, JuMP.index(vref))
            end
        end
    end
    # delete remaining measure information
    delete!(model.meas_in_objective, JuMP.index(mref))
    delete!(model.measures, JuMP.index(mref))
    delete!(model.meas_to_name, JuMP.index(mref))
    return
end
