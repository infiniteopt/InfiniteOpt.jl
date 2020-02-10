# Extend Base.copy for new variable types
Base.copy(v::MeasureRef, new_model::InfiniteModel) = MeasureRef(new_model,
                                                                v.index)

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

# Parse the string for displaying a measure
function _make_meas_name(meas::Measure)::String
    return string(meas.data.name, "(", JuMP.function_string(JuMP.REPLMode,
                                                            meas.func), ")")
end

## Internal functions for adding measure data supports to the parameter supports
# scalar pref
function _add_supports_to_parameters(pref::ParameterRef,
                                     supports::Vector{<:Number})
    add_supports(pref, supports)
    return
end

# array pref
function _add_supports_to_parameters(pref::JuMPC.SparseAxisArray{<:ParameterRef},
                                     supports::Array{<:JuMPC.SparseAxisArray{<:Number}})
    for i = eachindex(supports)
        for key in keys(pref.data)
            add_supports(pref.data[key], supports[i].data[key])
        end
    end
    return
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
                                    data::DiscreteMeasureData,
                                    mindex::Int)
    if haskey(model.param_to_meas, JuMP.index(data.parameter_ref))
        if !(mindex in model.param_to_meas[JuMP.index(data.parameter_ref)])
            push!(model.param_to_meas[JuMP.index(data.parameter_ref)], mindex)
        end
    else
        model.param_to_meas[JuMP.index(data.parameter_ref)] = [mindex]
    end
    return
end

# MultiDiscreteMeasureData
function _update_param_data_mapping(model::InfiniteModel,
                                    data::MultiDiscreteMeasureData,
                                    mindex::Int)
    for pref in data.parameter_ref
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

# Fallback
function _update_param_data_mapping(model::InfiniteModel, data::T,
                                    mindex::Int) where {T <: AbstractMeasureData}
    @warn "Unable to map parameter dependence for measure data type $T. " *
          "Parameter deletion methods should not be used."
    return
end

"""
    add_measure(model::InfiniteModel, meas::Measure)::MeasureRef

Add a measure to `model` and return the corresponding measure reference. This
operates in a manner similar to [`JuMP.add_variable`](@ref).
"""
function add_measure(model::InfiniteModel, meas::Measure)::MeasureRef
    model.next_meas_index += 1
    index = model.next_meas_index
    JuMP.check_belongs_to_model(meas.func, model)
    _add_supports_to_parameters(meas.data.parameter_ref, meas.data.supports)
    vrefs = _all_function_variables(meas.func)
    _update_var_meas_mapping(vrefs, index)
    _update_param_data_mapping(model, meas.data, index)
    mref = MeasureRef(model, model.next_meas_index)
    model.measures[mref.index] = meas
    JuMP.set_name(mref, _make_meas_name(meas))
    model.meas_in_objective[index] = false
    return mref
end

# Set a default weight function
_w(t) = 1

"""
    DiscreteMeasureData(parameter_ref::ParameterRef,
                        coefficients::Vector{<:Number},
                        supports::Vector{<:Number}; name::String = "measure",
                        weight_function::Function = w(t) = 1)::DiscreteMeasureData

Returns a `DiscreteMeasureData` object that can be utilized to define
measures using [`measure`](@ref). This accepts input for a scalar (single)
parameter. Note that `name` is used for printing purposes and a description of
the other arguments is provided in the documentation for
[`DiscreteMeasureData`](@ref). Errors if supports are out bounds or an unequal
number of supports and coefficients are given.

**Example**
```julia-repl
julia> data = DiscreteMeasureData(pref, [0.5, 0.5], [1, 2], name = "example")
DiscreteMeasureData(pref, [0.5, 0.5], [1, 2], "example", InfiniteOpt._w)
```
"""
function DiscreteMeasureData(parameter_ref::ParameterRef,
                             coefficients::Vector{<:Number},
                             supports::Vector{<:Number};
                             name::String = "measure",
                             weight_function::Function = _w
                             )::DiscreteMeasureData
    is_finite_parameter(parameter_ref) && error("Measure parameter cannot be " *
                                                "finite.")
    return DiscreteMeasureData(parameter_ref, coefficients, supports, name,
                               weight_function)
end

"""
    DiscreteMeasureData(parameter_ref::AbstractArray{<:ParameterRef},
                        coefficients::Vector{<:Number},
                        supports::Vector{<:AbstractArray{<:Number}};
                        name::String = "measure",
                        weight_function::Function = w(t) = 1
                        )::MultiDiscreteMeasureData

Returns a `MultiDiscreteMeasureData` object that can be utilized to
define measures using [`measure`](@ref). This accepts input for an array (multi)
parameter. The inner arrays in the supports vector need to match the formatting
of the array used for `parameter_ref`. Note that `name` is used for printing
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
function DiscreteMeasureData(parameter_ref::AbstractArray{<:ParameterRef},
                             coefficients::Vector{<:Number},
                             supports::Vector{<:AbstractArray};
                             name::String = "measure",
                             weight_function::Function = _w
                             )::MultiDiscreteMeasureData
    any(is_finite_parameter.(parameter_ref)) && error("Measure parameter cannot " *
                                                      " be finite.")
    supports = [convert(JuMPC.SparseAxisArray, s) for s in supports]
    parameter_ref = convert(JuMPC.SparseAxisArray, parameter_ref)
    return MultiDiscreteMeasureData(parameter_ref, coefficients, supports, name,
                                    weight_function)
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
                              pref::JuMPC.SparseAxisArray{<:ParameterRef})
    for key in keys(pref.data)
        if !_has_parameter(vrefs, pref.data[key])
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

## Check that variables don't violate the parameter bounds
# GeneralVariableRef
function _check_var_bounds(vref::GeneralVariableRef, data::AbstractMeasureData)
    return
end

# HoldVariableRef (single parameter)
function _check_var_bounds(vref::HoldVariableRef, data::DiscreteMeasureData)
    bounds = parameter_bounds(vref)
    pref = data.parameter_ref
    supports = data.supports
    if haskey(bounds.intervals, pref)
        if bounds.intervals[pref].lower_bound > minimum(supports) ||
            bounds.intervals[pref].upper_bound < maximum(supports)
            error("Measure bounds violate hold variable bounds.")
        end
    end
    return
end

# HoldVariableRef (multiple parameters)
function _check_var_bounds(vref::HoldVariableRef, data::MultiDiscreteMeasureData)
    bounds = parameter_bounds(vref)
    prefs = data.parameter_ref
    supports = data.supports
    mins = minimum(supports)
    maxs = maximum(supports)
    for key in keys(prefs)
        if haskey(bounds.intervals, prefs[key])
            if bounds.intervals[prefs[key]].lower_bound > mins[key] ||
                bounds.intervals[prefs[key]].upper_bound < maxs[key]
                error("Measure bounds violate hold variable bounds.")
            end
        end
    end
    return
end

# HoldVariableRef (fallback)
function _check_var_bounds(vref::HoldVariableRef, data::AbstractMeasureData)
    type = typeof(data)
    @warn "Unable to check if hold variables bounds are valid in measure with" *
          " custom measure data type $type."
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
    if !isa(expr, Union{InfiniteExpr, MeasureExpr, ParameterExpr})
        error("Expression must contain infinite variables, infinite " *
              "parameters, or measure references")
    end
    vrefs = _all_function_variables(expr)
    model = _model_from_expr(vrefs)
    if model == nothing
        error("Expression contains no variables.")
    end
    pref = data.parameter_ref
    _check_has_parameter(vrefs, pref)
    if model.has_hold_bounds
        for vref in vrefs
            _check_var_bounds(vref, data)
        end
    end
    meas = Measure(expr, data)
    return add_measure(model, meas)
end

"""
    measure(expr::JuMP.AbstractJuMPScalar,
            [params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                          Nothing} = nothing,
            lb::Union{Number, AbstractArray{<:Number}, Nothing},
            ub::Union{Number, AbstractArray{<:Number}, Nothing}];
            [eval_method::Function = mc_sampling, num_supports::Int = 50,
            weight_func::Function = _w, name = "measure",
            use_existing_supports::Bool = false,
            call_from_expect::Bool = false])::MeasureRef

Returns a measure reference that evaluates `expr` without using an object of
[`AbstractMeasureData`](@ref) type. Similar to the main [`measure`](@ref)
method, this function aims to implement measures of the form:
``\\int_{p \\in P} expr(p) w(p) dp`` where ``p`` is an infinite parameter (scalar
or vector) and ``w`` is the weight function. This function serves as a flexible
interface where users only have to provide necessary data about the
integration. Instead of taking an [`AbstractMeasureData`](@ref) object as input,
this function constructs the [`AbstractMeasureData`](@ref) object using some
default numerical integration schemes. By default, the function will generate
points by Monte Carlo sampling from the interval if the parameter is in an
[`IntervalSet`](@ref), or from the underlying distribution if the parameter is
in a [`DistributionSet`](@ref). If the expression involves multiple groups of
parameters, the user needs to specify the parameter. The user can also specify
lower bounds and upper bounds for the parameters, number of supports to
generate, and the function to generate the supports with. The last option makes
the method extendible for more schemes to generate supports.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(seed = true))
julia> @infinite_parameter(model, x in [0., 1.])
x

julia> @infinite_variable(model, f(x))
f(x)

julia> meas = measure(f, num_supports = 5)
(f(x))

julia> expand(meas)
0.2 f(0.8236475079774124) + 0.2 f(0.9103565379264364) + 0.2 f(0.16456579813368521) + 0.2 f(0.17732884646626457) + 0.2 f(0.278880109331201)
```
"""
function measure(expr::JuMP.AbstractJuMPScalar,
                 params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                                                             Nothing} = nothing,
                 lb::Union{Number, AbstractArray{<:Number}, Nothing} = nothing,
                 ub::Union{Number, AbstractArray{<:Number}, Nothing} = nothing;
                 eval_method::Function = mc_sampling, num_supports::Int = 50,
                 weight_func::Function = _w, name = "measure",
                 use_existing_supports::Bool = false,
                 call_from_expect::Bool = false)::MeasureRef
    # Measure function that takes non-AbstractMeasureData types
    if isa(params, Nothing)
        if isa(expr, MeasureRef)
            error("Nested call of measure must specify parameters.")
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

    if isa(params, ParameterRef)
        num_params = 1
    else
        num_params = length(params)
        if num_params == 0
            error("No infinite parameter is provided.")
        elseif num_params == 1
            params = first(params)
        else
            ids = unique(group_id.(params))
            if length(ids) > 1
                error("Multiple groups of parameters are specified.")
            end
            params = convert(JuMPC.SparseAxisArray, params)
            if isa(lb, AbstractArray)
                lb = convert(JuMPC.SparseAxisArray, lb)
            end
            if isa(ub, AbstractArray)
                ub = convert(JuMPC.SparseAxisArray, ub)
            end
        end
    end

    set = _parameter_set(first(params))
    if num_params > 1
        if isa(lb, Number)
            lb = JuMPC.SparseAxisArray(Dict(k => lb for k in keys(params)))
        end
        if isa(ub, Number)
            ub = JuMPC.SparseAxisArray(Dict(k => ub for k in keys(params)))
        end
    end

    if isa(set, IntervalSet)
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
                error("Lower bound is not less than upper bound for some parameter." *
                      " Please check the input lower bounds and upper bounds.")
            end
        end
    end

    # construct data and return measure if we use existing supports
    if use_existing_supports
        support = supports(params)
        if !isa(lb, Nothing) && !isa(ub, Nothing)
            support = [i for i in support if all(i .>= lb) && all(i .<= ub)]
        end

        # TODO: generate reasonable coefficients for given supports
        len = length(support)
        if call_from_expect
            data = DiscreteMeasureData(params, ones(len) ./ len, support,
                                       name = name, weight_function = weight_func)
        else
            data = DiscreteMeasureData(params, ones(len), support,
                                       name = name, weight_function = weight_func)
        end
        return measure(expr, data)
    end

    # construct AbstractMeasureData as data
    data = generate_measure_data(params, num_supports, lb, ub, method = eval_method,
                                 name = name, weight_func = weight_func)

    # call measure function to construct the measure
    return measure(expr, data)
end

"""
    expect(expr::JuMP.AbstractJuMPScalar,
           [params::Union{ParameterRef, AbstractArray{<:ParameterRef},
                          Nothing} = nothing];
           [num_supports::Int = 50,
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
                num_supports::Int = 50,
                use_existing_supports::Bool = false)::MeasureRef
    # expectation measure
    return measure(expr, params, num_supports = num_supports,
                   name = "expect", use_existing_supports = use_existing_supports,
                   call_from_expect = true)
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
                     params::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}
                     = nothing)::MeasureRef
    # sum measure
    return measure(expr, params, use_existing_supports = true, name = "sum")
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
