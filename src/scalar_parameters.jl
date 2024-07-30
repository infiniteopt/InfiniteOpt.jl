################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::IndependentParameterIndex
    )
    return IndependentParameterRef(model, index)
end

function dispatch_variable_ref(
    model::InfiniteModel,
    index::FiniteParameterIndex
    )
    return FiniteParameterRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::ScalarParameterData{<:IndependentParameter}
    )
    index =  MOIUC.add_item(model.independent_params, object)
    push!(model.param_group_indices, index)
    return index
end

function _add_data_object(
    model::InfiniteModel,
    object::ScalarParameterData{<:FiniteParameter}
    )
    return MOIUC.add_item(model.finite_params, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(model::InfiniteModel, ::Type{IndependentParameter})
    return model.independent_params
end

function _data_dictionary(model::InfiniteModel, ::Type{FiniteParameter})
    return model.finite_params
end

# Extend _data_dictionary (ref based)
function _data_dictionary(pref::IndependentParameterRef)
    return JuMP.owner_model(pref).independent_params
end

function _data_dictionary(pref::FiniteParameterRef)
    return JuMP.owner_model(pref).finite_params
end

# Extend _data_object
function _data_object(pref::ScalarParameterRef)
    object = get(_data_dictionary(pref), JuMP.index(pref), nothing)
    if isnothing(object)
        error("Invalid scalar parameter reference, cannot find ",
              "corresponding parameter in the model. This is likely ",
              "caused by using the reference of a deleted parameter.")
    end
    return object
end

################################################################################
#                             CORE OBJECT METHODS
################################################################################
"""
    core_object(pref::IndependentParameterRef)::IndependentParameter

Retrieve the underlying core [`IndependentParameter`] object for `pref`. 
This is intended as an advanced method for developers.
"""
function core_object(pref::IndependentParameterRef)
    return _data_object(pref).parameter
end

"""
    core_object(pref::FiniteParameterRef)::FiniteParameter

Retrieve the underlying core [`FiniteParameter`] object for `pref`. 
This is intended as an advanced method for developers.
"""
function core_object(pref::FiniteParameterRef)
    return _data_object(pref).parameter
end

# Extend _parameter_number
function _parameter_number(pref::IndependentParameterRef)
    return _data_object(pref).parameter_num
end

# Extend _parameter_numbers
function _parameter_numbers(pref::IndependentParameterRef)
    return [_parameter_number(pref)]
end

"""
    `parameter_group_int_index(pref::IndependentParameterRef)::Int

Return the infinite parameter group integer index corresponding to `pref`.
"""
function parameter_group_int_index(pref::IndependentParameterRef)
    return _data_object(pref).group_int_idx
end

# Extend parameter_group_int_indices
function parameter_group_int_indices(pref::IndependentParameterRef)
    return [parameter_group_int_index(pref)]
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(
    pref::ScalarParameterRef, 
    param::P, 
    data::ScalarParameterData{P}
    ) where {P <: ScalarParameter}
    data.parameter = param
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(
    pref::ScalarParameterRef, 
    param::P1, 
    data::ScalarParameterData{P2}
    )  where {P1, P2}
    new_data = ScalarParameterData(param, data.group_int_idx, data.parameter_num, 
                                   data.name, data.parameter_func_indices,
                                   data.infinite_var_indices, 
                                   data.derivative_indices, data.measure_indices,
                                   data.constraint_indices, data.in_objective,
                                   data.generative_measures,
                                   data.has_internal_supports, 
                                   data.has_generative_supports,
                                   data.has_deriv_constrs)
    _data_dictionary(pref)[JuMP.index(pref)] = new_data
    return
end

# Extend _set_core_object for ScalarParameterRefs
function _set_core_object(
    pref::ScalarParameterRef, 
    param::ScalarParameter
    )
    _adaptive_data_update(pref, param, _data_object(pref))
    return
end

################################################################################
#                            PARAMETER DEFINITION
################################################################################
# Define the default derivative evaluation method 
const DefaultDerivativeMethod = FiniteDifference()

# Check that supports don't violate the domain bounds
function _check_supports_in_bounds(
    _error::Function,
    supports::Union{<:Real, Vector{<:Real}},
    domain::AbstractInfiniteDomain
    )
    if !supports_in_domain(supports, domain)
        _error("Supports violate the domain bounds.")
    end
    return
end

"""
    build_parameter(
        _error::Function, domain::InfiniteScalarDomain;
        [num_supports::Int = 0,
        supports::Union{Real, Vector{<:Real}} = Float64[],
        sig_digits::Int = DefaultSigDigits,
        derivative_method::AbstractDerivativeMethod = DefaultDerivativeMethod]
    )::IndependentParameter

Returns a [`IndependentParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. Errors if supports violate the
bounds associated with `domain`. This is meant to primarily serve as a
helper method for [`@infinite_parameter`](@ref). Here `derivative_method` 
specifies the numerical evalution method that will be applied to derivatives that 
are taken with respect to this infinite parameter.

**Example**
```julia-repl
julia> param = build_parameter(error, IntervalDomain(0, 3), supports = Vector(0:3));
```
"""
function build_parameter(
    _error::Function,
    domain::InfiniteScalarDomain;
    num_supports::Int = 0,
    supports::Union{Real, Vector{<:Real}} = Float64[],
    sig_digits::Int = DefaultSigDigits,
    derivative_method::AbstractDerivativeMethod = DefaultDerivativeMethod,
    extra_kwargs...
    )
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    domain = round_domain(domain, sig_digits)
    label = UserDefined
    length_supports = length(supports)
    if !isempty(supports)
        supports = round.(supports, sigdigits = sig_digits)
        _check_supports_in_bounds(_error, supports, domain)
        num_supports == 0 || @warn("Ignoring `num_supports` since `supports` is not empty.")
    elseif num_supports != 0
        supports, label = generate_support_values(domain, num_supports = num_supports,
                                                  sig_digits = sig_digits)
    end
    supports_dict = DataStructures.SortedDict{Float64, Set{DataType}}(
                                            i => Set([label]) for i in supports)
    if length_supports != 0 && (length(supports_dict) != length_supports)
        @warn("Support points are not unique, eliminating redundant points.")
    end
    return IndependentParameter(domain, supports_dict, sig_digits, derivative_method,
                                generative_support_info(derivative_method))
end

# Fallback for bad domain types 
function build_parameter(
    _error::Function, 
    domain::AbstractInfiniteDomain,
    kwargs...
    )
    _error("Expected scalar infinite domain for each independent parameter, ",
           "but got a domain of type `$(domain)`. If you are trying to use an ",
           "`InfiniteArrayDomain`, try setting `independent = false`.")
end

"""
    build_parameter(_error::Function, value::Real)::FiniteParameter

Returns a [`FiniteParameter`](@ref) given the appropriate information.
This is analagous to `JuMP.build_variable`. This is meant to primarily serve as
a helper method for [`@finite_parameter`](@ref).

**Example**
```jldoctest; setup = :(using InfiniteOpt)
julia> build_parameter(error, 1)
FiniteParameter(1.0)
```
"""
function build_parameter(
    _error::Function, 
    value::Real;
    extra_kwargs...
    )
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    return FiniteParameter(value)
end

# Generic fallback
function build_parameter(_error::Function, arg, kwargs...)
    _error("Unexpected input given for the value of the parameter.")
end

"""
    add_parameter(model::InfiniteModel, p::IndependentParameter,
                  [name::String = ""])::GeneralVariableRef

Returns a [`GeneralVariableRef`](@ref) associated with the parameter `p` that is added 
to `model`. This adds a parameter to the model in a manner similar to 
`JuMP.add_variable`. This is used to add parameters with the use of 
[`@infinite_parameter`](@ref). 
[`build_parameter`](@ref build_parameter(::Function, ::InfiniteScalarDomain)) 
should be used to construct `p`.

**Example**
```julia-repl
julia> p = build_parameter(error, IntervalDomain(0, 3), supports = Vector(0:3));

julia> param_ref = add_parameter(model, p, "name")
name
```
"""
function add_parameter(
    model::InfiniteModel, 
    p::IndependentParameter,
    name::String = ""
    )
    group_int_idx = length(parameter_group_indices(model)) + 1
    param_num = model.last_param_num += 1
    data_object = ScalarParameterData(p, group_int_idx, param_num, name)
    obj_index = _add_data_object(model, data_object)
    model.name_to_param = nothing
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

"""
    add_parameter(model::InfiniteModel, p::FiniteParameter,
                  [name::String = ""])::GeneralVariableRef

Returns a [`GeneralVariableRef`](@ref) associated with the parameter `p` that is 
added to `model`. This adds a parameter to the model in a manner similar to
`JuMP.add_variable`. This is to add parameters with the use of 
[`@finite_parameter`](@ref). 
[`build_parameter`](@ref build_parameter(::Function, ::Real)) should be used to
construct `p`.

**Example**
```julia-repl
julia> p = build_parameter(error, 42);

julia> param_ref = add_parameter(model, p, "name")
name
```
"""
function add_parameter(
    model::InfiniteModel, 
    p::FiniteParameter,
    name::String = ""
    )
    data_object = ScalarParameterData(p, -1, -1, name)
    obj_index = _add_data_object(model, data_object)
    model.name_to_param = nothing
    return GeneralVariableRef(model, obj_index.value, typeof(obj_index))
end

################################################################################
#                           PARAMETER DEPENDENCIES
################################################################################
# Extend _infinite_variable_dependencies
function _infinite_variable_dependencies(pref::ScalarParameterRef)
    return _data_object(pref).infinite_var_indices
end

# Extend _parameter_function_dependencies
function _parameter_function_dependencies(pref::ScalarParameterRef)
    return _data_object(pref).parameter_func_indices
end

# Extend _derivative_dependencies
function _derivative_dependencies(pref::ScalarParameterRef)
    return _data_object(pref).derivative_indices
end

# Extend _measure_dependencies
function _measure_dependencies(pref::ScalarParameterRef)
    return _data_object(pref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(pref::ScalarParameterRef)
    return _data_object(pref).constraint_indices
end

# Extend _generative_measures
function _generative_measures(pref::ScalarParameterRef)
    return _data_object(pref).generative_measures
end

################################################################################
#                             USED_BY FUNCTIONS
################################################################################
"""
    used_by_infinite_variable(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by an infinite variable or false otherwise.

**Example**
```julia-repl
julia> used_by_infinite_variable(t)
true
```
"""
function used_by_infinite_variable(pref::IndependentParameterRef)
    return !isempty(_infinite_variable_dependencies(pref))
end

# FiniteParameter
used_by_infinite_variable(pref::FiniteParameterRef)::Bool = false

"""
    used_by_parameter_function(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by an infinite parameter function or false otherwise.

**Example**
```julia-repl
julia> used_by_parameter_function(t)
false
```
"""
function used_by_parameter_function(pref::IndependentParameterRef)
    return !isempty(_parameter_function_dependencies(pref))
end

# FiniteParameter
used_by_parameter_function(pref::FiniteParameterRef) = false

"""
    used_by_measure(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used by a measure or false otherwise.

**Example**
```julia-repl
julia> used_by_measure(t)
false
```
"""
function used_by_measure(pref::ScalarParameterRef)
    return !isempty(_measure_dependencies(pref))
end

"""
    used_by_constraint(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used by a constraint or false otherwise.

**Example**
```julia-repl
julia> used_by_constraint(t)
true
```
"""
function used_by_constraint(pref::ScalarParameterRef)
    return !isempty(_constraint_dependencies(pref))
end

"""
    used_by_objective(pref::FiniteParameterRef)::Bool

Return true if `pref` is used by the objective function.

**Example**
```julia-repl
```
"""
function used_by_objective(pref::FiniteParameterRef)
    return _data_object(pref).in_objective
end

# IndependentParameter
used_by_objective(::IndependentParameterRef) = false

"""
    used_by_derivative(pref::IndependentParameterRef)::Bool

Return true if `pref` is used by a derivative or false otherwise.

**Example**
```julia-repl
julia> used_by_derivative(t)
false
```
"""
function used_by_derivative(pref::IndependentParameterRef)
    return !isempty(_derivative_dependencies(pref))
end

# FiniteParameter
used_by_derivative(::FiniteParameterRef) = false

"""
    is_used(pref::Union{IndependentParameterRef, FiniteParameterRef})::Bool

Return true if `pref` is used in the model or false otherwise.

**Example**
```julia-repl
julia> is_used(t)
true
```
"""
function is_used(pref::ScalarParameterRef)
    return used_by_measure(pref) || used_by_constraint(pref) ||
           used_by_infinite_variable(pref) || used_by_objective(pref) || 
           used_by_derivative(pref) || used_by_parameter_function(pref)
end

################################################################################
#                              NAME METHODS
################################################################################

"""
    JuMP.name(pref::Union{IndependentParameterRef, FiniteParameterRef})::String

Extend the `JuMP.name` function to accomodate infinite parameters. Returns the 
name string associated with `pref`.

**Example**
```julia-repl
julia> name(t)
"t"
```
"""
function JuMP.name(pref::ScalarParameterRef)
    object = get(_data_dictionary(pref), JuMP.index(pref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(pref::ScalarParameterRef, name::String)

Extend the `JuMP.set_name` function to accomodate infinite parameters. Set a new 
base name to be associated with `pref`.

**Example**
```julia-repl
julia> set_name(t, "time")

julia> name(t)
"time"
```
"""
function JuMP.set_name(pref::ScalarParameterRef, name::String)
    _data_object(pref).name = name
    JuMP.owner_model(pref).name_to_param = nothing
    return
end

# Get the name_to_param Dictionary
function _param_name_dict(model::InfiniteModel)
    return model.name_to_param
end

# Update name_to_param
function _update_param_name_dict(
    model::InfiniteModel,
    param_dict::MOIUC.CleverDict{K, V}
    ) where {K, V <: ScalarParameterData}
    name_dict = _param_name_dict(model)
    for (index, data_object) in param_dict
        param_name = data_object.name
        if haskey(name_dict, param_name)
            # IndependentParameterIndex(-1) is a special value that means
            # this string does not map to a unique variable name.
            name_dict[param_name] = IndependentParameterIndex(-1)
        else
            name_dict[param_name] = index
        end
    end
    model.name_to_param = name_dict
    return
end

function _update_param_name_dict(
    model::InfiniteModel,
    param_dict::MOIUC.CleverDict{K, V}
    ) where {K, V <: MultiParameterData}
    name_dict = _param_name_dict(model)
    for (index, data_object) in param_dict
        param_nums = data_object.parameter_nums
        for i in eachindex(param_nums)
            name = data_object.names[i]
            if haskey(name_dict, name)
                # IndependentParameterIndex(-1) is a special value that means
                # this string does not map to a unique variable name.
                name_dict[name] = IndependentParameterIndex(-1)
            else
                individual_index = DependentParameterIndex(index, i)
                name_dict[name] = individual_index
            end
         end
    end
    model.name_to_param = name_dict
    return
end

"""
    parameter_by_name(model::InfiniteModel,
                      name::String)::Union{GeneralVariableRef, Nothing}

Return the parameter reference assoociated with a parameter name. Errors if
multiple parameters have the same name. Returns nothing if no such name exists.

**Example**
```julia-repl
julia> parameter_by_name(model, "t")
t
```
"""
function parameter_by_name(model::InfiniteModel, name::String)
    if isnothing(_param_name_dict(model))
        model.name_to_param = Dict{String, AbstractInfOptIndex}()
        _update_param_name_dict(model, model.independent_params)
        _update_param_name_dict(model, model.dependent_params)
        _update_param_name_dict(model, model.finite_params)
    end
    index = get(_param_name_dict(model), name, nothing)
    if isnothing(index)
        return nothing
    elseif index == IndependentParameterIndex(-1)
        error("Multiple parameters have the name $name.")
    else
        return GeneralVariableRef(model, index)
    end
end

################################################################################
#                       GENERATIVE SUPPORT FUNCTIONS
################################################################################
# Extend copy for NoGenerativeSupports
function Base.copy(d::NoGenerativeSupports)
    return NoGenerativeSupports()
end

# Extend copy for UniformGenerativeInfo
function Base.copy(d::UniformGenerativeInfo)
    return UniformGenerativeInfo(copy(d.support_basis), d.label)
end

"""
    support_label(info::AbstractGenerativeInfo)::DataType 

Return the support label to be associated with generative supports produced in 
accordance with `info`. This is intended an internal method that should be 
extended for user defined types of [`AbstractGenerativeInfo`](@ref).
"""
function support_label(info::AbstractGenerativeInfo)
    error("`support_label` not defined for generative support info type " *
          "$(typeof(info)).")
end

# UniformGenerativeInfo
function support_label(info::UniformGenerativeInfo)
    return info.label
end

# NoGenerativeSupports
function support_label(info::NoGenerativeSupports)
    return _NoLabel
end

"""
    generative_support_info(pref::IndependentParameterRef)::AbstractGenerativeInfo

Return the generative support information associated with `pref`.
"""
function generative_support_info(pref::IndependentParameterRef)
    return core_object(pref).generative_supp_info
end

"""
    has_generative_supports(pref::IndependentParameterRef)::Bool

Return whether generative supports have been added to `pref` in accordance 
with its generative support info.
"""
function has_generative_supports(pref::IndependentParameterRef)
    return _data_object(pref).has_generative_supports
end

# Specify if a parameter has generative supports
function _set_has_generative_supports(pref::IndependentParameterRef, status::Bool)
    _data_object(pref).has_generative_supports = status 
    return
end

# Reset (remove) the generative supports if needed 
function _reset_generative_supports(pref::IndependentParameterRef)
    if has_generative_supports(pref)
        label = support_label(generative_support_info(pref))
        delete_supports(pref, label = label) # this also calls _set_has_generative_supports
    end
    return
end

# Specify the generative_support_info
function _set_generative_support_info(
    pref::IndependentParameterRef,
    info::AbstractGenerativeInfo
    )
    sig_digits = significant_digits(pref)
    method = derivative_method(pref)
    domain = _parameter_domain(pref)
    supps = _parameter_supports(pref)
    new_param = IndependentParameter(domain, supps, sig_digits, method, info)
    _reset_generative_supports(pref)
    _set_core_object(pref, new_param)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    make_generative_supports(info::AbstractGenerativeInfo,
                             pref::IndependentParameterRef,
                             existing_supps::Vector{Float64}
                             )::Vector{Float64}

Generate the generative supports for `pref` in accordance with `info` and the 
`existing_supps` that `pref` has. The returned supports should not include 
`existing_supps`. This is intended as internal method to enable 
[`add_generative_supports`](@ref) and should be extended for any user defined 
`info` types that are created to enable new measure and/or derivative evaluation 
techniques that require the creation of generative supports.
"""
function make_generative_supports(info::AbstractGenerativeInfo, pref, supps)
    error("`make_generative_supports` is not defined for generative support " * 
          "info of type $(typeof(info)).")
end

# UniformGenerativeInfo
function make_generative_supports(info::UniformGenerativeInfo, pref, supps)
    # collect the preliminaries
    basis = info.support_basis
    num_internal = length(basis)
    num_existing = length(supps)
    num_existing <= 1 && error("$(pref) does not have enough supports for " *
                                "creating generative supports.")
    internal_nodes = Vector{Float64}(undef, num_internal * (num_existing - 1))
    # generate the internal node supports
    for i in Iterators.take(eachindex(supps), num_existing - 1)
        lb = supps[i]
        ub = supps[i+1]
        internal_nodes[(i-1)*num_internal+1:i*num_internal] = basis * (ub - lb) .+ lb
    end
    return internal_nodes
end

## Define internal dispatch methods for adding generative supports
# AbstractGenerativeInfo
function _add_generative_supports(pref, info::AbstractGenerativeInfo)
    if !has_generative_supports(pref)
        existing_supps = supports(pref, label = All)
        supps = make_generative_supports(info, pref, existing_supps)
        add_supports(pref, supps, label = support_label(info))
        _set_has_generative_supports(pref, true)
    end
    return
end

# NoGenerativeSupports
function _add_generative_supports(pref, info::NoGenerativeSupports)
    return
end

"""
    add_generative_supports(pref::IndependentParameterRef)::Nothing

Create generative supports for `pref` if needed in accordance with its 
generative support info using [`make_generative_supports`](@ref) and add them to 
`pref`. This is intended as an internal function, but can be useful user defined 
transformation backend extensions that utlize our support system.
"""
function add_generative_supports(pref::IndependentParameterRef)
    info = generative_support_info(pref)
    _add_generative_supports(pref, info)
    return
end

################################################################################
#                        DERIVATIVE METHOD FUNCTIONS
################################################################################
# Determine if any derivatives have derivative constraints
function has_derivative_constraints(pref::IndependentParameterRef)
    return _data_object(pref).has_deriv_constrs
end

# Make update function for whether it has derivative supports 
function _set_has_derivative_constraints(pref::IndependentParameterRef, status::Bool)
    _data_object(pref).has_deriv_constrs = status
    return
end

"""
    derivative_method(pref::IndependentParameterRef)::AbstractDerivativeMethod

Returns the numerical derivative evaluation method employed with `pref` when it 
is used as an operator parameter in a derivative.

**Example**
```julia-repl
julia> derivative_method(pref) 
FiniteDifference(Backward, true)
```
"""
function derivative_method(pref::IndependentParameterRef)
    return core_object(pref).derivative_method
end

# Make method to reset derivative constraints (supports are handled separately)
function _reset_derivative_constraints(
    pref::Union{IndependentParameterRef, DependentParameterRef}
    )
    if has_derivative_constraints(pref)
        @warn("Support/method changes will invalidate existing derivative evaluation " *
              "constraints that have been added to the InfiniteModel. Thus, " *
              "these are being deleted.")
        for idx in _derivative_dependencies(pref)
            delete_derivative_constraints(DerivativeRef(JuMP.owner_model(pref), idx))
        end
        _set_has_derivative_constraints(pref, false)
    end
    return
end

"""
    set_derivative_method(pref::IndependentParameterRef, 
                          method::AbstractDerivativeMethod)::Nothing

Specfies the desired derivative evaluation method `method` for derivatives that are 
taken with respect to `pref`. Any internal supports exclusively associated with 
the previous method will be deleted. Also, if any derivatives were evaluated 
manually, the associated derivative evaluation constraints will be deleted. Errors 
if new derivative method generates supports that are incompatible with existing 
measures.

**Example**
```julia-repl
julia> set_derivative_method(d, OrthogonalCollocation(2))

```
"""
function set_derivative_method(
    pref::IndependentParameterRef,
    method::NonGenerativeDerivativeMethod
    )
    old_param = core_object(pref)
    domain = _parameter_domain(pref)
    supps = _parameter_supports(pref)
    sig_figs = significant_digits(pref)
    if isempty(_generative_measures(pref))
        _reset_generative_supports(pref)
        new_param = IndependentParameter(domain, supps, sig_figs, method, 
                                         NoGenerativeSupports())
    else
        info = generative_support_info(pref)
        new_param = IndependentParameter(domain, supps, sig_figs, method, info)
    end
    _reset_derivative_constraints(pref)
    _set_core_object(pref, new_param)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

# GenerativeDerivativeMethod
function set_derivative_method(
    pref::IndependentParameterRef, 
    method::GenerativeDerivativeMethod
    )
    new_info = generative_support_info(method)
    old_info = generative_support_info(pref)
    if !isempty(_generative_measures(pref)) && new_info != old_info 
        error("Generative derivative method conflicts with existing generative " *
              "measures.")
    end
    old_param = core_object(pref)
    domain = _parameter_domain(pref)
    supps = _parameter_supports(pref)
    sig_figs = significant_digits(pref)
    new_param = IndependentParameter(domain, supps, sig_figs, method, new_info)
    _reset_derivative_constraints(pref)
    _reset_generative_supports(pref)
    _set_core_object(pref, new_param)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

################################################################################
#                            INFINITE DOMAIN FUNCTIONS
################################################################################
# Internal functions
function _parameter_domain(pref::IndependentParameterRef)
    return core_object(pref).domain
end
function _update_parameter_domain(
    pref::IndependentParameterRef,
    domain::AbstractInfiniteDomain
    )
    # old supports will always be discarded
    sig_digits = significant_digits(pref)
    method = derivative_method(pref)
    info = generative_support_info(pref)
    new_param = IndependentParameter(domain, DataStructures.SortedDict{Float64, Set{DataType}}(),
                                     sig_digits, method, info)
    _set_core_object(pref, new_param)
    _reset_derivative_constraints(pref)
    _set_has_generative_supports(pref, false)
    _set_has_internal_supports(pref, false)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    infinite_domain(pref::IndependentParameterRef)::InfiniteScalarDomain

Return the infinite domain associated with `pref`.

**Example**
```julia-repl
julia> infinite_domain(t)
[0, 1]
```
"""
function infinite_domain(pref::IndependentParameterRef)
    return _parameter_domain(pref)
end

"""
    set_infinite_domain(pref::IndependentParameterRef,
                     domain::InfiniteScalarDomain)::Nothing

Reset the infinite domain of `pref` with another `InfiniteScalarDomain`. An error will 
be thrown if `pref` is being used by some measure.

**Example**
```julia-repl
julia> set_infinite_domain(t, IntervalDomain(0, 2))

julia> infinite_domain(t)
[0, 2]
```
"""
function set_infinite_domain(
    pref::IndependentParameterRef,
    domain::InfiniteScalarDomain
    )
    if used_by_measure(pref)
        error("$pref is used by a measure so changing its " *
              "infinite domain is not allowed.")
    end
    _update_parameter_domain(pref, domain)
    return
end

################################################################################
#                        LOWER/UPPER BOUND FUNCTIONS
################################################################################
"""
    JuMP.has_lower_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_lower_bound` function to accomodate infinite parameters.
Return true if the domain associated with `pref` has a defined lower bound or if a
lower bound can be found. Extensions with user-defined infinite domain types
should extend `JuMP.has_lower_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_lower_bound(t)
true
```
"""
function JuMP.has_lower_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    return JuMP.has_lower_bound(domain)
end

"""
    JuMP.lower_bound(pref::IndependentParameterRef)::Real

Extend the `JuMP.lower_bound` function to accomodate infinite parameters.
Returns the lower bound associated with the infinite domain. Errors if such a bound
is not well-defined.

**Example**
```julia-repl
julia> lower_bound(t)
0.0
```
"""
function JuMP.lower_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    return JuMP.lower_bound(domain)
end

"""
    JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_lower_bound` function to accomodate infinite parameters.
Updates the infinite domain lower bound if such an operation is supported. Set
extensions that seek to employ this should extend
`JuMP.set_lower_bound(domain::NewType, lower::Number)`.

**Example**
```julia-repl
julia> set_lower_bound(t, -1)

julia> lower_bound(t)
-1.0
```
"""
function JuMP.set_lower_bound(pref::IndependentParameterRef, lower::Real)
    domain = _parameter_domain(pref)
    new_domain = JuMP.set_lower_bound(domain, lower)
    _update_parameter_domain(pref, new_domain)
    return
end

"""
    JuMP.has_upper_bound(pref::IndependentParameterRef)::Bool

Extend the `JuMP.has_upper_bound` function to accomodate infinite parameters.
Return true if the domain associated with `pref` has a defined upper bound or if a
upper bound can be found. Extensions with user-defined domains should extend
`JuMP.has_upper_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_upper_bound(t)
true
```
"""
function JuMP.has_upper_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    return JuMP.has_upper_bound(domain)
end

"""
    JuMP.upper_bound(pref::IndependentParameterRef)::Real

Extend the `JuMP.upper_bound` function to accomodate infinite parameters.
Returns the upper bound associated with the infinite domain. Errors if such a bound
is not well-defined. Extensions with user-defined domain types should extend
`JuMP.has_upper_bound(domain::NewType)` and `JuMP.upper_bound(domain::NewType)` if
appropriate.

**Example**
```julia-repl
julia> upper_bound(t)
1.0
```
"""
function JuMP.upper_bound(pref::IndependentParameterRef)
    domain = _parameter_domain(pref)
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    return JuMP.upper_bound(domain)
end

"""
    JuMP.set_upper_bound(pref::IndependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_upper_bound` function to accomodate infinite parameters.
Updates the infinite domain upper bound if and only if it is an IntervalDomain. Errors
otherwise. Extensions with user-defined infinite domains should extend
`JuMP.set_upper_bound(domain::NewType, upper::Number)` if appropriate.

**Example**
```julia-repl
julia> set_upper_bound(t, 2)

julia> upper_bound(t)
2.0
```
"""
function JuMP.set_upper_bound(pref::IndependentParameterRef, upper::Real)
    domain = _parameter_domain(pref)
    new_domain = JuMP.set_upper_bound(domain, upper)
    _update_parameter_domain(pref, new_domain)
    return
end

################################################################################
#                               SUPPORT FUNCTIONS
################################################################################
# Internal functions
function _parameter_supports(pref::IndependentParameterRef)
    return core_object(pref).supports
end
function _parameter_support_values(pref::IndependentParameterRef)
    return collect(keys(_parameter_supports(pref)))
end
function _update_parameter_supports(
    pref::IndependentParameterRef,
    supports::DataStructures.SortedDict{Float64, Set{DataType}}
    )
    domain = _parameter_domain(pref)
    method = derivative_method(pref)
    sig_figs = significant_digits(pref)
    info = generative_support_info(pref)
    new_param = IndependentParameter(domain, supports, sig_figs, method, info)
    _set_core_object(pref, new_param)
    _reset_derivative_constraints(pref)
    _set_has_generative_supports(pref, false)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    has_internal_supports(pref::Union{IndependentParameterRef, DependentParameterRef})::Bool

Indicate if `pref` has internal supports that will be hidden from the user by 
default. 
"""
function has_internal_supports(
    pref::Union{IndependentParameterRef, DependentParameterRef}
    )
    return _data_object(pref).has_internal_supports
end

# update has internal supports 
function _set_has_internal_supports(
    pref::Union{IndependentParameterRef, DependentParameterRef}, 
    status::Bool
    )
    _data_object(pref).has_internal_supports = status
    return
end

"""
    significant_digits(pref::IndependentParameterRef)::Int

Return the number of significant digits enforced on the supports of `pref`.

**Example**
```julia-repl
julia> significant_digits(t)
12
```
"""
function significant_digits(pref::IndependentParameterRef)
    return core_object(pref).sig_digits
end

"""
    num_supports(pref::IndependentParameterRef; 
                 [label::Type{<:AbstractSupportLabel} = PublicLabel])::Int

Return the number of support points associated with `pref`. By default, only the 
number of public supports are counted. The full amount can be determined by setting 
`label = All`. Moreover, the amount of labels that satisfy `label` is obtained 
using an [`AbstractSupportLabel`](@ref).

**Example**
```julia-repl
julia> num_supports(t)
2
```
"""
function num_supports(
    pref::IndependentParameterRef; 
    label::Type{<:AbstractSupportLabel} = PublicLabel
    )
    supports_dict = _parameter_supports(pref)
    if label == All || (!has_internal_supports(pref) && label == PublicLabel)
        return length(supports_dict)
    else
        return count(p -> any(v -> v <: label, p[2]), supports_dict)
    end
end

"""
    has_supports(pref::IndependentParameterRef)::Bool

Return true if `pref` has supports or false otherwise.

**Example**
```julia-repl
julia> has_supports(t)
true
```
"""
has_supports(pref::IndependentParameterRef) = !isempty(_parameter_supports(pref))

"""
    supports(pref::IndependentParameterRef; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel])::Vector{Float64}

Return the support points associated with `pref`. Errors if there are no
supports. Users can query just support points generated by a certain method
using the keyword argument `label`. By default, the function returns all public
support points regardless of the associated label. The full collection is given by setting 
`label = All`. Moreover, the amount of labels that satisfy `label` is obtained 
using an [`AbstractSupportLabel`](@ref).

**Example**
```julia-repl
julia> supports(t)
2-element Array{Float64,1}:
 0.0
 1.0
```
"""
function supports(
    pref::IndependentParameterRef; 
    label::Type{<:AbstractSupportLabel} = PublicLabel
    )
    if label == All || (!has_internal_supports(pref) && label == PublicLabel)
        return _parameter_support_values(pref)
    else
        return findall(x -> any(v -> v <: label, x), _parameter_supports(pref))
    end
end

# Return a matrix os supports when given a vector of IndependentParameterRefs (for measures)
function supports(
    prefs::Vector{IndependentParameterRef};
    label::Type{<:AbstractSupportLabel} = PublicLabel,
    use_combinatorics::Bool = true
    )
    # generate the support matrix considering all the unique combinations
    if use_combinatorics 
        supp_list = Tuple(supports(p, label = label) for p in prefs)
        inds = CartesianIndices(ntuple(i -> 1:length(supp_list[i]), length(prefs)))
        supps = Matrix{Float64}(undef, length(prefs), length(inds))
        for (k, idx) in enumerate(inds) 
            supps[:, k] = [supp_list[i][j] for (i, j) in enumerate(idx.I)]
        end
        return supps
    # generate the support matrix while negating the unique combinations
    else 
        num_supps = num_supports(first(prefs), label = label)
        trans_supps = Matrix{Float64}(undef, num_supps, length(prefs))
        for i in eachindex(prefs)
            supp = supports(prefs[i], label = label)
            if length(supp) != num_supps
                error("Cannot simultaneously query the supports of multiple " *
                      "independent parameters if the support dimensions do not match " *
                      "while ignoring the combinatorics. Try setting `use_combinatorics = true`.")
            else
                @inbounds trans_supps[:, i] = supp
            end
        end
        return permutedims(trans_supps)
    end
end

"""
    set_supports(pref::IndependentParameterRef, supports::Vector{<:Real};
                 [force::Bool = false,
                 label::Type{<:AbstractSupportLabel} = UserDefined]
                 )::Nothing

Specify the support points for `pref`. Errors if the supports violate the bounds
associated with the infinite domain. Warns if the points are not unique. If `force`
this will overwrite exisiting supports otherwise it will error if there are
existing supports.

**Example**
```julia-repl
julia> set_supports(t, [0, 1])

julia> supports(t)
2-element Array{Int64,1}:
 0
 1
```
"""
function set_supports(
    pref::IndependentParameterRef,
    supports::Vector{<:Real};
    force::Bool = false,
    label::Type{<:AbstractSupportLabel} = UserDefined
    )
    if has_supports(pref) && !force
        error("Unable set supports for $pref since it already has supports." *
              " Consider using `add_supports` or use `force = true` to " *
              "overwrite the existing supports.")
    end
    domain = _parameter_domain(pref)
    supports = round.(supports, sigdigits = significant_digits(pref))
    _check_supports_in_bounds(error, supports, domain)
    supports_dict = DataStructures.SortedDict{Float64, Set{DataType}}(
                                            i => Set([label]) for i in supports)
    if length(supports_dict) != length(supports)
        @warn("Support points are not unique, eliminating redundant points.")
    end
    _update_parameter_supports(pref, supports_dict)
    _set_has_internal_supports(pref, label <: InternalLabel)
    return
end

"""
    add_supports(pref::IndependentParameterRef,
                 supports::Union{Real, Vector{<:Real}};
                 [label::Type{<:AbstractSupportLabel} = UserDefined])::Nothing

Add additional support points for `pref` with identifying label `label`.

**Example**
```julia-repl
julia> add_supports(t, 0.5)

julia> supports(t)
3-element Array{Float64,1}:
 0.0
 0.5
 1.0

julia> add_supports(t, [0.25, 1])

julia> supports(t)
4-element Array{Float64,1}:
 0.0
 0.25
 0.5
 1.0
```
"""
function add_supports(
    pref::IndependentParameterRef,
    supports::Union{Real, Vector{<:Real}};
    label::Type{<:AbstractSupportLabel} = UserDefined,
    check::Bool = true)
    domain = infinite_domain(pref)
    supports = round.(supports, sigdigits = significant_digits(pref))
    check && _check_supports_in_bounds(error, supports, domain)
    supports_dict = _parameter_supports(pref)
    added_new_support = false
    for s in supports
        if haskey(supports_dict, s)
            push!(supports_dict[s], label)
        else
            supports_dict[s] = Set([label])
            added_new_support = true
        end
    end
    if label <: InternalLabel
        _set_has_internal_supports(pref, true)
    end
    if added_new_support
        _reset_derivative_constraints(pref)
        _reset_generative_supports(pref)
        if is_used(pref)
            set_transformation_backend_ready(JuMP.owner_model(pref), false)
        end
    end
    return
end

"""
    delete_supports(pref::IndependentParameterRef; 
                    [label::Type{<:AbstractSupportLabel} = All])::Nothing

Delete the support points for `pref`. If `label != All` then delete `label` and 
any supports that solely depend on it.

**Example**
```julia-repl
julia> delete_supports(t)

julia> supports(t)
ERROR: Parameter t does not have supports.
```
"""
function delete_supports(
    pref::IndependentParameterRef; 
    label::Type{<:AbstractSupportLabel} = All
    )
    supp_dict = _parameter_supports(pref)
    if has_derivative_constraints(pref)
        @warn("Deleting supports invalidated derivative evaluations. Thus, these " * 
              "are being deleted as well.")
        for idx in _derivative_dependencies(pref)
            delete_derivative_constraints(DerivativeRef(JuMP.owner_model(pref), idx))
        end
        _set_has_derivative_constraints(pref, false)
    end
    if label == All
        if used_by_measure(pref)
            error("Cannot delete the supports of $pref since it is used by " *
                  "a measure.")
        end
        empty!(supp_dict)
        _set_has_generative_supports(pref, false)
        _set_has_internal_supports(pref, false)
    else
        if has_generative_supports(pref) && support_label(generative_support_info(pref)) != label
            label = Union{label, support_label(generative_support_info(pref))}
        end
        _set_has_generative_supports(pref, false)
        filter!(p -> !all(v -> v <: label, p[2]), supp_dict)
        for (k, v) in supp_dict 
            filter!(l -> !(l <: label), v)
        end
        if has_internal_supports(pref) && num_supports(pref, label = InternalLabel) == 0
            _set_has_internal_supports(pref, false)
        end
    end
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

# Make dispatch for an array of parameters 
function delete_supports(
    prefs::AbstractArray{<:IndependentParameterRef}; 
    label::Type{<:AbstractSupportLabel} = All
    )
    delete_supports.(prefs, label = label)
    return
end

"""
    parameter_value(pref::FiniteParameterRef)::Float64

Return the value of a finite parameter reference `pref`. Errors if it is
an infinite parameter.

**Example**
```julia-repl
julia> value(cost)
42.0
```
"""
function parameter_value(pref::FiniteParameterRef)
    return core_object(pref).value
end

"""
    JuMP.set_value(pref::FiniteParameterRef, value::Real)::Nothing

Set the value of `pref` so long as it is a finite parameter. Errors if it is
an infinite parameter.

**Example**
```julia-repl
julia> set_value(cost, 27)

julia> value(cost)
27.0
```
"""
function JuMP.set_value(pref::FiniteParameterRef, value::Real)
    _data_object(pref).parameter = FiniteParameter(value)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    fill_in_supports!(pref::IndependentParameterRef;
                      [num_supports::Int = DefaultNumSupports])::Nothing

Automatically generate support points for a particular independent parameter `pref`.
Generating `num_supports` for the parameter. The supports are generated uniformly
if the underlying infinite domain is an `IntervalDomain` or they are generating randomly
accordingly to the distribution if the domain is a `UniDistributionDomain`.
Will add nothing if there are supports
and `modify = false`. Extensions that use user defined domain types should extend
[`generate_and_add_supports!`](@ref) and/or [`generate_support_values`](@ref)
as needed. Errors if the infinite domain type is not recognized.

**Example**
```julia-repl
julia> fill_in_supports!(x, num_supports = 4)

julia> supports(x)
4-element Array{Number,1}:
 0.0
 0.333
 0.667
 1.0

```
"""
function fill_in_supports!(
    pref::IndependentParameterRef;
    num_supports::Int = DefaultNumSupports,
    modify::Bool = true
    )
    domain = infinite_domain(pref)
    current_amount = length(_parameter_supports(pref))
    if (modify || current_amount == 0) && current_amount < num_supports
        generate_and_add_supports!(pref, domain,
                                   num_supports = num_supports - current_amount,
                                   adding_extra = (current_amount > 0))
    end
    return
end

"""
    generate_and_add_supports!(pref::IndependentParameterRef,
                               domain::AbstractInfiniteDomain,
                               [method::Type{<:AbstractSupportLabel}];
                               [num_supports::Int = DefaultNumSupports])::Nothing

Generate supports for independent parameter `pref` via [`generate_support_values`](@ref)
and add them to `pref`. This is intended as an extendable internal method for
[`fill_in_supports!`](@ref fill_in_supports!(::IndependentParameterRef)).
Most extensions that empoy user-defined infinite domains can typically enable this
by extending [`generate_support_values`](@ref). Errors if the infinite domain type
is not recognized.
"""
function generate_and_add_supports!(
    pref::IndependentParameterRef,
    domain::AbstractInfiniteDomain;
    num_supports::Int = DefaultNumSupports,
    adding_extra::Bool = false)
    sig_digits = significant_digits(pref)
    if isa(domain, IntervalDomain) && adding_extra
        supports, label = generate_support_values(domain, MCSample,
                                                  num_supports = num_supports,
                                                  sig_digits = sig_digits)
    else
        supports, label = generate_supports(domain,
                                            num_supports = num_supports,
                                            sig_digits = sig_digits)
    end
    add_supports(pref, supports, label = label)
    return
end

# Dispatch with method 
function generate_and_add_supports!(
    pref::IndependentParameterRef,
    domain::AbstractInfiniteDomain,
    method::Type{<:AbstractSupportLabel};
    num_supports::Int = DefaultNumSupports,
    adding_extra::Bool = false)
    sig_digits = significant_digits(pref)
    supports, label = generate_supports(domain, method,
                                        num_supports = num_supports,
                                        sig_digits = sig_digits)
    add_supports(pref, supports, label = label)
    return
end

################################################################################
#                               DELETE FUNCTIONS
################################################################################
# Check if parameter is used by measure data and error if it is to prevent bad
# deleting behavior
function _check_param_in_data(pref::GeneralVariableRef, data::AbstractMeasureData)
    prefs = parameter_refs(data)
    if isequal(pref, prefs) || any(isequal(pref), prefs)
        error("Unable to delete `$pref` since it is used to evaluate measures.")
    end
    return
end

# Update the dependent measures
function _update_measures(model::InfiniteModel, pref::GeneralVariableRef)
    for mindex in _measure_dependencies(pref)
        mref = dispatch_variable_ref(model, mindex)
        func = measure_function(mref)
        if func isa GeneralVariableRef
            data = measure_data(mref)
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_meas = Measure(new_func, data, Int[], Int[], true)
            _set_core_object(mref, new_meas)
        else
            _remove_variable(func, pref)
        end
    end
    return
end

# Update the dependent constraints
function _update_constraints(model::InfiniteModel, pref::GeneralVariableRef)
    for cindex in copy(_constraint_dependencies(pref))
        cref = InfOptConstraintRef(model, cindex)
        func = JuMP.jump_function(JuMP.constraint_object(cref))
        if func isa GeneralVariableRef
            set = JuMP.moi_set(JuMP.constraint_object(cref))
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            new_constr = JuMP.ScalarConstraint(new_func, set)
            _set_core_object(cref, new_constr)
            empty!(parameter_group_int_indices(cref))
        elseif func isa AbstractArray && any(isequal(pref), func)
            JuMP.delete(model, cref)
        else
            _remove_variable(func, pref)
        end
    end
    return
end

# Remove given object/parameter number and update the list
function _update_number_list(nums::Vector{Int}, list::Vector{Int})
    filter!(e -> !(e in nums), list)
    max_num = maximum(nums)
    for i in eachindex(list)
        if list[i] > max_num
            list[i] -= length(nums)
        end
    end
    return
end

# Update the model with the removed parameter number(s) and parameter group index
function _update_model_numbers(
    model::InfiniteModel,
    group_int_idx::Int,
    param_nums::Vector{Int}
    )
    # update the independent parameters
    for (_, object) in _data_dictionary(model, IndependentParameter)
        if object.group_int_idx > group_int_idx
            object.group_int_idx -= 1
            object.parameter_num -= length(param_nums)
        end
    end
    # update the dependent parameters
    for (_, object) in _data_dictionary(model, DependentParameters)
        if object.group_int_idx > group_int_idx
            object.group_int_idx -= 1
            object.parameter_nums = object.parameter_nums .- length(param_nums)
        end
    end
    # update the infinite parameter functions
    for (_, object) in model.param_functions
        _update_number_list([group_int_idx], object.func.group_int_idxs)
        _update_number_list(param_nums, object.func.parameter_nums)
    end
    # update the infinite variables
    for vref in JuMP.all_variables(model, InfiniteVariable)
        _update_number_list([group_int_idx], parameter_group_int_indices(vref))
        _update_number_list(param_nums, _parameter_numbers(vref))
    end
    # update the semi-infinite variables
    for vref in JuMP.all_variables(model, SemiInfiniteVariable)
        _update_number_list([group_int_idx], parameter_group_int_indices(vref))
        _update_number_list(param_nums, _parameter_numbers(vref))
    end
    # update the measures
    for mref in all_measures(model)
        _update_number_list([group_int_idx], parameter_group_int_indices(mref))
        _update_number_list(param_nums, _parameter_numbers(mref))
    end
    # update the constraints
    for (_, object) in model.constraints
        _update_number_list([group_int_idx], object.group_int_idxs)
    end
    # update the central info
    deleteat!(parameter_group_indices(model), group_int_idx)
    model.last_param_num -= length(param_nums)
    return
end

"""
    JuMP.delete(model::InfiniteModel, pref::ScalarParameterRef)::Nothing

Extend `JuMP.delete` to delete
scalar parameters and their dependencies. All variables, constraints, and
measure functions that depend on `pref` are updated to exclude it. Errors if the
parameter is used by an infinite variable or if it is contained in an 
`AbstractMeasureData` DataType that is employed by
a measure since the measure becomes invalid otherwise. Thus, measures that
contain this dependency must be deleted first. Note that
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) needs to be
extended to allow deletion of parameters when custom `AbstractMeasureData`
datatypes are used.

**Example**
```julia-repl
julia> delete(model, x)
```
"""
function JuMP.delete(model::InfiniteModel, pref::IndependentParameterRef)
    @assert JuMP.is_valid(model, pref) "Parameter reference is invalid."
    gvref = GeneralVariableRef(JuMP.owner_model(pref), JuMP.index(pref))
    # ensure deletion is okay (pref isn't used by measure data)
    for mindex in _measure_dependencies(pref)
        data = measure_data(dispatch_variable_ref(model, mindex))
        _check_param_in_data(gvref, data)
    end
    # ensure pref is not used by an infinite variable
    if used_by_infinite_variable(pref)
        error("Cannot delete `$pref` since it is used by an infinite ",
              "variable(s).")
    end
    # ensure pref is not used by a parameter function 
    if used_by_parameter_function(pref)
        error("Cannot delete `$pref` since it is used by an parameter ",
              "function(s).")
    end
    # update transformation backend status
    if is_used(pref)
        set_transformation_backend_ready(model, false)
    end
    # delete dependence of measures on pref
    _update_measures(model, gvref)
    # delete any derivatives that use pref 
    for index in _derivative_dependencies(pref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # update constraints in mapping to remove the parameter
    _update_constraints(model, gvref)
    # delete parameter information stored in model
    group_int_idx = parameter_group_int_index(pref)
    param_nums = _parameter_numbers(pref)
    _delete_data_object(pref)
    # update the parameter group integer indices and parameter numbers
    _update_model_numbers(model, group_int_idx, param_nums)
    return
end

# FiniteParameterRef
function JuMP.delete(model::InfiniteModel, pref::FiniteParameterRef)
    @assert JuMP.is_valid(model, pref) "Parameter reference is invalid."
    # update transformation backend status
    if is_used(pref)
        set_transformation_backend_ready(model, false)
    end
    gvref = GeneralVariableRef(model, JuMP.index(pref))
    # delete dependence of measures on pref
    _update_measures(model, gvref)
    # update constraints in mapping to remove the parameter
    _update_constraints(model, gvref)
    # update the objective if necessary
    if used_by_objective(pref)
        func = JuMP.objective_function(model)
        if func isa GeneralVariableRef
            new_func = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
            JuMP.set_objective_function(model, new_func)
            JuMP.set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        else
            _remove_variable(func, gvref)
        end
    end
    # delete parameter information stored in model
    _delete_data_object(pref)
    return
end
