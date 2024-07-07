################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::DependentParameterIndex
    )
    return DependentParameterRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::MultiParameterData
    )
    index = MOIUC.add_item(model.dependent_params, object)
    push!(model.param_object_indices, index)
    return index
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel,
    ::Type{DependentParameters}
    )
    return model.dependent_params
end

# Extend _data_dictionary (ref based)
function _data_dictionary(pref::DependentParameterRef)
    return JuMP.owner_model(pref).dependent_params
end

# Extend _data_object
function _data_object(pref::DependentParameterRef)
    object = get(_data_dictionary(pref), JuMP.index(pref).object_index, nothing)
    if isnothing(object) 
        error("Invalid dependent parameter reference, cannot find ",
              "corresponding parameter in the model. This is likely ",
              "caused by using the reference of a deleted parameter.")
    end
    return object
end

# Extend _core_variable_object
function _core_variable_object(pref::DependentParameterRef)
    return _data_object(pref).parameters
end

# Return the number of dependent parameters involved
function _num_parameters(pref::DependentParameterRef)
    return length(_data_object(pref).names)
end

# Extend _delete_data_object
function _delete_data_object(vref::DependentParameterRef)
    delete!(_data_dictionary(vref), JuMP.index(vref).object_index)
    return
end

################################################################################
#                             PARAMETER DEFINITION
################################################################################
# Check that multi-dimensional domains are all the same
function _check_same_domain(_error::Function, domains)
    if !_allequal(domains)
        _error("Conflicting infinite domains. Only one multi-dimensional ",
               "can be specified. Otherwise, scalar domains can be used ",
               "element-wise.")
    end
    return 
end

## Use domain type dispatch to make the proper InfiniteArrayDomain
# MultiDistributionDomain with non SparseAxisArray
function _make_array_domain(
    _error::Function,
    domains::Vector{T},
    inds::Collections.ContainerIndices
    ) where {T <: MultiDistributionDomain}
    _check_same_domain(_error, domains)
    dist = first(domains).distribution
    if size(dist) != size(inds)
        _error("The dimensions of the parameters are incompatible with the ",
               "specified multi-dimensional distribution $(dist).")
    end
    return first(domains)
end

# MultiDistributionDomain with SparseAxisArray
function _make_array_domain(
    _error::Function,
    domains::Vector{<:MultiDistributionDomain},
    inds::Collections.ContainerIndices{1, <:Vector}
    )
    _error("Cannot specify multiple-dimensional distribution domain with a ",
           "`SparseAxisArray` of dependent infinite parameters.")
end

# CollectionDomain with non SparseAxisArray
function _make_array_domain(
    _error::Function,
    domains::Vector{T},
    inds::Collections.ContainerIndices
    ) where {T <: CollectionDomain}
    _check_same_domain(_error, domains)
    if length(collection_domains(first(domains))) != length(inds)
        _error("The dimensions of the parameters and the specified ",
               "`CollectionDomain` do not match.")
    end
    return first(domains)
end

# CollectionDomain with SparseAxisArray
function _make_array_domain(
    _error::Function,
    domains::Vector{<:CollectionDomain},
    inds::Collections.ContainerIndices{1, <:Vector}
    )
    _error("Cannot specify a `CollectionDomain` with a `SparseAxisArray` ",
           "of dependent infinite parameters, consider instead specifying ",
           "the `InfiniteScalarDomain` for each parameter using the `domain` ",
           "keyword and the appropriate indices.")
end

# Fallback for other InfiniteArrayDomains
function _make_array_domain(
    _error::Function,
    domains::Vector{T},
    inds::Collections.ContainerIndices
    ) where {T <: InfiniteArrayDomain}
    _check_same_domain(_error, domains)
    return first(domains)
end

# InfiniteScalarDomains
function _make_array_domain(
    _error::Function,
    domains::Vector{T},
    inds::Collections.ContainerIndices
    ) where {T <: InfiniteScalarDomain}
    return CollectionDomain(domains)
end

# Generic fallback
function _make_array_domain(_error::Function, domains, inds)
    _error("Unrecognized infinite domain format.")
end

## Process the supports format via dispatch
# Vector{<:Real}
function _process_supports(
    _error::Function, 
    supps::Vector{<:Real},
    domain,
    sig_digits
    )
    if !supports_in_domain(reshape(supps, length(supps), 1), domain) 
        _error("Support violates the infinite domain.")
    end
    supps = round.(supps, sigdigits = sig_digits)
    return Dict{Vector{Float64}, Set{DataType}}(supps => Set([UserDefined]))
end

# Vector{Vector{<:Real}}
function _process_supports(
    _error::Function, 
    vect_supps::Vector{<:Vector{<:Real}},
    domain,
    sig_digits
    )
    len = length(first(vect_supps))
    if any(length(s) != len for s in vect_supps)
        _error("Inconsistent support dimensions.")
    end
    supps = permutedims(reduce(hcat, vect_supps))
    if !supports_in_domain(supps, domain) 
        _error("Supports violate the infinite domain.")
    end
    supps = round.(supps, sigdigits = sig_digits)
    return Dict{Vector{Float64}, Set{DataType}}(s => Set([UserDefined]) 
                                                for s in eachcol(supps))
end

## Use dispatch to make the formatting of the derivative method vector 
# Valid vector 
function _process_derivative_methods(
    _error::Function,
    methods::V,
    domains
    ) where {V <: Vector{<:NonGenerativeDerivativeMethod}}
    return methods
end

# Default or invalid vector 
function _process_derivative_methods(
    _error::Function,
    methods::Vector,
    domains
    )
    if isempty(methods)
        return map(i -> DefaultDerivativeMethod, domains)
    end
    _error("Can only subtypes of `NonGenerativeDerivativeMethod` can be ",
           "used with the `derivative_method` keywords argument.") 
end

# A helper function for dependent parameter builds in @infinite_parameter
function _build_parameters(
    _error::Function,
    domains::Vector,
    orig_inds::Collections.ContainerIndices;
    num_supports::Int = 0,
    sig_digits::Int = DefaultSigDigits,
    supports::Union{Vector{<:Real}, Vector{<:Vector{<:Real}}} = Float64[],
    derivative_method::Vector = [],
    extra_kwargs...
    )
    # error with extra keywords
    for (kwarg, _) in extra_kwargs
       _error("Unrecognized keyword argument $kwarg")
    end
    # process the infinite domain
    domain = _make_array_domain(_error, domains, orig_inds)
    domain = round_domain(domain, sig_digits)
    # we have supports
    if !isempty(supports)
        supp_dict = _process_supports(_error, supports, domain, sig_digits)
    # we want to generate supports
    elseif !iszero(num_supports)
        supps, label = generate_support_values(domain, 
                                               num_supports = num_supports,
                                               sig_digits = sig_digits)
        supp_dict = Dict{Vector{Float64}, Set{DataType}}(s => Set([label]) 
                                                         for s in eachcol(supps))
    # no supports are specified
    else
        supp_dict = Dict{Vector{Float64}, Set{DataType}}()
    end
    # check the derivative methods
    methods = _process_derivative_methods(_error, derivative_method, domains)
    # make the parameter object
    return DependentParameters(domain, supp_dict, sig_digits, methods)
end

"""
    add_parameters(model::InfiniteModel,
                   params::DependentParameters,
                   names::Vector{String}
                   )::Vector{GeneralVariableRef}

Add `params` to `model` and return an appropriate container of the dependent
infinite parameter references. This is intended as an internal method for use
with [`@infinite_parameter`](@ref). However, if desired users can use this
to add a `DependentParameters` object to `model`. Here, `names` denote the name 
of each parameter. 

**Example**
```julia-repl
julia> using Distributions

julia> dist = MvNormal(ones(3)); # 3 dimensional

julia> domain = MultiDistributionDomain(dist); # 3 dimensional

julia> params = DependentParameters(domain, Dict{Vector{Float64}, Set{DatatType}}(), 10);

julia> prefs = add_parameters(model, params, ["par1", "par2", "par3"])
3-element Array{GeneralVariableRef,1}:
 par1
 par2
 par3
```
"""
function add_parameters(
    model::InfiniteModel,
    params::DependentParameters,
    names::Vector{String}
    )
    # get the number of parameters
    num_params = length(params.domain)
    # process the names
    if length(names) != num_params
        error("The amounts of names and dependent parameters do not match.")
    end
    # make the parameter model object
    group_int_idx = length(_param_object_indices(model)) + 1
    first_param_num = model.last_param_num + 1
    last_param_num = model.last_param_num += num_params
    param_nums = first_param_num:last_param_num
    data_object = MultiParameterData(params, group_int_idx, param_nums, names)
    # add the data object to the model and make the references
    obj_index = _add_data_object(model, data_object)
    # reset the name dictionary 
    model.name_to_param = nothing
    return [GeneralVariableRef(model, obj_index.value, DependentParameterIndex, i)
            for i in 1:num_params]
end

################################################################################
#                                  NAMING
################################################################################
# Get the parameter index in the DependentParameters object
_param_index(pref::DependentParameterRef) = JuMP.index(pref).param_index

"""
    JuMP.name(pref::DependentParameterRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to return the names of
infinite dependent parameters.

**Example**
```julia-repl
julia> name(pref)
"par_name"
```
"""
function JuMP.name(pref::DependentParameterRef)
    object = get(_data_dictionary(pref), JuMP.index(pref).object_index, nothing)
    return isnothing(object) ? "" : object.names[_param_index(pref)]
end

"""
    JuMP.set_name(pref::DependentParameterRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of dependent infinite parameters.

**Example**
```julia-repl
julia> set_name(vref, "par_name")

julia> name(vref)
"para_name"
```
"""
function JuMP.set_name(pref::DependentParameterRef, name::String)
    _data_object(pref).names[_param_index(pref)] = name
    JuMP.owner_model(pref).name_to_param = nothing
    return
end

################################################################################
#                           PARAMETER DEPENDENCIES
################################################################################
# Extend _infinite_variable_dependencies
function _infinite_variable_dependencies(pref::DependentParameterRef
    )::Vector{InfiniteVariableIndex}
    return _data_object(pref).infinite_var_indices
end

# Extend _parameter_function_dependencies
function _parameter_function_dependencies(pref::DependentParameterRef
    )::Vector{ParameterFunctionIndex}
    return _data_object(pref).parameter_func_indices
end

# Extend _measure_dependencies
function _measure_dependencies(pref::DependentParameterRef
    )::Vector{MeasureIndex}
    return _data_object(pref).measure_indices[_param_index(pref)]
end

# Extend _constraint_dependencies
function _constraint_dependencies(pref::DependentParameterRef
    )::Vector{InfOptConstraintIndex}
    return _data_object(pref).constraint_indices[_param_index(pref)]
end

# Extend _derivative_dependencies
function _derivative_dependencies(pref::DependentParameterRef
    )::Vector{DerivativeIndex}
    return _data_object(pref).derivative_indices[_param_index(pref)]
end

"""
    used_by_infinite_variable(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
an infinite variable.

**Example**
```julia-repl
julia> used_by_infinite_variable(pref)
true
```
"""
function used_by_infinite_variable(pref::DependentParameterRef)
    return !isempty(_infinite_variable_dependencies(pref))
end

"""
    used_by_parameter_function(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
an infinite parameter function.

**Example**
```julia-repl
julia> used_by_parameter_function(pref)
true
```
"""
function used_by_parameter_function(pref::DependentParameterRef)
    return !isempty(_parameter_function_dependencies(pref))
end

"""
    used_by_measure(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used
by a measure.

**Example**
```julia-repl
julia> used_by_measure(pref)
true
```
"""
function used_by_measure(pref::DependentParameterRef)
    return !isempty(_measure_dependencies(pref))
end

"""
    used_by_constraint(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
a constraint.

**Example**
```julia-repl
julia> used_by_constraint(pref)
false
```
"""
function used_by_constraint(pref::DependentParameterRef)
    return !isempty(_constraint_dependencies(pref))
end

"""
    used_by_derivative(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used by
a derivative.

**Example**
```julia-repl
julia> used_by_derivative(pref)
false
```
"""
function used_by_derivative(pref::DependentParameterRef)
    return !isempty(_derivative_dependencies(pref))
end

# Extend used by objective
used_by_objective(pref::DependentParameterRef) = false

"""
    is_used(pref::DependentParameterRef)::Bool

Return a `Bool` indicating if the dependent infinite parameter `pref` is used in
the model.

**Example**
```julia-repl
julia> is_used(pref)
true
```
"""
function is_used(pref::DependentParameterRef)
    return used_by_measure(pref) || used_by_constraint(pref) ||
           used_by_infinite_variable(pref) || used_by_derivative(pref) ||
           used_by_parameter_function(pref)
end

################################################################################
#                          PARAMETER OBJECT METHODS
################################################################################
# Extend _parameter_number
function _parameter_number(pref::DependentParameterRef)
    return _data_object(pref).parameter_nums[_param_index(pref)]
end

# Extend _parameter_numbers
function _parameter_numbers(pref::DependentParameterRef)
    return [_parameter_number(pref)]
end

"""
    parameter_group_int_index(pref::DependentParameterRef)::Int

Return the infinite parameter group integer index that corresponds to `pref`.
"""
function parameter_group_int_index(pref::DependentParameterRef)
    return _data_object(pref).group_int_idx
end

# Extend parameter_group_int_indices
function parameter_group_int_indices(pref::DependentParameterRef)
    return [parameter_group_int_index(pref)]
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(pref::DependentParameterRef, params::P, 
    data::MultiParameterData{P}) where {P <: DependentParameters}
    data.parameters = params
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(pref::DependentParameterRef, params::P1, 
    data::MultiParameterData{P2})  where {P1, P2}
    new_data = MultiParameterData(params, data.object_num, data.parameter_nums, 
                                  data.names, data.parameter_func_indices,
                                  data.infinite_var_indices, 
                                  data.derivative_indices, data.measure_indices,
                                  data.constraint_indices,
                                  data.has_internal_supports, 
                                  data.has_deriv_constrs)
    _data_dictionary(pref)[JuMP.index(pref).object_index] = new_data
    return
end

# Extend _set_core_variable_object
function _set_core_variable_object(pref::DependentParameterRef,
                                   params::DependentParameters)
    _adaptive_data_update(pref, params, _data_object(pref))
    return
end

################################################################################
#                        DERIVATIVE METHOD FUNCTIONS
################################################################################
# Extend fallback for dependent parameters
function has_generative_supports(pref::DependentParameterRef)::Bool
    return false
end

# Extend fallback for dependent parameters
function _set_has_generative_supports(pref::DependentParameterRef, 
                                      status::Bool)::Nothing
    return
end

# Extend add_generative_supports
function add_generative_supports(pref::DependentParameterRef)::Nothing
    return
end

# Extend has derivative constraints 
function has_derivative_constraints(pref::DependentParameterRef)::Bool 
    return _data_object(pref).has_deriv_constrs[_param_index(pref)]
end

# Extend setting if has derivative constraints 
function _set_has_derivative_constraints(pref::DependentParameterRef, 
                                         status::Bool)::Nothing 
    _data_object(pref).has_deriv_constrs[_param_index(pref)] = status
    return
end

# Get the raw derivative method vector
function _derivative_methods(pref::DependentParameterRef)
    return _core_variable_object(pref).derivative_methods
end

"""
    derivative_method(pref::DependentParameterRef)::NonGenerativeDerivativeMethod

Returns the numerical derivative evaluation method employed with `pref` when it 
is used as an operator parameter in a derivative.

**Example**
```julia-repl
julia> derivative_method(pref) 
FiniteDifference
```
"""
function derivative_method(pref::DependentParameterRef)::NonGenerativeDerivativeMethod
    return _derivative_methods(pref)[_param_index(pref)]
end

## Define helper methods for setting the derivative method efficiently
# Compatible with vector type 
function _adaptive_method_update(pref, 
    p::DependentParameters{S, M1}, 
    method::M2
    )::Nothing where {S, M1 <: NonGenerativeDerivativeMethod, M2 <: M1}
    p.derivative_methods[_param_index(pref)] = method
    return
end

# Not compatible
function _adaptive_method_update(pref, 
    p::DependentParameters{S, M1}, 
    method::M2
    )::Nothing where {S, M1, M2}
    methods = p.derivative_methods
    new_methods = [i == _param_index(pref) ? method : m 
                   for (i, m) in enumerate(methods)]
    new_params = DependentParameters(p.domain, p.supports, p.sig_digits, new_methods)
    _set_core_variable_object(pref, new_params)
    return
end

"""
    set_derivative_method(pref::DependentParameterRef, 
                          method::NonGenerativeDerivativeMethod)::Nothing

Specfies the desired derivative evaluation method `method` for derivatives that are 
taken with respect to `pref`. Errors if `method` is generative (i.e., it requires 
the definition of additional supports)

**Example**
```julia-repl
julia> set_derivative_method(d, FiniteDifference())

```
"""
function set_derivative_method(pref::DependentParameterRef, 
    method::AbstractDerivativeMethod
    )::Nothing
    if !(method isa NonGenerativeDerivativeMethod)
        error("Must specify a subtype of `NonGenerativeDerivativeMethod` for " *
              "for a dependent parameter.")
    end
    _adaptive_method_update(pref, _core_variable_object(pref), method)
    _reset_derivative_constraints(pref)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    set_all_derivative_methods(model::InfiniteModel, 
                               method::AbstractDerivativeMethod)::Nothing

Sets the desired evaluation method `method` for all the derivatives currently added 
to `model`. Note that this is done with respect to the infinite parameters. Errors 
if a generative method is specified and the model contains dependent parameters.

**Example**
```julia-repl
julia> set_all_derivative_methods(model, OrthogonalCollocation(2))

```
"""
function set_all_derivative_methods(model::InfiniteModel, 
    method::AbstractDerivativeMethod
    )::Nothing
    for pref in all_parameters(model, InfiniteParameter)
        set_derivative_method(pref, method)
    end
    return
end

################################################################################
#                           INFINITE DOMAIN METHODS
################################################################################
## Get the individual infinite domain if possible
# raw_domain
function _parameter_domain(pref::DependentParameterRef)::InfiniteArrayDomain
    return _core_variable_object(pref).domain
end

# CollectionDomain
function _parameter_domain(domain::CollectionDomain{S},
                        pref::DependentParameterRef
                        )::S where {S <: InfiniteScalarDomain}
    return collection_domains(domain)[_param_index(pref)]
end

# InfiniteArrayDomain (Fallback)
function _parameter_domain(domain::InfiniteArrayDomain, pref::DependentParameterRef)
    error("An individual infinite domain is not well-defined for $pref which " *
          "is part of a group of dependent infinite parameters that correspond " *
          "to an multi-dimensional infinite domain of type `$(typeof(domain))`.")
end

"""
    infinite_domain(pref::DependentParameterRef)::InfiniteScalarDomain

Return the infinite domain associated with the particular infinite dependent
parameter `pref` if valid. Errors if the underlying [`DependentParameters`](@ref)
object does not use a [`CollectionDomain`](@ref).

**Example**
```julia-repl
julia> infinite_domain(x[1])
[-1, 1]
```
"""
function infinite_domain(pref::DependentParameterRef)::InfiniteScalarDomain
    return _parameter_domain(_parameter_domain(pref), pref)
end

# Check that prefs are complete
function _check_complete_param_array(
    prefs::AbstractArray{<:DependentParameterRef}
    )::Nothing
    if length(prefs) != _num_parameters(first(prefs))
        error("Dimensions of parameter container and the infinite domain do not " *
              "match, ensure all related dependent parameters are included.")
    end
    return
end

"""
    infinite_domain(prefs::AbstractArray{<:DependentParameterRef})::InfiniteArrayDomain

Return the infinite domain associated with the container of infinite dependent
parameters `prefs`. Errors if the container `prefs` is incomplete.

**Example**
```julia-repl
julia> infinite_domain(x)
ZeroMeanDiagNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 1.0]
)
```
"""
function infinite_domain(prefs::AbstractArray{<:DependentParameterRef}
                      )::InfiniteArrayDomain
    _check_complete_param_array(prefs)
    return _parameter_domain(first(prefs))
end

# Update the underlying domain and delete the supports
function _update_parameter_domain(pref::DependentParameterRef,
                               new_domain::InfiniteArrayDomain)::Nothing
    old_params = _core_variable_object(pref)
    new_supports = Dict{Vector{Float64}, Set{DataType}}()
    sig_figs = significant_digits(pref)
    methods = _derivative_methods(pref)
    new_params = DependentParameters(new_domain, new_supports, sig_figs, methods)
    _set_core_variable_object(pref, new_params)
    for i in 1:length(new_domain)
        idx = DependentParameterIndex(JuMP.index(pref).object_index, i)
        p = DependentParameterRef(JuMP.owner_model(pref), idx)
        _reset_derivative_constraints(p)
    end
    _set_has_internal_supports(pref, false)
    if is_used(pref)
        set_transformation_backend_ready(JuMP.owner_model(pref), false)
    end
    return
end

"""
    set_infinite_domain(pref::DependentParameterRef,
                     domain::InfiniteScalarDomain)::Nothing

Specify the scalar infinite domain of the dependent infinite parameter `pref` to
`domain` if `pref` is part of a [`CollectionDomain`](@ref), otherwise an error is
thrown. Note this will reset/delete all the supports contained in the
underlying [`DependentParameters`](@ref) object. Also, errors if `pref` is used
by a measure.

**Example**
```julia-repl
julia> set_infinite_domain(x[1], IntervalDomain(0, 2))

julia> infinite_domain(x[1])
[0, 2]
```
"""
function set_infinite_domain(pref::DependentParameterRef,
                          domain::InfiniteScalarDomain)::Nothing
    old_domain = _parameter_domain(pref)
    if !(old_domain isa CollectionDomain)
        error("Cannot set the individual infinite domain of $pref if the " *
              "underlying domain is not a CollectionDomain.")
    elseif used_by_measure(pref)
        error("Cannot override the infinite domain of $pref since it is used by " *
              "a measure.")
    end
    param_idx = _param_index(pref)
    new_domain = CollectionDomain([i != param_idx ? collection_domains(old_domain)[i] : domain
                             for i in eachindex(collection_domains(old_domain))])
    _update_parameter_domain(pref, new_domain)
    return
end

"""
    set_infinite_domain(prefs::AbstractArray{<:DependentParameterRef},
                     domain::InfiniteArrayDomain)::Nothing

Specify the multi-dimensional infinite domain of the dependent infinite parameters
`prefs` to `domain`. Note this will reset/delete all the supports contained in the
underlying [`DependentParameters`](@ref) object. This will error if the not all
of the dependent infinite parameters are included, if any of them are used by
measures.

**Example**
```julia-repl
julia> set_infinite_domain(x, CollectionDomain([IntervalDomain(0, 1), IntervalDomain(0, 2)]))
```
"""
function set_infinite_domain(prefs::AbstractArray{<:DependentParameterRef},
                          domain::InfiniteArrayDomain)::Nothing
    if any(used_by_measure(pref) for pref in prefs)
        error("Cannot override the infinite domain of $prefs since it is used by " *
              "a measure.")
    end
    _check_complete_param_array(prefs)
    _update_parameter_domain(first(prefs), domain)
    return
end

"""
    JuMP.has_lower_bound(pref::DependentParameterRef)::Bool

Extend the `JuMP.has_lower_bound` function to accomodate a single dependent
infinite parameter.
Return true if the domain associated with `pref` has a defined lower bound or if a
lower bound can be found. Extensions with user-defined scalar infinite domain types
should extend `JuMP.has_lower_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_lower_bound(x[1])
true
```
"""
function JuMP.has_lower_bound(pref::DependentParameterRef)::Bool
    domain = _parameter_domain(pref)
    if domain isa CollectionDomain
        return JuMP.has_lower_bound(collection_domains(domain)[_param_index(pref)])
    else
        return false
    end
end

"""
    JuMP.lower_bound(pref::DependentParameterRef)::Number

Extend the `JuMP.lower_bound` function to accomodate a single dependent infinite
parameter. Returns the lower bound associated with the infinite domain. Errors if
such a bound is not well-defined.

**Example**
```julia-repl
julia> lower_bound(x[1])
0.0
```
"""
function JuMP.lower_bound(pref::DependentParameterRef)::Number
    if !JuMP.has_lower_bound(pref)
        error("Parameter $(pref) does not have a lower bound.")
    end
    return JuMP.lower_bound(infinite_domain(pref))
end

"""
    JuMP.set_lower_bound(pref::DependentParameterRef, lower::Real)::Nothing

Extend the `JuMP.set_lower_bound` function to accomodate a single dependent
infinite parameter. Updates the infinite domain lower bound if such an operation
is supported. Infinite scalar domain extensions that seek to employ this should extend
`JuMP.set_lower_bound(domain::NewType, lower::Number)`. This will call
[`set_infinite_domain`](@ref) and will error if this is not well-defined. Note
that existing supports will be deleted.

**Example**
```julia-repl
julia> set_lower_bound(t, -1)

julia> lower_bound(t)
-1.0
```
"""
function JuMP.set_lower_bound(pref::DependentParameterRef, lower::Real)::Nothing
    domain = infinite_domain(pref)
    new_domain = JuMP.set_lower_bound(domain, lower)
    set_infinite_domain(pref, new_domain)
    return
end

"""
    JuMP.has_upper_bound(pref::DependentParameterRef)::Bool

Extend the `JuMP.has_upper_bound` function to accomodate a single dependent
infinite parameter.
Return true if the domain associated with `pref` has a defined upper bound or if a
upper bound can be found. Extensions with user-defined scalar infinite domain types
should extend `JuMP.has_upper_bound(domain::NewType)`.

**Example**
```julia-repl
julia> has_upper_bound(x[1])
true
```
"""
function JuMP.has_upper_bound(pref::DependentParameterRef)::Bool
    domain = _core_variable_object(pref).domain
    if domain isa CollectionDomain
        return JuMP.has_upper_bound(collection_domains(domain)[_param_index(pref)])
    else
        return false
    end
end

"""
    JuMP.upper_bound(pref::DependentParameterRef)::Number

Extend the `JuMP.upper_bound` function to accomodate a single dependent infinite
parameter. Returns the upper bound associated with the infinite domain. Errors if
such a bound is not well-defined.

**Example**
```julia-repl
julia> upper_bound(x[1])
0.0
```
"""
function JuMP.upper_bound(pref::DependentParameterRef)::Number
    if !JuMP.has_upper_bound(pref)
        error("Parameter $(pref) does not have a upper bound.")
    end
    return JuMP.upper_bound(infinite_domain(pref))
end

"""
    JuMP.set_upper_bound(pref::DependentParameterRef, upper::Real)::Nothing

Extend the `JuMP.set_upper_bound` function to accomodate a single dependent
infinite parameter. Updates the infinite domain upper bound if such an operation
is supported. Infinite scalar domain extensions that seek to employ this should extend
`JuMP.set_upper_bound(domain::NewType, upper::Number)`. This will call
[`set_infinite_domain`](@ref) and will error if this is not well-defined. Note
that existing supports will be deleted.

**Example**
```julia-repl
julia> set_upper_bound(t, -1)

julia> upper_bound(t)
-1.0
```
"""
function JuMP.set_upper_bound(pref::DependentParameterRef, upper::Real)::Nothing
    domain = infinite_domain(pref)
    new_domain = JuMP.set_upper_bound(domain, upper)
    set_infinite_domain(pref, new_domain)
    return
end

################################################################################
#                              SUPPORT METHODS
################################################################################
# Get the raw supports
function _parameter_supports(pref::DependentParameterRef
                             )::Dict{Vector{Float64}, Set{DataType}}
    return _core_variable_object(pref).supports
end

"""
    significant_digits(pref::DependentParameterRef)::Int

Return the number of significant digits enforced on the supports of `pref`.

**Example**
```julia-repl
julia> significant_digits(x[1])
12
```
"""
function significant_digits(pref::DependentParameterRef)::Int
    return _core_variable_object(pref).sig_digits
end

"""
    num_supports(pref::DependentParameterRef; 
                 [label::Type{<:AbstractSupportLabel} = PublicLabel])::Int

Return the number of support points associated with a single dependent infinite
parameter `pref`. Specify a subset of supports via `label` to only count the
supports with `label`. By default only the amount of public supports are given, but 
the full amount is obtained via `label == All`.

**Example**
```julia-repl
julia> num_supports(x[1])
2

julia> num_supports(x[1], label = MCSample)
0
```
"""
function num_supports(pref::DependentParameterRef; 
                      label::Type{<:AbstractSupportLabel} = PublicLabel)::Int
    supp_dict = _parameter_supports(pref)
    if label == All || (!has_internal_supports(pref) && label == PublicLabel)
        return length(supp_dict)
    else
        return count(p -> any(v -> v <: label, p[2]), supp_dict)
    end
end

"""
    num_supports(prefs::AbstractArray{<:DependentParameterRef};
                 [label::Type{<:AbstractSupportLabel} = PublicLabel])::Int

Return the number of support points associated with dependent infinite
parameters `prefs`. Errors if not all from the same underlying object.
Specify a subset of supports via `label` to only count the supports with `label`.
By default only the amount of public supports are given, but the full amount is 
obtained via `label == All`.

**Example**
```julia-repl
julia> num_supports(x)
2
```
"""
function num_supports(prefs::AbstractArray{<:DependentParameterRef};
    label::Type{<:AbstractSupportLabel} = PublicLabel
    )::Int
    _check_complete_param_array(prefs)
    return num_supports(first(prefs), label = label)
end

"""
    has_supports(pref::DependentParameterRef)::Bool

Return true if `pref` has supports or false otherwise.

**Example**
```julia-repl
julia> has_supports(x[1])
true
```
"""
has_supports(pref::DependentParameterRef)::Bool = !isempty(_parameter_supports(pref))

"""
    has_supports(prefs::AbstractArray{<:DependentParameterRef})::Bool

Return true if `prefs` have supports or false otherwise. Errors if not all of the
infinite dependent parameters are from the same object.

**Example**
```julia-repl
julia> has_supports(x)
true
```
"""
function has_supports(prefs::AbstractArray{<:DependentParameterRef})::Bool
    _check_complete_param_array(prefs)
    return has_supports(first(prefs))
end

"""
    supports(pref::DependentParameterRef; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel])::Vector{Float64}

Return the support points associated with `pref`. A subset of supports can be
returned via `label` to return just the supports associated with `label`. By 
default only the public supports are given, but the full set is 
obtained via `label == All`.

**Example**
```julia-repl
julia> supports(x[1])
2-element Array{Float64,1}:
 0.0
 1.0
```
"""
function supports(pref::DependentParameterRef;
                  label::Type{<:AbstractSupportLabel} = PublicLabel)::Vector{Float64}
    supp_dict = _parameter_supports(pref)
    pindex = _param_index(pref)
    if label == All || (!has_internal_supports(pref) && label == PublicLabel)
        return Float64[supp[pindex] for supp in keys(supp_dict)]
    else
        reduced_supps = findall(e -> any(v -> v <: label, e), supp_dict)
        return Float64[supp[pindex] for supp in reduced_supps]
    end
end

"""
    supports(prefs::AbstractArray{<:DependentParameterRef};
             [label::Type{<:AbstractSupportLabel} = PublicLabel]
             )::Union{Vector{<:AbstractArray{<:Real}}, Array{Float64, 2}}

Return the support points associated with `prefs`. Errors if not all of the
infinite dependent parameters are from the same object. This will return a
matrix if `prefs` is `Vector`, otherwise a vector of arrays is returned where each 
array is a support point matching the format of `prefs`. A subset of supports can be
returned via `label` to return just the supports associated with `label`. By 
default only the public supports are given, but the full set is  obtained via 
`label == All`.

**Example**
```julia-repl
julia> supports(x) # columns are supports
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0
```
"""
function supports(prefs::AbstractArray{<:DependentParameterRef};
                  label::Type{<:AbstractSupportLabel} = PublicLabel
                  )::Vector{<:AbstractArray{<:Real}}
    _check_complete_param_array(prefs)
    inds = Collections.indices(prefs)
    supp_dict = _parameter_supports(first(prefs))
    if label == All || (!has_internal_supports(first(prefs)) && label == PublicLabel)
        raw_supps = collect(keys(supp_dict))
    else
        raw_supps = findall(e -> any(v -> v <: label, e), supp_dict)
    end
    return map(s -> Collections.unvectorize(s, inds), raw_supps)
end

# Dispatch for Vectors to make predictable matrix outputs
function supports(prefs::Vector{DependentParameterRef};
                  label::Type{<:AbstractSupportLabel} = PublicLabel
                  )::Array{Float64, 2}
    if !has_supports(prefs)
        return zeros(Float64, _num_parameters(first(prefs)), 0)
    elseif label == All || (!has_internal_supports(first(prefs)) && label == PublicLabel)
        raw_supps = keys(_parameter_supports(first(prefs)))
        if length(raw_supps) == 1
            return reduce(hcat, collect(raw_supps))
        else
            return reduce(hcat, raw_supps)
        end
    else
        raw_supps = findall(e -> any(v -> v <: label, e), 
                            _parameter_supports(first(prefs)))
        if isempty(raw_supps)
            return zeros(Float64, _num_parameters(first(prefs)), 0)
        else
            return reduce(hcat, raw_supps)
        end
    end
end

# Define method for overriding the current supports
function _update_parameter_supports(prefs::AbstractArray{<:DependentParameterRef},
                                    supports::Array{<:Real, 2},
                                    label::Type{<:AbstractSupportLabel})::Nothing
    domain = _parameter_domain(first(prefs))
    new_supps = Dict{Vector{Float64}, Set{DataType}}(s => Set([label]) 
                                                     for s in eachcol(supports))
    sig_figs = significant_digits(first(prefs))
    methods = _derivative_methods(first(prefs))
    new_params = DependentParameters(domain, new_supps, sig_figs, methods)
    _set_core_variable_object(first(prefs), new_params)
    _set_has_internal_supports(first(prefs), label <: InternalLabel)
    for pref in prefs 
        _reset_derivative_constraints(pref)
    end
    if any(is_used(pref) for pref in prefs)
        set_transformation_backend_ready(JuMP.owner_model(first(prefs)), false)
    end
    return
end

# Process an array of vectors into a support matrix
function _make_support_matrix(
    supports::Vector{<:AbstractArray{<:Real}},
    inds::Collections.ContainerIndices    
    )::Array{Float64, 2}
    supp_inds = Collections.indices(first(supports))
    supp_inds == inds || error("Inconsistent support indices")
    lens = [length(supp) for supp in supports]
    _allequal(lens) || error("Inconsistent support dimensions.")
    supps = Array{Float64}(undef, length(inds), length(supports))
    for i in eachindex(supports)
        @inbounds supps[:, i] = Collections.vectorize(supports[i], supp_inds)
    end
    return supps
end

"""
    set_supports(prefs::AbstractArray{<:DependentParameterRef},
                 supports::Vector{<:AbstractArray{<:Real}};
                 [force::Bool = false,
                 label::Type{<:AbstractSupportLabel} = UserDefined])::Nothing

Specify the support points for `prefs`. Errors if the supports violate the domain
of the infinite domain, if the dimensions don't match up properly,
if `prefs` and `supports` have different indices, not all of the `prefs` are
from the same dependent infinite parameter container, there are existing
supports and `force = false`. Note that it is strongly preferred to use
`add_supports` if possible to avoid destroying measure dependencies.

```julia
    set_supports(prefs::Vector{DependentParameterRef},
                 supports::Array{<:Real, 2};
                 [force::Bool = false,
                 label::Type{<:AbstractSupportLabel} = UserDefined])::Nothing
```
Specify the supports for a vector `prefs` of dependent infinite parameters.
Here rows of `supports` correspond to `prefs` and the columns correspond to the
supports. This is more efficient than the above method and will error for the
same reasons.

**Example**
```julia-repl
julia> set_supports(y, [[0, 1], [0, 1]])

julia> set_supports(x, [0 1; 0 1])

julia> supports(x)
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0
```
"""
function set_supports(
    prefs::AbstractArray{<:DependentParameterRef},
    supports::Vector{<:AbstractArray{<:Real}};
    force::Bool = false,
    label::Type{<:AbstractSupportLabel} = UserDefined
    )::Nothing
    inds = Collections.indices(prefs)
    supps = _make_support_matrix(supports, inds)
    set_supports(Collections.vectorize(prefs, inds), supps, force = force, 
                 label = label)
    return
end

# Efficient method for vector prefs and matrix of supports
function set_supports(
    prefs::Vector{DependentParameterRef},
    supports::Array{<:Real, 2};
    force::Bool = false,
    label::Type{<:AbstractSupportLabel} = UserDefined
    )::Nothing
    domain = infinite_domain(prefs) # this does a check on prefs
    if has_supports(prefs) && !force
        error("Unable set supports for $prefs since they already have supports." *
              " Consider using `add_supports` or use `force = true` to " *
              "overwrite the existing supports.")
    elseif !supports_in_domain(supports, domain)
        error("Supports violate the domain of the infinite domain.")
    end
    supports = round.(supports, sigdigits = significant_digits(first(prefs)))
    _update_parameter_supports(prefs, supports, label)
    return
end

# Error for single dependent parameters
function set_supports(pref::DependentParameterRef, supports; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    add_supports(prefs::AbstractArray{<:DependentParameterRef},
                 supports::Vector{<:AbstractArray{<:Real}};
                 [label::Type{<:AbstractSupportLabel} = UserDefined])::Nothing

Add additional support points for `prefs`. Errors if the supports violate the domain
of the infinite domain, if the dimensions don't match up properly,
if `prefs` and `supports` have different indices, or not all of the `prefs` are
from the same dependent infinite parameter container.

```julia
    add_supports(prefs::Vector{DependentParameterRef},
                 supports::Array{<:Real, 2};
                 [label::Type{<:AbstractSupportLabel} = UserDefined])::Nothing
```
Specify the supports for a vector `prefs` of dependent infinite parameters.
Here rows of `supports` correspond to `prefs` and the columns correspond to the
supports. This is more efficient than the above method and will error for the
same reasons.

**Example**
```julia-repl
julia> add_supports(x, [[1], [1]])

julia> supports(x)
2×2 Array{Float64,2}:
 0.0  1.0
 0.0  1.0

julia> add_supports(x, ones(2, 1) * 0.5)

julia> supports(t)
2×3 Array{Float64,2}:
 0.0  1.0  0.5
 0.0  1.0  0.5
```
"""
function add_supports(
    prefs::AbstractArray{<:DependentParameterRef},
    supports::Vector{<:AbstractArray{<:Real}};
    label::Type{<:AbstractSupportLabel} = UserDefined, # interal keyword args
    check::Bool = true
    )::Nothing
    inds = Collections.indices(prefs)
    supps = _make_support_matrix(supports, inds)
    add_supports(Collections.vectorize(prefs, inds), supps, label = label, 
                 check = check)
    return
end

# More efficient version for supports in the correct format
function add_supports(
    prefs::Vector{DependentParameterRef},
    supports::Array{<:Real, 2};
    label::Type{<:AbstractSupportLabel} = UserDefined, # internal keyword args
    check::Bool = true
    )::Nothing
    domain = infinite_domain(prefs) # this does a check on prefs
    if check && !supports_in_domain(supports, domain)
        error("Supports violate the domain of the infinite domain.")
    end
    supports = round.(supports, sigdigits = significant_digits(first(prefs)))
    current_supports = _parameter_supports(first(prefs))
    added_new_support = false
    for i in 1:size(supports, 2)
        s = @view(supports[:, i])
        if haskey(current_supports, s)
            push!(current_supports[s], label)
        else
            current_supports[s] = Set([label])
            added_new_support = true
        end
    end
    if label <: InternalLabel
        _set_has_internal_supports(first(prefs), true)
    end
    if added_new_support
        for pref in prefs
            _reset_derivative_constraints(pref)
        end
        if any(is_used(pref) for pref in prefs)
            set_transformation_backend_ready(JuMP.owner_model(first(prefs)), false)
        end
    end
    return
end

# Error for single dependent parameters
function add_supports(pref::DependentParameterRef, supports; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    delete_supports(prefs::AbstractArray{<:DependentParameterRef};
                    [label::Type{<:AbstractSupportLabel} = All])::Nothing

Delete the support points for `prefs`. Errors if any of the parameters are
used by a measure or if not all belong to the same set of dependent parameters.
If `label != All` then that label is removed along with any supports that solely 
contain that label.

**Example**
```julia-repl
julia> delete_supports(w)

```
"""
function delete_supports(
    prefs::AbstractArray{<:DependentParameterRef};
    label::Type{<:AbstractSupportLabel} = All
    )::Nothing
    _check_complete_param_array(prefs)
    supp_dict = _parameter_supports(first(prefs))
    for pref in prefs
        _reset_derivative_constraints(pref)
    end
    if label == All
        if any(used_by_measure(pref) for pref in prefs)
            error("Cannot delete supports with measure dependencies.")
        end
        empty!(supp_dict)
        _set_has_internal_supports(first(prefs), false)
    else 
        filter!(p -> !all(v -> v <: label, p[2]), supp_dict)
        for (k, v) in supp_dict 
            filter!(l -> !(l <: label), v)
        end
        pref1 = first(prefs)
        if has_internal_supports(pref1) && num_supports(pref1, label = InternalLabel) == 0
            _set_has_internal_supports(pref1, false)
        end
    end
    if any(is_used(pref) for pref in prefs)
        set_transformation_backend_ready(JuMP.owner_model(first(prefs)), false)
    end
    return
end

# Error for single dependent parameters
function delete_supports(pref::DependentParameterRef)
    error("Cannot delete the supports of a single dependent infinite parameter.")
end

# TODO resolve case that there are existing UniformGrid supports
"""
    generate_and_add_supports!(prefs::AbstractArray{<:DependentParameterRef},
                               domain::InfiniteArrayDomain,
                               [method::Type{<:AbstractSupportLabel}];
                               [num_supports::Int = DefaultNumSupports])::Nothing

Generate supports for `prefs` via [`generate_support_values`](@ref) and add them
to `pref`. This is intended as an extendable internal method for
[`fill_in_supports!`](@ref fill_in_supports!(::AbstractArray{<:DependentParameterRef})).
Most extensions that employ user-defined infinite domains can typically enable this
by extending [`generate_support_values`](@ref). However, in some cases it may be
necessary to extend this when more complex operations need to take place then just
adding supports to a set of infinite parameters. Errors if the
infinite domain type is not recognized.
"""
function generate_and_add_supports!(
    prefs::AbstractArray{<:DependentParameterRef},
    domain::InfiniteArrayDomain;
    num_supports::Int = DefaultNumSupports
    )::Nothing
    new_supps, label = generate_supports(domain,
                                         num_supports = num_supports,
                                    sig_digits = significant_digits(first(prefs)))
    add_supports(Collections.vectorize(prefs), new_supps, check = false, 
                 label = label)
    return
end

# Method dispatch
function generate_and_add_supports!(
    prefs::AbstractArray{<:DependentParameterRef},
    domain::InfiniteArrayDomain,
    method::Type{<:AbstractSupportLabel};
    num_supports::Int = DefaultNumSupports
    )::Nothing
    new_supps, label = generate_supports(domain, method,
                                         num_supports = num_supports,
                                    sig_digits = significant_digits(first(prefs)))
    add_supports(Collections.vectorize(prefs), new_supps, check = false, 
                 label = label)
    return
end

"""
    fill_in_supports!(prefs::AbstractArray{<:DependentParameterRef};
                      [num_supports::Int = DefaultNumSupports,
                       modify::Bool = true])::Nothing

Automatically generate support points for a container of dependent infinite
parameters `prefs`. Generating up to `num_supports` for the parameters in accordance
with `generate_and_add_supports!`. Will add nothing if there are supports and
`modify = false`. Extensions that use user defined
domain types should extend [`generate_and_add_supports!`](@ref) and/or
[`generate_support_values`](@ref) as needed. Errors if the infinite domain type is
not recognized.

**Example**
```julia-repl
julia> fill_in_supports!(x, num_supports = 4)

julia> supports(x)
2×4 Array{Float64,2}:
 0.0  0.333  0.667  1.0
 0.0  0.333  0.667  1.0
```
"""
function fill_in_supports!(
    prefs::AbstractArray{<:DependentParameterRef};
    num_supports::Int = DefaultNumSupports,
    modify::Bool = true
    )::Nothing
    domain = infinite_domain(prefs) # does check for bad container
    current_amount = InfiniteOpt.num_supports(first(prefs))
    if (modify || current_amount == 0) && (current_amount < num_supports)
        generate_and_add_supports!(prefs, domain,
                                   num_supports = num_supports - current_amount)
    end
    return
end

# Error for single dependent parameters
function fill_in_supports!(pref::DependentParameterRef; kwargs...)
    error("Cannot modify the supports of a single dependent infinite parameter.")
end

"""
    fill_in_supports!(model::InfiniteModel; [num_supports::Int = DefaultNumSupports,
                      modify::Bool = true])::Nothing

Automatically generate support points for all infinite parameters in model.
This calls `fill_in_supports!` for each parameter in the model.
See [`fill_in_supports!`](@ref)
for more information. Errors if one of the infinite domain types is unrecognized.
Note that no supports will be added to a particular parameter if it already has
some and `modify = false`.

**Example**
```julia-repl
julia> fill_in_supports!(model, num_supports = 4)

julia> supports(t)
4-element Array{Float64,1}:
 0.0
 0.333
 0.667
 1.0
```
"""
function fill_in_supports!(
    model::InfiniteModel; 
    num_supports::Int = DefaultNumSupports,
    modify::Bool = true
    )::Nothing
    # fill in the the supports of each independent parameter
    for (key, _) in model.independent_params
        pref = dispatch_variable_ref(model, key)
        fill_in_supports!(pref, num_supports = num_supports, modify = modify)
    end
    # fill in the supports of each dependent parameter domain
    for (key, data_object) in model.dependent_params
        prefs = [dispatch_variable_ref(model, DependentParameterIndex(key, i))
                 for i in 1:length(data_object.names)]
        fill_in_supports!(prefs, num_supports = num_supports, modify = modify)
    end
    return
end

################################################################################
#                          MODEL PARAMETER QUERIES
################################################################################
"""
    num_parameters(model::InfiniteModel,
                   [type::Type{InfOptParameter} = InfOptParameter])::Int

Return the number of `InfiniteOpt` parameters assigned to `model`. By default,
the total number of infinite and finite parameters is returned. The amount
of a particular type is obtained by specifying the concrete parameter type
of [`InfOptParameter`](@ref) via `type`. Type options include:
 - `InfOptParameter`: all parameters
 - `ScalarParameter`: all scalar parameters
 - `InfiniteParameter`: all infinite parameters
 - `FiniteParameter`: all finite parameters
 - `IndependentParameter`: all independent infinite parameters
 - `DependentParameters`: all dependent infinite parameters

**Example**
```julia-repl
julia> num_parameters(model)
3

julia> num_parameters(model, IndependentParameter)
2
```
"""
function num_parameters(
    model::InfiniteModel,
    type::Type{InfOptParameter} = InfOptParameter
    )::Int
    num_pars = num_parameters(model, IndependentParameter)
    num_pars += num_parameters(model, FiniteParameter)
    num_pars += num_parameters(model, DependentParameters)
    return num_pars
end

# Particular scalar parameter types
function num_parameters(
    model::InfiniteModel,
    type::Type{C}
    )::Int where {C <: ScalarParameter}
    return length(_data_dictionary(model, type))
end

# ScalarParameter
function num_parameters(
    model::InfiniteModel,
    type::Type{ScalarParameter}
    )::Int
    num_pars = num_parameters(model, FiniteParameter)
    num_pars += num_parameters(model, IndependentParameter)
    return num_pars
end

# DependentParameters
function num_parameters(
    model::InfiniteModel,
    type::Type{DependentParameters}
    )::Int
    num_pars = 0
    for (_, object) in _data_dictionary(model, type)
        num_pars += length(object.names)
    end
    return num_pars
end

# InfiniteParameter
function num_parameters(
    model::InfiniteModel,
    type::Type{InfiniteParameter}
    )::Int
    num_pars = num_parameters(model, IndependentParameter)
    num_pars += num_parameters(model, DependentParameters)
    return num_pars
end

"""
    all_parameters(model::InfiniteModel,
                   type::Type{InfOptParameter} = InfOptParameter
                   )::Vector{GeneralVariableRef}

Return a list of all the `InfiniteOpt` parameters assigned to `model`. By default,
all of the infinite and finite parameters is returned. The search is reduced to
a particular type is obtained by specifying the concrete parameter type
of [`InfOptParameter`](@ref) via `type`. Type options include:
- `InfOptParameter`: all parameters
- `ScalarParameter`: all scalar parameters
- `InfiniteParameter`: all infinite parameters
- `FiniteParameter`: all finite parameters
- `IndependentParameter`: all independent infinite parameters
- `DependentParameters`: all dependent infinite parameters

**Examples**
```julia-repl
julia> all_parameters(model)
4-element Array{GeneralVariableRef,1}:
 t
 x[1]
 x[2]
 alpha

julia> all_parameters(model, FiniteParameter)
1-element Array{GeneralVariableRef,1}:
 alpha
```
"""
function all_parameters(
    model::InfiniteModel,
    type::Type{InfOptParameter} = InfOptParameter
    )::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, IndependentParameter)
    append!(prefs_list, all_parameters(model, DependentParameters))
    append!(prefs_list, all_parameters(model, FiniteParameter))
    return prefs_list
end

# Particular scalar parameter types
function all_parameters(
    model::InfiniteModel,
    type::Type{C}
    )::Vector{GeneralVariableRef} where {C <: InfOptParameter}
    prefs_list = Vector{GeneralVariableRef}(undef, num_parameters(model, type))
    for (i, (index, _)) in enumerate(_data_dictionary(model, type))
        prefs_list[i] = _make_parameter_ref(model, index)
    end
    return prefs_list
end

# ScalarParameter
function all_parameters(
    model::InfiniteModel,
    type::Type{ScalarParameter}
    )::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, IndependentParameter)
    append!(prefs_list, all_parameters(model, FiniteParameter))
    return prefs_list
end

# DependentParameters
function all_parameters(
    model::InfiniteModel,
    type::Type{DependentParameters}
    )::Vector{GeneralVariableRef}
    prefs_list = Vector{GeneralVariableRef}(undef, num_parameters(model, type))
    counter = 1
    for (index, object) in _data_dictionary(model, type)
        for i in eachindex(object.names)
            dep_idx = DependentParameterIndex(index, i)
            prefs_list[counter] = _make_parameter_ref(model, dep_idx)
            counter += 1
        end
    end
    return prefs_list
end

# InfiniteParameter
function all_parameters(
    model::InfiniteModel,
    type::Type{InfiniteParameter}
    )::Vector{GeneralVariableRef}
    prefs_list = all_parameters(model, DependentParameters)
    append!(prefs_list, all_parameters(model, IndependentParameter))
    return prefs_list
end

################################################################################
#                                 DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel,
                prefs::AbstractArray{<:DependentParameterRef})::Nothing

Extend `JuMP.delete` to delete
dependent infinite parameters and their dependencies. All variables, constraints, and
measure functions that depend on `prefs` are updated to exclude them. Errors if the
parameters are contained in an `AbstractMeasureData` datatype that is employed by
a measure since the measure becomes invalid otherwise. Thus, measures that
contain this dependency must be deleted first. Note that
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) needs to be
extended to allow deletion of parameters when custom `AbstractMeasureData`
datatypes are used. Note that any dependent infinite variables will have their
start values reset to the default via [`reset_start_value_function`](@ref).

**Example**
```julia-repl
julia> print(model)
Min measure(g(t, x)*t + x) + z
Subject to
 z ≥ 0.0
 g(t, x) + z ≥ 42.0, ∀ t ∈ [0, 6], x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
 g(0.5, x) = 0, x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]

julia> delete(model, x)

julia> print(model)
Min measure(g(t)*t) + z
Subject to
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]
 g(0.5) = 0
```
"""
function JuMP.delete(
    model::InfiniteModel,
    prefs::AbstractArray{<:DependentParameterRef}
    )::Nothing
    @assert JuMP.is_valid(model, first(prefs)) "Parameter references are invalid."
    _check_complete_param_array(prefs)
    gvrefs = [_make_parameter_ref(model, JuMP.index(pref)) for pref in prefs]
    # ensure deletion is okay (prefs are not used by measure data)
    for pref in gvrefs
        for mindex in _measure_dependencies(pref)
            data = measure_data(dispatch_variable_ref(model, mindex))
            _check_param_in_data(pref, data)
        end
    end
    # make sure it isn't used by infinite variable
    if used_by_infinite_variable(first(prefs))
        error("Cannot delete `$prefs` since they are used by an infinite " * 
              "variable(s).")
    end
    # make sure it isn't used by parameter function 
    if used_by_parameter_function(first(prefs))
        error("Cannot delete `$prefs` since they are used by an infinite " * 
              "parameter function(s).")
    end
    # update transformation backend status
    if any(is_used(pref) for pref in prefs)
        set_transformation_backend_ready(model, false)
    end
    # delete dependence of measures and constraints on prefs
    for pref in gvrefs
        _update_measures(model, pref)
        _update_constraints(model, pref)
    end
    # get the object and parameter numbers
    group_int_idx = parameter_group_int_index(first(prefs))
    param_nums = collect(_data_object(first(prefs)).parameter_nums)
    # delete derivatives that depend on any of these parameters 
    for pref in gvrefs 
        for index in _derivative_dependencies(pref)
            JuMP.delete(model, dispatch_variable_ref(model, index))
        end
    end
    # delete parameter information stored in model
    _delete_data_object(first(prefs))
    # update the parameter group integer indices and parameter numbers
    _update_model_numbers(model, group_int_idx, param_nums)
    return
end
