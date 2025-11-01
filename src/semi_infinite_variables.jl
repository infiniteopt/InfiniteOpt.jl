################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::SemiInfiniteVariableIndex
    )
    return SemiInfiniteVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:SemiInfiniteVariable}
    )
    return MOIUC.add_item(model.semi_infinite_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel, 
    ::Type{SemiInfiniteVariable}
    )
    return model.semi_infinite_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(vref::SemiInfiniteVariableRef)
    return JuMP.owner_model(vref).semi_infinite_vars
end

# Extend _data_object
function _data_object(vref::SemiInfiniteVariableRef)
    object = get(_data_dictionary(vref), JuMP.index(vref), nothing)
    if isnothing(object) 
        error("Invalid point variable reference, cannot find ",
              "corresponding variable in the model. This is likely ",
              "caused by using the reference of a deleted variable.")
    end
    return object
end

"""
    internal_semi_infinite_variable(
        vref::SemiInfiniteVariableRef,
        backend::AbstractTransformationBackend
        )::SemiInfiniteVariable

Return the semi-infinite variable object of `vref` assuming it is an internal variable
made during measure expansion within a transformation backend. This will apply to
transformation backend extensions that utilize `add_measure_variable` in combination
with `expand_measure`.
"""
function internal_semi_infinite_variable end

"""
    core_object(vref::SemiInfiniteVariableRef)::SemiInfiniteVariable

Retrieve the underlying core [`SemiInfiniteVariable`] object for `vref`. 
This is intended as an advanced method for developers.
"""
function core_object(vref::SemiInfiniteVariableRef)
    if !haskey(_data_dictionary(vref), JuMP.index(vref))
        model = JuMP.owner_model(vref)
        return internal_semi_infinite_variable(vref, model.backend)
    else
        return _data_object(vref).variable
    end
end

"""
    parameter_group_int_indices(vref::SemiInfiniteVariableRef)::Vector{Int}

Return the list of infinite parameter group integer indices used by `vref`.
"""
function parameter_group_int_indices(vref::SemiInfiniteVariableRef)
    return core_object(vref).group_int_idxs
end

# Extend _parameter_numbers
function _parameter_numbers(vref::SemiInfiniteVariableRef)
    return core_object(vref).parameter_nums
end

################################################################################
#                             DEFINITION METHODS
################################################################################
## Process the value arguments to be formatted correctly 
## (avoid AffExpr array promotions)
# Not an AffExpr
function _process_value(v)
    return v
end

# An AffExpr
function _process_value(v::JuMP.GenericAffExpr)
    if isempty(v.terms)
        return v.constant
    elseif length(v.terms) == 1 && isone(collect(values(v.terms))[1]) && 
           iszero(v.constant)
        return first(keys(v.terms))
    else
        return v
    end
end

# Convenient dispatch for restrict syntax
function JuMP.build_variable(
    _error::Function,  
    ivref::GeneralVariableRef,
    raw_params::Collections.VectorTuple,
    restricted_info::RestrictedDomainInfo
    )
    # check that the values are the same format as the infinite variable 
    ivref_prefs = raw_parameter_refs(ivref)
    if !all(v isa Union{Real, GeneralVariableRef} for v in raw_params)
        _error("Unexpected inputs given, expected them to be parameter references ",
               "and real numbers.")
    elseif !Collections.same_structure(raw_params, ivref_prefs)
        _error("The parameter reference/value input format $(raw_params) does ",
               "not match that of the infinite variable $(ivref).")
    end
    # make the evaluation supports
    eval_support = Float64[p isa Real ? p : NaN for p in raw_params]
    # check that no dependent parameters are partially transcribed
    for r in ivref_prefs.ranges
        real_inds = map(s -> !isnan(s), @view(eval_support[r]))
        if any(real_inds) != all(real_inds)
            _error("Cannot partially evaluate a multi-dimensional input of an " *
                   "infinite variable.")
        end
    end
    # dispatch to the main builder 
    return JuMP.build_variable(_error, ivref, eval_support, restricted_info)
end

"""
    SemiInfinite{T} <: InfOptVariableType 

A `DataType` to assist in making semi-infinite variables. This can be passed as an 
extra argument to `@variable` to make such a variable: 
```julia 
@variable(model, var_expr, SemiInfinite(inf_var, parameter_values...), kwargs...)
```
Here `parameter_values` must match the format of the infinite parameter 
references associated with the infinite variable `inf_var` and can be comprised 
of both real valued supports and/or infinite parameters.

**Fields**
- `infinite_variable_ref::GeneralVariableRef`
- `parameter_values::VectorTuple{T}`: The infinite parameters and/or infinite 
   parameter support values the variable will depend on.
"""
struct SemiInfinite{T} <: InfOptVariableType 
    infinite_variable_ref::GeneralVariableRef
    parameter_values::Collections.VectorTuple{T}
    function SemiInfinite(
        ivref::GeneralVariableRef,
        vt::Collections.VectorTuple{T}
        ) where {T}
        processed_vals = _process_value.(vt.values)
        if !isequal(processed_vals, vt.values)
            vt2 = Collections.VectorTuple(processed_vals, vt.ranges, vt.dimensions, vt.num_columns)
            return new{only(typeof(vt2).parameters)}(ivref, vt2)
        end
        return new{T}(ivref, vt)
    end
end
function SemiInfinite(ivref::GeneralVariableRef, vals...)
    vt = Collections.VectorTuple(vals)
    return SemiInfinite(ivref, vt)
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, 
                        var_type::SemiInfinite)::SemiInfiniteVariable{GeneralVariableRef}

Build and return a semi-infinite variable based on `info` and `var_type`. Errors 
if the information stored in `var_type` is invalid. See [`SemiInfinite`](@ref) 
for more information.

**Example**
```julia-repl
julia> y
y(t, x)

julia> info = VariableInfo(false, 0, false, 0, false, 0, true, 0, false, false);

julia> semi_inf_var = build_variable(error, info, SemiInfinite(y, t, 0));
```
"""
function JuMP.build_variable(
    _error::Function, 
    info::JuMP.VariableInfo, 
    var_type::SemiInfinite; 
    extra_kwargs...
    )
    if !isempty(extra_kwargs)
        _error("Keyword argument $(first(keys(extra_kwargs))) is not for use with semi-infinite ",
               "variables.")
    end
    restricted_info = _process_restricted_info(_error, info)
    ivref = var_type.infinite_variable_ref
    raw_params = var_type.parameter_values
    return JuMP.build_variable(_error, ivref, raw_params, restricted_info)
end

"""
    JuMP.build_variable(
        _error::Function, 
        ivref::GeneralVariableRef,
        eval_support::Vector{Float64},
        [restricted_info::RestrictedDomainInfo]; 
        [check::Bool = true]
        )::SemiInfiniteVariable{GeneralVariableRef}

Extend the `JuMP.build_variable` function to build a semi-infinite variable
based on the infinite variable/derivative/parameter function `ivref` with 
reduction support `eval_support`. Will check that input is appropriate if 
`check = true`. Errors if `ivref` is not an infinite variable, `eval_support` 
violate infinite parameter domains, or if the support dimensions don't match the 
infinite parameter dimensions of `ivref`. This is intended an internal method for 
use in evaluating measures.
"""
function JuMP.build_variable(
    _error::Function, 
    ivref::GeneralVariableRef,
    eval_support::Vector{Float64},
    restricted_info::RestrictedDomainInfo;
    check::Bool = true
    )
    # check the inputs
    dvref = dispatch_variable_ref(ivref)
    if check && !(dvref isa Union{InfiniteVariableRef, DerivativeRef, ParameterFunctionRef})
         _error("Must specify an infinite variable/derivative/parameter function dependency.")
    elseif check && length(eval_support) != length(parameter_list(dvref))
        _error("Support evaluation dictionary indices do not match the infinite " *
               "parameter dependencies of $(ivref).")
    end
    prefs = parameter_list(dvref)
    if check
        for (pref, value) in zip(prefs, eval_support)
            if !isnan(value) && JuMP.has_lower_bound(pref) && !supports_in_domain(value, infinite_domain(pref))
                _error("Evaluation support violates infinite parameter domain(s).")
            end
        end
    end
    # get the parameter group integer indices of the dependencies
    raw_prefs = raw_parameter_refs(dvref)
    group_int_idxs = [
        idx
        for (i, idx) in enumerate(parameter_group_int_indices(dvref))
        if isnan(eval_support[raw_prefs.ranges[i].start])
    ]
    # get the parameter numbers
    orig_nums = _parameter_numbers(ivref)
    param_nums = [orig_nums[i] for i in eachindex(orig_nums)
                  if isnan(eval_support[i])]
    # round the support values in accordance with the significant digits
    for i in eachindex(eval_support)
        eval_support[i] = round(
            eval_support[i], 
            sigdigits = significant_digits(prefs[i])
        )
    end
    # build the variable
    return SemiInfiniteVariable(
        restricted_info,
        ivref,
        eval_support,
        param_nums, 
        group_int_idxs
    )
end
function JuMP.build_variable(
    _error::Function, 
    ivref::GeneralVariableRef,
    eval_support::Vector{Float64};
    check::Bool = true
    )
    return JuMP.build_variable(
        _error,
        ivref,
        eval_support,
        RestrictedDomainInfo(),
        check = check
    )
end

## Make dispatch functions to add semi-infinite variable based supports
# Independent parameter
function _add_semi_infinite_support(
    prefs::Vector{IndependentParameterRef}, 
    eval_supp::Vector{Float64}
    )
    add_supports(
        only(prefs),
        only(eval_supp), 
        check = false, 
        label = UserDefined
    )
    return
end

# Dependent parameters
function _add_semi_infinite_support(
    prefs::Vector{DependentParameterRef}, 
    eval_supp::Vector{Float64}
    )
    supp = reshape(eval_supp, length(prefs), 1)
    add_supports(prefs, supp, check = false, label = UserDefined)
    return
end

"""
    JuMP.add_variable(
        model::InfiniteModel,
        var::SemiInfiniteVariable,
        [name::String = ""]
    )::GeneralVariableRef

Extend the `JuMP.add_variable` function to accomodate `InfiniteOpt` 
semi-infinite variable types. Adds `var` to the infinite model `model` and 
returns a [`GeneralVariableRef`](@ref). Primarily intended to be an internal 
function used in evaluating measures.
"""
function JuMP.add_variable(
    model::InfiniteModel, 
    var::SemiInfiniteVariable,
    name::String = "";
    add_support = true
    )
    ivref = var.infinite_variable_ref
    divref = dispatch_variable_ref(ivref)
    eval_supp = var.eval_support
    info = var.info
    JuMP.check_belongs_to_model(divref, model)
    existing_index = get(model.semi_lookup, (ivref, eval_supp), nothing)
    is_active_info = info.active_lower_bound_info || info.active_upper_bound_info || 
                    info.active_fix_info || info.active_start_info
    if isnothing(existing_index)
        data_object = VariableData(var, name)
        vindex = _add_data_object(model, data_object)
        push!(_semi_infinite_variable_dependencies(divref), vindex)
        if add_support
            prefs = raw_parameter_refs(divref)
            for r in prefs.ranges
                if !isnan(eval_supp[r.start])
                    dprefs = dispatch_variable_ref.(prefs[r])
                    _add_semi_infinite_support(dprefs, eval_supp[r])
                end
            end
        end
        model.semi_lookup[(ivref, eval_supp)] = vindex
        vref = SemiInfiniteVariableRef(model, vindex)
    else
        vindex = existing_index
        if !isempty(name)
            model.semi_infinite_vars[vindex].name = name
        end
        vref = SemiInfiniteVariableRef(model, vindex)
        if is_active_info
            _delete_info_constraints(vref)
            _set_core_object(vref, var)
        end
    end
    gvref = GeneralVariableRef(vref)
    if is_active_info
        _set_info_constraints(info, gvref, vref)
    end
    return gvref
end

################################################################################
#                          PARAMETER REFERENCE METHODS
################################################################################
"""
    infinite_variable_ref(vref::SemiInfiniteVariableRef)::GeneralVariableRef

Return the infinite variable/derivative/parameter function reference associated 
with the semi-infinite variable `vref`.

**Example**
```julia-repl
julia> infinite_variable_ref(vref)
g(t, x)
```
"""
function infinite_variable_ref(vref::SemiInfiniteVariableRef)
    return core_object(vref).infinite_variable_ref
end

"""
    eval_support(vref::SemiInfiniteVariableRef)::Vector{Float64}

Return the evaluation supports associated with the semi-infinite variable
`vref`. Any element corresponding to an infinite parameter that is not
evaluated is filled with a `NaN`.

**Example**
```julia-repl
julia> supp = eval_support(vref)
[0.0, NaN]
```
"""
function eval_support(vref::SemiInfiniteVariableRef)
    return core_object(vref).eval_support
end

# helper version of raw_parameter_refs
function raw_parameter_refs(var::SemiInfiniteVariable)
    orig_prefs = raw_parameter_refs(var.infinite_variable_ref)
    delete_indices = [
        isnan(var.eval_support[orig_prefs.ranges[i].start])
        for i in 1:size(orig_prefs, 1)
    ]
    return Collections.restricted_copy(orig_prefs, delete_indices)
end

"""
    raw_parameter_refs(vref::Union{SemiInfiniteVariableRef, SemiInfiniteVariable})::VectorTuple
Return the raw [`VectorTuple`](@ref InfiniteOpt.Collections.VectorTuple) of the 
parameter references that `vref` depends on. This is primarily an internal method 
where [`parameter_refs`](@ref parameter_refs(vref::SemiInfiniteVariableRef)) 
is intended as the preferred user function.
"""
function raw_parameter_refs(vref::SemiInfiniteVariableRef)
    return raw_parameter_refs(core_object(vref))
end

"""
    parameter_refs(vref::SemiInfiniteVariableRef)::Tuple

Return the infinite parameter references associated with the semi-infinite variable
`vref`. This is formatted as a `Tuple` of containing the parameter references as
they were inputted to define the untranscripted infinite variable except, the
evaluated parameters are excluded.

**Example**
```julia-repl
julia> parameter_refs(vref)
(t, [x[1], x[2]])
```
"""
function parameter_refs(vref::SemiInfiniteVariableRef)
    return Tuple(raw_parameter_refs(vref))
end

"""
    parameter_list(vref::SemiInfiniteVariableRef)::Vector{GeneralVariableRef}

Return a vector of the parameter references that `vref` depends on. This is
primarily an internal method where [`parameter_refs`](@ref parameter_refs(vref::SemiInfiniteVariableRef))
is intended as the preferred user function.
"""
function parameter_list(vref::SemiInfiniteVariableRef)
    orig_prefs = parameter_list(infinite_variable_ref(vref))
    eval_supp = eval_support(vref)
    return [orig_prefs[i] for i in eachindex(orig_prefs) if isnan(eval_supp[i])]
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _derivative_dependencies
function _derivative_dependencies(vref::SemiInfiniteVariableRef)
    return _data_object(vref).derivative_indices
end

"""
    used_by_derivative(vref::SemiInfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a derivative.

**Example**
```julia-repl
julia> used_by_derivative(vref)
true
```
"""
function used_by_derivative(vref::SemiInfiniteVariableRef)
    return !isempty(_derivative_dependencies(vref))
end

"""
    is_used(vref::SemiInfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```julia-repl
julia> is_used(vref)
true
```
"""
function is_used(vref::SemiInfiniteVariableRef)
    if used_by_measure(vref) || used_by_constraint(vref)
        return true
    end
    for dindex in _derivative_dependencies(vref)
        if is_used(DerivativeRef(JuMP.owner_model(vref), dindex))
            return true
        end
    end
    return false
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Set info for point variables
function _update_variable_info(
    vref::SemiInfiniteVariableRef,
    info::RestrictedDomainInfo
    )
    var = core_object(vref)
    new_var = SemiInfiniteVariable(
        info,
        var.infinite_variable_ref,
        var.eval_support,
        var.parameter_nums,
        var.group_int_idxs
        )
    _set_core_object(vref, new_var)
    return
end

################################################################################
#                                  DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::SemiInfiniteVariableRef)
    # remove mapping to infinite variable
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref), _semi_infinite_variable_dependencies(ivref))
    # delete associated derivative variables and mapping 
    model = JuMP.owner_model(vref)
    for index in copy(_derivative_dependencies(vref))
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # remove the lookup entry
    delete!(JuMP.owner_model(vref).semi_lookup, (ivref, eval_support(vref)))
    return
end
