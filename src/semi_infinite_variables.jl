################################################################################
#                   CORE DISPATCHVARIABLEREF METHOD EXTENSIONS
################################################################################
# Extend dispatch_variable_ref
function dispatch_variable_ref(
    model::InfiniteModel,
    index::SemiInfiniteVariableIndex
    )::SemiInfiniteVariableRef
    return SemiInfiniteVariableRef(model, index)
end

# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::VariableData{<:SemiInfiniteVariable}
    )::SemiInfiniteVariableIndex
    return MOIUC.add_item(model.semi_infinite_vars, object)
end

# Extend _data_dictionary (type based)
function _data_dictionary(
    model::InfiniteModel, 
    ::Type{SemiInfiniteVariable}
    )::MOIUC.CleverDict{SemiInfiniteVariableIndex, VariableData{SemiInfiniteVariable{GeneralVariableRef}}}
    return model.semi_infinite_vars
end

# Extend _data_dictionary (reference based)
function _data_dictionary(
    vref::SemiInfiniteVariableRef
    )::MOIUC.CleverDict{SemiInfiniteVariableIndex, VariableData{SemiInfiniteVariable{GeneralVariableRef}}}
    return JuMP.owner_model(vref).semi_infinite_vars
end

# Extend _data_object
function _data_object(
    vref::SemiInfiniteVariableRef
    )::VariableData{SemiInfiniteVariable{GeneralVariableRef}}
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

"""
    SemiInfinite{V, VT <: VectorTuple} <: InfOptVariableType 

A `DataType` to assist in making semi-infinite variables. This can be passed as an 
extra argument to `@variable` to make such a variable: 
```julia 
@variable(model, var_expr, SemiInfinite(inf_var, parameter_values...), kwargs...)
```
Here `parameter_values` must match the format of the infinite parameter 
references associated with the infinite variable `inf_var` and can be comprised 
of both real valued supports and/or infinite parameters.

**Fields**
- `infinite_variable_ref::V`
- `parameter_values::VT`: The infinite parameters and/or infinite 
   parameter support values the variable will depend on.
"""
struct SemiInfinite{V, VT <: Collections.VectorTuple} <: InfOptVariableType 
    infinite_variable_ref::V 
    parameter_values::VT
    function SemiInfinite(ivref::V, vt::VT) where {V, VT <: Collections.VectorTuple}
        processed_vals = _process_value.(vt.values)
        if !isequal(processed_vals, vt.values)
            vt2 = Collections.VectorTuple(processed_vals, vt.ranges, vt.indices)
            return new{V, typeof(vt2)}(ivref, vt2)
        end
        return new{V, VT}(ivref, vt)
    end
end
function SemiInfinite(ivref, vals...)
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
    )::SemiInfiniteVariable{GeneralVariableRef}
    # check for unneeded keywords
    for (kwarg, _) in extra_kwargs
        _error("Keyword argument $kwarg is not for use with semi-infinite ",
               "variables.")
    end
    # ensure that all the info was left alone (TODO relax this)
    if info.has_lb || info.has_ub || info.has_fix || info.has_start || 
       info.integer || info.binary
        _error("Semi-infinite variables do not currently support the ",
               "specification of variable information fields (e.g., lower ",
               "bounds, upper bounds, start values, etc). These are simply ",
               "inherited from the corresponding infinite variable.")
    end
    # check the infinite variable reference
    ivref = var_type.infinite_variable_ref
    if !(ivref isa GeneralVariableRef)
        _error("Expected an infinite variable/derivative reference dependency ",
               "of type `GeneralVariableRef`, but got an argument of type ",
               "$(typeof(ivref)).")
    end
    # check that the values are the same format as the infinite variable 
    raw_params = var_type.parameter_values
    if !all(v isa Union{Real, GeneralVariableRef} for v in raw_params)
        _error("Unexpected inputs given, expected them to be parameter references ",
               "and real numbers.")
    elseif !Collections.same_structure(raw_params, raw_parameter_refs(ivref))
        _error("The parameter reference/value input format $(raw_params) does ",
               "not match that of the infinite variable $(ivref).")
    end
    # make the evaluation support dictionary 
    inds = findall(p -> p isa Real, raw_params)
    eval_supports = Dict{Int, Float64}(i => raw_params[i] for i in inds)
    # dispatch to the main builder 
    return JuMP.build_variable(_error, ivref, eval_supports)
end

"""
    JuMP.build_variable(_error::Function, ivref::GeneralVariableRef,
                        eval_supports::Dict{Int, Float64}; [check::Bool = true]
                        )::SemiInfiniteVariable{GeneralVariableRef}

Extend the `JuMP.build_variable` function to build a semi-infinite variable
based on the infinite variable/derivative/parameter function `ivref` with 
reduction support `eval_supports`. Will check that input is appropriate if 
`check = true`. Errors if `ivref` is not an infinite variable, `eval_supports` 
violate infinite parameter domains, or if the support dimensions don't match the 
infinite parameter dimensions of `ivref`. This is intended an internal method for 
use in evaluating measures.
"""
function JuMP.build_variable(
    _error::Function, 
    ivref::GeneralVariableRef,
    eval_supports::Dict{Int, Float64};
    check::Bool = true
    )::SemiInfiniteVariable{GeneralVariableRef}
    # check the inputs
    dvref = dispatch_variable_ref(ivref)
    if check && !(dvref isa Union{InfiniteVariableRef, DerivativeRef, ParameterFunctionRef})
         _error("Must specify an infinite variable/derivative/parameter function dependency.")
    elseif check && maximum(keys(eval_supports)) > length(parameter_list(dvref))
        _error("Support evaluation dictionary indices do not match the infinite " *
               "parameter dependencies of $(ivref).")
    end
    prefs = parameter_list(dvref)
    if check
        for (index, value) in eval_supports
            pref = prefs[index]
            if JuMP.has_lower_bound(pref) && !supports_in_domain(value, infinite_domain(pref))
                _error("Evaluation support violates infinite parameter domain(s).")
            end
        end
    end
    # get the parameter group integer indices of the dependencies
    group_int_idxs = Int[]
    for i in eachindex(prefs)
        if !haskey(eval_supports, i)
            union!(group_int_idxs, parameter_group_int_index(prefs[i]))
        end
    end
    # get the parameter numbers
    orig_nums = _parameter_numbers(ivref)
    param_nums = [orig_nums[i] for i in eachindex(orig_nums)
                  if !haskey(eval_supports, i)]
    # round the support values in accordance with the significant digits
    for (k, v) in eval_supports
        eval_supports[k] = round(v, sigdigits = significant_digits(prefs[k]))
    end
    # build the variable
    return SemiInfiniteVariable(ivref, eval_supports, param_nums, group_int_idxs)
end

## Make dispatch functions to add semi-infinite variable based supports
# Independent parameters
function _add_semi_infinite_support(
    prefs::Vector{IndependentParameterRef}, 
    eval_supps::Dict{Int, Float64}, 
    counter::Int
    )::Int
    for pref in prefs
        if haskey(eval_supps, counter)
            add_supports(pref, eval_supps[counter], check = false, 
                        label = UserDefined)
        end
        counter += 1
    end
    return counter
end

# Dependent parameters
function _add_semi_infinite_support(
    prefs::Vector{DependentParameterRef}, 
    eval_supps::Dict{Int, Float64}, 
    counter::Int
    )::Int
    len = length(prefs)
    if all(haskey(eval_supps, i) for i in counter:counter+len-1)
        supp = Matrix{Float64}(undef, len, 1)
        for i in eachindex(supp)
            supp[i] = eval_supps[counter + i - 1]
        end
        add_supports(prefs, supp, check = false, label = UserDefined)
    end
    return counter + len
end

# Update the parameter supports as needed
function _update_param_supports(
    ivref::Union{InfiniteVariableRef, DerivativeRef}, 
    eval_supps::Dict{Int, Float64}
    )::Nothing
    raw_prefs = raw_parameter_refs(ivref)
    counter = 1
    for i in 1:size(raw_prefs, 1)
        prefs = dispatch_variable_ref.(raw_prefs[i, :])
        counter = _add_semi_infinite_support(prefs, eval_supps, counter)
    end
end

# Update the parameter supports as needed
function _update_param_supports(
    ivref::ParameterFunctionRef, 
    eval_supps::Dict{Int, Float64}
    )::Nothing
    return
end

"""
    JuMP.add_variable(model::InfiniteModel, var::SemiInfiniteVariable,
                      [name::String = ""])::GeneralVariableRef

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
    )::GeneralVariableRef
    ivref = var.infinite_variable_ref
    divref = dispatch_variable_ref(ivref)
    eval_supps = var.eval_supports
    JuMP.check_belongs_to_model(divref, model)
    existing_index = get(model.semi_lookup, (ivref, eval_supps), nothing)
    if isnothing(existing_index)
        data_object = VariableData(var, name)
        vindex = _add_data_object(model, data_object)
        push!(_semi_infinite_variable_dependencies(divref), vindex)
        if add_support 
            _update_param_supports(divref, eval_supps)
        end
        model.semi_lookup[(ivref, eval_supps)] = vindex
    else
        vindex = existing_index
        if !isempty(name)
            model.semi_infinite_vars[vindex].name = name
        end
    end
    return GeneralVariableRef(model, vindex)
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
function infinite_variable_ref(vref::SemiInfiniteVariableRef)::GeneralVariableRef
    return core_object(vref).infinite_variable_ref
end

"""
    eval_supports(vref::SemiInfiniteVariableRef)::Dict{Int, Float64}

Return the evaluation supports associated with the semi-infinite variable
`vref`.

**Example**
```julia-repl
julia> eval_supports(vref)
Dict{Int64,Float64} with 1 entry:
  1 => 0.5
```
"""
function eval_supports(vref::SemiInfiniteVariableRef)::Dict{Int, Float64}
    return core_object(vref).eval_supports
end

# helper version of raw_parameter_refs
function raw_parameter_refs(var::SemiInfiniteVariable)
    orig_prefs = raw_parameter_refs(var.infinite_variable_ref)
    eval_supps = var.eval_supports
    delete_indices = [!haskey(eval_supps, i) for i in eachindex(orig_prefs)]
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
function parameter_list(vref::SemiInfiniteVariableRef)::Vector{GeneralVariableRef}
    orig_prefs = parameter_list(infinite_variable_ref(vref))
    eval_supps = eval_supports(vref)
    return [orig_prefs[i] for i in eachindex(orig_prefs) if !haskey(eval_supps, i)]
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _derivative_dependencies
function _derivative_dependencies(
    vref::SemiInfiniteVariableRef
    )::Vector{DerivativeIndex}
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
function used_by_derivative(vref::SemiInfiniteVariableRef)::Bool
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
function is_used(vref::SemiInfiniteVariableRef)::Bool
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
#                            VARIABLE INFO METHODS
################################################################################
"""
    JuMP.has_lower_bound(vref::SemiInfiniteVariableRef)::Bool

Extend `JuMP.has_lower_bound` to return a `Bool` whether the original
infinite variable of `vref` has a lower bound.

**Example**
```julia-repl
julia> has_lower_bound(vref)
true
```
"""
function JuMP.has_lower_bound(vref::SemiInfiniteVariableRef)::Bool
    return JuMP.has_lower_bound(infinite_variable_ref(vref))
end

"""
    JuMP.lower_bound(vref::SemiInfiniteVariableRef)::Float64

Extend `JuMP.lower_bound` to return the lower bound of the original
infinite variable of `vref`. Errors if `vref` doesn't have a lower bound.

**Example**
```julia-repl
julia> lower_bound(vref)
0.0
```
"""
function JuMP.lower_bound(vref::SemiInfiniteVariableRef)::Float64
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.has_lower_bound(ivref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return JuMP.lower_bound(ivref)
end

# Extend to return the index of the lower bound constraint associated with the
# original infinite variable of `vref`.
function InfiniteOpt._lower_bound_index(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintIndex
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.has_lower_bound(ivref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return InfiniteOpt._lower_bound_index(ivref)
end

"""
    JuMP.LowerBoundRef(vref::SemiInfiniteVariableRef)::InfOptConstraintRef

Extend `JuMP.LowerBoundRef` to extract a constraint reference for the
lower bound of the original infinite variable of `vref`.

**Example**
```julia-repl
julia> cref = LowerBoundRef(vref)
var >= 0.0
```
"""
function JuMP.LowerBoundRef(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintRef
    return JuMP.LowerBoundRef(infinite_variable_ref(vref))
end

"""
    JuMP.has_upper_bound(vref::SemiInfiniteVariableRef)::Bool

Extend `JuMP.has_upper_bound` to return a `Bool` whether the original
infinite variable of `vref` has an upper bound.

**Example**
```julia-repl
julia> has_upper_bound(vref)
true
```
"""
function JuMP.has_upper_bound(vref::SemiInfiniteVariableRef)::Bool
    return JuMP.has_upper_bound(infinite_variable_ref(vref))
end

"""
    JuMP.upper_bound(vref::SemiInfiniteVariableRef)::Float64

Extend `JuMP.upper_bound` to return the upper bound of the original
infinite variable of `vref`. Errors if `vref` doesn't have a upper bound.

**Example**
```julia-repl
julia> upper_bound(vref)
0.0
```
"""
function JuMP.upper_bound(vref::SemiInfiniteVariableRef)::Float64
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.has_upper_bound(ivref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return JuMP.upper_bound(ivref)
end

# Extend to return the index of the upper bound constraint associated with the
# original infinite variable of `vref`.
function InfiniteOpt._upper_bound_index(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintIndex
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.has_upper_bound(ivref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return InfiniteOpt._upper_bound_index(ivref)
end

"""
    JuMP.UpperBoundRef(vref::SemiInfiniteVariableRef)::InfOptConstraintRef

Extend `JuMP.UpperBoundRef` to extract a constraint reference for the
upper bound of the original infinite variable of `vref`.

**Example**
```julia-repl
julia> cref = UpperBoundRef(vref)
var <= 1.0
```
"""
function JuMP.UpperBoundRef(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintRef
    return JuMP.UpperBoundRef(infinite_variable_ref(vref))
end

"""
    JuMP.is_fixed(vref::SemiInfiniteVariableRef)::Bool

Extend `JuMP.is_fixed` to return `Bool` whether the original infinite
variable of `vref` is fixed.

**Example**
```julia-repl
julia> is_fixed(vref)
true
```
"""
function JuMP.is_fixed(vref::SemiInfiniteVariableRef)::Bool
    return JuMP.is_fixed(infinite_variable_ref(vref))
end

"""
    JuMP.fix_value(vref::SemiInfiniteVariableRef)::Float64

Extend `JuMP.fix_value` to return the fix value of the original infinite
variable of `vref`. Errors if variable is not fixed.

**Example**
```julia-repl
julia> fix_value(vref)
0.0
```
"""
function JuMP.fix_value(vref::SemiInfiniteVariableRef)::Float64
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.is_fixed(ivref)
        error("Variable $(vref) is not fixed.")
    end
    return JuMP.fix_value(ivref)
end

# Extend to return the index of the fix constraint associated with the original
# infinite variable of `vref`.
function InfiniteOpt._fix_index(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintIndex
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.is_fixed(ivref)
        error("Variable $(vref) is not fixed.")
    end
    return InfiniteOpt._fix_index(ivref)
end

"""
    JuMP.FixRef(vref::SemiInfiniteVariableRef)::InfOptConstraintRef

Extend `JuMP.FixRef` to return the constraint reference of the fix
constraint associated with the original infinite variable of `vref`. Errors
`vref` is not fixed.

**Examples**
```julia-repl
julia> cref = FixRef(vref)
var == 1.0
```
"""
function JuMP.FixRef(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintRef
    return JuMP.FixRef(infinite_variable_ref(vref))
end

# Fallback dispatch
function JuMP.start_value(vref::SemiInfiniteVariableRef)
    return JuMP.start_value(infinite_variable_ref(vref))
end

"""
    start_value_function(vref::SemiInfiniteVariableRef)::Union{Nothing, Function}

Return the function that is used to generate the start values of `vref` for
particular support values. Returns `nothing` if no start behavior has been
specified.

**Example**
```julia-repl
julia> start_value_func(vref)
my_func
```
"""
function start_value_function(vref::SemiInfiniteVariableRef)
    return start_value_function(infinite_variable_ref(vref))
end

"""
    JuMP.is_binary(vref::SemiInfiniteVariableRef)::Bool

Extend `JuMP.is_binary` to return `Bool` whether the original infinite
variable of `vref` is binary.

**Example**
```julia-repl
julia> is_binary(vref)
true
```
"""
function JuMP.is_binary(vref::SemiInfiniteVariableRef)::Bool
    return JuMP.is_binary(infinite_variable_ref(vref))
end

# Extend to return the index of the binary constraint associated with the
# original infinite variable of `vref`.
function InfiniteOpt._binary_index(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintIndex
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.is_binary(ivref)
        error("Variable $(vref) is not binary.")
    end
    return InfiniteOpt._binary_index(ivref)
end

"""
    JuMP.BinaryRef(vref::SemiInfiniteVariableRef)::InfOptConstraintRef

Extend `JuMP.BinaryRef` to return a constraint reference to the
constraint constrainting the original infinite variable of `vref` to be binary.
Errors if one does not exist.

**Example**
```julia-repl
julia> cref = BinaryRef(vref)
var binary
```
"""
function JuMP.BinaryRef(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintRef
    return JuMP.BinaryRef(infinite_variable_ref(vref))
end

"""
    JuMP.is_integer(vref::SemiInfiniteVariableRef)::Bool

Extend `JuMP.is_integer` to return `Bool` whether the original infinite
variable of `vref` is integer.

**Example**
```julia-repl
julia> is_integer(vref)
true
```
"""
function JuMP.is_integer(vref::SemiInfiniteVariableRef)::Bool
    return JuMP.is_integer(infinite_variable_ref(vref))
end

# Extend to return the index of the integer constraint associated with the
# original infinite variable of `vref`.
function InfiniteOpt._integer_index(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintIndex
    ivref = dispatch_variable_ref(infinite_variable_ref(vref))
    if !JuMP.is_integer(ivref)
        error("Variable $(vref) is not an integer.")
    end
    return InfiniteOpt._integer_index(ivref)
end

"""
    JuMP.IntegerRef(vref::SemiInfiniteVariableRef)::InfOptConstraintRef

Extend `JuMP.IntegerRef` to return a constraint reference to the
constraint constrainting the original infinite variable of `vref` to be integer.
Errors if one does not exist.

**Example**
```julia-repl
julia> cref = IntegerRef(vref)
var integer
```
"""
function JuMP.IntegerRef(
    vref::SemiInfiniteVariableRef
    )::InfOptConstraintRef
    return JuMP.IntegerRef(infinite_variable_ref(vref))
end

################################################################################
#                                  DELETION
################################################################################
# Extend _delete_variable_dependencies (for use with JuMP.delete)
function _delete_variable_dependencies(vref::SemiInfiniteVariableRef)::Nothing
    # remove mapping to infinite variable
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref), _semi_infinite_variable_dependencies(ivref))
    # delete associated derivative variables and mapping 
    model = JuMP.owner_model(vref)
    for index in _derivative_dependencies(vref)
        JuMP.delete(model, dispatch_variable_ref(model, index))
    end
    # remove the lookup entry
    delete!(JuMP.owner_model(vref).semi_lookup, (ivref, eval_supports(vref)))
    return
end
