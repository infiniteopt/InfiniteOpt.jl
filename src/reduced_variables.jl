# Helper function to get reduced variable info
function _reduced_info(vref::ReducedInfiniteVariableRef)::ReducedInfiniteInfo
    return JuMP.owner_model(vref).reduced_info[JuMP.index(vref)]
end

"""
    infinite_variable_ref(vref::ReducedInfiniteVariableRef)::InfiniteVariableRef

Return the `InfiniteVariableRef` associated with the reduced infinite variable
`vref`.

**Example**
```julia
julia> infinite_variable_ref(vref)
g(t, x)
```
"""
function infinite_variable_ref(vref::ReducedInfiniteVariableRef)::InfiniteVariableRef
    return _reduced_info(vref).infinite_variable_ref
end

"""
    eval_supports(vref::ReducedInfiniteVariableRef)::Dict

Return the evaluation supports associated with the reduced infinite variable
`vref`.

**Example**
```julia
julia> eval_supports(vref)
Dict{Int64,Float64} with 1 entry:
  1 => 0.5
```
"""
function eval_supports(vref::ReducedInfiniteVariableRef)::Dict
    return _reduced_info(vref).eval_supports
end

"""
    parameter_refs(vref::ReducedInfiniteVariableRef)::Tuple

Return the `ParameterRef`(s) associated with the reduced infinite variable
`vref`. This is formatted as a Tuple of containing the parameter references as
they were inputted to define the untracripted infinite variable except, the
evaluated parameters are excluded.

**Example**
```julia
julia> parameter_refs(vref)
(t,   [2]  =  x[2]
  [1]  =  x[1])
```
"""
function parameter_refs(vref::ReducedInfiniteVariableRef)
    orig_prefs = parameter_refs(infinite_variable_ref(vref))
    prefs = Tuple(orig_prefs[i] for i = 1:length(orig_prefs) if !haskey(eval_supports(vref), i))
    return prefs
end

"""
    JuMP.name(vref::ReducedInfiniteVariableRef)::String

Extend `JuMP.name` to return name of reduced infinite variable references. This
is used when displaying measure expansions that contain such variables.

**Exanple**
```julia
julia> name(rvref)
g(1.25, x)
```
"""
function JuMP.name(vref::ReducedInfiniteVariableRef)::String
    root_name = _root_name(infinite_variable_ref(vref))
    prefs = parameter_refs(infinite_variable_ref(vref))
    param_names = [_root_name(first(pref)) for pref in prefs]
    for (k, v) in eval_supports(vref)
        param_names[k] = string(v)
    end
    param_name_tuple = "("
    for i = 1:length(param_names)
        if i != length(param_names)
            param_name_tuple *= string(param_names[i], ", ")
        else
            param_name_tuple *= string(param_names[i])
        end
    end
    param_name_tuple *= ")"
    return string(root_name, param_name_tuple)
end

"""
    JuMP.has_lower_bound(vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.has_lower_bound`](@ref) to return a `Bool` whether the original
infinite variable of `vref` has a lower bound.

**Example**
```julia
julia> has_lower_bound(vref)
true
```
"""
function JuMP.has_lower_bound(vref::ReducedInfiniteVariableRef)::Bool
    return JuMP.has_lower_bound(infinite_variable_ref(vref))
end

"""
    JuMP.lower_bound(vref::ReducedInfiniteVariableRef::Float64

Extend [`JuMP.lower_bound`](@ref) to return the lower bound of the original
infinite variable of `vref`. Errors if `vref` doesn't have a lower bound.

**Example**
```julia
julia> lower_bound(vref)
0.0
```
"""
function JuMP.lower_bound(vref::ReducedInfiniteVariableRef)::Float64
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return JuMP.lower_bound(infinite_variable_ref(vref))
end

# Extend to return the index of the lower bound constraint associated with the
# original infinite variable of `vref`.
function JuMP._lower_bound_index(vref::ReducedInfiniteVariableRef)::Int
    if !JuMP.has_lower_bound(vref)
        error("Variable $(vref) does not have a lower bound.")
    end
    return JuMP._lower_bound_index(infinite_variable_ref(vref))
end

"""
    JuMP.LowerBoundRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef

Extend [`JuMP.LowerBoundRef`](@ref) to extract a constraint reference for the
lower bound of the original infinite variable of `vref`.

**Example**
```julia
julia> cref = LowerBoundRef(vref)
var >= 0.0
```
"""
function JuMP.LowerBoundRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef
    return JuMP.LowerBoundRef(infinite_variable_ref(vref))
end

"""
    JuMP.has_upper_bound(vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.has_upper_bound`](@ref) to return a `Bool` whether the original
infinite variable of `vref` has an upper bound.

**Example**
```julia
julia> has_upper_bound(vref)
true
```
"""
function JuMP.has_upper_bound(vref::ReducedInfiniteVariableRef)::Bool
    return JuMP.has_upper_bound(infinite_variable_ref(vref))
end

"""
    JuMP.upper_bound(vref::ReducedInfiniteVariableRef)::Float64

Extend [`JuMP.upper_bound`](@ref) to return the upper bound of the original
infinite variable of `vref`. Errors if `vref` doesn't have a upper bound.

**Example**
```julia
julia> upper_bound(vref)
0.0
```
"""
function JuMP.upper_bound(vref::ReducedInfiniteVariableRef)::Float64
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return JuMP.upper_bound(infinite_variable_ref(vref))
end

# Extend to return the index of the upper bound constraint associated with the
# original infinite variable of `vref`.
function JuMP._upper_bound_index(vref::ReducedInfiniteVariableRef)::Int
    if !JuMP.has_upper_bound(vref)
        error("Variable $(vref) does not have a upper bound.")
    end
    return JuMP._upper_bound_index(infinite_variable_ref(vref))
end

"""
    JuMP.UpperBoundRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef

Extend [`JuMP.UpperBoundRef`](@ref) to extract a constraint reference for the
upper bound of the original infinite variable of `vref`.

**Example**
```julia
julia> cref = UpperBoundRef(vref)
var <= 1.0
```
"""
function JuMP.UpperBoundRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef
    return JuMP.UpperBoundRef(infinite_variable_ref(vref))
end

"""
    JuMP.is_fixed(vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.is_fixed`](@ref) to return `Bool` whether the original infinite
variable of `vref` is fixed.

**Example**
```julia
julia> is_fixed(vref)
true
```
"""
function JuMP.is_fixed(vref::ReducedInfiniteVariableRef)::Bool
    return JuMP.is_fixed(infinite_variable_ref(vref))
end

"""
    JuMP.fix_value(vref::ReducedInfiniteVariableRef)::Float64

Extend [`JuMP.fix_value`](@ref) to return the fix value of the original infinite
variable of `vref`. Errors if variable is not fixed.

**Example**
```julia
julia> fix_value(vref)
0.0
```
"""
function JuMP.fix_value(vref::ReducedInfiniteVariableRef)::Float64
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return JuMP.fix_value(infinite_variable_ref(vref))
end

# Extend to return the index of the fix constraint associated with the original
# infinite variable of `vref`.
function JuMP._fix_index(vref::ReducedInfiniteVariableRef)::Int
    if !JuMP.is_fixed(vref)
        error("Variable $(vref) is not fixed.")
    end
    return JuMP._fix_index(infinite_variable_ref(vref))
end

"""
    JuMP.FixRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef

Extend [`JuMP.FixRef`](@ref) to return the constraint reference of the fix
constraint associated with the original infinite variable of `vref`. Errors
`vref` is not fixed.

**Examples**
```julia
julia> cref = FixRef(vref)
var == 1.0
```
"""
function JuMP.FixRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef
    return JuMP.FixRef(infinite_variable_ref(vref))
end

"""
    JuMP.start_value(vref::ReducedInfiniteVariableRef)::Union{Nothing, Float64}

Extend [`JuMP.start_value`](@ref) to return starting value of the original
infinite variable of `vref` if it has one. Returns `nothing` otherwise.

**Example**
```julia
julia> start_value(vref)
0.0
```
"""
function JuMP.start_value(vref::ReducedInfiniteVariableRef)::Union{Nothing, Float64}
    return JuMP.start_value(infinite_variable_ref(vref))
end

"""
    JuMP.is_binary(vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.is_binary`](@ref) to return `Bool` whether the original infinite
variable of `vref` is binary.

**Example**
```julia
julia> is_binary(vref)
true
```
"""
function JuMP.is_binary(vref::ReducedInfiniteVariableRef)::Bool
    return JuMP.is_binary(infinite_variable_ref(vref))
end

# Extend to return the index of the binary constraint associated with the
# original infinite variable of `vref`.
function JuMP._binary_index(vref::ReducedInfiniteVariableRef)::Int
    if !JuMP.is_binary(vref)
        error("Variable $(vref) is not binary.")
    end
    return JuMP._binary_index(infinite_variable_ref(vref))
end

"""
    JuMP.BinaryRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef

Extend [`JuMP.BinaryRef`](@ref) to return a constraint reference to the
constraint constrainting the original infinite variable of `vref` to be binary.
Errors if one does not exist.

**Example**
```julia
julia> cref = BinaryRef(vref)
var binary
```
"""
function JuMP.BinaryRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef
    return JuMP.BinaryRef(infinite_variable_ref(vref))
end

"""
    JuMP.is_integer(vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.is_integer`](@ref) to return `Bool` whether the original infinite
variable of `vref` is integer.

**Example**
```julia
julia> is_integer(vref)
true
```
"""
function JuMP.is_integer(vref::ReducedInfiniteVariableRef)::Bool
    return JuMP.is_integer(infinite_variable_ref(vref))
end

# Extend to return the index of the integer constraint associated with the
# original infinite variable of `vref`.
function JuMP._integer_index(vref::ReducedInfiniteVariableRef)::Int
    if !JuMP.is_integer(vref)
        error("Variable $(vref) is not an integer.")
    end
    return JuMP._integer_index(infinite_variable_ref(vref))
end

"""
    JuMP.IntegerRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef

Extend [`JuMP.IntegerRef`](@ref) to return a constraint reference to the
constraint constrainting the original infinite variable of `vref` to be integer.
Errors if one does not exist.

**Example**
```julia
julia> cref = IntegerRef(vref)
var integer
```
"""
function JuMP.IntegerRef(vref::ReducedInfiniteVariableRef)::GeneralConstraintRef
    return JuMP.IntegerRef(infinite_variable_ref(vref))
end

"""
    used_by_constraint(vref::ReducedInfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a constraint.

**Example**
```julia
julia> used_by_constraint(vref)
false
```
"""
function used_by_constraint(vref::ReducedInfiniteVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).reduced_to_constrs, JuMP.index(vref))
end

"""
    used_by_measure(vref::ReducedInfiniteVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a measure.

**Example**
```julia
julia> used_by_measure(vref)
true
```
"""
function used_by_measure(vref::ReducedInfiniteVariableRef)::Bool
    return haskey(JuMP.owner_model(vref).reduced_to_meas, JuMP.index(vref))
end

"""
    JuMP.is_valid(model::InfiniteModel, vref::ReducedInfiniteVariableRef)::Bool

Extend [`JuMP.is_valid`](@ref) to accomodate reduced infinite variables.

**Example**
```julia
julia> is_valid(model, vref)
true
```
"""
function JuMP.is_valid(model::InfiniteModel,
                       vref::ReducedInfiniteVariableRef)::Bool
    return (model === JuMP.owner_model(vref) && JuMP.index(vref) in keys(model.reduced_info))
end

"""
    JuMP.delete(model::InfiniteModel, vref::ReducedInfiniteVariableRef)

Extend [`JuMP.delete`](@ref) to delete reduced infinite variables and its
dependencies. Errors if `vref` is invalid, meaning it has already been deleted
or it belongs to another model.

**Example**
``julia
julia> print(model)
Min measure(g(0, t)*t + g(1, t)*t) + z
Subject to
 z >= 0.0
 g(0, t) + g(1, t) == 0
 g(x, t) + z >= 42.0
 g(0.5, 0.5) == 0
 t in [0, 6]
 x in [0, 1]

julia> delete(model, rvref1)

julia> print(model)
Min measure(t + g(1, t)*t) + z
Subject to
 z >= 0.0
 g(1, t) == 0
 g(x, t) + z >= 42.0
 g(0.5, 0.5) == 0
 t in [0, 6]
 x in [0, 1]
```
"""
function JuMP.delete(model::InfiniteModel, vref::ReducedInfiniteVariableRef)
    # check valid reference
    @assert JuMP.is_valid(model, vref) "Invalid variable reference."
    # remove from measures if used
    if used_by_measure(vref)
        for mindex in model.reduced_to_meas[JuMP.index(vref)]
            if isa(model.measures[mindex].func, ReducedInfiniteVariableRef)
                model.measures[mindex] = Measure(zero(JuMP.AffExpr),
                                                 model.measures[mindex].data)
            else
                _remove_variable(model.measures[mindex].func, vref)
            end
            JuMP.set_name(MeasureRef(model, mindex),
                          _make_meas_name(model.measures[mindex]))
        end
        # delete mapping
        delete!(model.reduced_to_meas, JuMP.index(vref))
    end
    # remove from constraints if used
    if used_by_constraint(vref)
        for cindex in model.reduced_to_constrs[JuMP.index(vref)]
            if isa(model.constrs[cindex].func, ReducedInfiniteVariableRef)
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr),
                                                      model.constrs[cindex].set)
            else
                _remove_variable(model.constrs[cindex].func, vref)
            end
        end
        # delete mapping
        delete!(model.reduced_to_constrs, JuMP.index(vref))
    end
    # remove mapping to infinite variable
    ivref = infinite_variable_ref(vref)
    filter!(e -> e != JuMP.index(vref),
            model.infinite_to_reduced[JuMP.index(ivref)])
    if length(model.infinite_to_reduced[JuMP.index(ivref)]) == 0
        delete!(model.infinite_to_reduced, JuMP.index(ivref))
    end
    # delete the info
    delete!(model.reduced_info, JuMP.index(vref))
    return
end
