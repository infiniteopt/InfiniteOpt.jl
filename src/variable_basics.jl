################################################################################
#                               VARIABLE DEFINITION
################################################################################
# Define symbol inputs for different variable types
const Infinite = :Infinite
const Point = :Point
const Hold = :Hold

# Fallback _make_variable (the methods are defined in the variable files)
function _make_variable(_error::Function, info::JuMP.VariableInfo, type;
                        extra_kw_args...)
    _error("Unrecognized variable type $type, should be Infinite, " *
           "Point, or Hold.")
end

"""
    JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                        var_type::Symbol;
                        [parameter_refs::Union{GeneralVariableRef,
                                              AbstractArray{<:GeneralVariableRef},
                                              Tuple, Nothing} = nothing,
                        infinite_variable_ref::Union{GeneralVariableRef,
                                                     Nothing} = nothing,
                        parameter_values::Union{Number, AbstractArray{<:Real},
                                                Tuple, Nothing} = nothing,
                        parameter_bounds::Union{ParameterBounds{GeneralVariableRef},
                                                Nothing} = nothing]
                        )::InfOptVariable

Extend the `JuMP.build_variable` function to accomodate `InfiniteOpt`
variable types. Returns the appropriate variable Datatype (i.e.,
[`InfiniteVariable`](@ref), [`PointVariable`](@ref), and
[`HoldVariable`](@ref)). Primarily, this method is to be used internally by the
appropriate constructor macros [`@infinite_variable`](@ref),
[`@point_variable`](@ref), and [`@hold_variable`](@ref). However, it can be
called manually to build `InfiniteOpt` variables. Errors if an unneeded keyword
argument is given or if the keywoard arguments are formatted incorrectly (e.g.,
`parameter_refs` contains repeated parameter references when an infinite variable
is defined). Also errors if needed keyword arguments are negated.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel())
julia> @independent_parameter(m, t in [0, 1])
t

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = t)
InfiniteVariable{Int64,Int64,Int64,Int64,GeneralVariableRef}(VariableInfo{Int64,Int64,Int64,Int64}(false, 0, false, 0, false, 0, false, 0, false, false), (t,), [], [])

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, info, Point, infinite_variable_ref = ivref,
                               parameter_values = 0.5)
PointVariable{Int64,Int64,Int64,Float64,GeneralVariableRef}(VariableInfo{Int64,Int64,Int64,Float64}(false, 0, false, 0, false, 0, true, 0.0, false, false), var_name(t), [0.5])

julia> hd_var = build_variable(error, info, Hold)
HoldVariable{Int64,Int64,Int64,Int64,GeneralVariableRef}(VariableInfo{Int64,Int64,Int64,Int64}(false, 0, false, 0, false, 0, false, 0, false, false), Subdomain bounds (0): )
```
"""
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo,
                             var_type::Symbol;
                             macro_error::Union{Function, Nothing} = nothing,
                             kw_args...)::InfOptVariable
    if macro_error != nothing
        _error = macro_error # replace with macro error function
    end
    # make the variable and conduct necessary checks
    return _make_variable(_error, info, Val(var_type); kw_args...)
end

# Fallback
function _check_and_make_variable_ref(model::InfiniteModel, v::T) where {T}
    throw(ArgumentError("Invalid variable object type `$T`."))
end

"""
    JuMP.add_variable(model::InfiniteModel, var::InfOptVariable,
                      [name::String = ""])::GeneralVariableRef

Extend the [`JuMP.add_variable`](@ref JuMP.add_variable(::JuMP.Model, ::JuMP.ScalarVariable, ::String))
function to accomodate `InfiniteOpt` variable types. Adds a variable to an
infinite model `model` and returns a [`GeneralVariableRef`](@ref).
Primarily intended to be an internal function of the
constructor macros [`@infinite_variable`](@ref), [`@point_variable`](@ref), and
[`@hold_variable`](@ref). However, it can be used in combination with
[`JuMP.build_variable`](@ref) to add variables to an infinite model object.
Errors if invalid parameters reference(s) or an invalid infinite variable
reference is included in `var`.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel())
julia> @infinite_parameter(m, t in [0, 10]);

julia> info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false);

julia> inf_var = build_variable(error, info, Infinite, parameter_refs = t);

julia> ivref = add_variable(m, inf_var, "var_name")
var_name(t)

julia> pt_var = build_variable(error, info, Point, infinite_variable_ref = ivref,
                               parameter_values = 0.5);

julia> pvref = add_variable(m, pt_var, "var_alias")
var_alias

julia> hd_var = build_variable(error, info, Hold);

julia> hvref = add_variable(m, hd_var, "var_name")
var_name
```
"""
function JuMP.add_variable(model::InfiniteModel, var::InfOptVariable,
                           name::String = "")::GeneralVariableRef
    dvref = _check_and_make_variable_ref(model, var)
    JuMP.set_name(dvref, name)
    vindex = JuMP.index(dvref)
    gvref = GeneralVariableRef(model, vindex.value, typeof(vindex))
    if var.info.has_lb
        newset = MOI.GreaterThan(convert(Float64, var.info.lower_bound))
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_lower_bound_index(dvref, JuMP.index(cref))
    end
    if var.info.has_ub
        newset = MOI.LessThan(convert(Float64, var.info.upper_bound))
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_upper_bound_index(dvref, JuMP.index(cref))
    end
    if var.info.has_fix
        newset = MOI.EqualTo(convert(Float64, var.info.fixed_value))
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, newset),
                                   is_info_constr = true)
        _set_fix_index(dvref, JuMP.index(cref))
    end
    if var.info.binary
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.ZeroOne()),
                                   is_info_constr = true)
        _set_binary_index(dvref, JuMP.index(cref))
    elseif var.info.integer
        cref = JuMP.add_constraint(model, JuMP.ScalarConstraint(gvref, MOI.Integer()),
                                   is_info_constr = true)
        _set_integer_index(dvref, JuMP.index(cref))
    end
    return gvref
end

################################################################################
#                            VARIABLE DEPENDENCIES
################################################################################
# Extend _measure_dependencies
function _measure_dependencies(vref::DecisionVariableRef)::Vector{MeasureIndex}
    return _data_object(vref).measure_indices
end

# Extend _constraint_dependencies
function _constraint_dependencies(vref::DecisionVariableRef)::Vector{ConstraintIndex}
    return _data_object(vref).constraint_indices
end

"""
    used_by_measure(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a measure.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref); m.var_to_meas[1] = [1])
julia> used_by_measure(vref)
true
```
"""
function used_by_measure(vref::DecisionVariableRef)::Bool
    return !isempty(_measure_dependencies(vref))
end

"""
    used_by_constraint(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by a constraint.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref))
julia> used_by_constraint(vref)
false
```
"""
function used_by_constraint(vref::DecisionVariableRef)::Bool
    return !isempty(_constraint_dependencies(vref))
end

"""
    used_by_objective(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used by the objective.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref); m.var_in_objective[1] = true)
julia> used_by_objective(vref)
true
```
"""
function used_by_objective(vref::DecisionVariableRef)::Bool
    return _data_object(vref).in_objective
end

"""
    is_used(vref::DecisionVariableRef)::Bool

Return a `Bool` indicating if `vref` is used in the model.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, 0 <= vref))
julia> is_used(vref)
true
```
"""
function is_used(vref::DecisionVariableRef)::Bool
    return used_by_measure(vref) || used_by_constraint(vref) || used_by_objective(vref)
end

################################################################################
#                                VARIABLE NAMING
################################################################################
"""
    JuMP.name(vref::DecisionVariableRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.VariableRef)) to return the names of
`InfiniteOpt` variables.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref, base_name = "var_name"))
julia> name(vref)
"var_name"
```
"""
function JuMP.name(vref::DecisionVariableRef)::String
    return _data_object(vref).name
end

"""
    JuMP.set_name(vref::DecisionVariableRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.VariableRef, ::String)) to set
names of decision variables.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, vref))
julia> set_name(vref, "var_name")

julia> name(vref)
"var_name"
```
"""
function JuMP.set_name(vref::DecisionVariableRef, name::String)::Nothing
    _data_object(vref).name = name
    JuMP.owner_model(vref).name_to_var = nothing
    return
end

# Make a variable reference
function _make_variable_ref(model::InfiniteModel, index::ObjectIndex)::GeneralVariableRef
    return GeneralVariableRef(model, index.value, typeof(index))
end

# Get the name_to_var Dictionary
function _var_name_dict(model::InfiniteModel)::Union{Nothing, Dict{String, ObjectIndex}}
    return model.name_to_var
end

# Update name_to_var
function _update_var_name_dict(model::InfiniteModel, var_dict::MOIUC.CleverDict)::Nothing
    for (index, data_object) in var_dict
        var_name = data_object.name
        if haskey(_var_name_dict(model), var_name)
            # -1 is a special value that means this string does not map to
            # a unique variable name.
            _var_name_dict(model)[var_name] = -1
        else
            _var_name_dict(model)[var_name] = index
        end
    end
    return
end

"""
    JuMP.variable_by_name(model::InfiniteModel,
                          name::String)::Union{GeneralVariableRef, Nothing}

Extend [`JuMP.variable_by_name`](@ref JuMP.variable_by_name(::JuMP.Model, ::String))
for `InfiniteModel` objects. Return the variable reference assoociated with a
variable name. Errors if multiple variables have the same name. Returns nothing
if no such name exists.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; m = InfiniteModel(); @hold_variable(m, base_name = "var_name"))
julia> variable_by_name(m, "var_name")
var_name

julia> variable_by_name(m, "fake_name")

```
"""
function JuMP.variable_by_name(model::InfiniteModel,
                               name::String)::Union{GeneralVariableRef, Nothing}
    if _var_name_dict(model) === nothing
        model.name_to_var = Dict{String, ObjectIndex}()
        _update_name_dict(model, model.infinite_vars)
        _update_name_dict(model, model.reduced_vars)
        _update_name_dict(model, model.point_vars)
        _update_name_dict(model, model.hold_vars)
    end
    index = get(_var_name_dict(model), name, nothing)
    if index isa Nothing
        return nothing
    elseif index == -1
        error("Multiple variables have the name $name.")
    else
        return _make_variable_ref(model, index)
    end
end

################################################################################
#                         VARIABLE OBJECT MODIFICATION
################################################################################
# Extend _set_core_variable_object
function _set_core_variable_object(vref::DecisionVariableRef,
                                   object::VariableData)::Nothing
    _data_object(vref).variable = object
    return
end

################################################################################
#                           VARIABLE INFO METHODS
################################################################################
# Include all the extension functions for manipulating the properties associated
# with VariableInfo
# include("variable_info.jl")


################################################################################
#                          MODEL VARIABLE QUERIES
################################################################################
# TODO make like constraint queries to search based on type
#=
"""
    JuMP.num_variables(model::InfiniteModel)::Int

Extend [`JuMP.num_variables`](@ref JuMP.num_variables(::JuMP.Model)) to return the
number of `InfiniteOpt` variables assigned to `model`.

**Example**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel(); @hold_variable(model, vref[1:3]))
julia> num_variables(model)
3
```
"""
JuMP.num_variables(model::InfiniteModel)::Int = length(model.vars)

"""
    JuMP.all_variables(model::InfiniteModel)::Vector{GeneralVariableRef}

Extend [`JuMP.all_variables`](@ref JuMP.all_variables(::JuMP.Model)) to return a
list of all the variable references associated with `model`.

**Examples**
```julia-repl
julia> all_variables(model)
4-element Array{GeneralVariableRef,1}:
 y(t)
 w(t, x)
 y(0)
 z
```
"""
function JuMP.all_variables(model::InfiniteModel)::Vector{GeneralVariableRef}
    vrefs_list = Vector{GeneralVariableRef}(undef, JuMP.num_variables(model))
    indexes = sort([index for index in keys(model.vars)])
    counter = 1
    for index in indexes
        vrefs_list[counter] = _make_variable_ref(model, index)
        counter += 1
    end
    return vrefs_list
end

################################################################################
#                                 DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.VariableRef)) to delete
`InfiniteOpt` variables and their dependencies. Errors if variable is invalid,
meaning it has already been deleted or it belongs to another model.

**Example**
```julia-repl
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z ≥ 0.0
 g(t) + z ≥ 42.0
 g(0.5) = 0
 t ∈ [0, 6]

julia> delete(model, g)

julia> print(model)
Min measure(t) + z
Subject to
 z ≥ 0.0
 z ≥ 42.0
 t ∈ [0, 6]
```
"""
function JuMP.delete(model::InfiniteModel, vref::InfOptVariableRef)
    @assert JuMP.is_valid(model, vref) "Variable is invalid."
    # update the optimizer model status
    if is_used(vref)
        set_optimizer_model_ready(model, false)
    end
    # remove variable info constraints associated with vref
    if JuMP.has_lower_bound(vref)
        JuMP.delete_lower_bound(vref)
    end
    if JuMP.has_upper_bound(vref)
        JuMP.delete_upper_bound(vref)
    end
    if JuMP.is_fixed(vref)
        JuMP.unfix(vref)
    end
    if JuMP.is_binary(vref)
        JuMP.unset_binary(vref)
    elseif JuMP.is_integer(vref)
        JuMP.unset_integer(vref)
    end
    # remove dependencies from measures and update them
    if used_by_measure(vref)
        for mindex in model.var_to_meas[JuMP.index(vref)]
            if isa(model.measures[mindex].func, InfOptVariableRef)
                model.measures[mindex] = Measure(zero(JuMP.AffExpr),
                                                 model.measures[mindex].data)
            else
                _remove_variable(model.measures[mindex].func, vref)
            end
            JuMP.set_name(MeasureRef(model, mindex),
                           _make_meas_name(model.measures[mindex]))
        end
        # delete mapping
        delete!(model.var_to_meas, JuMP.index(vref))
    end
    # remove dependencies from measures and update them
    if used_by_constraint(vref)
        for cindex in model.var_to_constrs[JuMP.index(vref)]
            if isa(model.constrs[cindex].func, InfOptVariableRef)
                model.constrs[cindex] = JuMP.ScalarConstraint(zero(JuMP.AffExpr),
                                                      model.constrs[cindex].set)
            else
                _remove_variable(model.constrs[cindex].func, vref)
            end
        end
        # delete mapping
        delete!(model.var_to_constrs, JuMP.index(vref))
    end
    # remove from objective if vref is in it
    if used_by_objective(vref)
        if isa(model.objective_function, InfOptVariableRef)
            model.objective_function = zero(JuMP.AffExpr)
        else
            _remove_variable(model.objective_function, vref)
        end
    end
    # do specific updates if vref is infinite
    if isa(vref, InfiniteVariableRef)
        # update parameter mapping
        all_prefs = parameter_list(vref)
        for pref in all_prefs
            filter!(e -> e != JuMP.index(vref),
                    model.param_to_vars[JuMP.index(pref)])
            if length(model.param_to_vars[JuMP.index(pref)]) == 0
                delete!(model.param_to_vars, JuMP.index(pref))
            end
        end
        # delete associated point variables and mapping
        if used_by_point_variable(vref)
            for index in model.infinite_to_points[JuMP.index(vref)]
                JuMP.delete(model, PointVariableRef(model, index))
            end
            delete!(model.infinite_to_points, JuMP.index(vref))
        end
        # delete associated reduced variables and mapping
        if used_by_reduced_variable(vref)
            for index in model.infinite_to_reduced[JuMP.index(vref)]
                JuMP.delete(model, ReducedInfiniteVariableRef(model, index))
            end
            delete!(model.infinite_to_reduced, JuMP.index(vref))
        end
    end
    # update mappings if is point variable
    if isa(vref, PointVariableRef)
        ivref = infinite_variable_ref(vref)
        filter!(e -> e != JuMP.index(vref),
                model.infinite_to_points[JuMP.index(ivref)])
        if length(model.infinite_to_points[JuMP.index(ivref)]) == 0
            delete!(model.infinite_to_points, JuMP.index(ivref))
        end
    end
    # delete the variable information
    delete!(model.var_in_objective, JuMP.index(vref))
    delete!(model.vars, JuMP.index(vref))
    delete!(model.var_to_name, JuMP.index(vref))
    return
end
=#
