"""
    JuMP.owner_model(cref::GeneralConstraintRef)::InfiniteModel

Extend [`JuMP.owner_model`](@ref) to return the infinite model associated with
`cref`.

**Example**
```julia
julia> model = owner_model(cref)
An InfiniteOpt Model
Minimization problem with:
Variables: 3
Objective function type: GlobalVariableRef
`GenericAffExpr{Float64,FiniteVariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
Names registered in the model: g, t, h, x
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
JuMP.owner_model(cref::GeneralConstraintRef)::InfiniteModel = cref.model

"""
    JuMP.index(cref::GeneralConstraintRef)::Int

Extend [`JuMP.index`](@ref) to return the index of an `InfiniteOpt` constraint
`cref`.

**Example**
```julia
julia> index(cref)
2
```
"""
JuMP.index(cref::GeneralConstraintRef)::Int = cref.index

# Extend Base and JuMP functions
function Base.:(==)(v::GeneralConstraintRef, w::GeneralConstraintRef)::Bool
    return v.model === w.model && v.index == w.index && v.shape == w.shape && typeof(v) == typeof(w)
end
Base.broadcastable(cref::GeneralConstraintRef) = Ref(cref)
JuMP.constraint_type(m::InfiniteModel) = GeneralConstraintRef

# Set default bounded set
const default_bounds = Dict{ParameterRef, IntervalSet}()

# This might not be necessary...
function JuMP.build_constraint(_error::Function,
                               v::Union{InfiniteVariableRef, MeasureRef},
                               set::MOI.AbstractScalarSet;
                               parameter_bounds::Dict{ParameterRef,
                                                  IntervalSet} = default_bounds)
    if length(parameter_bounds) != 0
        return BoundedScalarConstraint(v, set, parameter_bounds)
    else
        return JuMP.ScalarConstraint(v, set)
    end
end

"""
    JuMP.build_constraint(_error::Function, expr::InfiniteExpr,
                          set::MOI.AbstractScalarSet;
                          parameter_bounds::Dict{ParameterRef, IntervalSet} = Dict())

Extend [`JuMP.build_constraint`](@ref) to accept the parameter_bounds argument
and return a [`BoundedScalarConstraint`](@ref) if the `parameter_bounds` keyword
argument is specifed or return a [`ScalarConstraint`](@ref) otherwise. This is
primarily intended to work as an internal function for constraint macros.

**Example**
```julia
julia> constr = build_constraint(error, g + x, MOI.EqualTo(42.0),
                               parameter_bounds = Dict(t => IntervalSet(0, 1)));

julia> isa(constr, BoundedScalarConstraint)
true
```
"""
function JuMP.build_constraint(_error::Function,
                               expr::Union{InfiniteExpr, MeasureExpr},
                               set::MOI.AbstractScalarSet;
                               parameter_bounds::Dict{ParameterRef,
                                                  IntervalSet} = default_bounds)
    offset = JuMP.constant(expr)
    JuMP.add_to_expression!(expr, -offset)
    if length(parameter_bounds) != 0
        return BoundedScalarConstraint(expr, MOIU.shift_constant(set, -offset),
                                       parameter_bounds)
    else
        return JuMP.ScalarConstraint(expr, MOIU.shift_constant(set, -offset))
    end
end

# Used to update the model.var_to_constrs field
function _update_var_constr_mapping(vrefs::Vector{<:GeneralVariableRef},
                                    cindex::Int)
    for vref in vrefs
        model = JuMP.owner_model(vref)
        if isa(vref, InfOptVariableRef)
            if haskey(model.var_to_constrs, JuMP.index(vref))
                push!(model.var_to_constrs[JuMP.index(vref)], cindex)
            else
                model.var_to_constrs[JuMP.index(vref)] = [cindex]
            end
        elseif isa(vref, ParameterRef)
            if haskey(model.param_to_constrs, JuMP.index(vref))
                push!(model.param_to_constrs[JuMP.index(vref)], cindex)
            else
                model.param_to_constrs[JuMP.index(vref)] = [cindex]
            end
        elseif isa(vref, MeasureRef)
            if haskey(model.meas_to_constrs, JuMP.index(vref))
                push!(model.meas_to_constrs[JuMP.index(vref)], cindex)
            else
                model.meas_to_constrs[JuMP.index(vref)] = [cindex]
            end
        end
    end
    return
end

# Check that parameter_bounds argument is valid
function _check_bounds(model::InfiniteModel, bounds::Dict)
    prefs = collect(keys(bounds))
    for pref in prefs
        !JuMP.is_valid(model, pref) && error("Parameter bound reference " *
                                             "is invalid.")
        if JuMP.has_lower_bound(pref)
            if bounds[pref].lower_bound < JuMP.lower_bound(pref)
                error("Specified parameter lower bound exceeds that defined " *
                      "for $pref.")
            end
        end
        if JuMP.has_upper_bound(pref)
            if bounds[pref].upper_bound > JuMP.upper_bound(pref)
                error("Specified parameter upper bound exceeds that defined " *
                      "for $pref.")
            end
        end
    end
    return
end

# Extend functions for bounded constraints
JuMP.shape(c::BoundedScalarConstraint) = JuMP.shape(JuMP.ScalarConstraint(c.func, c.set))
JuMP.jump_function(c::BoundedScalarConstraint) = c.func
JuMP.moi_set(c::BoundedScalarConstraint) = c.set

"""
    JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint,
                        name::String = "")

Extend [`JuMP.add_constraint`](@ref) to add a constraint `c` to an infinite model
`model` with name `name`. Returns an appropriate constraint reference whose type
depends on what variables are used to define the constraint. Errors if a vector
constraint is used, the constraint only constains parameters, or if any
variables do not belong to `model`. This is primarily used as an internal
method for the cosntraint macros.

**Example**
```julia
julia> constr = build_constraint(error, g + x, MOI.EqualTo(42));

julia> cref = add_constraint(model, constr, "name")
name : g(t) + x == 42.0
```
"""
function JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint,
                             name::String = "")
    isa(c, JuMP.VectorConstraint) && error("Vector constraints not supported.")
    vrefs = _all_function_variables(c.func)
    isa(vrefs, Vector{ParameterRef}) && error("Constraints cannot contain " *
                                              "only parameters.")
    for vref in vrefs
        JuMP.owner_model(vref) != model && error("Variable $vref does not " *
                                                 "belong to model.")
    end
    if isa(c, BoundedScalarConstraint)
        _check_bounds(model, c.bounds)
    end
    model.next_constr_index += 1
    index = model.next_constr_index
    if length(vrefs) != 0
        _update_var_constr_mapping(vrefs, index)
    end
    if c.func isa InfiniteExpr
        cref = InfiniteConstraintRef(model, index, JuMP.shape(c))
    elseif c.func isa MeasureExpr
        cref = MeasureConstraintRef(model, index, JuMP.shape(c))
    else
        cref = FiniteConstraintRef(model, index, JuMP.shape(c))
    end
    model.constrs[index] = c
    JuMP.set_name(cref, name)
    model.constr_in_var_info[index] = false
    set_optimizer_model_ready(model, false)
    return cref
end

"""
    JuMP.delete(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.delete` function to accomodate our new constraint types.
"""
function JuMP.delete(model::InfiniteModel, cref::GeneralConstraintRef)
    @assert JuMP.is_valid(model, cref)
    all_vrefs = _all_function_variables(model.constrs[JuMP.index(cref)].func)
    for vref in all_vrefs
        if isa(vref, InfOptVariableRef)
            filter!(e -> e != JuMP.index(cref), model.var_to_constrs[JuMP.index(vref)])
            if length(model.var_to_constrs[JuMP.index(vref)]) == 0
                delete!(model.var_to_constrs, JuMP.index(vref))
            end
        elseif isa(vref, ParameterRef)
            filter!(e -> e != JuMP.index(cref), model.param_to_constrs[JuMP.index(vref)])
            if length(model.param_to_constrs[JuMP.index(vref)]) == 0
                delete!(model.param_to_constrs, JuMP.index(vref))
            end
        elseif isa(vref, MeasureRef)
            filter!(e -> e != JuMP.index(cref), model.meas_to_constrs[JuMP.index(vref)])
            if length(model.meas_to_constrs[JuMP.index(vref)]) == 0
                delete!(model.meas_to_constrs, JuMP.index(vref))
            end
        end
    end
    delete!(model.constrs, JuMP.index(cref))
    delete!(model.constr_to_name, JuMP.index(cref))
    delete!(model.constr_in_var_info, JuMP.index(cref))
    set_optimizer_model_ready(model, false)
    return
end

"""
    JuMP.is_valid(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.is_valid` function to accomodate our new constraint types.
"""
function JuMP.is_valid(model::InfiniteModel, cref::GeneralConstraintRef)
    return (model === JuMP.owner_model(cref) && JuMP.index(cref) in keys(model.constrs))
end

"""
    JuMP.constraint_object(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.constraint_object` function to accomodate our new constraint types.
"""
function JuMP.constraint_object(cref::GeneralConstraintRef)
    return JuMP.owner_model(cref).constrs[JuMP.index(cref)]
end

"""
    JuMP.name(cref::GeneralConstraintRef)::String

Extend [`JuMP.name`](@ref) to return the name of an `InfiniteOpt` constraint.

**Example**
```julia
julia> name(cref)
constr_name
```
"""
function JuMP.name(cref::GeneralConstraintRef)::String
    return JuMP.owner_model(cref).constr_to_name[JuMP.index(cref)]
end

"""
    JuMP.set_name(cref::GeneralConstraintRef, name::String)

Extend [`JuMP.set_name`](@ref) to specify the name of a constraint `cref`.

**Example**
```julia
julia> set_name(cref, "new_name")

julia> name(cref)
new_name
```
"""
function JuMP.set_name(cref::GeneralConstraintRef, name::String)
    JuMP.owner_model(cref).constr_to_name[JuMP.index(cref)] = name
    JuMP.owner_model(cref).name_to_constr = nothing
    return
end

# Return the appropriate constraint reference given the index and model
function _make_constraint_ref(model::InfiniteModel,
                              index::Int)::GeneralConstraintRef
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index,
                                     JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index,
                                    JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index,
                                   JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.constraint_by_name(model::InfiniteModel,
                            name::String)::Union{GeneralConstraintRef, Nothing}

Extend [`JuMP.constraint_by_name`](@ref) to return the constraint reference
associated with `name` if one exists or returns nothing. Errors if more than
one constraint uses the same name.

**Example**
```julia
julia> constraint_by_name(model, "constr_name")
constr_name : x + pt == 3.0
```
"""
function JuMP.constraint_by_name(model::InfiniteModel, name::String)
    if model.name_to_constr === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_constr = Dict{String, Int}()
        for (constr, constr_name) in model.constr_to_name
            if haskey(model.name_to_constr, constr_name)
                # -1 is a special value that means this string does not map to
                # a unique constraint name.
                model.name_to_constr[constr_name] = -1
            else
                model.name_to_constr[constr_name] = constr
            end
        end
    end
    index = get(model.name_to_constr, name, nothing)
    if index isa Nothing
        return nothing
    elseif index == -1
        error("Multiple constraints have the name $name.")
    else
        return _make_constraint_ref(model, index)
    end
end

"""
    JuMP.num_constraints(model::InfiniteModel, function_type, set_type)::Int64
Extend the `JuMP.num_constraints` function to accomodate our new constraint types.
"""
function JuMP.num_constraints(model::InfiniteModel, function_type::Type{<:JuMP.AbstractJuMPScalar},
                              set_type::Type{<:MOI.AbstractSet})::Int64
    counter = 0
    for k in keys(model.constrs)
        if typeof(model.constrs[k].func) <: function_type && typeof(model.constrs[k].set) <: set_type
            counter += 1
        end
    end
    return counter
end

"""
    JuMP.num_constraints(model::InfiniteModel, function_type)::Int64
Extend the `JuMP.num_constraints` function to search by function types for all MOI sets.
"""
function JuMP.num_constraints(model::InfiniteModel, function_type::Type{<:JuMP.AbstractJuMPScalar})
    return JuMP.num_constraints(model, function_type, MOI.AbstractSet)
end

"""
    JuMP.num_constraints(model::InfiniteModel, function_type)::Int64
Extend the `JuMP.num_constraints` function to search by MOI set type for all function types.
"""
function JuMP.num_constraints(model::InfiniteModel, set_type::Type{<:MOI.AbstractSet})
    return JuMP.num_constraints(model, JuMP.AbstractJuMPScalar, set_type)
end

"""
    JuMP.num_constraints(model::InfiniteModel)::Int64
Extend the `JuMP.num_constraints` function to return the total number of constraints.
"""
function JuMP.num_constraints(model::InfiniteModel)
    return length(keys(model.constrs))
end

"""
    JuMP.all_constraints(model::InfiniteModel, function_type, set_type)::Vector{<:GeneralConstraintRef}
Extend the `JuMP.all_constraints` function to accomodate our new constraint types.
"""
function JuMP.all_constraints(model::InfiniteModel, function_type::Type{<:JuMP.AbstractJuMPScalar},
                              set_type::Type{<:MOI.AbstractSet})::Vector{<:GeneralConstraintRef}
    constr_list = Vector{GeneralConstraintRef}(undef, JuMP.num_constraints(model, function_type, set_type))
    indexes = sort([index for index in keys(model.constrs)])
    counter = 1
    for index in indexes
        if typeof(model.constrs[index].func) <: function_type && typeof(model.constrs[index].set) <: set_type
            constr_list[counter] = _make_constraint_ref(model, index)
            counter += 1
        end
    end
    return constr_list
end

"""
    JuMP.all_constraints(model::InfiniteModel, function_type)::Vector{<:GeneralConstraintRef}
Extend the `JuMP.all_constraints` function to search by function types for all MOI sets.
"""
function JuMP.all_constraints(model::InfiniteModel, function_type::Type{<:JuMP.AbstractJuMPScalar})
    return JuMP.all_constraints(model, function_type, MOI.AbstractSet)
end

"""
    JuMP.all_constraints(model::InfiniteModel, function_type)::Vector{<:GeneralConstraintRef}
Extend the `JuMP.all_constraints` function to search by MOI set type for all function types.
"""
function JuMP.all_constraints(model::InfiniteModel, set_type::Type{<:MOI.AbstractSet})
    return JuMP.all_constraints(model, JuMP.AbstractJuMPScalar, set_type)
end

"""
    JuMP.all_constraints(model::InfiniteModel)::Vector{<:GeneralConstraintRef}
Extend the `JuMP.all_constraints` function to return the total number of constraints.
"""
function JuMP.all_constraints(model::InfiniteModel)
    return JuMP.all_constraints(model, JuMP.AbstractJuMPScalar, MOI.AbstractSet)
end

"""
    JuMP.list_of_constraint_types(model::InfiniteModel)
Extend the `JuMP.list_of_constraint_types` function to accomodate our new constraint types.
"""
function JuMP.list_of_constraint_types(model::InfiniteModel)
    type_list = Vector{Tuple{DataType, DataType}}(undef, JuMP.num_constraints(model))
    counter = 1
    for k in keys(model.constrs)
        type_list[counter] = (typeof(model.constrs[k].func), typeof(model.constrs[k].set))
        counter += 1
    end
    return unique(type_list)
end
