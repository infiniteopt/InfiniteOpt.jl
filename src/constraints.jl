################################################################################
#                             BASIC REFERENCE METHODS
################################################################################
"""
    JuMP.owner_model(cref::InfOptConstraintRef)::InfiniteModel

Extend [`JuMP.owner_model`](@ref JuMP.owner_model(::JuMP.ConstraintRef)) to
return the infinite model associated with `cref`.

**Example**
```julia-repl
julia> model = owner_model(cref)
An InfiniteOpt Model
Minimization problem with:
Finite Parameters: 0
Infinite Parameters: 3
Variables: 3
Measures: 0
Objective function type: GeneralVariableRef
`GenericAffExpr{Float64,GeneralVariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
Names registered in the model: g, t, h, x
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
JuMP.owner_model(cref::InfOptConstraintRef)::InfiniteModel = cref.model

"""
    JuMP.index(cref::InfOptConstraintRef)::ConstraintIndex

Extend [`JuMP.index`](@ref JuMP.index(::JuMP.ConstraintRef)) to return the index
of an `InfiniteOpt` constraint `cref`.

**Example**
```julia-repl
julia> index(cref)
ConstrainIndex(2)
```
"""
JuMP.index(cref::InfOptConstraintRef)::ConstraintIndex = cref.index

# Extend Base and JuMP functions
function Base.:(==)(v::InfOptConstraintRef, w::InfOptConstraintRef)::Bool
    return v.model === w.model && v.index == w.index && v.shape == w.shape &&
           typeof(v) == typeof(w)
end
Base.broadcastable(cref::InfOptConstraintRef) = Ref(cref)
JuMP.constraint_type(m::InfiniteModel) = InfOptConstraintRef

################################################################################
#                             CORE OBJECT METHODS
################################################################################
# Extend _add_data_object
function _add_data_object(model::InfiniteModel,
                          object::ConstraintData)::ConstraintIndex
    return MOIUC.add_item(model.constraints, object)
end

# Extend _data_dictionary
function _data_dictionary(cref::InfOptConstraintRef
    )::MOIUC.CleverDict{ConstraintIndex, ConstraintData}
    return JuMP.owner_model(cref).constraints
end

# Extend _data_object
function _data_object(cref::InfOptConstraintRef)::ConstraintData
  object = _get(_data_dictionary(cref), JuMP.index(cref), nothing)
  isnothing(object) && error("Invalid constraint reference, cannot find " *
                        "corresponding constraint in the model. This is likely " *
                        "caused by using the reference of a deleted constraint.")
  return object
end

# Return the core constraint object
function _core_constraint_object(cref::InfOptConstraintRef)::JuMP.AbstractConstraint
    return _data_object(cref).constraint
end

# Update the core constraint object
function _set_core_constraint_object(cref::InfOptConstraintRef,
                                     constr::JuMP.AbstractConstraint)::Nothing
    _data_object(cref).constraint = constr
    return
end

# Extend _object_numbers
function _object_numbers(cref::InfOptConstraintRef)::Vector{Int}
    return _data_object(cref).object_nums
end

# Extend _measure_dependencies
function _measure_dependencies(cref::InfOptConstraintRef)::Vector{MeasureIndex}
    return _data_object(cref).measure_indices
end

# Return if this constraint is an info constraint
function _is_info_constraint(cref::InfOptConstraintRef)::Bool
    return _data_object(cref).is_info_constraint
end

# Extend _delete_data_object
function _delete_data_object(cref::InfOptConstraintRef)::Nothing
    delete!(_data_dictionary(cref), JuMP.index(cref))
    return
end

################################################################################
#                             DEFINITION METHODS
################################################################################
# This might not be necessary...
function JuMP.build_constraint(_error::Function,
    v::GeneralVariableRef,
    set::MOI.AbstractScalarSet;
    parameter_bounds::ParameterBounds{GeneralVariableRef} = ParameterBounds()
    )::JuMP.AbstractConstraint
    # make the constraint
    if length(parameter_bounds) != 0
        _check_bounds(parameter_bounds)
        return BoundedScalarConstraint(v, set, parameter_bounds,
                                       copy(parameter_bounds))
    else
        return JuMP.ScalarConstraint(v, set)
    end
end

"""
    JuMP.build_constraint(_error::Function, expr::JuMP.AbstractJuMPScalar,
                          set::MOI.AbstractScalarSet;
                          [parameter_bounds::ParameterBounds = ParameterBounds()]
                          )::JuMP.AbstractConstraint

Extend `JuMP.build_constraint` to accept the ```parameter_bounds``` argument
and return a [`BoundedScalarConstraint`](@ref) if the ```parameter_bounds``` keyword
argument is specifed or return a [`JuMP.ScalarConstraint`](@ref) otherwise. This is
primarily intended to work as an internal function for constraint macros.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_variable(model, g(t));

julia> @hold_variable(model, x);

julia> constr = build_constraint(error, g + x, MOI.EqualTo(42.0),
              parameter_bounds = ParameterBounds((t => IntervalSet(0, 1),)));

julia> isa(constr, BoundedScalarConstraint)
true
```
"""
function JuMP.build_constraint(_error::Function,
    expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}},
    set::MOI.AbstractScalarSet;
    parameter_bounds::ParameterBounds{GeneralVariableRef} = ParameterBounds()
    )::JuMP.AbstractConstraint where {C, V <: GeneralVariableRef}
    # make the constraint
    offset = JuMP.constant(expr)
    JuMP.add_to_expression!(expr, -offset)
    if !isempty(parameter_bounds)
        _check_bounds(parameter_bounds)
        return BoundedScalarConstraint(expr, MOIU.shift_constant(set, -offset),
                                       parameter_bounds, copy(parameter_bounds))
    else
        return JuMP.ScalarConstraint(expr, MOIU.shift_constant(set, -offset))
    end
end

## Perfrom bound checks and update them if needed, then return the updated constraint
# BoundedScalarConstraint
function _check_and_update_bounds(model::InfiniteModel,
                                  c::BoundedScalarConstraint,
                                  vrefs::Vector{GeneralVariableRef}
                                  )::JuMP.AbstractConstraint
    # look for bounded hold variables and update bounds
    for vref in vrefs
        _update_var_bounds(vref, c.bounds)
    end
    # now validate and return
    _validate_bounds(model, c.bounds)
    # TODO should we check that bounds don't violate point variables?
    return c
end

# ScalarConstraint
function _check_and_update_bounds(model::InfiniteModel,
                                  c::JuMP.ScalarConstraint,
                                  vrefs::Vector{GeneralVariableRef}
                                  )::JuMP.AbstractConstraint
    bounds = ParameterBounds()
    # check for bounded hold variables and build the intersection of the bounds
    for vref in vrefs
        _update_var_bounds(vref, bounds)
    end
    # if we added bounds, change to a bounded constraint and validate
    if !isempty(bounds)
        new_c = BoundedScalarConstraint(c.func, c.set, bounds, ParameterBounds())
        _validate_bounds(model, new_c.bounds)
        return new_c
    end
    # TODO should we check that bounds don't violate point variables?
    return c
end

# Used to update the model.var_to_constrs field
function _update_var_constr_mapping(vrefs::Vector{GeneralVariableRef},
                                    cref::InfOptConstraintRef)::Nothing
    for vref in vrefs
        dvref = dispatch_variable_ref(vref)
        push!(_constraint_dependencies(dvref), JuMP.index(cref))
        if dvref isa MeasureRef
            push!(_measure_dependencies(cref), JuMP.index(dvref))
        end
    end
    return
end

# Extend functions for bounded constraints
JuMP.shape(c::BoundedScalarConstraint) = JuMP.shape(JuMP.ScalarConstraint(c.func, c.set))
JuMP.jump_function(c::BoundedScalarConstraint) = c.func
JuMP.moi_set(c::BoundedScalarConstraint) = c.set
parameter_bounds(c::BoundedScalarConstraint) = c.bounds
original_parameter_bounds(c::BoundedScalarConstraint) = c.orig_bounds

"""
    JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint,
                        [name::String = ""])::InfOptConstraintRef

Extend [`JuMP.add_constraint`](@ref JuMP.add_constraint(::JuMP.Model, ::JuMP.AbstractConstraint, ::String))
to add a constraint `c` to an infinite model
`model` with name `name`. Returns an appropriate constraint reference whose type
depends on what variables are used to define the constraint. Errors if a vector
constraint is used, the constraint only constains parameters, or if any
variables do not belong to `model`. This is primarily used as an internal
method for the cosntraint macros.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10]);

julia> @infinite_variable(model, g(t));

julia> @hold_variable(model, x);

julia> constr = build_constraint(error, g + x, MOI.EqualTo(42));

julia> cref = add_constraint(model, constr, "name")
name : g(t) + x = 42.0, ∀ t ∈ [0, 10]
```
"""
function JuMP.add_constraint(model::InfiniteModel,
                             c::JuMP.AbstractConstraint,
                             name::String = "";
                             is_info_constr::Bool = false
                             )::InfOptConstraintRef
    # gather the unique list of variable references for testing and mapping
    vrefs = _all_function_variables(JuMP.jump_function(c))
    # test in the model
    for vref in vrefs
        JuMP.check_belongs_to_model(vref, model)
    end
    # test parameter bounds and update with hold variable bounds
    if model.has_hold_bounds
        c = _check_and_update_bounds(model, c, vrefs)
    elseif c isa BoundedScalarConstraint
        _validate_bounds(model, parameter_bounds(c))
    end
    # get the parameter object numbers
    object_nums = sort!(_object_numbers(vrefs))
    # add the constaint to the model
    constr_object = ConstraintData(c, object_nums, name, MeasureIndex[],
                                   is_info_constr)
    cindex = _add_data_object(model, constr_object)
    cref = InfOptConstraintRef(model, cindex, JuMP.shape(c))
    # update the variable mappings and model status
    _update_var_constr_mapping(vrefs, cref)
    set_optimizer_model_ready(model, false)
    return cref
end

# Fallback for vector constraints
function JuMP.add_constraint(model::InfiniteModel, c::JuMP.VectorConstraint,
                             name::String = "";
                             is_info_constr::Bool = false)
    error("Vector constraints not supported in InfiniteOpt.")
end

################################################################################
#                                 JuMP METHODS
################################################################################
"""
    JuMP.is_valid(model::InfiniteModel, cref::InfOptConstraintRef)::Bool

Extend [`JuMP.is_valid`](@ref JuMP.is_valid(::JuMP.Model, ::JuMP.ConstraintRef{JuMP.Model}))
to return `Bool` whether an `InfiniteOpt` constraint reference is valid.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel(); cref = @constraint(model, 2 == 2))
julia> is_valid(model, cref)
true
```
"""
function JuMP.is_valid(model::InfiniteModel, cref::InfOptConstraintRef)::Bool
    return (model === JuMP.owner_model(cref) &&
            JuMP.index(cref) in keys(_data_dictionary(cref)))
end

"""
    JuMP.constraint_object(cref::InfOptConstraintRef)::JuMP.AbstractConstraint

Extend [`JuMP.constraint_object`](@ref JuMP.constraint_object)
to return the constraint object associated with `cref`.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10]);

julia> @hold_variable(model, x <= 1);

julia> cref = UpperBoundRef(x);

julia> obj = constraint_object(cref)
ScalarConstraint{GeneralVariableRef,MathOptInterface.LessThan{Float64}}(x,
MathOptInterface.LessThan{Float64}(1.0))
```
"""
function JuMP.constraint_object(cref::InfOptConstraintRef)::JuMP.AbstractConstraint
    return _core_constraint_object(cref)
end

"""
    JuMP.name(cref::InfOptConstraintRef)::String

Extend [`JuMP.name`](@ref JuMP.name(::JuMP.ConstraintRef{JuMP.Model,<:JuMP._MOICON})
to return the name of an `InfiniteOpt` constraint.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel(); cref = @constraint(model, constr_name, 2 == 2))
julia> name(cref)
"constr_name"
```
"""
function JuMP.name(cref::InfOptConstraintRef)::String
    object = _get(_data_dictionary(cref), JuMP.index(cref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(cref::InfOptConstraintRef, name::String)::Nothing

Extend [`JuMP.set_name`](@ref JuMP.set_name(::JuMP.ConstraintRef{JuMP.Model,<:JuMP._MOICON}, ::String))
to specify the name of a constraint `cref`.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel(); cref = @constraint(model, 2 == 2))
julia> set_name(cref, "new_name")

julia> name(cref)
"new_name"
```
"""
function JuMP.set_name(cref::InfOptConstraintRef, name::String)::Nothing
    _data_object(cref).name = name
    JuMP.owner_model(cref).name_to_constr = nothing
    return
end

# Return a constraint set with an updated value
function _set_set_value(set::S, value::Real) where {T, S <: Union{MOI.LessThan{T},
                                            MOI.GreaterThan{T}, MOI.EqualTo{T}}}
    return S(convert(T, value))
end

"""
    JuMP.set_normalized_rhs(cref::InfOptConstraintRef, value::Real)::Nothing

Set the right-hand side term of `constraint` to `value`.
Note that prior to this step, JuMP will aggregate all constant terms onto the
right-hand side of the constraint. For example, given a constraint `2x + 1 <=
2`, `set_normalized_rhs(con, 4)` will create the constraint `2x <= 4`, not `2x +
1 <= 4`.

```julia-repl
julia> @constraint(model, con, 2x + 1 <= 2)
con : 2 x ≤ 1.0

julia> set_normalized_rhs(con, 4)

julia> con
con : 2 x ≤ 4.0
```
"""
function JuMP.set_normalized_rhs(cref::InfOptConstraintRef, value::Real)::Nothing
    old_constr = JuMP.constraint_object(cref)
    new_set = _set_set_value(JuMP.moi_set(old_constr), value)
    if old_constr isa BoundedScalarConstraint
        new_constr = BoundedScalarConstraint(JuMP.jump_function(old_constr),
                        new_set, parameter_bounds(old_constr),
                        original_parameter_bounds(old_constr))
    else
        new_constr = JuMP.ScalarConstraint(JuMP.jump_function(old_constr),
                                           new_set)
    end
    _set_core_constraint_object(cref, new_constr)
    return
end

"""
    JuMP.normalized_rhs(cref::InfOptConstraintRef)::Float64

Return the right-hand side term of `cref` after JuMP has converted the
constraint into its normalized form.
"""
function JuMP.normalized_rhs(cref::InfOptConstraintRef)::Float64
    constr = JuMP.constraint_object(cref)
    return MOI.constant(JuMP.moi_set(constr))
end

"""
    JuMP.add_to_function_constant(cref::InfOptConstraintRef, value::Real)::Nothing

Add `value` to the function constant term.
Note that for scalar constraints, JuMP will aggregate all constant terms onto the
right-hand side of the constraint so instead of modifying the function, the set
will be translated by `-value`. For example, given a constraint `2x <=
3`, `add_to_function_constant(c, 4)` will modify it to `2x <= -1`.
```
"""
function JuMP.add_to_function_constant(cref::InfOptConstraintRef,
                                       value::Real)::Nothing
    current_value = JuMP.normalized_rhs(cref)
    JuMP.set_normalized_rhs(cref, current_value - value)
    return
end

"""
    JuMP.set_normalized_coefficient(cref::InfOptConstraintRef,
                                    variable::GeneralVariableRef,
                                    value::Real)::Nothing

Set the coefficient of `variable` in the constraint `constraint` to `value`.
Note that prior to this step, JuMP will aggregate multiple terms containing the
same variable. For example, given a constraint `2x + 3x <= 2`,
`set_normalized_coefficient(con, x, 4)` will create the constraint `4x <= 2`.

```julia-repl
julia> con
con : 5 x ≤ 2.0

julia> set_normalized_coefficient(con, x, 4)

julia> con
con : 4 x ≤ 2.0
```
"""
function JuMP.set_normalized_coefficient(cref::InfOptConstraintRef,
                                         variable::GeneralVariableRef,
                                         value::Real)::Nothing
    # update the constraint expression and update the constraint
    old_constr = JuMP.constraint_object(cref)
    new_expr = _set_variable_coefficient!(JuMP.jump_function(old_constr),
                                          variable, value)
    if old_constr isa BoundedScalarConstraint
        new_constr = BoundedScalarConstraint(new_expr, JuMP.moi_set(old_constr),
                        parameter_bounds(old_constr),
                        original_parameter_bounds(old_constr))
    else
        new_constr = JuMP.ScalarConstraint(new_expr, JuMP.moi_set(old_constr))
    end
    _set_core_constraint_object(cref, new_constr)
    return
end

"""
    JuMP.normalized_coefficient(cref::InfOptConstraintRef,
                                variable::GeneralVariableRef)::Float64

Return the coefficient associated with `variable` in `constraint` after JuMP has
normalized the constraint into its standard form.
"""
function JuMP.normalized_coefficient(cref::InfOptConstraintRef,
                                     variable::GeneralVariableRef)::Float64
    constr = JuMP.constraint_object(cref)
    expr = JuMP.jump_function(constr)
    if expr isa GeneralVariableRef && expr == variable
        return 1.0
    elseif expr isa GeneralVariableRef
        return 0.0
    else
        return JuMP._affine_coefficient(expr, variable)
    end
end

# Return the appropriate constraint reference given the index and model
function _make_constraint_ref(model::InfiniteModel,
                              index::ConstraintIndex)::InfOptConstraintRef
    constr = model.constraints[index].constraint
    return InfOptConstraintRef(model, index, JuMP.shape(constr))
end

"""
    JuMP.constraint_by_name(model::InfiniteModel,
                            name::String)::Union{InfOptConstraintRef, Nothing}

Extend [`JuMP.constraint_by_name`](@ref JuMP.constraint_by_name)
to return the constraint reference
associated with `name` if one exists or returns nothing. Errors if more than
one constraint uses the same name.

**Example**
```julia-repl
julia> constraint_by_name(model, "constr_name")
constr_name : x + pt = 3.0
```
"""
function JuMP.constraint_by_name(model::InfiniteModel,
    name::String)::Union{InfOptConstraintRef, Nothing}
    if model.name_to_constr === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_constr = Dict{String, Int}()
        for (index, data_object) in model.constraints
            constr_name = data_object.name
            if haskey(model.name_to_constr, constr_name)
                model.name_to_constr[constr_name] = ConstraintIndex(-1)
            else
                model.name_to_constr[constr_name] = index
            end
        end
    end
    index = get(model.name_to_constr, name, nothing)
    if index === nothing
        return nothing
    elseif index == ConstraintIndex(-1)
        error("Multiple constraints have the name $name.")
    else
        return _make_constraint_ref(model, index)
    end
end

"""
    JuMP.num_constraints(model::InfiniteModel,
                         function_type::Type{<:JuMP.AbstractJuMPScalar},
                         set_type::Type{<:MOI.AbstractSet})::Int

Extend [`JuMP.num_constraints`](@ref JuMP.num_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to return the number of constraints
with a partiuclar function type and set type.

**Example**
```julia-repl
julia> num_constraints(model, HoldVariableRef, MOI.LessThan)
1
```
"""
function JuMP.num_constraints(model::InfiniteModel,
                              function_type::Type{<:JuMP.AbstractJuMPScalar},
                              set_type::Type{<:MOI.AbstractSet})::Int
    counter = 0
    for (index, data_object) in model.constraints
        if isa(JuMP.jump_function(data_object.constraint), function_type) &&
           isa(JuMP.moi_set(data_object.constraint), set_type)
            counter += 1
        end
    end
    return counter
end

"""
    JuMP.num_constraints(model::InfiniteModel,
                         function_type::Type{<:JuMP.AbstractJuMPScalar})::Int

Extend [`JuMP.num_constraints`](@ref JuMP.num_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to search by function types for all MOI
sets and return the total number of constraints with a particular function type.

```julia-repl
julia> num_constraints(model, HoldVariableRef)
3
```
"""
function JuMP.num_constraints(model::InfiniteModel,
                            function_type::Type{<:JuMP.AbstractJuMPScalar})::Int
    return JuMP.num_constraints(model, function_type, MOI.AbstractSet)
end

"""
    JuMP.num_constraints(model::InfiniteModel,
                         function_type::Type{<:MOI.AbstractSet})::Int

Extend [`JuMP.num_constraints`](@ref JuMP.num_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to search by MOI set type for all function
types and return the total number of constraints that use a particular MOI set
type.

```julia-repl
julia> num_constraints(model, MOI.LessThan)
2
```
"""
function JuMP.num_constraints(model::InfiniteModel,
                              set_type::Type{<:MOI.AbstractSet})::Int
    return JuMP.num_constraints(model, JuMP.AbstractJuMPScalar, set_type)
end

"""
    JuMP.num_constraints(model::InfiniteModel)::Int

Extend [`JuMP.num_constraints`](@ref JuMP.num_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to return the total number of constraints
in an infinite model `model`.

```julia-repl
julia> num_constraints(model)
4
```
"""
function JuMP.num_constraints(model::InfiniteModel)::Int
    return length(model.constraints)
end

"""
    JuMP.all_constraints(model::InfiniteModel,
                         function_type::Type{<:JuMP.AbstractJuMPScalar},
                         set_type::Type{<:MOI.AbstractSet}
                         )::Vector{InfOptConstraintRef}

Extend [`JuMP.all_constraints`](@ref JuMP.all_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to return a list of all the constraints with a particular function type and set type.

```julia-repl
julia> all_constraints(model, GeneralVariableRef, MOI.LessThan)
1-element Array{InfOptConstraintRef,1}:
 x ≤ 1.0
```
"""
function JuMP.all_constraints(model::InfiniteModel,
                              function_type::Type{<:JuMP.AbstractJuMPScalar},
                              set_type::Type{<:MOI.AbstractSet}
                              )::Vector{InfOptConstraintRef}
    constr_list = Vector{InfOptConstraintRef}(undef,
                           JuMP.num_constraints(model, function_type, set_type))
    counter = 1
    for (index, object) in model.constraints
        if isa(JuMP.jump_function(object.constraint), function_type) &&
           isa(JuMP.moi_set(object.constraint), set_type)
            constr_list[counter] = _make_constraint_ref(model, index)
            counter += 1
        end
    end
    return constr_list
end

"""
    JuMP.all_constraints(model::InfiniteModel,
                         function_type::Type{<:JuMP.AbstractJuMPScalar}
                         )::Vector{InfOptConstraintRef}

Extend [`JuMP.all_constraints`](@ref JuMP.all_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to search by function types for all MOI
sets and return a list of all constraints use a particular function type.

```julia-repl
julia> all_constraints(model, GeneralVariableRef)
3-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 x ≤ 3.0
 x integer
```
"""
function JuMP.all_constraints(model::InfiniteModel,
                              function_type::Type{<:JuMP.AbstractJuMPScalar}
                              )::Vector{InfOptConstraintRef}
    return JuMP.all_constraints(model, function_type, MOI.AbstractSet)
end

"""
    JuMP.all_constraints(model::InfiniteModel,
                         set_type::Type{<:MOI.AbstractSet}
                         )::Vector{InfOptConstraintRef}

Extend [`JuMP.all_constraints`](@ref JuMP.all_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to search by MOI set type for all function
types and return a list of all constraints that use a particular set type.

```julia-repl
julia> all_constraints(model, MOI.GreaterThan)
3-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 g(t) ≥ 0.0, ∀ t ∈ [0, 6]
 g(0.5) ≥ 0.0
```
"""
function JuMP.all_constraints(model::InfiniteModel,
                              set_type::Type{<:MOI.AbstractSet}
                              )::Vector{InfOptConstraintRef}
    return JuMP.all_constraints(model, JuMP.AbstractJuMPScalar, set_type)
end

"""
    JuMP.all_constraints(model::InfiniteModel)::Vector{InfOptConstraintRef}

Extend [`JuMP.all_constraints`](@ref JuMP.all_constraints(::JuMP.Model, ::Type{<:Union{JuMP.AbstractJuMPScalar, Vector{<:JuMP.AbstractJuMPScalar}}}, ::Type{<:MOI.AbstractSet}))
to return all a list of all the constraints
in an infinite model `model`.

```julia-repl
julia> all_constraints(model)
5-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 x ≤ 3.0
 x integer
 g(t) ≥ 0.0, ∀ t ∈ [0, 6]
 g(0.5) ≥ 0.0
```
"""
function JuMP.all_constraints(model::InfiniteModel)::Vector{InfOptConstraintRef}
    return [_make_constraint_ref(model, index)
            for (index, object) in model.constraints]
end

"""
    JuMP.list_of_constraint_types(model::InfiniteModel)::Vector{Tuple{DataType, DataType}}

Extend [`JuMP.list_of_constraint_types`](@ref JuMP.list_of_constraint_types(::JuMP.Model))
to return a list of tuples that
contain all the used combinations of function types and set types in the model.

```julia-repl
julia> all_constraints(model)
3-element Array{Tuple{DataType,DataType},1}:
 (GeneralVariableRef, MathOptInterface.LessThan{Float64})
 (GeneralVariableRef, MathOptInterface.GreaterThan{Float64})
 (GeneralVariableRef, MathOptInterface.Integer)
```
"""
function JuMP.list_of_constraint_types(model::InfiniteModel
    )::Vector{Tuple{DataType, DataType}}
    type_set = Set{Tuple{DataType, DataType}}()
    for (index, object) in model.constraints
        push!(type_set, (typeof(JuMP.jump_function(object.constraint)),
                          typeof(JuMP.moi_set(object.constraint))))
    end
    return collect(type_set)
end

################################################################################
#                         PARAMETER REFERENCE METHODS
################################################################################
"""
    parameter_refs(cref::InfOptConstraintRef)::Tuple

Return the tuple of infinite parameter references that determine the infinite
dependencies of `cref`.

**Example**
```julia-repl
julia> parameter_refs(cref)
(t,)
```
"""
function parameter_refs(cref::InfOptConstraintRef)::Tuple
    model = JuMP.owner_model(cref)
    obj_indices = _param_object_indices(model)[_object_numbers(cref)]
    return Tuple(_make_param_tuple_element(model, idx) for idx in obj_indices)
end

################################################################################
#                           PARAMETER BOUND METHODS
################################################################################
"""
    has_parameter_bounds(cref::InfOptConstraintRef)::Bool

Return a `Bool` indicating if `cref` is limited to a sub-domain as defined
by a [`ParameterBounds`](@ref) object.

**Example**
```julia-repl
julia> has_parameter_bounds(cref)
true
```
"""
function has_parameter_bounds(cref::InfOptConstraintRef)::Bool
    if JuMP.constraint_object(cref) isa BoundedScalarConstraint
        return !isempty(JuMP.constraint_object(cref).bounds)
    else
        return false
    end
end

"""
    parameter_bounds(cref::InfOptConstraintRef)::ParameterBounds{GeneralVariableRef}

Return the [`ParameterBounds`](@ref) object associated with the constraint
`cref`. Errors if `cref` does not have parameter bounds.

**Example**
```julia-repl
julia> parameter_bounds(cref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function parameter_bounds(cref::InfOptConstraintRef)::ParameterBounds{GeneralVariableRef}
    !has_parameter_bounds(cref) && error("$cref does not have parameter bounds.")
    return JuMP.constraint_object(cref).bounds
end

# Internal function used to change the parameter bounds of a constraint
function _update_constr_param_bounds(cref::InfOptConstraintRef,
                                     bounds::ParameterBounds{GeneralVariableRef},
                                     orig_bounds::ParameterBounds{GeneralVariableRef}
                                     )::Nothing
    old_constr = JuMP.constraint_object(cref)
    func = JuMP.jump_function(old_constr)
    set = JuMP.moi_set(old_constr)
    if !isempty(bounds)
        new_constr = BoundedScalarConstraint(func, set, bounds, orig_bounds)
    else
        new_constr = JuMP.ScalarConstraint(func, set)
    end
    _set_core_constraint_object(cref, new_constr)
    return
end

"""
    set_parameter_bounds(cref::InfOptConstraintRef,
                         bounds:ParameterBounds{GeneralVariableRef};
                         [force::Bool = false])::Nothing

Specify a new [`ParameterBounds`](@ref) object `bounds` for the constraint `cref`.
This is meant to be primarily used by [`@set_parameter_bounds`](@ref) which
provides a more intuitive syntax.

**Example**
```julia-repl
julia> set_parameter_bounds(cref, ParameterBounds((t => IntervalSet(0, 2),)))

julia> parameter_bounds(cref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function set_parameter_bounds(cref::InfOptConstraintRef,
                              bounds::ParameterBounds{GeneralVariableRef};
                              force::Bool = false, _error::Function = error
                              )::Nothing
    if has_parameter_bounds(cref) && !force
        _error("$cref already has parameter bounds. Consider adding more using " *
               "`add_parameter_bounds` or overwriting them by setting " *
               "the keyword argument `force = true`")
    else
        # check that bounds are valid and add support(s) if necessary
        _check_bounds(bounds, _error = _error)
        _validate_bounds(JuMP.owner_model(cref), bounds, _error = _error)
        # consider hold variables
        orig_bounds = copy(bounds)
        vrefs = _all_function_variables(JuMP.jump_function(JuMP.constraint_object(cref)))
        for vref in vrefs
            _update_var_bounds(vref, bounds)
        end
        # set the new bounds
        _update_constr_param_bounds(cref, bounds, orig_bounds)
        # update status
        set_optimizer_model_ready(JuMP.owner_model(cref), false)
    end
    return
end

"""
    add_parameter_bounds(cref::InfOptConstraintRef,
                         new_bounds::ParameterBounds{GeneralVariableRef}
                         )::Nothing

Add an additional parameter bound to `cref` such that it is defined over the
sub-domain based on `pref` from `lower` to `upper`. This is primarily meant to be
used by [`@add_parameter_bounds`](@ref).

```julia-repl
julia> add_parameter_bounds(cref, ParameterBounds((t => IntervalSet(0, 2),))

julia> parameter_bounds(cref)
Subdomain bounds (1): t ∈ [0, 2]
```
"""
function add_parameter_bounds(cref::InfOptConstraintRef,
                              new_bounds::ParameterBounds{GeneralVariableRef};
                              _error::Function = error)::Nothing
    # check the new bounds
    _check_bounds(new_bounds, _error = _error)
    _validate_bounds(JuMP.owner_model(cref), new_bounds, _error = _error)
    # add the bounds
    if JuMP.constraint_object(cref) isa BoundedScalarConstraint
        _update_bounds(parameter_bounds(cref), new_bounds, _error = _error)
        _update_bounds(original_parameter_bounds(JuMP.constraint_object(cref)),
                       new_bounds, _error = _error)
    else
        _update_constr_param_bounds(cref, new_bounds, copy(new_bounds))
    end
    # update the optimizer model status
    set_optimizer_model_ready(JuMP.owner_model(cref), false)
    return
end

"""
    delete_parameter_bounds(cref::InfOptConstraintRef)::Nothing

Delete all the parameter bounds of the constraint `cref`. Note any bounds that
are needed for hold variables inside in `cref` will be unaffected.

**Example**
```julia-repl
julia> @BDconstraint(model, c1(x == 0), y <= 42)
c1 : y(x) ≤ 42, ∀ x[1] = 0, x[2] = 0

julia> delete_parameter_bounds(c1)

julia> c1
c1 : y(x) ≤ 42, ∀ x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
```
"""
function delete_parameter_bounds(cref::InfOptConstraintRef)::Nothing
    # check if has bounds on pref from the constraint
    constr = JuMP.constraint_object(cref)
    if has_parameter_bounds(cref) && !isempty(original_parameter_bounds(constr))
        # consider hold variables
        new_bounds = ParameterBounds()
        vrefs = _all_function_variables(JuMP.jump_function(constr))
        for vref in vrefs
            _update_var_bounds(vref, new_bounds)
        end
        # set the new bounds
        _update_constr_param_bounds(cref, new_bounds, ParameterBounds())
    end
    return
end

################################################################################
#                                 DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel, cref::InfOptConstraintRef)::Nothing

Extend [`JuMP.delete`](@ref JuMP.delete(::JuMP.Model, ::JuMP.ConstraintRef{JuMP.Model}))
to delete an `InfiniteOpt` constraint and all associated information. Errors
if `cref` is invalid.

**Example**
```julia-repl
julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z ≥ 0.0
 g(t) + z ≥ 42.0, ∀ t ∈ [0, 6]

julia> delete(model, cref)

julia> print(model)
Min measure(g(t)*t) + z
Subject to
 z ≥ 0.0
```
"""
function JuMP.delete(model::InfiniteModel, cref::InfOptConstraintRef)::Nothing
    # check valid reference
    @assert JuMP.is_valid(model, cref) "Invalid constraint reference."
    # update variable dependencies
    constr = _core_constraint_object(cref)
    all_vrefs = _all_function_variables(JuMP.jump_function(constr))
    for vref in all_vrefs
        filter!(e -> e != JuMP.index(cref), _constraint_dependencies(vref))
    end
    # delete constraint information
    _delete_data_object(cref)
    # reset optimizer model status
    set_optimizer_model_ready(model, false)
    return
end
