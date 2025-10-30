################################################################################
#                             BASIC REFERENCE METHODS
################################################################################
"""
    JuMP.owner_model(cref::InfOptConstraintRef)::InfiniteModel

Extend `JuMP.owner_model` to return the infinite model associated with `cref`.

**Example**
```julia-repl
julia> model = owner_model(cref)
An InfiniteOpt Model
Minimization problem with:
  Finite parameters: 0
  Infinite parameter: 1
  Variables: 3
  Derivatives: 0
  Measures: 0
  Objective function type: GenericAffExpr{Float64, GeneralVariableRef}
  `GenericAffExpr{Float64, GeneralVariableRef}`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
  `GenericAffExpr{Float64, GeneralVariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
  `GeneralVariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 3 constraints
Names registered in the model: c1, c2, t, y, z
Transformation backend information:
  Backend type: TranscriptionBackend
  `t` transcribed over 10 supports
  Solver: Ipopt
  Transformation built and up-to-date: true
```
"""
JuMP.owner_model(cref::InfOptConstraintRef) = cref.model

"""
    JuMP.index(cref::InfOptConstraintRef)::InfOptConstraintIndex

Extend `JuMP.index` to return the index of an `InfiniteOpt` constraint `cref`.

**Example**
```julia-repl
julia> index(cref)
InfOptConstraintIndex(2)
```
"""
JuMP.index(cref::InfOptConstraintRef) = cref.index

# Extend Base and JuMP functions
function Base.:(==)(v::InfOptConstraintRef, w::InfOptConstraintRef)
    return v.model === w.model && v.index == w.index
end
Base.broadcastable(cref::InfOptConstraintRef) = Ref(cref)

################################################################################
#                             CORE OBJECT METHODS
################################################################################
# Extend _add_data_object
function _add_data_object(
    model::InfiniteModel,
    object::ConstraintData
    )
    return MOIUC.add_item(model.constraints, object)
end

# Extend _data_dictionary
function _data_dictionary(cref::InfOptConstraintRef)
    return JuMP.owner_model(cref).constraints
end

# Extend _data_object
function _data_object(cref::InfOptConstraintRef)
    object = get(_data_dictionary(cref), JuMP.index(cref), nothing)
    if isnothing(object)
        error("Invalid constraint reference, cannot find corresponding ", 
              "constraint in the model. This is likely caused by using the ",
              "reference of a deleted constraint.")
    end
    return object
end

"""
    JuMP.constraint_object(cref::InfOptConstraintRef)::JuMP.AbstractConstraint

Extend `JuMP.constraint_object` to return the constraint object associated with 
`cref`.

**Example**
```julia-repl
julia> @infinite_parameter(model, t in [0, 10]);

julia> @variable(model, x <= 1);

julia> cref = UpperBoundRef(x);

julia> obj = constraint_object(cref)
ScalarConstraint{GeneralVariableRef,MathOptInterface.LessThan{Float64}}(x,
MathOptInterface.LessThan{Float64}(1.0))
```
"""
function JuMP.constraint_object(cref::InfOptConstraintRef)
    return _data_object(cref).constraint
end

"""
    core_object(cref::InfOptConstraintRef)::JuMP.AbstractConstraint

Return the core underlying constraint object for `cref`.
This is intended for the developer API. For general usage,
[`JuMP.constraint_object`](@ref) should be used instead.
"""
function core_object(cref::InfOptConstraintRef)
    return JuMP.constraint_object(cref::InfOptConstraintRef)
end

## Set helper methods for adapting data_objects with parametric changes 
# No change needed 
function _adaptive_data_update(
    cref::InfOptConstraintRef, 
    c::C, 
    data::ConstraintData{C}
    ) where {C <: JuMP.AbstractConstraint}
    data.constraint = c
    return
end

# Reconstruction is necessary 
function _adaptive_data_update(
    cref::InfOptConstraintRef, 
    c::C1, 
    data::ConstraintData{C2}
    )  where {C1, C2}
    new_data = ConstraintData(c, data.group_int_idxs, data.name, 
                              data.measure_indices, data.is_info_constraint)
    _data_dictionary(cref)[JuMP.index(cref)] = new_data
    return
end

# Update the core constraint object
function _set_core_object(
    cref::InfOptConstraintRef,
    constr::JuMP.AbstractConstraint
    )
    _adaptive_data_update(cref, constr, _data_object(cref))
    set_transformation_backend_ready(JuMP.owner_model(cref), false)
    return
end

"""
    parameter_group_int_indices(cref::InfOptConstraintRef)::Vector{Int}

Return the list of infinite parameter group integer indices used by `cref`.
"""
function parameter_group_int_indices(cref::InfOptConstraintRef)
    return _data_object(cref).group_int_idxs
end

# Extend _measure_dependencies
function _measure_dependencies(cref::InfOptConstraintRef)
    return _data_object(cref).measure_indices
end

"""
    is_variable_domain_constraint(cref::InfOptConstraintRef)::Bool

Returns a `Bool` whether `cref` was created based on a variable's
domain. For instance, it could be the upper bound of a variable
`y(t)` which is normally queried via `UpperBoundRef`. This is
intended as a helper function for developers of new
transformation backends which typically ignore these constraints,
since they are taken care of when the variables are processed.
"""
function is_variable_domain_constraint(cref::InfOptConstraintRef)
    return _data_object(cref).is_info_constraint
end

# Extend _delete_data_object
function _delete_data_object(cref::InfOptConstraintRef)
    delete!(_data_dictionary(cref), JuMP.index(cref))
    return
end

################################################################################
#                             DEFINITION METHODS
################################################################################
"""
    DomainRestriction{F <: Function}

Tag for specifiying restricted domains for constraints. These are created via
```julia
DomainRestriction(restriction_func::Function, parameter_refs...)
```
where `restrict_func` is a function that accepts a support that follows the
formatting of `parameter_refs` and returns `Bool` on whether that support should
be included in the domain of the constraint.

**Example**
```julia
@infinite_parameter(model, t in [0, 1])
@infinite_parameter(model, x in [-1, 1])
@variable(model, y, Infinite(t, x))

restrict_func(t_s, x_s) = (0 <= t_s <= 0.5) && (x_s < 0)
restriction = DomainRestriction(restrict_func, t, x)
@constraint(model, y^2 >= 42, restriction)
```
"""
struct DomainRestriction{F}
    restriction_func::F
    parameter_refs::Collections.VectorTuple{GeneralVariableRef}
    function DomainRestriction(restrict_func::F, prefs...) where {F}
        return new{F}(restrict_func, Collections.VectorTuple(prefs))
    end
end

"""
    JuMP.build_constraint(
        _error::Function,
        func,
        set::MOI.AbstractSet,
        restriction::DomainRestriction
    )::DomainRestrictedConstraint

Extend `JuMP.buid_constraint` to handle including a `restriction` to its inherit 
infinite parameter domain in addition to the traditional `func` in `set` setup. 
This returns a `DomainRestrictedConstraint` that can then 
be added via `JuMP.add_constraint`. Errors if the restriction is not compatible 
with infinite parameter domains. 

**Example**
```julia
@infinite_parameter(model, t in [0, 1])
@infinite_parameter(model, x in [-1, 1])
@variable(model, y, Infinite(t, x))

restrict_func(t_s, x_s) = (0 <= t_s <= 0.5) && (x_s < 0)
restriction = DomainRestriction(restrict_func, t, x)
con = build_constraint(error, y + 2, MOI.LessThan(0.0), restriction);
```
"""
function JuMP.build_constraint(
    _error::Function,
    func,
    set::MOI.AbstractSet,
    restriction::DomainRestriction
    )
    # make the constraint and check the domain restrictions
    constr = JuMP.build_constraint(_error, func, set)
    pfunc_restrict = build_parameter_function(
        _error,
        restriction.restriction_func,
        restriction.parameter_refs
    )
    return DomainRestrictedConstraint(constr, pfunc_restrict)
end

# Used to update the model.var_to_constrs field
function _update_var_constr_mapping(
    vrefs::Vector{GeneralVariableRef},
    cref::InfOptConstraintRef
    )
    for vref in vrefs
        dvref = dispatch_variable_ref(vref)
        push!(_constraint_dependencies(dvref), JuMP.index(cref))
        if dvref isa MeasureRef
            push!(_measure_dependencies(cref), JuMP.index(dvref))
        end
        # TODO maybe make mapping for derivatives...
    end
    return
end

"""
    JuMP.add_constraint(
        model::InfiniteModel,
        c::JuMP.AbstractConstraint,
        [name::String = ""]
    )::InfOptConstraintRef

Extend `JuMP.add_constraint` to add a constraint `c` to an infinite model
`model` with name `name`. Returns an appropriate constraint reference whose type
depends on what variables are used to define the constraint. Errors if any 
variables do not belong to `model`. This is primarily used as an internal
method for the constraint macros.

**Example**
```julia-repl
julia> @infinite_parameter(model, t in [0, 10]);

julia> @variable(model, g, Infinite(t));

julia> @variable(model, x);

julia> constr = build_constraint(error, g + x, MOI.EqualTo(42));

julia> cref = add_constraint(model, constr, "name")
name : g(t) + x = 42.0, ∀ t ∈ [0, 10]
```
"""
function JuMP.add_constraint(
    model::InfiniteModel,
    c::JuMP.AbstractConstraint,
    name::String = "";
    is_info_constr::Bool = false
    )
    # gather the unique list of variable references for testing and mapping
    vrefs = all_expression_variables(JuMP.jump_function(c))
    # test in the model
    for vref in vrefs
        JuMP.check_belongs_to_model(vref, model)
    end
    # get the parameter group integer indices
    group_int_idxs = sort!(parameter_group_int_indices(vrefs))
    # add the constaint to the model
    constr_object = ConstraintData(
        c,
        group_int_idxs,
        name,
        MeasureIndex[],
        is_info_constr
        )
    cindex = _add_data_object(model, constr_object)
    cref = InfOptConstraintRef(model, cindex)
    # update the variable mappings and model status
    _update_var_constr_mapping(vrefs, cref)
    set_transformation_backend_ready(model, false)
    # clear out the name dictionary 
    model.name_to_constr = nothing
    # return the constraint reference
    return cref
end

################################################################################
#                                 JuMP METHODS
################################################################################
"""
    JuMP.is_valid(model::InfiniteModel, cref::InfOptConstraintRef)::Bool

Extend `JuMP.is_valid` to return `Bool` whether an `InfiniteOpt` constraint 
reference is valid.

**Example**
```julia-repl
julia> is_valid(model, cref)
true
```
"""
function JuMP.is_valid(
    model::InfiniteModel, 
    cref::InfOptConstraintRef
    )
    return (model === JuMP.owner_model(cref) &&
            JuMP.index(cref) in keys(_data_dictionary(cref)))
end

"""
    JuMP.name(cref::InfOptConstraintRef)::String

Extend `JuMP.name` to return the name of an `InfiniteOpt` constraint.

**Example**
```julia-repl
julia> name(cref)
"constr_name"
```
"""
function JuMP.name(cref::InfOptConstraintRef)
    object = get(_data_dictionary(cref), JuMP.index(cref), nothing)
    return isnothing(object) ? "" : object.name
end

"""
    JuMP.set_name(cref::InfOptConstraintRef, name::String)::Nothing

Extend `JuMP.set_name` to specify the name of a constraint `cref`.

**Example**
```julia-repl
julia> set_name(cref, "new_name")

julia> name(cref)
"new_name"
```
"""
function JuMP.set_name(
    cref::InfOptConstraintRef, 
    name::String
    )
    _data_object(cref).name = name
    JuMP.owner_model(cref).name_to_constr = nothing
    return
end

# Extend basic access methods for DomainRestrictedConstraint
function JuMP.jump_function(c::DomainRestrictedConstraint) 
    return JuMP.jump_function(c.constraint)
end
function JuMP.moi_set(c::DomainRestrictedConstraint)
    return JuMP.moi_set(c.constraint)
end
function _update_constraint(old_c::JuMP.ScalarConstraint, func, set)
    return JuMP.ScalarConstraint(func, set)
end
function _update_constraint(old_c::DomainRestrictedConstraint, func, set)
    return DomainRestrictedConstraint(
        JuMP.ScalarConstraint(func, set),
        old_c.restriction
    )
end

# Return a constraint set with an updated value
function _set_set_value(
    set::S, 
    value::Real
    ) where {T, S <: Union{MOI.LessThan{T}, MOI.GreaterThan{T}, MOI.EqualTo{T}}}
    return S(convert(T, value))
end

# Enforce that the MOI set is a traditional scalar one 
function _enforce_rhs_set(
    set::Union{MOI.LessThan{T}, MOI.GreaterThan{T}, MOI.EqualTo{T}}
    ) where {T}
    return 
end
function _enforce_rhs_set(set)
    error("The right hand side value is not well defined for constraints ", 
          "with sets of type `$(typeof(set))`.")
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
function JuMP.set_normalized_rhs(
    cref::InfOptConstraintRef, 
    value::Real
    )
    old_constr = JuMP.constraint_object(cref)
    set = JuMP.moi_set(old_constr)
    func = JuMP.jump_function(old_constr)
    _enforce_rhs_set(set)
    new_set = _set_set_value(set, value)
    new_constr = _update_constraint(old_constr, func, new_set)
    _set_core_object(cref, new_constr)
    return
end

"""
    JuMP.normalized_rhs(cref::InfOptConstraintRef)::Float64

Return the right-hand side term of `cref` after JuMP has converted the
constraint into its normalized form.
"""
function JuMP.normalized_rhs(cref::InfOptConstraintRef)
    constr = JuMP.constraint_object(cref)
    set = JuMP.moi_set(constr)
    _enforce_rhs_set(set)
    return MOI.constant(set)
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
function JuMP.add_to_function_constant(
    cref::InfOptConstraintRef,
    value::Real
    )
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
function JuMP.set_normalized_coefficient(
    cref::InfOptConstraintRef,
    variable::GeneralVariableRef,
    value::Real
    )
    # update the constraint expression and update the constraint
    old_constr = JuMP.constraint_object(cref)
    func = JuMP.jump_function(old_constr)
    set = JuMP.moi_set(old_constr)
    new_func = _set_variable_coefficient!(func, variable, value) # checks valid
    new_constr = _update_constraint(old_constr, new_func, set)
    _set_core_object(cref, new_constr)
    return
end

"""
    JuMP.normalized_coefficient(cref::InfOptConstraintRef,
                                variable::GeneralVariableRef)::Float64

Return the coefficient associated with `variable` in `constraint` after JuMP has
normalized the constraint into its standard form.
"""
function JuMP.normalized_coefficient(
    cref::InfOptConstraintRef,
    variable::GeneralVariableRef
    )
    constr = JuMP.constraint_object(cref)
    func = JuMP.jump_function(constr)
    return _affine_coefficient(func, variable) # checks valid
end

"""
    JuMP.constraint_by_name(
        model::InfiniteModel,
        name::String
    )::Union{InfOptConstraintRef, Nothing}

Extend `JuMP.constraint_by_name` to return the constraint reference
associated with `name` if one exists or returns nothing. Errors if more than
one constraint uses the same name.

**Example**
```julia-repl
julia> constraint_by_name(model, "constr_name")
constr_name : x + pt = 3.0
```
"""
function JuMP.constraint_by_name(
    model::InfiniteModel,
    name::String
    )
    if isnothing(model.name_to_constr)
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_constr = Dict{String, Int}()
        for (index, data_object) in model.constraints
            constr_name = data_object.name
            if haskey(model.name_to_constr, constr_name)
                model.name_to_constr[constr_name] = InfOptConstraintIndex(-1)
            else
                model.name_to_constr[constr_name] = index
            end
        end
    end
    index = get(model.name_to_constr, name, nothing)
    if isnothing(index)
        return nothing
    elseif index == InfOptConstraintIndex(-1)
        error("Multiple constraints have the name $name.")
    else
        return InfOptConstraintRef(model, index)
    end
end

"""
    JuMP.num_constraints(model::InfiniteModel, [function_type], [set_type])::Int

Extend `JuMP.num_constraints` to return the number of constraints with a 
partiuclar function type and set type.

**Example**
```julia-repl
julia> num_constraints(model, FiniteVariableRef, MOI.LessThan)
1

julia> num_constraints(model, FiniteVariableRef)
3

julia> num_constraints(model, MOI.LessThan)
2

julia> num_constraints(model)
4
```
"""
function JuMP.num_constraints(
    model::InfiniteModel,
    function_type,
    set_type
    )
    counter = 0
    for (index, data_object) in model.constraints
        if isa(JuMP.jump_function(data_object.constraint), function_type) &&
           isa(JuMP.moi_set(data_object.constraint), set_type)
            counter += 1
        end
    end
    return counter
end

# Function type only
function JuMP.num_constraints(
    model::InfiniteModel,
    function_type
    )
    return JuMP.num_constraints(model, function_type, MOI.AbstractSet)
end

# Set type only
function JuMP.num_constraints(
    model::InfiniteModel,
    set_type::Type{<:MOI.AbstractSet}
    )
    return JuMP.num_constraints(model, Any, set_type)
end

# All the constraints
JuMP.num_constraints(model::InfiniteModel) = length(model.constraints)

"""
    JuMP.all_constraints(
        model::InfiniteModel,
        [function_type], 
        [set_type]
    )::Vector{InfOptConstraintRef}

Extend `JuMP.all_constraints` to return a list of all the constraints with a 
particular function type and set type.

**Example**
```julia-repl
julia> all_constraints(model, GeneralVariableRef, MOI.LessThan)
1-element Array{InfOptConstraintRef,1}:
 x ≤ 1.0

julia> all_constraints(model, GeneralVariableRef)
3-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 x ≤ 3.0
 x integer

julia> all_constraints(model, MOI.GreaterThan)
3-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 g(t) ≥ 0.0, ∀ t ∈ [0, 6]
 g(0.5) ≥ 0.0

julia> all_constraints(model)
5-element Array{InfOptConstraintRef,1}:
 x ≥ 0.0
 x ≤ 3.0
 x integer
 g(t) ≥ 0.0, ∀ t ∈ [0, 6]
 g(0.5) ≥ 0.0
```
"""
function JuMP.all_constraints(
    model::InfiniteModel,
    function_type,
    set_type
    )
    num_constrs = JuMP.num_constraints(model, function_type, set_type)
    constr_list = Vector{InfOptConstraintRef}(undef, num_constrs)
    counter = 1
    for (index, object) in model.constraints
        if isa(JuMP.jump_function(object.constraint), function_type) &&
           isa(JuMP.moi_set(object.constraint), set_type)
            constr_list[counter] = InfOptConstraintRef(model, index)
            counter += 1
        end
    end
    return constr_list
end

# Function type only
function JuMP.all_constraints(model::InfiniteModel, function_type)
    return JuMP.all_constraints(model, function_type, MOI.AbstractSet)
end

# Set type only
function JuMP.all_constraints(
    model::InfiniteModel,
    set_type::Type{<:MOI.AbstractSet}
    )
    return JuMP.all_constraints(model, JuMP.AbstractJuMPScalar, set_type)
end

# All the constraints
function JuMP.all_constraints(model::InfiniteModel)
    return [InfOptConstraintRef(model, idx) for (idx, _) in model.constraints]
end

"""
    JuMP.list_of_constraint_types(model::InfiniteModel)::Vector{Tuple{DataType, DataType}}

Extend `JuMP.list_of_constraint_types` to return a list of tuples that contain 
all the used combinations of function types and set types in the model.

**Example**
```julia-repl
julia> all_constraints(model)
3-element Array{Tuple{DataType,DataType},1}:
 (GeneralVariableRef, MathOptInterface.LessThan{Float64})
 (GeneralVariableRef, MathOptInterface.GreaterThan{Float64})
 (GeneralVariableRef, MathOptInterface.Integer)
```
"""
function JuMP.list_of_constraint_types(model::InfiniteModel)
    type_set = Set{Tuple{DataType, DataType}}()
    for (_, object) in model.constraints
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
function parameter_refs(cref::InfOptConstraintRef)
    model = JuMP.owner_model(cref)
    prefs = parameter_refs(model)
    group_int_idxs = parameter_group_int_indices(cref)
    length(prefs) == length(group_int_idxs) && return prefs
    return prefs[group_int_idxs]
end

################################################################################
#                         DOMAIN RESTRICTION METHODS
################################################################################
# Check if domains have domain restriction
has_domain_restriction(::DomainRestrictedConstraint) = true
has_domain_restriction(::JuMP.AbstractConstraint) = false

"""
    has_domain_restriction(cref::InfOptConstraintRef)::Bool

Return a `Bool` indicating if `cref` is limited to a sub-domain as defined
by a [`DomainRestriction`](@ref) object.

**Example**
```julia-repl
julia> has_domain_restriction(cref)
true
```
"""
function has_domain_restriction(cref::InfOptConstraintRef)
    return has_domain_restriction(JuMP.constraint_object(cref))
end

"""
    domain_restriction(cref::InfOptConstraintRef)::ParameterFunction

Return the domain restriction formated as a [`ParameterFunction`](@ref) 
associated with the constraint`cref`. Errors if it does not have a
domain restriction.

**Example**
```julia-repl
julia> domain_restriction(cref)

```
"""
function domain_restriction(cref::InfOptConstraintRef)
    constr = JuMP.constraint_object(cref)
    if !has_domain_restriction(constr) 
        return error("`$cref` doesn't have a domain restriction.")
    end
    return constr.restriction
end

# jump constraint object getter
_jump_constraint(c::JuMP.AbstractConstraint) = c
_jump_constraint(c::DomainRestrictedConstraint) = c.constraint 

"""
    set_domain_restriction(
        cref::InfOptConstraintRef,
        restriction:DomainRestriction
    )::Nothing

Specify a new [`DomainRestriction`](@ref) object `restriction` for the 
constraint `cref`

**Example**
```julia-repl
julia> restrict_func(t_s) = 0 <= t_s <= 0.5;

julia> set_domain_restriction(cref, DomainRestrictions(restrict_func, t))

julia> domain_restriction(cref)

```
"""
function set_domain_restriction(
    cref::InfOptConstraintRef,
    restriction::DomainRestriction
    )
    pfunc = build_parameter_function(
        error,
        restriction.restriction_func,
        restriction.parameter_refs
    )
    constr = JuMP.constraint_object(cref)
    jump_constr = _jump_constraint(constr)
    _set_core_object(cref, DomainRestrictedConstraint(jump_constr, pfunc))
    set_transformation_backend_ready(JuMP.owner_model(cref), false)
    return
end

"""
    delete_domain_restriction(cref::InfOptConstraintRef)::Nothing

Delete the domain restriction of the constraint `cref`.

**Example**
```julia-repl
julia> delete_domain_restriction(c1)
```
"""
function delete_domain_restrictions(cref::InfOptConstraintRef)
    constr = JuMP.constraint_object(cref)
    if has_domain_restriction(constr)
        _set_core_object(cref, constr.constraint)
        set_transformation_backend_ready(JuMP.owner_model(cref), false)
    end
    return
end

################################################################################
#                                 DELETION
################################################################################
"""
    JuMP.delete(model::InfiniteModel, cref::InfOptConstraintRef)::Nothing

Extend `JuMP.delete` to delete an `InfiniteOpt` constraint and all associated 
information. Errors if `cref` is invalid.

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
function JuMP.delete(model::InfiniteModel, cref::InfOptConstraintRef)
    # check valid reference
    @assert JuMP.is_valid(model, cref) "Invalid constraint reference."
    # update variable dependencies
    constr = JuMP.constraint_object(cref)
    all_vrefs = all_expression_variables(JuMP.jump_function(constr))
    for vref in all_vrefs
        filter!(e -> e != JuMP.index(cref), _constraint_dependencies(vref))
    end
    # delete constraint information
    _delete_data_object(cref)
    # reset transformation backend status
    set_transformation_backend_ready(model, false)
    return
end
