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
Finite Parameters: 0
Infinite Parameters: 3
Variables: 3
Derivatives: 0
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
    object = Base.get(_data_dictionary(cref), JuMP.index(cref), nothing)
    if isnothing(object)
        error("Invalid constraint reference, cannot find corresponding ", 
              "constraint in the model. This is likely caused by using the ",
              "reference of a deleted constraint.")
    end
    return object
end

# Return the core constraint object
function _core_constraint_object(cref::InfOptConstraintRef)
    return _data_object(cref).constraint
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
    new_data = ConstraintData(c, data.object_nums, data.name, 
                              data.measure_indices, data.is_info_constraint)
    _data_dictionary(cref)[JuMP.index(cref)] = new_data
    return
end

# Update the core constraint object
function _set_core_constraint_object(
    cref::InfOptConstraintRef,
    constr::JuMP.AbstractConstraint
    )
    _adaptive_data_update(cref, constr, _data_object(cref))
    set_optimizer_model_ready(JuMP.owner_model(cref), false)
    return
end

# Extend _object_numbers
function _object_numbers(cref::InfOptConstraintRef)
    return _data_object(cref).object_nums
end

# Extend _measure_dependencies
function _measure_dependencies(cref::InfOptConstraintRef)
    return _data_object(cref).measure_indices
end

# Return if this constraint is an info constraint
function _is_info_constraint(cref::InfOptConstraintRef)
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
# Check that domain_restrictions argument is valid
function _check_restrictions(
    restrictions::DomainRestrictions{GeneralVariableRef};
    _error = error
    )
    depend_counter = Dict{DependentParameterRef, Int}()
    for (pref, domain) in restrictions
        # check that pref is an infinite parameter
        if !(_index_type(pref) <: InfiniteParameterIndex)
            _error("Can only specify infinite parameters for domain restrictions.")
        end
        # check that respects lower bound
        if JuMP.has_lower_bound(pref) && (JuMP.lower_bound(domain) < JuMP.lower_bound(pref))
                _error("Specified parameter lower bound exceeds that defined " *
                       "for $pref.")
        end
        # check that respects upper bound
        if JuMP.has_upper_bound(pref) && (JuMP.upper_bound(domain) > JuMP.upper_bound(pref))
                _error("Specified parameter upper bound exceeds that defined " *
                       "for $pref.")
        end
        # keep track of dependent parameters with equality restrictions to ensure completeness
        if (_index_type(pref) == DependentParameterIndex) &&
           (JuMP.lower_bound(domain) == JuMP.upper_bound(domain))
            idx = DependentParameterIndex(DependentParametersIndex(_raw_index(pref)), 1)
            dumby_pref = dispatch_variable_ref(pref.model, idx)
            if haskey(depend_counter, dumby_pref)
                depend_counter[dumby_pref] += 1
            else
                depend_counter[dumby_pref] = 1
            end
        end
    end
    # check dimensions of dependent parameter equalities (if any are provided)
    for (pref, amt) in depend_counter
        if _num_parameters(pref) != amt
            _error("Cannot specify equality domain restrictions for a subset of ",
                   "dependent infinite parameters.")
        end
    end
    return
end

"""
    JuMP.build_constraint(_error::Function, func, set,
                          restrictions::DomainRestrictions{GeneralVariableRef}
                          )::DomainRestrictedConstraint

Extend `JuMP.buid_constraint` to handle including `restrictions` to its inherit 
infinite parameter domains in addition to the traditional `func` in `set` setup. 
This returns a `DomainRestrictedConstraint` that can then 
be added via `JuMP.add_constraint`. Errors if the restrictions are incompadible 
with infinite parameter domains. 

**Example**
```julia-repl
julia> restrictions = DomainRestrictions(t => 0)
Subdomain restrictions (1): t = 0

julia> con = build_constraint(error, y + 2, MOI.LessThan(0.0), restrictions);
```
"""
function JuMP.build_constraint(
    _error::Function,
    func,
    set::MOI.AbstractSet,
    restrictions::DomainRestrictions
    )
    # make the constraint and check the domain restrictions
    constr = JuMP.build_constraint(_error, func, set)
    _check_restrictions(restrictions, _error = _error)
    return DomainRestrictedConstraint(constr, restrictions)
end

# Validate domain restrictions and add support(s) if needed
function _validate_restrictions(
    model::InfiniteModel,
    restrictions::DomainRestrictions
    )
    depend_supps = Dict{DependentParametersIndex, Matrix{Float64}}()
    for (pref, domain) in restrictions
        # check validity
        JuMP.check_belongs_to_model(pref, model)
        # ensure has a support if a point constraint was given
        if (JuMP.lower_bound(domain) == JuMP.upper_bound(domain))
            if _index_type(pref) == IndependentParameterIndex
                # label will be UserDefined
                add_supports(pref, JuMP.lower_bound(domain), check = false)
            else
                idx = DependentParametersIndex(_raw_index(pref))
                if !haskey(depend_supps, idx)
                    dumby_pref = dispatch_variable_ref(model, DependentParameterIndex(idx, 1))
                    depend_supps[idx] = Matrix{Float64}(undef, _num_parameters(dumby_pref), 1)
                end
                depend_supps[idx][_param_index(pref)] = JuMP.lower_bound(domain)
            end
        end
    end
    # add dependent supports if any are given
    for (idx, supp) in depend_supps
        prefs = [dispatch_variable_ref(model, DependentParameterIndex(idx, i))
                 for i in 1:length(supp)]
        add_supports(prefs, supp, check = false) # label will be UserDefined
    end
    return
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
    JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint,
                        [name::String = ""])::InfOptConstraintRef

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
    vrefs = _all_function_variables(JuMP.jump_function(c))
    # test in the model
    for vref in vrefs
        JuMP.check_belongs_to_model(vref, model)
    end
    # get the parameter object numbers
    object_nums = sort!(_object_numbers(vrefs))
    # add the constaint to the model
    constr_object = ConstraintData(c, object_nums, name, MeasureIndex[],
                                   is_info_constr)
    cindex = _add_data_object(model, constr_object)
    cref = InfOptConstraintRef(model, cindex)
    # update the variable mappings and model status
    _update_var_constr_mapping(vrefs, cref)
    set_optimizer_model_ready(model, false)
    # clear out the name dictionary 
    model.name_to_constr = nothing
    # return the constraint reference
    return cref
end

# Domain restricted constraints
function JuMP.add_constraint(
    model::InfiniteModel,
    c::DomainRestrictedConstraint,
    name::String = ""
    )
    # test domain restrictions and add needed supports
    _validate_restrictions(model, c.restrictions)
    # add the underlying constraint 
    cref = JuMP.add_constraint(model, c.constraint, name)
    # add the domain restrictions 
    model.constraint_restrictions[JuMP.index(cref)] = c.restrictions
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
function JuMP.constraint_object(
    cref::InfOptConstraintRef
    )
    return _core_constraint_object(cref)
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
    object = Base.get(_data_dictionary(cref), JuMP.index(cref), nothing)
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
    _enforce_rhs_set(JuMP.moi_set(old_constr))
    new_set = _set_set_value(JuMP.moi_set(old_constr), value)
    new_constr = JuMP.ScalarConstraint(JuMP.jump_function(old_constr), new_set)
    _set_core_constraint_object(cref, new_constr)
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
    new_func = _set_variable_coefficient!(JuMP.jump_function(old_constr),
                                          variable, value) # checks valid
    new_constr = JuMP.ScalarConstraint(new_func, JuMP.moi_set(old_constr))
    _set_core_constraint_object(cref, new_constr)
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

# Return the appropriate constraint reference given the index and model
function _make_constraint_ref(model::InfiniteModel,
    index::InfOptConstraintIndex
    )
    return InfOptConstraintRef(model, index)
end

"""
    JuMP.constraint_by_name(model::InfiniteModel,
                            name::String)::Union{InfOptConstraintRef, Nothing}

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
    index = Base.get(model.name_to_constr, name, nothing)
    if isnothing(index)
        return nothing
    elseif index == InfOptConstraintIndex(-1)
        error("Multiple constraints have the name $name.")
    else
        return _make_constraint_ref(model, index)
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
function JuMP.num_constraints(model::InfiniteModel)
    return length(model.constraints)
end

"""
    JuMP.all_constraints(model::InfiniteModel, [function_type], [set_type])::Vector{InfOptConstraintRef}

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
    return [_make_constraint_ref(model, idx) for (idx, _) in model.constraints]
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
    obj_indices = _param_object_indices(model)[_object_numbers(cref)]
    return Tuple(_make_param_tuple_element(model, idx) for idx in obj_indices)
end

################################################################################
#                         DOMAIN RESTRICTION METHODS
################################################################################
"""
    has_domain_restrictions(cref::InfOptConstraintRef)::Bool

Return a `Bool` indicating if `cref` is limited to a sub-domain as defined
by a [`DomainRestrictions`](@ref) object.

**Example**
```julia-repl
julia> has_domain_restrictions(cref)
true
```
"""
function has_domain_restrictions(cref::InfOptConstraintRef)
    return !isempty(domain_restrictions(cref))
end

"""
    domain_restrictions(cref::InfOptConstraintRef)::DomainRestrictions{GeneralVariableRef}

Return the [`DomainRestrictions`](@ref) object associated with the constraint
`cref`.

**Example**
```julia-repl
julia> domain_restrictions(cref)
Subdomain restrictions (1): t ∈ [0, 2]
```
"""
function domain_restrictions(cref::InfOptConstraintRef)
    return Base.get(JuMP.owner_model(cref).constraint_restrictions, JuMP.index(cref), 
               DomainRestrictions())
end

"""
    set_domain_restrictions(cref::InfOptConstraintRef,
                         restrictions:DomainRestrictions{GeneralVariableRef};
                         [force::Bool = false])::Nothing

Specify a new [`DomainRestrictions`](@ref) object `restrictions` for the 
constraint `cref`. Errors if `cref` already has restrictions and `force = false`. 
Where possible it is recommended to use [`add_domain_restrictions`](@ref) instead.  

**Example**
```julia-repl
julia> set_domain_restrictions(cref, DomainRestrictions(t => [0, 2]))

julia> domain_restrictions(cref)
Subdomain restrictions (1): t ∈ [0, 2]
```
"""
function set_domain_restrictions(
    cref::InfOptConstraintRef,
    restrictions::DomainRestrictions{GeneralVariableRef};
    force::Bool = false
    )
    if has_domain_restrictions(cref) && !force
        error("$cref already has domain restrictions. Consider adding more using " *
               "`add_domain_restrictions` or overwriting them by setting " *
               "the keyword argument `force = true`")
    end
    # check that restrictions are valid and add support(s) if necessary
    _check_restrictions(restrictions)
    model = JuMP.owner_model(cref)
    _validate_restrictions(model, restrictions)
    # set the new restrictions
    model.constraint_restrictions[JuMP.index(cref)] = restrictions
    # update status
    set_optimizer_model_ready(JuMP.owner_model(cref), false)
    return
end

# Update the current restrictions to overlap with the new ones if possible
function _update_restrictions(
    old::DomainRestrictions{GeneralVariableRef},
    new::DomainRestrictions{GeneralVariableRef}
    )
    # check each new restriction
    for (pref, domain) in new
        # we have a new restriction
        if !haskey(old, pref)
            old[pref] = domain
        # the previous domain and the new one do not overlap
        elseif (JuMP.lower_bound(domain) > JuMP.upper_bound(old[pref])) ||
               (JuMP.upper_bound(domain) < JuMP.lower_bound(old[pref]))
            error("The new domain restrictions are incompatible with the ", 
                  "existing ones.")
        # we have an existing restriction
        else
            # we have a new stricter lower bound to update with
            if JuMP.lower_bound(domain) > JuMP.lower_bound(old[pref])
                old[pref] = IntervalDomain(JuMP.lower_bound(domain),
                                            JuMP.upper_bound(old[pref]))
            end
            # we have a new stricter upper bound to update with
            if JuMP.upper_bound(domain) < JuMP.upper_bound(old[pref])
                old[pref] = IntervalDomain(JuMP.lower_bound(old[pref]),
                                            JuMP.upper_bound(domain))
            end
        end
    end
    return
end

"""
    add_domain_restrictions(cref::InfOptConstraintRef,
                         new_restrictions::DomainRestrictions{GeneralVariableRef}
                         )::Nothing

Add additional domain restrictions to `cref` such that it is defined over the
sub-domain based on `pref` from `lower` to `upper`.

```julia-repl
julia> add_domain_restrictions(cref, DomainRestrictions(t => [0, 2]))

julia> domain_restrictions(cref)
Subdomain restrictions (1): t ∈ [0, 2]
```
"""
function add_domain_restrictions(
    cref::InfOptConstraintRef,
    new_restrictions::DomainRestrictions{GeneralVariableRef}
    )
    # check the new restrictions
    _check_restrictions(new_restrictions)
    model = JuMP.owner_model(cref)
    _validate_restrictions(model, new_restrictions)
    # add the restrictions
    if has_domain_restrictions(cref)
        _update_restrictions(domain_restrictions(cref), new_restrictions)
    else 
        model.constraint_restrictions[JuMP.index(cref)] = new_restrictions
    end
    # update the optimizer model status
    set_optimizer_model_ready(JuMP.owner_model(cref), false)
    return
end

"""
    delete_domain_restrictions(cref::InfOptConstraintRef)::Nothing

Delete all the domain restrictions of the constraint `cref`. Note any restrictions that
are needed for finite variables inside in `cref` will be unaffected.

**Example**
```julia-repl
julia> @constraint(model, c1, y <= 42, DomainRestrictions(x => 0))
c1 : y(x) ≤ 42, ∀ x[1] = 0, x[2] = 0

julia> delete_domain_restrictions(c1)

julia> c1
c1 : y(x) ≤ 42, ∀ x[1] ∈ [-1, 1], x[2] ∈ [-1, 1]
```
"""
function delete_domain_restrictions(cref::InfOptConstraintRef)
    # delete the restrictions if there are any
    delete!(JuMP.owner_model(cref).constraint_restrictions, JuMP.index(cref))
    # update status
    set_optimizer_model_ready(JuMP.owner_model(cref), false)
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
    constr = _core_constraint_object(cref)
    all_vrefs = _all_function_variables(JuMP.jump_function(constr))
    for vref in all_vrefs
        filter!(e -> e != JuMP.index(cref), _constraint_dependencies(vref))
    end
    # delete any restrictions there are 
    delete_domain_restrictions(cref)
    # delete constraint information
    _delete_data_object(cref)
    # reset optimizer model status
    set_optimizer_model_ready(model, false)
    return
end
