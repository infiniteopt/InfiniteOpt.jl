"""
    JuMP.owner_model(cref::GeneralConstraintRef)
Extends the method `JuMP.owner_model` for our new constraints.
"""
JuMP.owner_model(cref::GeneralConstraintRef) = cref.model

"""
    JuMP.index(cref::GeneralConstraintRef)::Int
Extent `JuMP.index` to return the index of a InfOpt constraint.
"""
JuMP.index(cref::GeneralConstraintRef) = cref.index

# Extend Base and JuMP functions
Base.:(==)(v::GeneralConstraintRef, w::GeneralConstraintRef) = v.model === w.model && v.index == w.index && v.shape == w.shape
Base.broadcastable(cref::GeneralConstraintRef) = Ref(cref)
JuMP.constraint_type(m::InfiniteModel) = GeneralConstraintRef

# This might not be necessary...
function JuMP.build_constraint(_error::Function, v::InfiniteVariableRef,
                               set::MOI.AbstractScalarSet;
                               parameter_bounds::Dict{ParameterRef, IntervalSet} = Dict{ParameterRef, IntervalSet}())
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
Extent `JuMP.build_constraint` to accept the parameter_bounds argument.
"""
function JuMP.build_constraint(_error::Function, expr::InfiniteExpr,
                               set::MOI.AbstractScalarSet;
                               parameter_bounds::Dict{ParameterRef, IntervalSet} = Dict{ParameterRef, IntervalSet}())
    offset = JuMP.constant(expr)
    JuMP.add_to_expression!(expr, -offset)
    if length(parameter_bounds) != 0
        return BoundedScalarConstraint(expr, MOIU.shift_constant(set, -offset), parameter_bounds)
    else
        return JuMP.ScalarConstraint(expr, MOIU.shift_constant(set, -offset))
    end
end

# Used to update the model.var_to_constrs field
function _update_var_constr_mapping(vrefs::Vector{<:GeneralVariableRef}, cindex::Int)
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
        end
    end
    return
end

# Extend functions for bounded constraints
JuMP.shape(c::BoundedScalarConstraint) = JuMP.shape(JuMP.ScalarConstraint(c.func, c.set))
JuMP.jump_function(c::BoundedScalarConstraint) = c.func
JuMP.moi_set(c::BoundedScalarConstraint) = c.set

"""
    JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint, name::String="")
Extend the `JuMP.add_constraint` function to accomodate the `InfiniteModel` object.
"""
function JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint, name::String="")
    isa(c, JuMP.VectorConstraint) && error("Vector constraints not supported.")
    vrefs = _all_function_variables(c.func)
    isa(vrefs, Vector{ParameterRef}) && error("Constraints cannot contain only parameters.")
    for vref in vrefs
        JuMP.owner_model(vref) != model && error("Variable $vref does not belong to model.")
    end
    # TODO add checks for bounded constraints, like check prefs are valid
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
        # else
        #     JuMP.delete(model, vref) --> Might not want to delete measures
        end
    end
    delete!(model.constrs, JuMP.index(cref))
    delete!(model.constr_to_name, JuMP.index(cref))
    delete!(model.constr_in_var_info, JuMP.index(cref))
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
    JuMP.name(cref::GeneralConstraintRef)
Extend the `JuMP.name` function to accomodate our new constraint types.
"""
JuMP.name(cref::GeneralConstraintRef) = JuMP.owner_model(cref).constr_to_name[JuMP.index(cref)]

"""
    JuMP.set_name(cref::GeneralConstraintRef, name::String)
Extend the `JuMP.set_name` function to accomodate our new constraint types.
"""
function JuMP.set_name(cref::GeneralConstraintRef, name::String)
    JuMP.owner_model(cref).constr_to_name[JuMP.index(cref)] = name
    JuMP.owner_model(cref).name_to_constr = nothing
    return
end

# Return the appropriate constraint reference given the index and model
function _make_constraint_ref(model::InfiniteModel, index::Int)
    if model.constrs[index].func isa InfiniteExpr
        return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    elseif model.constrs[index].func isa MeasureExpr
        return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    else
        return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
    end
end

"""
    JuMP.constraint_by_name(model::InfiniteModel, name::String)
Extend the `JuMP.constraint_by_name` function to accomodate our new constraint types.
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
