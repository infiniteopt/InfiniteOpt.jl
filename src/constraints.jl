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

"""
    JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint, name::String="")
Extend the `JuMP.add_constraint` function to accomodate the `InfiniteModel` object.
"""
function JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint, name::String="")
    model.next_constr_index += 1
    index = model.next_constr_index
    if c.func isa InfiniteExpr
        cref = InfiniteConstraintRef(model, index, JuMP.shape(c))
    elseif c.func isa MeasureExpr
        cref = MeasureConstraintRef(model, index, JuMP.shape(c))
    else
        cref = FiniteConstraintRef(model, index, JuMP.shape(c))
    end
    model.constrs[index] = c
    JuMP.set_name(cref, name)
    return cref
end

"""
    JuMP.delete(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.delete` function to accomodate our new constraint types.
"""
function JuMP.delete(model::InfiniteModel, cref::GeneralConstraintRef)
    @assert JuMP.is_valid(model, cref)
    delete!(model.constrs, JuMP.index(cref))
    delete!(model.constr_to_name, JuMP.index(cref))
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
        if model.constrs[index].func isa InfiniteExpr
            return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
        elseif model.constrs[index].func isa MeasureExpr
            return MeasureConstraintRef(model, index, JuMP.shape(model.constrs[index]))
        else
            return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
        end
    end
end

# TODO Implement and extend the num_constraint, all_constraint, and list_of_constraint_types methods
