"""
GeneralConstraintRef
An abstract type to define new variable types.
"""
abstract type GeneralConstraintRef end

"""
InfiniteConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain infinite variables
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct InfiniteConstraintRef <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::JuMP.AbstractShape
end

"""
FiniteConstraintRef <: GeneralConstraintRef
A DataType for constraints that contain finite variables
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of constraint in model.
- `shape::JuMP.AbstractShape` Shape of constraint
"""
struct FiniteConstraintRef <: GeneralConstraintRef
    model::InfiniteModel
    index::Int
    shape::JuMP.AbstractShape
end

"""
    owner_model(cref::GeneralConstraintRef)
Extends the method `JuMP.owner_model` for our new constraints.
"""
owner_model(cref::GeneralConstraintRef) = cref.model

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
    if c.func isa InfiniteAffExpr || c.func isa InfiniteQuadExpr
        cref = InfiniteConstraintRef(model, index, JuMP.shape(c))
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
    delete!(model.constrs, cref.index)
    delete!(model.constr_to_name, cref.index)
end

"""
    JuMP.is_valid(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.is_valid` function to accomodate our new constraint types.
"""
function JuMP.is_valid(model::InfiniteModel, cref::GeneralConstraintRef)
    return (model === cref.model &&
            cref.index in keys(model.constrs))
end

"""
    JuMP.constraint_object(model::InfiniteModel, cref::GeneralConstraintRef)
Extend the `JuMP.constraint_object` function to accomodate our new constraint types.
"""
function JuMP.constraint_object(cref::GeneralConstraintRef)
    return cref.model.constrs[cref.index]
end

"""
    JuMP.name(cref::GeneralConstraintRef)
Extend the `JuMP.name` function to accomodate our new constraint types.
"""
JuMP.name(cref::GeneralConstraintRef) = cref.model.constr_to_name[cref.index]

"""
    JuMP.set_name(cref::GeneralConstraintRef, name::String)
Extend the `JuMP.set_name` function to accomodate our new constraint types.
"""
function JuMP.set_name(cref::GeneralConstraintRef, name::String)
    cref.model.constr_to_name[cref.index] = name
    cref.model.name_to_constr = nothing
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
        if model.constrs[index].func isa InfiniteAffExpr || model.constrs[index].func isa InfiniteQuadExpr
            return InfiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
        else
            return FiniteConstraintRef(model, index, JuMP.shape(model.constrs[index]))
        end
    end
end
