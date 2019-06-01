const MyConstraintRef = JuMP.ConstraintRef{InfiniteModel, Int}
JuMP.constraint_type(::InfiniteModel) = MyConstraintRef
function JuMP.add_constraint(model::InfiniteModel, c::JuMP.AbstractConstraint,
                             name::String="")
    model.next_constr_index += 1
    index = model.next_constr_index
    cref = JuMP.ConstraintRef(model, index, JuMP.shape(c))
    model.constraints[index] = c
    JuMP.set_name(cref, name)
    return cref
end
function JuMP.delete(model::InfiniteModel, constraint_ref::MyConstraintRef)
    @assert JuMP.is_valid(model, constraint_ref)
    delete!(model.constraints, constraint_ref.index)
    delete!(model.constr_to_name, constraint_ref.index)
end
function JuMP.is_valid(model::InfiniteModel, constraint_ref::MyConstraintRef)
    return (model === constraint_ref.model &&
            constraint_ref.index in keys(model.constraints))
end
function JuMP.constraint_object(cref::MyConstraintRef)
    return cref.model.constraints[cref.index]
end

JuMP.name(cref::MyConstraintRef) = cref.model.constr_to_name[cref.index]
function JuMP.set_name(cref::MyConstraintRef, name::String)
    cref.model.constr_to_name[cref.index] = name
    cref.model.name_to_constr = nothing
end
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
    elseif index.value == -1
        error("Multiple constraints have the name $name.")
    else
        # We have no information on whether this is a vector constraint
        # or a scalar constraint
        return JuMP.ConstraintRef(model, index, JuMP.ScalarShape())
    end
end
