function JuMP.show_backend_summary(io::IO, model::InfiniteModel) end

function JuMP.show_objective_function_summary(io::IO, model::InfiniteModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(model))
    return
end

function JuMP.objective_function_string(print_mode, model::InfiniteModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end

_plural(n) = (isone(n) ? "" : "s")

function JuMP.show_constraints_summary(io::IO, model::InfiniteModel)
    n = length(model.constrs)
    print(io, "Constraint", _plural(n), ": ", n)
    return
end

function JuMP.constraints_string(print_mode, model::InfiniteModel)
    strings = String[]
    # Sort by creation order
    constraints = sort(collect(model.constrs), by = c -> c.first)
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    # TODO print the infinite parameters with their sets
    # params = sort(collect(model.params), by = c -> c.first)
    # for (index, param) in params
    #     push!(strings, string(model.param_to_name[index], " ",
    #                           _math_symbol(print_mode, :in), " ",
    #                           )
    # end
    return strings
end

function Base.show(io::IO, ref::GeneralConstraintRef)
    print(io, JuMP.constraint_string(JuMP.REPLMode, ref))
    return
end

function Base.show(io::IO, ::MIME"text/latex", ref::GeneralConstraintRef)
    print(io, JuMP.constraint_string(JuMP.IJuliaMode, ref))
    return
end

function JuMP.constraint_string(print_mode, ref::GeneralConstraintRef)
    return JuMP.constraint_string(print_mode, JuMP.name(ref), JuMP.constraint_object(ref))
end
