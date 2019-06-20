function Base.show(io::IO, model::InfiniteModel)
    println(io, "An InfiniteOpt Model")
    sense = JuMP.objective_sense(model)
    if sense == MOI.MAX_SENSE
        print(io, "Maximization")
    elseif sense == MOI.MIN_SENSE
        print(io, "Minimization")
    else
        print(io, "Feasibility")
    end
    println(io, " problem with:")
    println(io, "Variable", _plural(JuMP.num_variables(model)), ": ",
            JuMP.num_variables(model))
    if sense != MOI.FEASIBILITY_SENSE
        JuMP.show_objective_function_summary(io, model)
    end
    JuMP.show_constraints_summary(io, model)
    JuMP.show_backend_summary(io, model)
    names_in_scope = sort(collect(keys(JuMP.object_dictionary(model))))
    if !isempty(names_in_scope)
        println(io)
        print(io, "Names registered in the model: ",
              join(string.(names_in_scope), ", "))
    end
end

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
    # TODO print multivariate distributions with variables grouped together
    params = sort(collect(model.params), by = c -> c.first)
    for (index, param) in params
        pref = ParameterRef(model, index)
        if is_used(pref)
            push!(strings, string(JuMP.function_string(print_mode, pref), " ",
                                  JuMP.in_set_string(print_mode, param.set)))
        end
    end
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
    return JuMP.constraint_string(print_mode, JuMP.name(ref),
           JuMP.constraint_object(ref))
end

function JuMP.in_set_string(print_mode, set::IntervalSet)
    return string(JuMP._math_symbol(print_mode, :in), " [", set.lower_bound,
                  ", ", set.upper_bound, "]")
end

function JuMP.in_set_string(print_mode, set::DistributionSet)
    d_string = string(set.distribution)
    first_bracket = findfirst(isequal('{'), d_string)
    last_bracket = findall(isequal('}'), d_string)[end]
    d_string = d_string[1:first_bracket-1] * d_string[last_bracket+1:end]
    new_lines = findall(isequal('\n'), d_string)
    for i = 1:length(new_lines)
        if new_lines[1] == length(d_string)
            d_string = d_string[1:end-1]
        elseif d_string[new_lines[1]-1] == '(' || d_string[new_lines[1]+1] == ')'
            d_string = d_string[1:new_lines[1]-1] * d_string[new_lines[1]+1:end]
        else
            d_string = d_string[1:new_lines[1]-1] * ", " * d_string[new_lines[1]+1:end]
        end
        new_lines = findall(isequal('\n'), d_string)
    end
    return string(JuMP._math_symbol(print_mode, :in), " ", d_string)
end
