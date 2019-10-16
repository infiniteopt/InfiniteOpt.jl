# Extend to return of in set string for interval sets
function JuMP.in_set_string(print_mode, set::IntervalSet)::String
    if set.lower_bound != set.upper_bound
        return string(JuMP._math_symbol(print_mode, :in), " [",
                      JuMP._string_round(set.lower_bound), ", ",
                      JuMP._string_round(set.upper_bound), "]")
    else
        return string(JuMP._math_symbol(print_mode, :eq), " ",
                      JuMP._string_round(set.lower_bound))
    end
end

# Extend to return of in set string for distribution sets
function JuMP.in_set_string(print_mode, set::DistributionSet)::String
    # get distribution string
    d_string = string(set.distribution)
    # remove number type
    first_bracket = findfirst(isequal('{'), d_string)
    if first_bracket != nothing
        last_bracket = findall(isequal('}'), d_string)[end]
        d_string = d_string[1:first_bracket-1] * d_string[last_bracket+1:end]
    end
    # remove new line characters
    new_lines = findall(isequal('\n'), d_string)
    for i = 1:length(new_lines)
        if new_lines[1] == length(d_string)
            d_string = d_string[1:end-1]
        elseif d_string[new_lines[1]-1] == '(' || d_string[new_lines[1]+1] == ')'
            d_string = d_string[1:new_lines[1]-1] * d_string[new_lines[1]+1:end]
        else
            d_string = d_string[1:new_lines[1]-1] * ", " *
                       d_string[new_lines[1]+1:end]
        end
        new_lines = findall(isequal('\n'), d_string)
    end
    return string(JuMP._math_symbol(print_mode, :in), " ", d_string)
end

# Return contraint parameter bounds as a string
function bound_string(print_mode, bounds::Dict)::String
    string_list = JuMP._math_symbol(print_mode, :for_all) * " "
    for (pref, set) in bounds
        string_list *= string(JuMP.function_string(print_mode, pref), " ",
                              JuMP.in_set_string(print_mode, set), ", ")
    end
    return string_list[1:end-2]
end

# Return string of a bounded constraint
function JuMP.constraint_string(print_mode,
                                constraint_object::BoundedScalarConstraint)
    # get needed strings
    func_str = JuMP.function_string(print_mode, constraint_object)
    in_set_str = JuMP.in_set_string(print_mode, constraint_object)
    bound_str = bound_string(print_mode, constraint_object.bounds)
    # modify as needed
    if print_mode == JuMP.REPLMode
        lines = split(func_str, '\n')
        lines[1 + div(length(lines), 2)] *= " " * in_set_str * ", " * bound_str
        return join(lines, '\n')
    else
        return func_str * " " * in_set_str * ", " * bound_str
    end
end

# Return a string for a constraint reference
function JuMP.constraint_string(print_mode, ref::GeneralConstraintRef)
    return JuMP.constraint_string(print_mode, JuMP.name(ref),
           JuMP.constraint_object(ref))
end

# return list of constraint
function JuMP.constraints_string(print_mode, model::InfiniteModel)
    strings = String[]
    # Sort by creation order
    constraints = sort(collect(model.constrs), by = c -> c.first)
    # produce a string for each constraint
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    # produce the parameters that are used along with their sets
    params = sort(collect(model.params), by = c -> c.first)
    ignore_groups = [] # used to indicate which groups should no longer be used
    for (index, param) in params
        pref = ParameterRef(model, index)
        if is_used(pref) && !(group_id(pref) in ignore_groups)
            if isa(infinite_set(pref), DistributionSet{<:Distributions.MultivariateDistribution})
                # print multivariate distributions with variables grouped together
                push!(ignore_groups, group_id(pref))
                push!(strings, string(_root_name(pref), " ",
                                JuMP.in_set_string(print_mode, param.set)))
            else
                # print parameter individually
                push!(strings, string(JuMP.function_string(print_mode, pref), " ",
                                      JuMP.in_set_string(print_mode, param.set)))
            end
        end
    end
    return strings
end

# Return the objective string
function JuMP.objective_function_string(print_mode, model::InfiniteModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end

# TODO show hold variables better

# Show constraint in REPLMode
function Base.show(io::IO, ref::GeneralConstraintRef)
    print(io, JuMP.constraint_string(JuMP.REPLMode, ref))
    return
end

# Show constraint in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", ref::GeneralConstraintRef)
    print(io, JuMP.constraint_string(JuMP.IJuliaMode, ref))
    return
end

# Show the backend information associated with the optimizer model
function JuMP.show_backend_summary(io::IO, model::InfiniteModel)
    println(io, "Optimizer model backend information: ")
    JuMP.show_backend_summary(io, optimizer_model(model))
    return
end

# Show the objective function type
function JuMP.show_objective_function_summary(io::IO, model::InfiniteModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(model))
    return
end

# Return "s" if n is greater than one
_plural(n) = (isone(n) ? "" : "s")

# Return constraint summary in JuMP like manner
function JuMP.show_constraints_summary(io::IO, model::InfiniteModel)
    for (F, S) in JuMP.list_of_constraint_types(model)
        n_constraints = JuMP.num_constraints(model, F, S)
        println(io, "`$F`-in-`$S`: $n_constraints constraint",
                _plural(n_constraints))
    end
    return
end

# Have the infinitemodel object show similar to Model
function Base.show(io::IO, model::InfiniteModel)
    # indicate is an infinite model
    println(io, "An InfiniteOpt Model")
    # show objective ssnse info
    sense = JuMP.objective_sense(model)
    if sense == MOI.MAX_SENSE
        print(io, "Maximization")
    elseif sense == MOI.MIN_SENSE
        print(io, "Minimization")
    else
        print(io, "Feasibility")
    end
    println(io, " problem with:")
    # show variable info
    println(io, "Variable", _plural(JuMP.num_variables(model)), ": ",
            JuMP.num_variables(model))
    # show obejctive function info
    if sense != MOI.FEASIBILITY_SENSE
        JuMP.show_objective_function_summary(io, model)
    end
    # show constraint info
    JuMP.show_constraints_summary(io, model)
    # show other info
    names_in_scope = sort(collect(keys(JuMP.object_dictionary(model))))
    if !isempty(names_in_scope)
        print(io, "Names registered in the model: ",
              join(string.(names_in_scope), ", "), "\n")
    end
    JuMP.show_backend_summary(io, model)
end
