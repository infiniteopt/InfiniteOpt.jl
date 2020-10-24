################################################################################
#                            MATH SYMBOL METHODS
################################################################################
# Support additional math symbols beyond what JuMP does
function _infopt_math_symbol(::Type{JuMP.REPLMode}, name::Symbol)::String
    if name == :intersect
        return Sys.iswindows() ? "and" : "∩"
    elseif name == :prop
        return "~"
    elseif name == :partial 
        return Sys.iswindows() ? "d" : "∂"
    elseif name == :expect 
        return Sys.iswindows() ? "E" : "𝔼"
    elseif name == :integral 
        return Sys.iswindows() ? "integral" : "∫"
    else
        return JuMP._math_symbol(JuMP.REPLMode, name)
    end
end

# Support additional math symbols beyond what JuMP does
function _infopt_math_symbol(::Type{JuMP.IJuliaMode}, name::Symbol)::String
    if name == :intersect
        return "\\cap"
    elseif name == :prop
        return "\\sim"
    elseif name == :partial 
        return "\\partial"
    elseif name == :open_rng 
        return "\\left["
    elseif name == :close_rng 
        return "\\right]"
    elseif name == :expect 
        return "\\mathbb{E}"
    elseif name == :integral 
        return "\\int"
    else
        return JuMP._math_symbol(JuMP.IJuliaMode, name)
    end
end

################################################################################
#                            INFINITE SET METHODS
################################################################################
# Return "s" if n is greater than one
_plural(n) = (isone(n) ? "" : "s")

## Return the string of an infinite set
# IntervalSet
function set_string(print_mode, set::IntervalSet)::String
    return string("[", JuMP._string_round(JuMP.lower_bound(set)), ", ",
                  JuMP._string_round(JuMP.upper_bound(set)), "]")
end

# DistributionSet
function set_string(print_mode, set::DistributionSet)::String
    return string(set.distribution)
end

# CollectionSet
function set_string(print_mode, set::CollectionSet)::String
    sets = collection_sets(set)
    num_sets = length(sets)
    cs_string = string("CollectionSet with ", num_sets, " set",
                       _plural(num_sets), ":")
    for s in sets
        cs_string *= string("\n ", set_string(print_mode, s))
    end
    return cs_string
end

# Fallback
function set_string(print_mode, set::S) where {S <: AbstractInfiniteSet}
    # hack fix because just calling `string(set)` --> StackOverflow Error
    values = map(f -> getfield(set, f), fieldnames(S))
    return string(S, "(", join(values, ", "), ")")
end

## Return in set strings of infinite sets
# Extend to return of in set string for interval sets
function JuMP.in_set_string(print_mode, set::IntervalSet)::String
    if JuMP.lower_bound(set) != JuMP.upper_bound(set)
        return string(_infopt_math_symbol(print_mode, :in), " ",
                      set_string(print_mode, set))
    else
        return string(_infopt_math_symbol(print_mode, :eq), " ",
                      JuMP._string_round(JuMP.lower_bound(set)))
    end
end

# Extend to return of in set string for distribution sets
function JuMP.in_set_string(print_mode, set::DistributionSet)::String
    dist = set.distribution
    name = string(typeof(dist))
    bracket_index = findfirst(isequal('{'), name)
    if bracket_index !== nothing
        name = name[1:bracket_index-1]
    end
    if !isempty(size(dist))
        dims = size(dist)
        name *= string("(dim", _plural(length(dims)), ": (", join(dims, ", "), "))")
    end
    return string(_infopt_math_symbol(print_mode, :prop), " ", name)
end

# Extend to return of in set string of other sets
function JuMP.in_set_string(print_mode, set::AbstractInfiniteSet)::String
    return string(_infopt_math_symbol(print_mode, :in), " ",
                  set_string(print_mode, set))
end

## Extend in set string to consider parameter bounds
# IntervalSet
function JuMP.in_set_string(print_mode,
    pref::GeneralVariableRef,
    set::IntervalSet,
    bounds::ParameterBounds{GeneralVariableRef})::String
    # determine if in bounds
    in_bounds = haskey(bounds, pref)
    # make the string
    interval = in_bounds ? bounds[pref] : set
    return JuMP.in_set_string(print_mode, interval)
end

# InfiniteScalarSet
function JuMP.in_set_string(print_mode,
    pref::GeneralVariableRef,
    set::InfiniteScalarSet,
    bounds::ParameterBounds{GeneralVariableRef})::String
    # determine if in bounds
    if haskey(bounds, pref)
        bound_set = bounds[pref]
        if JuMP.lower_bound(bound_set) == JuMP.upper_bound(bound_set)
            return JuMP.in_set_string(print_mode, bound_set)
        else
            return  string(JuMP.in_set_string(print_mode, set), " ",
                           _infopt_math_symbol(print_mode, :intersect),
                           " ", set_string(print_mode, bound_set))
        end
    else
        return JuMP.in_set_string(print_mode, set)
    end
end

################################################################################
#                           MEASURE STRING METHODS
################################################################################
## Convert measure data into a useable string for measure printing
# 1-D DiscreteMeasureData/FunctionalDiscreteMeasureData
function measure_data_string(print_mode,
    data::Union{DiscreteMeasureData{GeneralVariableRef},
                FunctionalDiscreteMeasureData{GeneralVariableRef}}
    )::String
    pref = parameter_refs(data)
    lb = JuMP.lower_bound(data)
    ub = JuMP.upper_bound(data)
    nan_bound = isnan(lb) || isnan(ub)
    if nan_bound
        return JuMP.function_string(print_mode, pref)
    else
        set = IntervalSet(lb, ub)
        return string(JuMP.function_string(print_mode, pref), " ",
                      JuMP.in_set_string(print_mode, set))
    end
end

# Multi-D DiscreteMeasureData/FunctionalDiscreteMeasureData
function measure_data_string(print_mode,
    data::Union{DiscreteMeasureData{Vector{GeneralVariableRef}},
                FunctionalDiscreteMeasureData{Vector{GeneralVariableRef}}}
    )::String
    prefs = parameter_refs(data)
    lbs = JuMP.lower_bound(data)
    ubs = JuMP.upper_bound(data)
    has_bounds = !isnan(first(lbs)) && !isnan(first(ubs))
    homo_bounds = has_bounds && _allequal(lbs) && _allequal(ubs)
    names = map(p -> _remove_name_index(p), prefs)
    homo_names = _allequal(names)
    num_prefs = length(prefs)
    if homo_names && homo_bounds
        set = IntervalSet(first(lbs), first(ubs))
        return string(first(names), " ", JuMP.in_set_string(print_mode, set),
                      "^", num_prefs)
    elseif has_bounds
        str_list = [JuMP.function_string(print_mode, prefs[i]) * " " *
                    JuMP.in_set_string(print_mode, IntervalSet(lbs[i], ubs[i]))
                    for i in eachindex(prefs)]
        return _make_str_value(str_list)[2:end-1]
    elseif homo_names
        return first(names)
    else
        return _make_str_value(prefs)
    end
end

# extract the most compact parameter name possible
function _get_root_parameter_name(data::AbstractMeasureData)::String 
    prefs = parameter_refs(data)
    names = map(p -> _remove_name_index(p), prefs)
    if _allequal(names)
        return first(names)
    else
        return _make_str_value(prefs)
    end
end 

# Fallback for measure_data_string
function measure_data_string(print_mode, data::AbstractMeasureData)::String
    return _get_root_parameter_name(data)
end

# Make strings to represent measures in REPLMode
function variable_string(m::Type{JuMP.REPLMode}, mref::MeasureRef)::String
    data = measure_data(mref)
    data_str = measure_data_string(m, data)
    func_str = JuMP.function_string(m, measure_function(mref))
    name = JuMP.name(mref)
    if name == "integral"
        name = _infopt_math_symbol(m, :integral)
    elseif name == "expect"
        name = _infopt_math_symbol(m, :expect)
    end
    return string(name, "{", data_str, "}[", func_str, "]")
end

# Make strings to represent measures in IJuliaMode
function variable_string(m::Type{JuMP.IJuliaMode}, mref::MeasureRef)::String
    data = measure_data(mref)
    data_str = measure_data_string(m, data)
    func_str = JuMP.function_string(m, measure_function(mref))
    name = JuMP.name(mref) 
    if name == "integral" 
        param_name = _get_root_parameter_name(data)
        end_str = length(param_name) == 1 ? "d$param_name" : "d($param_name)"
        return string("\\int_{", data_str, "}", func_str, end_str)
    else 
        name = name == "expect" ? "\\mathbb{E}" : string("\\text{", name, "}")
        return string(name, "_{", data_str, "}", 
                  InfiniteOpt._infopt_math_symbol(m, :open_rng), func_str, 
                  InfiniteOpt._infopt_math_symbol(m, :close_rng))
    end
end

################################################################################
#                          VARIABLE STRING METHODS
################################################################################
## helper function for getting the variable names
# REPLMode
function _get_base_name(::Type{JuMP.REPLMode}, vref)::String
    var_name = JuMP.name(vref)
    if !isempty(var_name)
        return var_name
    else
        return "noname"
    end
end

# IJuliaMode
function _get_base_name(::Type{JuMP.IJuliaMode}, vref)::String
    var_name = JuMP.name(vref)
    if !isempty(var_name)
        # TODO: This is wrong if variable name constains extra "]"
        return replace(replace(var_name, "[" => "_{", count = 1), "]" => "}")
    else
        return "noname"
    end
end

# Helper method for infinite variable construction 
function _add_on_parameter_refs(base_name::String, prefs::VectorTuple)::String 
    param_name_tuple = "("
    for i in 1:size(prefs, 1)
        element_prefs = prefs[i, :]
        type = _index_type(first(element_prefs))
        if type == DependentParameterIndex
            param_name = _remove_name_index(first(element_prefs))
        elseif length(element_prefs) == 1
            param_name = JuMP.name(first(element_prefs))
        else
            # TODO this isn't quite right with a subset of an independent container
            names = map(p -> _remove_name_index(p), element_prefs)
            if _allequal(names)
                param_name = first(names)
            else
                param_name = string("[", join(element_prefs, ", "), "]")
            end
        end
        if i != size(prefs, 1)
            param_name_tuple *= string(param_name, ", ")
        else
            param_name_tuple *= string(param_name, ")")
        end
    end
    return string(base_name, param_name_tuple)
end

# Make a string for InfiniteVariableRef
function variable_string(print_mode, vref::InfiniteVariableRef)::String
    base_name = _get_base_name(print_mode, vref)
    if !haskey(_data_dictionary(vref), JuMP.index(vref))
        return base_name
    else
        prefs = raw_parameter_refs(vref)
        return _add_on_parameter_refs(base_name, prefs)
    end
end

## Make helper function for making derivative operators 
# REPL 
function _deriv_operator(m::Type{JuMP.REPLMode}, pref)::String
    return string(_infopt_math_symbol(m, :partial), "/", 
                  _infopt_math_symbol(m, :partial), variable_string(m, pref))
end

# IJulia 
function _deriv_operator(m::Type{JuMP.IJuliaMode}, pref)::String
    return string("\\frac{", _infopt_math_symbol(m, :partial), "}{", 
                  _infopt_math_symbol(m, :partial), " ", 
                  variable_string(m, pref), "}")
end

# TODO implement more intelligent naming for nested derivatives (i.e., use exponents)
# TODO account for container naming when variable macro is used (maybe deal with this at the macro end)
# Make a string for DerivativeRef 
function variable_string(print_mode, dref::DerivativeRef)::String
    base_name = _get_base_name(print_mode, dref)
    if !haskey(_data_dictionary(dref), JuMP.index(dref))
        return base_name
    elseif base_name != "noname"
        prefs = raw_parameter_refs(dref)
        return _add_on_parameter_refs(base_name, prefs)
    else
        vref = dispatch_variable_ref(derivative_argument(dref))
        pref = operator_parameter(dref)
        return string(_deriv_operator(print_mode, pref), 
                      _infopt_math_symbol(print_mode, :open_rng), 
                      variable_string(print_mode, vref), 
                      _infopt_math_symbol(print_mode, :close_rng))
    end
end

## Return the parameter value as an appropriate string
# Number
function _make_str_value(value)::String
    return JuMP._string_round(value)
end

# Array{<:Number}
function _make_str_value(values::Array)::String
    if length(values) == 1
        return _make_str_value(first(values))
    end
    if length(values) <= 4
        str_value = "["
        for i in eachindex(values)
            if i != length(values)
                str_value *= JuMP._string_round(values[i]) * ", "
            else
                str_value *= JuMP._string_round(values[i]) * "]"
            end
        end
        return str_value
    else
        return string("[", JuMP._string_round(first(values)), ", ..., ",
                      JuMP._string_round(last(values)), "]")
    end
end

# Make a string for PointVariableRef
# TODO improve so numerator of derivative contains the point
function variable_string(print_mode, vref::PointVariableRef)::String
    if !haskey(_data_dictionary(vref), JuMP.index(vref)) || !isempty(JuMP.name(vref))
        return _get_base_name(print_mode, vref)
    else
        ivref = dispatch_variable_ref(infinite_variable_ref(vref))
        if ivref isa InfiniteVariableRef || !isempty(JuMP.name(ivref))
            base_name = _get_base_name(print_mode, ivref)
        else 
            base_name = variable_string(print_mode, ivref) # we have a derivative
        end
        prefs = raw_parameter_refs(ivref)
        values = raw_parameter_values(vref)
        name = string(base_name, "(")
        for i in 1:size(prefs, 1)
            if i != size(prefs, 1)
                name *= string(_make_str_value(values[prefs.ranges[i]]), ", ")
            else
                name *= string(_make_str_value(values[prefs.ranges[i]]), ")")
            end
        end
        return name
    end
end

# Make a string for ReducedVariableRef
# TODO improve so numerator of derivative contains the parameter tuple
function variable_string(print_mode, vref::ReducedVariableRef)::String
    if !haskey(_data_dictionary(vref), JuMP.index(vref)) || !isempty(JuMP.name(vref))
        return _get_base_name(print_mode, vref)
    else
        ivref = dispatch_variable_ref(infinite_variable_ref(vref))
        if ivref isa InfiniteVariableRef || !isempty(JuMP.name(ivref))
            base_name = _get_base_name(print_mode, ivref)
        else 
            base_name = variable_string(print_mode, ivref) # we have a derivative
        end
        prefs = raw_parameter_refs(ivref)
        eval_supps = eval_supports(vref)
        raw_list = [i in keys(eval_supps) ? eval_supps[i] : prefs[i]
                    for i in eachindex(prefs)]
        param_name_tuple = "("
        for i in 1:size(prefs, 1)
            value = raw_list[prefs.ranges[i]]
            if i != size(prefs, 1)
                param_name_tuple *= string(_make_str_value(value), ", ")
            else
                param_name_tuple *= string(_make_str_value(value), ")")
            end
        end
        return string(base_name, param_name_tuple)
    end
end

# Fallback
function variable_string(print_mode, vref::JuMP.AbstractVariableRef)::String
    return _get_base_name(print_mode, vref)
end

# Extend function string for DispatchVariableRefs (REPL)
function JuMP.function_string(::Type{JuMP.REPLMode},
                              vref::DispatchVariableRef)::String
    return variable_string(JuMP.REPLMode, vref)
end

# Extend function string for DispatchVariableRefs (IJulia)
function JuMP.function_string(::Type{JuMP.IJuliaMode},
                              vref::DispatchVariableRef)::String
    return variable_string(JuMP.IJuliaMode, vref)
end

# Extend function string for GeneralVariableRefs (REPL)
function JuMP.function_string(::Type{JuMP.REPLMode},
                              vref::GeneralVariableRef)::String
    return variable_string(JuMP.REPLMode, dispatch_variable_ref(vref))
end

# Extend function string for GeneralVariableRefs (IJulia)
function JuMP.function_string(::Type{JuMP.IJuliaMode},
                              vref::GeneralVariableRef)::String
    return variable_string(JuMP.IJuliaMode, dispatch_variable_ref(vref))
end

################################################################################
#                         CONSTRAINT STRING METHODS
################################################################################
# Return parameter bound list as a string
function bound_string(print_mode, bounds::ParameterBounds{GeneralVariableRef}
                      )::String
    string_list = ""
    for (pref, set) in bounds
        string_list *= string(JuMP.function_string(print_mode, pref), " ",
                              JuMP.in_set_string(print_mode, set), ", ")
    end
    return string_list[1:end-2]
end

## Return the parameter domain string given the object index
# IndependentParameter
function _param_domain_string(print_mode, model::InfiniteModel,
                              index::IndependentParameterIndex,
                              bounds::ParameterBounds{GeneralVariableRef}
                              )::String
    pref = dispatch_variable_ref(model, index)
    set = infinite_set(pref)
    gvref = GeneralVariableRef(model, MOIUC.key_to_index(index), typeof(index))
    return string(JuMP.function_string(print_mode, pref), " ",
                  JuMP.in_set_string(print_mode, gvref, set, bounds))
end

# DependentParameters
function _param_domain_string(print_mode, model::InfiniteModel,
                              index::DependentParametersIndex,
                              bounds::ParameterBounds{GeneralVariableRef}
                              )::String
    # parse the infinite set
    first_gvref = GeneralVariableRef(model, MOIUC.key_to_index(index),
                                     DependentParameterIndex, 1)
    set = _parameter_set(dispatch_variable_ref(first_gvref))
    # iterate over the individual scalar sets if is a collection set
    if set isa CollectionSet
        domain_str = ""
        for i in eachindex(collection_sets(set))
            sset = collection_sets(set)[i]
            gvref = GeneralVariableRef(model, MOIUC.key_to_index(index),
                                       DependentParameterIndex, i)
            domain_str *= string(JuMP.function_string(print_mode, gvref), " ",
                        JuMP.in_set_string(print_mode, gvref, sset, bounds), ", ")
        end
        domain_str = domain_str[1:end-2]
    else
        # determine if bounds contain equalities and filter to the bounds of interest
        is_eq = haskey(bounds, first_gvref) && JuMP.lower_bound(bounds[first_gvref]) == JuMP.upper_bound(bounds[first_gvref])
        filtered_bounds = filter(e -> e[1].raw_index == MOIUC.key_to_index(index) &&
                                 e[1].index_type == DependentParameterIndex, bounds)
        # build the domain string
        if is_eq
            domain_str = bound_string(print_mode, filtered_bounds)
        else
            domain_str = string(_remove_name_index(first_gvref), " ",
                            JuMP.in_set_string(print_mode, set))
            if !isempty(filtered_bounds)
                domain_str *= string(" ", _infopt_math_symbol(print_mode, :intersect),
                                     " (", bound_string(print_mode, filtered_bounds), ")")
            end
        end
    end
    return domain_str
end

# Return string of an infinite constraint
function JuMP.constraint_string(print_mode, cref::InfOptConstraintRef;
                                in_math_mode = false)::String
    # get the function and set strings
    func_str = JuMP.function_string(print_mode, _core_constraint_object(cref))
    in_set_str = JuMP.in_set_string(print_mode, _core_constraint_object(cref))
    # check if constraint if finite
    obj_nums = _object_numbers(cref)
    if isempty(obj_nums)
        bound_str = ""
    else
        # get the parameter bounds if there are any
        if has_parameter_bounds(cref)
            bounds = parameter_bounds(cref)
        else
            bounds = ParameterBounds()
        end
        # prepare the parameter domains
        model = JuMP.owner_model(cref)
        bound_str = string(", ", JuMP._math_symbol(print_mode, :for_all), " ")
        for index in _param_object_indices(model)[_object_numbers(cref)]
            bound_str *= string(_param_domain_string(print_mode, model, index, bounds),
                                ", ")
        end
        bound_str = bound_str[1:end-2]
    end
    # form the constraint string
    if print_mode == JuMP.REPLMode
        lines = split(func_str, '\n')
        lines[1 + div(length(lines), 2)] *= " " * in_set_str * bound_str
        constr_str = join(lines, '\n')
    else
        constr_str = string(func_str, " ", in_set_str, bound_str)
    end
    # format for IJulia
    if print_mode == JuMP.IJuliaMode && !in_math_mode
        constr_str = JuMP._wrap_in_inline_math_mode(constr_str)
    end
    # add name if it has one
    name = JuMP.name(cref)
    if isempty(name) || (print_mode == JuMP.IJuliaMode && in_math_mode)
        return constr_str
    else
        return string(name, " : ", constr_str)
    end
end

# return list of constraint
function JuMP.constraints_string(print_mode, model::InfiniteModel)::Vector{String}
    # allocate the string vector
    strings = Vector{String}(undef, JuMP.num_constraints(model))
    # produce a string for each constraint
    counter = 1
    for (index, object) in model.constraints
        cref = _make_constraint_ref(model, index)
        strings[counter] = JuMP.constraint_string(print_mode, cref,
                                                  in_math_mode = true)
        counter += 1
    end
    return strings
end

################################################################################
#                           OBJECTIVE STRING METHODS
################################################################################
# Return the objective string
function JuMP.objective_function_string(print_mode, model::InfiniteModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end

################################################################################
#                           DATATYPE SHOWING METHODS
################################################################################
# Show infinite sets in REPLMode
function Base.show(io::IO, set::AbstractInfiniteSet)
    print(io, set_string(JuMP.REPLMode, set))
    return
end

# Show infinite sets in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", set::AbstractInfiniteSet)
    print(io, set_string(JuMP.IJuliaMode, set))
end

# TODO use ... when necessary --> Subdomain bounds: t in [0, 1], x[1] == 0, x[2] == 0, ... x[6] == 0, a in [2, 3]
# Show ParameterBounds in REPLMode
function Base.show(io::IO, bounds::ParameterBounds)
    print(io, "Subdomain bounds (", length(bounds), "): ",
          bound_string(JuMP.REPLMode, bounds))
    return
end

# Show ParameterBounds in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", bounds::ParameterBounds)
    print(io, "Subdomain bounds (", length(bounds), "): ",
          bound_string(JuMP.IJuliaMode, bounds))
end

# Show constraint in REPLMode
function Base.show(io::IO, ref::InfOptConstraintRef)
    print(io, JuMP.constraint_string(JuMP.REPLMode, ref))
    return
end

# Show constraint in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", ref::InfOptConstraintRef)
    print(io, JuMP.constraint_string(JuMP.IJuliaMode, ref))
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
    # show finite parameter info
    num_finite_params = num_parameters(model, FiniteParameter)
    println(io, "Finite Parameter", _plural(num_finite_params), ": ",
            num_finite_params)
    # show infinite parameter info
    num_infinite_params = num_parameters(model, InfiniteParameter)
    println(io, "Infinite Parameter", _plural(num_infinite_params), ": ",
            num_infinite_params)
    # show variable info
    num_vars = JuMP.num_variables(model)
    println(io, "Variable", _plural(num_vars), ": ", num_vars)
    # show the derivative info 
    num_derivs = num_derivatives(model)
    println(io, "Derivative", _plural(num_derivs), ": ", num_derivs)
    # show measure info
    num_meas = num_measures(model)
    println(io, "Measure", _plural(num_meas), ": ", num_meas)
    # show objective function info
    if sense != MOI.FEASIBILITY_SENSE
        JuMP.show_objective_function_summary(io, model)
    end
    # show constraint info
    JuMP.show_constraints_summary(io, model)
    # show other info
    names_in_scope = sort!(collect(keys(JuMP.object_dictionary(model))))
    if !isempty(names_in_scope)
        print(io, "Names registered in the model: ",
              join(string.(names_in_scope), ", "), "\n")
    end
    JuMP.show_backend_summary(io, model)
end
