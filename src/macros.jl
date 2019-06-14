using Base.Meta
using JuMP: _valid_model, _error_if_cannot_register, object_dictionary, variable_type

function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:<=}, Val{:≤}},
    upper)
    JuMP._set_upper_bound_or_error(_error, infoexpr, upper)
end

function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:>=}, Val{:≥}},
    lower)
    JuMP._set_lower_bound_or_error(_error, infoexpr, lower)
end

function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:in}, Val{:∈}}, value)
    _dist_or_error(_error, infoexpr, value)
end

function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Val{S}, value) where S
    _error("Unknown sense $S.")
end

# In that case we assume the variable is the lhs.
function _parse_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                        sense::Symbol, var, value)
    _parse_one_operator_parameter(_error, infoexpr, Val(sense),
                                JuMP._esc_non_constant(value))
    return var
end

# If the lhs is a number and not the rhs, we can deduce that the rhs is
# the variable.
function _parse_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                        sense::Symbol, value::Number, var)
    _parse_one_operator_parameter(_error, infoexpr, JuMP.reverse_sense(Val(sense)),
                                JuMP._esc_non_constant(value))
    return var
end

function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                 ::Union{Val{:<=}, Val{:≤}}, lower,
                                 ::Union{Val{:<=}, Val{:≤}}, upper)
    JuMP._set_lower_bound_or_error(_error, infoexpr, lower)
    JuMP._set_upper_bound_or_error(_error, infoexpr, upper)
end

function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                 ::Union{Val{:>=}, Val{:≥}}, upper,
                                 ::Union{Val{:>=}, Val{:≥}}, lower)
    _parse_ternary_parameter(_error, infoexpr, Val(:≤), lower, Val(:≤), upper)
end

function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                 ::Val, lvalue, ::Val, rvalue)
    _error("Use the form lb <= ... <= ub.")
end

function _parse_parameter(_error::Function, infoexpr::_ParameterInfoExpr, lvalue,
                        lsign::Symbol, var, rsign::Symbol, rvalue)
    # lvalue lsign var rsign rvalue
    _parse_ternary_parameter(_error, infoexpr, Val(lsign),
                           JuMP._esc_non_constant(lvalue), Val(rsign),
                           JuMP._esc_non_constant(rvalue))
    return var
end

# TODO Check for dimensionality with multivariate distribution
"""
    @infinite_parameter(model, args)
A macro for the defining parameters of type `ParameterRef`.
"""
macro infinite_parameter(model, args...)
    _error(str...) = JuMP._macro_error(:infinite_parameter, (model, args...), str...)

    esc_model = esc(model)

    extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(extra) == 0
        x = gensym()
        anon_singleton = true
    else
        x = popfirst!(extra)
        anon_singleton = false
    end

    info_kw_args = filter(_is_param_keyword, kw_args)
    extra_kw_args = filter(kw -> kw.args[1] != :base_name && !InfOpt._is_param_keyword(kw), kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    infoexpr = InfOpt._ParameterInfoExpr(; JuMP._keywordify.(info_kw_args)...)

    # There are four cases to consider:
    # x                                         | type of x | x.head
    # ------------------------------------------+-----------+------------
    # param                                     | Symbol    | NA
    # param[1:2]                                | Expr      | :ref
    # param <= ub or var[1:2] <= ub             | Expr      | :call
    # lb <= param <= ub or lb <= var[1:2] <= ub | Expr      | :comparison
    # In the two last cases, we call parse_variable
    explicit_comparison = isexpr(x, :comparison) || isexpr(x, :call)
    if explicit_comparison
        param = InfOpt._parse_parameter(_error, infoexpr, x.args...)
    else
        param = x
    end

    anonvar = isexpr(param, :vect) || isexpr(param, :vcat) || anon_singleton
    anonvar && explicit_comparison && _error("Cannot use explicit bounds via >=, <= with an anonymous parameter")
    parameter = gensym()
    name = JuMP._get_name(param)
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        _error("Expression $name should not be used as a parameter name. Use the \"anonymous\" syntax $name = @infinite_parameter(model, ...) instead.")
    end

    set = InfOpt._constructor_set(_error, infoexpr)
    if isa(param, Symbol)
        # Easy case - a single variable
        buildcall = :( build_parameter($_error, $set, $(extra...)) )
        JuMP._add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $base_name) )
        # The looped code is trivial here since there is a single variable
        creationcode = :($parameter = $parametercall)
        final_parameter = parameter
    else
        isa(param, Expr) || _error("Expected $param to be a parameter name")
        # We now build the code to generate the variables (and possibly the
        # SparseAxisArray to contain them)
        refcall, idxparams, idxsets, condition = JuMP._build_ref_sets(param, parameter)
        clear_dependencies(i) = (JuMP.Containers.is_dependent(idxparams, idxsets[i], i) ? () : idxsets[i])

        # Code to be used to create each variable of the container.
        buildcall = :( build_parameter($_error, $set, $(extra...)) )
        JuMP._add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $(JuMP._name_call(base_name, idxparams))) )
        code = :( $(refcall) = $parametercall )
        # Determine the return type of add_variable. This is needed to create the container holding them.
        vartype = :( variable_type($esc_model, Parameter) )
        creationcode = JuMP._get_looped_code(parameter, code, condition, idxparams, idxsets, vartype, requestedcontainer)
        final_parameter = parameter
    end
    if anonvar
        # Anonymous variable, no need to register it in the model-level
        # dictionary nor to assign it to a variable in the user scope.
        # We simply return the variable
        macro_code = JuMP._macro_return(creationcode, final_parameter)
    else
        # We register the variable reference to its name and
        # we assign it to a variable in the local scope of this name
        macro_code = JuMP._macro_assign_and_return(creationcode, parameter, name,
                                              final_variable=final_parameter,
                                              model_for_registering = esc_model)
    end
    return JuMP._assert_valid_model(esc_model, macro_code)
end

# TODO Enable expression parsing of the form var(params)
# """
#     @infinite_variable(model, param_refs, args...)
# A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
# it defines variables of type `InfiniteVariableRef`.
# """
# macro infinite_variable(model, param_refs, args...)
#     # TODO properly implement error messages
#     code = quote
#         @assert isa($model, InfiniteModel)
#         if $param_refs isa Tuple
#             InfOpt._check_parameter_tuple($param_refs)
#             InfOpt._check_tuple_names($param_refs)
#         else
#             @assert typeof($param_refs) <: Union{ParameterRef, AbstractArray{<:ParameterRef}}
#             @assert InfOpt._only_one_name($param_refs)
#         end
#         JuMP.@variable($model, ($(args...)), variable_type = Infinite, param_refs = $param_refs)
#     end
#     return esc(code)
# end

# Assume variable on lhs
function _less_than_parse(arg1, arg2)
    if isexpr(arg1, :call)
        return Expr(:call, :<=, arg1.args[1], arg2), Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :<=, arg1, arg2), nothing
    end
end

# Assume variable on lhs
function _greater_than_parse(arg1, arg2)
    if isexpr(arg1, :call)
        return Expr(:call, :>=, arg1.args[1], arg2), Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :>=, arg1, arg2), nothing
    end
end

# The variable is on the rhs
function _less_than_parse(arg1::Number, arg2)
    return _greater_than_parse(arg2, arg1)
end

# The variable is on the rhs
function _greater_than_parse(arg1::Number, arg2)
    return _less_than_parse(arg2, arg1)
end

# Assume variable on lhs
function _equal_to_parse(arg1, arg2)
    if isexpr(arg1, :call)
        return Expr(:call, :(==), arg1.args[1], arg2), Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :(==), arg1, arg2), nothing
    end
end

# The variable is on the rhs
function _equal_to_parse(arg1::Number, arg2)
    return _equal_to_parse(arg2, arg1)
end

function _parse_parameters(_error::Function, head::Val{:call}, ::Union{Val{:<=}, Val{:≤}}, args)
    return _less_than_parse(args[2], args[3])
end

function _parse_parameters(_error::Function, head::Val{:call}, ::Union{Val{:>=}, Val{:≥}}, args)
    return _greater_than_parse(args[2], args[3])
end

function _parse_parameters(_error::Function, head::Val{:call}, ::Val{:(==)}, args)
    return _equal_to_parse(args[2], args[3])
end

function _parse_parameters(_error::Function, head::Val{:call}, first, args)
    if !(typeof(first) <: Union{Val{:in}, Val{:∈}})
        return args[1], Expr(:tuple, args[2:end]...)
    else
        _error("Bad in operator $first.")
    end
end

# TODO Make work with anonymous syntax, comparison heads, indexed bounds, and add checks
macro infinite_variable(model, expr, args...)
    _error(str...) = JuMP._macro_error(:infinite_parameter, (model, args...), str...)

    esc_model = esc(model)

    # extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)

    # Check for easy case if it is anonymous single variable
    # if length(extra) == 0
    #     code = quote
    #         @assert isa($esc_model, InfiniteModel)
    #         JuMP.@variable($esc_model, ($(args...)), variable_type = Infinite)
    #     end
    # else
        # x = popfirst!(extra)

    # There are several cases to consider:
    # x                                                     | type of x | x.head
    # ------------------------------------------------------+-----------+------------
    # var                                                   | Symbol    | NA
    # var[1:2]                                              | Expr      | :ref
    # var(x, y)                                             | Expr      | :call
    # var[1:2](x, y)                                        | Expr      | :call
    # var <= ub or var[1:2] <= ub                           | Expr      | :call
    # var(x, y) <= ub or var[1:2](x, y) <= ub               | Expr      | :call
    # lb <= var <= ub or lb <= var[1:2] <= ub               | Expr      | :comparison
    # lb <= var(x, y) <= ub or lb <= var[1:2](x, y) <= ub   | Expr      | :comparison
    x = expr
    if isexpr(x, :comparison) || isexpr(x, :call)
        inf_expr, params = InfOpt._parse_parameters(_error, Val(x.head), Val(x.args[1]), x.args)
    else
        x = inf_expr
        params = nothing
    end

    # new_args = (,)
    # for i = 2:length(args)
    #     new_args = (new_args..., args[i])
    # end
    # TODO make sure param_refs are not double specified
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(inf_expr)), ($(args...)), variable_type = Infinite, param_refs = $params)
    end
    # end
    return esc(code)
end

# TODO Streamline inputs and do checks, perhaps implement expression parsing
"""
    @point_variable(model, args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `PointVariableRef`.
"""
macro point_variable(model, inf_var, index, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Point, inf_var_ref = $inf_var, param_values = $index)
    end
    return esc(code)
end

"""
    @global_variable(model, args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `GlobalVariableRef`.
"""
macro global_variable(model, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Global)
    end
    return esc(code)
end
