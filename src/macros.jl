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

function _assert_valid_model_call(m, macrocode)
    # assumes m is already escaped
    quote
        JuMP._valid_model($m, $(quot(m.args[1])))
        $(m).next_param_id += 1
        $macrocode
    end
end

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

    info_kw_args = filter(_is_set_keyword, kw_args)
    extra_kw_args = filter(kw -> kw.args[1] != :base_name && !InfOpt._is_set_keyword(kw), kw_args)
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
        buildcall = :( build_parameter($_error, $set, 1, $(extra...)) )
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
        vartype = :( variable_type($esc_model, Parameter) )
        container_code, = JuMP.Containers.generate_container(vartype, idxparams, idxsets, requestedcontainer)
        buildcall = :( build_parameter($_error, $set, length($container_code), $(extra...)) )
        JuMP._add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $(JuMP._name_call(base_name, idxparams))) )
        code = :( $(refcall) = $parametercall )
        # Determine the return type of add_variable. This is needed to create the container holding them.
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
                                                   final_variable = final_parameter,
                                                   model_for_registering = esc_model)
    end
    return _assert_valid_model_call(esc_model, macro_code)
end

# Check rhs to to ensure is not a variable
function _check_rhs(arg1, arg2)
    if isexpr(arg2, :ref)
        if isexpr(arg2, :kw)
            temp = arg2
            arg2 = arg1
            arg1 = temp
        end
    elseif isexpr(arg2, :call)
        temp = arg2
        arg2 = arg1
        arg1 = temp
    end
    return arg1, arg2
end

# Assume variable on lhs
function _less_than_parse(arg1, arg2)
    arg1, arg2 = _check_rhs(arg1, arg2)
    if isexpr(arg1, :call)
        return Expr(:call, :<=, arg1.args[1], arg2), Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :<=, arg1, arg2), nothing
    end
end

# Assume variable on lhs
function _greater_than_parse(arg1, arg2)
    arg1, arg2 = _check_rhs(arg1, arg2)
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
    arg1, arg2 = _check_rhs(arg1, arg2)
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

function _parse_parameters(_error::Function, head::Val{:call}, args)
    if args[1] in [:<=, :≤]
        return _less_than_parse(args[2], args[3])
    elseif args[1] in [:>=, :≥]
        return _greater_than_parse(args[2], args[3])
    elseif args[1] == :(==)
        return _equal_to_parse(args[2], args[3])
    elseif !(args[1] in [:in, :∈])
        return args[1], Expr(:tuple, args[2:end]...)
    else
        first = args[1]
        _error("Invalid operator $first.")
    end
end

function _parse_parameters(_error::Function, head::Val{:comparison}, args)
    if isexpr(args[3], :call)
        return Expr(:comparison, args[1:2]..., args[3].args[1], args[4:5]...), Expr(:tuple, args[3].args[2:end]...)
    else
        return Expr(:comparison, args...), nothing
    end
end

"""
    @infinite_variable(model, args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `InfiniteVariableRef` Support syntax
@infinite_variable(m, x(params...), args...) or @infinite_variable(m, parameter_refs= (params...), kwargs...).
"""
macro infinite_variable(model, args...)
    _error(str...) = JuMP._macro_error(:infinite_parameter, (model, args...), str...)

    extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)
    param_kw_args = filter(kw -> kw.args[1] == :parameter_refs, kw_args)

    # Check for easy case if it is anonymous single variable
    if length(extra) == 0
        code = quote
            @assert isa($model, InfiniteModel)
            JuMP.@variable($model, ($(args...)), variable_type = Infinite, error = $_error)
        end
    else
        x = popfirst!(extra)

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
        if isexpr(x, :comparison) || isexpr(x, :call)
            inf_expr, params = InfOpt._parse_parameters(_error, Val(x.head), x.args)
        else
            inf_expr = x
            params = nothing
        end

        if length(param_kw_args) != 0 && params != nothing
            _error("Cannot specify double specify the infinite parameter references.")
        end

        if length(args) == 1
            if params == nothing
                code = quote
                    @assert isa($model, InfiniteModel)
                    JuMP.@variable($model, ($(inf_expr)), variable_type = Infinite, error = $_error)
                end
            else
                code = quote
                    @assert isa($model, InfiniteModel)
                    JuMP.@variable($model, ($(inf_expr)), variable_type = Infinite, parameter_refs = $params, error = $_error)
                end
            end
        else
            rest_args = [args[i] for i = 2:length(args)]
            if params == nothing
                code = quote
                    @assert isa($model, InfiniteModel)
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)), variable_type = Infinite, error = $_error)
                end
            else
                code = quote
                    @assert isa($model, InfiniteModel)
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)), variable_type = Infinite, parameter_refs = $params, error = $_error)
                end
            end
        end
    end
    return esc(code)
end

"""
#     @point_variable(model, args...)
# A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
# it defines variables of type `PointVariableRef` Support syntax
@point_variable(m, inf_var(param_vals...), args...) or @point_variable(m, infinite_variable_ref = inf_var, parameter_values = param_vals, kw_args...).
# """
macro point_variable(model, args...)
    _error(str...) = JuMP._macro_error(:point_parameter, (model, args...), str...)

    extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)
    param_kw_args = filter(kw -> kw.args[1] == :infinite_variable_ref || kw.args[1] == :parameter_values, kw_args)

    # Check for easy case if it is anonymous single variable
    if length(extra) == 0
        code = quote
            @assert isa($model, InfiniteModel)
            JuMP.@variable($model, ($(args...)), variable_type = Point, error = $_error)
        end
    else
        x = popfirst!(extra)

        if isexpr(x, :call)
            if length(param_kw_args) != 0
                _error("Cannot double specify the infinite variable reference and/or its paramter values.")
            elseif x.args[1] in [:in, :<=, :>=, :(==), :≥, :≤, ∈]
                _error("Invalid input syntax.")
            end
            # TODO handle vectorized point variables with vector infinite vars
            rest_args = [args[i] for i = 2:length(args)]
            inf_var = x.args[1]
            param_vals = Expr(:tuple, x.args[2:end]...)
            code = quote
                @assert isa($model, InfiniteModel)
                JuMP.@variable($model, ($(rest_args...)), variable_type = Point, infinite_variable_ref = $inf_var, parameter_values = $param_vals, error = $_error)
            end
        elseif isexpr(x, :vect) && length(extra) == 0
            code = quote
                @assert isa($model, InfiniteModel)
                JuMP.@variable($model, ($(args...)), variable_type = Point, error = $_error)
            end
        else
            _error("Invalid input syntax.")
        end
    end
    return esc(code)
end

"""
    @global_variable(model, args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `GlobalVariableRef`.
"""
macro global_variable(model, args...)
    _error(str...) = JuMP._macro_error(:global_parameter, (model, args...), str...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Global, error = $_error)
    end
    return esc(code)
end

# TODO make custom constraint macro with better subset syntax (i.e. using form
# @constraint(m, name[i = 1:2](t in [0,1], x[1] in [-1, 0], expr)))
# TODO Add bridge constraints for chance constraints and derivatives.
