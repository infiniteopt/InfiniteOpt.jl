using Base.Meta
using JuMP: _valid_model, _error_if_cannot_register, object_dictionary, variable_type

# Parse raw input to define the upper bound for an interval set
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:<=}, Val{:≤}},
    upper)
    JuMP._set_upper_bound_or_error(_error, infoexpr, upper)
end

# Parse raw input to define the lower bound for an interval set
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:>=}, Val{:≥}},
    lower)
    JuMP._set_lower_bound_or_error(_error, infoexpr, lower)
end

# Parse raw input to define the distribution for a dist set
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:in}, Val{:∈}},
    value)
    _dist_or_error(_error, infoexpr, value)
end

# Default function for unrecognized operators
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
    _parse_one_operator_parameter(_error, infoexpr,
                                  JuMP.reverse_sense(Val(sense)),
                                  JuMP._esc_non_constant(value))
    return var
end

# Parse raw input to define the upper and lower bounds for an interval set
function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                  ::Union{Val{:<=}, Val{:≤}}, lower,
                                  ::Union{Val{:<=}, Val{:≤}}, upper)
    JuMP._set_lower_bound_or_error(_error, infoexpr, lower)
    JuMP._set_upper_bound_or_error(_error, infoexpr, upper)
end

# Parse raw input to define the upper and lower bounds for an interval set
function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                  ::Union{Val{:>=}, Val{:≥}}, upper,
                                  ::Union{Val{:>=}, Val{:≥}}, lower)
    _parse_ternary_parameter(_error, infoexpr, Val(:≤), lower, Val(:≤), upper)
end

# Default for unrecognized operators
function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                  ::Val, lvalue, ::Val, rvalue)
    _error("Use the form lb <= ... <= ub.")
end

# Interpret a param_expr from macro input
function _parse_parameter(_error::Function, infoexpr::_ParameterInfoExpr, lvalue,
                          lsign::Symbol, var, rsign::Symbol, rvalue)
    # lvalue lsign var rsign rvalue
    _parse_ternary_parameter(_error, infoexpr, Val(lsign),
                             JuMP._esc_non_constant(lvalue), Val(rsign),
                             JuMP._esc_non_constant(rvalue))
    return var
end

# Extend to increase the param group id
function _assert_valid_model_call(m, macrocode)
    # assumes m is already escaped
    quote
        JuMP._valid_model($m, $(quot(m.args[1])))
        $(m).next_param_id += 1
        $macrocode
    end
end

"""
    @infinite_parameter(model, kw_args...)

Add an *anonymous* infinite parameter to the model `model` described by the
keyword arguments `kw_args` and returns the parameter reference.

    @infinite_parameter(model, expr, kw_args...)

Add a parameter to the model `model` described by the expression `expr`, the
positional arguments `args` and the keyword arguments `kw_args`. (note that in
the following the symbol `<=` can be used instead of `≤`, the symbol `>=`can
be used instead of `≥`, and the symbo `in` can be used instead of `∈`) The
expression `expr` can be of the form:
- `paramexpr` creating parameters described by `paramexpr`.
- `lb ≤ varexpr ≤ ub` creating parameters described by `paramexpr` characterized
   by a continuous interval set with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ dist` creating parameters described by `paramexpr` characterized
   by the `Distributions.jl` distribution object `dist`.

The expression `varexpr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The recognized keyword arguments in `kw_args` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `lower_bound`: Sets the value of the parameter lower bound for an interval set.
- `upper_bound`: Sets the value of the parameter upper bound for an interval set.
- `set`: The `InfiniteSet` characterizing the parameters see [`IntervalSet`](@ref)
   and [`DistributionSet`](@ref).
- `distribution`: Sets the `Distributions.jl` distribution object that characterizes
  the parameters.
- `supports`: Sets the support points for the parameters.
- `independent`: Specifies if the each parameter is independent from each other
  or not.
- `container`: Specify the container type.

 **Examples**
 ```julia
julia> @infinite_parameter(m, 0 <= x <= 1)
x

julia> supps = [[0, 1, 2], [-1, 1]];

julia> @infinite_parameter(m, y[i = 1:2] in Normal(), supports = supps[i])
2-element Array{ParameterRef,1}:
 y[1]
 y[2]

julia> z = @infinite_parameter(m, ["a", "b"], distribution = Uniform(), independent = true)
2-dimensional DenseAxisArray{ParameterRef,2,...} with index sets:
    Dimension 1, "a"
    Dimension 2, "b"
And data, a 1×1 Array{ParameterRef,2}:
 noname
 ```
"""
macro infinite_parameter(model, args...)
    _error(str...) = JuMP._macro_error(:infinite_parameter, (model, args...),
                                       str...)

    esc_model = esc(model)

    extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model
    if length(extra) == 0
        x = gensym()
        anon_singleton = true
    elseif length(extra) == 1
        x = popfirst!(extra)
        anon_singleton = false
    else
        x = popfirst!(extra)
        arg = popfirst!(extra)
        _error("Unrecognized argument $arg provided.")
    end

    info_kw_args = filter(_is_set_keyword, kw_args)
    extra_kw_args = filter(kw -> kw.args[1] != :base_name && !InfiniteOpt._is_set_keyword(kw), kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    infoexpr = InfiniteOpt._ParameterInfoExpr(; JuMP._keywordify.(info_kw_args)...)

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
        param = InfiniteOpt._parse_parameter(_error, infoexpr, x.args...)
    else
        param = x
    end

    anonvar = isexpr(param, :vect) || isexpr(param, :vcat) || anon_singleton
    anonvar && explicit_comparison && _error("Cannot use explicit bounds via " *
                                             ">=, <= with an anonymous parameter")
    parameter = gensym()
    name = JuMP._get_name(param)
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        _error("Expression $name should not be used as a parameter name. Use " *
               "the \"anonymous\" syntax $name = @infinite_parameter(model, " *
               "...) instead.")
    end

    set = InfiniteOpt._constructor_set(_error, infoexpr)
    if isa(param, Symbol)
        # Easy case - a single variable
        buildcall = :( build_parameter($_error, $set, 1, $(extra...)) )
        JuMP._add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $base_name,
                                         macro_call = true) )
        # The looped code is trivial here since there is a single variable
        creationcode = :($parameter = $parametercall)
        final_parameter = parameter
    else
        # isa(param, Expr) || _error("Expected $param to be a parameter name") --> not needed... I think
        # We now build the code to generate the variables (and possibly the
        # SparseAxisArray to contain them)
        refcall, idxparams, idxsets, condition = JuMP._build_ref_sets(param,
                                                                      parameter)
        clear_dependencies(i) = (JuMP.Containers.is_dependent(idxparams,
                                               idxsets[i], i) ? () : idxsets[i])

        # Code to be used to create each variable of the container.
        vartype = :( variable_type($esc_model, Parameter) )
        container_code, = JuMP.Containers.generate_container(vartype, idxparams,
                                                    idxsets, requestedcontainer)
        buildcall = :( build_parameter($_error, $set, length($container_code),
                                       $(extra...)) )
        JuMP._add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall,
                                         $(JuMP._name_call(base_name, idxparams)),
                                         macro_call = true) )
        code = :( $(refcall) = $parametercall )
        # Determine the return type of add_variable. This is needed to create
        # the container holding them.
        creationcode = JuMP._get_looped_code(parameter, code, condition,
                                             idxparams, idxsets, vartype,
                                             requestedcontainer)
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

## Define helper functions needed to parse infinite variable expressions
# Check rhs to to ensure is not a variable as best we can
function _check_rhs(arg1, arg2)
    if isexpr(arg2, :ref)
        if isexpr(arg2.args[2], :kw)
            temp = arg2
            arg2 = arg1
            arg1 = temp
        end
    # elseif isexpr(arg2, :call) --> infinite loop if both sides are :calls
    #     temp = arg2
    #     arg2 = arg1
    #     arg1 = temp
    end
    return arg1, arg2
end

# Assume variable on lhs, but check if swapped
function _less_than_parse(arg1, arg2)
    # check if swapped
    arg1_orig = arg1
    arg1, arg2 = _check_rhs(arg1, arg2)
    if arg1 != arg1_orig
        return _greater_than_parse(arg1, arg2)
    end
    # check if has parameter tuple
    if isexpr(arg1, :call)
        return Expr(:call, :<=, arg1.args[1], arg2),
                    Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :<=, arg1, arg2), nothing
    end
end

# Assume variable on lhs, but check if swapped
function _greater_than_parse(arg1, arg2)
    # check if swapped
    arg1_orig = arg1
    arg1, arg2 = _check_rhs(arg1, arg2)
    if arg1 != arg1_orig
        return _less_than_parse(arg1, arg2)
    end
    # check if has parameter tuple
    if isexpr(arg1, :call)
        return Expr(:call, :>=, arg1.args[1], arg2),
                    Expr(:tuple, arg1.args[2:end]...)
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
        return Expr(:call, :(==), arg1.args[1], arg2),
                    Expr(:tuple, arg1.args[2:end]...)
    else
        return Expr(:call, :(==), arg1, arg2), nothing
    end
end

# The variable is on the rhs
function _equal_to_parse(arg1::Number, arg2)
    return _equal_to_parse(arg2, arg1)
end

# Return the tuple of form (var_expr, param_tuple_expr)
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

# Return the tuple of form (var_expr, param_tuple_expr)
function _parse_parameters(_error::Function, head::Val{:comparison}, args)
    if isexpr(args[3], :call)
        return Expr(:comparison, args[1:2]..., args[3].args[1], args[4:5]...),
                    Expr(:tuple, args[3].args[2:end]...)
    else
        return Expr(:comparison, args...), nothing
    end
end

"""
    @infinte_variable(model, kw_args...)

Add an *anonymous* infinite variable to the model `model` described by the
keyword arguments `kw_args` and returns the variable reference. Note that the
`parameter_refs` keyword is required in this case.

    @infinite_variable(model, varexpr, args..., kw_args...)

Add an infinite variable to `model` described by the expression `var_expr`, the
positional arguments `args` and the keyword arguments `kw_args`. The expression
`varexpr` can either be (note that in the following the symbol `<=` can be used
instead of `≤` and the symbol `>=`can be used instead of `≥`) of the form:

- `varexpr` creating variables described by `varexpr`
- `varexpr ≤ ub` (resp. `varexpr ≥ lb`) creating variables described by
  `varexpr` with upper bounds given by `ub` (resp. lower bounds given by `lb`)
- `varexpr == value` creating variables described by `varexpr` with fixed values
   given by `value`
- `lb ≤ varexpr ≤ ub` or `ub ≥ varexpr ≥ lb` creating variables described by
  `varexpr` with lower bounds given by `lb` and upper bounds given by `ub`

The expression `varexpr` can be of the form:

- `varname` creating a scalar real variable of name `varname`
- `varname(params...)` creating a scalar real variable of name `varname` with
  infinite parameters `params...` see `parameter_refs` for format.
- `varname[...]` or `[...]` creating a container of variables.
- `varname[...](params...)` or `[...]` creating a container of variables with
  infinite parameters `params...` in the first case.

The recognized positional arguments in `args` are the following:

- `Bin`: Sets the variable to be binary, i.e. either 0 or 1.
- `Int`: Sets the variable to be integer, i.e. one of ..., -2, -1, 0, 1, 2, ...

The recognized keyword arguments in `kw_args` are the following:

- `parameter_refs`: This is mandatory if not specified in `varexpr`. Can be a
  single parameter reference, a single parameter array with parameters defined
  in the same call of [`@infinite_parameter`](@ref) (i.e., have same group ID),
  or a tuple where each element is either of the first two options listed.
- `base_name`: Sets the name prefix used to generate variable names. It
  corresponds to the variable name for scalar variable, otherwise, the
  variable names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `lower_bound`: Sets the value of the variable lower bound.
- `upper_bound`: Sets the value of the variable upper bound.
- `start`: Sets the variable starting value used as initial guess in optimization.
- `binary`: Sets whether the variable is binary or not.
- `integer`: Sets whether the variable is integer or not.
- `container`: Specify the container type.
"""
macro infinite_variable(model, args...)
    _error(str...) = JuMP._macro_error(:infinite_parameter, (model, args...),
                                       str...)

    extra, kw_args, requestedcontainer = JuMP._extract_kw_args(args)
    param_kw_args = filter(kw -> kw.args[1] == :parameter_refs, kw_args)

    # Check for easy case if it is anonymous single variable
    if length(extra) == 0
        code = quote
            @assert isa($model, InfiniteModel) "Model must be an `InfiniteModel`."
            JuMP.@variable($model, ($(args...)), variable_type = Infinite,
                           error = $_error)
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
            inf_expr, params = InfiniteOpt._parse_parameters(_error, Val(x.head),
                                                             x.args)
        else
            inf_expr = x
            params = nothing
        end

        if length(param_kw_args) != 0 && params != nothing
            _error("Cannot specify double specify the infinite parameter " *
                   "references.")
        end

        if length(args) == 1
            if params == nothing
                code = quote
                    @assert isa($model, InfiniteModel) "Model must be an " *
                                                       "`InfiniteModel`."
                    JuMP.@variable($model, ($(inf_expr)),
                                   variable_type = Infinite, error = $_error)
                end
            else
                code = quote
                    @assert isa($model, InfiniteModel) "Model must be an " *
                                                       "`InfiniteModel`."
                    JuMP.@variable($model, ($(inf_expr)),
                                   variable_type = Infinite,
                                   parameter_refs = $params, error = $_error)
                end
            end
        else
            rest_args = [args[i] for i = 2:length(args)]
            if params == nothing
                code = quote
                    @assert isa($model, InfiniteModel) "Model must be an " *
                                                       "`InfiniteModel`."
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)),
                                   variable_type = Infinite, error = $_error)
                end
            else
                code = quote
                    @assert isa($model, InfiniteModel) "Model must be an " *
                                                       "`InfiniteModel`."
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)),
                                   variable_type = Infinite,
                                   parameter_refs = $params, error = $_error)
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
