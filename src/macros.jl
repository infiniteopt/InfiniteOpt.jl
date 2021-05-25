const Containers = JuMP.Containers

# Parse raw input to define the upper bound for an interval domain
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:<=}, Val{:≤}},
    upper)
    _set_upper_bound_or_error(_error, infoexpr, upper)
end

# Parse raw input to define the lower bound for an interval domain
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:>=}, Val{:≥}},
    lower)
    _set_lower_bound_or_error(_error, infoexpr, lower)
end

# Parse raw input to define the distribution for a dist domain
function _parse_one_operator_parameter(
    _error::Function, infoexpr::_ParameterInfoExpr, ::Union{Val{:in}, Val{:∈}},
    value)
    # check if interval domain
    if isexpr(value.args[1], :vect)
        _set_lower_bound_or_error(_error, infoexpr,
                                  _esc_non_constant(value.args[1].args[1]))
        _set_upper_bound_or_error(_error, infoexpr,
                                  _esc_non_constant(value.args[1].args[2]))
    else
        _domain_or_error(_error, infoexpr, value)
    end
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
                                _esc_non_constant(value))
    return var
end

# If the lhs is a number and not the rhs, we can deduce that the rhs is
# the variable.
function _parse_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                          sense::Symbol, value::Number, var)
    _parse_one_operator_parameter(_error, infoexpr,
                                  JuMP.reverse_sense(Val(sense)),
                                  _esc_non_constant(value))
    return var
end

# Parse raw input to define the upper and lower bounds for an interval domain
function _parse_ternary_parameter(_error::Function, infoexpr::_ParameterInfoExpr,
                                  ::Union{Val{:<=}, Val{:≤}}, lower,
                                  ::Union{Val{:<=}, Val{:≤}}, upper)
    _set_lower_bound_or_error(_error, infoexpr, lower)
    _set_upper_bound_or_error(_error, infoexpr, upper)
end

# Parse raw input to define the upper and lower bounds for an interval domain
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
                             _esc_non_constant(lvalue), Val(rsign),
                             _esc_non_constant(rvalue))
    return var
end

"""
    @independent_parameter(model::InfiniteModel, kw_args...)::GeneralVariableRef

Add an *anonymous* infinite parameter to the model `model` described by the
keyword arguments `kw_args` and returns the parameter reference.

```julia
@independent_parameter(model::InfiniteModel, expr, kw_args...)::GeneralVariableRef
```

Add a parameter to the model `model` described by the expression `expr`
and the keyword arguments `kw_args`. (note that in
the following the symbol `<=` can be used instead of `≤`, the symbol `>=`can
be used instead of `≥`, and the symbol `in` can be used instead of `∈`) The
expression `expr` can be of the form:
- `paramexpr` creating parameters described by `paramexpr`.
- `lb ≤ paramexpr ≤ ub` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ [lb, ub]` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ dist` creating parameters described by `paramexpr` characterized
   by the `Distributions.jl` distribution object `dist`.

The expression `paramexpr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The recognized keyword arguments in `kw_args` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `lower_bound`: Sets the value of the parameter lower bound for an interval domain.
- `upper_bound`: Sets the value of the parameter upper bound for an interval domain.
- `domain`: The `InfiniteDomain` characterizing the parameters see [`IntervalDomain`](@ref)
   and [`UniDistributionDomain`](@ref).
- `distribution`: Sets the `Distributions.jl` distribution object that characterizes
  the parameters.
- `supports`: Sets the support points for the parameters.
- `num_supports`: Specifies the number of supports to be automatically generated.
                  Note that `supports` takes precedence. Defaults to `DefaultNumSupports`.
- `derivative_method`: Specifies the numerical method used to evaluate derivatives that 
                       are taken with respect to the parameter.
- `sig_digits`: Specifies the number of significant digits that should be used
              in automatic support generation. Defaults to `DefaultSigDigits`.
- `container`: Specify the container type. Defaults to `automatic`

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP, Distributions; m = InfiniteModel())
julia> @independent_parameter(m, x in [0, 1])
x

julia> @independent_parameter(m, y[i = 1:2] in Normal(), num_supports = 10)
2-element Array{GeneralVariableRef,1}:
 y[1]
 y[2]

julia> z = @independent_parameter(m, [["a", "b"]], lower_bound = 0,
                                  upper_bound = 1, supports = [0, 0.5, 1])
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
   Dimension 1, ["a", "b"]
And data, a 2-element Array{GeneralVariableRef,1}:
noname[a]
noname[b]
```
"""
macro independent_parameter(model, args...)
    esc_model = esc(model)

    extra, kw_args, requestedcontainer = _extract_kw_args(args)

    # check to see if error_args are provided define error function
    error_kwargs = filter(kw -> kw.args[1] == :error_func, kw_args)
    if isempty(error_kwargs)
        _error(str...) = _macro_error(:independent_parameter, (model, args...), 
                                      __source__, str...)
    else
        _error = error_kwargs[1].args[2]
    end

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

    info_kw_args = filter(InfiniteOpt._is_domain_keyword, kw_args)
    extra_kw_args = filter(kw -> kw.args[1] != :base_name && 
                           !InfiniteOpt._is_domain_keyword(kw) && 
                           kw.args[1] != :error_func, kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    infoexpr = InfiniteOpt._ParameterInfoExpr(; InfiniteOpt._keywordify.(info_kw_args)...)

    # There are five cases to consider:
    # x                                         | type of x | x.head
    # ------------------------------------------+-----------+------------
    # param                                       | Symbol    | NA
    # param in [lb, ub]                           | Expr      | :call
    # lb <= param <= ub                           | Expr      | :comparison
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
    name = _get_name(param)
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        _error("Expression $name should not be used as a parameter name. Use " *
               "the \"anonymous\" syntax $name = @independent_parameter(model, " *
               "...) instead.")
    end

    domain = InfiniteOpt._constructor_domain(_error, infoexpr)
    if isa(param, Symbol)
        # Easy case - a single variable
        buildcall = :( build_parameter($_error, $domain) )
        _add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $base_name) )
        creationcode = :($parameter = $parametercall)
    else
        # isa(param, Expr) || _error("Expected $param to be a parameter name") --> not needed... I think
        # We now build the code to generate the variables (and possibly the
        # SparseAxisArray to contain them)
        idxvars, indices = JuMPC._build_ref_sets(_error, param)
        buildcall = :( build_parameter($_error, $domain) )
        _add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall,
                                         $(_name_call(base_name, idxvars))) )
        creationcode = JuMPC.container_code(idxvars, indices, parametercall,
                                             requestedcontainer)
    end
    if anonvar
        # Anonymous variable, no need to register it in the model-level
        # dictionary nor to assign it to a variable in the user scope.
        # We simply return the variable
        macro_code = creationcode
    else
        # We register the variable reference to its name and
        # we assign it to a variable in the local scope of this name
        macro_code = _macro_assign_and_return(creationcode, parameter, name,
                                              model_for_registering = esc_model)
    end
    return _finalize_macro(esc_model, macro_code, __source__)
end

"""
    @finite_parameter(model::InfiniteModel, value)::GeneralVariableRef

Define and add an anonymous finite parameter to `model` and return its
parameter reference. Its value is equal to `value`.

```julia
    @finite_parameter(model::InfiniteModel, param_expr,
                      value_expr)::GeneralVariableRef
```
Define and add a finite parameter(s) to `model` and return appropriate parameter
reference(s). The parameter(s) has/have value(s) indicated by the `value_expr`.
The expression `param_expr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The expression `value_expr` simply expresses the value of the parameter(s). This
is typically a number but could be an array indexed using an index defined in
`param_expr`.

The recognized keyword arguments in `kw_args` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `container`: Specify the container type.

**Examples**
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> par = @finite_parameter(model, 2)
noname

julia> vals = [3, 2];

julia> pars = @finite_parameter(model, [i = 1:2], vals[i], base_name = "par")
2-element Array{ParameterRef,1}:
 par[1]
 par[2]

julia> @finite_parameter(model, par2, 42)
par2
```
"""
macro finite_parameter(model, args...)
    esc_model = esc(model)

    extra, kw_args, requestedcontainer = _extract_kw_args(args)

    macro_name = :finite_parameter
    _error(str...) = _macro_error(macro_name, (model, args...), __source__, str...)

    # if there is only a single non-keyword argument, this is an anonymous
    # variable spec and the one non-kwarg is the model

    if length(extra) == 1
        param = gensym()
        value = popfirst!(extra)
        anon_singleton = true
    elseif length(extra) == 2
        param = popfirst!(extra)
        value = popfirst!(extra)
        anon_singleton = false
    else
        _error("Incorrect number of arguments. Must be of form " *
               "@finite_parameter(model, name_expr, value).")
    end
    value = _esc_non_constant(value)

    extra_kw_args = filter(kw -> kw.args[1] != :base_name, kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)

    anonvar = isexpr(param, :vect) || isexpr(param, :vcat) || anon_singleton

    parameter = gensym()
    name = _get_name(param)
    if isempty(base_name_kw_args)
        base_name = anonvar ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end

    if !isa(name, Symbol) && !anonvar
        _error("Expression $name should not be used as a parameter name. Use " *
               "the \"anonymous\" syntax $name = @finite_parameter(model, " *
               "...) instead.")
    end

    if isa(param, Symbol)
        buildcall = :( build_parameter($_error, $value) )
        _add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall, $base_name) )
        creationcode = :($parameter = $parametercall)
    else
        idxvars, indices = JuMPC._build_ref_sets(_error, param)
        buildcall = :( build_parameter($_error, $value) )
        _add_kw_args(buildcall, extra_kw_args)
        parametercall = :( add_parameter($esc_model, $buildcall,
                                         $(_name_call(base_name, idxvars))) )
        creationcode = JuMPC.container_code(idxvars, indices, parametercall,
                                             requestedcontainer)
    end
    if anonvar
        macro_code = creationcode
    else
        macro_code = _macro_assign_and_return(creationcode, parameter, name,
                                              model_for_registering = esc_model)
    end
    return _finalize_macro(esc_model, macro_code, __source__)
end

"""
    @dependent_parameters(model::InfiniteModel, kw_args...)::GeneralVariableRef

Add *anonymous* dependent infinite parameters to the model `model` described by
the keyword arguments `kw_args` and returns the container of parameter references.

```julia
@dependent_parameters(model::InfiniteModel, expr, kw_args...)::GeneralVariableRef
```

Add a container of dependent infinite parameters to the model `model` described
by the expression `expr`, and the keyword arguments `kw_args`. (note that in
the following the symbol `<=` can be used instead of `≤`, the symbol `>=`can
be used instead of `≥`, and the symbol `in` can be used instead of `∈`) The
expression `expr` can be of the form:
- `paramexpr` creating parameters described by `paramexpr`.
- `lb ≤ paramexpr ≤ ub` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ [lb, ub]` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ dist` creating parameters described by `paramexpr` characterized
   by the `Distributions.jl` distribution object `dist`.
- `paramexpr ∈ domain` creating parameters described by `paramexpr` characterized
   by the `AbstractInfiniteDomain` object `domain`.

The expression `paramexpr` must be of the form:
- `paramname[...]` or `[...]` creating a container of parameters

The recognized keyword arguments in `kw_args` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `lower_bound`: Sets the value of the parameter lower bound for an interval domain.
- `upper_bound`: Sets the value of the parameter upper bound for an interval domain.
- `domain`: The `InfiniteDomain` characterizing the parameters that are subtypes of
         [`AbstractInfiniteDomain`](@ref).
- `distribution`: Sets the `Distributions.jl` distribution object that characterizes
  the parameters.
- `supports`: Sets the support points for the parameters.
- `num_supports`: Specifies the number of supports to be automatically generated.
                  Note that `supports` takes precedence. Defaults to `DefaultNumSupports`.
- `derivative_method`: Specifies the numerical derivative method used to evaluate 
                       derivatives that depend on a particular dependent parameter.
- `sig_digits`: Specifies the number of significant digits that should be used
              in automatic support generation. Defaults to `DefaultSigDigits`.
- `container`: Specify the container type. Defaults to `automatic`.

**Examples**
```julia-repl
julia> @dependent_parameters(m, x[1:2] in [0, 1])
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]

julia> @dependent_parameters(m, y[i = 1:2] in MvNormal(ones(2)), num_supports = 10)
2-element Array{GeneralVariableRef,1}:
 y[1]
 y[2]

julia> z = @dependent_parameters(m, [i = ["a", "b"], j = 1:2],
                                 distribution = MatrixBeta(2, 2, 2))
2-dimensional DenseAxisArray{GeneralVariableRef,2,...} with index sets:
    Dimension 1, ["a", "b"]
    Dimension 2, Base.OneTo(2)
And data, a 2×2 Array{GeneralVariableRef,2}:
 noname[a,1]  noname[a,2]
 noname[b,1]  noname[b,2]
```
"""
macro dependent_parameters(model, args...)
    # extract the raw arguments
    esc_model = esc(model)
    extra, kw_args, requestedcontainer = _extract_kw_args(args)

    # check to see if error_args are provided define error function
    error_kwargs = filter(kw -> kw.args[1] == :error_func, kw_args)
    if isempty(error_kwargs)
        _error(str...) = _macro_error(:dependent_parameters, (model, args...), 
                                      __source__, str...)
    else
        _error = error_kwargs[1].args[2]
    end

    # there must be one extra positional argument to specify the dimensions
    if length(extra) == 1
        x = popfirst!(extra) 
    elseif length(extra) == 0
        _error("Must specify more than one dependent parameter.")
    else
        x = popfirst!(extra)
        arg = popfirst!(extra)
        _error("Unrecognized argument $arg provided.")
    end

    # parse the keyword arguments
    info_kw_args = filter(_is_domain_keyword, kw_args)
    supports_kw_arg = filter(kw -> kw.args[1] == :supports, kw_args)
    method_kw_arg = filter(kw -> kw.args[1] == :derivative_method, kw_args)
    extra_kw_args = filter(kw -> kw.args[1] != :base_name &&
                           !InfiniteOpt._is_domain_keyword(kw) &&
                           kw.args[1] != :error_func &&
                           kw.args[1] != :supports &&
                           kw.args[1] != :derivative_method, kw_args)
    base_name_kw_args = filter(kw -> kw.args[1] == :base_name, kw_args)
    infoexpr = InfiniteOpt._ParameterInfoExpr(; JuMP._keywordify.(info_kw_args)...)

    # There are six cases to consider:
    # x                                         | type of x | x.head
    # ------------------------------------------+-----------+------------
    # param                                       | Symbol    | NA
    # [1:2]                                       | Expr      | :vect
    # param[1:2]                                  | Expr      | :ref
    # param[1:2] <= ub                            | Expr      | :call
    # param[1:2] in [lb, ub] (other domains)         | Expr      | :call
    # lb <= param[1:2] <= ub                      | Expr      | :comparison
    # In the three last cases, we call parse_variable
    symbolic_domain = isexpr(x, :comparison) || isexpr(x, :call)
    if symbolic_domain
        param = InfiniteOpt._parse_parameter(_error, infoexpr, x.args...)
    elseif isa(x, Symbol)
        _error("Must specify more than one dependent parameter.")
    else
        param = x
    end

    # determine it is an anonymous call
    anonparam = isexpr(param, :vcat) || isexpr(param, :vect)
    anonparam && symbolic_domain && _error("Cannot use symbolic infinite defintion " *
                                        "with an anonymous parameter")
    # process the parameter name
    parameter = gensym()
    name = _get_name(param)
    if isempty(base_name_kw_args)
        base_name = anonparam ? "" : string(name)
    else
        base_name = esc(base_name_kw_args[1].args[2])
    end
    if !isa(name, Symbol) && !anonparam
        _error("Expression $name should not be used as a parameter name. Use " *
               "the \"anonymous\" syntax $name = @dependent_parameters(model, " *
               "...) instead.")
    end

    # process the supports
    if isempty(supports_kw_arg)
        supports = :(Float64[])
    else
        supports = esc(supports_kw_arg[1].args[2])
    end

    # process the derivative methods
    if isempty(method_kw_arg)
        method = :(InfiniteOpt.DefaultDerivativeMethod)
    else
        method = esc(method_kw_arg[1].args[2])
    end

    # parse the infinite domain(s)
    domain = InfiniteOpt._construct_array_domain(_error, infoexpr)

    # make code to build the DependentParameters object and the references
    idxvars, indices = JuMPC._build_ref_sets(_error, param)
    namecode = _name_call(base_name, idxvars)
    parambuildcall = :( InfiniteOpt._DependentParameter($domain, $supports, $namecode, $method) )
    param_container_call = JuMPC.container_code(idxvars, indices, parambuildcall,
                                                requestedcontainer)
    buildcall = :( InfiniteOpt._build_parameters($_error, $param_container_call) )
    _add_kw_args(buildcall, extra_kw_args)
    creationcode = :( add_parameters($esc_model, $(buildcall)...) )

    # make the final return code based on if it is anonymous or not
    if anonparam
        macro_code = creationcode
    else
        macro_code = _macro_assign_and_return(creationcode, parameter, name,
                                              model_for_registering = esc_model)
    end
    return _finalize_macro(esc_model, macro_code, __source__)
end

"""
    @infinite_parameter(model::InfiniteModel, kw_args...)::GeneralVariableRef

Add *anonymous* infinite parameter to the model `model` described by
the keyword arguments `kw_args` and returns the parameter reference.

```julia
@infinite_parameter(model::InfiniteModel, expr, kw_args...)::GeneralVariableRef
```

Add an infinite parameter to the model `model` described by the expression `expr`,
and the keyword arguments `kw_args`. This is just a wrapper macro that will make
the appropriate call to either [`@independent_parameter`](@ref) or
[`@dependent_parameters`](@ref). (Note that in
the following the symbol `<=` can be used instead of `≤`, the symbol `>=`can
be used instead of `≥`, and the symbol `in` can be used instead of `∈`) The
expression `expr` can be of the form:
- `paramexpr` creating parameters described by `paramexpr`.
- `lb ≤ paramexpr ≤ ub` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ [lb, ub]` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ∈ dist` creating parameters described by `paramexpr` characterized
   by the `Distributions.jl` distribution object `dist`.
- `paramexpr ∈ domain` creating parameters described by `paramexpr` characterized
  by the `AbstractInfiniteDomain` object `domain`.

The expression `paramexpr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The recognized keyword arguments in `kw_args` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `lower_bound`: Sets the value of the parameter lower bound for an interval domain.
- `upper_bound`: Sets the value of the parameter upper bound for an interval domain.
- `domain`: The `InfiniteDomain` characterizing the parameters see subtypes of
         [`AbstractInfiniteDomain`](@ref).
- `distribution`: Sets the `Distributions.jl` distribution object that characterizes
  the parameters.
- `supports`: Sets the support points for the parameters.
- `num_supports`: Specifies the number of supports to be automatically generated.
                  Note that `supports` takes precedence. Defaults to `DefaultNumSupports`.
- `derivative_method`: Specify the numerical method to evaluate derivatives that are
                       taken with respect to the parameter.
- `sig_digits`: Specifies the number of significant digits that should be used
              in automatic support generation. Defaults to `DefaultSigDigits`.
- `independent`: Specifies if the each parameter is independent from each other
  or not. Defaults to false.
- `container`: Specify the container type. Defaults to `automatic`

**Examples**
```julia-repl
julia> @infinite_parameter(m, x in [0, 1])
x

julia> @infinite_parameter(m, y[i = 1:2] in MvNormal(ones(2)), num_supports = 10)
2-element Array{GeneralVariableRef,1}:
 y[1]
 y[2]

julia> z = @infinite_parameter(m, [["a", "b"]], distribution = Normal(),
                               independent = true)
1-dimensional DenseAxisArray{GeneralVariableRef,1,...} with index sets:
    Dimension 1, ["a", "b"]
And data, a 2-element Array{GeneralVariableRef,1}:
 noname[a]
 noname[b]
```
"""
macro infinite_parameter(model, args...)
    # define error message function
    _error(str...) = _macro_error(:infinite_parameter, (model, args...),
                                  __source__, str...)
    # parse the arguments
    extra, kw_args, requestedcontainer = _extract_kw_args(args)
    indep_kwarg = filter(kw -> kw.args[1] == :independent, kw_args)
    new_kwargs = filter(kw -> kw.args[1] != :independent, kw_args)
    if isempty(indep_kwarg)
        independent = false
    else
        independent = indep_kwarg[1].args[2]
    end

    # determine if there is only one parameter
    if length(extra) == 0
        code = :( @independent_parameter($model, $(extra...), $(new_kwargs...),
                                         container = $requestedcontainer,
                                         error_func = $(_error)) )
    else
        # get the param expression and check it is an array --> TODO make more efficient
        x = first(extra)
        infoexpr = InfiniteOpt._ParameterInfoExpr()
        symbolic_domain = isexpr(x, :comparison) || isexpr(x, :call)
        if symbolic_domain
            param = InfiniteOpt._parse_parameter(_error, infoexpr, x.args...)
        else
            param = x
        end
        if param isa Symbol
            code = :( @independent_parameter($model, $(extra...), $(new_kwargs...),
                                             container = $requestedcontainer,
                                             error_func = $(_error)) )
        else
            code = quote
                if $independent
                    @independent_parameter($model, $(extra...), $(new_kwargs...),
                                           container = $requestedcontainer,
                                           error_func = $(_error))
                else
                    @dependent_parameters($model, $(extra...), $(new_kwargs...),
                                          container = $requestedcontainer,
                                          error_func = $(_error))
                end
            end
        end
    end
    return esc(code)
end

# Helper method for checking that function call is not from an operator
function _ensure_not_operator(_error::Function, func)
    if func in [:in, :<=, :>=, :(==), :≥, :≤, ∈]
        _error("Invalid input syntax.")
    end
    return
end

# Helper method to process parameter function expressions 
function _process_func_expr(_error::Function, raw_expr, pref_expr)
    if isexpr(raw_expr, :call)
        # we have symbolic arguments given 
        _ensure_not_operator(_error, raw_expr.args[1])
        func_expr = _esc_non_constant(raw_expr.args[1])
        # restructure to expected in case keywords were given using ;
        if length(raw_expr.args) > 1 && isexpr(raw_expr.args[2], :parameters)
            append!(raw_expr.args, raw_expr.args[2].args)
            deleteat!(raw_expr.args, 2)
        end
        # check that the prefs are given first and match the pref_expr
        num_pref_args = length(pref_expr.args)
        if length(raw_expr.args) <= num_pref_args || pref_expr.args != raw_expr.args[2:num_pref_args+1]
            _error("Inconsistent infinite parameter reference dependencies.")
        end 
        # extract the extra arguments given symbolically
        extra_func_args = raw_expr.args[num_pref_args+2:end]
        raw_fargs = filter(e -> !isexpr(e, :kw), extra_func_args)
        if isempty(raw_fargs) 
            fargs = nothing
        else
            fargs = esc(Expr(:tuple, raw_fargs...))
        end
        # extract the keyword arguments given symbolically
        raw_fkwargs = filter(e -> isexpr(e, :kw), extra_func_args)
        if isempty(raw_fkwargs)
            fkwargs = nothing
        else
            fkwargs = esc(Expr(:tuple, map(e -> Expr(:(=), e.args...), raw_fkwargs)...))
        end
    else 
        # easy case where only the function name was given
        func_expr = _esc_non_constant(raw_expr)
        fargs = nothing
        fkwargs = nothing
    end
    return func_expr, fargs, fkwargs
end

"""
    @parameter_function(model::InfiniteModel, kw_args...)::GeneralVariableRef

Add an *anonymous* parameter function to the model `model` described by the
keyword arguments `kw_args` and returns the object reference. Note that the
`func` and `parameter_refs` keywords are required in this
case.

```julia
@parameter_function(model::InfiniteModel, varexpr, funcexpr, kw_args...)::GeneralVariableRef

@parameter_function(model::InfiniteModel, varexpr == funcexpr, kw_args...)::GeneralVariableRef
```

Add a parameter function to `model` described by the expression `varexpr`, the
function expression `funcexpr`, and the keyword arguments `kw_args`. The expression 
`varexpr` is used to define the parameter function references and determine the 
infinite parameters that are intended for this function. Here the accepted forms 
are:
- `varname(params...)` creating a scalar parameter function with alias name `varname` 
  that depends on the parameter references `params...`
- `varname[...](params...)` or `[...](params...)` creating a container of parameter 
  functions with parameter references given in `paramexpr`.
The expression `params` can match that of the `parameter_refs` keyword argument 
(see the keyword argument description for info on accepted formats). 

The expression `funcexpr` determines the concrete Julia function that defines the 
behavior of parameter function and allows us to specify the `func`, `func_args`, 
and/or `func_kwargs` keyword arguments. The accepted forms are:
- `funcname`: Give the name of the concrete Julia function or give an anonymous function
- `funcname(params..., func_args..., func_kwargs...)`: Provide the concrete Julia 
  function name and provide the addtional positional/keyword arguments as needed.
  Note the at the format of `params` must match that provided in `varexpr`. 

The recognized keyword arguments in `kw_args` are the following:
- `parameter_refs`: This is mandatory if not specified in `varexpr`. Can be a
  single parameter reference, a single parameter array with parameters defined
  in the same call of [`@infinite_parameter`](@ref),
  or a tuple where each element is either of the first two options listed.
- `func`: A concrete Julia function of the form 
  `func(paramvals..., func_args...; func_kwargs...)::Float64` where the format of 
  `param_vals` matches that of `parameter_refs` but accepts numeric supports. We could 
  also instead give an anonymous function.
- `func_args`: A `Tuple` of additional positional arguments for `func`
- `func_kwargs`: A `NameTuple` of keyword arguments for `func`
- `base_name`: Sets the name prefix used to generate object names. It
  corresponds to the object name for scalar parameter function, otherwise, the
  object names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `container`: Specify the container type.

**Examples**
```julia-repl 
julia> @parameter_function(model, func = sin, parameter_refs = t, base_name = "sin")
sin(t)

julia> @parameter_function(model, [i = 1:2], func = [sin, cos][i], parameter_refs = t)
2-element Array{GeneralVariableRef,1}:
 noname(t)
 noname(t)

julia> f(t_val, x_vals) = t_val + sum(x_vals)
f (generic function with 1 method)

julia> @parameter_function(model, pf(t, x) == f)
pf(t, x)

julia> g(t_val, a; b = 0) = t_val + a + b
g (generic function with 1 method)

julia> @parameter_function(model, pf2[i = 1:2](t) == g(t, i, b = 2 * i ))
2-element Array{GeneralVariableRef,1}:
 pf2[1](t)
 pf2[2](t)

julia> @parameter_function(model, pf2_alt[i = 1:2](t) == t -> g(t, i, b = 2 * i ))
2-element Array{GeneralVariableRef,1}:
 pf2_alt[1](t)
 pf2_alt[2](t)
```
"""
macro parameter_function(model, args...)
    # prepare the model 
    esc_model = esc(model)

    # define error message function
    _error(str...) = _macro_error(:parameter_function, (model, args...),
                                  __source__, str...)

    # parse the arguments
    extra, kwargs, requestedcontainer = _extract_kw_args(args)
    func_kwarg = filter(kw -> kw.args[1] == :func, kwargs)
    prefs_kwarg = filter(kw -> kw.args[1] == :parameter_refs, kwargs)
    fargs_kwarg = filter(kw -> kw.args[1] == :func_args, kwargs)
    fkwargs_kwarg = filter(kw -> kw.args[1] == :func_kwargs, kwargs)
    base_name_kwarg = filter(kw -> kw.args[1] == :base_name, kwargs)
    filter!(kw -> kw.args[1] != :func && kw.args[1] != :parameter_refs &&
            kw.args[1] != :base_name && kw.args[1] == :func_args &&
            kw.args[1] == :func_kwargs, kwargs)

    # if an equality syntax is given convert it into 2 separate arguments
    if length(extra) == 1 && isexpr(extra[1], :call) && extra[1].args[1] == :(==)
        extra = [extra[1].args[2], extra[1].args[3]]
    end

    # process and extract the symbolic inputs
    fargs = nothing
    fkwargs = nothing
    anon_kwargs = !isempty(func_kwarg) && !isempty(prefs_kwarg)
    if isempty(extra) && anon_kwargs
        # single anonymous definition
        var_expr = gensym()
        func_expr = _esc_non_constant(func_kwarg[1].args[2])
        pref_expr = _esc_non_constant(prefs_kwarg[1].args[2])
        anon = true
    elseif length(extra) == 1 && (isexpr(extra[1], :vect) || isexpr(extra[1], :vcat)) && anon_kwargs
        # container definition of anonymous parameter functions
        var_expr = extra[1]
        func_expr = _esc_non_constant(func_kwarg[1].args[2])
        pref_expr = _esc_non_constant(prefs_kwarg[1].args[2])
        anon = true
    elseif isempty(extra) || (length(extra) == 1 && (isexpr(extra[1], :vect) || isexpr(extra[1], :vcat)))
        # we have one of the above situations, but not enough information
        _error("Anonymous calls must specify both the `func` and " *
               "`parameter_refs` keyword arguments.")
    elseif !isempty(func_kwarg) || !isempty(prefs_kwarg)
        # anonymous keyword args cannot be mixed with symbolic syntax
        _error("Keyword arguments `func` and `parameter_refs` are only for " * 
               "anonymous parameter function definitions.")
    elseif length(extra) == 2 && isexpr(extra[1], :call)
        # 2 argument symbolic definition which may or may not be anonymous
        # anonymous if the first argument is of form [idx_expr](pref_expr)
        # otherwise the first argument should be name[idx_expr](pref_expr) or name(pref_expr)
        _ensure_not_operator(_error, extra[1].args[1])
        var_expr = extra[1].args[1]
        pref_expr = Expr(:tuple, extra[1].args[2:end]...)
        func_expr, fargs, fkwargs = _process_func_expr(_error, extra[2], pref_expr)
        pref_expr = esc(pref_expr)
        anon = isexpr(var_expr, :vect) || isexpr(var_expr, :vcat)
    elseif length(extra) == 2 
        # some other 2 argument syntax
        _error("Invalid explicit syntax, should be of form ",
               "`@parameter_function(model, nameexpr(params...), funcexpr, kwargs...)")
    elseif length(extra) > 2
        # cannot have more than 2 positional arguments
        _error("Too many positional arguments.")
    else
        # something else was given incorrectly with one positional argument
        _error("Invalid syntax.")
    end

    # parse the base name appropriately
    name = _get_name(var_expr)
    if isempty(base_name_kwarg)
        base_name = anon ? "" : string(name) # TODO perhaps change default behavior
    else
        base_name = esc(base_name_kwarg[1].args[2])
    end
    if !isa(name, Symbol) && !anon
        _error("Expression $name should not be used as a parameter function name. Use " *
               "the \"anonymous\" syntax $name = @parameter_function(model, " *
               "...) instead.")
    end

    # update and check fargs
    if !isempty(fargs_kwarg) && fargs === nothing 
        fargs = _esc_non_constant(fargs_kwarg[1].args[2])
    elseif !isempty(fargs_kwarg)
        _error("Cannot double specify the `func_args` keyword argument.")
    end

    # update and check fkwargs
    if !isempty(fkwargs_kwarg) && fkwargs === nothing 
        fkwargs = _esc_non_constant(fkwargs_kwarg[1].args[2])
    elseif !isempty(fkwargs_kwarg)
        _error("Cannot double specify the `func_kwargs` keyword argument.")
    end

    # generate the needed code 
    pfunc = gensym()
    buildcall = :( build_parameter_function($_error, $func_expr, $pref_expr, 
                                                func_args = $fargs, 
                                                func_kwargs = $fkwargs) )
    _add_kw_args(buildcall, kwargs)
    if isa(var_expr, Symbol)
        # create a single parameter function 
        pfunctioncall = :( add_parameter_function($esc_model, $buildcall, $base_name) )
        creationcode = :( $pfunc = $pfunctioncall )
    else
        # create code for a container 
        idxvars, indices = JuMPC._build_ref_sets(_error, var_expr)
        pfunctioncall = :( add_parameter_function($esc_model, $buildcall,
                                         $(_name_call(base_name, idxvars))) )
        creationcode = JuMPC.container_code(idxvars, indices, pfunctioncall,
                                            requestedcontainer)
    end
    if anon
        # anonymous creation
        macro_code = creationcode
    else
        # want to register the object/container reference
        macro_code = _macro_assign_and_return(creationcode, pfunc, name,
                                              model_for_registering = esc_model)
    end
    return _finalize_macro(esc_model, macro_code, __source__)
end
