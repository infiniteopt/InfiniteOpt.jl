################################################################################
#                                BASIC HELPERS
################################################################################
# Macro error function
# inspired from https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/macros.jl#L923-L928
function _macro_error(macroname, args, source, str...)
    error("At $(source.file):$(source.line): `@$macroname($(join(args, ", ")))`: ", 
          str...)
end

# Escape when needed
# taken from https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/macros.jl#L895-L897
_esc_non_constant(x::Number) = x
_esc_non_constant(x::Expr) = isexpr(x,:quote) ? x : esc(x)
_esc_non_constant(x) = esc(x)

# Extract the name from a macro expression 
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/45ce630b51fb1d72f1ff8fed35a887d84ef3edf7/src/Containers/macro.jl#L8-L17
_get_name(c::Symbol) = c
_get_name(c::Nothing) = ()
_get_name(c::AbstractString) = c
function _get_name(c::Expr)
    if isexpr(c, :string)
        return c
    else
        return c.args[1]
    end
end

# Given a base_name and idxvars, returns an expression that constructs the name
# of the object.
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/macros.jl#L930-L946
function _name_call(base_name, idxvars)
    if isempty(idxvars) || base_name == ""
        return base_name
    end
    ex = Expr(:call, :string, base_name, "[")
    for i in eachindex(idxvars)
        # Converting the arguments to strings before concatenating is faster:
        # https://github.com/JuliaLang/julia/issues/29550.
        esc_idxvar = esc(idxvars[i])
        push!(ex.args, :(string($esc_idxvar)))
        i < length(idxvars) && push!(ex.args, ",")
    end
    push!(ex.args, "]")
    return ex
end

# Process macro arugments 
function _extract_kwargs(args)
    arg_list = collect(args)
    if !isempty(args) && isexpr(args[1], :parameters)
        p = popfirst!(arg_list)
        append!(arg_list, p.args)
    end
    extra_kwargs = filter(x -> isexpr(x, :(=)) && x.args[1] != :container &&
                          x.args[1] != :base_name, arg_list)
    container_type = :Auto
    base_name = nothing
    for kw in arg_list
        if isexpr(kw, :(=)) && kw.args[1] == :container
            container_type = kw.args[2]
        elseif isexpr(kw, :(=)) && kw.args[1] == :base_name
            base_name = esc(kw.args[2])
        end
    end
    pos_args = filter!(x -> !isexpr(x, :(=)), arg_list)
    return pos_args, extra_kwargs, container_type, base_name
end

# Add on keyword arguments to a function call expression and escape the expressions
# Adapted from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L11-L36
function _add_kwargs(call, kwargs)::Nothing
    for kw in kwargs
        push!(call.args, esc(Expr(:kw, kw.args...)))
    end
    return
end

# Determine if an expression contains any index variable symbols
function _has_idxvars(expr, idxvars)::Bool
    expr in idxvars && return true
    if expr isa Expr
        return any(_has_idxvars(a, idxvars) for a in expr.args)
    end
    return false
end

# Ensure a model argument is valid
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L38-L44
_valid_model(_error::Function, model::InfiniteModel, name) = nothing
function _valid_model(_error::Function, model, name)
    _error("Expected $name to be an `InfiniteModel`, but it has type ", 
           typeof(model))
end

# Check if a macro julia variable can be registered 
# Adapted from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L66-L86
function _error_if_cannot_register(
    _error::Function, 
    model::InfiniteModel, 
    name::Symbol
    )::Nothing
    if haskey(JuMP.object_dictionary(model), name)
        _error("An object of name $name is already attached to this model. If ",
               "this is intended, consider using the anonymous construction ",
               "syntax, e.g., `x = @macro_name(model, ...)` where the name ",
               "of the object does not appear inside the macro. Alternatively, ",
               "use `unregister(model, :$(name))` to first unregister the ",
               "existing name from the model. Note that this will not delete ",
               "the object; it will just remove the reference at ",
               "`model[:$(name)]`")
    end
    return
end
function _error_if_cannot_register(_error::Function, model::InfiniteModel, name)
    return _error("Invalid name $name.")
end

# Update the creation code to register and assign the object to the name
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L88-L120
function _macro_assign_and_return(_error, code, name, model)
    return quote
        _error_if_cannot_register($_error, $model, $(quot(name)))
        $(esc(name)) = $code
        $model[$(quot(name))] = $(esc(name))
    end
end

# Wrap the macro generated code for better stacttraces (assumes model is escaped)
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L46-L64
function _finalize_macro(_error, model, code, source::LineNumberNode)
    return Expr(:block, source, 
                :(_valid_model($_error, $model, $(quot(model.args[1])))), code)
end

################################################################################
#                          INFINITE PARAMETER MACRO
################################################################################
## Process a distribution input from a macro into a distribution domain
# Univariate Distribution
function _distribution_or_error(
    _error::Function, 
    dist::D
    )::UniDistributionDomain{D} where {D <: Distributions.UnivariateDistribution}
    return UniDistributionDomain(dist)
end

# Multivariate Distribution
function _distribution_or_error(
    _error::Function, 
    dist::D
    )::MultiDistributionDomain{D} where {D <: NonUnivariateDistribution}
    return MultiDistributionDomain(dist)
end

# Fallback
function _distribution_or_error(_error::Function, dist)
    _error("Expected distribution from `Distributions.jl`, but got input of ",
           "type `$(typeof(dist))`.")
end

## Process an infinite domain macro input and make sure it is correct 
# AbstractInfiniteDomain
function _domain_or_error(
    _error::Function, 
    domain::D
    )::D where {D <: AbstractInfiniteDomain}
    return domain
end

# Vector{Real}
function _domain_or_error(
    _error::Function, 
    vect::Vector{<:Real}
    )::IntervalDomain
    if length(vect) == 2
        return IntervalDomain(vect...)
    else
        _error("Expected interval domain of format `[lb, ub]`, but got `$vect`.")
    end
end

# Distribution 
function _domain_or_error(_error::Function, dist::Distributions.Distribution)
    _error("Distribution `$dist` was specified as a domain, but should " * 
          "be given as a distribution via the syntax: `@infinite_parameter( " * 
          "model, param_expr ~ distribution, kwargs...)`.")
end

# Fallback
function _domain_or_error(_error::Function, domain)
    _error("Expected an infinite domain, but got input of ",
           "type `$(typeof(domain))`.")
end

# Process an expr as a distribution
function _make_distribution_call(_error, expr)
    return :( _distribution_or_error($_error, $(_esc_non_constant(expr))) )
end

# Process an expr as an infinite domain
function _make_domain_call(_error, expr)
    return :( _domain_or_error($_error, $(_esc_non_constant(expr))) )
end

# Parse a distribution input
function _parse_parameter(_error::Function, ::Val{:~}, param, dist)
    return param, _make_distribution_call(_error, dist)
end

# Parse a domain input
function _parse_parameter(
    _error::Function, 
    ::Union{Val{:in}, Val{:∈}}, 
    param,
    domain
    )
    return param, _make_domain_call(_error, domain)
end

# Fallback
function _parse_parameter(_error::Function, ::Val{S}, param, rhs) where S
    _error("Unexpected operator $S.")
end

# Define the keyword arg types that can vary for individual dependent parameters
const _ArrayKwargNames = (:supports, :derivative_method)

"""
    @infinite_parameter(model::InfiniteModel, kwargs...)

Add *anonymous* infinite parameter to the model `model` described by
the keyword arguments `kw_args` and returns the parameter reference.

```julia
@infinite_parameter(model::InfiniteModel, expr, kwargs...)
```

Add an infinite parameter to the model `model` described by the expression `expr`,
and the keyword arguments `kw_args`. (Note that in the following the symbol `in` 
can be used instead of `∈`) The expression `expr` can be of the form:
- `paramexpr` creating parameters described by `paramexpr`.
- `paramexpr ∈ [lb, ub]` creating parameters described by `paramexpr` characterized
   by a continuous interval domain with lower bound `lb` and upper bound `ub`.
- `paramexpr ~ dist` creating parameters described by `paramexpr` characterized
   by the `Distributions.jl` distribution object `dist`.
- `paramexpr ∈ domain` creating parameters described by `paramexpr` characterized
   by the `AbstractInfiniteDomain` object `domain`.

The expression `paramexpr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The recognized keyword arguments in `kwargs` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `domain`: The `InfiniteDomain` characterizing the parameters see subtypes of
            [`AbstractInfiniteDomain`](@ref).
- `distribution`: Sets the `Distributions.jl` distribution object that 
   characterizes the parameters (specified instead of a domain).
- `supports`: Sets the support points for the parameters.
- `num_supports`: Specifies the number of supports to be automatically generated.
                  Note that `supports` takes precedence. Defaults to 0.
- `derivative_method`: Specify the numerical method to evaluate derivatives that 
                       are taken with respect to the parameter.
- `sig_digits`: Specifies the number of significant digits that should be used
                in automatic support generation. Defaults to `DefaultSigDigits`.
- `independent`: Specifies if the each parameter is independent from each other
  or not. Defaults to false.
- `container`: Specify the container type. Defaults to `Auto`

**Examples**
```julia-repl
julia> @infinite_parameter(m, x in [0, 1])
x

julia> @infinite_parameter(m, y[i = 1:2] ~ MvNormal(ones(2)), num_supports = 10)
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
    pos_args, extra_kwargs, container_type, base_name = _extract_kwargs(args)

    # process the keyword arguments
    domain_kwarg = filter(kw -> kw.args[1] == :domain, extra_kwargs)
    dist_kwarg = filter(kw -> kw.args[1] == :distribution, extra_kwargs)
    indep_kwarg = filter(kw -> kw.args[1] == :independent, extra_kwargs)
    filter!(kw -> !(kw.args[1] in (:domain, :independent, :distribution)), 
            extra_kwargs)

    # process the positional arguments
    if isempty(pos_args)
        p = gensym()
        is_anon_single = true
    elseif length(pos_args) == 1
        p = first(pos_args)
        is_anon_single = false
    else
        _error("Too many positional arguments. See docs for correct syntax.")
    end

    # extract the domain if given symbolically 
    # There are 9 cases to consider:
    # p                                  | type of p | p.head
    # -----------------------------------+-----------+------------
    # name                               | Symbol    | NA
    # name[1:2]                          | Expr      | :ref
    # name[i = 1:2, j = 1:2; i + j >= 3] | Expr      | :typed_vcat
    # [1:2]                              | Expr      | :vect
    # [i = 1:2, j = 1:2; i + j >= 3]     | Expr      | :vcat
    # name in domain                     | Expr      | :call
    # name[...] in domain                | Expr      | :call
    # name ~ dist                        | Expr      | :call
    # name[...] ~ dist                   | Expr      | :call
    # In the 4 last cases, we call _parse_parameter
    if isexpr(p, :call)
        param, domain = _parse_parameter(_error, Val(p.args[1]), p.args[2:end]...)
        if !isempty(domain_kwarg) || !isempty(dist_kwarg)
            _error("Cannot double specify the infinite domain and/or ",
                   "distribution.")
        end
    elseif !isempty(domain_kwarg) && isempty(dist_kwarg)
        param = p 
        domain = _make_domain_call(_error, domain_kwarg[1].args[2])
    elseif isempty(domain_kwarg) && !isempty(dist_kwarg)
        param = p 
        domain = _make_distribution_call(_error, dist_kwarg[1].args[2])
    elseif !isempty(domain_kwarg) && !isempty(dist_kwarg)
        _error("Cannot specify both a domain and a distribution.")
    else
        _error("Must specify the infinite domain and/or the given syntax is ",
               "unrecognized. See docs for accepted forms.")
    end

    # determine if it is a single parameter
    is_single = param isa Symbol 

    # determine if independent 
    is_independent = false
    if is_single 
        is_independent = true
    elseif !isempty(indep_kwarg)
        if !isa(indep_kwarg[1].args[2], Bool)
            _error("Can only specify boolean literals (`false` or `true`) ", 
                   "with keyword `independent`.")
        end
        is_independent = indep_kwarg[1].args[2] 
    end

    # determine if anonymous 
    is_anon = isexpr(param, (:vect, :vcat)) || is_anon_single

    # process the naming 
    name = _get_name(param)
    if !isa(name, Symbol) && !is_anon 
        _error("Expression $name should not be used as a parameter name. Use ",
               "the \"anonymous\" syntax ",
               "`$name = @infinite_parameter(model, kwargs...)`` instead.")
    end
    if base_name === nothing
        base_name = is_anon ? "" : string(name)
    end

    # prepare the name code if it is multi-dimensional
    if !is_single
        idxvars, inds = JuMPC.build_ref_sets(_error, param)
        name_code = _name_call(base_name, idxvars)
        if model in idxvars
            _error("Index $(model) is the same symbol as the model. Use a ",
                   "different name for the index.")
        end
    end

    # make the build call
    if is_independent
        # we only have a single parameter or an array of independent parameters
        build_call = :( build_parameter($_error, $domain) )
    else
        # we have multi-dimensional dependent parameters
        # let's first process the kwargs accordingly (array ones are made into 
        # containers)
        array_kwargs = filter(kw -> kw.args[1] in _ArrayKwargNames, extra_kwargs)
        filter!(kw -> !(kw.args[1] in _ArrayKwargNames), extra_kwargs)
        # check that the non array keywords do not contain index variables used 
        # for building containers
        for kw in extra_kwargs
            if _has_idxvars(kw.args[2], idxvars)
                _error("Cannot use index variables with `$(kw.args[1])` keyword.")
            end
        end
        # now let's make the build call
        inds_var = gensym() # used as a placeholder variable for the ContainerIndices used to vectorize
        domain_var = gensym() # used as a placeholder for the domain container
        vect_domain_call = :( Collections.vectorize($domain_var, $inds_var) )
        build_call = :( _build_parameters($_error, $vect_domain_call, $inds_var) )
        for kw in array_kwargs
            code = JuMPC.container_code(idxvars, inds, 
                                        _esc_non_constant(kw.args[2]), 
                                        container_type)
            vect_code = :( Collections.vectorize($code, $inds_var) )
            push!(build_call.args, Expr(:kw, kw.args[1], vect_code))
        end
    end
    _add_kwargs(build_call, extra_kwargs)

    # make the creation code
    esc_model = esc(model)
    if is_single 
        # only single independent parameter
        creation_code = :( add_parameter($esc_model, $build_call, $base_name) )
    elseif is_independent
        # multi-dimensional independent parameters 
        parameter_call = :( add_parameter($esc_model, $build_call, $name_code) )
        creation_code = JuMPC.container_code(idxvars, inds, parameter_call,
                                             container_type)
    else
        # multi-dimensional dependent parameters
        domain_call = JuMPC.container_code(idxvars, inds, domain, container_type)
        name_code = JuMPC.container_code(idxvars, inds, name_code, container_type)
        name_code = :( Collections.vectorize($name_code, $inds_var) )
        vect_add_call = :( add_parameters($esc_model, $build_call, $name_code) )
        creation_code = quote 
            $domain_var = $domain_call
            $inds_var = Collections.indices($domain_var)
            Collections.unvectorize($vect_add_call, $inds_var)
        end
    end

    # finalize the code and register the parameter name if needed
    if is_anon
        macro_code = creation_code
    else
        macro_code = _macro_assign_and_return(_error, creation_code, name,
                                              esc_model)
    end
    return _finalize_macro(_error, esc_model, macro_code, __source__)
end

################################################################################
#                           FINITE PARAMETER MACRO
################################################################################

"""
    @finite_parameter(model::InfiniteModel, value, kwargs...)

Define and add an anonymous finite parameter to `model` and return its
parameter reference. Its value is equal to `value`.

```julia
@finite_parameter(model::InfiniteModel, param_expr == value_expr, kwargs...)
```
Define and add a finite parameter(s) to `model` and return appropriate parameter
reference(s). The parameter(s) has/have value(s) indicated by the `value_expr`.
The expression `param_expr` can be of the form:
- `paramname` creating a scalar parameter of name `paramname`
- `paramname[...]` or `[...]` creating a container of parameters

The expression `value_expr` simply expresses the value of the parameter(s). This
is typically a number but could be an array indexed using an index defined in
`param_expr`.

The recognized keyword arguments in `kwargs` are the following:
- `base_name`: Sets the name prefix used to generate parameter names. It
  corresponds to the parameter name for scalar parameter, otherwise, the
  parameter names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `container`: Specify the container type, defaults to `Auto`.

**Examples**
```julia-repl
julia> par = @finite_parameter(model, 2)
noname

julia> vals = [3, 2];

julia> pars = @finite_parameter(model, [i = 1:2] == vals[i], base_name = "par")
2-element Array{ParameterRef,1}:
 par[1]
 par[2]

julia> @finite_parameter(model, par2 == 42)
par2
```
"""
macro finite_parameter(model, args...)
    # process the inputs
    esc_model = esc(model)
    pos_args, kwargs, container_type, base_name = _extract_kwargs(args)

    # make an error function
    _error(str...) = _macro_error(:finite_parameter, (model, args...), 
                                  __source__, str...)

    # process the positional arguments
    if length(pos_args) == 1
        expr = popfirst!(pos_args)
    else
        _error("Incorrect number of arguments. Must be of form ",
               "@finite_parameter(model, name_expr == value).")
    end
    if isexpr(expr, :call)
        if expr.args[1] !== :(==)
            _error("Unrecognized operator. Must be of form ",
                   "@finite_parameter(model, name_expr == value).")
        end
        param = expr.args[2]
        value = _esc_non_constant(expr.args[3])
        is_anon = isexpr(param, (:vect, :vcat))
    else
        param = gensym()
        value = _esc_non_constant(expr)
        is_anon = true
    end

    # process the name
    name = _get_name(param)
    if base_name === nothing
        base_name = is_anon ? "" : string(name)
    end
    if !isa(name, Symbol) && !is_anon
        _error("Expression $name should not be used as a parameter name. Use " *
               "the \"anonymous\" syntax $name = @finite_parameter(model, " *
               "...) instead.")
    end

    # make the build call 
    build_call = :( build_parameter($_error, $value) )
    _add_kwargs(build_call, kwargs)

    # make the creation code
    if isa(param, Symbol)
        # easy case with single parameter
        creation_code = :( add_parameter($esc_model, $build_call, $base_name) )
    else
        # we have a container of parameters
        idxvars, inds = JuMPC.build_ref_sets(_error, param)
        if model in idxvars
            _error("Index $(model) is the same symbol as the model. Use a ",
                   "different name for the index.")
        end
        name_code = _name_call(base_name, idxvars)
        parameter_call = :( add_parameter($esc_model, $build_call, $name_code) )
        creation_code = JuMPC.container_code(idxvars, inds, parameter_call,
                                             container_type)
    end
    
    # finalize the macro
    if is_anon
        macro_code = creation_code
    else
        macro_code = _macro_assign_and_return(_error, creation_code, name,
                                              esc_model)
    end
    return _finalize_macro(_error, esc_model, macro_code, __source__)
end

################################################################################
#                         PARAMETER FUNCTION MACRO
################################################################################
const _BadOperators = (:in, :<=, :>=, :(==), :≥, :≤, ∈)

## Make methods to find and replace a part of an expression
# Expr
function _expr_replace!(ex::Expr, old, new)
    for (i, a) in enumerate(ex.args)
        if a == old
            ex.args[i] = new
        elseif a isa Expr
            _expr_replace!(a, old, new)
        end
    end
    return ex
end

# Symbol
function _expr_replace!(ex::Symbol, old, new)
    return ex == old ? new : ex
end

## Safely extract the parameter references (copy as needed)
# Symbol 
function _extract_parameters(ex::Symbol)
    return esc(ex)
end

# Expr 
function _extract_parameters(ex::Expr) 
    return esc(copy(ex))
end

# Helper method to process parameter function expressions 
function _process_func_expr(_error::Function, raw_expr)
    if isexpr(raw_expr, :call)
        # check that the call is not some operator
        if raw_expr.args[1] in _BadOperators
            _error("Invalid input syntax.")
        end
        func_expr = esc(raw_expr.args[1])
        is_anon = false
        # check for keywords
        if isexpr(raw_expr.args[2], :parameters) || 
            any(isexpr(a, :kw) for a in raw_expr.args[2:end])
            _error("Cannot specify keyword arguements directly, try using an ",
                   "anonymous function.")
        end
        # extract the parameter inputs
        pref_expr = esc(Expr(:tuple, raw_expr.args[2:end]...)) 
    elseif isexpr(raw_expr, :(->))
        # extract the parameter inputs
        pref_expr = _extract_parameters(raw_expr.args[1]) # this will create a copy if needed
        # fix the function if the parameter arguments were given as references
        if isexpr(raw_expr.args[1], :ref) 
            raw_expr = _expr_replace!(raw_expr, raw_expr.args[1], gensym())
        elseif isexpr(raw_expr.args[1], :tuple)
            for a in raw_expr.args[1].args
                if isexpr(a, :ref)
                    raw_expr = _expr_replace!(raw_expr, a, gensym())
                end
            end
        end
        # extract the function expression
        func_expr = esc(raw_expr)
        is_anon = true
    else 
        _error("Unrecognized syntax.")
    end
    return func_expr, pref_expr, is_anon
end

"""
    @parameter_function(model::InfiniteModel, func_expr, kwargs...)

Add an *anonymous* parameter function to the model `model` described by the
keyword arguments `kw_args` and returns the object reference.

```julia
@parameter_function(model::InfiniteModel, var_expr == func_expr, kwargs...)
```

Add a parameter function to `model` described by the expression `var_expr`, the
function expression `func_expr`, and the keyword arguments `kwargs`. The 
expression `var_expr` is used to define the parameter function references of the 
form `varname[...]` where the indexing matches the container syntax of other 
macros.

The expression `func_expr` determines the concrete Julia function that defines the 
behavior of parameter function and also specifies the infinite parameters it 
depends on. The accepted forms are:
- `func(params...)`: where `func` is the function that takes supports of the 
  infinite parameters `params` as input and outputs a scalar value.
- `(params...) -> my_func_expr`: where `params` are the infinite parameters and 
  `my_func_expr` is the source code of the anonymous function.

The recognized keyword arguments in `kwargs` are the following:
- `base_name`: Sets the name prefix used to generate object names. It
  corresponds to the object name for scalar parameter function, otherwise, the
  object names are set to `base_name[...]` for each index `...` of the axes
  `axes`.
- `container`: Specify the container type. Defaults to `:Auto`.

**Examples**
```julia-repl 
julia> @parameter_function(model, sin(t))
sin(t)

julia> func_vect = [sin, cos];

julia> @parameter_function(model, [i = 1:2] == func_vect[i](t))
2-element Array{GeneralVariableRef,1}:
 sin(t)
 cos(t)

julia> f(t_val, x_vals) = t_val + sum(x_vals)
f (generic function with 1 method)

julia> @parameter_function(model, pf == f(t, x))
pf(t, x)

julia> g(t_val, a; b = 0) = t_val + a + b
g (generic function with 1 method)

julia> @parameter_function(model, pf2[i = 1:2] == t -> g(t, i, b = 2 * i ))
2-element Array{GeneralVariableRef,1}:
 pf2[1](t)
 pf2[2](t)
```
"""
macro parameter_function(model, args...)
    # prepare the model 
    esc_model = esc(model)

    # define error message function
    _error(str...) = _macro_error(:parameter_function, (model, args...),
                                  __source__, str...)

    # parse the arguments 
    pos_args, kwargs, container_type, base_name = _extract_kwargs(args)

    # process the positional arguements 
    if length(pos_args) == 1
        expr = popfirst!(pos_args)
    else
        _error("Incorrect amount of positional arguments given.")
    end
    if isexpr(expr, :call) && expr.args[1] === :(==)
        var = expr.args[2]
        func, prefs, is_anon_func = _process_func_expr(_error, expr.args[3])
        is_anon = isexpr(var, (:vect, :vcat))
    elseif isexpr(expr, (:call, :(->)))
        var = gensym()
        func, prefs, is_anon_func = _process_func_expr(_error, expr)
        is_anon = true
    else
        _error("Unrecognized syntax.")
    end

    # parse the name appropriately
    name = _get_name(var)
    if !isa(name, Symbol) && !is_anon
        _error("Expression $name should not be used as a parameter function name. Use " *
               "the \"anonymous\" syntax $name = @parameter_function(model, " *
               "...) instead.")
    end
    idxvars, inds = JuMPC.build_ref_sets(_error, var)
    if model in idxvars
        _error("Index $(model) is the same symbol as the model. Use a ",
               "different name for the index.")
    end
    if base_name === nothing && is_anon && !is_anon_func
        name_code = :( string(nameof($func)) )
    elseif base_name === nothing && is_anon
        name_code = _name_call("", idxvars)
    elseif base_name === nothing
        name_code = _name_call(string(name), idxvars)
    else
        name_code = _name_call(base_name, idxvars)
    end

    # make the build call
    build_call = :( build_parameter_function($_error, $func, $prefs) )
    _add_kwargs(build_call, kwargs)

    # make the creation code
    if isa(var, Symbol)
        # create a single parameter function 
        creation_code = :( add_parameter_function($esc_model, $build_call, $name_code) )
    else
        # create code for a container 
        pfunction_call = :( add_parameter_function($esc_model, $build_call,
                                                   $name_code) )
        creation_code = JuMPC.container_code(idxvars, inds, pfunction_call,
                                             container_type)
    end

    # finalize the macro
    if is_anon
        # anonymous creation
        macro_code = creation_code
    else
        # want to register the object/container reference
        macro_code = _macro_assign_and_return(_error, creation_code, name,
                                              esc_model)
    end
    return _finalize_macro(_error, esc_model, macro_code, __source__)
end
