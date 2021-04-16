################################################################################
#                                BASIC METHODS
################################################################################
using Base.Meta

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

# Properly convert macro keywords into a keyword argument format 
# Taken from https://github.com/jump-dev/JuMP.jl/blob/45ce630b51fb1d72f1ff8fed35a887d84ef3edf7/src/variables.jl#L66-L70
function _keywordify(kw::Expr)
    return (kw.args[1], _esc_non_constant(kw.args[2]))
end

# Extract the name from a macro expression 
# Taken from https://github.com/jump-dev/JuMP.jl/blob/45ce630b51fb1d72f1ff8fed35a887d84ef3edf7/src/Containers/macro.jl#L8-L17
_get_name(c::Symbol) = c
_get_name(c::Nothing) = ()
_get_name(c::AbstractString) = c
function _get_name(c::Expr)
    if c.head == :string
        return c
    else
        return c.args[1]
    end
end

# Given a base_name and idxvars, returns an expression that constructs the name
# of the object.
# Taken from https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/macros.jl#L930-L946
function _name_call(base_name, idxvars)
    if isempty(idxvars) || base_name == ""
        return base_name
    end
    ex = Expr(:call, :string, base_name, "[")
    for i in 1:length(idxvars)
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
# Taken from https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/Containers/macro.jl#L26-L36
function _extract_kw_args(args)
    kw_args = filter(x -> isexpr(x, :(=)) && x.args[1] != :container , collect(args))
    flat_args = filter(x->!isexpr(x, :(=)), collect(args))
    requested_container = :Auto
    for kw in args
        if isexpr(kw, :(=)) && kw.args[1] == :container
            requested_container = kw.args[2]
        end
    end
    return flat_args, kw_args, requested_container
end

# Add on keyword arguments to a function call expression and escape the expressions
# Taken from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L11-L36
function _add_kw_args(call, kw_args)
    for kw in kw_args
        @assert isexpr(kw, :(=)) "Unrecognized keyword argument format."
        push!(call.args, esc(Expr(:kw, kw.args...)))
    end
end

# Ensure a model argument is valid
# Inspired from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L38-L44
_valid_model(m::InfiniteModel, name) = nothing
function _valid_model(m, name)
    error("Expected $name to be an `InfiniteModel`, but it has type ", typeof(m))
end

# Check if a macro julia variable can be registered 
# Adapted from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L66-L86
function _error_if_cannot_register(model::InfiniteModel, name::Symbol)
    obj_dict = JuMP.object_dictionary(model)
    if haskey(obj_dict, name)
        error(
            """An object of name $name is already attached to this model. If this
          is intended, consider using the anonymous construction syntax, e.g.,
          `x = @finite_variable(model, [1:N], ...)` where the name of the object does
          not appear inside the macro.
          Alternatively, use `unregister(model, :$(name))` to first unregister
          the existing name from the model. Note that this will not delete the
          object; it will just remove the reference at `model[:$(name)]`.
      """,
        )
    end
    return
end
function _error_if_cannot_register(model::InfiniteModel, name)
    return error("Invalid name $name.")
end

# Update the creation code to register and assign the object to the name
# Taken from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L88-L120
function _macro_assign_and_return(code, variable, name; 
                                  model_for_registering = nothing)
    return quote
        $(
            if model_for_registering !== nothing
                :(_error_if_cannot_register(
                    $model_for_registering,
                    $(quot(name)),
                ))
            end
        )
        $variable = $code
        $(
            if model_for_registering !== nothing
                :($model_for_registering[$(quot(name))] = $variable)
            end
        )
        # This assignment should be in the scope calling the macro
        $(esc(name)) = $variable
    end
end

# Wrap the macro generated code for better stacttraces (assumes model is escaped)
# Taken from https://github.com/jump-dev/JuMP.jl/blob/d9cd5fb16c2d0a7e1c06aa9941923492fc9a28b5/src/macros.jl#L46-L64
function _finalize_macro(model, code, source::LineNumberNode)
    return Expr(:block, source, :(_valid_model($model, $(quot(model.args[1])))),
                code)
end
