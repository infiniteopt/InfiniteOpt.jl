# TODO add deprecation errors for old methods

################################################################################
#                               USER FUNCTIONS
################################################################################
# Keep track of the predefined functions in MOI
const _NativeNLPFunctions = append!(copy(MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS), 
                                    MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS)
append!(_NativeNLPFunctions, (:&&, :||, :<=, :(==), :>=, :<, :>))

# TODO expand details in docstring
"""
    JuMP.add_user_defined_function(
        model::InfiniteModel,
        op::Symbol,
        dim::Int,
        f::Function,
        [∇f::Function,]
        [∇²f::Function]
    )

Extend `JuMP.add_user_defined_function` for `InfiniteModel`s. 

Add a user-defined function with `dim` input arguments to `model` and associate
it with the operator `op`.

The function `f` evaluates the function. The optional function `∇f` evaluates
the first derivative, and the optional function `∇²f` evaluates the second
derivative. `∇²f` may be provided only if `∇f` is also provided.
"""
function JuMP.add_user_defined_function(
    model::InfiniteModel,
    op::Symbol,
    dim::Int,
    funcs...
    )
    if isempty(funcs)
        error("Tried to register `$op`, but no evaluation function was given.")
    elseif !all(f -> f isa Function, funcs) 
        error("The gradient and/or hessian must be functions, but got argument(s) " *
              "of type `" * join(Tuple(typeof(f) for f in funcs if !(f isa Function)), "`, `") *
              "`.")
    elseif op in _NativeNLPFunctions || op in keys(model.func_lookup)
        error("A function with name `$op` arguments is already " *
               "registered. Please use a function with a different name.")
    elseif !hasmethod(funcs[1], NTuple{dim, Real})
        error("The function `$op` is not defined for arguments of type `Real`.")
    end
    push!(model.registrations, RegisteredFunction(op, dim, funcs...))
    model.func_lookup[op] = (funcs[1], dim)
    # TODO should we set the optimizer model to be out of date?
    return JuMP.UserDefinedFunction(op)
end

"""
    name_to_function(model::InfiniteModel, op::Symbol)::Union{Function, Nothing} 

Return the registered function that corresponds to `op`. 
Returns `nothing` if no such registered function exists. This helps retrieve the 
functions of user-defined functions.

!!! warning
    Currently, this does not return functions for default operators.
"""
function name_to_function(model::InfiniteModel, op::Symbol)
    haskey(model.func_lookup, op) && return model.func_lookup[op][1]
    return
end

"""
    all_registered_functions(model::InfiniteModel)::Vector{Symbol}

Retrieve all the functions that are currently registered to `model`.
"""
function all_registered_functions(model::InfiniteModel) 
    return append!(copy(_NativeNLPFunctions), map(v -> Symbol(first(v)), values(model.func_lookup)))
end

"""
    user_registered_functions(model::InfiniteModel)::Vector{RegisteredFunction}

Return all the functions (and their associated information) that the user has 
registered to `model`. Each is stored as a [`RegisteredFunction`](@ref).
"""
function user_registered_functions(model::InfiniteModel)
    return model.registrations
end

## Define helper function to add registered functions to JuMP
# No gradient or hessian
function _add_func_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredFunction{F, Nothing, Nothing}
    ) where {F <: Function}
    JuMP.add_user_defined_function(model, data.op, data.dim, data.f)
    return
end

# Only gradient information
function _add_func_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredFunction{F, G, Nothing}
    ) where {F <: Function, G <: Function}
    JuMP.add_user_defined_function(model, data.op, data.dim, data.f, data.∇f)
    return
end

# Gradient and hessian information
function _add_func_data_to_jump(model::JuMP.Model, data::RegisteredFunction)
    JuMP.add_user_defined_function(model, data.op, data.dim, data.f, data.∇f, data.∇²f)
    return
end

"""
    add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)::Nothing

Add the user registered functions in `inf_model` to a `JuMP` model `opt_model`. 
This is intended as an internal method, but it is provided for developers that 
extend `InfiniteOpt` to use other optimizer models.
"""
function add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)
    for data in user_registered_functions(inf_model)
        _add_func_data_to_jump(opt_model, data)
    end
    return
end
