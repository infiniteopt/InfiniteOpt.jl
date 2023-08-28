# TODO add deprecation errors for old methods

################################################################################
#                               USER OPERATORS
################################################################################
# Keep track of the predefined functions in MOI
const _NativeNLPOperators = append!(copy(MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS), 
                                    MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS)
append!(_NativeNLPOperators, (:&&, :||, :<=, :(==), :>=, :<, :>))

# TODO expand details in docstring
"""
    JuMP.register_nonlinear_operator(
        model::InfiniteModel,
        dim::Int,
        f::Function,
        [∇f::Function,]
        [∇²f::Function];
        [name::Symbol = Symbol(f)]
    )

Extend `register_nonlinear_operator` for `InfiniteModel`s. 

Add a new nonlinear operator with `dim` input arguments to `model` and associate
it with the name `name`.

The function `f` evaluates the operator. The optional function `∇f` evaluates
the first derivative, and the optional function `∇²f` evaluates the second
derivative. `∇²f` may be provided only if `∇f` is also provided.
"""
function JuMP.register_nonlinear_operator(
    model::InfiniteModel,
    dim::Int,
    funcs::Vararg{Function, N};
    name::Symbol = Symbol(f)
    ) where {N}
    if isempty(funcs)
        error("Tried to register `$name`, but no evaluation function was given.")
    elseif !all(f -> f isa Function, funcs) 
        error("The gradient and/or hessian must be functions, but got argument(s) " *
              "of type `" * join(Tuple(typeof(f) for f in funcs if !(f isa Function)), "`, `") *
              "`.")
    elseif name in _NativeNLPOperators || name in keys(model.op_lookup)
        error("An operator with name `$op` arguments is already " *
               "registered. Please use a operator with a different name.")
    elseif !hasmethod(funcs[1], NTuple{dim, Float64})
        error("The operator `$op` is not defined for arguments of type `Float64`.")
    end
    push!(model.registrations, RegisteredOperator(name, dim, funcs...))
    model.op_lookup[name] = (funcs[1], dim)
    # TODO should we set the optimizer model to be out of date?
    return JuMP.NonlinearOperator(name, funcs[1])
end

"""
    name_to_operator(model::InfiniteModel, name::Symbol)::Union{Function, Nothing} 

Return the registered operator that corresponds to `name`. 
Returns `nothing` if no such registered operator exists. This helps retrieve the 
functions of user-defined nonlinear operators.

!!! warning
    Currently, this does not return functions for default operators.
"""
function name_to_operator(model::InfiniteModel, name::Symbol)
    haskey(model.op_lookup, name) && return model.op_lookup[op][1]
    return
end

"""
    all_registered_operators(model::InfiniteModel)::Vector{Symbol}

Retrieve all the operators that are currently registered to `model`.
"""
function all_registered_operators(model::InfiniteModel) 
    return append!(copy(_NativeNLPOperators), map(v -> Symbol(first(v)), values(model.op_lookup)))
end

"""
    user_registered_operators(model::InfiniteModel)::Vector{RegisteredOperator}

Return all the operators (and their associated information) that the user has 
registered to `model`. Each is stored as a [`RegisteredOperator`](@ref).
"""
function user_registered_operators(model::InfiniteModel)
    return model.registrations
end

## Define helper function to add registered operators to JuMP
# No gradient or hessian
function _add_op_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredOperator{F, Nothing, Nothing}
    ) where {F <: Function}
    JuMP.register_nonlinear_operator(model, data.dim, data.f, name = data.name)
    return
end

# Only gradient information
function _add_op_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredOperator{F, G, Nothing}
    ) where {F <: Function, G <: Function}
    JuMP.register_nonlinear_operator(model, data.dim, data.f, data.∇f, name = data.name)
    return
end

# Gradient and hessian information
function _add_op_data_to_jump(model::JuMP.Model, data::RegisteredOperator)
    JuMP.register_nonlinear_operator(model, data.dim, data.f, data.∇f, data.∇²f, name = data.name)
    return
end

"""
    add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)::Nothing

Add the user registered nonlinear operators in `inf_model` to a `JuMP` model `opt_model`. 
This is intended as an internal method, but it is provided for developers that 
extend `InfiniteOpt` to use other optimizer models.
"""
function add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)
    for data in user_registered_operators(inf_model)
        _add_op_data_to_jump(opt_model, data)
    end
    return
end
