################################################################################
#                               USER OPERATORS
################################################################################
# Keep track of the predefined functions in MOI
const _NativeNLPOperators = append!(copy(MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS),
                                    MOI.Nonlinear.DEFAULT_MULTIVARIATE_OPERATORS)
append!(_NativeNLPOperators, (:&&, :||, :<=, :(==), :>=, :<, :>))

"""
    JuMP.add_nonlinear_operator(
        model::InfiniteModel,
        dim::Int,
        f::Function,
        [∇f::Function,]
        [∇²f::Function];
        [name::Symbol = Symbol(f)]
    )

Extend `add_nonlinear_operator` for `InfiniteModel`s.

Add a new nonlinear operator with `dim` input arguments to `model` and associate
it with the name `name`. Alternatively, [`@operator`](https://jump.dev/JuMP.jl/v1/api/JuMP/#@operator)
can be used for a more convenient syntax.

The function `f` evaluates the operator. The optional function `∇f` evaluates
the first derivative, and the optional function `∇²f` evaluates the second
derivative. `∇²f` may be provided only if `∇f` is also provided.

```julia-repl
julia> @variable(model, y);

julia> g(x) = x^2;

julia> new_op = add_nonlinear_operator(model, 1, g)
NonlinearOperator(g, :g)

julia> @expression(model, new_op(y))
g(y)
```
"""
function JuMP.add_nonlinear_operator(
    model::InfiniteModel,
    dim::Int,
    f::Function,
    funcs::Vararg{Function, N};
    name::Symbol = Symbol(f)
    ) where {N}
    if name in _NativeNLPOperators || name in keys(model.op_lookup)
        error("An operator with name `$name` arguments is already " *
               "added. Please use a operator with a different name.")
    elseif !hasmethod(f, NTuple{dim, Float64})
        error("The operator evaluation function `$f` is not defined for arguments of type `Float64`.")
    end
    push!(model.operators, NLPOperator(name, dim, f, funcs...))
    model.op_lookup[name] = (f, dim)
    # TODO should we set the transformation backend to be out of date?
    return JuMP.NonlinearOperator(f, name)
end

"""
    name_to_operator(model::InfiniteModel, name::Symbol)::Union{Function, Nothing}

Return the nonlinear operator that corresponds to `name`.
Returns `nothing` if no such operator exists.

!!! warning
    Currently, this does not return functions for default operators.
"""
function name_to_operator(model::InfiniteModel, name::Symbol)
    haskey(model.op_lookup, name) && return model.op_lookup[name][1]
    return
end

"""
    all_nonlinear_operators(model::InfiniteModel)::Vector{Symbol}

Retrieve all the operators that are currently added to `model`.
"""
function all_nonlinear_operators(model::InfiniteModel)
    return append!(copy(_NativeNLPOperators), map(v -> Symbol(first(v)), values(model.op_lookup)))
end

"""
    user_defined_operators(model::InfiniteModel)::Vector{NLPOperator}

Return all the operators (and their associated information) that the user has
added to `model`. Each is stored as a [`NLPOperator`](@ref).
"""
function added_nonlinear_operators(model::InfiniteModel)
    return model.operators
end

## Define helper function to add nonlinear operators to JuMP
# No gradient or hessian
function _add_op_data_to_jump(
    model::JuMP.Model,
    data::NLPOperator{F, Nothing, Nothing}
    ) where {F <: Function}
    JuMP.add_nonlinear_operator(model, data.dim, data.f, name = data.name)
    return
end

# Only gradient information
function _add_op_data_to_jump(
    model::JuMP.Model,
    data::NLPOperator{F, G, Nothing}
    ) where {F <: Function, G <: Function}
    JuMP.add_nonlinear_operator(model, data.dim, data.f, data.∇f, name = data.name)
    return
end

# Gradient and hessian information
function _add_op_data_to_jump(model::JuMP.Model, data::NLPOperator)
    JuMP.add_nonlinear_operator(model, data.dim, data.f, data.∇f, data.∇²f, name = data.name)
    return
end

"""
    add_operators_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)::Nothing

Add the additional nonlinear operators in `inf_model` to a `JuMP` model `opt_model`.
This is intended as an internal method, but it is provided for developers that
extend `InfiniteOpt` to use new [`JuMPBackend`](@ref)s.
"""
function add_operators_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)
    for data in added_nonlinear_operators(inf_model)
        _add_op_data_to_jump(opt_model, data)
    end
    return
end
