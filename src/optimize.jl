################################################################################
#                              CORE BACKEND API
################################################################################
"""
    transformation_backend_ready(model::InfiniteModel)::Bool

Return `Bool` if the transformation backend model is up-to-date with `model` and 
ready to be optimized.
**Example**
```julia-repl
julia> transformation_backend_ready(model)
false
```
"""
transformation_backend_ready(model::InfiniteModel) = model.ready_to_optimize

"""
    set_transformation_backend_ready(model::InfiniteModel, status::Bool)

Set the status of the transformation backend model to whether it is up-to-date or 
not. Note is more intended as an internal function, but is useful for extensions.

**Example**
```julia-repl
julia> set_transformation_backend_ready(model, true)

julia> transformation_backend_ready(model)
true
```
"""
function set_transformation_backend_ready(model::InfiniteModel, status::Bool)
     model.ready_to_optimize = status
     return
end

"""

"""
function transform_model(backend::AbstractTransformationBackend)
    error("")
end

"""

"""
transform_model(model::InfiniteModel) = backend_model(model.backend)

"""

"""
function set_transformation_backend(
    model::InfiniteModel, 
    backend::AbstractTransformationBackend
    )
    model.backend = backend
    set_transformation_backend_ready(model, false)
    return
end

"""

"""
function JuMP.get_attribute(
    backend::AbstractTransformationBackend,
    attr
    )
    error("")
end

"""

"""
function JuMP.get_attribute(model::InfiniteModel, attr)
    return JuMP.get_attribute(model.backend, attr)
end

"""

"""
function JuMP.set_attribute(
    backend::AbstractTransformationBackend,
    attr,
    value
    )
    error("")
end

"""

"""
function JuMP.set_attribute(model::InfiniteModel, attr, value)
    return JuMP.set_attribute(model.backend, attr, value)
end

"""

"""
function build_transformation_model!(
    model::InfiniteModel, 
    backend::AbstractTransformationBackend;
    kwargs...
    )
    return error("")
end

"""

"""
function build_transformation_model!(model::InfiniteModel; kwargs...)
    if num_parameters(model, InfiniteParameter) == 0
        @warn("Finite models (i.e., `InfiniteModel`s with no infinite " * 
              "parameters) should be modeled directly via a `Model` in JuMP.jl.")
    end
    build_optimizer_model!(model, model.backend; kwargs...)
    set_transformation_backend_ready(model, true)
    return
end

"""
    JuMP.set_optimize_hook(
        model::InfiniteModel, 
        hook::Union{Function, Nothing}
        )::Nothing

Set the function `hook` as the optimize hook for `model` where `hook` should 
have be of the form `hook(model::InfiniteModel; hook_specfic_kwargs..., kwargs...)`. 
The `kwargs` are those passed to [`optimize!`](@ref). The `hook_specifc_kwargs` 
are passed as additional keywords by the user when they call [`optimize!`](@ref).

## Notes

* The optimize hook should generally modify the model, or some external state
in some way, and then call `optimize!(model; ignore_optimize_hook = true)` to
optimize the problem, bypassing the hook.
* Use `set_optimize_hook(model, nothing)` to unset an optimize hook.
"""
function JuMP.set_optimize_hook(
    model::InfiniteModel, 
    hook::Union{Function, Nothing}
    )
    model.optimize_hook = hook
    set_optimizer_model_ready(model, false)
    return
end

"""

"""
function JuMP.optimize!(backend::AbstractTransformationBackend)
    return error("")
end


"""
    JuMP.optimize!(model::InfiniteModel; [kwargs...])

Extend `JuMP.optimize!` to optimize infinite models using the internal
optimizer model. Calls [`build_optimizer_model!`](@ref) if the optimizer
model isn't up to date. The `kwargs` correspond to keyword arguments passed to
[`build_optimizer_model!`](@ref) if any are defined. The `kwargs` can also 
include arguments that are passed to an optimize hook if one was set with 
[`JuMP.set_optimize_hook`](@ref). 

**Example**
```julia-repl
julia> optimize!(model)

julia> has_values(model)
true
```
"""
function JuMP.optimize!(
    model::InfiniteModel; 
    ignore_optimize_hook = isnothing(model.optimize_hook),
    kwargs...)
    if !ignore_optimize_hook
        return model.optimize_hook(model; kwargs...)
    end
    if !optimizer_model_ready(model)
        build_optimizer_model!(model; kwargs...)
    end
    JuMP.optimize!(model.backend)
    return
end

"""

"""
function optimizer_model_variable(
    vref::GeneralVariableRef, 
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    error("")
end

"""

"""
function optimizer_model_variable(vref::GeneralVariableRef; kwargs...)
    model = JuMP.owner_model(vref)
    return optimizer_model_variable(vref, model.backend; kwargs...)
end

