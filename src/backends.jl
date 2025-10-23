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
    set_transformation_backend_ready(model::InfiniteModel, status::Bool)::Nothing

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
    transformation_model(backend::AbstractTransformationBackend)

Return the underlying model used by the `backend`. This serves as an 
extension point for new backend types. No extension is needed for 
[`JuMPBackend`](@ref)s.
"""
function transformation_model(backend::AbstractTransformationBackend)
    error("`transformation_model` not implemented for transformation backends " * 
          "of type `$(typeof(backend))`.")
end

"""
    transformation_model(model::InfiniteModel)

Return the underlying model used by the transformation backend.

**Example**
```julia-repl
julia> trans_model = transformation_model(model)
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```
"""
transformation_model(model::InfiniteModel) = transformation_model(model.backend)

"""
    transformation_data(backend::AbstractTransformationBackend)

Return the underlying data (typically mapping data) used by the `backend`. 
This serves as an extension point for new backend types. No extension 
is needed for [`JuMPBackend`](@ref)s.
"""
function transformation_data(backend::AbstractTransformationBackend)
    error("`transformation_data` not implemented for transformation backends " * 
          "of type `$(typeof(backend))`.")
end

"""
    transformation_data(model::InfiniteModel)

Return the underlying data (typically mapping data) used by the 
transformation backend. 

**Example**
```julia-repl
julia> mapping_data = transformation_data(model);
```
"""
function transformation_data(model::InfiniteModel)
    return transformation_data(model.backend)
end

"""
    set_transformation_backend(
        model::InfiniteModel, 
        backend::AbstractTransformationBackend
        )::Nothing

Specify a new transformation backend `backend` for the `model`. Note
that all data/settings/results associated with the previous backend 
will be removed.

**Example**
```julia-repl
julia> transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> set_transformation_backend(model, TranscriptionBackend(Ipopt.Optimizer))

julia> transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: Ipopt
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```
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
    transformation_backend(
        model::InfiniteModel
        )::AbstractTransformationBackend

Retrieve the transformation backend used by the `model`.

**Example**
```julia-repl
julia> transformation_backend(model)
A TranscriptionBackend that uses a
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```
"""
function transformation_backend(model::InfiniteModel)
    return model.backend
end

"""
    JuMP.get_attribute(backend::AbstractTransformationBackend, attr)

Retrieve some attribute `attr` from the `backend`. This is a general 
purpose method typically used to query optimizer related information. 
This serves as an extension point for new backend types. New backends
should include extensions for the following attributes as appropriate:
    - `MOI.Silent`
    - `MOI.TimeLimitSec`
    - `MOI.RawOptimizerAttribute`
    - `MOI.SolverName`
    - `MOI.TerminationStatus`
    - `MOI.RawStatusString`
    - `MOI.PrimalStatus`
    - `MOI.DualStatus`
    - `MOI.SolveTimeSec`
    - `MOI.ResultCount`
    - `MOI.SimplexIterations`
    - `MOI.BarrierIterations`
    - `MOI.NodeCount`
    - `MOI.ObjectiveBound`
    - `MOI.RelativeGap`
    - `MOI.ObjectiveValue`
    - `MOI.DualObjectiveValue`
    
No extension is needed for [`JuMPBackend`](@ref)s. 
"""
function JuMP.get_attribute(
    backend::AbstractTransformationBackend,
    attr
    )
    error("`JuMP.get_attribute` not implemented for transformation backends " * 
          "of type `$(typeof(backend))` with attribute `$attr`.")
end

"""
    JuMP.get_attribute(model::InfiniteModel, attr)

Retrieve an attribute `attr` from the transformation backend of 
`model`. Typically, this corresponds to `MOI.AbstractOptimizerAttribute`s.

**Example**
```julia-repl
julia> get_attribute(model, MOI.TimeLimitSec())
60.0
```
"""
function JuMP.get_attribute(model::InfiniteModel, attr)
    return JuMP.get_attribute(model.backend, attr)
end
function JuMP.get_attribute(model::InfiniteModel, attr::String)
    return JuMP.get_attribute(model.backend, MOI.RawOptimizerAttribute(attr))
end
function JuMP.get_optimizer_attribute(model::InfiniteModel, attr)
    return JuMP.get_attribute(model, attr)
end

"""
    JuMP.set_attribute(backend::AbstractTransformationBackend, attr, value)::Nothing

Specify some attribute `attr` to the `backend`. This is a general 
purpose method typically used to set optimizer related information. 
This serves as an extension point for new backend types. New backends
should include extensions for attributes of type:
    - `MOI.Silent`
    - `MOI.TimeLimitSec`
    - `MOI.RawOptimizerAttribute`
    
No extension is needed for [`JuMPBackend`](@ref)s. 
"""
function JuMP.set_attribute(
    backend::AbstractTransformationBackend,
    attr,
    value
    )
    error("`JuMP.set_attribute` not implemented for transformation backends " * 
          "of type `$(typeof(backend))` with attribute `$attr`.")
end

"""
    JuMP.set_attribute(model::InfiniteModel, attr, value)::Nothing

Specify an attribute `attr` to the transformation backend of 
`model`. Typically, this corresponds to `MOI.AbstractOptimizerAttribute`s.

**Example**
```julia-repl
julia> set_attribute(model, MOI.TimeLimitSec(), 42.0)
```
"""
function JuMP.set_attribute(model::InfiniteModel, attr, value)
    return JuMP.set_attribute(model.backend, attr, value)
end
function JuMP.set_attribute(model::InfiniteModel, attr::String, value)
    return JuMP.set_attribute(model.backend, MOI.RawOptimizerAttribute(attr), value)
end
function JuMP.set_optimizer_attribute(model::InfiniteModel, attr, value)
    return JuMP.set_attribute(model, attr, value)
end

"""
    JuMP.set_attributes(model::InfiniteModel, pairs::Pair...)::Nothing

Specify multiple optimizer transformation backend attributes as `Pair`s
of the form `attr => value` which are used for `set_attribute(model, attr, value)`.

**Example**
```julia-repl
julia> set_attributes(model, "tol" => 1e-4, "max_iter" => 100)
```
"""
function JuMP.set_attributes(model::InfiniteModel, pairs::Pair...)
    for (attr, value) in pairs
        JuMP.set_attribute(model.backend, attr, value)
    end
    return
end

"""
    Base.empty!(backend::AbstractTransformationBackend)

Empty `backend` of all its contents. For new backend types, this needs to 
be defined such that `empty!(model::InfiniteModel)` works. For 
[`JuMPBackend`](@ref)s this defaults to
```julia
empty!(transformation_model(backend))
empty!(transformation_data(backend))
```
"""
function Base.empty!(backend::AbstractTransformationBackend)
    error("`empty!` not implemented for transformation backends of type " * 
          "`$(typeof(backend))`.")
end

"""
    build_transformation_backend!(
        model::InfiniteModel, 
        backend::AbstractTransformationBackend;
        [kwargs...]
        )::Nothing

Given `model`, transform it into the representation used by `backend`. 
Once completed, `backend` should be ready to be solved. This serves as 
an extension point for new types of backends. If needed, keyword arguments
can be added. Typically, this should clear out the backend before reconstructing
it.
"""
function build_transformation_backend!(
    model::InfiniteModel, 
    backend::AbstractTransformationBackend;
    kwargs...
    )
    error("`build_transformation_backend!` not implemented for transformation backends " * 
          "of type `$(typeof(backend))`.")
end

"""
    build_transformation_backend!(model::InfiniteModel; [kwargs...])::Nothing

Build the model used by the underlying transformation backend stored in `model` such 
that it is ready to solve. Specifically, translate the InfiniteOpt formulation 
stored in `model` into (typically an appoximate) formulation that is compatible 
with the backend. This is called automatically by `optimize!`; however, it this 
method can be used to build the transformation backend without solving it.

**Example**
```julia-repl
julia> build_transformation_backend!(model)

julia> transformation_backend_ready(model)
true
```
"""
function build_transformation_backend!(model::InfiniteModel; kwargs...)
    if num_parameters(model, InfiniteParameter) == 0
        @warn("Finite models (i.e., `InfiniteModel`s with no infinite " * 
              "parameters) should be modeled directly via a `Model` in JuMP.jl.")
    end
    build_transformation_backend!(model, model.backend; kwargs...)
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
    set_transformation_backend_ready(model, false)
    return
end

"""
    JuMP.optimize!(backend::AbstractTransformationBackend)

Invoke the relevant routines to solve the underlying model used by 
`backend`. Note that [`build_transformation_backend!`](@ref) will be 
called before this method is. This needs to be extended for new 
backend types, but no extension is needed for [`JuMPBackend`](@ref)s.
Optionally, information can be returned if desired.
"""
function JuMP.optimize!(backend::AbstractTransformationBackend)
    error("`JuMP.optimize!` not implemented for transformation backends " * 
          "of type `$(typeof(backend))`.")
end


"""
    JuMP.optimize!(model::InfiniteModel; [kwargs...])

Extend `JuMP.optimize!` to optimize infinite models using the internal
transformation backend. Calls [`build_transformation_backend!`](@ref) if the optimizer
model isn't up-to-date. The `kwargs` correspond to keyword arguments passed to
[`build_transformation_backend!`](@ref) if any are defined. The `kwargs` can also 
include arguments that are passed to an optimize hook if one was set with 
[`JuMP.set_optimize_hook`](@ref). Typically, this returns `nothing`, but 
certain backends may return something.

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
    kwargs...
    )
    if !ignore_optimize_hook
        return model.optimize_hook(model; kwargs...)
    end
    if !transformation_backend_ready(model)
        build_transformation_backend!(model; kwargs...)
    end
    return JuMP.optimize!(model.backend)
end

################################################################################
#                            JUMP-BASED OPTIMIZER API
################################################################################
# 1 arg setters
for (func, Attr, val) in (
    (:set_silent, :Silent, true),
    (:unset_silent, :Silent, false),
    (:unset_time_limit_sec, :TimeLimitSec, nothing)
    )
    @eval begin
        @doc """
            JuMP.$($func)(model::InfiniteModel)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func))
        to accept `InfiniteModel`s. This relies on the underlying transformation 
        backend supporting `JuMP.set_attribute` with attribute `$(MOI.$Attr)()`.
        """
        function JuMP.$func(model::InfiniteModel)
            return JuMP.set_attribute(model.backend, MOI.$Attr(), $val)
        end
    end
end

# 1 arg getters
for (func, Attr) in ((:time_limit_sec, :TimeLimitSec), (:solver_name, :SolverName))
    @eval begin
        @doc """
            JuMP.$($func)(model::InfiniteModel)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func))
        to accept `InfiniteModel`s. This relies on the underlying transformation 
        backend supporting `JuMP.get_attribute` with attribute `$(MOI.$Attr)()`.
        """
        function JuMP.$func(model::InfiniteModel)
            return JuMP.get_attribute(model.backend, MOI.$Attr())
        end
    end
end

"""
    JuMP.set_time_limit_sec(model::InfiniteModel, value::Real)

Extend [`JuMP.set_time_limit_sec`](https://jump.dev/JuMP.jl/v1/api/JuMP/#set_time_limit_sec)
to accept `InfiniteModel`s. This relies on the underlying transformation 
backend supporting `JuMP.set_attribute` with attribute `MOI.TimeLimitSec()`.
"""
function JuMP.set_time_limit_sec(model::InfiniteModel, value::Real)
    return JuMP.set_attribute(model.backend, MOI.TimeLimitSec(), Float64(value))
end

# Single argument methods that don't rely on `[get/set]_attribute`
for func in (:bridge_constraints, :backend, :mode, :unsafe_backend, 
             :compute_conflict!, :copy_conflict)
    @eval begin
        @doc """
            JuMP.$($func)(backend::AbstractTransformationBackend)

        Implement [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func))
        for transformation backends. If applicable, this should be extended for 
        new backend types. No extension is needed for [`JuMPBackend`](@ref)s.
        """
        function JuMP.$func(backend::AbstractTransformationBackend)
            error("`JuMP.$($func)` not defined for backends of type " *
                "`$(typeof(backend))`.")
        end

        # Define for JuMPBackend
        function JuMP.$func(backend::JuMPBackend)
            return JuMP.$func(backend.model)
        end

        @doc """
            JuMP.$($func)(model::InfiniteModel)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func))
        to accept `InfiniteModel`s. This relies on the underlying transformation 
        backend supporting `JuMP.$($func)`.
        """
        function JuMP.$func(model::InfiniteModel)
            return JuMP.$func(model.backend)
        end 
    end
end


"""
    JuMP.add_bridge(backend::AbstractTransformationBackend, value)

Implement [`JuMP.add_bridge`](https://jump.dev/JuMP.jl/v1/api/JuMP/#add_bridge)
for transformation backends. If applicable, this should be extended for 
new backend types. No extension is needed for [`JuMPBackend`](@ref)s.
"""
function JuMP.add_bridge(backend::AbstractTransformationBackend, value)
    error("`JuMP.add_bridge` not defined for backends of type " *
        "`$(typeof(backend))`.")
end

# Define for JuMPBackend
function JuMP.add_bridge(backend::JuMPBackend, value)
    return JuMP.add_bridge(backend.model, value)
end

"""
    JuMP.add_bridge(model::InfiniteModel, value)

Extend [`JuMP.add_bridge`](https://jump.dev/JuMP.jl/v1/api/JuMP/#add_bridge)
to accept `InfiniteModel`s. This relies on the underlying transformation 
backend supporting `JuMP.add_bridge`.
"""
function JuMP.add_bridge(model::InfiniteModel, value)
    return JuMP.add_bridge(model.backend, value)
end 

"""
    JuMP.print_active_bridges(
        io::IO, 
        backend::AbstractTransformationBackend, 
        args...
        )

Implment `JuMP.print_active_bridges` for transformation backends. If applicable, this 
should be extended for new backend types. No extension is needed for 
[`JuMPBackend`](@ref)s. Here, `args` can be one of the following:
- empty (print all the bridges)
- the objective type (print the objective bridges)
- a function type and set type from a constraint
- a constraint set type
"""
function JuMP.print_active_bridges(
    io::IO, 
    backend::AbstractTransformationBackend, 
    args...
    )
    error("`JuMP.print_active_bridges` not defined for backends of type " *
            "`$(typeof(backend))`.")
end

# Define for JuMPBackend
function JuMP.print_active_bridges(io::IO, backend::JuMPBackend, args...)
    return JuMP.print_active_bridges(io, backend.model, args...)
end

"""
    JuMP.print_active_bridges([io::IO = stdout,] model::InfiniteModel)

    JuMP.print_active_bridges([io::IO = stdout,] model::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar})

    JuMP.print_active_bridges([io::IO = stdout,] model::InfiniteModel, ::Type{<:JuMP.AbstractJuMPScalar}, ::Type{<:MOI.AbstractSet})

    JuMP.print_active_bridges([io::IO = stdout,] model::InfiniteModel, ::Type{<:MOI.AbstractSet})

Extend [`JuMP.print_active_bridges`](https://jump.dev/JuMP.jl/v1/api/JuMP/#print_active_bridges)
to accept `InfiniteModel`s. This relies on the underlying transformation 
backend supporting `JuMP.print_active_bridges`.
"""
function JuMP.print_active_bridges(io::IO, model::InfiniteModel, args...)
    return JuMP.print_active_bridges(io, model.backend, args...)
end 
function JuMP.print_active_bridges(model::InfiniteModel, args...)
    return JuMP.print_active_bridges(Base.stdout, model.backend, args...)
end

"""
    JuMP.print_bridge_graph(
        io::IO, 
        backend::AbstractTransformationBackend
        )

Implment `JuMP.print_bridge_graph` for transformation backends. If applicable, this 
should be extended for new backend types. No extension is needed for 
[`JuMPBackend`](@ref)s.
"""
function JuMP.print_bridge_graph(
    io::IO, 
    backend::AbstractTransformationBackend,
    )
    error("`JuMP.print_bridge_graph` not defined for backends of type " *
            "`$(typeof(backend))`.")
end

# Define for JuMPBackend
function JuMP.print_bridge_graph(io::IO, backend::JuMPBackend)
    return JuMP.print_bridge_graph(io, backend.model)
end

"""
    JuMP.print_bridge_graph([io::IO = stdout,] model::InfiniteModel)

Extend [`JuMP.print_bridge_graph`](https://jump.dev/JuMP.jl/v1/api/JuMP/#print_bridge_graph)
to accept `InfiniteModel`s. This relies on the underlying transformation 
backend supporting `JuMP.print_bridge_graph`.
"""
function JuMP.print_bridge_graph(io::IO, model::InfiniteModel)
    return JuMP.print_bridge_graph(io, model.backend)
end 
function JuMP.print_bridge_graph(model::InfiniteModel)
    return JuMP.print_bridge_graph(Base.stdout, model.backend)
end

"""
    JuMP.set_optimizer(
        backend::AbstractTransformationBackend,
        optimizer_constructor;
        [kwargs...]
        )::Nothing

Specify the optimizer `optimizer_constructor` that should be used by `backend`. 
This is intended as an extension point for new transformation backend types. 
Keyword arguments can be added as needed. No extension is necessary for 
[`JuMPBackend`](@ref)s.
"""
function JuMP.set_optimizer(
    backend::AbstractTransformationBackend,
    optimizer_constructor;
    kwargs...
    )
    error("`JuMP.set_optimizer` not defined for transformation backends " *
          "of type `$(typeof(backend))` with optimizer input " * 
          "`$(optimizer_constructor)`.")
end

# JuMPBackend
function JuMP.set_optimizer(
    backend::JuMPBackend,
    optimizer_constructor;
    add_bridges::Bool = true
    )
    JuMP.set_optimizer(backend.model, optimizer_constructor,
                       add_bridges = add_bridges)
    return
end

"""
    JuMP.set_optimizer(
        model::InfiniteModel,
        [optimizer_constructor;
        add_bridges::Bool = true, 
        kwargs...]
        )

Extend `JuMP.set_optimizer` to set optimizer used by the underlying 
transformation backend associated with `model`. If a backend uses 
`JuMP` then `add_bridges` can be used as a keyword argument.

**Example**
```julia-repl
julia> set_optimizer(model, HiGHS.Optimizer)

julia> transformation_model(model)
A JuMP Model
├ solver: HiGHS
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none
```
"""
function JuMP.set_optimizer(
    model::InfiniteModel,
    optimizer_constructor;
    kwargs...
    )
    return JuMP.set_optimizer(model.backend, optimizer_constructor; kwargs...)
end

# JuMPBackend dispatches
transformation_model(backend::JuMPBackend) = backend.model
transformation_data(backend::JuMPBackend) = backend.data
function JuMP.get_attribute(backend::JuMPBackend, attr)
    return JuMP.get_attribute(backend.model, attr)
end
function JuMP.set_attribute(backend::JuMPBackend, attr, val)
    return JuMP.set_attribute(backend.model, attr, val)
end
function Base.empty!(backend::JuMPBackend)
    empty!(transformation_model(backend))
    empty!(transformation_data(backend))
    return backend
end
function JuMP.optimize!(backend::JuMPBackend)
    return JuMP.optimize!(backend.model)
end

################################################################################
#                             VARIABLE MAPPING API
################################################################################
"""
    transformation_variable(
        vref::GeneralVariableRef, 
        backend::AbstractTransformationBackend; 
        [kwargs...]
        )

Return the variable(s) that map to `vref` used by `backend`. This serves as an 
extension point for new backend types. If needed, keywords arguments can be 
added.
"""
function transformation_variable(
    vref::GeneralVariableRef, 
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    error("`transformation_variable` not defined for backends of type " *
          "`$(typeof(backend))`.")
end

"""
    transformation_variable(vref::GeneralVariableRef; [kwargs...])

Returns the variable(s) used by the transformation backend to represent `vref`. 
Certain backends may also allow the use of keyward arguments. 

The default backend `TranscriptionOpt` uses the keyword arguments:
- `label::Type{<:AbstractSupportLabel} = PublicLabel`
By default only variables corresponding to public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the transcripted
variables of infinite variables are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on.

**Example**
```julia-repl
julia> transformation_variable(x) # infinite variable
2-element Array{VariableRef,1}:
 x(0.0)
 x(1.0)

julia> transformation_variable(z) # finite variable
z
```
"""
function transformation_variable(vref::GeneralVariableRef; kwargs...)
    model = JuMP.owner_model(vref)
    return transformation_variable(vref, model.backend; kwargs...)
end

"""
    variable_supports(
        vref::DecisionVariableRef,
        backend::AbstractTransformationBackend;
        [kwargs...]
        )

Return the supports associated with the mappings of `vref` in `backend`.
This dispatches off of `backend` which permits transformation backend extensions. This
should throw an error if `vref` is not associated with the variable mappings
stored in `backend`. Keyword arguments can be added as needed. Note that
no extension is necessary for point or finite variables. 
"""
function variable_supports(vref, backend::AbstractTransformationBackend; kwargs...)
    error("`variable_supports` not implemented for transformation backend of type " *
          "`$(typeof(backend))` and/or variable type $(typeof(vref)).")
end

# FiniteRef
function variable_supports(
    vref::FiniteRef,
    backend::AbstractTransformationBackend;
    kwargs...
    )
    return ()
end

"""
    supports(
        vref::DecisionVariableRef; 
        [label::Type{<:AbstractSupportLabel} = PublicLabel, 
        kwargs...]
        )::Vector{<:Tuple}

Return the supports associated with `vref` in the transformation
model. Errors if [`InfiniteOpt.variable_supports`](@ref) has not been extended for the
transformation backend type or if `vref` is not reformulated in the transformation backend.

The keyword argument `label` is what `TranscriptionOpt` employs
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `variable_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the supports
of infinite variables are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on.

**Example**
```julia-repl
julia> supports(vref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(
    vref::Union{DecisionVariableRef, MeasureRef, ParameterFunctionRef}; 
    kwargs...
    )
    backend = JuMP.owner_model(vref).backend
    return variable_supports(vref, backend; kwargs...)
end

################################################################################
#                             EXPRESSION MAPPING API
################################################################################
"""
    transformation_expression(expr, backend::AbstractTransformationBackend; [kwargs...])

Return the reformulation expression(s) stored in the transformation backend that correspond
to `expr`. This needs to be defined for extensions that implement a new 
[`AbstractTransformationBackend`](@ref). Keyword arguments can be added as needed.
Note that if `expr` is a `GeneralVariableRef` this just dispatches to
`transformation_variable`.
"""
function transformation_expression(
    expr, 
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    error("`transformation_expression` not defined for transformation backends " *
          "of type `$(typeof(backend))` and expression type `$(typeof(expr))`.")
end

# Define for variable reference expressions
function transformation_expression(
    expr::GeneralVariableRef, 
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    return transformation_variable(expr, backend; kwargs...)
end

"""
    transformation_expression(
        expr::JuMP.AbstractJuMPScalar; 
        [label::Type{<:AbstractSupportLabel} = PublicLabel,
        kwargs...]
        )

Return the reformulation expression(s) stored in the transformation backend that correspond
to `expr`. Also errors if no such expression can be found in
the transformation backend (meaning one or more of the underlying variables have not
been transformed).

The keyword argument `label` is what `TranscriptionOpt` employs
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of [`transformation_expression`](@ref). Errors if such an
extension has not been written. 

By default only the expressions associated with public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the transformed
expressions are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on. The corresponding supports are obtained via 
`supports` using the same keyword arguments.

**Example**
```julia-repl
julia> transformation_expression(my_expr) # finite expression
x(0.0) - y
```
"""
function transformation_expression(expr::JuMP.AbstractJuMPScalar; kwargs...)
    model = JuMP.owner_model(expr)
    if isnothing(model)
        return zero(JuMP.AffExpr) + JuMP.constant(expr)
    else
        return transformation_expression(expr, model.backend; kwargs...)
    end
end

"""
    expression_supports(
        expr,
        backend::AbstractTransformationBackend;
        [kwargs...]
        )

Return the supports associated with the mappings of `expr` in `backend`.
This should throw an error if `expr` is not associated with the variable mappings
stored in `backend`. Keyword arguments can be added as needed. Note that
if `expr` is a `GeneralVariableRef` this just dispatches to `variable_supports`.
"""
function expression_supports(expr, backend::AbstractTransformationBackend; kwargs...)
  error("`expression_supports` not implemented for transformation backend of type " *
        "`$(typeof(backend))` and/or expressions of type `$(typeof(expr))`.")
end

# Variable reference expressions
function expression_supports(
    vref::GeneralVariableRef, 
    backend::AbstractTransformationBackend;
    kwargs...
    )
    return variable_supports(dispatch_variable_ref(vref), backend; kwargs...)
end

"""
    supports(
        expr::JuMP.AbstractJuMPScalar; 
        [label::Type{<:AbstractSupportLabel} = PublicLabel,
        kwargs...]
        )

Return the support associated with `expr`. Errors if `expr` is
not associated with the constraint mappings stored in the transformation backend.

The keyword arguments `label` is what `TranscriptionOpt` employs
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `expression_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the supports
of an infinite expression are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on.

**Example**
```julia-repl
julia> supports(cref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(expr::JuMP.AbstractJuMPScalar; kwargs...)
    model = JuMP.owner_model(expr)
    if isnothing(model)
        return ()
    else
        return expression_supports(expr, model.backend; kwargs...)
    end
end

################################################################################
#                             CONSTRAINT MAPPING API
################################################################################
"""
    transformation_constraint(
        cref::InfOptConstraintRef,
        backend::AbstractTransformationBackend; 
        [kwargs...]
        )

Return the reformulation constraint(s) stored in the transformation backend 
that correspond to `cref`. This needs to be defined for extensions that 
implement a custom transformation backend type. Keyword arguments can be 
added as needed.
"""
function transformation_constraint(
    cref::InfOptConstraintRef,
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    error("`transformation_constraint` not implemented for " * 
          "transformation backends of type `$(typeof(backend))`.")
end

"""
    transformation_constraint(
        cref::InfOptConstraintRef; 
        [label::Type{<:AbstractSupportLabel} = PublicLabel, 
        kwargs...]
        )

Return the reformulation constraint(s) stored in the transformation backend that 
correspond to `cref`. Errors if no such constraint can be found in
the transformation backend.

The keyword argument `label` is what `TranscriptionOpt` employs
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of [`transformation_constraint`](@ref). Errors if such an
extension has not been written. 

By default only the constraints associated with public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the transformed
cosntraints are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on. The corresponding supports are obtained via 
`supports` using the same keyword arguments.

**Example**
```julia-repl
julia> transformation_constraint(c1) # finite constraint
c1 : x(0.0) - y <= 3.0
```
"""
function transformation_constraint(
    cref::InfOptConstraintRef; 
    kwargs...
    )
    backend = JuMP.owner_model(cref).backend
    return transformation_constraint(cref, backend; kwargs...)
end

"""
    constraint_supports(
        cref::InfOptConstraintRef
        backend::AbstractTransformationBackend; 
        [kwargs...]
        )

Return the supports associated with the mappings of `cref` in `backend`.
This should throw an error if `cref` is not associated with the variable mappings
stored in `backend`. Keyword arguments can be added as needed.
"""
function constraint_supports(
    cref::InfOptConstraintRef,
    backend::AbstractTransformationBackend; 
    kwargs...
    )
    error("`constraint_supports` not implemented for transformation backends " *
          "of type `$(typeof(backend))`.")
end

"""
    supports(cref::InfOptConstraintRef; 
             [label::Type{<:AbstractSupportLabel} = PublicLabel,
             kwargs...])

Return the support associated with `cref`. Errors if `cref` is
not associated with the constraint mappings stored in the transformation backend.

The keyword argument `label` is what `TranscriptionOpt` employs
and `kwargs` denote extra ones that user extensions may employ in accordance with
their implementation of `constraint_supports`. Errors if such an
extension has not been written. 

By default only the public supports are returned, the 
full set can be accessed via `label = All`. Where possible, all the supports
of the constraint are returned as an n-dimensional array 
where each dimension is determined by the each independent group of
infinite parameters it depends on.

**Example**
```julia-repl
julia> supports(cref)
2-element Array{Tuple{Float64},1}:
 (0.0,)
 (1.0,)
```
"""
function supports(cref::InfOptConstraintRef; kwargs...)
    model = JuMP.owner_model(cref)
    return constraint_supports(cref, model.backend; kwargs...)
end

################################################################################
#                                 UPDATING API
################################################################################
"""
    update_parameter_value(
        backend::AbstractTransformationBackend,
        ref::Union{FiniteParameterRef, ParameterFunctionRef},
        value
    )::Bool

If `backend`` is built, then this method updates what `ref` corresponds to in the
`backend` to `value`, then it returns a `Bool` on whether the update was successful.
This is intended as an extension point for new `AbstractTransformationBackend`s to
more efficiently handle parameter updates for resolves. This defaults to `false`,
meaning that no update occurs (forcing the backend to be rebuilt). Users should
use [`JuMP.set_parameter_value`](@ref) rather than call this method directly.
"""
update_parameter_value(backend::AbstractTransformationBackend, ref, value) = false
