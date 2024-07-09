################################################################################
#                                  MODEL QUERIES
################################################################################
# Ensure results are reliable
function _check_result_is_current(model::InfiniteModel, func)
    if JuMP.termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        error("Unable to query `$func` since up-to-date solution results " *
              "are not available. This is likely because the transformation " *
              "backend has not been built and optimized yet, or the model " *
              "has been modified since the backend was optimized.")
    end
    return
end

# Make `get_attribute` check for reliability
function JuMP.get_attribute(model::InfiniteModel, attr::MOI.TerminationStatus)
    if !transformation_backend_ready(model)
        return MOI.OPTIMIZE_NOT_CALLED
    end
    return JuMP.get_attribute(model.backend, attr)
end
function JuMP.get_attribute(model::InfiniteModel, attr::MOI.ResultCount)
    if JuMP.termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        return 0
    end
    return JuMP.get_attribute(model.backend, attr)
end
function JuMP.get_attribute(model::InfiniteModel, attr::MOI.RawStatusString)
    if JuMP.termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        return "optimize not called"
    end
    return JuMP.get_attribute(model.backend, attr)
end
function JuMP.get_attribute(
    model::InfiniteModel,
    attr::Union{MOI.PrimalStatus, MOI.DualStatus}
    )
    if JuMP.termination_status(model) == MOI.OPTIMIZE_NOT_CALLED
        return MOI.NO_SOLUTION
    end
    return JuMP.get_attribute(model.backend, attr)
end
function JuMP.get_attribute(
    model::InfiniteModel,
    attr::MOI.AbstractModelAttribute
    )
    if MOI.is_set_by_optimize(attr)
        _check_result_is_current(model::InfiniteModel, attr)
    end
    return JuMP.get_attribute(model.backend, attr)
end

# Simple model queries
for (func, Attr) in (
    (:termination_status, :TerminationStatus),
    (:raw_status, :RawStatusString),
    (:solve_time, :SolveTimeSec),
    (:simplex_iterations, :SimplexIterations),
    (:barrier_iterations, :BarrierIterations),
    (:node_count, :NodeCount),
    (:objective_bound, :ObjectiveBound),
    (:relative_gap, :RelativeGap),
    (:result_count, :ResultCount)
    )
    @eval begin
        @doc """
            JuMP.$($func)(model::InfiniteModel)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func)) 
        for `InfiniteModel`s in accordance with that reported by its 
        transformation backend. Errors if such a query is not supported or if 
        the transformation backend hasn't been solved.
        """
        function JuMP.$func(model::InfiniteModel)
            JuMP.get_attribute(model, MOI.$Attr())
        end
    end
end

# Simple result dependent model queries
for (func, Attr) in (
    (:primal_status, :PrimalStatus),
    (:dual_status, :DualStatus),
    (:objective_value, :ObjectiveValue),
    (:dual_objective_value, :DualObjectiveValue)
    )
    @eval begin 
        @doc """
            JuMP.$($func)(model::InfiniteModel; result::Int = 1)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func)) 
        for `InfiniteModel`s in accordance with that reported by its 
        transformation backend. Errors if such a query is not supported or if the 
        transformation backend hasn't be solved. Accepts keyword `result`
        to access the solution index of interest (if the solver/backend supports 
        multiple solutions).
        """
        function JuMP.$func(model::InfiniteModel; result::Int = 1)
            return JuMP.get_attribute(model, MOI.$Attr(result))
        end
    end
end

# Simple result dependent model queries
for (func, Attr) in ((:has_values, :PrimalStatus), (:has_duals, :DualStatus))
    @eval begin 
        @doc """
            JuMP.$($func)(model::InfiniteModel; result::Int = 1)

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func)) 
        for `InfiniteModel`s in accordance with that reported by its 
        transformation backend. Errors if such a query is not supported or if the 
        transformation backend hasn't be solved. Accepts keyword `result`
        to access the solution index of interest (if the solver/backend supports 
        multiple solutions).
        """
        function JuMP.$func(model::InfiniteModel; result::Int = 1)
            return JuMP.get_attribute(model, MOI.$Attr(result)) != MOI.NO_SOLUTION
        end
    end
end

"""
    JuMP.is_solved_and_feasible(
        model::InfiniteModel;
        [dual::Bool = false,
        allow_local::Bool = true,
        allow_almost::Bool = false,
        result::Int = 1]
        )::Bool

Extend [`JuMP.is_solved_and_feasible`](https://jump.dev/JuMP.jl/v1/api/JuMP/#is_solved_and_feasible))
for `model`. See the JuMP docs details. 
For new transformation backend types, this relies on `JuMP.termination_status`,
`JuMP.primal_status`, and `JuMP.dual_status`.
"""
function JuMP.is_solved_and_feasible(
    model::InfiniteModel;
    dual::Bool = false,
    allow_local::Bool = true,
    allow_almost::Bool = false,
    result::Int = 1
    )
    status = JuMP.termination_status(model)
    ret =
        (status == MOI.OPTIMAL) ||
        (allow_local && (status == MOI.LOCALLY_SOLVED)) ||
        (allow_almost && (status == MOI.ALMOST_OPTIMAL)) ||
        (allow_almost && allow_local && (status == MOI.ALMOST_LOCALLY_SOLVED))
    if ret
        primal = JuMP.primal_status(model; result)
        ret &=
            (primal == MOI.FEASIBLE_POINT) ||
            (allow_almost && (primal == MOI.NEARLY_FEASIBLE_POINT))
    end
    if ret && dual
        dual_stat = JuMP.dual_status(model; result)
        ret &=
            (dual_stat == MOI.FEASIBLE_POINT) ||
            (allow_almost && (dual_stat == MOI.NEARLY_FEASIBLE_POINT))
    end
    return ret
end

################################################################################
#                                 VALUE QUERIES
################################################################################
"""
    map_value([ref/expr], backend::AbstractTransformationBackend; [kwargs...])

Map the value(s) of `ref` to its counterpart in the `backend`.
Here `ref` need refer to methods for both variable references and constraint
references. No extension is needed for [`JuMPBackend`](@ref)s that support
`transformation_variable`, `transformation_expression`, and 
`transformation_constraint`. In this case, `transformation_variable`, 
`transformation_expression`, and `transformation_constraint` are 
used to make these mappings by default where `kwargs` are passed on these functions. 
For mapping the values of infinite parameters, refer to 
[`map_infinite_parameter_value`](@ref).
"""
function map_value(ref, backend::AbstractTransformationBackend; kwargs...)
    error("Value queries are not supported for `$(typeof(ref))`s with a " *
          "transformation backend of type `$(typeof(backend))`. If you are " * 
          "writing an extension be sure to extend `map_value`.")
end

# Dispatch to deal with what is returned by parameter functions
_get_jump_value(v, result) = JuMP.value(v, result = result)
_get_jump_value(v::Real, result) = v

# Default method that depends on transformation_variable --> making extensions easier
function map_value(
    vref::GeneralVariableRef,
    backend::JuMPBackend;
    result::Int = 1,
    kwargs...
    )
    opt_vref = transformation_variable(vref, backend; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> _get_jump_value(v, result), opt_vref)
    else
        return _get_jump_value(opt_vref, result)
    end
end

# Default method that depends on transformation_expression --> making extensions easier
function map_value(
    expr::JuMP.AbstractJuMPScalar,
    backend::JuMPBackend;
    result::Int = 1,
    kwargs...
    )
    opt_expr = transformation_expression(expr, backend; kwargs...)
    if opt_expr isa AbstractArray
        return map(v -> _get_jump_value(v, result), opt_expr)
    else
        return _get_jump_value(opt_expr, result)
    end
end

# Default method that depends on transformation_constraint --> making extensions easier
function map_value(
    cref::InfOptConstraintRef,
    backend::JuMPBackend;
    result::Int = 1,
    kwargs...
    )
    opt_cref = transformation_constraint(cref, backend; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> _get_jump_value(c, result), opt_cref)
    else
        return _get_jump_value(opt_cref, result)
    end
end

"""
    map_infinite_parameter_value(
        pref::GeneralVariableRef, 
        backend::AbstractTransformationBackend;
        [kwargs...]
        )

Return the mapped value of the infinite parameter `pref` according to the 
`backend`. This serves as an optional extension point for new type of 
backends that do not rely on using supports. Otherwise, this defaults 
to:
```julia
map_infinite_parameter_value(pref; [label = PublicLabel]) = supports(pref, label = label)
```
"""
function map_infinite_parameter_value(
    pref::GeneralVariableRef, 
    backend::AbstractTransformationBackend; 
    label = PublicLabel
    )
    return supports(pref, label = label)
end

## Define dispatch methods to collect value of parameters 
# InfiniteParameter 
function _get_value(pref, ::Type{<:InfiniteParameterIndex}; kwargs...)
    backend = JuMP.owner_model(pref).backend
    return map_infinite_parameter_value(pref, backend; kwargs...)
end

# FiniteParameter 
function _get_value(pref, ::Type{FiniteParameterIndex}; kwargs...)
    return parameter_value(pref)
end

# Others 
function _get_value(vref, index_type; kwargs...)
    return map_value(vref, JuMP.owner_model(vref).backend; kwargs...)
end

"""
    JuMP.value(vref::GeneralVariableRef; [kwargs...])

Extend `JuMP.value` to return the value(s) of `vref` in accordance with its 
reformulation variable(s) stored in the transformation backend. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
whether a result exists before checking the values. 
    
Thw keyword arguments `kwargs` depend on the transformation backend that is 
being used. The default backend `TranscriptionOpt` uses the keyword 
arguments:
- `result::Int = 1`: indexes the solution result to be queried
- `label::Type{<:AbstractSupportLabel} = PublicLabel`: the label of supports to be returned
- `ndarray::Bool = false`: indicates whether the output should be formatted as an array
By default only the values associated with public supports (i.e., `PublicLabel`s) 
are returned, the full set can be accessed via `label = All`. Moreover, the values 
of infinite variables are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the variable has multiple 
infinite parameter dependencies.

To provide context for the values, it may be helpful to also query the variable's 
`parameter_refs` and `supports` which will have a one-to-one correspondence with 
the value(s). It may also be helpful to query via [`transformation_variable`](@ref) 
to retrieve the variables(s) that these values are based on. These functions should 
all be called with the same keyword arguments for consistency.

For extensions, this only works if 
[`transformation_variable`](@ref) has been extended correctly and/or 
[`map_value`](@ref) has been extended for variables.

**Example**
```julia-repl
julia> value(z)
42.0
```
"""
function JuMP.value(vref::GeneralVariableRef; kwargs...)
    _check_result_is_current(JuMP.owner_model(vref), JuMP.value)
    return _get_value(vref, _index_type(vref); kwargs...)
end

"""
    JuMP.value(expr::JuMP.AbstractJuMPScalar; [kwargs...])

Extend `JuMP.value` to return the value(s) of `vref` in accordance with its 
reformulation expression(s) stored in the transformation backend. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
whether a result exists before checking the values. 
    
Thw keyword arguments `kwargs` depend on the transformation backend that is 
being used. The default backend `TranscriptionOpt` uses the keyword 
arguments:
- `result::Int = 1`: indexes the solution result to be queried
- `label::Type{<:AbstractSupportLabel} = PublicLabel`: the label of supports to be returned
- `ndarray::Bool = false`: indicates whether the output should be formatted as an array
By default only the values associated with public supports (i.e., `PublicLabel`s) 
are returned, the full set can be accessed via `label = All`. Moreover, the values 
of infinite expressions are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the expression has multiple 
infinite parameter dependencies.

To provide context for the values, it may be helpful to also query the expression's 
`parameter_refs` and `supports` which will have a one-to-one correspondence with 
the value(s). It may also be helpful to query via [`transformation_expression`](@ref) 
to retrieve the expression(s) that these values are based on. These functions should 
all be called with the same keyword arguments for consistency.

For extensions, this only works if 
[`transformation_expression`](@ref) has been extended correctly and/or 
[`map_value`](@ref) has been extended for expressions.

**Example**
```julia-repl
julia> value(my_finite_expr)
23.34

julia> value(my_infinite_expr)
4-element Array{Float64,1}:
 -0.0
 20.9
 20.9
 20.9
```
"""
function JuMP.value(
    expr::Union{
        JuMP.GenericAffExpr{Float64, GeneralVariableRef}, 
        JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
        JuMP.GenericNonlinearExpr{GeneralVariableRef}
        };
    kwargs...
    )
    # get the model
    model = JuMP.owner_model(expr)
    # if no model then the expression only contains a constant
    if isnothing(model)
        expr isa JuMP.GenericNonlinearExpr && return JuMP.value(identity, expr)
        return JuMP.constant(expr)
    # otherwise let's call map_value
    else
        _check_result_is_current(model, JuMP.value)
        return map_value(expr, model.backend; kwargs...)
    end
end

"""
    JuMP.value(cref::InfOptConstraintRef; [kwargs...])

Extend `JuMP.value` to return the value(s) of `cref` in accordance with its 
reformulation constraint(s) stored in the transformation backend. Use
[`JuMP.has_values`](@ref JuMP.has_values(::InfiniteModel)) to check
whether a result exists before checking the values. 
    
Thw keyword arguments `kwargs` depend on the transformation backend that is 
being used. The default backend `TranscriptionOpt` uses the keyword 
arguments:
- `result::Int = 1`: indexes the solution result to be queried
- `label::Type{<:AbstractSupportLabel} = PublicLabel`: the label of supports to be returned
- `ndarray::Bool = false`: indicates whether the output should be formatted as an array
By default only the values associated with public supports (i.e., `PublicLabel`s) 
are returned, the full set can be accessed via `label = All`. Moreover, the values 
of infinite constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

To provide context for the values, it may be helpful to also query the constraint's 
`parameter_refs` and `supports` which will have a one-to-one correspondence with 
the value(s). It may also be helpful to query via [`transformation_constraint`](@ref) 
to retrieve the constraint(s) that these values are based on. These functions should 
all be called with the same keyword arguments for consistency.

For extensions, this only works if 
[`transformation_constraint`](@ref) has been extended correctly and/or 
[`map_value`](@ref) has been extended for constraints.

**Example**
```julia-repl
julia> value(c1)
4-element Array{Float64,1}:
 -0.0
 20.9
 20.9
 20.9
```
"""
function JuMP.value(cref::InfOptConstraintRef; kwargs...)
    model = JuMP.owner_model(cref)
    _check_result_is_current(model, JuMP.value)
    return map_value(cref, model.backend; kwargs...)
end

################################################################################
#                            BOILERPLATE REF QUERIES
################################################################################
for (Ref, func, mapper) in (
    (:GeneralVariableRef, :reduced_cost, :transformation_variable), 
    (:GeneralVariableRef, :optimizer_index, :transformation_variable),
    (:InfOptConstraintRef, :optimizer_index, :transformation_constraint),
    (:InfOptConstraintRef, :shadow_price, :transformation_constraint)
    )
    @eval begin 
        @doc """
            map_$($func)(
                ref::$($Ref),
                backend::AbstractTransformationBackend;
                [kwargs...]
                )

        Map `JuMP.$($func)` of `ref` to its counterpart in the `backend`.
        No extension is needed for [`JuMPBackend`](@ref)s that support
        `$($mapper)`, in which case, `$($mapper)` is used to make these 
        mappings using `kwargs`.
        """
        function $(Symbol(string("map_", func)))(
            ref::$Ref,
            backend::AbstractTransformationBackend;
            kwargs...
            )
            error("`$($func)` queries are not supported for a " *
                  "transformation backend of type `$(typeof(backend))`. If you " * 
                  "are writing an extension be sure to extend `map_$($func)`.")
        end

        # JuMPBackend
        function $(Symbol(string("map_", func)))(
            ref::$Ref,
            backend::JuMPBackend;
            kwargs...
            )
            opt_ref = $mapper(ref, backend; kwargs...)
            if opt_ref isa AbstractArray
                return map(r -> JuMP.$func(r), opt_ref)
            else
                return JuMP.$func(opt_ref)
            end
        end

        @doc """
            JuMP.$($func)(ref::$($Ref); [kwargs...])

        Extend [`JuMP.$($func)`](https://jump.dev/JuMP.jl/v1/api/JuMP/#$($func))
        for `ref`s in InfiniteModel. The exact format of output will depend 
        on the transformation backend that is being used.

        Thw keyword arguments `kwargs` depend on the transformation backend that is 
        being used. The default backend `TranscriptionOpt` uses the keyword 
        arguments:
        - `label::Type{<:AbstractSupportLabel} = PublicLabel`: the label of supports to be returned
        - `ndarray::Bool = false`: indicates whether the output should be formatted as an array
        By default only the values associated with public supports (i.e., `PublicLabel`s) 
        are returned, the full set can be accessed via `label = All`. Moreover, the values 
        of infinite variables/constraints are returned as a list. However, a n-dimensional array 
        can be obtained via `ndarray = true` which is handy when the constraint has multiple 
        infinite parameter dependencies.

        To provide context for the values, it may be helpful to also query the
        `parameter_refs` and `supports` which will have a one-to-one correspondence with 
        the output(s) of this function. These functions should 
        all be called with the same keyword arguments for consistency.
        """
        function JuMP.$func(ref::$Ref; kwargs...)
            model = JuMP.owner_model(ref)
            backend = model.backend
            _check_result_is_current(model, JuMP.$func)
            return $(Symbol(string("map_", func)))(ref, backend; kwargs...)
        end
    end
end

################################################################################
#                                  DUAL QUERIES
################################################################################
"""
    map_dual(
        cref::InfOptConstraintRef,
        backend::AbstractTransformationBackend;
        [kwargs...]
        )

Map the dual(s) of `cref` to its counterpart in the `backend`.
No extension is needed for [`JuMPBackend`](@ref)s that support
`transformation_constraint`. In this case, `transformation_constraint` 
are used to make these mappings by default where `kwargs` are passed on these 
functions.
"""
function map_dual(
    cref::InfOptConstraintRef,
    backend::AbstractTransformationBackend;
    kwargs...
    )
    error("Dual queries are not supported for a " *
          "transformation backend of type `$(typeof(backend))`. If you are " * 
          "writing an extension be sure to extend `map_dual`.")
end

# JuMPBackend default
function map_dual(
    cref::InfOptConstraintRef,
    backend::JuMPBackend;
    result::Int = 1,
    kwargs...
    )
    opt_cref = transformation_constraint(cref, backend; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.dual(c, result = result), opt_cref)
    else
        return JuMP.dual(opt_cref, result = result)
    end
end

"""
    JuMP.dual(cref::InfOptConstraintRef; [kwargs...])

Extend `JuMP.dual` to return the dual(s) of `cref` in accordance with its 
reformulation constraint(s) stored in the transformation backend. Use
[`JuMP.has_duals`](@ref JuMP.has_duals(::InfiniteModel)) to check
whether a result exists before checking the duals. 
    
Thw keyword arguments `kwargs` depend on the transformation backend that is 
being used. The default backend `TranscriptionOpt` uses the keyword 
arguments:
- `result::Int = 1`: indexes the solution result to be queried
- `label::Type{<:AbstractSupportLabel} = PublicLabel`: the label of supports to be returned
- `ndarray::Bool = false`: indicates whether the output should be formatted as an array
By default only the values associated with public supports (i.e., `PublicLabel`s) 
are returned, the full set can be accessed via `label = All`. Moreover, the duals 
of infinite constraints are returned as a list. However, a n-dimensional array 
can be obtained via `ndarray = true` which is handy when the constraint has multiple 
infinite parameter dependencies.

To provide context for the duals, it may be helpful to also query the constraint's 
`parameter_refs` and `supports` which will have a one-to-one correspondence with 
the value(s). It may also be helpful to query via [`transformation_constraint`](@ref) 
to retrieve the constraint(s) that these values are based on. These functions should 
all be called with the same keyword arguments for consistency.

For extensions, this only works if 
[`transformation_constraint`](@ref) has been extended correctly and/or 
[`map_dual`](@ref) has been extended for constraints.

**Example**
```julia-repl
julia> dual(c1)
4-element Array{Float64,1}:
 -42.0
 -42.0
 32.3
 0.0
```
"""
function JuMP.dual(cref::InfOptConstraintRef; kwargs...)
    model = JuMP.owner_model(cref)
    _check_result_is_current(model, JuMP.dual)
    return map_dual(cref, model.backend; kwargs...)
end

# Error redirect for variable call
function JuMP.dual(vref::GeneralVariableRef; kwargs...)
    return JuMP.dual(JuMP.VariableRef(JuMP.Model(), MOI.VariableIndex(1)))
end

################################################################################
#                           LP SENSITIVITY ANALYSIS
################################################################################
"""
    InfOptSensitivityReport

A wrapper `DataType` for `JuMP.SensitivityReport`s in `InfiniteOpt`. 
These are generated based on the transformation backend and should be made via 
the use of [`lp_sensitivity_report`](@ref JuMP.lp_sensitivity_report(::InfiniteModel)). 
Once made these can be indexed to get the sensitivies with respect to variables and/or 
constraints. The indexing syntax for these is: 
```julia
report[ref::[GeneralVariableRef/InfOptConstraintRef]; 
       [label::Type{<:AbstractSupportLabel} = PublicLabel,
       ndarray::Bool = false, kwargs...]]
```

This is enabled for new transformation backends by appropriately 
extending [`transformation_variable`](@ref) and 
[`transformation_constraint`](@ref).

**Fields**
- `opt_report::JuMP.SensitivityReport`: The LP sensitivity captured from the backend.
"""
struct InfOptSensitivityReport 
    opt_report::JuMP.SensitivityReport
end

# Extend Base.getindex for variables on InfOptSensitivityReport
function Base.getindex(s::InfOptSensitivityReport, v::GeneralVariableRef; kwargs...)
    backend = JuMP.owner_model(v).backend
    opt_vref = transformation_variable(v, backend; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> s.opt_report[v], opt_vref)
    else
        return s.opt_report[opt_vref]
    end
end

# Extend Base.getindex for constraints on InfOptSensitivityReport
function Base.getindex(s::InfOptSensitivityReport, c::InfOptConstraintRef; kwargs...)
    backend = JuMP.owner_model(c).backend
    opt_cref = transformation_constraint(c, backend; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> s.opt_report[c], opt_cref)
    else
        return s.opt_report[opt_cref]
    end
end

"""
    JuMP.lp_sensitivity_report(
        backend::AbstractTransformationBackend;
        [atol::Float64 = 1e-8]
        )::InfOptSensitivityReport

Extend `JuMP.lp_sensitivity_report` as appropriate for `backend`. This is
intended as an extension point. For [`JuMPBackend`](@ref)s, this simply 
calls `JuMP.lp_sensitivity_report` on the underlying JuMP model. 
"""
function JuMP.lp_sensitivity_report(
    backend::AbstractTransformationBackend;
    kwargs...
    )
    error("`JuMP.lp_sensitivity_report` not defined for backends of type " *
          "`$(typeof(backend))`.")
end

# JuMPBackend
function JuMP.lp_sensitivity_report(backend::JuMPBackend; atol::Float64 = 1e-8)
    report = JuMP.lp_sensitivity_report(backend.model, atol = atol)
    return InfOptSensitivityReport(report)
end

"""
    JuMP.lp_sensitivity_report(
        model::InfiniteModel;
        [atol::Float64 = 1e-8]
        )::InfOptSensitivityReport

Extends `JuMP.lp_sensitivity_report` to generate and return an LP sensitivity 
report in accordance with the transformation backend. See 
[`InfOptSensitivityReport`](@ref) for syntax details on how to query it. `atol` 
denotes the optimality tolerance and should match that used by the solver to 
compute the basis. Please refer to `JuMP`'s documentation for more technical 
information on interpretting the output of the report.

**Example**
```julia-repl 
julia> report = lp_sensitivity_report(model);

julia> report[x]
(0.0, 0.5)
```
"""
function JuMP.lp_sensitivity_report(model::InfiniteModel; atol::Float64 = 1e-8)
    _check_result_is_current(model, JuMP.lp_sensitivity_report)
    return JuMP.lp_sensitivity_report(model.backend; atol = atol)
end
