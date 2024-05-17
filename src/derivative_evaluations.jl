################################################################################
#                           DERIVATIVE METHOD DATATYPES
################################################################################
"""
    OrthogonalCollocation{Q <: MeasureToolbox.AbstractUnivariateMethod
                          } <: GenerativeDerivativeMethod 

A `DataType` for storing information about orthogonal collocation over finite 
elements to approximate derivatives. The constructor is of the form:
```
    OrthogonalCollocation(num_nodes::Int, 
                          [quad::AbstractUnivariateMethod = GaussLobatto])
```
Note that if `num_nodes > 2` and the problem contains control variables, then 
[`constant_over_collocation`](@ref constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)) 
should probably be used on those variables.

Note that only 1st order derivatives are supported. If higher order derivatives 
are given, then they are recursively reformulated into a system of 1st order 
differential equations which are then discretized using orthogonal collocation. 

**Fields**
- `num_nodes::Int`: The number of collocation points (nodes) per finite element.
- `quadrature_method::Q`: The quadrature method uses to choose the collocation points.
"""
struct OrthogonalCollocation{Q <: MeasureToolbox.AbstractUnivariateMethod
                             } <: GenerativeDerivativeMethod 
    num_nodes::Int
    quadrature_method::Q
    function OrthogonalCollocation{Q}(
        num_nodes::Int, 
        quad::Q
        ) where {Q <: MeasureToolbox.AbstractUnivariateMethod}
        return new{Q}(num_nodes, quad)
    end
end

# Define Lobatto constructor
function OrthogonalCollocation(
    num_nodes::Int, 
    quad::Q
    ) where {Q <: MeasureToolbox.GaussLobatto}
    num_nodes >= 2 || error("Must specify at least 2 collocation points (i.e., " *
                            "the bounds of each support interval with no internal " * 
                            "support nodes).")
    return OrthogonalCollocation{Q}(num_nodes, quad)
end

# Define default constructor
function OrthogonalCollocation(num_nodes::Int)
    return OrthogonalCollocation(num_nodes, MeasureToolbox.GaussLobatto())
end

"""
    allows_high_order_derivatives(method::AbstractDerivativeMethod)::Bool

Returns a `Bool` on whether `method` supports derivatives with orders higher than 
one. If `false` is returned, then higher order derivatives will be automatically 
reformulated into a system of first order derivatives. This is intended as an 
internal method and needs to be extended for any new `AbstractDerivativeMethod`.
"""
function allows_high_order_derivatives(method::AbstractDerivativeMethod)
    error("`allows_high_order_derivatives` not implemented for derivative method " * 
          "of type $(typeof(method)). If you are writing an extension, you need " * 
          "to extend `allows_high_order_derivatives` appropriately.")
end
allows_high_order_derivatives(::FiniteDifference) = true
allows_high_order_derivatives(::OrthogonalCollocation{MeasureToolbox.GaussLobatto}) = false

################################################################################
#                                HELPER METHODS
################################################################################
"""
    generative_support_info(method::AbstractDerivativeMethod)::AbstractGenerativeInfo

Return the [`AbstractGenerativeInfo`](@ref) associated with `method`. This is 
intended as an internal method and should be extended for user-defined derivative 
methods are [`GenerativeDerivativeMethod`](@ref)s.
"""
function generative_support_info(method::AbstractDerivativeMethod)
    error("`generative_support_info` is undefined for derivative method type " * 
          "`$(typeof(method))`.")
end

# NonGenerativeDerivativeMethods
function generative_support_info(method::NonGenerativeDerivativeMethod)
    return NoGenerativeSupports()
end

# OrthogonalCollocation with GaussLabatto
function generative_support_info(
    method::OrthogonalCollocation{MeasureToolbox.GaussLobatto}
    )
    num_nodes = method.num_nodes
    if num_nodes == 2
        return NoGenerativeSupports()
    else
        (nodes, _) = FastGaussQuadrature.gausslobatto(method.num_nodes)
        deleteat!(nodes, length(nodes))
        deleteat!(nodes, 1)
        return UniformGenerativeInfo(nodes, MeasureToolbox.InternalGaussLobatto, -1, 1)
    end
end

# OrthogonalCollocation fallback
function generative_support_info(::OrthogonalCollocation{Q}) where {Q}
    error("Invalid `OrthogonalCollocation` with unknown quadrature method `$Q`.")
end

"""
    support_label(method::GenerativeDerivativeMethod)

Return the support label associated with `method` if there is one, errors otherwise. 
This depends on [`generative_support_info`](@ref generative_support_info(::AbstractDerivativeMethod)) 
being defined for the type of `method`.
"""
function support_label(method::AbstractDerivativeMethod)
    error("`support_label` not defined for derivative methods of type " *
          "`$(typeof(method))`.")
end

# Extend support_label for GenerativeDerivativeMethods
function support_label(method::GenerativeDerivativeMethod)
    return support_label(generative_support_info(method))
end

"""
    make_reduced_expr(vref::GeneralVariableRef, pref::GeneralVariableRef, 
                      support::Float64, write_model::Union{InfiniteModel, JuMP.Model})

Given the argument variable `vref` and the operator parameter `pref` from a 
derivative, build and return the reduced expression in accordance to the support 
`support` with respect to `pref`. New point/semi-infinite variables will be written to 
`write_model`. This is solely intended as a helper function for derivative 
evaluation.
"""
function make_reduced_expr(
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef, 
    support::Float64, 
    write_model::JuMP.AbstractModel
    )
    return make_reduced_expr(vref, _index_type(vref), pref, support, write_model)
end

# MeasureIndex
function make_reduced_expr(
    mref, 
    ::Type{MeasureIndex}, 
    pref, 
    support, 
    write_model
    )
    data = DiscreteMeasureData(pref, [1], [support], InternalLabel, # NOTE the label will not be used
                               default_weight, NaN, NaN, false)
    return expand_measure(mref, data, write_model)
end

# InfiniteVariableIndex/DerivativeIndex
function make_reduced_expr(
    vref, 
    ::Union{Type{InfiniteVariableIndex}, Type{DerivativeIndex}, Type{ParameterFunctionIndex}}, 
    pref, 
    support, 
    write_model
    )
    prefs = parameter_list(vref)
    # only one parameter so we have to make a point variable (we know that this must be pref)
    if length(prefs) == 1
        return make_point_variable_ref(write_model, vref, [support])
    # there are other parameters so make semi-infinite variable
    else 
        pindex = findfirst(isequal(pref), prefs)
        return make_semi_infinite_variable_ref(write_model, vref, [pindex], 
                                               [support])
    end
end

# SemiInfiniteVariableIndex
function make_reduced_expr(
    vref, 
    ::Type{SemiInfiniteVariableIndex}, 
    pref, 
    support, 
    write_model
    )
    # get the preliminary info
    dvref = dispatch_variable_ref(vref)
    ivref = infinite_variable_ref(vref)
    var_prefs = parameter_list(dvref)
    orig_prefs = parameter_list(ivref)
    eval_supps = eval_supports(dvref)
    pindex = findfirst(isequal(pref), orig_prefs)
    # we only have 1 parameter so we need to make a point variable
    if length(var_prefs) == 1
        processed_support = _make_point_support(orig_prefs, eval_supps, pindex, 
                                                support)
        return make_point_variable_ref(write_model, ivref, processed_support)
    # otherwise we need to make another semi-infinite variable
    else
        indices = collect(keys(eval_supps))
        vals = map(k -> eval_supps[k], indices)
        push!(indices, pindex)
        push!(vals, support)
        return make_semi_infinite_variable_ref(write_model, ivref, indices, vals)
    end
end

################################################################################
#                        EVALUATE_DERIVIATIVE DEFINITIONS
################################################################################
"""
    evaluate_derivative(dref::GeneralVariableRef, 
                        vref::GeneralVariableRef,
                        method::AbstractDerivativeMethod,
                        write_model::JuMP.AbstractModel)::Vector{JuMP.AbstractJuMPScalar}

Build expressions for derivative `dref` evaluated in accordance with `method`. Here, 
`vref` is normally the [`derivative_argument`](@ref), but for derivative methods that 
do not support higher order derivatives, `vref` will be substituted with appropriate placeholder 
variables such that `dref` can be reformulated as a first derivative.
The expressions are of the form `lhs - rhs`, where `lhs` is a function of derivatives evaluated at some supports
for certain infinite parameter, and `rhs` is a function of the derivative arguments
evaluated at some supports for certain infinite parameter. For example, for finite difference
methods at point `t = 1`, `lhs` is `Δt * ∂/∂t[T(1)]`, and `rhs` could be `T(1+Δt) - T(1)` in
case of forward difference mode. This is intended as a helper function for `evaluate`, which 
will take the the expressions generated by this method and generate constraints that approximate
the derivative values by setting the expressions as 0. However, one can extend this function 
to encode custom methods for approximating derivatives. This should invoke 
`add_derivative_supports` if the method is generative and users will likely find 
it convenient to use `make_reduced_expr`.
"""
function evaluate_derivative(
    dref::GeneralVariableRef,  
    vref::GeneralVariableRef,                      
    method::AbstractDerivativeMethod,                         
    write_model::JuMP.AbstractModel)
    error("`evaluate_derivative` not defined for derivative method of type " * 
          "$(typeof(method)).")
end

## Define helper methods for finite difference 
# Forward
function _make_difference_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef,
    order::Int,
    expr_idx::Int, 
    supps::Vector{Float64}, # ordered
    write_model::JuMP.AbstractModel, 
    type::Forward
    )
    idx = expr_idx + 1 # support index of appropriate derivative variable
    supp_product = prod(supps[idx + k] - supps[idx + k - 1] for k in 1:order)
    return @_expr(
        make_reduced_expr(dref, pref, supps[idx], write_model) * supp_product - 
        sum((-1)^(order-k) * binomial(order, k) * make_reduced_expr(vref, pref, supps[idx + k], write_model) for k in 0:order)
        )
end

# Central
function _make_difference_expr(
    dref::GeneralVariableRef, 
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef,
    order::Int,
    expr_idx::Int, 
    supps::Vector{Float64}, # ordered
    write_model::JuMP.AbstractModel, 
    type::Central
    )
    if order == 1
        idx = expr_idx + 1 # support index of appropriate derivative variable
        return @_expr(
            make_reduced_expr(dref, pref, supps[idx], write_model) * (supps[idx+1] - supps[idx-1]) -
            make_reduced_expr(vref, pref, supps[idx+1], write_model) +
            make_reduced_expr(vref, pref, supps[idx-1], write_model)
            )
    elseif iseven(order)
        offset = order ÷ 2
        idx = expr_idx + offset # support index of appropriate derivative variable
        supp_product = prod(supps[idx - k + 1 + offset] - supps[idx - k + offset] for k in 1:order)
        return @_expr(
            make_reduced_expr(dref, pref, supps[idx], write_model) * supp_product - 
            sum((-1)^(k) * binomial(order, k) * make_reduced_expr(vref, pref, supps[idx - k + offset], write_model) for k in 0:order)
            )
    else
        error("Central difference with odd derivative orders other than 1 is not supported.")
    end
end

# Backward
function _make_difference_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef,
    order::Int,
    expr_idx::Int, 
    supps::Vector{Float64}, # ordered
    write_model::JuMP.AbstractModel,
    type::Backward
    )
    idx = expr_idx + order # support index of appropriate derivative variable
    supp_product = prod(supps[idx - k + 1] - supps[idx - k] for k in 1:order)
    return @_expr(
        make_reduced_expr(dref, pref, supps[idx], write_model) * supp_product - 
        sum((-1)^(k) * binomial(order, k) * make_reduced_expr(vref, pref, supps[idx - k], write_model) for k in 0:order)
        )
end

# Fallback
function _make_difference_expr(dref, vref, pref, order, idx, supps, model, type)
    error("Undefined `FiniteDifference` technique `$(typeof(type))`.")
end

# evaluate_derivative for FiniteDifference
function evaluate_derivative(
    dref::GeneralVariableRef, 
    vref::GeneralVariableRef,
    method::FiniteDifference, 
    write_model::JuMP.AbstractModel
    )
    # gather the arugment and parameter 
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    # get the supports and check validity
    ordered_supps = sort!(supports(pref, label = All))
    n_supps = length(ordered_supps)
    if method.technique isa Central && iseven(order)
        n_exprs = n_supps - order
    else
        n_exprs = n_supps - order - 1
    end
    n_exprs < 1 && error("$(pref) does not have enough supports for derivative evaluation of $(dref).")
    # make the expressions
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, n_exprs)
    for i in eachindex(exprs)
        exprs[i] = _make_difference_expr(dref, vref, pref, order, i, ordered_supps, 
                                         write_model, method.technique)
    end
    if method.add_boundary_constraint && method.technique isa Forward
        push!(exprs, _make_difference_expr(dref, vref, pref, order, 0, ordered_supps, 
                                           write_model, Forward()))
    elseif method.add_boundary_constraint && method.technique isa Backward
        push!(exprs, _make_difference_expr(dref, vref, pref, order, n_supps - order, ordered_supps, 
                                           write_model, Backward()))
    end
    return exprs
end

# Evaluate_derivative for OrthogonalCollocation with Labatto
# This is based on the collocation scheme presented in "Nonlinear Modeling, 
# Estimation and Predictive Control in APMonitor"
function evaluate_derivative(
    dref::GeneralVariableRef, 
    vref::GeneralVariableRef,
    method::OrthogonalCollocation{MeasureToolbox.GaussLobatto}, 
    write_model::JuMP.AbstractModel
    )
    # gather the arugment and parameter 
    pref = operator_parameter(dref)
    # ensure that the internal supports are added in case they haven't already 
    add_generative_supports(pref) # this checks for enough supports
    # gather the supports and indices of the internal supports
    all_supps = supports(pref, label = All) # sorted by virtue of SortedDict
    bool_inds = map(s -> !(support_label(method) in s), 
                    values(_parameter_supports(dispatch_variable_ref(pref))))
    interval_inds = findall(bool_inds)
    n_nodes = method.num_nodes - 2 # we want the # of internal nodes
    # compute M matrix based on the node values and generate expressions using 
    # the M matrix and `make_reduced_expr`
    exprs = Vector{JuMP.AbstractJuMPScalar}(undef, length(all_supps) - 1)
    counter = 1
    for i in Iterators.take(eachindex(interval_inds), length(interval_inds) - 1)
        # collect the supports
        lb = all_supps[interval_inds[i]]
        i_nodes = all_supps[interval_inds[i]+1:interval_inds[i+1]] .- lb
        # build the matrices in memory order (column-wise using the transpose form)
        M1t = Matrix{Float64}(undef, n_nodes+1, n_nodes+1)
        M2t = Matrix{Float64}(undef, n_nodes+1, n_nodes+1)
        for j in eachindex(i_nodes) # support index
            for k in eachindex(i_nodes) # polynomial index
                M1t[k, j] = k * i_nodes[j]^(k-1)
                M2t[k, j] = i_nodes[j]^k
            end
        end
        Minvt = M1t \ M2t
        for j in eachindex(i_nodes)
            exprs[counter] = @_expr(
                sum(Minvt[k, j] *  make_reduced_expr(dref, pref, i_nodes[k] + lb, write_model) for k in eachindex(i_nodes)) - 
                make_reduced_expr(vref, pref, i_nodes[j] + lb, write_model) + 
                make_reduced_expr(vref, pref, lb, write_model)
                )
            counter += 1
        end
    end
    return exprs
end

################################################################################
#                              EVALUATION METHODS
################################################################################
"""
    evaluate(dref::DerivativeRef)::Nothing

Numerically evaluate `dref` by computing its auxiliary derivative constraints 
(e.g., collocation equations) and add them to the model. For normal usage, it is 
recommended that this method not be called directly and instead have TranscriptionOpt 
handle these equations. Errors if `evaluate_derivative` is not 
defined for the derivative method employed.

The resulting constraints can be accessed via `derivative_constraints`.

**Example**
```julia-repl
julia> m = InfiniteModel(); @infinite_parameter(m, t in [0,2]); @variable(m, T, Infinite(t));

julia> dref = @deriv(T,t)
∂/∂t[T(t)]

julia> add_supports(t, [0, 0.5, 1, 1.5, 2])

julia> evaluate(dref)

julia> derivative_constraints(dref)
Feasibility
4-element Array{InfOptConstraintRef,1}:
 0.5 ∂/∂t[T(t)](0.5) - T(0.5) + T(0) = 0.0
 0.5 ∂/∂t[T(t)](1) - T(1) + T(0.5) = 0.0
 0.5 ∂/∂t[T(t)](1.5) - T(1.5) + T(1) = 0.0
 0.5 ∂/∂t[T(t)](2) - T(2) + T(1.5) = 0.0
```
"""
function evaluate(dref::DerivativeRef)
    if !has_derivative_constraints(dref)
        # collect the basic info
        method = derivative_method(dref)
        model = JuMP.owner_model(dref)
        constr_list = _derivative_constraint_dependencies(dref)
        gvref = GeneralVariableRef(model, JuMP.index(dref).value, DerivativeIndex)
        # recursively build 1st order derivatives if necessary
        vref = derivative_argument(dref)
        if !allows_high_order_derivatives(method) && derivative_order(dref) > 1 
            pref = operator_parameter(dref)
            vref = _build_add_derivative(vref, pref, derivative_order(dref) - 1)
            evaluate(vref)
        end
        # get the expressions 
        exprs = evaluate_derivative(gvref, vref, method, model)
        # add the constraints
        for expr in exprs
            cref = JuMP.@constraint(model, expr == 0)
            push!(constr_list, JuMP.index(cref))
        end
        # update the parameter status 
        pref = dispatch_variable_ref(operator_parameter(dref))
        _set_has_derivative_constraints(pref, true)
    end
    return
end

"""
    evaluate_all_derivatives!(model::InfiniteModel)::Nothing

Evaluate all the derivatives in `model` by adding the corresponding auxiliary 
equations to `model`. See [`evaluate`](@ref) for more information.

**Example**
```julia-repl
julia> m = InfiniteModel();

julia> @infinite_parameter(m, t in [0,2], supports = [0, 1, 2]);

julia> @infinite_parameter(m, x in [0,1], supports = [0, 0.5, 1]);

julia> @variable(m, T, Infinite(x, t));

julia> dref1 = @deriv(T, t); dref2 = @deriv(T, x^2);

julia> evaluate_all_derivatives!(m)

julia> print(m)
Feasibility
Subject to
 ∂/∂t[T(x, t)](x, 1) - T(x, 1) + T(x, 0) = 0.0, ∀ x ∈ [0, 1]
 ∂/∂t[T(x, t)](x, 2) - T(x, 2) + T(x, 1) = 0.0, ∀ x ∈ [0, 1]
 0.25 ∂²/∂x²[T(x, t)](1, t) - T(1, t) + 2 T(0.5, t) - T(0, t) = 0.0, ∀ t ∈ [0, 2]
```
"""
function evaluate_all_derivatives!(model::InfiniteModel)
    for dref in all_derivatives(model)
        evaluate(dref)
    end
    return
end
