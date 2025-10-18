################################################################################
#                           DERIVATIVE METHOD DATATYPES
################################################################################
"""
    OrthogonalCollocation{
        Q <: MeasureToolbox.AbstractUnivariateMethod
    } <: GenerativeDerivativeMethod 

A `DataType` for storing information about orthogonal collocation over finite 
elements to approximate derivatives. The constructor is of the form:
```
    OrthogonalCollocation(
        num_nodes::Int, 
        [quad::AbstractUnivariateMethod = GaussLobatto]
    )
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
struct OrthogonalCollocation{
    Q <: MeasureToolbox.AbstractUnivariateMethod
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
    make_reduced_expr(
        vref::GeneralVariableRef, 
        pref::GeneralVariableRef, 
        support::Float64, 
        write_model::Union{InfiniteModel, AbstractTransformationBackend}
        )

    make_reduced_expr(
        vref::GeneralVariableRef, 
        pref::GeneralVariableRef, 
        supports::Vector{Float64}, 
        idx::Int, 
        write_model::Union{InfiniteModel, AbstractTransformationBackend}
        )

Given the argument variable `vref` and the operator parameter `pref` from a 
derivative, build and return the reduced expression in accordance to the support 
`support` with respect to `pref`. New point/semi-infinite variables will be written to 
`write_model`. This is solely intended as a helper function for derivative 
evaluation. Instead of providing the support directly, one can also provide the 
vector of supports `supports` and the `index` that retreives the support of interest; 
this is useful an an extension point for backends that only want the index. 
"""
function make_reduced_expr(
    vref::GeneralVariableRef, 
    pref::GeneralVariableRef, 
    support::Float64, 
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    return make_reduced_expr(vref, _index_type(vref), pref, support, write_model)
end

# unindexed supports (fallback) --> useful for backends that just want the index of the support
function make_reduced_expr(
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    supps::Vector{Float64},
    idx,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    return make_reduced_expr(vref, pref, supps[idx], write_model)
end

# MeasureIndex
function make_reduced_expr(
    mref, 
    ::Type{MeasureIndex}, 
    pref, 
    support, 
    write_model
    )
    data = DiscreteMeasureData(
        pref,
        [1],
        [support],
        InternalLabel, # NOTE the label will not be used
        default_weight, 
        NaN, 
        NaN, 
        false
    )
    return expand_measure(mref, data, write_model)
end

# InfiniteVariableIndex/DerivativeIndex/ParameterFunctionIndex
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
        idx = findfirst(isequal(pref), prefs)
        supp = fill(NaN, length(prefs))
        supp[idx] = support
        return make_semi_infinite_variable_ref(write_model, vref, supp)
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
    eval_supp = eval_support(dvref)
    # we only have 1 parameter so we need to make a point variable
    if length(var_prefs) == 1
        supp = [isnan(s) ? support : s for s in eval_supp]
        return make_point_variable_ref(write_model, ivref, supp)
    # otherwise we need to make another semi-infinite variable
    else
        idx = findfirst(isequal(pref), orig_prefs)
        new_eval_support = copy(eval_supp)
        new_eval_support[idx] = support
        return make_semi_infinite_variable_ref(
            write_model,
            ivref,
            new_eval_support
        )
    end
end

################################################################################
#                        EVALUATE_DERIVIATIVE DEFINITIONS
################################################################################
"""
    make_indexed_derivative_expr(
        dref::GeneralVariableRef,
        vref::GeneralVariableRef,
        pref::GeneralVariableRef,
        order::Int,
        idx,
        supps::Vector{Float64},
        write_model::Union{InfiniteModel, AbstractTransformationBackend},
        method::AbstractDerivativeMethod,
        expr_params...
        )

Produce an expression that numerically approvides the derivative `dref` at a 
particular support value (`supps[idx]`) where `supps` are the ordered support 
values of the derivative parameter `pref`. This is an extension point for 
developers adding a new `method` type. It should be defined in connection with 
an extension to [`derivative_expr_data`](@ref). [`evaluate_derivative`](@ref) 
will then use these to generate derivative approximation equations. 

This considers a derivative of 
order `order` on `vref` taken with respect to `pref`. The resulting expression 
is intended for `write_model` and the method [`make_reduced_expr`](@ref) 
can be helpful for generating semi-infinite and point variables of `vref`.
Additionally, parameter values needed for the approximation (e.g., Δt) can 
be provided as extra arguments (note that these must be scalars). The 
iteration data used to provide `idx` and `expr_params` is created using 
[`derivative_expr_data`](@ref). 
"""
function make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64},
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    method::AbstractDerivativeMethod,
    constants...
    )
    error("`make_indexed_derivative_expr` not defined for derivative method of type " * 
          "`$(typeof(method))` with extra arguments `$(string(constants)[2:end-1])`." * 
          "Developers that encounter this error should extend `make_indexed_derivative_expr` " *
          "and `derivative_expr_data`.")
end

"""
    derivative_expr_data(
        dref::GeneralVariableRef, 
        order::Int,
        supps::Vector{Float64},
        method::AbstractDerivativeMethod,
        )::Tuple

Produce the data needed to generated derivative approximation equations for 
the derivative `dref` in accordance with [`make_indexed_derivative_expr`](@ref). 
This is intended as an extension point when defining new derivative approximation 
methods. It should return a tuple of interators needed to generate the constraint 
expressions with `make_indexed_derivative_expr` and it is used within 
[`evaluate_derivative`](@ref) (given `dref`, `vref`, `pref`, `order`, 
`supps`, `write_model`, and `method`) to generate a constraint expression 
`constr_expr` for each index of `idxs`:
```julia
idxs, extra_itrs... = derivative_expr_data(dref, order, supps, method)
for (idx, extra_data...) in zip(idxs, extra_itrs...)
    constr_expr = make_indexed_derivative_expr(
        dref, 
        vref, 
        pref, 
        order, 
        idx, 
        supps, 
        write_model, 
        method, 
        extra_data...
        )
end
```
"""
function derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    method::AbstractDerivativeMethod,
    )
    error("`derivative_expr_data` not defined for derivative method of type " * 
          "`$(typeof(method))`. It might be because the transformation backend doesn't " *
          "support this particular derivative method. If you are extending InfiniteOpt" * 
          "to have a new derivative method, please extend `derivative_expr_data` if " *
          "possible. See the documentation for details.")
end

# Forward Difference
function make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    ::FiniteDifference{Forward},
    supp_product
    )
    return @_expr(
        make_reduced_expr(dref, pref, supps, idx, write_model) * supp_product - 
        sum((-1)^(order-k) * binomial(order, k) * make_reduced_expr(vref, pref, supps, idx + k, write_model) for k in 0:order)
        )
end
function derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    method::FiniteDifference{Forward},
    )
    # check the number of supports
    n_exprs = length(supps) - order - 1
    n_exprs < 1 && error("$(operator_parameter(dref)) does not have enough supports for derivative evaluation of $(dref).")
    # determine the derivative support indices
    idxs = method.add_boundary_constraint ? (1:n_exprs+1) : (2:n_exprs+1)
    # determine the support products (e.g., Δt) used by the approximation
    supp_products = (prod(supps[idx + k] - supps[idx + k - 1] for k in 1:order) for idx in idxs)
    return idxs, supp_products
end

# Central Difference
function make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    ::FiniteDifference{Central},
    supp_product
    )
    if isone(order)
        return @_expr(
            make_reduced_expr(dref, pref, supps, idx, write_model) * supp_product -
            make_reduced_expr(vref, pref, supps, idx+1, write_model) +
            make_reduced_expr(vref, pref, supps, idx-1, write_model)
            )
    else
        offset = order ÷ 2
        return @_expr(
            make_reduced_expr(dref, pref, supps, idx, write_model) * supp_product - 
            sum((-1)^(k) * binomial(order, k) * make_reduced_expr(vref, pref, supps, idx - k + offset, write_model) for k in 0:order)
            )
    end
end
function derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    ::FiniteDifference{Central},
    )
    # check the number of supports
    if iseven(order)
        n_exprs = length(supps) - order
    else
        n_exprs = length(supps) - order - 1
    end
    n_exprs < 1 && error("$(operator_parameter(dref)) does not have enough supports for derivative evaluation of $(dref).")
    # determine the derivative support indices and the support products
    if isone(order)
        idxs = 2:n_exprs+1
        supp_products = (supps[idx+1] - supps[idx-1] for idx in idxs)
    elseif iseven(order)
        offset = order ÷ 2
        idxs = 1+offset:n_exprs+offset
        supp_products = (prod(supps[idx - k + 1 + offset] - supps[idx - k + offset] for k in 1:order) for idx in idxs)
    else
        error("Central difference with odd derivative orders other than 1 is not supported.")
    end
    return idxs, supp_products
end

# Backward Difference
function make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    ::FiniteDifference{Backward},
    supp_product
    )
    return @_expr(
        make_reduced_expr(dref, pref, supps, idx, write_model) * supp_product - 
        sum((-1)^(k) * binomial(order, k) * make_reduced_expr(vref, pref, supps, idx - k, write_model) for k in 0:order)
        )
end
function derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    method::FiniteDifference{Backward},
    )
    # check the number of supports
    n_exprs = length(supps) - order - 1
    n_exprs < 1 && error("$(operator_parameter(dref)) does not have enough supports for derivative evaluation of $(dref).")
    # determine the derivative support indices
    idxs = method.add_boundary_constraint ? (1+order:n_exprs+order+1) : (1+order:n_exprs+order)
    # determine the support products (e.g., Δt) used by the approximation
    supp_products = (prod(supps[idx - k + 1] - supps[idx - k] for k in 1:order) for idx in idxs)
    return idxs, supp_products
end

# OrthogonalCollocation (Lobatto)
function make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    ::OrthogonalCollocation{MeasureToolbox.GaussLobatto},
    lb_idx,
    coeffs...
    )
    n_coeffs = length(coeffs)
    return @_expr(
        sum(coeffs[k] * make_reduced_expr(dref, pref, supps, lb_idx + k, write_model) for k in 1:n_coeffs) - 
        make_reduced_expr(vref, pref, supps, idx + 1, write_model) + 
        make_reduced_expr(vref, pref, supps, lb_idx, write_model)
        )
end
function derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    method::OrthogonalCollocation{MeasureToolbox.GaussLobatto},
    )
    # determine basic finite element info
    n_coeffs = method.num_nodes - 1 # dimension of M matrix
    # prepare the indices for expression generation
    idxs = 1:length(supps)-1
    lb_idxs = (((idx - 1) ÷ n_coeffs) * n_coeffs + 1 for idx in idxs)
    # compute M matrix for each finite element and gather the coefficients for each equation
    coeffs_list = sizehint!(SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}[], length(supps) - 1)
    for i in 1:n_coeffs:length(supps)-1
        # collect the supports
        lb = supps[i]
        i_nodes = @view(supps[i+1:i+n_coeffs])
        # build the matrices in memory order (column-wise using the transpose form)
        M1t = Matrix{Float64}(undef, n_coeffs, n_coeffs)
        M2t = Matrix{Float64}(undef, n_coeffs, n_coeffs)
        for j in eachindex(i_nodes) # support index
            for k in eachindex(i_nodes) # polynomial index
                @inbounds M1t[k, j] = k * (i_nodes[j] - lb)^(k-1)
                @inbounds M2t[k, j] = (i_nodes[j] - lb)^k
            end
        end
        Minvt = M1t \ M2t
        append!(coeffs_list, eachcol(Minvt))
    end
    coeff_itrs = ((cs[j] for cs in coeffs_list) for j in 1:n_coeffs)
    return idxs, lb_idxs, coeff_itrs...
end

"""
    evaluate_derivative(
        dref::GeneralVariableRef, 
        vref::GeneralVariableRef,
        method::AbstractDerivativeMethod,
        write_model::Union{InfiniteModel, AbstractTransformationBackend}
        )::Vector{JuMP.AbstractJuMPScalar}

Build expressions for derivative `dref` evaluated in accordance with 
`method`. Here, `vref` is normally the [`derivative_argument`](@ref), but 
for derivative methods that do not support higher order derivatives, `vref` 
will be substituted with appropriate placeholder variables such that `dref` 
can be reformulated as a first derivative. The expressions are of the form 
`lhs - rhs`, where `lhs` is a function of derivatives evaluated at some 
supports for certain infinite parameter, and `rhs` is a function of the 
derivative arguments evaluated at some supports for certain infinite parameter. 
For example, for finite difference methods at point `t = 1`, `lhs` is 
`Δt * ∂/∂t[T(1)]`, and `rhs` could be `T(1+Δt) - T(1)` in case of forward 
difference mode. This is intended as a helper function for `evaluate`, which 
will take the the expressions generated by this method and generate 
constraints that approximate the derivative values by setting the expressions 
as 0. 

For developers defining a new type of derivative method, they should extend 
[`derivative_expr_data`](@ref) and [`make_indexed_derivative_expr`](@ref) if at 
all possible. However, if this is not possible, then one can extend 
[`evaluate_derivative`](@ref) directly to encode custom methods for approximating 
derivatives. This should invoke `add_derivative_supports` if the method is 
generative and users will likely find it convenient to use `make_reduced_expr`.
"""
function evaluate_derivative(
    dref::GeneralVariableRef,  
    vref::GeneralVariableRef,                          
    method::AbstractDerivativeMethod,                         
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
    # gather the arugment and parameter 
    pref = operator_parameter(dref)
    order = derivative_order(dref)
    # get and sort supports
    add_generative_supports(pref)
    supps = sort!(supports(pref, label = All))
    # get the necessary data and make the expressions
    idxs, arg_itrs... = derivative_expr_data(dref, order, supps, method)
    return [make_indexed_derivative_expr(dref, vref, pref, order, idx, supps, write_model, method, args...) 
            for (idx, args...) in zip(idxs, arg_itrs...)]
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
