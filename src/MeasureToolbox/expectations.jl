################################################################################
#                             EXPECTATION METHODS
################################################################################
"""
    generate_expect_data(domain::AbstractInfiniteDomain, 
                         prefs::Union{GeneralVariableRef, Vector{GeneralVariableRef}}, 
                         num_supports::Int; 
                         [kwargs...])::AbstractMeasureData

Generate a concrete instance of `AbstractMeasureData` in accordance with the 
`domain` and infinite parameter(s) `prefs` given for computing the expectation. 
This is intended as an internal method, but should be extended for user defined 
infinite domain types.
"""
function generate_expect_data end 

# Fallback
function generate_expect_data(domain, prefs, num_supports; kwargs...)
    error("Unsupported infinite domain and/or parameter(s) combo for `expect`. ",
          "If this syntax is desired consider using `support_sum` or `integral` ",
          "to carry out the measure.\n\n",
          "If you are extending InfiniteOpt with a ",
          "new infinite domain type, you'll need to extend ",
          "`InfiniteOpt.MeasureToolbox.generate_expect_data`.")
end

# Define coefficient function for expectations
_expect_coeffs(x::Array) = ones(size(x)[end]) ./ size(x)[end]

# Univariate distribution domain 
function generate_expect_data(domain::InfiniteOpt.DistributionDomain, 
    pref::InfiniteOpt.GeneralVariableRef, 
    num_supports;
    kwargs...
    )
    for (k, _) in kwargs
        error("Keyword argument `$k` not supported for expectations over ",
              "distribution domains.")
    end
    return InfiniteOpt.FunctionalDiscreteMeasureData(pref, _expect_coeffs,
                                                     num_supports, 
                                                     InfiniteOpt.WeightedSample,
                                                     InfiniteOpt.NoGenerativeSupports(),
                                                     InfiniteOpt.default_weight,
                                                     NaN, NaN, true)
end

# Multivariate distribution domain 
function generate_expect_data(domain::InfiniteOpt.MultiDistributionDomain, 
    prefs::Vector{InfiniteOpt.GeneralVariableRef}, 
    num_supports;
    kwargs...
    )
    for (k, _) in kwargs
        error("Keyword argument `$k` not supported for expectations over ",
              "distribution domains.")
    end
    bounds = map(e -> NaN, prefs)
    return InfiniteOpt.FunctionalDiscreteMeasureData(prefs, _expect_coeffs,
                                                     num_supports, 
                                                     InfiniteOpt.WeightedSample,
                                                     InfiniteOpt.NoGenerativeSupports(),
                                                     InfiniteOpt.default_weight,
                                                     bounds, bounds, true)
end

# Collection domain 
function generate_expect_data(domain::InfiniteOpt.CollectionDomain, 
    prefs::Vector{InfiniteOpt.GeneralVariableRef}, 
    num_supports;
    kwargs...
    )
    domains = InfiniteOpt.collection_domains(domain)
    if !all(d isa InfiniteOpt.UniDistributionDomain for d in domains)
        error("`Expect` not supported for dependent infinite parameters with ",
              "heterogeneous infinite domains containing non-distribution ",
              "domains.")
    end
    for (k, _) in kwargs
        error("Keyword argument `$k` not supported for expectations over ",
              "distribution domains.")
    end
    bounds = map(e -> NaN, prefs)
    return InfiniteOpt.FunctionalDiscreteMeasureData(prefs, _expect_coeffs,
                                                     num_supports, 
                                                     InfiniteOpt.WeightedSample,
                                                     InfiniteOpt.NoGenerativeSupports(),
                                                     InfiniteOpt.default_weight,
                                                     bounds, bounds, true)
end

# Collection domain of intervals
function generate_expect_data(
    domain::InfiniteOpt.CollectionDomain{InfiniteOpt.IntervalDomain}, 
    pref::InfiniteOpt.GeneralVariableRef, 
    num_supports;
    kwargs...
    )
    domains = InfiniteOpt.collection_domains(domain)
    domain = domains[InfiniteOpt._param_index(pref)]
    return generate_expect_data(domain, pref, num_supports; kwargs...)
end

# Interval domain coefficient function 
function _default_pdf(supp, lb, ub)
    return 1 / (ub - lb)
end

# Univariate interval domain
function generate_expect_data(domain::InfiniteOpt.IntervalDomain, 
    pref::InfiniteOpt.GeneralVariableRef, 
    num_supports;
    pdf::Function = InfiniteOpt.default_weight,
    kwargs...
    )
    for (k, _) in kwargs
        error("Keyword argument `$k` not supported for expectations over ",
              "interval domains.")
    end
    lb = JuMP.lower_bound(domain)
    ub = JuMP.upper_bound(domain)
    if isinf(lb) || isinf(ub)
        error("(Semi-)infinite interval domains are not currently supported ",
              "for expectations.")
    elseif pdf == InfiniteOpt.default_weight
        pdf = (supp) -> _default_pdf(supp, lb, ub)
    end
    return generate_integral_data(pref, lb, ub, UniTrapezoid(), 
                                  num_supports = num_supports,
                                  weight_func = pdf)
end

"""
    expect(expr::JuMP.AbstractJuMPScalar,
           prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef};
           [num_supports::Int = DefaultNumSupports])::GeneralVariableRef

Makes a measure for `expr` based on its expectation with respect to `prefs`. For 
`prefs` with distribution domains this is essentially equivalent to 
```julia
1/total_num_supports * support_sum(expr, prefs, label = WeightedSample)
```
Thus, for these domain types it only considers supports that are added to `prefs` 
via generation on creation (i.e., specifying the `num_supports` keyword when 
creating `prefs`). For incorporating other supports consider 
calling [`integral`](@ref) and using the `weight_func` argument to specify the 
probability density function.

For a single infinite parameter defined over a bounded interval domain the syntax 
becomes:
```julia 
    expect(expr::JuMP.AbstractJuMPScalar,
           prefs::GeneralVariableRef;
           [num_supports::Int = DefaultNumSupports,
           pdf::Function = (supp) -> 1 / (ub - lb)])::GeneralVariableRef
```
The behavior with the default `pdf` is equivalent to evaluating the mean value 
theorem for integrals for `expr` with respect to `pref` using 
[`UniTrapezoid`](@ref). Other density functions can be given via `pdf`. Errors 
if the interval domain is not bounded.

Note that num_supports should be 0 if a single dependent parameter is given.
Also, note that it is preferred to call [`@expect`](@ref) when `expr` is not
just a single variable reference.

**Example**
```julia-repl
julia> @infinite_parameter(model, x ~ Normal(), num_supports = 2)
x

julia> @variable(model, f, Infinite(x))
f(x)

julia> meas = expect(f, x)
𝔼{x}[f(x)]

julia> expand(meas)
0.5 f(0.6791074260357777) + 0.5 f(0.8284134829000359)
```
"""
function expect(expr::JuMP.AbstractJuMPScalar,
    prefs::Union{InfiniteOpt.GeneralVariableRef, AbstractArray{InfiniteOpt.GeneralVariableRef}};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    kwargs...
    )::InfiniteOpt.GeneralVariableRef
    # check the inputs
    InfiniteOpt._check_params(prefs)
    dpref = InfiniteOpt.dispatch_variable_ref(first(prefs))
    # update properly for single dependent parameters
    if prefs isa InfiniteOpt.GeneralVariableRef && dpref isa InfiniteOpt.DependentParameterRef
        if !(num_supports in [0, InfiniteOpt.DefaultNumSupports])
            @warn("Cannot specify a nonzero `num_supports` for individual " *
                  "dependent parameters.")
        end
        num_supports = 0
    end
    # prepare the data
    domain = InfiniteOpt._parameter_domain(dpref)
    vect_prefs = InfiniteOpt.Collections.vectorize(prefs)
    data = generate_expect_data(domain, vect_prefs, num_supports; kwargs...)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "expect")
end

"""
    @expect(expr::JuMP.AbstractJuMPScalar,
            prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef};
            [num_supports::Int = DefaultNumSupports, 
            kwargs...]
            )::GeneralVariableRef

An efficient wrapper for [`expect`](@ref). Please see its doc string more
information.
"""
macro expect(expr, prefs, args...)
    _error(str...) = InfiniteOpt._macro_error(:expect, (expr, prefs, args...), 
                                              str...)
    extra, kwargs, _, _ = InfiniteOpt._extract_kwargs(args)
    if length(extra) > 0
        _error("Unexpected positional arguments." *
               "Must be of form @expect(expr, prefs, kwargs...).")
    end
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( expect($expression, $prefs; ($(kwargs...))) )
    return esc(mref)
end

"""
    𝔼(expr::JuMP.AbstractJuMPScalar,
      prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}};
      [num_supports::Int = DefaultNumSupports, 
      kwargs...]
      )::GeneralVariableRef)

A convenient wrapper for [`expect`](@ref). The unicode symbol `𝔼` is produced by 
`\\bbE`.
"""
function 𝔼(expr::JuMP.AbstractJuMPScalar,
    prefs::Union{InfiniteOpt.GeneralVariableRef, AbstractArray{InfiniteOpt.GeneralVariableRef}};
    num_supports::Int = InfiniteOpt.DefaultNumSupports
    )::InfiniteOpt.GeneralVariableRef
    return expect(expr, prefs, num_supports = num_supports)
end

"""
    @𝔼(expr::JuMP.AbstractJuMPScalar,
       prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef};
       [num_supports::Int = DefaultNumSupports],
       kwargs...)::GeneralVariableRef

A convenient wrapper for [`@expect`](@ref). The unicode symbol `𝔼` is produced by 
`\\bbE`.
"""
macro 𝔼(expr, prefs, args...)
    _error(str...) = InfiniteOpt._macro_error(:𝔼, (expr, prefs, args...), 
                                              str...)
    extra, kwargs, _, _ = InfiniteOpt._extract_kwargs(args)
    if length(extra) > 0
        _error("Unexpected positional arguments." *
               "Must be of form @𝔼(expr, prefs, kwargs...).")
    end
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( expect($expression, $prefs; ($(kwargs...))) )
    return esc(mref)
end
