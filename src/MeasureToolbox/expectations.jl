################################################################################
#                             EXPECTATION METHODS
################################################################################
## Determine if parameter(s) use a distrubution set
# IndependentParameter and UniDistributionSet
function _has_distribution_set(pref::InfiniteOpt.IndependentParameterRef,
                               set::InfiniteOpt.UniDistributionSet)::Bool
    return true
end

# DependentParameter and MultiDistributionSet
function _has_distribution_set(pref::InfiniteOpt.DependentParameterRef,
                               set::InfiniteOpt.MultiDistributionSet)::Bool
    return true
end

# DependentParameter and CollectionSet
function _has_distribution_set(pref::InfiniteOpt.DependentParameterRef,
                               set::InfiniteOpt.CollectionSet)::Bool
    sets = InfiniteOpt.collection_sets(set)
    return sets[InfiniteOpt._param_index(pref)] isa InfiniteOpt.UniDistributionSet
end

# Fallback
function _has_distribution_set(pref, set)::Bool
    return false
end

# Define coefficient function for expectations
_expect_coeffs(x::Array) = ones(size(x)[end]) ./ size(x)[end]

"""
    expect(expr::JuMP.AbstractJuMPScalar,
           prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef};
           [min_num_supports::Int = DefaultNumSupports])::GeneralVariableRef

Creates a measure that represents the expected value of an expression based on
`prefs`. If `prefs` are not random parameters then this will be equivalent
to the following call:
```julia
1/total_num_supports * support_sum(expr, prefs)
```
Note that min_num_supports should be 0 if a single dependent parameter is given.
Also, note that it is preferred to call [`@expect`](@ref) when `expr` is not
just a single variable reference.

**Example**
```julia-repl
julia> @infinite_parameter(model, x in Normal())
x

julia> @infinite_variable(model, f(x))
f(x)

julia> meas = expect(f, min_num_supports = 2)
expect{x}[f(x)]

julia> expand(meas)
0.5 f(0.6791074260357777) + 0.5 f(0.8284134829000359)
```
"""
function expect(expr::JuMP.AbstractJuMPScalar,
    prefs::Union{InfiniteOpt.GeneralVariableRef, AbstractArray{InfiniteOpt.GeneralVariableRef}};
    min_num_supports::Int = InfiniteOpt.DefaultNumSupports
    )::InfiniteOpt.GeneralVariableRef
    # check the inputs
    InfiniteOpt._check_params(prefs)
    dpref = InfiniteOpt.dispatch_variable_ref(first(prefs))
    # update properly for single dependent parameters
    if prefs isa InfiniteOpt.GeneralVariableRef && dpref isa InfiniteOpt.DependentParameterRef
        if !(min_num_supports in [0, InfiniteOpt.DefaultNumSupports])
            @warn("Cannot specify a nonzero `min_num_supports` for individual " *
                  "dependent parameters.")
        end
        min_num_supports = 0
    end
    # prepare the label
    set = InfiniteOpt._parameter_set(dpref)
    if _has_distribution_set(dpref, set)
        label = InfiniteOpt.WeightedSample
        is_expect = true
    else
        label = InfiniteOpt.All # TODO maybe do something more rigorous
        is_expect = false
    end
    ordered_prefs = InfiniteOpt._make_ordered_vector(prefs)
    length(ordered_prefs) == 1 ? bounds = NaN : bounds = map(e -> NaN, ordered_prefs)
    # make the data
    data = InfiniteOpt.FunctionalDiscreteMeasureData(ordered_prefs, _expect_coeffs,
                                                     min_num_supports, label,
                                                     InfiniteOpt.default_weight,
                                                     bounds, bounds, is_expect)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "expect")
end

"""
    @expect(expr::JuMP.AbstractJuMPScalar,
            prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef};
            [min_num_supports::Int = DefaultNumSupports])::GeneralVariableRef

An efficient wrapper for [`expect`](@ref). Please see its doc string more
information.
"""
macro expect(expr, prefs, args...)
    _error(str...) = InfiniteOpt._macro_error(:integral, (expr, prefs, args...), 
                                              str...)
    extra, kw_args, requestedcontainer = InfiniteOpt._extract_kw_args(args)
    if length(extra) > 0
        _error("Unexpected positional arguments." *
               "Must be of form @expect(expr, prefs, min_num_supports = some_integer).")
    end
    if !isempty(filter(kw -> kw.args[1] != :min_num_supports, kw_args))
        _error("Unexpected keyword arugments.")
    end
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( expect($expression, $prefs; ($(kw_args...))) )
    return esc(mref)
end
