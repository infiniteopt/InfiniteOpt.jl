################################################################################
#                             SUPPORT SUM METHODS
################################################################################
# Define coefficient function for summing over supports
_support_sum_coeffs(x::Array) = ones(size(x)[end])

"""
    support_sum(expr::JuMP.AbstractJuMPScalar,
                params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}};
                label = All
                )::GeneralVariableRef

Creates a measure that represents the sum of the expression over a parameter(s)
using all of its supports corresponding to `label`. Also, note that it is 
preferred to call [`@support_sum`](@ref) when `expr` is not just a 
single variable reference.

**Example**
```julia-repl
julia> @infinite_parameter(model, x in [0, 1], supports = [0.3, 0.7])
x

julia> @variable(model, f, Infinite(x))
f(x)

julia> meas = support_sum(f, x)
support_sum{x}[f(x)]

julia> expand(meas)
f(0.3) + f(0.7)
```
"""
function support_sum(
    expr::JuMP.AbstractJuMPScalar,
    prefs::Union{InfiniteOpt.GeneralVariableRef, AbstractArray{InfiniteOpt.GeneralVariableRef}};
    label = InfiniteOpt.All
    )::InfiniteOpt.GeneralVariableRef
    # make the data
    vect_prefs = InfiniteOpt.Collections.vectorize(prefs)
    length(prefs) == 1 ? bounds = NaN : bounds = map(e -> NaN, vect_prefs)
    data = InfiniteOpt.FunctionalDiscreteMeasureData(vect_prefs, _support_sum_coeffs,
                                                     0, label,
                                                     InfiniteOpt.NoGenerativeSupports(),
                                                     InfiniteOpt.default_weight,
                                                     bounds, bounds, false)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "support_sum")
end

"""
    @support_sum(expr::JuMP.AbstractJuMPScalar,
                 prefs::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}};
                 label = All
                 )::GeneralVariableRef

An efficient wrapper for [`support_sum`](@ref) please see its doc string for
more information.
"""
macro support_sum(expr, prefs, args...)
    _error(str...) = InfiniteOpt._macro_error(:support_sum, (expr, prefs, args...), 
                                              str...)
    extra, kwargs, _, _ = InfiniteOpt._extract_kwargs(args)
    if length(extra) > 0
        _error("Unexpected positional arguments." *
               "Must be of form @support_sum(expr, prefs, kwargs...).")
    end
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( support_sum($expression, $prefs; ($(kwargs...))) )
    return esc(mref)
end
