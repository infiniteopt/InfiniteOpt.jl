################################################################################
#                             SUPPORT SUM METHODS
################################################################################
# Define coefficient function for summing over supports
_support_sum_coeffs(x::Array) = ones(size(x)[end])

"""
    support_sum(expr::JuMP.AbstractJuMPScalar,
                params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}}
                )::GeneralVariableRef

Creates a measure that represents the sum of the expression over a parameter(s)
using all of its supports. Also, note that it is preferred to call
[`@support_sum`](@ref) when `expr` is not just a single variable reference.

**Example**
```julia-repl
julia> @infinite_parameter(model, x in [0, 1], supports = [0.3, 0.7])
x

julia> @infinite_variable(model, f(x))
f(x)

julia> meas = support_sum(f, x)
support_sum{x}[f(x)]

julia> expand(meas)
f(0.3) + f(0.7)
```
"""
function support_sum(expr::JuMP.AbstractJuMPScalar,
    prefs::Union{InfiniteOpt.GeneralVariableRef, AbstractArray{InfiniteOpt.GeneralVariableRef}}
    )::InfiniteOpt.GeneralVariableRef
    # make the data
    ordered_prefs = InfiniteOpt._make_ordered_vector(prefs)
    length(prefs) == 1 ? bounds = NaN : bounds = map(e -> NaN, ordered_prefs)
    data = InfiniteOpt.FunctionalDiscreteMeasureData(ordered_prefs, _support_sum_coeffs,
                                                     0, InfiniteOpt.All,
                                                     InfiniteOpt.default_weight,
                                                     bounds, bounds, false)
    # make the measure
    return InfiniteOpt.measure(expr, data, name = "support_sum")
end

"""
    @support_sum(expr::JuMP.AbstractJuMPScalar,
                 params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}}
                 )::GeneralVariableRef

An efficient wrapper for [`support_sum`](@ref) please see its doc string for
more information.
"""
macro support_sum(expr, params)
    expression = :( JuMP.@expression(InfiniteOpt._Model, $expr) )
    mref = :( support_sum($expression, $params) )
    return esc(mref)
end
