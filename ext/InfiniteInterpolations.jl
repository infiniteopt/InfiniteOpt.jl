module InfiniteInterpolations

import JuMP
import InfiniteOpt
import Interpolations as IP

# Define convenient type unions
const _IrregularGridMethods = Union{
    typeof(IP.linear_interpolation),
    typeof(IP.constant_interpolation)
}
const _RefTypes = Union{
    InfiniteOpt.GeneralVariableRef,
    JuMP.GenericAffExpr{Float64, InfiniteOpt.GeneralVariableRef},
    JuMP.GenericQuadExpr{Float64, InfiniteOpt.GeneralVariableRef},
    JuMP.GenericNonlinearExpr{InfiniteOpt.GeneralVariableRef},
    InfiniteOpt.InfOptConstraintRef
}

"""
    JuMP.value(
        ref::Union{GeneralVariableRef, JuMP.AbstractJuMPScalar, InfOptConstraintRef},
        method::Union{Interpolations.InterpolationType, Interpolations.Degree};
        [kwargs...]
    )::Interpolations.Extrapolation

Extend `JuMP.value` to return `ref` as an interpolation object from
Interpolations.jl, based on `method` which specifies the interpolation method. 
Currently supported method(s) are:
- `constant_interpolation` or `Constant()`
- `linear_interpolation` or `Linear()`
- `cubic_spline_interpolation` or `Cubic()`

All methods support equidistant grid points. However, nonequidistant discretization grids
are not compatible with cubic splines.

**Examples**
```julia-repl
julia> y_interp_func = value(y, cubic_spline_interpolation)

julia> y_interp_func(5.4);
42.0

julia> y_interp_func2 = value(y, Cubic())

julia> y_interp_func2(5.4)
42.0
```
"""
function JuMP.value(
    obj::_RefTypes,
    method::IP.InterpolationType; 
    kwargs...
    )
    error("Unsupported interpolation type: $(method).",
          "Supported interpolation types are: linear_interpolation, ",
          "constant_interpolation, and cubic_spline_interpolation.")
end
function JuMP.value(
    ref::_RefTypes,
    interp_method::typeof(IP.cubic_spline_interpolation);
    kwargs...
    )
    # get the variable values
    ref_values = JuMP.value(ref; kwargs...)
    isempty(InfiniteOpt.parameter_group_int_indices(ref)) && return ref_values
    # get infinite parameter references for which the variable depends on
    prefs = InfiniteOpt.parameter_refs(ref)
    supp_ranges = Tuple(
        begin
            if length(pref) != 1
                error("Cannot interpolate over dependent infinite parameters.")
            end
            s = JuMP.value(pref; kwargs...)
            if !all(isapprox(s[i+1] - s[i], s[2] - s[1]) for i in 1:length(s)-1)
                error("Cubic spline interpolation requires uniform grids for supports.")
            end
            LinRange(s[1], s[end], length(s))
        end
        for pref in prefs
    )
    # return the interpolation
    return interp_method(supp_ranges, ref_values)
end
function JuMP.value(
    ref::_RefTypes,
    interp_method::_IrregularGridMethods; 
    kwargs...
    )
    # get the variable values
    ref_values = JuMP.value(ref; kwargs...)
    isempty(InfiniteOpt.parameter_group_int_indices(ref)) && return ref_values
    # get infinite parameter references for which the variable depends on
    prefs = InfiniteOpt.parameter_refs(ref)
    supp_values = Tuple(
        begin
            if length(pref) != 1
                error("Cannot interpolate over dependent infinite parameters.")
            end
            JuMP.value(pref; kwargs...)
        end
        for pref in prefs
    )
    # return the interpolation
    return interp_method(supp_values, ref_values)
end

# Interpolation degrees
function JuMP.value(
    ref::_RefTypes,
    degree::IP.Degree; 
    kwargs...
    )
    error("Unsupported interpolation degree: $(degree). Supported ",
          "interpolation degrees are: Linear(), Constant(), and Cubic().")
end
for (degree, method) in zip(
    (:Linear, :Constant, :Cubic),
    (:linear_interpolation, :constant_interpolation, :cubic_spline_interpolation)
    )
    @eval begin
        function JuMP.value(
            ref::_RefTypes,
            degree::IP.$degree;
            kwargs...
            )
            return JuMP.value(ref, IP.$method; kwargs...)
        end
    end
end

# Provide fallback for map_value_to_start (kind of type piracy, but there isn't a great alternative)
function InfiniteOpt.map_value_to_start(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::InfiniteOpt.AbstractTransformationBackend;
    degree::IP.Degree = IP.Linear()
    )
    interp = JuMP.value(vref, degree)
    return (p...) -> interp(p...)
end

end # end of module
