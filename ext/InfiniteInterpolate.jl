module InfiniteInterpolate

import JuMP
import InfiniteOpt as infOpt
import Interpolations as IP

const _irregularGridMethods = Union{typeof(IP.linear_interpolation),
                            typeof(IP.constant_interpolation)}

const _convenienceConstructors = Union{typeof(IP.linear_interpolation),
                            typeof(IP.constant_interpolation),
                            typeof(IP.cubic_spline_interpolation)}

"""
    JuMP.value(vref::GeneralVariableRef, method::Union{typeof(Interpolations.constant_interpolation), typeof(Interpolations.cubic_spline_interpolation), typeof(Interpolations.linear_interpolation)}; [kwargs...])

Extend `JuMP.value` to return `vref` as an interpolation object from Interpolations.jl, where the `method` argument specifies the interpolation method to be used.
The currently supported methods are `linear_interpolation`, `constant_interpolation` and `cubic_spline_interpolation`.

**Example**
```julia-repl
julia> zFunc = value(z, cubic_spline_interpolation)
julia> zFunc(5.4)
42.0
```
"""
function JuMP.value(vref::infOpt.GeneralVariableRef, interpMethod::_convenienceConstructors; kwargs...)
    infOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)

    # Get infinite parameter references for which the variable depends on
    prefs = infOpt.parameter_refs(vref)

    if isempty(prefs)
        # If no infinite parameters, return the value directly
        return infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)
    else
        # Get the parameter supports
        Vparams = []
        for pref in prefs
            # If user defined irregular grid points for supports, throw an error if the specified method doesn't support it
            suppsLabel = first(infOpt.core_object(pref).supports)[2]
            if !(infOpt.UniformGrid in suppsLabel) && !(interpMethod isa _irregularGridMethods)
                throw(ArgumentError("Interpolation method $(interpMethod) does not support irregular grids for supports. Please specify a uniform grid or choose a different interpolation method."))
            end

            paramVals = infOpt._get_value(pref, infOpt._index_type(pref); kwargs...)
            numSupps = length(paramVals)
            if interpMethod isa _irregularGridMethods
                # Can directly pass in a vector of support values, which may be irregularly spaced
                push!(Vparams, paramVals)
            else
                # Create an equidistant range for support values
                paramRange = LinRange(paramVals[1], paramVals[end], numSupps)
                push!(Vparams, paramRange)
            end
        end

        # Ensure Vparams is a tuple for interpolation
        Vparams = Tuple(Vparams)

        # Get the variable supports
        Vsupps = infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)

        # Pass the parameter and variable values to interpolation function
        varFunc = interpMethod(Vparams, Vsupps)
        return varFunc
    end
end

"""
    JuMP.value(vref::GeneralVariableRef, degree::Interpolations.Degree; [kwargs...])

Extend `JuMP.value` to return `vref` as an interpolation object from Interpolations.jl, where the `degree` argument specifies the degree of interpolation. The corresponding interpolation method is called depending on the degree.
The currently supported degrees are `Linear()`, `Constant()` and `Cubic()`.

**Example**
```julia-repl
julia> zFunc = value(z, Cubic())
julia> zFunc(5.4)
42.0
```
"""
function JuMP.value(vref::infOpt.GeneralVariableRef, degree::IP.Degree; kwargs...)
    if degree == IP.Linear()
        return JuMP.value(vref, IP.linear_interpolation; kwargs...)
    elseif degree == IP.Constant()
        return JuMP.value(vref, IP.constant_interpolation; kwargs...)
    elseif degree == IP.Cubic()
        return JuMP.value(vref, IP.cubic_spline_interpolation; kwargs...)
    else
        throw(ArgumentError("Unsupported interpolation degree type: $(degree). Supported degree types are: Linear(), Constant(), and Cubic()."))
    end
end

# Fallback for unsupported interpolation types
function JuMP.value(vref::infOpt.GeneralVariableRef, method::IP.InterpolationType; kwargs...)
    throw(ArgumentError("Unsupported interpolation type: $(method). Supported interpolation types are: Linear(), Constant(), and Cubic()."))
end

end