module InfiniteInterpolate

import JuMP
import InfiniteOpt
import Interpolations as IP

const _irregularGridMethods = Union{typeof(IP.linear_interpolation),
                            typeof(IP.constant_interpolation)}

const _convenienceConstructors = Union{typeof(IP.linear_interpolation),
                        typeof(IP.constant_interpolation),
                        typeof(IP.cubic_spline_interpolation)}

"""
    JuMP.value(vref::GeneralVariableRef,
        method::Interpolations.InterpolationType);
        [kwargs...])

Extend `JuMP.value` to return `vref` as an interpolation object from
Interpolations.jl, based on `method` which specifies the interpolation method. 
Currently supported method(s) are:
- `linear_interpolation`
- `constant_interpolation`
- `cubic_spline_interpolation`

All methods support equidistant grid points. However, irregular grid points can 
only be used with `linear_interpolation` and `constant_interpolation`.

```julia
JuMP.value(vref::GeneralVariableRef, degree::Interpolations.Degree; kwargs...)
```
Extend `JuMP.value` to return `vref` as an interpolation object from 
Interpolations.jl, based on `degree` which specifies the degree of interpolation.
The currently supported degrees are:
- `Linear()`
- `Constant()`
- `Cubic()`

All methods support equidistant grid points. However, irregular grid points can 
only be used with `Linear()` and `Constant()`.

**Examples**
```julia-repl
julia> zFunc = value(z, cubic_spline_interpolation)

julia> zFunc(5.4);
42.0

julia> zFunc2 = value(z, Cubic())

julia> zFunc2(5.4)
42.0
```
"""
function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
                    method::IP.InterpolationType; kwargs...)
    throw(ArgumentError("Unsupported interpolation type: $(method). 
        Supported interpolation types are: linear_interpolation, 
        constant_interpolation, and cubic_spline_interpolation."))
end

function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
        interpMethod::typeof(IP.cubic_spline_interpolation);
        kwargs...)
    InfiniteOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)

    # Get infinite parameter references for which the variable depends on
    prefs = InfiniteOpt.parameter_refs(vref)

    if isempty(prefs)
        # If no infinite parameters, return the value directly
        return InfiniteOpt._get_value(vref,
                                      InfiniteOpt._index_type(vref);
                                      kwargs...)
    else
        # Get the parameter supports
        Vparams = []
        for pref in prefs
            # If user defined irregular grid points for supports, throw an error
            suppsLabel = first(InfiniteOpt.core_object(pref).supports)[2]
            if !(InfiniteOpt.UniformGrid in suppsLabel)
                throw(ArgumentError("Interpolation method $(interpMethod) does 
                not support irregular grids for supports. Please specify a 
                uniform grid or choose a different interpolation method."))
            end

            # Create an equidistant range for support values
            paramVals = InfiniteOpt._get_value(pref,
                                    InfiniteOpt._index_type(pref);
                                    kwargs...)
            numSupps = length(paramVals)
            paramRange = LinRange(paramVals[1], paramVals[end], numSupps)
            push!(Vparams, paramRange)
        end

        # Ensure Vparams is a tuple for interpolation
        Vparams = Tuple(Vparams)

        # Get the variable supports
        Vsupps = InfiniteOpt._get_value(vref,
                                        InfiniteOpt._index_type(vref);
                                        kwargs...)

        # Pass the parameter and variable values to interpolation function
        varFunc = interpMethod(Vparams, Vsupps)
        return varFunc
    end
end

function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
    interpMethod::_irregularGridMethods; kwargs...)
    InfiniteOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)

    # Get infinite parameter references for which the variable depends on
    prefs = InfiniteOpt.parameter_refs(vref)

    if isempty(prefs)
        # If no infinite parameters, return the value directly
        return InfiniteOpt._get_value(vref, InfiniteOpt._index_type(vref); kwargs...)
    else
        # Get the parameter supports
        Vparams = []
        for pref in prefs
            paramVals = InfiniteOpt._get_value(pref,
                                               InfiniteOpt._index_type(pref);
                                               kwargs...)
            # Directly pass in a vector of support values, which may be irregularly spaced
            push!(Vparams, paramVals)
        end

        # Ensure Vparams is a tuple for interpolation
        Vparams = Tuple(Vparams)

        # Get the variable supports
        Vsupps = InfiniteOpt._get_value(vref, 
                                        InfiniteOpt._index_type(vref); 
                                        kwargs...)

        # Pass the parameter and variable values to interpolation function
        varFunc = interpMethod(Vparams, Vsupps)
        return varFunc
    end
end

# Fallback for unsupported interpolation degrees
function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
    degree::IP.Degree; kwargs...)
    throw(ArgumentError("Unsupported interpolation degree: $(degree). Supported 
    interpolation degrees are: Linear(), Constant(), and Cubic()."))
end

function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
                    degree::IP.Linear;
                    kwargs...)
    return JuMP.value(vref, IP.linear_interpolation; kwargs...) 
end

function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
                    degree::IP.Constant;
                    kwargs...)
    return JuMP.value(vref, IP.constant_interpolation; kwargs...) 
end

function JuMP.value(vref::InfiniteOpt.GeneralVariableRef,
                    degree::IP.Cubic;
                    kwargs...)
    return JuMP.value(vref, IP.cubic_spline_interpolation; kwargs...) 
end
end