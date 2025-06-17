module InfiniteInterpolate

import JuMP
import InfiniteOpt as infOpt
import Interpolations as IP

"""
    JuMP.value(vref::GeneralVariableRef, method::Function; [kwargs...])

Extend `JuMP.value` to return `vref` as an interpolation object from Interpolations.jl, where the `method` argument specifies the interpolation method to be used.
The currently supported methods are `linear_interpolation`, `constant_interpolation` and `cubic_spline_interpolation`.

**Example**
```julia-repl
julia> zFunc = value(z, cubic_spline_interpolation)
julia> zFunc(5.4)
42.0
```
"""

# Extend value function to return interpolated infinite variables
function JuMP.value(vref::infOpt.GeneralVariableRef, method::Function; kwargs...)
    irregularGridMethod = Union{typeof(IP.linear_interpolation),
                                typeof(IP.constant_interpolation)}
    infOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)

    # Check if interpolation method is supported
    if !(method in (IP.constant_interpolation, IP.linear_interpolation, IP.cubic_spline_interpolation))
        throw(ArgumentError("Unsupported interpolation method: $(method). Supported methods are: constant_interpolation, linear_interpolation, cubic_spline_interpolation."))
    end

    # Get infinite parameter references for which the variable depends on
    prefs = infOpt.parameter_refs(vref)

    if isempty(prefs)
        # If no infinite parameters, return the value directly
        return infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)
    else
        # Get the parameter supports
        Vparams = []
        for pref in prefs
            JuMP.check_belongs_to_model(pref, JuMP.owner_model(vref))

            # If user defined irregular grid points for supports, throw an error if the specified method doesn't support it
            suppsLabel = first(infOpt.core_object(pref).supports)[2]
            if !(infOpt.UniformGrid in suppsLabel) && !(method isa irregularGridMethod)
                throw(ArgumentError("Interpolation method $(method) does not support irregular grids for supports. Please specify a uniform grid or choose a different interpolation method."))
            end

            paramVals = infOpt._get_value(pref, infOpt._index_type(pref); kwargs...)
            numSupps = infOpt.num_supports(pref)
            if method isa irregularGridMethod
                # Can directly pass in a vector of support values
                push!(Vparams, paramVals)
            else
                # Create range for support values
                paramRange = LinRange(paramVals[1], paramVals[end], numSupps)
                push!(Vparams, paramRange)
            end
        end

            # Ensure Vparams is a tuple for interpolation
            Vparams = Tuple(Vparams)

            # Get the variable supports
            Vsupps = infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)

            # Wrappers for interpolation method
            interpolateFunc(interpMethod::Function, params::Tuple, supps::Vector) = interpMethod(params, supps)

            # Wrapper for multi-parameter variables
            interpolateFunc(interpMethod::Function, params::Tuple, supps::Matrix) = interpMethod(params, supps)

            # Pass the parameter and variable values to interpolation function
            varFunc = interpolateFunc(method, Vparams, Vsupps)
            return varFunc
    end
end
end