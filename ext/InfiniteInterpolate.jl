module InfiniteInterpolate

import JuMP
import InfiniteOpt as infOpt
import Interpolations as IP

# Extend value function to return interpolated infinite variables
function JuMP.value(vref::infOpt.GeneralVariableRef, method::Function, interpolate::Bool; kwargs...)
    irregularGridMethod = Union{typeof(IP.linear_interpolation),
                                typeof(IP.constant_interpolation)}

    infOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)
    if !interpolate
        return infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)
    else
        # Get infinite parameter references for which the variable depends on
        prefs = infOpt.parameter_refs(vref)

        # Get the parameter supports
        Vparams = []
        for pref in prefs
            JuMP.check_belongs_to_model(pref, JuMP.owner_model(vref))
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
        Vparams = Tuple(Vparams)  # Used for constant or linear linear_interpolation

        # Get the variable supports
        Vsupps = infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)

        # TODO: What happens if a user manually defined their support when modelling but wants to use cubic spline interpolation when interpolating?

        # Wrappers for interpolation method
        interpolateFunc(interpMethod::Function, params::Tuple, supps::Vector) = interpMethod(params, supps)

        interpolateFunc(interpMethod::Function, params::Tuple, supps::Matrix) = interpMethod(params, supps)

        # Pass the parameter and variable values to interpolation function
        varFunc = interpolateFunc(method, Vparams, Vsupps)
        return varFunc
    end
end

# Extend for semi-infinite variables
# function JuMP.value(vref::SemiInfiniteVariableRef; interpolate::Bool=false, kwargs...)
#     _check_result_is_current(JuMP.owner_model(vref), JuMP.value)
#     if !interpolate
#         return _get_value(vref, _index_type(vref); kwargs...)
#     else
#         # Get the infinite parameter references that the variable depends on
#         prefs = parameter_refs(vref)
#         _check_parameters_valid(JuMP.owner_model(vref), prefs)

#         # Get the parameter values & supports
#         Vparams = _get_value.(pref, _index_type(vref); kwargs...)
#         Vsupps = supports(vref)

#         # Pass the supports and variable values to interpolation function
#         varFunc = IP.cubic_spline_interpolation(discV, supps)
#         return varFunc
# end
end