module InfiniteInterpolate

import JuMP
import InfiniteOpt as infOpt
import Interpolations as IP

# Import Revise.jl in the Julia session - allows me to immediately see any changes I make to the extension
# MAKE SURE to use the dev(elop) version of InfiniteOpt or you won't be able to update the package; should be on the how to contribute page
# Pkg.jl might have helpful resources too on how to check out a package

# Need to add InfiniteOpt as a prefix to the functions that are from InfiniteOpt (like I did with 'Interpolations as IP')
# Type piracy - need to work with data types I define myself VS data types from InfiniteOpt; otherwise there's direct conflict

# Might need to backtrack to _get_value potentially
# But _get_value is a private function, so I can't use it
# So we need to find/define an extension point that's public

# Can't dispatch off kwargs, only positional arguments before the semicolon
# Some sort of bool flag to turn interpolation on/off makes sense, but how to make it work?

# NEW IDEA
# GeneralVariableRef
# Extract the infinite parameter refs
# positional argument for interpolation method with the type being the constructors from Interpolations.jl
# Union of all the different types of interpolation methods - Union{Type1, Type2,...} - composite type

# concrete types inherited from abstract types
# Ex: Number is abstract, concrete types like Int64, Float64, etc...
# Function signature refers to the types of the arguments
# Julia will go with the most specific type
# Abstract - Function, concrete type is the specific name of the function
# Might have to wrap with Type{functionName} or smt similar to know we're working with the type of the function, not the function itself
# Try with just linear_interpolation and see if it works; can add more later after this is working
# value(x, linear_interpolation)

# Extend value function to return interpolated infinite variables
function JuMP.value(vref::infOpt.GeneralVariableRef, interpolate::Bool; kwargs...)
    infOpt._check_result_is_current(JuMP.owner_model(vref), JuMP.value)
    if !interpolate
        return infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)
    else
        # Get the infinite parameter references that the variable depends on
        prefs = infOpt.parameter_refs(vref)
        # infOpt._check_parameters_valid(JuMP.owner_model(vref), prefs)

        # Get the parameter values & supports
        Vparams = []
        for pref in prefs
            Vparam = infOpt._get_value(pref, infOpt._index_type(pref); kwargs...)
            push!(Vparams, Vparam)
        end
        
        Vparams = Tuple(Vparams)  # Convert to a tuple for interpolation
        Vsupps = infOpt._get_value(vref, infOpt._index_type(vref); kwargs...)
        numParams = length(Vparams)

        # Pass the supports and variable values to interpolation function
        varFunc = IP.linear_interpolation(Vparams, Vsupps)
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
#         # TODO: add interpolation method as argument?
#         varFunc = IP.cubic_spline_interpolation(discV, supps)
#         return varFunc
# end
end