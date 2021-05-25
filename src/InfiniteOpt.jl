module InfiniteOpt

# Import and export JuMP 
import Reexport 
Reexport.@reexport using JuMP

# Import the necessary packages.
import MathOptInterface
import Distributions
import DataStructures
import FastGaussQuadrature

# Make useful aliases (note we get MOI and MOIU from JuMP)
const JuMPC = JuMP.Containers
const MOIUC = MOIU.CleverDicts

# Import the Collections module
include("Collections/Collections.jl")
using .Collections: VectorTuple, same_structure, param_type

# Import methods
include("datatypes.jl")
include("infinite_domains.jl")
include("scalar_parameters.jl")
include("array_parameters.jl")
include("variable_basics.jl")
include("infinite_variables.jl")
include("semi_infinite_variables.jl")
include("point_variables.jl")
include("finite_variables.jl")
include("expressions.jl")
include("macro_utilities.jl")
include("measures.jl")

# Import and export MeasureToolbox
include("MeasureToolbox/MeasureToolbox.jl")
Reexport.@reexport using .MeasureToolbox

# Import more methods
include("derivatives.jl")
include("constraints.jl")
include("macros.jl")
include("objective.jl")
include("measure_expansions.jl")
include("derivative_evaluations.jl")
include("optimize.jl")
include("results.jl")
include("show.jl")
include("utilities.jl")
include("general_variables.jl")

# Import and export TranscriptionOpt
include("TranscriptionOpt/TranscriptionOpt.jl")
Reexport.@reexport using .TranscriptionOpt

# Deprecations introduced with v0.4.0
Base.@deprecate(ParameterBounds, DomainRestrictions)
Base.@deprecate(has_parameter_bounds, has_domain_restrictions)
Base.@deprecate(parameter_bounds, domain_restrictions)
Base.@deprecate(set_parameter_bounds, set_domain_restrictions)
Base.@deprecate(add_parameter_bounds, add_domain_restrictions)
Base.@deprecate(delete_parameter_bounds, delete_domain_restrictions)
Base.@deprecate(IntervalSet, IntervalDomain)
Base.@deprecate(UniDistributionSet, UniDistributionDomain)
Base.@deprecate(MultiDistributionSet, MultiDistributionDomain)
Base.@deprecate(CollectionSet, CollectionDomain)
Base.@deprecate(supports_in_set, supports_in_domain)
Base.@deprecate(collection_sets, collection_domains)
Base.@deprecate(infinite_set, infinite_domain)
Base.@deprecate(set_infinite_set, set_infinite_domain)
Base.@deprecate(InfiniteParameterFunction, ParameterFunction)

for op in (:has_parameter_bounds, :parameter_bounds, :set_parameter_bounds, 
           :add_parameter_bounds, :delete_parameter_bounds)
    @eval begin 
        function $op(::JuMP.AbstractVariableRef, args...; kwargs...)
            error("The use of parameter bounds (now called domain ",
                  "domain restrictions) has been discontinued for finite ",
                  "variables. The preferred syntax is to specify ",
                  "`DomainRestrictions` in connection with constraints.")
        end
    end
end

macro infinite_variable(model, args...)
    error("`@infinite_variable` has been deprecated in favor of `@variable`. ", 
          "\n\nOld Syntax: `@infinite_variable(model, var[idxs...](params...), ",
          "args..., kwargs...)`\nNew Syntax: `@variable(model, var[idxs...], ", 
          "Infinite(params...), args..., kwargs...)`.")
end
macro point_variable(model, args...)
    error("`@point_variable` has been deprecated in favor of `@variable`. ", 
          "\n\nOld Syntax: `@point_variable(model, ivref(param_values...), ",
          "var_expr, args..., kwargs...)`\nNew Syntax: `@variable(model, ",
          "var_expr, Point(ivref, param_values...), args..., kwargs...)`.")
end
macro hold_variable(model, args...)
    error("`@hold_variable` has been deprecated in favor of `@variable`.",
          "\n\nOld Syntax: `@hold_variable(model, var_expr, args..., ",
          "kwargs...)`\nNew Syntax: `@variable(model, var_expr, args..., ",
          "kwargs...)`.")
end
macro derivative_variable(model, args...)
    error("`@derivative_variable` has been deprecated in favor of `@variable`. ", 
          "\n\nOld Syntax: `@derivative_variable(model, d(ivref)/d(pref), ",
          "var_expr, kwargs...)`\nNew Syntax: `@variable(model, var_expr, ",
          "Deriv(ivref, pref), kwargs...)`.")
end
macro BDconstraint(model, args...)
    error("`@BDconstraint` has been deprecated in favor of `@constraint`. ", 
          "\n\nOld Syntax: `@BDconstraint(model, ref_expr(conditions...), ",
          "constr_expr, kwargs...)`\nNew Syntax: `@constraint(model, ref_expr, ",
          "constr_expr, DomainRestrictions(conditions...), kwargs...)`.")
end
macro set_parameter_bounds(args...)
    error("`@set_parameter_bounds` has been deprecated in favor of ",
          "`set_domain_restrictions`. Macro based modification is no longer ",
          "supported.")
end
macro add_parameter_bounds(args...)
    error("`@add_parameter_bounds` has been deprecated in favor of ",
          "`add_domain_restrictions`. Macro based modification is no longer ",
          "supported.")
end

# Define additional stuff that should not be exported
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

# Following JuMP, export everything that doesn't start with a _ 
for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_") || startswith(sym_string, "@_")
        continue
    end
    if !(Base.isidentifier(sym) || (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end])))
        continue
    end
    @eval export $sym
end

end # module
