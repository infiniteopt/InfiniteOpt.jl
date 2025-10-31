module InfiniteOpt

# Import and export JuMP 
import Reexport 
Reexport.@reexport using JuMP

# Import the necessary packages.
import Distributions
import DataStructures
import FastGaussQuadrature
import LinearAlgebra
import MutableArithmetics
using Base.Meta

# Make useful aliases (note we get MOI and MOIU from JuMP)
const JuMPC = JuMP.Containers
const MOIUC = MOIU.CleverDicts
const _MA = MutableArithmetics
export JuMPC, MOIUC # this makes these accessible to the submodules

# Import the Collections module
include("Collections/Collections.jl")

# Import core methods
include("datatypes.jl")
include("infinite_domains.jl")
include("scalar_parameters.jl")
include("array_parameters.jl")
include("variable_basics.jl")
include("infinite_variables.jl")
include("semi_infinite_variables.jl")
include("point_variables.jl")
include("finite_variables.jl")
include("nlp.jl")
include("macros.jl")
include("parameter_functions.jl")
include("expressions.jl")
include("measures.jl")

# Import and export MeasureToolbox
include("MeasureToolbox/MeasureToolbox.jl")
Reexport.@reexport using .MeasureToolbox

# import more core methods
include("derivatives.jl")
include("constraints.jl")
include("objective.jl")
include("measure_expansions.jl")
include("derivative_evaluations.jl")
include("backends.jl")
include("results.jl")
include("show.jl")
include("utilities.jl")
include("general_variables.jl")

# Import and export TranscriptionOpt
include("TranscriptionOpt/TranscriptionOpt.jl")
Reexport.@reexport using .TranscriptionOpt

# Add deprecation and old syntax methods 
macro register(args...)
    error("`@register` has now been replaced with `@operator`, see ",
           "the nonlinear documenation page for details.")
end
Base.@deprecate map_nlp_to_ast(f, expr) map_expression_to_ast(f, expr)
Base.@deprecate optimizer_model_variable(v; kwargs...) transformation_variable(v; kwargs...)
Base.@deprecate optimizer_model_expression(e; kwargs...) transformation_expression(e; kwargs...)
Base.@deprecate optimizer_model_constraint(c; kwargs...) transformation_constraint(c; kwargs...)
Base.@deprecate optimizer_model(m) transformation_model(m)
Base.@deprecate set_value(v::GeneralVariableRef, val) set_parameter_value(v, val)
Base.@deprecate has_domain_restrictions(cref) has_domain_restriction(cref)
Base.@deprecate domain_restrictions(cref) domain_restriction(cref)
Base.@deprecate set_domain_restrictions(cref, restrict) set_domain_restriction(cref, restrict)
Base.@deprecate delete_domain_restrictions(cref) delete_domain_restriction(cref)
function DomainRestrictions(args...)
    error("`DomainRestrictions` is deprecated in favor of `DomainRestriction`, see docs for details.")
end
Base.@deprecate start_value_function(vref)  start_value(vref)
Base.@deprecate set_start_value_function(vref, f)  set_start_value(vref, f)

# Define additional stuff that should not be exported
const _EXCLUDE_SYMBOLS = [
    Symbol(@__MODULE__), 
    :eval, 
    :include,
    :derivative_expr_data,
    :make_indexed_derivative_expr,
    :allows_high_order_derivatives
    ]

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
