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
using .Collections: VectorTuple, same_structure

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
