# import Pkg
# Pkg.pkg"add JuMP#od/nlp-expr"

using InfiniteOpt: _domain_or_error
using Test: Error
# Load in the dependencies
using InfiniteOpt, Distributions, Random, FastGaussQuadrature, DataStructures, 
Suppressor, LinearAlgebra
import MutableArithmetics

# load the test module
using Test

# Define convenient aliases
const IC = InfiniteOpt.Collections
const MT = InfiniteOpt.MeasureToolbox
const IOTO = InfiniteOpt.TranscriptionOpt
const JuMPC = JuMP.Containers
const MOIUC = MOIU.CleverDicts
const FGQ = FastGaussQuadrature
const IOMT = InfiniteOpt.MeasureToolbox
const MA = MutableArithmetics

# Load in testing utilities
include("utilities.jl")

# Run unit tests
println("----------------------------------------------------------------------------")
println("---------------------------------UNIT TESTS---------------------------------")
println("----------------------------------------------------------------------------")
@time @testset "Collections" begin
    include("Collections/vectorize.jl")
    include("Collections/VectorTuple.jl")
end
println("")
@time @testset "Datatypes" begin include("datatypes.jl") end
println("")
@time @testset "Utilities" begin include("utility_tests.jl") end
println("")
@time @testset "Infinite Domain Methods" begin include("infinite_domains.jl") end
println("")
@time @testset "General Variable Methods" begin
    include("general_variables.jl")
end
println("")
@time @testset "Optimizer Setup Methods" begin include("optimizer_setup.jl") end
println("")
@time @testset "Macro Utilities" begin include("macro_utilities.jl") end
println("")
@time @testset "Parameter Methods" begin
   @testset "Scalar" begin include("scalar_parameters.jl") end
   @testset "Array" begin include("array_parameters.jl") end
end
println("")
@time @testset "Variable Methods" begin
   @testset "Infinite Variables" begin include("infinite_variables.jl") end
   @testset "Point Variables" begin include("point_variables.jl") end
   @testset "Finite Variables" begin include("finite_variables.jl") end
   @testset "Info Constraints" begin include("variable_info.jl") end
   @testset "Semi-Infinite Variables" begin include("semi_infinite_variables.jl") end
end
println("")
@time @testset "Derivative Methods" begin include("derivatives.jl") end
println("")
@time @testset "Nonlinear" begin include("nlp.jl") end
println("")
@time @testset "Operators" begin include("operators.jl") end
println("")
@time @testset "Expression Methods" begin include("expressions.jl") end
println("")
@time @testset "Macro Expressions" begin include("macro_expressions.jl") end
println("")
@time @testset "Measure Methods" begin include("measures.jl") end
println("")
@time @testset "Measure Toolbox Methods" begin
    @testset "Integrals" begin include("MeasureToolbox/integrals.jl") end
    @testset "Expectations" begin include("MeasureToolbox/expectations.jl") end
    @testset "Support Sums" begin include("MeasureToolbox/support_sums.jl") end
end
println("")
@time @testset "Objective Methods" begin include("objective.jl") end
println("")
@time @testset "Constraint Methods" begin include("constraints.jl") end
println("")
@time @testset "Printing Methods" begin include("show.jl") end
println("")
@time @testset "Deletion Methods" begin include("deletion.jl") end
println("")
@time @testset "Expansion Methods" begin include("measure_expansions.jl") end
println("")
@time @testset "Derivative Evaluation" begin include("derivative_evaluation.jl") end
println("")
@time @testset "TranscriptionOpt" begin
    @testset "Model" begin include("TranscriptionOpt/model.jl") end
    @testset "Measures" begin include("TranscriptionOpt/measure.jl") end
    @testset "Transcribe" begin include("TranscriptionOpt/transcribe.jl") end
    @testset "Optimize" begin include("TranscriptionOpt/optimize.jl") end
end
println("")
@time @testset "Solution Methods" begin include("optimizer.jl") end
println("")
@time @testset "Solution Queries" begin include("results.jl") end
println("")
@time @testset "Extensions" begin include("extensions.jl") end
println("")
println("----------------------------------------------------------------------------")
println("-----------------------------TESTING COMPLETE!------------------------------")
println("----------------------------------------------------------------------------")
