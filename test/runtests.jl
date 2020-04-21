# Load in the dependencies
using InfiniteOpt, JuMP, MathOptInterface, Distributions, Random,
FastGaussQuadrature, DataStructures

# load the test module
using Test

# Define convenient aliases
const IC = InfiniteOpt.Collections
# const MEM = InfiniteOpt.MeasureEvalMethods
const IOTO = InfiniteOpt.TranscriptionOpt
const JuMPC = JuMP.Containers
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const MOIUC = MathOptInterface.Utilities.CleverDicts
const FGQ = FastGaussQuadrature

# Load in testing utilities
include("utilities.jl")

# Run unit tests
println("-----------------------------------------")
println("----------------UNIT TESTS---------------")
println("-----------------------------------------")
@time @testset "Collections" begin
    include("Collections/VectorTuple.jl")
    include("Collections/DualDict.jl")
end
println("")
@time @testset "Datatypes" begin include("datatypes.jl") end
println("")
@time @testset "Utilities" begin include("utility_tests.jl") end
println("")
@time @testset "Infinite Set Methods" begin include("infinite_sets.jl") end
println("")
@time @testset "Common General Variable Methods" begin
    include("standalone_general_variable_methods.jl")
end
println("")
@time @testset "Parameter Methods" begin
    @testset "Scalar" begin include("scalar_parameters.jl") end
    @testset "Array" begin include("array_parameters.jl") end
end
println("")
# @time @testset "Optimizer Setup Methods" begin include("optimizer_setup.jl") end
# println("")
# @time @testset "Variable Methods" begin
#     @testset "Basic " begin include("variable_basics.jl") end
#     @testset "Info" begin include("variable_info.jl") end
#     @testset "Definition" begin include("variable_definition.jl") end
#     @testset "Macros" begin include("variable_macros.jl") end
#     @testset "Query" begin include("variable_queries.jl") end
#     @testset "Reduction Variables" begin include("reduction_variables.jl") end
# end
# println("")
# @time @testset "General Variable Methods" begin
#     include("general_variables.jl")
# end
# println("")
# @time @testset "Operators" begin include("operators.jl") end
# println("")
# @time @testset "Expression Methods" begin include("expressions.jl") end
# println("")
# @time @testset "Macro Expressions" begin include("macro_expressions.jl") end
# println("")
# @time @testset "Measure Evaluation Methods" begin
#     include("MeasureEvalMethods/methods.jl")
# end
# println("")
# @time @testset "Measure Methods" begin include("measures.jl") end
# println("")
# @time @testset "Objective Methods" begin include("objective.jl") end
# println("")
# @time @testset "Constraint Methods" begin include("constraints.jl") end
# println("")
# @time @testset "Deletion Methods" begin include("deletion.jl") end
# println("")
# @time @testset "Expansion Methods" begin include("measure_expansions.jl") end
# println("")
# @time @testset "Printing Methods" begin include("show.jl") end
# println("")
# @time @testset "TranscriptionOpt" begin
#     @testset "Model" begin include("TranscriptionOpt/model.jl") end
#     @testset "Measures" begin include("TranscriptionOpt/measure.jl") end
#     @testset "Transcribe" begin include("TranscriptionOpt/transcribe.jl") end
#     @testset "Optimize" begin include("TranscriptionOpt/optimize.jl") end
# end
# println("")
# @time @testset "Solution Methods" begin include("optimizer.jl") end
# println("")
# @time @testset "Solution Queries" begin include("results.jl") end
# println("")
# @time @testset "Extensions" begin include("extensions.jl") end
# println("")

# TODO add involved deletion tests
println("-----------------------------------------")
println("-------------TESTING COMPLETE------------")
println("-----------------------------------------")
