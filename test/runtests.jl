# include("C:/Users/pulsipher/.julia/dev/InfiniteOpt/test/runtests.jl")

using InfiniteOpt, JuMP, MathOptInterface, Distributions
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const JuMPC = JuMP.Containers
const IOTO = InfiniteOpt.TranscriptionOpt
using Test
using Random
const MEM = InfiniteOpt.MeasureEvalMethods
using FastGaussQuadrature
const FGQ = FastGaussQuadrature

include("utilities.jl")

# Run unit tests
println("-----------------------------------------")
println("----------------UNIT TESTS---------------")
println("-----------------------------------------")
@time @testset "Datatypes" begin include("datatypes.jl") end
println("")
@time @testset "Utilities" begin include("utility_tests.jl") end
println("")
@time @testset "Parameter Methods" begin include("parameters.jl") end
println("")
@time @testset "Optimizer Setup Methods" begin include("optimizer_setup.jl") end
println("")
@time @testset "Variable Methods" begin
    @testset "Basic " begin include("variable_basics.jl") end
    @testset "Info" begin include("variable_info.jl") end
    @testset "Definition" begin include("variable_definition.jl") end
    @testset "Macros" begin include("variable_macros.jl") end
    @testset "Query" begin include("variable_queries.jl") end
    @testset "Reduction Variables" begin include("reduction_variables.jl") end
end
println("")
@time @testset "Operators" begin include("operators.jl") end
println("")
@time @testset "Expression Methods" begin include("expressions.jl") end
println("")
@time @testset "Measure Evaluation Methods" begin
    include("MeasureEvalMethods/methods.jl")
end
println("")
@time @testset "Measure Methods" begin include("measures.jl") end
println("")
@time @testset "Objective Methods" begin include("objective.jl") end
println("")
@time @testset "Constraint Methods" begin include("constraints.jl") end
println("")
@time @testset "Deletion Methods" begin include("deletion.jl") end
println("")
@time @testset "Expansion Methods" begin include("measure_expansions.jl") end
println("")
@time @testset "Printing Methods" begin include("show.jl") end
println("")
@time @testset "TranscriptionOpt" begin
    @testset "Model" begin include("TranscriptionOpt/model.jl") end
    @testset "Transcribe" begin include("TranscriptionOpt/transcribe.jl") end
    @testset "Optimize" begin include("TranscriptionOpt/optimize.jl") end
end
println("")
@time @testset "Solution Methods" begin include("optimizer.jl") end
println("")
@time @testset "Solution Queries" begin include("results.jl") end
println("")

# TODO add extension tests
# TODO add involved deletion tests
println("-----------------------------------------")
println("-------------TESTING COMPLETE------------")
println("-----------------------------------------")
