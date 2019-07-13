# include("C:/Users/puls446/.julia/dev/InfiniteOpt/test/runtests.jl")

using InfiniteOpt, JuMP
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
using Distributions
# using Ipopt
using Test

include("Utilities.jl")

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
@time @testset "Basic Variable Methods" begin include("variable_basics.jl") end
println("")
@time @testset "Variable Info Methods" begin include("variable_info.jl") end
println("")
@time @testset "Variable Definition Methods" begin include("variable_definition.jl") end
println("")
@time @testset "Variable Macros" begin include("variable_macros.jl") end
println("")
@time @testset "Variable Query Methods" begin include("variable_queries.jl") end
println("")

# TODO add tests for deleting
