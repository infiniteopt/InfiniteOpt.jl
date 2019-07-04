# include("C:/Users/puls446/.julia/dev/InfiniteOpt/test/runtests.jl")

using InfiniteOpt, JuMP
using MathOptInterface
const MOI = MathOptInterface
using Distributions
using Ipopt
using Test

# Run unit tests
println("-----------------------------------------")
println("----------------UNIT TESTS---------------")
println("-----------------------------------------")
@time @testset "Datatypes" begin include("datatypes.jl") end
println("")
@time @testset "Parameter Methods" begin include("parameters.jl") end
println("")
