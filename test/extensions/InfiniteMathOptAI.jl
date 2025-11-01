using MathOptAI, Ipopt

@testset "InfiniteMathOptAI.jl" begin
    predictor = MathOptAI.Tanh()
    model = InfiniteModel(Ipopt.Optimizer)
    set_silent(model)
    @infinite_parameter(model, ξ ~ Uniform(0, 1))
    @variable(model, -1 <= x <= 3, Infinite(ξ))
    @variable(model, z)
    y, _ = add_predictor(model, predictor, [x])
    w, _ = add_predictor(model, predictor, [z])
    @test isempty(parameter_refs(only(w)))
    @test parameter_refs(only(y)) == (ξ,)
    @objective(model, Max, 𝔼(only(y), ξ))
    @constraint(model, x <= ξ)
    @constraint(model, only(w) <= 2)
    optimize!(model)
    # @test get_variable_bounds(y[1]) == (tanh(-1), tanh(3))
    y_v = value(only(y))
    @test isapprox(y_v, tanh.(value(x)); atol = 1e-5)
    @test isapprox(objective_value(model), sum(y_v) / length(y_v); atol = 1e-5)
end