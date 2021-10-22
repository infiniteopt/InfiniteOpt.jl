# Test transcription_model
@testset "transcription_model" begin
    # initialize model
    m = InfiniteModel()
    # test normal usage
    @test isa(transcription_model(m), Model)
    # test error
    set_optimizer_model(m, Model())
    @test_throws ErrorException transcription_model(m)
end

# Test build_optimizer_model!
@testset "build_optimizer_model!" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    @infinite_parameter(m, par in [0, 1], num_supports = 2)
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @variable(m, 1 >= x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = 0)
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    meas1 = support_sum(x - w, par)
    meas2 = integral(y, pars)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, DomainRestrictions(par => [0.5, 1]))
    @constraint(m, c5, meas2 == 0)
    @constraint(m, @deriv(x, par) == 0)
    @constraint(m, sin(w) + integral(x^3, par) == 0)
    @objective(m, Min, x0 + meas1)
    set_silent(m)
    set_time_limit_sec(m, 42)
    # test normal usage
    @test isa(build_optimizer_model!(m, Val(:TransData)), Nothing)
    @test optimizer_model_ready(m)
    @test num_variables(optimizer_model(m)) == 30
    @test time_limit_sec(optimizer_model(m)) == 42
end
