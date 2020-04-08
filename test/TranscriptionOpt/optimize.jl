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
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1])
    @infinite_variable(m, 1 >= x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @hold_variable(m, 0 <= z <= 1, Bin)
    @hold_variable(m, w == 1, Int, start = 1)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @BDconstraint(m, c4(par in [0.5, 1]), meas2 - 2y0 + x <= 1)
    @constraint(m, c5, meas2 == 0)
    @objective(m, Min, x0 + meas1)
    set_silent(m)
    set_time_limit_sec(m, 42)
    # test normal usage
    @test isa(build_optimizer_model!(m, Val(:TransData)), Nothing)
    @test optimizer_model_ready(m)
    @test num_variables(optimizer_model(m)) == 8
    @test time_limit_sec(optimizer_model(m)) == 42
end
