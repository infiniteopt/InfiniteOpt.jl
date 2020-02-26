# Test build_optimizer_model!
@testset "build_optimizer_model!" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
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
    # test normal usage
    @test isa(build_optimizer_model!(m), Nothing)
    @test optimizer_model_ready(m)
    @test num_variables(optimizer_model(m)) == 8
end

# Test optimizer_model_variable
@testset "optimizer_model_variable" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par) >= 0)
    @point_variable(m, x(0), x0)
    @hold_variable(m, z >= 0)
    build_optimizer_model!(m)
    # test normal usage
    @test optimizer_model_variable(x) == transcription_variable(x)
    @test optimizer_model_variable(x0) == transcription_variable(x0)
    @test optimizer_model_variable(z) == transcription_variable(z)
    # test fallback
    @test_throws ErrorException optimizer_model_variable(x, Val(:Bad))
end

# Test optimizer_model_constraint
@testset "optimizer_model_constraint" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par))
    @point_variable(m, x(0), x0)
    @hold_variable(m, z)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - z, data1)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    build_optimizer_model!(m)
    # test normal usage
    @test optimizer_model_constraint(c1) == transcription_constraint(c1)
    @test optimizer_model_constraint(c2) == transcription_constraint(c2)
    @test optimizer_model_constraint(c3) == transcription_constraint(c3)
    # test fallback
    @test_throws ErrorException optimizer_model_constraint(c1, Val(:Bad))
end

# Test optimize!
@testset "JuMP.optimize!" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
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
    # test normal usage
    @test isa(optimize!(m), Nothing)
    @test optimizer_model_ready(m)
    @test num_variables(optimizer_model(m)) == 8
end

# Test JuMP.result_count
@testset "JuMP.result_count" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    tm = optimizer_model(m)
    model = MOIU.Model{Float64}()
    JuMP.optimize!(tm)
    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.ResultCount(), 2)
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    @test result_count(m) == 2
end
