# Test the optimizer model methods
@testset "Optimizer Model" begin
    m = InfiniteModel()
    # optimizer_model
    @testset "optimizer_model" begin
        @test isa(optimizer_model(m), Model)
    end
    # optimizer_model_ready
    @testset "optimizer_model_ready" begin
        @test !optimizer_model_ready(m)
        m.ready_to_optimize = true
        @test optimizer_model_ready(m)
    end
    # set_optimizer_model_ready
    @testset "set_optimizer_model_ready" begin
        @test isa(set_optimizer_model_ready(m, false), Nothing)
        @test !optimizer_model_ready(m)
    end
    # set_optimizer_model
    @testset "set_optimizer_model" begin
        @test isa(set_optimizer_model(m, Model()), Nothing)
        @test length(optimizer_model(m).ext) == 0
    end
    # optimizer_model_key
    @testset "optimizer_model_key" begin
        m = InfiniteModel()
        @test optimizer_model_key(m) == :TransData
        optimizer_model(m).ext[:extra] = 42
        @test_throws ErrorException optimizer_model_key(m)
    end
end

# Test JuMP extensions
@testset "JuMP Extensions" begin
    m = InfiniteModel()
    mockoptimizer = with_optimizer(MOIU.MockOptimizer,
                                   MOIU.Model{Float64}(),
                                   eval_objective_value=false)
    # bridge_constraints
    @testset "JuMP.bridge_constraints" begin
        @test !bridge_constraints(m)
        set_optimizer(optimizer_model(m), mockoptimizer)
        @test bridge_constraints(m)
    end
    # add_bridge
    @testset "JuMP.add_bridge" begin
        # @test isa(add_bridge(m, TestBridge), Nothing)
        @test isa(add_bridge(m, MOI.Bridges.Variable.VectorizeBridge), Nothing)
    end
    # set_optimizer
    @testset "JuMP.set_optimizer" begin
        m2 = InfiniteModel()
        @test isa(set_optimizer(m2, mockoptimizer), Nothing)
    end
    # set_silent
    @testset "JuMP.set_silent" begin
        m2 = InfiniteModel()
        @test isa(set_silent(m2), Nothing)
    end
    # unset_silent
    @testset "JuMP.unset_silent" begin
        m2 = InfiniteModel()
        @test isa(unset_silent(m2), Nothing)
    end
    # solve
    @testset "solve" begin
        @test_throws ErrorException solve(m)
    end
end
