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
    # bridge_constraints
    @testset "JuMP.bridge_constraints" begin
        @test !bridge_constraints(m)
        set_optimizer(optimizer_model(m), with_optimizer(Ipopt.Optimizer))
        @test bridge_constraints(m)
    end
    # add_bridge
    @testset "JuMP.add_bridge" begin
        struct test_bridge{C} <: MOI.Bridges.AbstractBridge where {C} end
        @test isa(add_bridge(m, test_bridge), Nothing)
    end
    # set_optimizer
    @testset "JuMP.set_optimizer" begin
        m = InfiniteModel()
        @test isa(set_optimizer(m, with_optimizer(Ipopt.Optimizer)), Nothing)
    end
    # TODO add silent methods with implemented by JuMP
end
