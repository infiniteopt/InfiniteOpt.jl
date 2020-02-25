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
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
     # set_optimizer
     @testset "JuMP.set_optimizer" begin
         m2 = InfiniteModel()
         @test isa(set_optimizer(m2, mockoptimizer), Nothing)
         @test m2.optimizer_constructor == mockoptimizer
     end
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
    # set_silent
    @testset "JuMP.set_silent" begin
        @test set_silent(m)
    end
    # unset_silent
    @testset "JuMP.unset_silent" begin
        @test !unset_silent(m)
    end
    # set_time_limit_sec
    @testset "JuMP.set_time_limit_sec" begin
        @test set_time_limit_sec(m, 100) == 100
    end
    # unset_time_limit_sec
    @testset "JuMP.unset_time_limit_sec" begin
        @test isa(unset_time_limit_sec(m), Nothing)
    end
    # time_limit_sec
    @testset "JuMP.time_limit_sec" begin
        @test time_limit_sec(m) == nothing
    end
    # set_optimizer_attribute
    @testset "JuMP.set_optimizer_attribute (String)" begin
        @test set_optimizer_attribute(m, "mine", 42) == 42
    end
    # set_optimizer_attribute
    @testset "JuMP.set_optimizer_attribute (MOI)" begin
        @test set_optimizer_attribute(m, MOI.Silent(), true)
    end
    # set_optimizer_attributes
    @testset "JuMP.set_optimizer_attributes" begin
        @test isa(set_optimizer_attributes(m, MOI.Silent() => false, "mine" => 1), Nothing)
        @test !MOI.get(optimizer_model(m), MOI.Silent())
        @test MOI.get(optimizer_model(m), MOI.RawParameter("mine")) == 1
    end
    # get_optimizer_attribute
    @testset "JuMP.get_optimizer_attribute (String)" begin
        @test get_optimizer_attribute(m, "mine") == 1
    end
    # get_optimizer_attribute
    @testset "JuMP.get_optimizer_attribute (MOI)" begin
        @test !get_optimizer_attribute(m, MOI.Silent())
    end
    # solver_name
    @testset "JuMP.solver_name" begin
        @test solver_name(m) == "Mock"
    end
    # backend
    @testset "JuMP.backend" begin
        @test backend(m) == backend(optimizer_model(m))
    end
    # mode
    @testset "JuMP.mode" begin
        @test JuMP.mode(m) == JuMP.mode(optimizer_model(m))
    end
    # solve
    @testset "solve" begin
        @test_throws ErrorException solve(m)
    end
end
