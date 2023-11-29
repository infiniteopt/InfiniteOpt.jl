# Test add_infinite_model_optimizer
@testset "add_infinite_model_optimizer" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    set_silent(m)
    set_time_limit_sec(m, 42.)
    # test normal
    tm = Model()
    @test InfiniteOpt.add_infinite_model_optimizer(tm, m) isa Nothing
    @test time_limit_sec(tm) == 42
    @test get_optimizer_attribute(tm, MOI.Silent())
end

# Test the optimizer model methods
@testset "Optimizer Model" begin
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
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
        # test with inheritance
        set_time_limit_sec(optimizer_model(m), 42.)
        @test isa(set_optimizer_model(m, Model()), Nothing)
        @test length(optimizer_model(m).ext) == 0
        @test time_limit_sec(optimizer_model(m)) == 42
        @test optimize!(optimizer_model(m)) isa Nothing
        # test without inheritance
        @test isa(set_optimizer_model(m, Model(), inherit_optimizer = false), Nothing)
        @test length(optimizer_model(m).ext) == 0
        @test_throws(
            Union{ErrorException,MOI.GetAttributeNotAllowed{MOI.TimeLimitSec}},
            time_limit_sec(optimizer_model(m)),
        )
        @test_throws NoOptimizer optimize!(optimizer_model(m))
    end
    # optimizer_model_key (optimizer models)
    @testset "optimizer_model_key (Model)" begin
        m = InfiniteModel()
        @test optimizer_model_key(optimizer_model(m)) == :TransData
        optimizer_model(m).ext[:extra] = 42
        @test_throws ErrorException optimizer_model_key(optimizer_model(m))
    end
    # optimizer_model_key (InfiniteModel)
    @testset "optimizer_model_key (InfiniteModel)" begin
        m = InfiniteModel()
        @test optimizer_model_key(m) == :TransData
        optimizer_model(m).ext[:extra] = 42
        @test_throws ErrorException optimizer_model_key(m)
    end
    # clear_optimizer_model_build! (optimizer models)
    @testset "clear_optimizer_model_build! (Model)" begin
        # setup
        m = TranscriptionModel(mockoptimizer)
        set_time_limit_sec(m, 42.)
        @variable(m, t)
        # test
        @test clear_optimizer_model_build!(m) isa Model
        @test num_variables(m) == 0
        @test length(m.ext) == 1
        @test time_limit_sec(m) == 42
        # test add variable again
        @test @variable(m, t) isa VariableRef
    end
    # clear_optimizer_model_build! (InfiniteModel)
    @testset "clear_optimizer_model_build! (InfiniteModel)" begin
        # setup
        m = InfiniteModel(mockoptimizer)
        set_time_limit_sec(optimizer_model(m), 42.)
        @variable(optimizer_model(m), t)
        # test
        @test clear_optimizer_model_build!(m) isa Model
        @test num_variables(optimizer_model(m)) == 0
        @test length(optimizer_model(m).ext) == 1
        @test time_limit_sec(optimizer_model(m)) == 42
        # test add variable again
        @test @variable(optimizer_model(m), t) isa VariableRef
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
        @test set_silent(m) isa Nothing
    end
    # unset_silent
    @testset "JuMP.unset_silent" begin
        @test unset_silent(m) isa Nothing
    end
    # set_time_limit_sec
    @testset "JuMP.set_time_limit_sec" begin
        @test set_time_limit_sec(m, 100.) isa Nothing
        @test time_limit_sec(m) == 100
    end
    # unset_time_limit_sec
    @testset "JuMP.unset_time_limit_sec" begin
        @test isa(unset_time_limit_sec(m), Nothing)
    end
    # time_limit_sec
    @testset "JuMP.time_limit_sec" begin
        @test time_limit_sec(m) === nothing
    end
    # set_optimizer_attribute
    @testset "JuMP.set_optimizer_attribute (String)" begin
        @test set_optimizer_attribute(m, "mine", 42.) isa Nothing
    end
    # set_optimizer_attribute
    @testset "JuMP.set_optimizer_attribute (MOI)" begin
        @test set_optimizer_attribute(m, MOI.Silent(), true) isa Nothing
    end
    # set_optimizer_attributes
    @testset "JuMP.set_optimizer_attributes" begin
        @test isa(set_optimizer_attributes(m, MOI.Silent() => false, "mine" => 1), Nothing)
        @test !MOI.get(optimizer_model(m), MOI.Silent())
        @test MOI.get(optimizer_model(m), MOI.RawOptimizerAttribute("mine")) == 1
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
end
