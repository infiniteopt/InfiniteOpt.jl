# Test the backend methods
@testset "Transformation Backends" begin
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    bmodel = Model(mockoptimizer)
    jump_backend = JuMPBackend{TestJuMPTag}(bmodel, 42)
    # transformation_backend_ready
    @testset "transformation_backend_ready" begin
        @test !transformation_backend_ready(m)
        m.ready_to_optimize = true
        @test transformation_backend_ready(m)
    end
    # set_transformation_backend_ready
    @testset "set_transformation_backend_ready" begin
        @test isa(set_transformation_backend_ready(m, false), Nothing)
        @test !transformation_backend_ready(m)
    end
    # transformation_model
    @testset "transformation_model" begin
        @test isa(transformation_model(m), Model)
        @test transformation_model(jump_backend) == bmodel
        @test_throws ErrorException transformation_model(TestBackend())
        @test (@test_deprecated optimizer_model(m)) == transformation_model(m)
    end
    # transformation_data
    @testset "transformation_data" begin
        @test isa(transformation_data(m), IOTO.TranscriptionData)
        @test transformation_data(jump_backend) == 42
        @test_throws ErrorException transformation_data(TestBackend())
    end
    # transformation_backend
    @testset "transformation_backend" begin
        @test transformation_backend(m) isa TranscriptionBackend
    end
    # set_transformation_backend
    @testset "set_transformation_backend" begin
        current_backend = m.backend
        @test set_transformation_backend(m, jump_backend) isa Nothing
        @test transformation_model(m) == bmodel
        @test set_transformation_backend(m, current_backend) isa Nothing
        @test transformation_data(m) isa IOTO.TranscriptionData
    end
    # JuMP.get_attribute and JuMP.set_attribute
    @testset "JuMP.[get/set]_attribute" begin
        @test set_attribute(m, MOI.TimeLimitSec(), 10.) isa Nothing
        @test get_attribute(m, MOI.TimeLimitSec()) == 10
        @test set_attribute(jump_backend, MOI.TimeLimitSec(), 10.) isa Nothing
        @test get_attribute(jump_backend, MOI.TimeLimitSec()) == 10
        @test_throws ErrorException get_attribute(TestBackend(), MOI.TimeLimitSec())
        @test_throws ErrorException set_attribute(TestBackend(), MOI.TimeLimitSec(), 10.)
        @test set_optimizer_attribute(m, MOI.TimeLimitSec(), 12.) isa Nothing
        @test get_optimizer_attribute(m, MOI.TimeLimitSec()) == 12
        @test set_attributes(m, MOI.TimeLimitSec() => 10.) isa Nothing
        @test get_attribute(m, MOI.TimeLimitSec()) == 10
        @test set_attribute(m, "something", 42) isa Nothing
        @test get_attribute(m, "something") == 42
    end
    # Base.empty!
    @testset "Base.empty!" begin
        @test empty!(JuMPBackend{TestJuMPTag}(Model(), [42])).data == []
        @test_throws ErrorException empty!(TestBackend())
    end
end

# Test JuMP extensions
@testset "JuMP Extensions" begin
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    bmodel = Model(mockoptimizer)
    jump_backend = JuMPBackend{TestJuMPTag}(bmodel, 42)
    set_time_limit_sec(bmodel, 10)
    @testset "Single Argument Direct Methods" begin
        for f in (bridge_constraints, backend, JuMP.mode, unsafe_backend,
                  compute_conflict!, copy_conflict)
            @test_throws ErrorException f(TestBackend())
            if f != copy_conflict
                @test f(jump_backend) == f(bmodel)
                @test f(m) == f(m.backend)
            else
                @test f(jump_backend) isa Tuple
                @test f(m) isa Tuple
            end
        end
    end
    @testset "JuMP.set_silent" begin
        @test set_silent(m) isa Nothing
        @test get_attribute(m, MOI.Silent())
        @test set_silent(InfiniteModel(jump_backend)) isa Nothing
        @test get_attribute(jump_backend, MOI.Silent())
        @test_throws ErrorException set_silent(InfiniteModel(TestBackend()))
    end
    @testset "JuMP.unset_silent" begin
        @test unset_silent(m) isa Nothing
        @test !get_attribute(m, MOI.Silent())
        @test unset_silent(InfiniteModel(jump_backend)) isa Nothing
        @test !get_attribute(jump_backend, MOI.Silent())
        @test_throws ErrorException unset_silent(InfiniteModel(TestBackend()))
    end
    @testset "JuMP.time_limit_sec" begin
        @test time_limit_sec(m) isa Nothing
        @test time_limit_sec(InfiniteModel(jump_backend)) == 10
        @test_throws ErrorException time_limit_sec(InfiniteModel(TestBackend()))
    end
    @testset "JuMP.set_time_limit_sec" begin
        @test_throws ErrorException set_time_limit_sec(InfiniteModel(TestBackend()), 42)
        @test set_time_limit_sec(InfiniteModel(jump_backend), 42) isa Nothing
        @test time_limit_sec(bmodel) == 42
        @test set_time_limit_sec(m, 42) isa Nothing
        @test time_limit_sec(m) == 42
    end
    @testset "JuMP.unset_time_limit_sec" begin
        @test unset_time_limit_sec(m) isa Nothing
        @test time_limit_sec(m) isa Nothing
        @test unset_time_limit_sec(InfiniteModel(jump_backend)) isa Nothing
        @test time_limit_sec(bmodel) isa Nothing
        @test_throws ErrorException unset_time_limit_sec(InfiniteModel(TestBackend()))
    end
    @testset "JuMP.solver_name" begin
        @test solver_name(m) == get_attribute(m.backend, MOI.SolverName())
        @test solver_name(InfiniteModel(jump_backend)) == solver_name(bmodel)
        @test_throws ErrorException solver_name(InfiniteModel(TestBackend()))
    end
    @testset "JuMP.add_bridge" begin
        bridge = MOI.Bridges.Variable.VectorizeBridge
        @test_throws ErrorException add_bridge(TestBackend(), bridge)
        @test add_bridge(jump_backend, bridge) isa Nothing
        @test add_bridge(m, bridge) isa Nothing
    end
    @testset "JuMP.print_active_bridges" begin
        @test_throws ErrorException sprint(print_active_bridges, TestBackend())
        expected = " * Supported objective: MOI.ScalarAffineFunction{Float64}\n"
        @test sprint(print_active_bridges, jump_backend) == expected
        @test sprint(print_active_bridges, m) == expected
        stdout_test(print_active_bridges, expected, m)
    end
    @testset "print_bridge_graph" begin
        @test_throws ErrorException sprint(print_bridge_graph, TestBackend())
        expected = "Bridge graph with 0 variable nodes, 0 constraint nodes and 0 objective nodes.\n"
        @test sprint(print_bridge_graph, jump_backend) == expected
        @test sprint(print_bridge_graph, m) == expected
        stdout_test(print_bridge_graph, expected, m)
    end
    @testset "JuMP.set_optimizer" begin
        @test_throws ErrorException set_optimizer(TestBackend(), mockoptimizer)
        bmodel2 = Model()
        jump_backend2 = JuMPBackend{TestJuMPTag}(bmodel2, 42)
        @test set_optimizer(jump_backend2, mockoptimizer) isa Nothing
        @test solver_name(jump_backend2.model) == "Mock"
        m2 = InfiniteModel()
        @test set_optimizer(m2, mockoptimizer) isa Nothing
        @test solver_name(m) == "Mock"
    end
end
