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
    @testset "Single Argument Methods" begin
        for f in (set_silent, unset_silent, bridge_constraints, 
                  time_limit_sec, unset_time_limit_sec, solver_name, backend,
                  JuMP.mode, unsafe_backend, compute_conflict!, copy_conflict,
                  set_string_names_on_creation)
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
    @testset "JuMP.set_time_limit_sec" begin
        @test_throws ErrorException set_time_limit_sec(TestBackend(), 42)
        @test set_time_limit_sec(jump_backend, 42) isa Nothing
        @test time_limit_sec(jump_backend) == 42
        @test set_time_limit_sec(m, 42) isa Nothing
        @test time_limit_sec(m) == 42
    end
    @testset "JuMP.set_string_names_on_creation" begin
        @test_throws ErrorException set_string_names_on_creation(TestBackend(), false)
        @test set_string_names_on_creation(jump_backend, false) isa Nothing
        @test set_string_names_on_creation(jump_backend) == false
        @test set_string_names_on_creation(m, true) isa Nothing
        @test set_string_names_on_creation(m) == true
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
        @test solver_name(jump_backend2) == "Mock"
        m2 = InfiniteModel()
        @test set_optimizer(m2, mockoptimizer) isa Nothing
        @test solver_name(m) == "Mock"
    end
end
