# Test general queries
@testset "General Queries" begin
    # Setup the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @objective(tm, Min, x^2)
    @constraint(tm, c1, 2 * x * y <= 1)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [x, y]
    modelstring = """
    variables: x, y
    minobjective: 1*x*x
    c1: 2*x*y <= 1.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x","y"], ["c1"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))

    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -1.0)
    MOI.set(mockoptimizer, MOI.SolveTime(), 0.42)
    # test termination_status
    @testset "JuMP.termination_status" begin
        @test termination_status(m) == MOI.OPTIMAL
    end
    # test primal_status
    @testset "JuMP.primal_status" begin
        @test primal_status(m) == MOI.FEASIBLE_POINT
    end
    # test dual_status
    @testset "JuMP.dul_status" begin
        @test dual_status(m) == MOI.FEASIBLE_POINT
    end
    # test solve_time
    @testset "JuMP.solve_time" begin
        @test solve_time(m) == 0.42
    end
end

# Test objective queries
@testset "Objective Queries" begin
    # Setup the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @objective(tm, Min, x^2)
    @constraint(tm, c1, 2 * x * y <= 1)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [x, y]
    modelstring = """
    variables: x, y
    minobjective: 1*x*x
    c1: 2*x*y <= 1.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x","y"], ["c1"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))

    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 2.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -1.0)
    MOI.set(mockoptimizer, MOI.SolveTime(), 0.42)
    # test objective_bound
    @testset "JuMP.objective_bound" begin
        @test objective_bound(m) == 2.0
    end
    # test objective_bound
    @testset "JuMP.objective_value" begin
        @test objective_value(m) == -1.
    end
end

# Test variable queries
@testset "Variable Queries" begin
    # Setup the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @objective(tm, Min, x^2)
    @constraint(tm, c1, 2 * x * y <= 1)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [x, y]
    modelstring = """
    variables: x, y
    minobjective: 1*x*x
    c1: 2*x*y <= 1.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x","y"], ["c1"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))

    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -1.0)
    MOI.set(mockoptimizer, MOI.SolveTime(), 0.42)
    # test has_values
    @testset "JuMP.has_values" begin
        @test has_values(m)
    end
    # test value
    @testset "JuMP.value" begin
        @test value(inf) == [1., 0.]
        @test value(g) == 1.
    end
end
