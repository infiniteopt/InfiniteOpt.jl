# Test map_value
@testset "map_value" begin
    # initialize the infinite model
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
    # test FiniteVariableRef
    @testset "FiniteVariableRef" begin
        @test map_value(g, Val(:TransData)) == 1.
    end
    # test InfiniteVariableRef
    @testset "InfiniteVariableRef" begin
        @test map_value(inf, Val(:TransData)) == [1., 0.]
    end
end
