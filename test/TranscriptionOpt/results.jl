# Test map_value
@testset "map_value" begin
    # initialize the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    @constraint(m, c1, g <= 0)
    @constraint(m, c2, inf >= 0)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @variable(tm, z)
    @objective(tm, Min, x^2)
    @constraint(tm, c3, x <= 0)
    @constraint(tm, c4, y >= 0)
    @constraint(tm, c5, z >= 0)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [y, z]
    tm.ext[:TransData].finite_to_constr[c1] = c3
    tm.ext[:TransData].infinite_to_constrs[c2] = [c4, c5]
    modelstring = """
    variables: x, y, z
    minobjective: 1*x*x
    c3: x <= 0.0
    c4: y >= 0.0
    c5: z >= 0.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))

    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(x), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(y), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(z), 1.0)
    # test FiniteVariableRef
    @testset "FiniteVariableRef" begin
        @test map_value(g, Val(:TransData)) == -1.
    end
    # test InfiniteVariableRef
    @testset "InfiniteVariableRef" begin
        @test map_value(inf, Val(:TransData)) == [0., 1.]
    end
    # test FiniteConstraintRef
    @testset "FiniteConstraintRef" begin
        @test map_value(c1, Val(:TransData)) == -1.
    end
    # test InfiniteConstraintRef
    @testset "InfiniteConstraintRef" begin
        @test map_value(c2, Val(:TransData)) == [0., 1.]
    end
end

# Test map_optimizer_index
@testset "map_optimizer_index" begin
    # initialize the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    @constraint(m, c1, g <= 0)
    @constraint(m, c2, inf >= 0)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @variable(tm, z)
    @objective(tm, Min, x^2)
    @constraint(tm, c3, x <= 0)
    @constraint(tm, c4, y >= 0)
    @constraint(tm, c5, z >= 0)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [y, z]
    tm.ext[:TransData].finite_to_constr[c1] = c3
    tm.ext[:TransData].infinite_to_constrs[c2] = [c4, c5]
    modelstring = """
    variables: x, y, z
    minobjective: 1*x*x
    c3: x <= 0.0
    c4: y >= 0.0
    c5: z >= 0.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    # MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x", "y", "z"],
    #                        ["c3", "c4", "c5"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))
    # test FiniteVariableRef
    @testset "FiniteVariableRef" begin
        @test isa(map_optimizer_index(g, Val(:TransData)), MOI.VariableIndex)
    end
    # test InfiniteVariableRef
    @testset "InfiniteVariableRef" begin
        @test isa(map_optimizer_index(inf, Val(:TransData)),
                  Vector{MOI.VariableIndex})
    end
    # test FiniteConstraintRef
    @testset "FiniteConstraintRef" begin
        @test isa(map_optimizer_index(c1, Val(:TransData)), MOI.ConstraintIndex)
    end
    # test InfiniteConstraintRef
    @testset "InfiniteConstraintRef" begin
        @test isa(map_optimizer_index(c2, Val(:TransData)),
                  Vector{<:MOI.ConstraintIndex})
    end
end

# Test map_dual
@testset "map_dual" begin
    # initialize the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    @constraint(m, c1, g <= 0)
    @constraint(m, c2, inf >= 0)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @variable(tm, z)
    @objective(tm, Min, x^2)
    @constraint(tm, c3, x <= 0)
    @constraint(tm, c4, y >= 0)
    @constraint(tm, c5, z >= 0)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [y, z]
    tm.ext[:TransData].finite_to_constr[c1] = c3
    tm.ext[:TransData].infinite_to_constrs[c2] = [c4, c5]
    modelstring = """
    variables: x, y, z
    minobjective: 1*x*x
    c3: x <= 0.0
    c4: y >= 0.0
    c5: z >= 0.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    # MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x", "y", "z"],
    #                        ["c3", "c4", "c5"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))
    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c3), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c5), 1.0)
    # test FiniteConstraintRef
    @testset "FiniteConstraintRef" begin
        @test map_dual(c1, Val(:TransData)) == -1.
    end
    # test InfiniteConstraintRef
    @testset "InfiniteConstraintRef" begin
        @test map_dual(c2, Val(:TransData)) == [0., 1.]
    end
end

# Test map_shadow_price
@testset "map_shadow_price" begin
    # initialize the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, inf(par))
    @global_variable(m, g)
    @constraint(m, c1, g <= 0)
    @constraint(m, c2, inf >= 0)
    # setup the mock solved transcription model
    tm = transcription_model(m)
    @variable(tm, x)
    @variable(tm, y)
    @variable(tm, z)
    @objective(tm, Min, x^2)
    @constraint(tm, c3, x <= 0)
    @constraint(tm, c4, y >= 0)
    @constraint(tm, c5, z >= 0)
    tm.ext[:TransData].global_to_var[g] = x
    tm.ext[:TransData].infinite_to_vars[inf] = [y, z]
    tm.ext[:TransData].finite_to_constr[c1] = c3
    tm.ext[:TransData].infinite_to_constrs[c2] = [c4, c5]
    modelstring = """
    variables: x, y, z
    minobjective: 1*x*x
    c3: x <= 0.0
    c4: y >= 0.0
    c5: z >= 0.0
    """
    model = MOIU.Model{Float64}()
    MOIU.loadfromstring!(model, modelstring)
    # MOIU.test_models_equal(JuMP.backend(tm).model_cache, model, ["x", "y", "z"],
    #                        ["c3", "c4", "c5"])
    JuMP.optimize!(tm, with_optimizer(MOIU.MockOptimizer,
                                     MOIU.Model{Float64}(),
                                     eval_objective_value=false))
    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c3), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c5), 1.0)
    # test FiniteConstraintRef
    @testset "FiniteConstraintRef" begin
        @test map_shadow_price(c1, Val(:TransData)) == -1.
    end
    # test InfiniteConstraintRef
    @testset "InfiniteConstraintRef" begin
        @test map_shadow_price(c2, Val(:TransData)) == [-0., -1.]
    end
end
