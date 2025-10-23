# Test general queries
@testset "General Queries" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, g)
    @objective(m, Min, g^2)
    @constraint(m, c1, 2 * inf * g <= 1)
    tb = m.backend

    # test not optimized yet
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not solved yet" begin 
        for f in (solve_time, simplex_iterations,
                  barrier_iterations, node_count)
            @test_throws ErrorException f(m)
        end
        @test termination_status(m) == MOI.OPTIMIZE_NOT_CALLED
        @test result_count(m) == 0
        @test raw_status(m) == "optimize not called"
        @test primal_status(m) == MOI.NO_SOLUTION
        @test dual_status(m) == MOI.NO_SOLUTION
    end

    # setup the results
    JuMP.optimize!(m)
    mockoptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ResultCount(), Int64(2))
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.PrimalStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.SimplexIterations(), Int64(4))
    MOI.set(mockoptimizer, MOI.BarrierIterations(), Int64(7))
    MOI.set(mockoptimizer, MOI.RelativeGap(), 9.0)
    MOI.set(mockoptimizer, MOI.NodeCount(), Int64(2))
    MOI.set(mockoptimizer, MOI.SolveTimeSec(), 0.42)

    # make model to set fallbacks
    test_model = InfiniteModel(TestBackend())
    set_transformation_backend_ready(test_model, true)

    # test termination_status
    @testset "JuMP.termination_status" begin
        @test termination_status(m) == MOI.OPTIMAL
    end
    # test is_solved_and_feasible
    @testset "JuMP.is_solved_and_feasible" begin
        @test is_solved_and_feasible(m)
        @test !is_solved_and_feasible(test_model)
        @test is_solved_and_feasible(test_model, allow_almost = true, dual = true)
    end
    # test primal_status
    @testset "JuMP.primal_status" begin
        @test primal_status(m) == MOI.FEASIBLE_POINT
        @test primal_status(m, result = 2) == MOI.FEASIBLE_POINT
    end
    # test dual_status
    @testset "JuMP.dual_status" begin
        @test dual_status(m) == MOI.FEASIBLE_POINT
        @test dual_status(m, result = 2) == MOI.FEASIBLE_POINT
    end
    # test solve_time
    @testset "JuMP.solve_time" begin
        @test solve_time(m) == 0.42
    end
    # test raw_status
    @testset "JuMP.raw_status" begin
        @test raw_status(m) == "solver specific string"
    end
    # test simplex_iterations
    @testset "JuMP.simplex_iterations" begin
        @test simplex_iterations(m) == 4
    end
    # test barrier_iterations
    @testset "JuMP.barrier_iterations" begin
        @test barrier_iterations(m) == 7
    end
    # test node_count
    @testset "JuMP.node_count" begin
        @test node_count(m) == 2
    end
    # test relative_gap
    @testset "JuMP.relative_gap" begin
        @test relative_gap(m) == 9
    end
    # test result_count
    @testset "JuMP.result_count" begin
        @test result_count(m) == 2
    end
    # test fallbacks
    @testset "Fallbacks" begin
        for f in (raw_status, solve_time, simplex_iterations,
                  barrier_iterations, node_count, result_count)
            @test_throws ErrorException f(test_model)
        end
    end
    # test model not up to date
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not up-to-date" begin 
        for f in (solve_time, simplex_iterations,
                  barrier_iterations, node_count)
            @test_throws ErrorException f(m)
        end
        @test raw_status(m) == "optimize not called"
        @test result_count(m) == 0
    end
end

# Test objective queries
@testset "Objective Queries" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, g)
    @objective(m, Min, g^2)
    @constraint(m, c1, 2 * g <= 1)
    tb = m.backend
    JuMP.optimize!(m)
    # setup the results
    mockoptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.DualObjectiveValue(), -2.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 2.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -2.0)
    # make model to set fallbacks
    test_model = InfiniteModel(TestBackend())
    # test objective_bound
    @testset "JuMP.objective_bound" begin
        @test objective_bound(m) == 2.0
    end
    # test objective_value
    @testset "JuMP.objective_value" begin
        @test objective_value(m) == -1.
    end
    # test dual_objective_bound
    @testset "JuMP.dual_objective_value" begin
        @test dual_objective_value(m) == -2.
    end
    # test fallbacks
    @testset "Fallbacks" begin
        for f in (objective_bound, objective_value, dual_objective_value)
            @test_throws ErrorException f(test_model)
        end
    end
    # test model not up to date
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not up-to-date" begin 
        for f in (objective_bound, objective_value, dual_objective_value)
            @test_throws ErrorException f(m)
        end
    end
end

# Test variable queries
@testset "Variable Queries" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], 
                        derivative_method = OrthogonalCollocation(3))
    @infinite_parameter(m, par2 in [0, 2], supports = [0, 2])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par2, par))
    @finite_parameter(m, fin == 42)
    @parameter_function(m, a == sin(par))
    @variable(m, g <=3)
    d1 = @deriv(inf, par)
    rv = inf2(par2, 0)
    @objective(m, Min, g^2)
    @constraint(m, c1, 2 * g <= 1)
    tb = m.backend
    JuMP.optimize!(m)
    inft = transformation_variable(inf, label = All)
    gt = transformation_variable(g)
    inf2t = transformation_variable(inf2, label = All)
    d1t = transformation_variable(d1, label = All)
    rvt = transformation_variable(rv, label = All)
    at = transformation_variable(a, label = All)
    cref = UpperBoundRef(g)
    creft = transformation_constraint(cref, label = All)
    # setup the optimizer
    mockoptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    for t in list_of_constraint_types(tb.model)
        for c in all_constraints(tb.model, t...)
            MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c), -12.0)
        end
    end
    # MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(creft), 7.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(gt), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inft[1]), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inft[2]), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inft[3]), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(rvt[1]), -2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(rvt[2]), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(d1t[1]), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(d1t[2]), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(d1t[3]), 2.0)
    # test has_values
    @testset "JuMP.has_values" begin
        @test has_values(m)
    end
    # test map_value
    @testset "map_value" begin
        @test InfiniteOpt.map_value(inf, tb, label = All) == [2., 1., 2.]
        @test InfiniteOpt.map_value(g, tb) == 1.
        @test InfiniteOpt.map_value(rv, tb, label = All) == [-2., -1.]
        @test_throws ErrorException InfiniteOpt.map_value(g, TestBackend())
    end
    # test _get_value 
    @testset "_get_value " begin
        @test InfiniteOpt._get_value(g, FiniteVariableIndex) == 1.
        @test InfiniteOpt._get_value(par, IndependentParameterIndex, label = All) == [0., 0.5, 1.]
        @test InfiniteOpt._get_value(fin, FiniteParameterIndex) == 42
    end
    # test value
    @testset "JuMP.value" begin
        @test value(inf) == [2., 2.]
        @test value(inf, label = All) == [2., 1., 2.]
        @test value(d1) == [2., 2.]
        @test value(d1, label = All) == [2., 1., 2.]
        @test value(g) == 1.
        @test value(rv) == [-2., -1.]
        @test value(par) == [0., 1.]
        @test value(par, label = All) == [0., 0.5, 1.]
        @test value(fin) == 42
        @test value(fin, label = All) == 42
        @test value(a) == sin.([0., 1.])
        @test value(a, label = All) == sin.([0., 0.5, 1.])
        @test value(sin(0.0)) == sin(0.0)
    end
    #test Reduced Cost
    @testset "map_reduced_cost" begin
        @test InfiniteOpt.map_reduced_cost(inf, tb, label = All) == [0.0, 0.0, 0.0]
        @test InfiniteOpt.map_reduced_cost(g, tb) == 26.0
        @test_throws ErrorException InfiniteOpt.map_reduced_cost(g, TestBackend())
    end
     #test Reduced Cost
     @testset "JuMP.reduced_cost" begin
        @test JuMP.reduced_cost(inf, label = All) == [0.0, 0.0, 0.0]
        @test JuMP.reduced_cost(g) == 26.0
    end
    # test map_optimizer_index
    @testset "map_optimizer_index" begin
        @test isa(InfiniteOpt.map_optimizer_index(g, tb), MOI.VariableIndex)
        @test isa(InfiniteOpt.map_optimizer_index(inf, tb), Vector{MOI.VariableIndex})
        @test isa(InfiniteOpt.map_optimizer_index(rv, tb), Vector{MOI.VariableIndex})
        @test_throws ErrorException InfiniteOpt.map_optimizer_index(inf, TestBackend())
    end
    # test optimizer_index
    @testset "JuMP.optimizer_index" begin
        @test isa(optimizer_index(g), MOI.VariableIndex)
        @test isa(optimizer_index(inf, label = InternalLabel), Vector{MOI.VariableIndex})
        @test isa(optimizer_index(rv), Vector{MOI.VariableIndex})
    end
    # test dual
    @testset "JuMP.dual" begin
        @test_throws ErrorException dual(g)
    end
    # test model not up to date
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not up-to-date" begin 
        @test_throws ErrorException value(inf)
        @test_throws ErrorException reduced_cost(inf)
        @test_throws ErrorException optimizer_index(inf)
    end
end

# Test expression/measure queries 
@testset "Expression/Measure Queries" begin 
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par2 in [0, 2], supports = [0, 2])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par2, par))
    @variable(m, g)
    meas1 = support_sum(2inf, par)
    meas2 = support_sum(inf2, par2)
    @objective(m, Min, g^2)
    @constraint(m, c1, 2 * g <= 1)
    tb = m.backend
    JuMP.optimize!(m)
    inft = transformation_variable(inf)
    gt = transformation_variable(g)
    inf2t = transformation_variable(inf2)
    # setup the optimizer
    mockoptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(gt), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inft[1]), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inft[2]), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inf2t[1]), 3.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inf2t[2]), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inf2t[3]), -3.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(inf2t[4]), -2.0)
    # test has_values
    @testset "JuMP.has_values" begin
        @test has_values(m)
    end
    # test map_value
    @testset "map_value" begin
        @test InfiniteOpt.map_value(meas1, tb) == 4.
        @test InfiniteOpt.map_value(meas2, tb) == [2., -5.]
        @test InfiniteOpt.map_value(3g - 1, tb) == 2.
        @test InfiniteOpt.map_value(inf^2 + g, tb) == [5., 1.]
        @test InfiniteOpt.map_value(zero(AffExpr) + 1, tb) == 1.
    end
    # test value
    @testset "JuMP.value" begin
        @test value(meas1, label = All) == 4.
        @test value(meas2, label = UserDefined) == [2., -5.]
        @test value(3g - 1) == 2.
        @test value(inf * inf + g - 2) == [3., -1.]
        @test value(zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) - 42) == -42.
        @test value(sin(g)) == sin(1)
        @test value(GenericNonlinearExpr{GeneralVariableRef}(:sin, Any[0])) == 0
    end
    # test dual
    @testset "JuMP.dual" begin
        @test_throws ErrorException dual(meas1)
    end
    # test model not up to date
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not up-to-date" begin 
        @test_throws ErrorException value(meas1)
        @test value(zero(AffExpr) + 1) == 1
    end
end

# Test constraint queries
@testset "Constraint Queries" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, g)
    @constraint(m, c1, g <= 0)
    @constraint(m, c2, inf >= 0)
    @constraint(m, c3, sin(g) == 0)
    @constraint(m, c4, sin(inf) == 0)
    @objective(m, Min, 42)
    tb = m.backend
    JuMP.optimize!(m)
    inft = transformation_variable(inf)
    gt = transformation_variable(g)
    c1t = transformation_constraint(c1)
    c2t = transformation_constraint(c2)
    c3t = transformation_constraint(c3)
    c4t = transformation_constraint(c4)
    # setup optimizer info
    mockoptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(inft[1]), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(inft[2]), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(gt), 1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1t), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c2t[1]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c2t[2]), 1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c3t), 4.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4t[1]), 2.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4t[2]), 3.0)
    # test map_value
    @testset "map_value" begin
        @test InfiniteOpt.map_value(c1, tb) == 1.
        @test InfiniteOpt.map_value(c2, tb) == [-1., 0.]
    end
    # test value
    @testset "JuMP.value" begin
        @test value(c1) == 1.
        @test value(c2, label = UserDefined) == [-1., 0.]
        @test value(c3) == sin(1)
        @test value(c4) == [sin(-1), sin(0)]
    end
    # test map_optimizer_index
    @testset "map_optimizer_index" begin
        @test isa(InfiniteOpt.map_optimizer_index(c1, tb), MOI.ConstraintIndex)
        @test isa(InfiniteOpt.map_optimizer_index(c2, tb), Vector{<:MOI.ConstraintIndex})
        @test_throws ErrorException InfiniteOpt.map_optimizer_index(c1, TestBackend())
    end
    # test optimizer_index
    @testset "JuMP.optimizer_index" begin
        @test isa(optimizer_index(c1), MOI.ConstraintIndex)
        @test isa(optimizer_index(c2, label = All), Vector{<:MOI.ConstraintIndex})
        @test isa(optimizer_index(c3), MOI.ConstraintIndex)
        @test isa(optimizer_index(c4, label = All), Vector{<:MOI.ConstraintIndex})
    end
    # test has_values
    @testset "JuMP.has_duals" begin
        @test has_duals(m)
    end
    # test map_dual
    @testset "map_dual" begin
        @test InfiniteOpt.map_dual(c1, tb) == -1.
        @test InfiniteOpt.map_dual(c2, tb) == [0., 1.]
        @test_throws ErrorException InfiniteOpt.map_dual(c1, TestBackend())
    end
    # test dual
    @testset "JuMP.dual" begin
        @test dual(c1) == -1.
        @test dual(c2, label = UserDefined) == [0., 1.]
        @test dual(c3) == 4
        @test dual(c4) == [2, 3]
    end
    # test shadow_price
    @testset "JuMP.shadow_price" begin
        @test shadow_price(c1) == -1.
        @test shadow_price(c2, label = PublicLabel) == [-0., -1.]
        @test shadow_price(c3) == -4
        @test shadow_price(c4) == [-2, -3]
        @test_throws ErrorException InfiniteOpt.map_shadow_price(c1, TestBackend())
    end
    # test model not up to date
    set_objective_sense(m, MOI.MAX_SENSE)
    @testset "Not up-to-date" begin 
        @test_throws ErrorException dual(c1)
        @test_throws ErrorException shadow_price(c1)
        @test_throws ErrorException optimizer_index(c1)
    end
end

# Test LP sensitivity methods
@testset "LP Sensitivities" begin
    # Setup the infinite model
    optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                         eval_objective_value=false)
    m = InfiniteModel(optimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, g)
    @constraint(m, c1, g >= 0)
    @constraint(m, c2, inf >= 0)
    @constraint(m, c3, g <= 1)
    @constraint(m, c4, inf <= 1)
    @objective(m, Max, 2g)
    tb = m.backend
    optimize!(m)
    inft = transformation_variable(inf)
    gt = transformation_variable(g)
    c1t = transformation_constraint(c1)
    c2t = transformation_constraint(c2)
    c3t = transformation_constraint(c3)
    c4t = transformation_constraint(c4)
    # setup the optimizer info
    mockoptimizer = JuMP.backend(tb).optimizer.model;
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(gt), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(inft[1]), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(), JuMP.optimizer_index(inft[2]), 0.0)
    MOI.set(mockoptimizer, MOI.VariableBasisStatus(), JuMP.optimizer_index(gt), MOI.NONBASIC)
    MOI.set(mockoptimizer, MOI.VariableBasisStatus(), JuMP.optimizer_index(inft[1]), MOI.NONBASIC)
    MOI.set(mockoptimizer, MOI.VariableBasisStatus(), JuMP.optimizer_index(inft[2]), MOI.NONBASIC)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1t), -2.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c2t[1]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c2t[2]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c3t), 2.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4t[1]), 1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c4t[2]), 1.0)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c1t), MOI.BASIC)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c2t[1]), MOI.BASIC)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c2t[2]), MOI.BASIC)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c3t), MOI.BASIC)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c4t[1]), MOI.BASIC)
    MOI.set(mockoptimizer, MOI.ConstraintBasisStatus(), JuMP.optimizer_index(c4t[2]), MOI.BASIC)
    # test making the report 
    @test lp_sensitivity_report(m, atol = 1e-6) isa InfOptSensitivityReport
    @test_throws ErrorException lp_sensitivity_report(TestBackend())
    # test constraint queries
    @test lp_sensitivity_report(m)[c1] == (-Inf, 0)
    @test lp_sensitivity_report(m)[c2, label = All] == [(-Inf, 0), (-Inf, 0)]
    # test variable queries
    @test lp_sensitivity_report(m)[g] == (0, 0)
    @test lp_sensitivity_report(m)[inf, label = UserDefined] == [(0, 0), (0, 0)]
    # test model not up to date
    set_objective_sense(m, MOI.MIN_SENSE)
    @testset "Not up-to-date" begin 
        @test_throws ErrorException lp_sensitivity_report(m)
    end
end