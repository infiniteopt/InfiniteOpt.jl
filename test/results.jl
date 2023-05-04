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
    tm = transcription_model(m)
    JuMP.optimize!(m)
    # setup the results
    mockoptimizer = JuMP.backend(tm).optimizer.model
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
    # test termination_status
    @testset "JuMP.termination_status" begin
        @test termination_status(m) == MOI.OPTIMAL
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
    tm = transcription_model(m)
    JuMP.optimize!(m)
    # setup the results
    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.DualObjectiveValue(), -2.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 2.0)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.ConstraintDual(), JuMP.optimizer_index(c1), -2.0)
    @testset "JuMP.objective_bound" begin
        @test objective_bound(m) == 2.0
    end
    # test objective_bound
    @testset "JuMP.objective_value" begin
        @test objective_value(m) == -1.
    end
    # test dual_objective_bound
    @testset "JuMP.dual_objective_value" begin
        @test dual_objective_value(m) == -2.
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
    var = build_variable(error, inf2, Dict{Int, Float64}(2 => 0))
    rv = add_variable(m, var)
    @objective(m, Min, g^2)
    @constraint(m, c1, 2 * g <= 1)
    tm = transcription_model(m)
    JuMP.optimize!(m)
    inft = transcription_variable(inf, label = All)
    gt = transcription_variable(g)
    inf2t = transcription_variable(inf2, label = All)
    d1t = transcription_variable(d1, label = All)
    rvt = transcription_variable(rv, label = All)
    cref = UpperBoundRef(g)
    creft = transcription_constraint(cref, label = All)
    # setup the optimizer
    mockoptimizer = JuMP.backend(tm).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.ResultCount(), 1)
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    for t in list_of_constraint_types(tm)
        for c in all_constraints(tm, t...)
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
        @test InfiniteOpt.map_value(inf, Val(:TransData), 1, label = All) == [2., 1., 2.]
        @test InfiniteOpt.map_value(g, Val(:TransData), 1) == 1.
        @test InfiniteOpt.map_value(rv, Val(:TransData), 1, label = All) == [-2., -1.]
    end
    # test _get_value 
    @testset "_get_value " begin
        @test InfiniteOpt._get_value(g, FiniteVariableIndex, 1) == 1.
        @test InfiniteOpt._get_value(par, IndependentParameterIndex, 1, label = All) == [0., 0.5, 1.]
        @test InfiniteOpt._get_value(fin, FiniteParameterIndex, 1) == 42
    end
    # test value
    @testset "JuMP.value" begin
        @test value(inf) == [2., 2.]
        @test value(inf, label = All) == [2., 1., 2.]
        @test value(inf, label = All, ndarray = true) == [2., 1., 2.]
        @test value(d1) == [2., 2.]
        @test value(d1, label = All) == [2., 1., 2.]
        @test value(g) == 1.
        @test value(rv) == [-2., -1.]
        @test value(par) == [0., 1.]
        @test value(par, label = All) == [0., 0.5, 1.]
        @test value(fin) == 42
        @test value(fin, label = All) == 42
        @test value(a) == [sin(0.), sin(1.)]
        @test value(a, label = All) == [sin(0.), sin(0.5), sin(1.)]
    end
    #test Reduced Cost
    @testset "map_reduced_cost" begin
    @test InfiniteOpt.map_reduced_cost(inf, Val(:TransData), label = All) == [0.0, 0.0, 0.0]
    @test InfiniteOpt.map_reduced_cost(g, Val(:TransData)) == 26.0
end
     #test Reduced Cost
     @testset "JuMP.reduced_cost" begin
     @test JuMP.reduced_cost(inf, label = All) == [0.0, 0.0, 0.0]
     @test JuMP.reduced_cost(g) == 26.0
 end
    # test map_optimizer_index
    @testset "map_optimizer_index" begin
        @test isa(InfiniteOpt.map_optimizer_index(g, Val(:TransData)), MOI.VariableIndex)
        @test isa(InfiniteOpt.map_optimizer_index(inf, Val(:TransData)), Vector{MOI.VariableIndex})
        @test isa(InfiniteOpt.map_optimizer_index(rv, Val(:TransData)), Vector{MOI.VariableIndex})
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
    tm = transcription_model(m)
    JuMP.optimize!(m)
    inft = transcription_variable(inf)
    gt = transcription_variable(g)
    inf2t = transcription_variable(inf2)
    # setup the optimizer
    mockoptimizer = JuMP.backend(tm).optimizer.model
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
        @test InfiniteOpt.map_value(meas1, Val(:TransData), 1) == 4.
        @test InfiniteOpt.map_value(meas2, Val(:TransData), 1) == [0., -3.]
        @test InfiniteOpt.map_value(3g - 1, Val(:TransData), 1) == 2.
        @test InfiniteOpt.map_value(inf^2 + g, Val(:TransData), 1) == [5., 1.]
        @test InfiniteOpt.map_value(zero(AffExpr) + 1, Val(:TransData), 1) == 1.
    end
    # test value
    @testset "JuMP.value" begin
        @test value(meas1, label = All) == 4.
        @test value(meas2, label = UserDefined) == [0., -3.]
        @test value(3g - 1) == 2.
        @test value(inf * inf + g - 2) == [3., -1.]
        @test value(inf * inf + g - 2, ndarray = true) == [3., -1.]
        @test value(zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) - 42) == -42.
        # @test value(sin(g)) == sin(1)
        # @test_throws ErrorException value(NonlinearExpr(:sin, Any[0]))
    end
    # test dual
    @testset "JuMP.dual" begin
        @test_throws ErrorException dual(meas1)
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
    tm = transcription_model(m)
    JuMP.optimize!(m)
    inft = transcription_variable(inf)
    gt = transcription_variable(g)
    c1t = transcription_constraint(c1)
    c2t = transcription_constraint(c2)
    c3t = transcription_constraint(c3)
    c4t = transcription_constraint(c4)
    # setup optimizer info
    mockoptimizer = JuMP.backend(tm).optimizer.model
    # block = MOI.get(tm, MOI.NLPBlock())
    # MOI.initialize(block.evaluator, Symbol[])
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
    # MOI.set(mockoptimizer, MOI.NLPBlockDual(1), [4.0, 2., 3.])
    # test map_value
    @testset "map_value" begin
        @test InfiniteOpt.map_value(c1, Val(:TransData), 1) == 1.
        @test InfiniteOpt.map_value(c2, Val(:TransData), 1) == [-1., 0.]
    end
    # test value
    @testset "JuMP.value" begin
        @test value(c1) == 1.
        @test value(c2, label = UserDefined) == [-1., 0.]
        @test value(c2, label = UserDefined, ndarray = true) == [-1., 0.]
        # @test value(c3) == sin(1)
        # @test value(c4) == [sin(-1), sin(0)]
    end
    # test map_optimizer_index
    @testset "map_optimizer_index" begin
        @test isa(InfiniteOpt.map_optimizer_index(c1, Val(:TransData)), MOI.ConstraintIndex)
        @test isa(InfiniteOpt.map_optimizer_index(c2, Val(:TransData)), Vector{<:MOI.ConstraintIndex})
    end
    # test optimizer_index
    @testset "JuMP.optimizer_index" begin
        @test isa(optimizer_index(c1), MOI.ConstraintIndex)
        @test isa(optimizer_index(c2, label = All), Vector{<:MOI.ConstraintIndex})
        @test isa(optimizer_index(c2, label = All, ndarray = true), Vector{<:MOI.ConstraintIndex})
        @test isa(optimizer_index(c3), MOI.ConstraintIndex)
        @test isa(optimizer_index(c4, label = All), Vector{<:MOI.ConstraintIndex})
    end
    # test has_values
    @testset "JuMP.has_duals" begin
        @test has_duals(m)
    end
    # test map_dual
    @testset "map_dual" begin
        @test InfiniteOpt.map_dual(c1, Val(:TransData), 1) == -1.
        @test InfiniteOpt.map_dual(c2, Val(:TransData), 1) == [0., 1.]
    end
    # test dual
    @testset "JuMP.dual" begin
        @test dual(c1) == -1.
        @test dual(c2, label = UserDefined) == [0., 1.]
        @test dual(c2, label = UserDefined, ndarray = true) == [0., 1.]
        @test dual(c3) == 4
        @test dual(c4) == [2, 3]
    end
    # test shadow_price
    @testset "JuMP.shadow_price" begin
        # test with FEASIBILITY_SENSE
        @test_throws ErrorException shadow_price(c1)
        @test_throws ErrorException shadow_price(c2)
        @test_throws ErrorException shadow_price(c3)
        # test with MIN_SENSE
        set_objective_sense(m, MOI.MIN_SENSE)
        @test shadow_price(c1) == -1.
        @test shadow_price(c1, ndarray = true) == [-1.]
        @test shadow_price(c2, label = PublicLabel) == [-0., -1.]
        @test shadow_price(c3) == -4
        @test shadow_price(c4) == [-2, -3]
        # test with MAX_SENSE
        set_objective_sense(m, MOI.MAX_SENSE)
        @test shadow_price(c1) == 1.
        @test shadow_price(c1, ndarray = true) == [1.]
        @test shadow_price(c2, label = PublicLabel) == [0., 1.]
        @test shadow_price(c3) == 4
        @test shadow_price(c4) == [2, 3]
        # test no duals error
        MOI.set(mockoptimizer, MOI.DualStatus(), MOI.NO_SOLUTION)
        @test_throws ErrorException shadow_price(c1)
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
    tm = transcription_model(m)
    optimize!(m)
    inft = transcription_variable(inf)
    gt = transcription_variable(g)
    c1t = transcription_constraint(c1)
    c2t = transcription_constraint(c2)
    c3t = transcription_constraint(c3)
    c4t = transcription_constraint(c4)
    # setup the optimizer info
    mockoptimizer = JuMP.backend(tm).optimizer.model;
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
    # test constraint queries
    @test lp_sensitivity_report(m)[c1] == (-Inf, 0)
    @test lp_sensitivity_report(m)[c2, label = All] == [(-Inf, 0), (-Inf, 0)]
    @test lp_sensitivity_report(m)[c2, ndarray = true] == [(-Inf, 0), (-Inf, 0)]
    # test variable queries
    @test lp_sensitivity_report(m)[g] == (0, 0)
    @test lp_sensitivity_report(m)[inf, label = UserDefined] == [(0, 0), (0, 0)]
    @test lp_sensitivity_report(m)[inf, ndarray = true] == [(0, 0), (0, 0)]
end
