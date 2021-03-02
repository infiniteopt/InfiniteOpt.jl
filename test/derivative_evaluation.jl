# Test all of the helper methods
@testset "Helper Methods" begin
    # Setup 
    m = InfiniteModel();
    @infinite_parameter(m, t in [0,1])
    @infinite_parameter(m, x[1:2] in [0,5], num_supports = 2)
    @infinite_variable(m, T(t, x))
    @infinite_variable(m, q(t))
    oc = OrthogonalCollocation(3, GaussLobatto())
    bad_oc = OrthogonalCollocation{GaussLegendre}(3, GaussLegendre())
    tref = dispatch_variable_ref(t)
    # test generative_support_info
    @testset "generative_support_info" begin 
        # test Fallback
        @test_throws ErrorException generative_support_info(TestGenMethod())
        # test NonGenerativeDerivativeMethod
        @test generative_support_info(TestMethod()) == NoGenerativeSupports()
        # test OrthogonalCollocation with GaussLobatto
        @test generative_support_info(oc) == UniformGenerativeInfo([0.5], InternalGaussLobatto)
        @test generative_support_info(OrthogonalCollocation(2)) == NoGenerativeSupports()
        # test bad OrthogonalCollocation
        @test_throws ErrorException generative_support_info(bad_oc)
    end
    # test support_label
    @testset "support_label" begin 
        # test Fallback
        @test_throws ErrorException support_label(TestMethod())
        # test GenerativeDerivativeMethod
        @test support_label(oc) == InternalGaussLobatto
    end
    # test make_reduced_expr 
    @testset "make_reduced_expr" begin 
        # test MeasureIndex based 
        meas = support_sum(T, x)
        pidx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, pidx, PointVariableIndex),
               GeneralVariableRef(m, pidx + 1, PointVariableIndex)]
        @test InfiniteOpt.make_reduced_expr(meas, t, 0.0, m) in [pts[1] + pts[2], pts[2] + pts[1]] # order is unknown
        @test parameter_values(pts[1]) in [(0., [0., 0.]), (0., [5., 5.])]
        @test parameter_values(pts[2]) in [(0., [0., 0.]), (0., [5., 5.])]
        # test InfiniteVariableIndex 
        pt = GeneralVariableRef(m, pidx + 2, PointVariableIndex)
        @test InfiniteOpt.make_reduced_expr(q, t, 0.0, m) == pt 
        @test parameter_values(pt) == (0.0,)
        ridx = 3 # 2 were previously made and deleted by the measure call
        rv = GeneralVariableRef(m, ridx, ReducedVariableIndex)
        @test InfiniteOpt.make_reduced_expr(T, t, 0.0, m) == rv
        @test parameter_list(rv) == [x[1], x[2]]
        @test eval_supports(rv)[1] == 0
        # test DerivativeIndex 
        dT = @deriv(T, x[1])
        dq = @deriv(q, t)
        pt = GeneralVariableRef(m, pidx + 3, PointVariableIndex)
        @test InfiniteOpt.make_reduced_expr(dq, t, 0.0, m) == pt 
        @test parameter_values(pt) == (0.0,)
        rv = GeneralVariableRef(m, ridx + 1, ReducedVariableIndex)
        @test InfiniteOpt.make_reduced_expr(dT, x[2], 0.0, m) == rv
        @test parameter_list(rv) == [t, x[1]]
        @test eval_supports(rv)[3] == 0
        # test ReducedVariableIndex
        rv2 = GeneralVariableRef(m, ridx + 2, ReducedVariableIndex)
        @test InfiniteOpt.make_reduced_expr(rv, x[1], 0.0, m) == rv2
        @test parameter_list(rv2) == [t]
        @test eval_supports(rv2)[3] == 0
        @test eval_supports(rv2)[2] == 0
        pt = GeneralVariableRef(m, pidx + 4, PointVariableIndex)
        @test InfiniteOpt.make_reduced_expr(rv2, t, 0.0, m) == pt 
        @test parameter_values(pt) == (0.0, [0.0, 0.0])
    end
end

# Test evaluate_derivative and its helpers 
@testset "evaluate_derivative" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10], num_supports = 3)
    @infinite_parameter(m, x[1:2] in [0, 1], num_supports = 2)
    @infinite_variable(m, y(t))
    @infinite_variable(m, q(t, x))
    meas = support_sum(q, t)
    d1 = @deriv(y, t)
    d2 = @deriv(q, t)
    d3 = @deriv(q, x[2])
    d4 = @deriv(meas, x[1])
    # test fallback 
    @testset "Fallback" begin 
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d1, TestMethod(), m)
    end
    # test _make_difference_expr for Forward 
    @testset "_make_difference_expr (Forward)" begin 
        supps = supports(t, label = All)
        pidx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, pidx, PointVariableIndex),
               GeneralVariableRef(m, pidx + 1, PointVariableIndex),
               GeneralVariableRef(m, pidx + 2, PointVariableIndex)]
        @test InfiniteOpt._make_difference_expr(d1, y, t, 1, supps, m, Forward()) == 5pts[1] - pts[2] + pts[3]
        @test parameter_values(pts[1]) == (0.,)
        @test parameter_values(pts[2]) == (5.,)
        @test parameter_values(pts[3]) == (0.,)
    end
    # test _make_difference_expr for Central
    @testset "_make_difference_expr (Central)" begin 
        supps = supports(t, label = All)
        ridx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rvs = [GeneralVariableRef(m, ridx, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 1, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 2, ReducedVariableIndex)]
        @test InfiniteOpt._make_difference_expr(d2, q, t, 2, supps, m, Central()) == 10rvs[1] - rvs[2] + rvs[3]
        @test eval_supports(rvs[1])[1] == 5
        @test eval_supports(rvs[2])[1] == 10
        @test eval_supports(rvs[3])[1] == 0
    end
    # test _make_difference_expr for FDBackward
    @testset "_make_difference_expr (Backward)" begin 
        supps = sort(supports(x[2], label = All))
        ridx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rvs = [GeneralVariableRef(m, ridx, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 1, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 2, ReducedVariableIndex)]
        @test InfiniteOpt._make_difference_expr(d3, q, x[2], 2, supps, m, Backward()) == rvs[1] - rvs[2] + rvs[3]
        @test eval_supports(rvs[1])[3] == 1
        @test eval_supports(rvs[2])[3] == 1
        @test eval_supports(rvs[3])[3] == 0
    end
    # test _make_difference_expr fallback
    @testset "_make_difference_expr (Fallback)" begin 
        @test_throws ErrorException InfiniteOpt._make_difference_expr(d3, q, x[2], 2, [1], m, BadFiniteTech())
    end
    # test FiniteDifference with evaluate_derivative
    @testset "FiniteDifference" begin 
        # test with independent parameter 
        method = FiniteDifference(Central())
        pidx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, pidx, PointVariableIndex),
               GeneralVariableRef(m, pidx + 1, PointVariableIndex),
               GeneralVariableRef(m, pidx + 2, PointVariableIndex)]
        exprs = [10pts[1] - pts[2] + pts[3]]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test with dependent parameter 
        method = FiniteDifference(Forward()) 
        ridx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rvs = [GeneralVariableRef(m, ridx, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 4, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 5, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 6, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 10, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 11, ReducedVariableIndex),
               GeneralVariableRef(m, ridx + 12, ReducedVariableIndex)]
        exprs = [rvs[1] - rvs[2] - rvs[3] - rvs[4] + rvs[5] + rvs[6] + rvs[7]] 
        @test InfiniteOpt.evaluate_derivative(d4, method, m) == exprs
        # test using Backward without boundary constraint
        method = FiniteDifference(Backward(), false)
        pidx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, pidx, PointVariableIndex),
               GeneralVariableRef(m, pidx + 1, PointVariableIndex),
               GeneralVariableRef(m, pidx + 2, PointVariableIndex)]
        exprs = [5pts[1] - pts[2] + pts[3]]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test Backward with boundary constraint 
        method = FiniteDifference(Backward(), true)
        pidx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, pidx, PointVariableIndex),
               GeneralVariableRef(m, pidx + 1, PointVariableIndex),
               GeneralVariableRef(m, pidx + 2, PointVariableIndex),
               GeneralVariableRef(m, pidx + 3, PointVariableIndex),
               GeneralVariableRef(m, pidx + 4, PointVariableIndex),
               GeneralVariableRef(m, pidx + 5, PointVariableIndex)]
        exprs = [5pts[1] - pts[2] + pts[3],
                 5pts[4] - pts[5] + pts[6]]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test without supports 
        delete_supports(x)
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d3, method, m)
        # test parameter function 
        f = parameter_function(sin, t)
        df = deriv(f, t)
        pts = [GeneralVariableRef(m, pidx + 6, PointVariableIndex),
               GeneralVariableRef(m, pidx + 7, PointVariableIndex)]
        @test InfiniteOpt.evaluate_derivative(df, method, m) == [5pts[1] - sin(5) + sin(0), 
                                                                 5pts[2] - sin(10) + sin(5)]
    end
    # test OrthogonalCollocation with evaluate_derivative
    @testset "OrthogonalCollocation" begin 
        # test first solve 
        method = OrthogonalCollocation(3)
        set_derivative_method(t, method)
        Mt = [1. 5.; 1. 10]' \ [2.5 2.5^2; 5. 25.]'
        M = Mt'
        ridx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 6
        rvs = [GeneralVariableRef(m, ridx + i, ReducedVariableIndex) for i = 1:16]
        exprs = [@expression(m, M[1, 1] * rvs[1] + M[1, 2] * rvs[2] - rvs[3] + rvs[4]),
                 @expression(m, M[2, 1] * rvs[5] + M[2, 2] * rvs[6] - rvs[7] + rvs[8]),
                 @expression(m, M[1, 1] * rvs[9] + M[1, 2] * rvs[10] - rvs[11] + rvs[12]),
                 @expression(m, M[2, 1] * rvs[13] + M[2, 2] * rvs[14] - rvs[15] + rvs[16])]
        @test InfiniteOpt.evaluate_derivative(d2, method, m) == exprs
        @test supports(t) == [0, 5, 10]
        @test supports(t, label = All) == [0, 2.5, 5, 7.5, 10]
        # test resolve 
        ridx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 6
        rvs = [GeneralVariableRef(m, ridx + i, ReducedVariableIndex) for i = 1:16]
        exprs = [@expression(m, M[1, 1] * rvs[1] + M[1, 2] * rvs[2] - rvs[3] + rvs[4]),
                 @expression(m, M[2, 1] * rvs[5] + M[2, 2] * rvs[6] - rvs[7] + rvs[8]),
                 @expression(m, M[1, 1] * rvs[9] + M[1, 2] * rvs[10] - rvs[11] + rvs[12]),
                 @expression(m, M[2, 1] * rvs[13] + M[2, 2] * rvs[14] - rvs[15] + rvs[16])]
        @test InfiniteOpt.evaluate_derivative(d2, method, m) == exprs
        @test supports(t) == [0, 5, 10]
        @test supports(t, label = All) == [0, 2.5, 5, 7.5, 10]
        @test has_generative_supports(t)
        @test has_internal_supports(t)
    end
end

# Test User Methods
@testset "User Methods" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10], supports = [0, 3, 10])
    @infinite_parameter(m, x in [-1, 1], supports = [-1, 1], 
                        derivative_method = OrthogonalCollocation(4))
    @infinite_variable(m, y(t))
    @infinite_variable(m, q(x))
    dy = @deriv(y, t)
    dq = @deriv(q, x)
    dy2 = @deriv(y, t^2)
    # test evaluate 
    @testset "evaluate" begin 
        # test evaluate first derivative 
        @test derivative_method(t) isa FiniteDifference
        @test evaluate(dy) isa Nothing 
        @test num_constraints(m) == 2
        @test !has_generative_supports(t)
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(dy)
        @test length(derivative_constraints(dy)) == 2
        # test evaluate second derivative 
        @test derivative_method(t) isa FiniteDifference
        @test evaluate(dy2) isa Nothing 
        @test num_constraints(m) == 4
        @test !has_generative_supports(t)
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(dy2)
        @test length(derivative_constraints(dy2)) == 2
    end
    # test evaluate_all_derivatives!
    @testset "evaluate_all_derivatives!" begin 
        @test derivative_method(x) isa OrthogonalCollocation
        @test evaluate_all_derivatives!(m) isa Nothing 
        @test num_constraints(m) == 7
        @test has_generative_supports(x)
        @test has_internal_supports(x)
        @test supports(x) == [-1, 1]
        @test num_supports(x, label = All) == 4
        @test supports(t) == [0, 3, 10]
        @test length(derivative_constraints(dy)) == 2
        @test length(derivative_constraints(dy2)) == 2
        @test length(derivative_constraints(dq)) == 3
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(x)
        @test has_derivative_constraints(dy)
        @test has_derivative_constraints(dy2)
        @test has_derivative_constraints(dq)
    end
end
