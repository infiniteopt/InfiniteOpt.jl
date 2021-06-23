# Test all of the helper methods
@testset "Helper Methods" begin
    # Setup 
    m = InfiniteModel();
    @infinite_parameter(m, t in [0,1])
    @infinite_parameter(m, x[1:2] in [0,5], num_supports = 2)
    @variable(m, T, Infinite(t, x))
    @variable(m, q, Infinite(t))
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
        expected = [T(0, [0, 0]) + T(0, [5, 5]), T(0, [5, 5]) + T(0, [0, 0])]
        @test InfiniteOpt.make_reduced_expr(meas, t, 0.0, m) in expected # order is unknown
        # test InfiniteVariableIndex 
        @test InfiniteOpt.make_reduced_expr(q, t, 0.0, m) == q(0)
        @test InfiniteOpt.make_reduced_expr(T, t, 0.0, m) == T(0, x)
        # test DerivativeIndex 
        dT = deriv(T, x[1])
        dq = deriv(q, t)
        @test InfiniteOpt.make_reduced_expr(dq, t, 0.0, m) == dq(0)
        v = dT(t, [x[1], 0])
        @test InfiniteOpt.make_reduced_expr(dT, x[2], 0.0, m) == v
        # test SemiInfiniteVariableIndex
        v2 = dT(t, [0, 0])
        @test InfiniteOpt.make_reduced_expr(v, x[1], 0.0, m) == v2
        @test InfiniteOpt.make_reduced_expr(v2, t, 0.0, m) == dT(0, [0, 0])
    end
end

# Test evaluate_derivative and its helpers 
@testset "evaluate_derivative" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10], num_supports = 3)
    @infinite_parameter(m, x[1:2] in [0, 1], num_supports = 2)
    @variable(m, y, Infinite(t))
    @variable(m, q, Infinite(t, x))
    meas = support_sum(q, t)
    d1 = deriv(y, t)
    d2 = deriv(q, t)
    d3 = deriv(q, x[2])
    d4 = deriv(meas, x[1])
    # test fallback 
    @testset "Fallback" begin 
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d1, TestMethod(), m)
    end
    # test _make_difference_expr for Forward 
    @testset "_make_difference_expr (Forward)" begin 
        supps = supports(t, label = All)
        expected = 5d1(0) - y(5) + y(0)
        @test InfiniteOpt._make_difference_expr(d1, y, t, 1, supps, m, Forward()) == expected
    end
    # test _make_difference_expr for Central
    @testset "_make_difference_expr (Central)" begin 
        supps = supports(t, label = All)
        expected = 10d2(5, x) - q(10, x) + q(0, x)
        @test InfiniteOpt._make_difference_expr(d2, q, t, 2, supps, m, Central()) == expected
    end
    # test _make_difference_expr for FDBackward
    @testset "_make_difference_expr (Backward)" begin 
        supps = sort(supports(x[2], label = All))
        expected = d3(t, [x[1], 1]) - q(t, [x[1], 1]) + q(t, [x[1], 0])
        @test InfiniteOpt._make_difference_expr(d3, q, x[2], 2, supps, m, Backward()) == expected
    end
    # test _make_difference_expr fallback
    @testset "_make_difference_expr (Fallback)" begin 
        @test_throws ErrorException InfiniteOpt._make_difference_expr(d3, q, x[2], 2, [1], m, BadFiniteTech())
    end
    # test FiniteDifference with evaluate_derivative
    @testset "FiniteDifference" begin 
        # test with independent parameter 
        method = FiniteDifference(Central())
        exprs = [10d1(5) - y(10) + y(0)]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test with dependent parameter 
        method = FiniteDifference(Forward()) 
        exprs = [d4([0, x[2]]) - q(0, [1, x[2]]) - q(5, [1, x[2]]) - 
                 q(10, [1, x[2]]) + q(0, [0, x[2]]) + q(5, [0, x[2]]) + 
                 q(10, [0, x[2]])] 
        @test InfiniteOpt.evaluate_derivative(d4, method, m) == exprs
        # test using Backward without boundary constraint
        method = FiniteDifference(Backward(), false)
        exprs = [5d1(5) - y(5) + y(0)]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test Backward with boundary constraint 
        method = FiniteDifference(Backward(), true)
        exprs = [5d1(5) - y(5) + y(0), 5d1(10) - y(10) + y(5)]
        @test InfiniteOpt.evaluate_derivative(d1, method, m) == exprs 
        # test without supports 
        delete_supports(x)
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d3, method, m)
        # test parameter function 
        f = parameter_function(sin, t)
        df = deriv(f, t)
        exprs = [5df(5) - sin(5) + sin(0), 5df(10) - sin(10) + sin(5)]
        @test InfiniteOpt.evaluate_derivative(df, method, m) == exprs
    end
    # test OrthogonalCollocation with evaluate_derivative
    @testset "OrthogonalCollocation" begin 
        # test first solve 
        method = OrthogonalCollocation(3)
        set_derivative_method(t, method)
        Mt = [1. 5.; 1. 10]' \ [2.5 2.5^2; 5. 25.]'
        M = Mt'
        exprs = [@expression(m, M[1, 1] * d2(2.5, x) + M[1, 2] * d2(5, x) - q(2.5, x) + q(0, x)),
                 @expression(m, M[2, 1] * d2(2.5, x) + M[2, 2] * d2(2.5, x) - q(5, x) + q(0, x)),
                 @expression(m, M[1, 1] * d2(7.5, x) + M[1, 2] * d2(10, x) - q(7.5, x) + q(5, x)),
                 @expression(m, M[2, 1] * d2(7.5, x) + M[2, 2] * d2(7.5, x) - q(10, x) + q(5, x))]
        delete_supports(t, label = UserDefined)
        function rm_zeros(exs)
            for e in exs 
                filter!((v, c) -> abs(c) > 1e-15, e)
            end
            return exs
        end
        @test rm_zeros(nfiniteOpt.evaluate_derivative(d2, method, m)) == exprs
        @test supports(t) == [0, 5, 10]
        @test supports(t, label = All) == [0, 2.5, 5, 7.5, 10]
        # test resolve 
        @test rm_zeros(nfiniteOpt.evaluate_derivative(d2, method, m)) == exprs
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
    @variable(m, y, Infinite(t))
    @variable(m, q, Infinite(x))
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
