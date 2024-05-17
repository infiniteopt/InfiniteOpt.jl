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
        expected = T(0, [5, 5]) + T(0, [0, 0])
        @test isequal_canonical(InfiniteOpt.make_reduced_expr(meas, t, 0.0, m), expected)
        # test InfiniteVariableIndex 
        @test isequal(InfiniteOpt.make_reduced_expr(q, t, 0.0, m), q(0))
        @test isequal(InfiniteOpt.make_reduced_expr(T, t, 0.0, m), T(0, x))
        # test DerivativeIndex 
        dT = deriv(T, x[1])
        dq = deriv(q, t)
        @test isequal(InfiniteOpt.make_reduced_expr(dq, t, 0.0, m), dq(0))
        v = dT(t, [x[1], 0])
        @test isequal(InfiniteOpt.make_reduced_expr(dT, x[2], 0.0, m), v)
        # test SemiInfiniteVariableIndex
        v2 = dT(t, [0, 0])
        @test isequal(InfiniteOpt.make_reduced_expr(v, x[1], 0.0, m), v2)
        @test isequal(InfiniteOpt.make_reduced_expr(v2, t, 0.0, m), dT(0, [0, 0]))
    end
     # test allows_high_order_derivatives
     @testset "allows_high_order_derivatives" begin 
        @test_throws ErrorException allows_high_order_derivatives(TestMethod())
        @test !allows_high_order_derivatives(OrthogonalCollocation(2))
        @test allows_high_order_derivatives(FiniteDifference(Central()))
    end
end

# Test evaluate_derivative and its helpers 
@testset "evaluate_derivative" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10], num_supports = 5)
    @infinite_parameter(m, x[1:2] in [0, 1], num_supports = 3)
    @variable(m, y, Infinite(t))
    @variable(m, q, Infinite(t, x))
    meas = support_sum(q, t)
    d1 = deriv(y, t)
    d2 = deriv(q, t)
    d3 = deriv(q, x[2])
    d4 = deriv(meas, x[1])
    d5 = deriv(y, t, t)
    d6 = deriv(y, t, t, t)
    d7 = @deriv(y, t^4)
    # test fallback 
    @testset "Fallback" begin 
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d1, y, TestMethod(), m)
    end
    # test _make_difference_expr for Forward 
    @testset "_make_difference_expr (Forward)" begin 
        supps = supports(t, label = All)
        # test 1st order
        expected = 2.5d1(0) - y(2.5) + y(0)
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d1, y, t, 1, 0, supps, m, Forward()), expected)
        # test 2nd order
        expected = 6.25d5(0) - y(5) + 2y(2.5) - y(0)
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d5, y, t, 2, 0, supps, m, Forward()), expected)
    end
    # test _make_difference_expr for Central
    @testset "_make_difference_expr (Central)" begin 
        supps = supports(t, label = All)
        # test 1st order
        expected = 5d2(2.5, x) - q(5, x) + q(0, x)
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d2, q, t, 1, 1, supps, m, Central()), expected)
        # test 2nd order
        expected = 6.25d5(2.5) - y(5) + 2y(2.5) - y(0)
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d5, y, t, 2, 1, supps, m, Central()), expected)
    end
    # test _make_difference_expr for Backward
    @testset "_make_difference_expr (Backward)" begin 
        supps = sort(supports(x[2], label = All))
        # test 1st order
        expected = 0.5d3(t, [x[1], 0.5]) - q(t, [x[1], 0.5]) + q(t, [x[1], 0])
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d3, q, x[2], 1, 1, supps, m, Backward()), expected)
        # test 2nd order
        supps = supports(t, label = All)
        expected = 6.25d5(5) - y(5) + 2y(2.5) - y(0)
        @test isequal_canonical(InfiniteOpt._make_difference_expr(d5, y, t, 2, 1, supps, m, Backward()), expected)
    end
    # test _make_difference_expr fallback
    @testset "_make_difference_expr (Fallback)" begin 
        @test_throws ErrorException InfiniteOpt._make_difference_expr(d3, q, x[2], 1, 2, [1], m, BadFiniteTech())
    end
    # test FiniteDifference with evaluate_derivative
    @testset "FiniteDifference" begin 
        # test with independent parameter 
        method = FiniteDifference(Central())
        exprs = [5d1(2.5) - y(5) + y(0), 5d1(5) - y(7.5) + y(2.5), 5d1(7.5) - y(10) + y(5)]
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(d1, y, method, m), exprs)
        exprs = [2.5^4 * d7(5) - y(10) + 4y(7.5) - 6y(5) + 4y(2.5) - y(0)]
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(d7, y, method, m), exprs)
        # test with dependent parameter 
        method = FiniteDifference(Forward()) 
        exprs = [0.5d4([s, x[2]]) - sum(q(i, [s+0.5, x[2]]) for i in supports(t)) + 
                 sum(q(i, [s, x[2]]) for i in supports(t)) for s in (0.5, 0)]  
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(d4, meas, method, m), exprs)
        # test using Backward without boundary constraint
        method = FiniteDifference(Backward(), false)
        exprs = [2.5d1(s) - y(s) + y(s-2.5) for s in (2.5, 5, 7.5)]
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(d1, y, method, m), exprs)
        # test Backward with boundary constraint 
        method = FiniteDifference(Backward(), true)
        exprs = [2.5d1(s) - y(s) + y(s-2.5) for s in (2.5, 5, 7.5, 10)]
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(d1, y, method, m), exprs)
        # test bad order for Central
        @test_throws ErrorException evaluate_derivative(d6, y, FiniteDifference(Central()), m)
        # test without supports 
        delete_supports(x)
        @test_throws ErrorException InfiniteOpt.evaluate_derivative(d3, q, method, m)
        # test parameter function 
        f = parameter_function(sin, t)
        df = deriv(f, t)
        exprs = [2.5df(s) - sin(s) + sin(s-2.5) for s in (2.5, 5, 7.5, 10)]
        @test isequal_canonical(InfiniteOpt.evaluate_derivative(df, f, method, m), exprs)
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
        set_supports(t, [0, 5, 10], force = true)
        function rm_zeros(exs)
            for e in exs 
                for (v, c) in e.terms
                    if abs(c) < 1e-15
                        delete!(e.terms, v)
                    end
                end
            end
            return exs
        end
        @test isequal_canonical(rm_zeros(evaluate_derivative(d2, q, method, m)), exprs)
        @test supports(t) == [0, 5, 10]
        @test supports(t, label = All) == [0, 2.5, 5, 7.5, 10]
        # test resolve 
        @test isequal_canonical(rm_zeros(evaluate_derivative(d2, q, method, m)), exprs)
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
    @infinite_parameter(m, t in [0, 10], supports = [0, 3, 10], 
                        derivative_method = FiniteDifference(Central()))
    @infinite_parameter(m, x in [-1, 1], supports = [-1, 1], 
                        derivative_method = OrthogonalCollocation(4))
    @variable(m, y, Infinite(t))
    @variable(m, q, Infinite(x))
    dy = @deriv(y, t)
    dq = @deriv(q, x)
    dy2 = @deriv(y, t^2)
    dq2 = @deriv(q, x^3)
    # test evaluate 
    @testset "evaluate" begin 
        # test evaluate first derivative 
        @test derivative_method(t) isa FiniteDifference
        @test evaluate(dy) isa Nothing 
        @test num_constraints(m) == 1
        @test !has_generative_supports(t)
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(dy)
        @test length(derivative_constraints(dy)) == 1
        # test evaluate second derivative 
        @test derivative_method(t) isa FiniteDifference
        @test evaluate(dy2) isa Nothing 
        @test num_constraints(m) == 2
        @test !has_generative_supports(t)
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(dy2)
        @test length(derivative_constraints(dy2)) == 1
        # test recursive 1st order derivative reformulation
        @test evaluate(dq2) isa Nothing
        @test num_derivatives(m) == 5
        @test num_constraints(m) == 11
        @test has_generative_supports(x)
        @test has_derivative_constraints(dq2)
        @test length(derivative_constraints(dq2)) == 3
    end
    # test evaluate_all_derivatives!
    @testset "evaluate_all_derivatives!" begin 
        @test derivative_method(x) isa OrthogonalCollocation
        @test evaluate_all_derivatives!(m) isa Nothing 
        @test num_constraints(m) == 11
        @test has_generative_supports(x)
        @test has_internal_supports(x)
        @test supports(x) == [-1, 1]
        @test num_supports(x, label = All) == 4
        @test supports(t) == [0, 3, 10]
        @test length(derivative_constraints(dy)) == 1
        @test length(derivative_constraints(dy2)) == 1
        @test length(derivative_constraints(dq)) == 3
        @test has_derivative_constraints(t)
        @test has_derivative_constraints(x)
        @test has_derivative_constraints(dy)
        @test has_derivative_constraints(dy2)
        @test has_derivative_constraints(dq)
    end
end
