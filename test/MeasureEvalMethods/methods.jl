m = InfiniteModel();
@infinite_parameter(m, t in [-Inf, Inf])
@infinite_parameter(m, x[1:2] in [-Inf, Inf])

@testset "Default Monte Carlo Sampling" begin
    # test univariate Monte Carlo sampling
    @testset "mc_sampling (univariate)" begin
        (supports, weights) = generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, 0., 2., Val(mc_sampling))
        @test length(supports) == 5
        @test length(weights) == 5
        @test weights == 0.4 * ones(5)
        @test all(supports .>= 0.)
        @test all(supports .<= 2.)
        (supports, weights) = generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, -Inf, Inf, Val(mc_sampling))
        @test length(supports) == 5
        @test length(weights) == 5
    end

    # test multivariate Monte Carlo sampling
    @testset "mc_sampling (multivariate)" begin
        lb = convert(JuMPC.SparseAxisArray, [0., 0.])
        ub = convert(JuMPC.SparseAxisArray, [1., 4.])
        (supports, weights) = generate_supports_and_coeffs(
                                  IntervalSet(-Inf, Inf), x, 5, lb, ub, Val(mc_sampling))
        @test length(supports) == 5
        @test length(weights) == 5
        @test weights == 0.8 * ones(5)
        @test all([i >= lb for i in supports])
        @test all([i <= ub for i in supports])
    end
end

@testset "Quadrature Methods" begin
    # test Gauss-Legendre method for bounded interval
    @testset "gauss_legendre" begin
        (supports, weights) = generate_supports_and_coeffs(
                               IntervalSet(-Inf, Inf), t, 5, 1., 5., Val(gauss_legendre))
        (expect_supports, expect_weights) = FGQ.gausslegendre(5)
        expect_supports = expect_supports * 2 .+ 3
        expect_weights = expect_weights * 2
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Hermite method for infinite interval
    # at the compiling step FGQ.gausshermite creates huge allocation
    @testset "gauss_hermite" begin
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, 0., Inf, Val(gauss_hermite))
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, -Inf, 0., Val(gauss_hermite))
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, Inf, Inf, Val(gauss_hermite))
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, 0., 1., Val(gauss_hermite))
        (supports, weights) = generate_supports_and_coeffs(
                             IntervalSet(-Inf, Inf), t, 5, -Inf, Inf, Val(gauss_hermite))
        @test length(supports) == 5
        (expect_supports, expect_weights) = FGQ.gausshermite(5)
        expect_weights = expect_weights .* exp.(expect_supports.^2)
        @test length(expect_weights) == 5
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Laguerre method for semi-infinite interval
    @testset "gauss_laguerre" begin
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, -Inf, Inf, Val(gauss_laguerre))
        @test_throws ErrorException generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, 0., 1., Val(gauss_laguerre))
        (supports, weights) = generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, 1., Inf, Val(gauss_laguerre))
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
        (supports, weights) = generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 5, -Inf, 1., Val(gauss_laguerre))
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = -expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
    end
end

@testset "trapezoid Rule" begin
    (supports, weights) = generate_supports_and_coeffs(IntervalSet(-Inf, Inf), t, 3, 0., 2., Val(trapezoid))
    @test length(supports) == 3
    @test length(weights) == 3
    @test weights == [0.5, 1., 0.5]
    @test supports == [0., 1., 2.]
end

# test the function that transform (semi-)infinite domain to finite domain
# with the default transform function
@testset "Infinite transform for (semi-)infinite domain" begin
    @testset "infinite_transform" begin
        @test_throws ErrorException infinite_transform(IntervalSet(-Inf, Inf), t, 5, 0., 1.)
        # Infinite domain
        (supports, weights) = infinite_transform(IntervalSet(-Inf, Inf), t, 5, -Inf, Inf,
                                                 Val(gauss_legendre))
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_weights = expect_weights .* MEM._default_dx.(expect_supports,
                                                            -Inf, Inf)
        expect_supports = MEM._default_x.(expect_supports, -Inf, Inf)
        @test supports == expect_supports
        @test weights == expect_weights
        # Lower bounded semi-infinite domain
        (supports, weights) = infinite_transform(IntervalSet(-Inf, Inf), t, 5, 0., Inf,
                                                 Val(gauss_legendre))
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_supports = (expect_supports .+ 1) ./ 2
        expect_weights = expect_weights ./ 2
        expect_weights = expect_weights .* MEM._default_dx.(expect_supports,
                                                            0., Inf)
        expect_supports = MEM._default_x.(expect_supports, 0., Inf)
        @test supports == expect_supports
        @test weights == expect_weights
        # Upper bounded semi-infinite domain
        (supports, weights) = infinite_transform(IntervalSet(-Inf, Inf), t, 5, -Inf, 0.,
                                                 Val(gauss_legendre))
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_supports = (expect_supports .+ 1) ./ 2
        expect_weights = expect_weights ./ 2
        expect_weights = expect_weights .* MEM._default_dx.(expect_supports,
                                                            -Inf, 0.)
        expect_supports = MEM._default_x.(expect_supports, -Inf, 0.)
        @test supports == expect_supports
        @test weights == expect_weights
    end
end

@testset "Data Generation" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in [0., 1.])
    @infinite_parameter(m, y[1:2] in [0., 1.])
    dist1 = Normal(0., 1.)
    dist2 = MvNormal([0., 0.], [1. 0.;0. 10.])
    @infinite_parameter(m, β in dist1)
    @infinite_parameter(m, ω[1:2] in dist2)
    @testset "generate_measure_data (univariate)" begin
        @test_throws ErrorException generate_measure_data(x, 10, 1., 2.)
        measure_data = generate_measure_data(x, 10, 0.3, 0.7, eval_method = gauss_legendre)
        (expect_supports, expect_weights) = generate_supports_and_coeffs(IntervalSet(0., 1.), x, 10, 0.3, 0.7, Val(gauss_legendre))
        @test measure_data.supports == expect_supports
        @test measure_data.coefficients == expect_weights
        measure_data = generate_measure_data(β, 10, -1., Inf)
        @test length(measure_data.supports) == 10
        @test all([i >= -1. for i in measure_data.supports])
        measure_data = generate_measure_data(β, 10, -Inf, 1.)
        @test length(measure_data.supports) == 10
        @test all([i <= 1. for i in measure_data.supports])
    end

    @testset "generate_measure_data (multivariate)" begin
        lb = convert(JuMPC.SparseAxisArray, [0.3, 0.2])
        ub = convert(JuMPC.SparseAxisArray, [0.5, 0.9])
        measure_data = generate_measure_data(y, 10, lb, ub)
        @test length(measure_data.supports) == 10
        @test length(measure_data.coefficients) == 10
        @test sum(measure_data.coefficients) == 0.14
        warn = "Truncated distribution for multivariate distribution is " *
               "not supported. Lower bounds and upper bounds are ignored."
        @test_logs (:warn, warn) generate_measure_data(ω, 10,
                                  convert(JuMPC.SparseAxisArray, [-1., -10.]),
                                  convert(JuMPC.SparseAxisArray, [1., 10.]))
    end

    @testset "generate_supports_and_coeffs for unextended types" begin
        @test_throws ErrorException generate_supports_and_coeffs(
                                    BadSet(), x, 10, 0., 1., Val(mc_sampling))
    end
end
