@testset "Default Monte Carlo Sampling" begin
    # test univariate Monte Carlo sampling
    @testset "mc_sampling (univariate)" begin
        (supports, weights) = mc_sampling(0., 2., 5)
        @test length(supports) == 5
        @test length(weights) == 5
        @test weights == 0.4 * ones(5)
        @test all(supports .>= 0.)
        @test all(supports .<= 2.)
        (supports, weights) = mc_sampling(-Inf, Inf, 5)
        @test length(supports) == 5
        @test length(weights) == 5
    end

    # test multivariate Monte Carlo sampling
    @testset "mc_sampling (multivariate)" begin
        lb = convert(JuMPC.SparseAxisArray, [0., 0.])
        ub = convert(JuMPC.SparseAxisArray, [1., 4.])
        (supports, weights) = mc_sampling(lb, ub, 5)
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
        (supports, weights) = gauss_legendre(1., 5., 5)
        (expect_supports, expect_weights) = FGQ.gausslegendre(5)
        expect_supports = expect_supports * 2 .+ 3
        expect_weights = expect_weights * 2
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Hermite method for infinite interval
    # at the compiling step FGQ.gausshermite creates huge allocation
    @testset "gauss_hermite" begin
        @test_throws ErrorException gauss_hermite(0., Inf, 5)
        @test_throws ErrorException gauss_hermite(-Inf, 0., 5)
        @test_throws ErrorException gauss_hermite(Inf, Inf, 5)
        @test_throws ErrorException gauss_hermite(0., 1., 5)
        (supports, weights) = gauss_hermite(-Inf, Inf, 5)
        @test length(supports) == 5
        (expect_supports, expect_weights) = FGQ.gausshermite(5)
        expect_weights = expect_weights .* exp.(expect_supports.^2)
        @test length(expect_weights) == 5
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Laguerre method for semi-infinite interval
    @testset "gauss_laguerre" begin
        @test_throws ErrorException gauss_laguerre(-Inf, Inf, 5)
        @test_throws ErrorException gauss_laguerre(0., 1., 5)
        (supports, weights) = gauss_laguerre(1., Inf, 5)
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
        (supports, weights) = gauss_laguerre(-Inf, 1., 5)
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = -expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
    end
end

# test the function that transform (semi-)infinite domain to finite domain
# with the default transform function
@testset "Infinite transform for (semi-)infinite domain" begin
    @testset "infinite_transform" begin
        @test_throws ErrorException infinite_transform(0., 1., 5)
        # Infinite domain
        (supports, weights) = infinite_transform(-Inf, Inf, 5,
                                                 sub_method = gauss_legendre)
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_weights = expect_weights .* MEM._default_dx.(expect_supports,
                                                            -Inf, Inf)
        expect_supports = MEM._default_x.(expect_supports, -Inf, Inf)
        @test supports == expect_supports
        @test weights == expect_weights
        # Lower bounded semi-infinite domain
        (supports, weights) = infinite_transform(0., Inf, 5,
                                                 sub_method = gauss_legendre)
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_supports = (expect_supports .+ 1) ./ 2
        expect_weights = expect_weights ./ 2
        expect_weights = expect_weights .* MEM._default_dx.(expect_supports,
                                                            0., Inf)
        expect_supports = MEM._default_x.(expect_supports, 0., Inf)
        @test supports == expect_supports
        @test weights == expect_weights
        # Upper bounded semi-infinite domain
        (supports, weights) = infinite_transform(-Inf, 0., 5,
                                                 sub_method = gauss_legendre)
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

@testset "Registry functions" begin
    struct NotASetType end
    struct NewSet <: AbstractInfiniteSet end
    function new_fn end
    m = InfiniteModel()
    @testset "_set_type" begin
        @test MEM._set_type(m, NewSet()) == nothing
        @test MEM._set_type(m, IntervalSet(0., 1.)) == IntervalSet
    end
    @testset "_set_method_check" begin
        @test_throws ErrorException MEM._set_method_check(m, NewSet(), new_fn)
        @test_throws ErrorException MEM._set_method_check(m, IntervalSet(0., 1.),
                                                          new_fn)
    end
    @testset "set_method_registry" begin
        @test_throws ErrorException set_method_registry(m, NotASetType,
                                                        mc_sampling)
        @test set_method_registry(m, NewSet, [new_fn, mc_sampling]) == nothing
        @test set_method_registry(m, DistributionSet, new_fn) == nothing
        @test m.meas_method_registry[DistributionSet] ==
              Set{Function}([new_fn, mc_sampling])
        @test m.meas_method_registry[NewSet] ==
              Set{Function}([new_fn, mc_sampling])
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
        (expect_supports, expect_weights) = gauss_legendre(0.3, 0.7, 10)
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
                                    BadSet(), x, 10, 0., 1., mc_sampling)
    end
end
