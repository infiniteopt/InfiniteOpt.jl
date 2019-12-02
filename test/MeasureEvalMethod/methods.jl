@testset "Default Monte Carlo Sampling" begin
    # test univariate Monte Carlo sampling
    @testset "MC_sampling (univariate)" begin
        (supports, weights) = MC_sampling(0., 2., 5)
        @test length(supports) == 5
        @test length(weights) == 5
        @test weights == 0.4 * ones(5)
        @test all(supports .>= 0.)
        @test all(supports .<= 2.)
    end

    # test multivariate Monte Carlo sampling
    @testset "MC_sampling (multivariate)" begin
        (supports, weights) = MC_sampling([0., 0.], [1., 4.], 5)
        @test size(supports)[2] == 5
        @test length(weights) == 5
        @test weights == 0.8 * ones(5)
        @test all(supports[1, :] .>= 0.)
        @test all(supports[1, :] .<= 1.)
        @test all(supports[2, :] .>= 0.)
        @test all(supports[2, :] .<= 4.)
    end
end

@testset "Default Quadrature Methods" begin
    # test Gauss-Legendre method for bounded interval
    @testset "Gauss_Legendre" begin
        (supports, weights) = Gauss_Legendre(1., 5., 5)
        (expect_supports, expect_weights) = FGQ.gausslegendre(5)
        expect_supports = expect_supports * 2 .+ 3
        expect_weights = expect_weights * 2
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Hermite method for infinite interval
    @testset "Gauss_Hermite" begin
        @test_throws ErrorException Gauss_Hermite(0., Inf, 5)
        @test_throws ErrorException Gauss_Hermite(-Inf, 0., 5)
        @test_throws ErrorException Gauss_Hermite(Inf, Inf, 5)
        @test_throws ErrorException Gauss_Hermite(0., 1., 5)
        (supports, weights) = Gauss_Hermite(-Inf, Inf, 5)
        @test length(supports) == 5
        (expect_supports, expect_weights) = gausshermite(5)
        expect_weights = expect_weights .* exp.(expect_supports.^2)
        @test length(expect_weights) == 5
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test Gauss-Laguerre method for semi-infinite interval
    @testset "Gauss_Laguerre" begin
        @test_throws ErrorException Gauss_Laguerre(-Inf, Inf, 5)
        @test_throws ErrorException Gauss_Laguerre(0., 1., 5)
        (supports, weights) = Gauss_Laguerre(1., Inf, 5)
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
        (supports, weights) = Gauss_Laguerre(-Inf, 1., 5)
        (expect_supports, expect_weights) = gausslaguerre(5)
        expect_weights = expect_weights .* exp.(expect_supports)
        expect_supports = -expect_supports .+ 1.
        @test supports == expect_supports
        @test weights == expect_weights
    end

    # test the function that transform (semi-)infinite domain to finite domain
    @testset "infinite_transform" begin
        @test_throws ErrorException infinite_transform(0., 1., 5)
        (supports, weights) = infinite_transform(-Inf, Inf, 5,
                                                 sub_method = Gauss_Legendre)
        (expect_supports, expect_weights) = gausslegendre(5)
        expect_supports = MEM._default_x.(expect_supports, -Inf, Inf)
        expect_weights = expect_weights .* MEM._default_dx.(expect_weights,
                                                            -Inf, Inf)
        @test supports == expect_supports
        @test weights == expect_weights
    end
end
