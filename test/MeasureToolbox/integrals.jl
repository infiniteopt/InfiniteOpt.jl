@testset "Integral Data Generation (Univariate)" begin
    m = InfiniteModel();
    @infinite_parameter(m, t in [-Inf, Inf], supports = [0, .5, 1])
    @infinite_parameter(m, x[1:2] in [-Inf, Inf])

    # test _trapezoid_coeff
    @testset "_trapezoid_coeff" begin
        @test InfiniteOpt.MeasureToolbox._trapezoid_coeff([1,2,3]) == [0.5, 1., 0.5]
    end
    # test generate_integral_data (trapezoid)
    @testset "generate_integral_data (trapezoid)" begin
        @test generate_integral_data(t, 0, 1, UniTrapezoid()) ==
                    FunctionalDiscreteMeasureData(t, IOMT._trapezoid_coeff, 0, All, 
                               NoGenerativeSupports(),default_weight, 0, 1, false)
    end
    # test _ensure_independent_param
    @testset "_ensure_independent_param" begin
        @test IOMT._ensure_independent_param(t, Automatic()) isa Nothing
        @test_throws ErrorException IOMT._ensure_independent_param(x[1], Automatic())
    end
    # test generate_integral_data (Gauss-Legendre)
    @testset "generate_integral_data (Gauss-Legendre)" begin
        @test generate_integral_data(t, 0, 1, GaussLegendre()) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, 1, GaussLegendre()).coefficients == FastGaussQuadrature.gausslegendre(3)[2] * 1/2
    end
    # test generate_integral_data (Gauss-Lobatto)
    @testset "generate_integral_data (Gauss-Lobatoo)" begin
        @test generate_integral_data(t, 0, 1, GaussLobatto()) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, 1, GaussLobatto()).coefficients == FastGaussQuadrature.gausslobatto(3)[2] * 1/2
    end

    # test generate_integral_data (FEGaussLobatto)
    @testset "generate_integral_data (FEGaussLobatto)" begin
        @test generate_integral_data(t, 0, 1, FEGaussLobatto()) isa FunctionalDiscreteMeasureData
        @test generate_integral_data(t, 0, 1, FEGaussLobatto()).label == All
        @test all(abs.(generate_integral_data(t, 0, 1, FEGaussLobatto()).coeff_function([0, 0.5, 1]) .- [0.08333333333333333,
                                                                                0.3333333333333333,
                                                                                0.16666666666666666,
                                                                                0.3333333333333333,
                                                                                0.08333333333333333]).<.00005) == true
        @infinite_variable(m, flob(t))
        integral(flob, t, 0, 1, eval_method = FEGaussLobatto(), num_nodes = 3)
        add_generative_supports(t)
        lobb_supps = supports(t, label = All)
        @test all(abs.(lobb_supps .- [0, .25, .5, .75, 1] .< .000005)) == true
        @test_throws ErrorException("Num supports must be greater than 2") integral(flob, t, 0, 1, eval_method = FEGaussLobatto(), num_nodes = 1)
    end
    # test generate_integral_data (Gauss-Radau)
    @testset "generate_integral_data (Gauss-Radau)" begin
        @test generate_integral_data(t, 0, 1, GaussRadau()) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, 1, GaussRadau()).coefficients == FastGaussQuadrature.gaussradau(3)[2] * 1/2
    end
    # test generate_integral_data (Gauss-Jacobi)
    @testset "generate_integral_data (Gauss-Jacobi)" begin
        @test generate_integral_data(t, 0, 1, GaussJacobi(5.0, 4.0)) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, 1, GaussJacobi(5.0, 4.0)).coefficients == FastGaussQuadrature.gaussjacobi(3, 5.0, 4.0)[2] * 1/2
        @test_throws ErrorException("α and β must be greater than -1 for GaussJacobi") generate_integral_data(t, 0, 1, GaussJacobi(-5.0, -5.0))
    end
     # test generate_integral_data (Gauss-Chebyshev)
    @testset "generate_integral_data (Gauss-Chebyshev)" begin
        @test generate_integral_data(t, 0, 1, GaussChebyshev(3)) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, 1, GaussChebyshev(3)).coefficients == FastGaussQuadrature.gausschebyshev(3, 3)[2] * 1/2
        @test_throws ErrorException("Order must be between 1 and 4 for GaussChebyshev") generate_integral_data(t, 0, 1, GaussChebyshev(5))
    end
    # test generate_integral_data (Gauss-Laguerre)
    @testset "generate_integral_data (Gauss-Laguerre)" begin
        @test generate_integral_data(t, -Inf, 0, GaussLaguerre()) isa DiscreteMeasureData
        @test generate_integral_data(t, 0, Inf, GaussLaguerre()) isa DiscreteMeasureData
    end
    # test generate_integral_data (Gauss-Hermite)
    @testset "generate_integral_data (Gauss-Hermite)" begin
        @test generate_integral_data(t, -Inf, Inf, GaussHermite()) isa DiscreteMeasureData
    end
    # test generate_integral_data (Quadrature)
    @testset "generate_integral_data (Quadrature)" begin
        @test generate_integral_data(t, 0, 1, Quadrature()) isa FunctionalDiscreteMeasureData
        @test generate_integral_data(t, -Inf, 0, Quadrature()) isa DiscreteMeasureData
        @test generate_integral_data(t, -Inf, Inf, Quadrature()) isa DiscreteMeasureData
    end
    # test handling of inappropriate quadrature method use
    @testset "generate_integral_data (improper quadrature method usage)" begin
        warn = "The `GaussLegendre` method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussLegendre())) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussLegendre())) isa DiscreteMeasureData
        warn = "The `GaussLobatto` method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussLobatto())) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussLobatto())) isa DiscreteMeasureData
        warn = "The Gauss Lobatto method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, FEGaussLobatto())) isa FunctionalDiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, FEGaussLobatto())) isa FunctionalDiscreteMeasureData
        warn = "The `GaussRadau` method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussRadau())) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussRadau())) isa DiscreteMeasureData
        warn = "The `GaussJacobi` method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussJacobi(5.0, 4.0))) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussJacobi(5.0, 4.0))) isa DiscreteMeasureData
        warn = "The `GaussChebyshev` method can only be applied on finite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussChebyshev(3))) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussChebyshev(3))) isa DiscreteMeasureData
        warn = "Gauss Laguerre quadrature can only be applied on semi-infinite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, 0, 1, GaussLaguerre())) isa DiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, Inf, GaussLaguerre())) isa DiscreteMeasureData
        warn = "Gauss Hermite quadrature can only be applied on infinite intervals, " *
               "switching to an appropriate method."
        @test (@test_logs (:warn, warn) generate_integral_data(t, 0, 1, GaussHermite())) isa FunctionalDiscreteMeasureData
        @test (@test_logs (:warn, warn) generate_integral_data(t, -Inf, 0, GaussHermite())) isa DiscreteMeasureData
    end
    # test generate_integral_data (uniform Monte Carlo)
    @testset "generate_integral_data (UniMCSampling())" begin
        @test generate_integral_data(t, 0, 1, UniMCSampling()) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException generate_integral_data(t, -Inf, 0, UniMCSampling())
        warn = "Cannot specify a nonzero minimum supports for an individual " *
               "dependent parameter. Setting `num_supports = 0`."
        @test_logs (:warn, warn) generate_integral_data(x[1], 0, 1, UniMCSampling(), num_supports = 1)
    end
    # test generate_integral_data (independent Monte Carlo)
    @testset "generate_integral_data (UniIndepMCSampling())" begin
        @test generate_integral_data(t, 0, 1, UniIndepMCSampling()) isa DiscreteMeasureData
        @test_throws ErrorException generate_integral_data(t, -Inf, 0, UniIndepMCSampling())
        @test_throws ErrorException generate_integral_data(t, -Inf, Inf, UniIndepMCSampling())
    end
    # test generate_integral_data (Automatic)
    @testset "generate_integral_data (Automatic)" begin
        @test generate_integral_data(t, 0, 1, Automatic()) isa FunctionalDiscreteMeasureData
        @test generate_integral_data(t, -Inf, 1, Automatic()) isa DiscreteMeasureData
        @test generate_integral_data(t, -Inf, Inf, Automatic()) isa DiscreteMeasureData
        @test generate_integral_data(x[1], 0, 1, Automatic()) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException generate_integral_data(x[1], -Inf, 0, Automatic())
        @test_throws ErrorException generate_integral_data(x[1], -Inf, Inf, Automatic())
    end
    # test fallback
    @testset "generate_integral_data (fallback)" begin
        @test_throws ErrorException generate_integral_data(x, 1, 2, Quadrature())
    end
end

@testset "Integral Data Generation (Multivariate)" begin
    m = InfiniteModel();
    @infinite_parameter(m, x[1:2] in [-Inf, Inf])
    # test generate_integral_data (MultiMCSampling)
    @testset "generate_integral_data (MultiMCSampling)" begin
        @test generate_integral_data(x, [0, 0], [1, 1], MultiMCSampling()) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException generate_integral_data(x, [0, 0], [Inf, 1], MultiMCSampling())
        @test_throws ErrorException("MC Sampling is not applicable to (semi-)infinite intervals.") generate_integral_data(x, [0, -Inf], [1, 1], MultiMCSampling())
    end
    # test _make_multi_mc_supports
    @testset "_make_multi_mc_supports" begin
        dx = dispatch_variable_ref.(x)
        @test IOMT._make_multi_mc_supports(dx, [0, 0], [1, 1], 5) isa Matrix{Float64}
        @test size(IOMT._make_multi_mc_supports(dx, [0, 0], [1, 1], 5)) == (2, 5)
    end
    # test generate_integral_data (MultiIndepMCSampling)
    @testset "generate_integral_data (MultiIndepMCSampling)" begin
        @test generate_integral_data(x, [0, 0], [1, 1], MultiIndepMCSampling()) isa DiscreteMeasureData
        @test_throws ErrorException generate_integral_data(x, [0, 0], [Inf, 1], MultiMCSampling())
        @test_throws ErrorException generate_integral_data(x, [0, -Inf], [1, 1], MultiMCSampling())
    end
    # test generate_integral_data (Automatic)
    @testset "generate_integral_data (Automatic)" begin
        @test generate_integral_data(x, [0, 0], [1, 1], Automatic()) isa FunctionalDiscreteMeasureData
    end
end

@testset "Univariate Integrals" begin
    m = InfiniteModel();
    @infinite_parameter(m, t in [0, 10])
    @infinite_variable(m, inf1(t))
    @testset "uni_integral_defaults" begin
        @test uni_integral_defaults() == IOMT.UniIntegralDefaults
    end
    @testset "set_uni_integral_defaults" begin
        @test uni_integral_defaults() == Dict{Symbol, Any}(:eval_method => Automatic())
        @test set_uni_integral_defaults(num_nodes = 5, new_kwarg = true) isa Nothing
        @test uni_integral_defaults()[:num_nodes] == 5
        @test uni_integral_defaults()[:new_kwarg]
        @test clear_uni_integral_defaults() isa Nothing
        @test uni_integral_defaults() == Dict{Symbol, Any}(:eval_method => Automatic())
    end
    @testset "integral" begin
        @test InfiniteOpt._index_type(integral(inf1, t)) == MeasureIndex
        @test_throws ErrorException integral(inf1, t, -1, 1)
        @test_throws ErrorException integral(inf1, t, 9, 11)
        @test_throws ErrorException integral(inf1, t, 3, 1)
        @test name(integral(inf1, t)) == "integral"
    end
    @testset "∫" begin
        @test InfiniteOpt._index_type(∫(inf1, t)) == MeasureIndex
        @test_throws ErrorException ∫(inf1, t, -1, 1)
        @test_throws ErrorException ∫(inf1, t, 9, 11)
        @test_throws ErrorException ∫(inf1, t, 3, 1)
        @test name(∫(inf1, t)) == "integral"
    end
end

@testset "Multivariate Integrals" begin
    m = InfiniteModel();
    @infinite_parameter(m, x[1:2] in [0, 10])
    @infinite_parameter(m, xi[1:2] in [0, 5])
    @infinite_parameter(m, y[2:3] in [0, 1])
    @infinite_variable(m, inf1(x))
    @infinite_variable(m, inf2(xi))
    @infinite_variable(m, inf3(y))
    @testset "multi_integral_defaults" begin
        @test multi_integral_defaults() == IOMT.MultiIntegralDefaults
    end
    @testset "set_multi_integral_defaults" begin
        @test multi_integral_defaults() == Dict{Symbol, Any}(:eval_method => Automatic())
        @test set_multi_integral_defaults(num_supports = 5, new_kwarg = true) isa Nothing
        @test multi_integral_defaults()[:num_supports] == 5
        @test multi_integral_defaults()[:new_kwarg]
        @test clear_multi_integral_defaults() isa Nothing
        @test multi_integral_defaults() == Dict{Symbol, Any}(:eval_method => Automatic())
    end
    @testset "integral" begin
        @test InfiniteOpt._index_type(integral(inf1, x)) == MeasureIndex
        @test InfiniteOpt._index_type(integral(inf2, xi)) == MeasureIndex
        @test InfiniteOpt._index_type(integral(inf1, x, 1, 2)) == MeasureIndex
        @test InfiniteOpt._index_type(integral(inf2, xi, 1, 2)) == MeasureIndex
        @test_throws ErrorException integral(inf3, y, [0, 0], [1, 1])
        @test_throws ErrorException integral(inf1, x, 3, 1)
        @test_throws ErrorException integral(inf1, x, -1, 1)
        @test_throws ErrorException integral(inf1, x, 9, 11)
        @test name(integral(inf1, x)) == "integral"
    end
    @testset "∫" begin
        @test InfiniteOpt._index_type(∫(inf1, x)) == MeasureIndex
        @test InfiniteOpt._index_type(∫(inf2, xi)) == MeasureIndex
        @test InfiniteOpt._index_type(∫(inf1, x, 1, 2)) == MeasureIndex
        @test InfiniteOpt._index_type(∫(inf2, xi, 1, 2)) == MeasureIndex
        @test_throws ErrorException ∫(inf3, y, [0, 0], [1, 1])
        @test_throws ErrorException ∫(inf1, x, 3, 1)
        @test_throws ErrorException ∫(inf1, x, -1, 1)
        @test_throws ErrorException ∫(inf1, x, 9, 11)
        @test name(∫(inf1, x)) == "integral"
    end
end

@testset "Macro" begin
    m = InfiniteModel();
    @infinite_parameter(m, t in [0, 10])
    @infinite_parameter(m, x[1:2] in [0, 10])
    @infinite_variable(m, inf1(t))
    @infinite_variable(m, inf2(x))
    @testset "@integral" begin
        @test InfiniteOpt._index_type(@integral(inf1, t)) == MeasureIndex
        @test InfiniteOpt._index_type(@integral(inf2, x)) == MeasureIndex
        @test InfiniteOpt._index_type(@integral(inf2, t)) == MeasureIndex
        @test InfiniteOpt._index_type(@integral(inf1, t, 3, 7)) == MeasureIndex
        @test InfiniteOpt._index_type(@integral(inf2, x, [0, 2], [3, 5])) == MeasureIndex
        @test_macro_throws ErrorException @integral(inf1, t, 1)
    end
    @testset "@∫" begin
        @test InfiniteOpt._index_type(@∫(inf1, t)) == MeasureIndex
        @test InfiniteOpt._index_type(@∫(inf2, x)) == MeasureIndex
        @test InfiniteOpt._index_type(@∫(inf2, t)) == MeasureIndex
        @test InfiniteOpt._index_type(@∫(inf1, t, 3, 7)) == MeasureIndex
        @test InfiniteOpt._index_type(@∫(inf2, x, [0, 2], [3, 5])) == MeasureIndex
        @test_macro_throws ErrorException @∫(inf1, t, 1)
    end
end
