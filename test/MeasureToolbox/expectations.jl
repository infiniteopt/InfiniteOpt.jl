@testset "Internal Methods" begin
    m = InfiniteModel()
    @infinite_parameter(m, x ~ Normal())
    @infinite_parameter(m, xi[1:2] ~ MvNormal(ones(2), [1. 0.; 0. 2.]))
    domain1 = UniDistributionDomain(Normal())
    domain2 = MultiDistributionDomain(MvNormal(ones(2), [1. 0.; 0. 2.]))
    domain3 = CollectionDomain([domain1, domain1])
    domain4 = BadDomain()
    domain5 = IntervalDomain(0, 1)
    domain6 = CollectionDomain([domain1, domain5])
    domain7 = IntervalDomain(0, Inf)
    domain8 = CollectionDomain([domain5, domain5])
    # test _expect_coeffs
    @testset "_expect_coeffs" begin
        @test IOMT._expect_coeffs([1, 2]) == [0.5, 0.5]
    end
    # test _default_pdf
    @testset "_default_pdf" begin
        @test IOMT._default_pdf(0, 1, 3) == 0.5 
    end
    # test generate_expect_data
    @testset "generate_expect_data" begin
        @test_throws ErrorException IOMT.generate_expect_data(domain4, x, 2)
        @test_throws ErrorException IOMT.generate_expect_data(domain1, x, 2, pdf = sin)
        @test IOMT.generate_expect_data(domain1, x, 2) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException IOMT.generate_expect_data(domain2, xi, 2, pdf = sin)
        @test IOMT.generate_expect_data(domain2, xi, 2) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException IOMT.generate_expect_data(domain6, xi, 2)
        @test_throws ErrorException IOMT.generate_expect_data(domain3, xi, 2, pdf = sin)
        @test IOMT.generate_expect_data(domain3, xi, 2) isa FunctionalDiscreteMeasureData
        @test_throws ErrorException IOMT.generate_expect_data(domain7, x, 2)
        @test_throws ErrorException IOMT.generate_expect_data(domain5, x, 2, d= 2)
        @test IOMT.generate_expect_data(domain5, x, 2) isa FunctionalDiscreteMeasureData
        @test IOMT.generate_expect_data(domain5, x, 2, pdf = (s) -> 0.5) isa FunctionalDiscreteMeasureData
        @test IOMT.generate_expect_data(domain8, xi[1], 2) isa FunctionalDiscreteMeasureData
    end
end

@testset "Expect" begin
    # Setup the model
    m = InfiniteModel()
    @infinite_parameter(m, x ~ Normal())
    @infinite_parameter(m, xi[1:2] ~ MvNormal(ones(2), [1. 0.; 0. 2.]))
    @infinite_parameter(m, y in [0, 1])
    @variable(m, inf, Infinite(x, xi, y))

    # Test normal usage
    @test InfiniteOpt._index_type(expect(inf, x)) == MeasureIndex
    warn = "Cannot specify a nonzero `num_supports` for individual " *
           "dependent parameters."
    @test_logs (:warn, warn) expect(inf, xi[1], num_supports = 3)
    @test !InfiniteOpt._is_expect(InfiniteOpt._core_variable_object(expect(inf, y)).data)
    @test 𝔼(inf, x) isa GeneralVariableRef

    # Test with JuMP container 
    @infinite_parameter(m, z[2:3] in [-1, 1])
    @test_throws ErrorException expect(z[2], z)
    @test expect(z[2], z[2]) isa GeneralVariableRef
end

@testset "Macro" begin
    # Setup the model
    m = InfiniteModel()
    @infinite_parameter(m, x ~ Normal())
    @infinite_parameter(m, xi[1:2] ~ MvNormal(ones(2), [1. 0.; 0. 2.]))
    @variable(m, inf, Infinite(x, xi))

    # Test normal usage
    @test InfiniteOpt._index_type(@expect(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(@expect(inf, xi, num_supports = 20)) == MeasureIndex
    @test_macro_throws ErrorException @expect(inf, x, 1)
    @test_macro_throws ErrorException @𝔼(inf, x, 5)
    @test @𝔼(inf, x) isa GeneralVariableRef
end
