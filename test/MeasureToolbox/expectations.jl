@testset "Internal Methods" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in Normal())
    @infinite_parameter(m, xi[1:2] in MvNormal(ones(2), [1. 0.; 0. 2.]))
    dx = dispatch_variable_ref(x)
    dxi = dispatch_variable_ref.(xi)
    set1 = UniDistributionSet(Normal())
    set2 = MultiDistributionSet(MvNormal(ones(2), [1. 0.; 0. 2.]))
    set3 = CollectionSet([set1, set1])
    # test _has_distribution_set
    @testset "_has_distribution_set" begin
        @test IOMT._has_distribution_set(dx, set1)
        @test IOMT._has_distribution_set(dxi[1], set2)
        @test IOMT._has_distribution_set(dxi[1], set3)
        @test !IOMT._has_distribution_set(dx, BadSet())
    end
    # test _expect_coeffs
    @testset "_expect_coeffs" begin
        @test IOMT._expect_coeffs([1, 2]) == [0.5, 0.5]
    end
end

# TODO test expect with JuMP container of infinite parameters

@testset "Expect" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in Normal())
    @infinite_parameter(m, xi[1:2] in MvNormal(ones(2), [1. 0.; 0. 2.]))
    @infinite_parameter(m, y in [0, 1])
    @infinite_variable(m, inf(x, xi, y))

    @test InfiniteOpt._index_type(expect(inf, x)) == MeasureIndex
    warn = "Cannot specify a nonzero `min_num_supports` for individual " *
           "dependent parameters."
    @test_logs (:warn, warn) expect(inf, xi[1], min_num_supports = 3)
    @test !InfiniteOpt._is_expect(InfiniteOpt._core_variable_object(expect(inf, y)).data)
end

@testset "Macro" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in Normal())
    @infinite_parameter(m, xi[1:2] in MvNormal(ones(2), [1. 0.; 0. 2.]))
    @infinite_variable(m, inf(x, xi))

    @test InfiniteOpt._index_type(@expect(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(@expect(inf, xi, min_num_supports = 20)) == MeasureIndex
    @test_macro_throws ErrorException @expect(inf, x, 1)
    @test_macro_throws ErrorException @expect(inf, x, bob = 35)
end
