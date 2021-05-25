@testset "Internal Methods" begin
    @testset "_support_sum_coeffs" begin
        @test IOMT._support_sum_coeffs([1, 2]) == [1, 1]
    end
end

@testset "Support Sum" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in [0, 1])
    @infinite_parameter(m, y[1:2] in [0, 1])
    @variable(m, inf, Infinite(x, y))

    @test InfiniteOpt._index_type(support_sum(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(support_sum(inf, y)) == MeasureIndex
    @test name(support_sum(inf, x)) == "support_sum"
    @test support_sum(inf, y, label = UserDefined) isa GeneralVariableRef
end

@testset "Macro" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in [0, 1])
    @infinite_parameter(m, y[1:2] in [0, 1])
    @variable(m, inf, Infinite(x, y))

    @test InfiniteOpt._index_type(@support_sum(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(@support_sum(inf, y)) == MeasureIndex
    @test_macro_throws ErrorException @support_sum(inf, y, 2)
end
