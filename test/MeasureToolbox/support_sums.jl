@testset "Internal Methods" begin
    @testset "_support_sum_coeffs" begin
        @test IOMT._support_sum_coeffs([1, 2]) == [1, 1]
    end
end

@testset "Support Sum" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in [0, 1])
    @infinite_parameter(m, y[1:2] in [0, 1])
    @infinite_variable(m, inf(x, y))

    @test InfiniteOpt._index_type(support_sum(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(support_sum(inf, y)) == MeasureIndex
    @test name(support_sum(inf, x)) == "support_sum"
end

@testset "Macro" begin
    m = InfiniteModel()
    @infinite_parameter(m, x in [0, 1])
    @infinite_parameter(m, y[1:2] in [0, 1])
    @infinite_variable(m, inf(x, y))

    @test InfiniteOpt._index_type(@support_sum(inf, x)) == MeasureIndex
    @test InfiniteOpt._index_type(@support_sum(inf, y)) == MeasureIndex
end
