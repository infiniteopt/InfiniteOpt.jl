# Test the core measure data structures
# TODO finish
@testset "Core Measure Datatypes" begin
    @test AbstractMeasureData isa DataType
    @test DiscreteMeasureData <: AbstractMeasureData
    # @test DiscreteMeasureData(pref, Int[], Int[], "", InfiniteOpt._w).name == ""
    @test_throws ErrorException DiscreteMeasureData(pref, Int[], Int[1], "",
                                                    InfiniteOpt._w)
    @test true
end
