# Test the core measure data structures
# TODO finish
@testset "Core Measure Datatypes" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    # test AbstractMeasureData
    @testset "AbstractMeasureData" begin
        @test AbstractMeasureData isa DataType
    end
    # test DiscreteMeasureData
    @testset "DiscreteMeasureData" begin
        @test DiscreteMeasureData <: AbstractMeasureData
        @test DiscreteMeasureData(par, [1], [1], "", InfiniteOpt._w).name == ""
        @test DiscreteMeasureData(par, [1], [1], "",
                                  InfiniteOpt._w).parameter_ref == par
        @test_throws ErrorException DiscreteMeasureData(par, Int[], [1], "",
                                                        InfiniteOpt._w)
        @test_throws ErrorException DiscreteMeasureData(par, [1], [2], "",
                                                        InfiniteOpt._w)
        @test_throws ErrorException DiscreteMeasureData(par, [1], [-1], "",
                                                        InfiniteOpt._w)
    end
    # test MultiDiscreteMeasureData
    @testset "MultiDiscreteMeasureData" begin
        @test MultiDiscreteMeasureData <: AbstractMeasureData
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0])
        @test MultiDiscreteMeasureData(pars, [1], [supp], "",
                                       InfiniteOpt._w).name == ""
        @test MultiDiscreteMeasureData(pars, [1], [supp], "",
                                       InfiniteOpt._w).parameter_ref == pars
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1],
                          JuMP.Containers.SparseAxisArray[], "", InfiniteOpt._w)
        @test_throws ErrorException MultiDiscreteMeasureData(pars, Int[],
                          JuMP.Containers.SparseAxisArray[], "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [-1, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [2, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
    end
    # test Measure
    @testset "Measure" begin

    end
end

# Test _has_variable
@testset "_has_variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    #TODO FINISH
end
