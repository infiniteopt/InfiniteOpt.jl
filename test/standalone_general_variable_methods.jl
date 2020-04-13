# test Base Extensions
@testset "Base Extensions" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, PointVariableIndex)
    gvref2 = GeneralVariableRef(m, 2, PointVariableIndex)
    ivref = InfiniteVariableRef(m, InfiniteVariableIndex(1))
    hvref = HoldVariableRef(m, HoldVariableIndex(1))
    # test Base.copy
    @testset "Base.copy" begin
        @test copy(gvref) == gvref
        @test copy(ivref) == ivref
    end
    # test Base.:(==)
    @testset "Base.:(==)" begin
        @test gvref == gvref
        @test ivref == ivref
        @test !(gvref == ivref)
        @test !(hvref == ivref)
        @test !(gvref == gvref2)
    end
    # test Base.broadcastable
    @testset "Base.broadcastable" begin
        @test Base.broadcastable(gvref) isa Base.RefValue{GeneralVariableRef}
        @test Base.broadcastable(ivref) isa Base.RefValue{InfiniteVariableRef}
    end
    # test Base.length
    @testset "Base.length" begin
        @test length(gvref) == 1
        @test length(ivref) == 1
    end
    # test JuMP.isequal_canonical
    @testset "JuMP.isequal_canonical" begin
        @test JuMP.isequal_canonical(gvref, gvref)
        @test JuMP.isequal_canonical(ivref, ivref)
    end
    # test JuMP.variable_type
    @testset "JuMP.variable_type" begin
        @test JuMP.variable_type(m) == GeneralVariableRef
    end
end

# test Reference Accessers
@testset "Reference Accessers" begin
    # Setup data
    m = InfiniteModel();
    gvref1 = GeneralVariableRef(m, 1, PointVariableIndex)
    gvref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
    ivref = InfiniteVariableRef(m, InfiniteVariableIndex(1))
    # test JuMP.index (GeneralVariableRef)
    @testset "JuMP.index (GeneralVariableRef)" begin
        @test index(gvref1) == PointVariableIndex(1)
        @test index(gvref2) == DependentParameterIndex(DependentParametersIndex(1), 2)
    end
    # test JuMP.index (DispatchVariableRef)
    @testset "JuMP.index (DispatchVariableRef)" begin
        @test index(ivref) == InfiniteVariableIndex(1)
    end
    # test JuMP.owner_model (GeneralVariableRef)
    @testset "JuMP.owner_model (GeneralVariableRef)" begin
        @test owner_model(gvref1) == m
        @test owner_model(gvref2) == m
    end
    # test JuMP.owner_model (DispatchVariableRef)
    @testset "JuMP.owner_model (DispatchVariableRef)" begin
        @test owner_model(ivref) == m
    end
end
