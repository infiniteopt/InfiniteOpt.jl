# Test the extensions of Base.convert
@testset "Convert Extensions" begin
    m = Model()
    @variable(m, x[1:2])
    # Array -> SparseAxisArray
    @testset "Array->Sparse" begin
        @test convert(JuMPC.SparseAxisArray,
                      x) isa JuMPC.SparseAxisArray
        @test length(convert(JuMPC.SparseAxisArray, x)) == length(x)
        inds = CartesianIndices(x)
        @test convert(JuMPC.SparseAxisArray,
                      x)[inds[1][1]] == x[inds[1]]
        @test convert(JuMPC.SparseAxisArray,
                      x)[inds[2][1]] == x[inds[2]]
    end
    # DenseAxisArray -> SparseAxisArray
    @testset "DenseArray->Sparse" begin
        @variable(m, y[2:3])
        @test convert(JuMPC.SparseAxisArray,
                      y) isa JuMPC.SparseAxisArray
        @test length(convert(JuMPC.SparseAxisArray, y)) == length(y)
        inds = collect(keys(y))
        @test convert(JuMPC.SparseAxisArray,
                      y)[inds[1][1]] == y[inds[1]]
        @test convert(JuMPC.SparseAxisArray,
                      y)[inds[2][1]] == y[inds[2]]
    end
    # Number -> AffExpr
    @testset "Number->AffExpr" begin
        expr = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) + 42
        @test isequal(convert(JuMP.AbstractJuMPScalar, 42), expr)
    end
end

# Test extension of keys
@testset "_keys" begin
    m = Model()
    @variable(m, x[1:2], container = SparseAxisArray)
    @test InfiniteOpt._keys(x) == keys(x.data)
    @test InfiniteOpt._keys(ones(2)) == keys(ones(2))
end

# Test allequal
@testset "_allequal" begin
    @test InfiniteOpt._allequal([1])
    @test !InfiniteOpt._allequal([1, 2])
    @test InfiniteOpt._allequal([1, 1])
    @test InfiniteOpt._allequal([1 1; 1 1])
end

# Test Base.:(==) for VariableInfo
@testset "VariableInfo Comparison" begin
    info1 = VariableInfo(true, 1, false, NaN, false, NaN, false, NaN, false, true)
    info2 = VariableInfo(true, 1., false, NaN, false, NaN, false, NaN, false, true)
    info3 = VariableInfo(true, 1, false, NaN, false, NaN, false, NaN, false, true)
    info4 = VariableInfo(true, 1, false, NaN, false, NaN, false, NaN, false, false)
    @test info1 == info3 
    @test info1 == info2
    @test info1 != info4
end

# Test _make_float_info
@testset "_make_float_info" begin
    info1 = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    num = Float64(0)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    @test InfiniteOpt._make_float_info(info1) isa VariableInfo{Float64, Float64, Float64, Float64}
    @test InfiniteOpt._make_float_info(info2) == info2
end

# Test _kwargs2string
@testset "_kwargs2string" begin
    nt1 = (d = 1,)
    nt2 = (d = 1, e = 2)
    @test InfiniteOpt._kwargs2string(nt1) == "d = 1"
    @test InfiniteOpt._kwargs2string(nt2) == "d = 1, e = 2"
end
