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
        @test convert(JuMP.AbstractJuMPScalar, 42) == expr
    end
end

# Test extension of keys
@testset "_keys" begin
    m = Model()
    @variable(m, x[1:2], container = SparseAxisArray)
    @test InfiniteOpt._keys(x) == keys(x.data)
    @test InfiniteOpt._keys(ones(2)) == keys(ones(2))
end

# Test extension of map
@testset "Base.map (DenseAxisArray)" begin
    JuMPC.@container(x[3:6], 2)
    @test map(i -> iseven(i), x) isa JuMPC.DenseAxisArray
    @test map(i -> iseven(i), x).data == [true, true, true, true]
end

# Test the _make_vector functions
@testset "_make_vector" begin
    x = [2, 4, 1, 8]
    # test with JuMP containers
    @testset "AbstractArray" begin
        y = convert(JuMPC.SparseAxisArray, x)
        @test sort(InfiniteOpt._make_vector(y)) == [1, 2, 4, 8]
        y = JuMPC.DenseAxisArray(x, [1, 2, 3, 4])
        @test InfiniteOpt._make_vector(y) == x
        @test InfiniteOpt._make_vector(zeros(2, 0)) == []
    end
    # test arrays
    @testset "Array" begin
        @test InfiniteOpt._make_vector(x) == x
        @test InfiniteOpt._make_vector(ones(2, 2)) == ones(4)
    end
end

# Test the _make_ordered_vector functions
@testset "_make_ordered_vector" begin
    x = [2, 4, 1, 8]
    # test with JuMP containers
    @testset "Vector" begin
        x = ones(2)
        @test InfiniteOpt._make_ordered_vector(x) === x
    end
    # test arrays
    @testset "Array" begin
        x = ones(2, 3)
        @test InfiniteOpt._make_ordered_vector(x) == ones(6)
        JuMPC.@container(x[3:6], 2)
        @test InfiniteOpt._make_ordered_vector(x) == ones(4) * 2
    end
    # test SparseAxisArray
    @testset "SparseAxisArray" begin
        JuMPC.@container(x[i = 3:6], i, container = SparseAxisArray)
        @test InfiniteOpt._make_ordered_vector(x) == [3, 4, 5, 6]
    end
end

# Test allequal
@testset "_allequal" begin
    @test InfiniteOpt._allequal([1])
    @test !InfiniteOpt._allequal([1, 2])
    @test InfiniteOpt._allequal([1, 1])
    @test InfiniteOpt._allequal([1 1; 1 1])
end

# Test _make_float_info
@testset "_make_float_info" begin
    info1 = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    num = Float64(0)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    @test InfiniteOpt._make_float_info(info1) isa VariableInfo{Float64, Float64, Float64, Float64}
    @test InfiniteOpt._make_float_info(info2) == info2
end
