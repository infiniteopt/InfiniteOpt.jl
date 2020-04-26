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
end

# Test extension of keys
@testset "Base.keys" begin
    m = Model()
    @variable(m, x[1:2], container = SparseAxisArray)
    @test keys(x) == keys(x.data)
end

# Test extension of isapprox
@testset "Base.isapprox" begin
    a = JuMPC.SparseAxisArray(Dict((1,)=>-0, (3,)=>1.1))
    b = JuMPC.SparseAxisArray(Dict((1,)=>-0, (3,)=>1.1 + 1e-30))
    @test isapprox(a, b)
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
    end
    # test arrays
    @testset "Array" begin
        @test InfiniteOpt._make_vector(x) == x
        @test InfiniteOpt._make_vector(ones(2, 2)) == ones(4)
    end
    # test other stuff
    @testset "Non-Array" begin
        @test InfiniteOpt._make_vector(2) == 2
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
