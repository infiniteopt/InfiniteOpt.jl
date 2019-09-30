# Test the extensions of Base.convert
@testset "Convert Extensions" begin
    m = Model()
    @variable(m, x[1:2])
    # GenericAffExpr -> GenericAffExpr
    @testset "Aff->Aff" begin
        aff = x[1] - 3 * x[2] - 2
        @test typeof(convert(GenericAffExpr{Float64, AbstractVariableRef},
                           aff)) == GenericAffExpr{Float64, AbstractVariableRef}
        @test convert(GenericAffExpr{Float64, AbstractVariableRef},
                      aff).terms == aff.terms
        @test convert(GenericAffExpr{Float64, AbstractVariableRef},
                      aff).constant == aff.constant
    end
    # GenericQuadExpr -> GenericQuadExpr
    @testset "Quad->Quad" begin
        quad = 2 * x[1] * x[2] + x[1] - 1
        @test typeof(convert(GenericQuadExpr{Float64, AbstractVariableRef},
                         quad)) == GenericQuadExpr{Float64, AbstractVariableRef}
        @test convert(GenericQuadExpr{Float64, AbstractVariableRef},
                      quad).aff.terms == quad.aff.terms
        @test convert(GenericQuadExpr{Float64, AbstractVariableRef},
                      quad).aff.constant == quad.aff.constant
        @test convert(GenericQuadExpr{Float64, AbstractVariableRef},
                      quad).terms == quad.terms
    end
    # UnorderedPair -> UnorderedPair
    @testset "UoPair->UoPair" begin
        pair = UnorderedPair(x[1], x[2])
        @test typeof(convert(UnorderedPair{AbstractVariableRef},
                                    pair)) == UnorderedPair{AbstractVariableRef}
        @test convert(UnorderedPair{AbstractVariableRef}, pair).a == pair.a
        @test convert(UnorderedPair{AbstractVariableRef}, pair).b == pair.b
    end
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

# Test the possible_convert functions
@testset "_possible_convert" begin
    m = Model()
    @variable(m, x[1:2])
    # GenericAffExpr conversions
    @testset "GenericAffExpr" begin
        aff = x[1] - 3 * x[2] - 2
        @test typeof(InfiniteOpt._possible_convert(AbstractVariableRef,
                           aff)) == GenericAffExpr{Float64, AbstractVariableRef}
        @test InfiniteOpt._possible_convert(AbstractVariableRef,
                                            aff).terms == aff.terms
        @test InfiniteOpt._possible_convert(AbstractVariableRef,
                                            aff).constant == aff.constant
        @test typeof(InfiniteOpt._possible_convert(GeneralVariableRef,
                                   aff)) == GenericAffExpr{Float64, VariableRef}
    end
    # GenericQuadExpr conversions
    @testset "GenericQuadExpr" begin
        quad = 2 * x[1] * x[2] + x[1] - 1
        @test typeof(InfiniteOpt._possible_convert(AbstractVariableRef,
                         quad)) == GenericQuadExpr{Float64, AbstractVariableRef}
        @test InfiniteOpt._possible_convert(AbstractVariableRef,
                                            quad).aff.terms == quad.aff.terms
        @test InfiniteOpt._possible_convert(AbstractVariableRef,
                                         quad).aff.constant == quad.aff.constant
        @test InfiniteOpt._possible_convert(AbstractVariableRef,
                                            quad).terms == quad.terms
        @test typeof(InfiniteOpt._possible_convert(GeneralVariableRef,
                                 quad)) == GenericQuadExpr{Float64, VariableRef}
    end
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
        @test InfiniteOpt._make_vector(ones(2, 2)) == ones(2, 2)
    end
    # test other stuff
    @testset "Non-Array" begin
        @test InfiniteOpt._make_vector(2) == 2
    end
end
