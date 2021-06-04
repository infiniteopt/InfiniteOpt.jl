@testset "Vectorizing Methods" begin 
    # Define useful data objects
    a = 42
    b = [4, -1, 2]
    b2 = JuMPC.DenseAxisArray(b, 3:4)
    c = [6 7; 5 1]
    c2 = JuMPC.DenseAxisArray(c, 3:4, [:q, :z])
    d = JuMPC.SparseAxisArray(Dict((2, 1) => 8, (2, 2) => 9, (1, 1) => -2))
    d_inds = [(1, 1), (2, 1), (2, 2)]

    # test ContainerIndices
    @testset "ContainerIndices" begin 
        # Julia arrays
        @test IC.ContainerIndices(b) isa IC.ContainerIndices{1, Nothing}
        @test IC.ContainerIndices(b) == IC.ContainerIndices(CartesianIndices(b), nothing)
        @test IC.ContainerIndices(c) isa IC.ContainerIndices{2, Nothing}
        @test IC.ContainerIndices(c) == IC.ContainerIndices(CartesianIndices(c), nothing)
        # DenseAxisArrays 
        @test IC.ContainerIndices(b2) isa IC.ContainerIndices{1, Tuple{UnitRange{Int64}}}
        @test IC.ContainerIndices(b2) == IC.ContainerIndices(CartesianIndices(b), (3:4,))
        @test IC.ContainerIndices(c2) isa IC.ContainerIndices{2, Tuple{UnitRange{Int64}, Array{Symbol,1}}}
        @test IC.ContainerIndices(c2) == IC.ContainerIndices(CartesianIndices(c), (3:4, [:q, :z]))
        # SparseAxisArrays
        @test IC.ContainerIndices(d) isa IC.ContainerIndices{1, Vector{Tuple{Int, Int}}}
        @test IC.ContainerIndices(d) == IC.ContainerIndices(CartesianIndices(1:3), d_inds)
    end

    # test basic extensions
    @testset "Base Extensions" begin 
        # ==
        @test IC.ContainerIndices(b) == IC.ContainerIndices(ones(3))
        @test IC.ContainerIndices(d) == IC.ContainerIndices(d)
        @test IC.ContainerIndices(c) != IC.ContainerIndices(c2)
        # size 
        @test size(IC.ContainerIndices(b)) == (3,)
        @test size(IC.ContainerIndices(c)) == (2,2)
        @test size(IC.ContainerIndices(c2)) == (2,2)
        # length
        @test length(IC.ContainerIndices(b)) == 3
        @test length(IC.ContainerIndices(c)) == 4
        @test length(IC.ContainerIndices(c2)) == 4
        @test length(IC.ContainerIndices(d)) == 3
    end

    # test indices
    @testset "indices" begin 
        @test IC.indices(a) isa Nothing 
        @test IC.indices(b) == IC.ContainerIndices(CartesianIndices(b), nothing)
        @test IC.indices(b2) == IC.ContainerIndices(CartesianIndices(b), (3:4,))
        @test IC.indices(c) == IC.ContainerIndices(CartesianIndices(c), nothing)
        @test IC.indices(c2) == IC.ContainerIndices(CartesianIndices(c), (3:4, [:q, :z]))
        @test IC.indices(d) == IC.ContainerIndices(CartesianIndices(1:3), d_inds)
    end

    # test vectorize
    @testset "vectorize" begin 
        # 1 argument
        @test IC.vectorize(a) == a
        @test IC.vectorize(b) == b
        @test IC.vectorize(b2) == b
        @test IC.vectorize(c) == [6, 5, 7, 1]
        @test IC.vectorize(c2) == [6, 5, 7, 1]
        @test IC.vectorize(d) == [-2, 8, 9]
        @test IC.vectorize([]) == []
        # 2 arguments 
        @test IC.vectorize(a, IC.indices(a)) == a
        @test IC.vectorize(b, IC.indices(b)) == b
        @test IC.vectorize(b2, IC.indices(b2)) == b
        @test IC.vectorize(c, IC.indices(c)) == [6, 5, 7, 1]
        @test IC.vectorize(c2, IC.indices(c2)) == [6, 5, 7, 1]
        @test IC.vectorize(d, IC.indices(d)) == [-2, 8, 9]
        @test IC.vectorize([], IC.indices([])) == []
    end

    # test unvectorize
    @testset "unvectorize" begin 
        @test IC.unvectorize(IC.vectorize(a), IC.indices(a)) == a
        @test IC.unvectorize(IC.vectorize(b), IC.indices(b)) == b
        @test IC.unvectorize(IC.vectorize(b2), IC.indices(b2)) == b2
        @test IC.unvectorize(IC.vectorize(c), IC.indices(c)) == c
        @test IC.unvectorize(IC.vectorize(c2), IC.indices(c2)) == c2
        @test IC.unvectorize(IC.vectorize(d), IC.indices(d)) == d
        @test IC.unvectorize(IC.vectorize([]), IC.indices([])) == []
    end
end
