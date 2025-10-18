@testset "Vectorizing Methods" begin 
    # Define useful data objects
    a = 42
    b = [4, -1, 2]
    b2 = JuMPC.DenseAxisArray(b, 3:4)
    c = [6 7; 5 1]
    c2 = JuMPC.DenseAxisArray(c, 3:4, [:q, :z])
    d = JuMPC.SparseAxisArray(Dict((2, 1) => 8, (2, 2) => 9, (1, 1) => -2))
    d_inds = [(1, 1), (2, 1), (2, 2)]
    # test vectorize
    @testset "vectorize" begin 
        @test IC.vectorize(a) == (a, (0, 1))
        @test IC.vectorize(b) == (b, (1, 1))
        @test_throws ErrorException IC.vectorize(b2)
        @test IC.vectorize(c) == ([6, 5, 7, 1], (2, 2))
        @test_throws ErrorException IC.vectorize(c2)
        @test_throws ErrorException IC.vectorize(d)
        @test IC.vectorize([]) == ([], (1, 1))
        @test IC.vectorize(zeros(2, 0)) == ([], (2, 0))
    end
    # test unvectorize
    @testset "unvectorize" begin 
        @test IC.unvectorize(IC.vectorize(a)...) == a
        @test IC.unvectorize(IC.vectorize(b)...) == b
        @test IC.unvectorize(IC.vectorize(c)...) == c
        @test IC.unvectorize(IC.vectorize([])...) == []
    end
end
