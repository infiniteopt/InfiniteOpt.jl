# Test VectorTuples
@testset "VectorTuple" begin
    # Define useful data objects
    a = 42
    b = [4, -1, 2]
    c = [6 7; 5 1]
    c2 = JuMPC.DenseAxisArray(c, 3:4, [:q, :z])

    @testset "Definition" begin
        # Existence check
        @testset "DataType" begin
            @test IC.VectorTuple isa UnionAll
        end
        # VectorTuple constructor from Tuple
        @testset "Constructor (Tuple)" begin
            # test simple
            @test IC.VectorTuple((a,)) isa IC.VectorTuple
            @test IC.VectorTuple((a,)).values == [a]
            @test IC.VectorTuple((a,)).ranges == [1:1]
            @test IC.VectorTuple((a,)).dimensions == [0]
            @test IC.VectorTuple((a,)).num_columns == [1]
            # test complicated
            tuple = (a, b, c)
            @test IC.VectorTuple(tuple).values == [(tuple...)...]
            @test IC.VectorTuple(tuple).ranges == [1:1, 2:4, 5:8]
            @test IC.VectorTuple(tuple).dimensions == [0, 1, 2]
            @test IC.VectorTuple(tuple).num_columns == [1, 1, 2]
            # test error
            @test_throws ErrorException IC.VectorTuple((c2,))
        end
        # VectorTuple constructor from args...
        @testset "Constructor (items...)" begin
            # test simple
            @test IC.VectorTuple(a) isa IC.VectorTuple
            @test IC.VectorTuple(a).values == [a]
            @test IC.VectorTuple(a).ranges == [1:1]
            @test IC.VectorTuple(a).dimensions == [0]
            @test IC.VectorTuple(a).num_columns == [1]
            # test complicated
            tuple = (a, b, c)
            @test IC.VectorTuple(tuple...).values == [(tuple...)...]
            @test IC.VectorTuple(tuple...).ranges == [1:1, 2:4, 5:8]
            @test IC.VectorTuple(tuple...).dimensions == [0, 1, 2]
            @test IC.VectorTuple(tuple...).num_columns == [1, 1, 2]
        end
        # Base.:(==)
        @testset "Base.:(==)" begin
            @test IC.VectorTuple(a) == IC.VectorTuple(a)
            @test IC.VectorTuple(a, b, c) == IC.VectorTuple(a, b, c)
            @test IC.VectorTuple(a, b) != IC.VectorTuple(b, c)
        end
        # test same_structure
        @testset "same_structure" begin
            vt1 = IC.VectorTuple(a, b, c)
            vt2 = IC.VectorTuple(a, b, b)
            vt3 = IC.VectorTuple('a', ['d', 't', 'q'], ['b' 'a'; 'd' 'f'])
            # test with indices
            @test IC.same_structure(vt1, vt3)
            @test !IC.same_structure(vt1, vt2)
            @test !IC.same_structure(vt2, vt3)
        end
        # VectorTuple copying
        @testset "Base.copy" begin
            vt1 = IC.VectorTuple(a, b, c)
            @test copy(vt1) == vt1
        end
    end

    @testset "Indexing" begin
        vt = IC.VectorTuple(b, a, c)
        values = copy(vt.values)
        # test size
        @testset "Base.size" begin
            @test size(vt, 1) == 3
            @test_throws ArgumentError size(vt, 2)
        end
        # test firstindex
        @testset "Base.firstindex" begin
            @test firstindex(vt) == 1
            @test firstindex(vt, 1) == 1
            @test firstindex(vt, 2) == 1
        end
        # test lastindex
        @testset "Base.lastindex" begin
            @test lastindex(vt) == 8
            @test lastindex(vt, 1) == 3
            @test_throws ArgumentError lastindex(vt, 2)
        end
        # test getindex (linear indexing)
        @testset "Base.getindex (Linear)" begin
            # test scalar indexing
            @test vt[2] == -1
            @test vt[end] == 1
            # test slicing
            @test vt[1:4] == [b; a]
            @test vt[5:end] == [c...]
            @test vt[:] == vt.values
            # test out of bounds
            @test_throws BoundsError vt[1:end+1]
            @test_throws BoundsError vt[13]
        end
        # test getindex (tuple indexing)
        @testset "Base.getindex (Tuple)" begin
            # test scalar
            @test vt[1, 2] == b[2]
            @test vt[2, 1] == a
            @test vt[end, 3] == c[3]
            # test slicing over 2nd dimension
            @test vt[3, :] == [c...]
            @test vt[2, :] == [a]
            @test vt[1, 2:3] == b[2:3]
            @test vt[3, 1:4] == [c...]
            # test slicing over first dimension
            @test vt[:, :] == [b, [a], [c...]]
            @test vt[:, 1] == [first(b), a, first(c)]
            @test vt[2:3, 1] == [a, first(c)]
            @test vt[[1, 3], 2:3] == [b[2:3], c[2:3]]
            # test out of bounds errors
            @test_throws BoundsError vt[:, 2]
            @test_throws BoundsError vt[4, 1]
            @test_throws BoundsError vt[2, 2]
            @test_throws BoundsError vt[2:3, 1:2]
        end
        # test setindex! (linear indexing)
        @testset "Base.setindex! (Linear)" begin
            # test scalar indexing
            @test (vt[2] = 1) == 1
            @test (vt[2] = -1) == -1
            @test (vt[end] = 2) == 2
            @test (vt[end] = 1) == 1
            # test slicing
            @test (vt[1:4] = collect(1:4)) == 1:4
            @test (vt[1:4] = [b; a]) == [b; a]
            @test (vt[6:7] = [3, 4]) == [3, 4]
            @test (vt[6:7] = c[2:3]) == c[2:3]
            # test everything reset
            @test vt.values == values
            # test out of bounds
            @test_throws BoundsError (vt[13] = 42)
            # test size mismatch
            @test_throws DimensionMismatch (vt[1:2] = [1])
        end
        # test setindex! (tuple indexing)
        @testset "Base.setindex! (Tuple)" begin
            # test scalar
            @test (vt[1, 2] = 42) == 42
            @test (vt[1, 2] = b[2]) == b[2]
            @test (vt[end, 3] = 3) == 3
            @test (vt[end, 3] = 9) == 9
            # test slicing over 2nd dimension
            @test (vt[3, :] = [1, 2, 5, 1]) == [1, 2, 5, 1]
            @test (vt[3, :] = [-2, 8, 9, 6]) == [-2, 8, 9, 6]
            @test (vt[1, 2:3] = [2, 6]) == [2, 6]
            @test (vt[1, 2:3] = b[2:3]) == b[2:3]
            # test slicing over first dimension
            @test (vt[:, 1] = [1, 3, 4]) == [1, 3, 4]
            @test (vt[:, 1] = [first(b), a, first(c)]) == [first(b), a, first(c)]
            @test (vt[2:3, 1] = [1, 2]) == [1, 2]
            @test (vt[[1, 3], 2:3] = [[1, 1], [2, 2]]) == [[1, 1], [2, 2]]
            @test (vt[[1, 3], 2:3] = [b[2:3], [8, 9]]) == [b[2:3], [8, 9]]
            @test (vt[2, 1] = a) == a
            @test (vt[3, :] = [c...]) == [c...] 
            # test more colon combos
            vt2 = IC.VectorTuple(b, c)
            values2 = copy(vt2.values)
            @test (vt2[:, 2:3] = [[1, 1], [3, 3]]) == [[1, 1], [3, 3]]
            @test (vt2[:, 2:3] = [b[2:3], [8, 9]]) == [b[2:3], [8, 9]]
            @test (vt2[:, :] = [[1, 1, 1], [3, 3, 3, 3]]) == [[1, 1, 1], [3, 3, 3, 3]]
            @test (vt2[:, :] = [b, [c...]]) == [b, [c...]]
            # test everything reset
            @test vt.values == values
            @test vt2.values == values2
            # test out of bounds errors
            @test_throws BoundsError (vt[2, 2] = 42)
            @test_throws BoundsError (vt[3, 1:5] = [1, 2, 3, 4, 5])
            @test_throws BoundsError (vt[1:2, 1:2] = [[1, 1], [2, 2]])
            # test dimension errors
            @test_throws AssertionError (vt[[1, 3], 1] = 1)
            @test_throws DimensionMismatch (vt[3, 2:3] = [1, 4, 2])
            @test_throws DimensionMismatch (vt[:, 1] = [1, 2])
        end
    end

    @testset "Basic Queries" begin
        vt = IC.VectorTuple(a, b, c)
        # test length
        @testset "Base.length" begin
            @test length(vt) == length(a) + length(b) + length(c)
        end
        # test isempty
        @testset "Base.isempty" begin
            @test !isempty(vt)
            @test isempty(IC.VectorTuple())
        end
        # test eachindex
        @testset "Base.eachindex" begin
            @test eachindex(vt) == eachindex(vt.values)
        end
        # test keys
        @testset "Base.keys" begin
            @test keys(vt) == keys(vt.values)
        end
        # test findfirst
        @testset "Base.findfirst" begin
            @test findfirst(isequal(-1), vt) == 3
        end
        # test findall
        @testset "Base.findall" begin
            @test findall(isequal(42), vt) == [1]
        end
        # test in
        @testset "Base.in" begin
            @test -1 in vt
            @test !(92 in vt)
        end
        # test iterate
        @testset "Base.iterate" begin
            @test [i for i in vt] == vt.values
        end
    end

    @testset "Copy Modification" begin
        # test restricted_copy
        @testset "restricted_copy" begin
            vt = IC.VectorTuple(b, a)
            @test isempty(IC.restricted_copy(vt, [false, false]))
            @test IC.restricted_copy(vt, [false, true]) == IC.VectorTuple(a)
            @test IC.restricted_copy(vt, [true, false]) == IC.VectorTuple(b)
            @test IC.restricted_copy(vt, [true, true]) == IC.VectorTuple(b, a)
            vt = IC.VectorTuple(c, a)
            @test IC.restricted_copy(vt, [true, false]) == IC.VectorTuple(c)
            @test IC.restricted_copy(vt, [false, true]) == IC.VectorTuple(a)
        end
    end

    @testset "Tuple Creation" begin
        # test Tuple
        @testset "Base.Tuple" begin
            vt = IC.VectorTuple(a, c, b)
            @test Tuple(vt) == (a, c, b)
        end
        # test Tuple
        @testset "Base.Tuple (other values)" begin
            vt = IC.VectorTuple(a, c, b)
            values = [:a for i in 1:length(vt)]
            @test Tuple(values, vt) == (:a, [:a :a; :a :a], [:a, :a, :a])
        end
    end

    @testset "Printing Methods" begin
        vt = IC.VectorTuple(a, c)
        @testset "Base.string" begin
            @test string(vt) == "(42, [6 7; 5 1])"
        end
        # test Base.show (REPL)
        @testset "Base.show (REPL)" begin
            show_test(MIME("text/plain"), vt, "(42, [6 7; 5 1])")
        end
        # test Base.show (IJulia)
        @testset "Base.show (IJulia)" begin
            show_test(MIME("text/latex"), vt, "(42, [6 7; 5 1])")
        end
    end
end
