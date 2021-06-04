# Test VectorTuples
@testset "VectorTuple" begin
    # Define useful data objects
    a = 42
    b = [4, -1, 2]
    c = [6 7; 5 1]
    c2 = JuMPC.DenseAxisArray(c, 3:4, [:q, :z])
    d = JuMPC.SparseAxisArray(Dict((2, 1) => 8, (2, 2) => 9, (1, 1) => -2))

    @testset "Definition" begin
        # Existence check
        @testset "DataType" begin
            @test IC.VectorTuple isa UnionAll
        end
        # _make_ordered (SparseAxisArray)
        @testset "_make_ordered (SparseAxisArray)" begin
            @test IC._make_ordered(d, IC.indices(d)) == [-2, 8, 9]
        end
        # _make_ordered (Fallback)
        @testset "_make_ordered (Fallback)" begin
            @test IC._make_ordered(a, nothing) == a
        end
        # VectorTuple constructor from Tuple
        @testset "Constructor (Tuple)" begin
            # test simple
            @test IC.VectorTuple((a,)) isa IC.VectorTuple
            @test IC.VectorTuple((a,)).values == [a]
            @test IC.VectorTuple((a,)).ranges == [1:1]
            @test IC.VectorTuple((a,)).indices == (nothing,)
            # test complicated
            tuple = (a, b, c, c2, d)
            @test sort(IC.VectorTuple(tuple).values) == sort([(tuple...)...])
            @test IC.VectorTuple(tuple).ranges == [1:1, 2:4, 5:8, 9:12, 13:15]
            @test IC.VectorTuple(tuple).indices == IC.indices.(tuple)
        end
        # VectorTuple constructor from args...
        @testset "Constructor (items...)" begin
            # test simple
            @test IC.VectorTuple(a) isa IC.VectorTuple
            @test IC.VectorTuple(a).values == [a]
            @test IC.VectorTuple(a).ranges == [1:1]
            @test IC.VectorTuple(a).indices == (nothing,)
            # test complicated
            tuple = (a, b, c, c2, d)
            @test sort(IC.VectorTuple(tuple...).values) == sort([(tuple...)...])
            @test IC.VectorTuple(tuple...).ranges == [1:1, 2:4, 5:8, 9:12, 13:15]
            @test IC.VectorTuple(tuple...).indices == IC.indices.(tuple)
        end
        # Base.:(==)
        @testset "Base.:(==)" begin
            @test IC.VectorTuple(a) == IC.VectorTuple(a)
            @test IC.VectorTuple(a, b, c, c2, d) == IC.VectorTuple(a, b, c, c2, d)
            @test IC.VectorTuple(a, b) != IC.VectorTuple(b, c)
        end
        # test same_structure
        @testset "same_structure" begin
            vt1 = IC.VectorTuple(a, b, c2)
            vt2 = IC.VectorTuple(a, b, c)
            vt3 = IC.VectorTuple(a, b, d)
            vt4 = IC.VectorTuple(a, b, b)
            vt5 = IC.VectorTuple('a', ['d', 't', 'q'], ['b' 'a'; 'd' 'f'])
            # test with indices
            @test IC.same_structure(vt2, vt5)
            @test !IC.same_structure(vt1, vt2)
            @test !IC.same_structure(vt2, vt3)
            @test !IC.same_structure(vt3, vt4)
            @test !IC.same_structure(vt1, vt5)
            # test without indices
            @test IC.same_structure(vt2, vt5, use_indices = false)
            @test IC.same_structure(vt1, vt2, use_indices = false)
            @test IC.same_structure(vt3, vt4, use_indices = false)
            @test !IC.same_structure(vt2, vt3, use_indices = false)
            @test !IC.same_structure(vt4, vt5, use_indices = false)
            @test IC.same_structure(vt1, vt5, use_indices = false)
        end
        # VectorTuple copying
        @testset "Base.copy" begin
            vt1 = IC.VectorTuple(a, b, c2)
            @test copy(vt1) == vt1
        end
    end

    @testset "Indexing" begin
        vt = IC.VectorTuple(b, a, c2, d)
        values = copy(vt.values)
        # test size
        @testset "Base.size" begin
            @test size(vt, 1) == 4
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
            @test lastindex(vt) == 11
            @test lastindex(vt, 1) == 4
            @test_throws ArgumentError lastindex(vt, 2)
        end
        # test getindex (linear indexing)
        @testset "Base.getindex (Linear)" begin
            # test scalar indexing
            @test vt[2] == -1
            @test vt[end] == 9
            # test slicing
            @test vt[1:4] == [b; a]
            @test vt[5:end] == [c2..., -2, 8, 9]
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
            @test vt[end, 3] == 9
            # test slicing over 2nd dimension
            @test vt[4, :] == [-2, 8, 9]
            @test vt[2, :] == [a]
            @test vt[1, 2:3] == b[2:3]
            @test vt[3, 1:4] == [c2...]
            # test slicing over first dimension
            @test vt[:, :] == [b, [a], [c2...], [-2, 8, 9]]
            @test vt[:, 1] == [first(b), a, first(c2), -2]
            @test vt[3:4, 2:3] == [c[2:3], [8, 9]]
            @test vt[[1, 4], 2:3] == [b[2:3], [8, 9]]
            # test out of bounds errors
            @test_throws BoundsError vt[:, 2]
            @test_throws BoundsError vt[5, 1]
            @test_throws BoundsError vt[2, 2]
            @test_throws BoundsError vt[3:4, 2:4]
        end
        # test setindex! (linear indexing)
        @testset "Base.setindex! (Linear)" begin
            # test scalar indexing
            @test (vt[2] = 1) == 1
            @test (vt[2] = -1) == -1
            @test (vt[end] = 2) == 2
            @test (vt[end] = 9) == 9
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
            @test (vt[4, :] = [1, 2, 5]) == [1, 2, 5]
            @test (vt[4, :] = [-2, 8, 9]) == [-2, 8, 9]
            @test (vt[1, 2:3] = [2, 6]) == [2, 6]
            @test (vt[1, 2:3] = b[2:3]) == b[2:3]
            # test slicing over first dimension
            @test (vt[:, 1] = [1, 3, 4, 5]) == [1, 3, 4, 5]
            @test (vt[:, 1] = [first(b), a, first(c2), -2]) == [first(b), a, first(c2), -2]
            @test (vt[3:4, 2:3] = [[1, 1], [2, 2]]) == [[1, 1], [2, 2]]
            @test (vt[3:4, 2:3] = [c[2:3], [8, 9]]) == [c[2:3], [8, 9]]
            @test (vt[[1, 4], 2:3] = [[1, 1], [2, 2]]) == [[1, 1], [2, 2]]
            @test (vt[[1, 4], 2:3] = [b[2:3], [8, 9]]) == [b[2:3], [8, 9]]
            # test more colon combos
            vt2 = IC.VectorTuple(b, d)
            values2 = copy(vt2.values)
            @test (vt2[:, 2:3] = [[1, 1], [3, 3]]) == [[1, 1], [3, 3]]
            @test (vt2[:, 2:3] = [b[2:3], [8, 9]]) == [b[2:3], [8, 9]]
            @test (vt2[:, :] = [[1, 1, 1], [3, 3, 3]]) == [[1, 1, 1], [3, 3, 3]]
            @test (vt2[:, :] = [b, [-2, 8, 9]]) == [b, [-2, 8, 9]]
            # test everything reset
            @test vt.values == values
            @test vt2.values == values2
            # test out of bounds errors
            @test_throws BoundsError (vt[2, 2] = 42)
            @test_throws BoundsError (vt[4, 1:4] = [1, 2, 3, 4])
            @test_throws BoundsError (vt[1:2, 1:2] = [[1, 1], [2, 2]])
            # test dimension errors
            @test_throws AssertionError (vt[[1, 3], 1] = 1)
            @test_throws DimensionMismatch (vt[3, 2:3] = [1, 4, 2])
            @test_throws DimensionMismatch (vt[:, 1] = [1, 4, 2])
        end
    end

    @testset "Basic Queries" begin
        vt = IC.VectorTuple(a, b, c2)
        # test length
        @testset "Base.length" begin
            @test length(vt) == length(a) + length(b) + length(c2)
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
        # test _update_indices (SparseAxisArray)
        @testset "_update_indices (SparseAxisArray)" begin
            inds = IC.indices(d)
            new_axes = [(1, 1), (2, 2)]
            @test IC._update_indices(inds, 2).axes == new_axes
            inds = IC.ContainerIndices(CartesianIndices(new_axes), new_axes)
            @test IC._update_indices(inds, 2) isa Nothing
        end
        # test _update_indices (DenseAxisArray)
        @testset "_update_indices (DenseAxisArray)" begin
            inds = IC.indices(c2)
            @test IC._update_indices(inds, 2).axes == [(3, :q), (3, :z), (4, :z)]
        end
        # test _update_indices (Array)
        @testset "_update_indices (Array)" begin
            inds = IC.indices(c)
            @test IC._update_indices(inds, 2).axes == [(1, 1), (1, 2), (2, 2)]
        end
        # test restricted_copy
        @testset "restricted_copy" begin
            # test with linear indexing with simple results
            vt = IC.VectorTuple(b, a)
            @test isempty(IC.restricted_copy(vt, [false, false, false, false]))
            vt = IC.VectorTuple(b, a)
            @test IC.restricted_copy(vt, [false, false, false, true]) == IC.VectorTuple(a)
            vt = IC.VectorTuple(c2, a)
            @test IC.restricted_copy(vt, [true, true, true, true, false]) == IC.VectorTuple(c2)
            # test with linear indexing and array reduction
            vt = IC.VectorTuple(d, a, c)
            new_d = JuMPC.SparseAxisArray(delete!(copy(d.data), (1, 1)))
            inds = CartesianIndices(c)
            cs = JuMPC.SparseAxisArray(Dict(Tuple(inds[k]) => c[inds[k]] for k = [1, 4]))
            copy_inds = [false, true, true, true, true, false, false, true]
            @test IC.restricted_copy(vt, copy_inds) == IC.VectorTuple(new_d, a, cs)
            @test vt == IC.VectorTuple(d, a, c)
            # test with linear indexing and array reduction again
            vt = IC.VectorTuple(a, c, d)
            new_d = JuMPC.SparseAxisArray(delete!(copy(d.data), (2, 1)))
            inds = CartesianIndices(c)
            cs = JuMPC.SparseAxisArray(Dict(Tuple(inds[k]) => c[inds[k]] for k = [2, 4]))
            copy_inds = [false, false, true, false, true, true, false, true]
            @test IC.restricted_copy(vt, copy_inds) == IC.VectorTuple(cs, new_d)
            @test vt == IC.VectorTuple(a, c, d)
        end
    end

    @testset "Tuple Creation" begin
        # test _make_array (Array)
        @testset "_make_array (Array)" begin
            vt = IC.VectorTuple(b, a, c)
            @test IC._make_array(vt.values[1:3], vt.indices[1]) == b
            @test IC._make_array(vt.values[5:8], vt.indices[3]) == c
        end
        # test _make_array (DenseAxisArray)
        @testset "_make_array (DenseAxisArray)" begin
            vt = IC.VectorTuple(c2, a)
            @test IC._make_array(vt.values[1:4], vt.indices[1]) == c2
        end
        # test _make_array (SparseAxisArray)
        @testset "_make_array (SparseAxisArray)" begin
            vt = IC.VectorTuple(a, d)
            @test IC._make_array(vt.values[2:4], vt.indices[2]) == d
        end
        # test _make_array (Fallback)
        @testset "_make_array (Fallback)" begin
            vt = IC.VectorTuple(a, d)
            @test IC._make_array(vt.values[1:1], vt.indices[1]) == a
            vt = IC.VectorTuple(b, (1, 2))
            @test IC._make_array(vt.values[4:5], vt.indices[2]) == [1, 2]
        end
        # test Tuple
        @testset "Base.Tuple" begin
            # test reproducing original
            vt = IC.VectorTuple(a, d, c2, b)
            @test Tuple(vt) == (a, d, c2, b)
            # test vector mode
            vt = IC.VectorTuple(a, c)
            @test Tuple(vt, use_indices = false) == (a, [c...])
        end
        # test Tuple
        @testset "Base.Tuple (other values)" begin
            # test reproducing original
            vt = IC.VectorTuple(a, c, b)
            values = [:a for i in 1:length(vt)]
            @test Tuple(values, vt) == (:a, [:a :a; :a :a], [:a, :a, :a])
            # test vector mode
            vt = IC.VectorTuple(a, c)
            values = [:a for i in 1:length(vt)]
            @test Tuple(values, vt, use_indices = false) == (:a, [:a, :a, :a, :a])
        end
    end

    @testset "Printing Methods" begin
        vt = IC.VectorTuple(a, c)
        @testset "Base.string" begin
            @test string(vt) == "(42, [6, 5, 7, 1])"
        end
        # test Base.show (REPL)
        @testset "Base.show (REPL)" begin
            show_test(REPLMode, vt, "(42, [6, 5, 7, 1])")
        end
        # test Base.show (IJulia)
        @testset "Base.show (IJulia)" begin
            show_test(IJuliaMode, vt, "(42, [6, 5, 7, 1])")
        end
    end
end
