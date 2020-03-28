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
            @test VectorTuple isa UnionAll
        end
        # _get_indices (Array)
        @testset "_get_indices (Array)" begin
            @test IC._get_indices(b) == CartesianIndices(b)
            @test IC._get_indices(c) == CartesianIndices(c)
        end
        # _get_indices (DenseAxisArray)
        @testset "_get_indices (DenseAxisArray)" begin
            @test IC._get_indices(c2) == (CartesianIndices(c2), axes(c2))
        end
        # _get_indices (SparseAxisArray)
        @testset "_get_indices (SparseAxisArray)" begin
            @test IC._get_indices(d) == [(1, 1), (2, 1), (2, 2)]
        end
        # _get_indices (Fallback)
        @testset "_get_indices (Fallback)" begin
            @test IC._get_indices(a) isa Nothing
        end
        # _make_ordered (SparseAxisArray)
        @testset "_make_ordered (SparseAxisArray)" begin
            @test IC._make_ordered(d, IC._get_indices(d)) == [-2, 8, 9]
        end
        # _make_ordered (Fallback)
        @testset "_make_ordered (Fallback)" begin
            @test IC._make_ordered(a, keys(d)) == a
        end
        # VectorTuple constructor from Tuple
        @testset "Constructor (Tuple)" begin
            # test empty
            @test VectorTuple(()) isa VectorTuple
            @test VectorTuple(()).values == []
            # test simple
            @test VectorTuple((a,)) isa VectorTuple
            @test VectorTuple((a,)).values == [a]
            @test VectorTuple((a,)).starts == [1]
            @test VectorTuple((a,)).ends == [1]
            @test VectorTuple((a,)).indices == [nothing]
            # test complicated
            tuple = (a, b, c, c2, d)
            @test sort(VectorTuple(tuple).values) == sort([(tuple...)...])
            @test VectorTuple(tuple).starts == [1, 2, 5, 9, 13]
            @test VectorTuple(tuple).ends == [1, 4, 8, 12, 15]
            expected = [nothing, CartesianIndices(b), CartesianIndices(c),
                        (CartesianIndices(c2), axes(c2)), [(1, 1), (2, 1), (2, 2)]]
            @test VectorTuple(tuple).indices == expected
        end
        # VectorTuple constructor from args...
        @testset "Constructor (items...)" begin
            # test simple
            @test VectorTuple(a) isa VectorTuple
            @test VectorTuple(a).values == [a]
            @test VectorTuple(a).starts == [1]
            @test VectorTuple(a).ends == [1]
            @test VectorTuple(a).indices == [nothing]
            # test complicated
            tuple = (a, b, c, c2, d)
            @test sort(VectorTuple(tuple...).values) == sort([(tuple...)...])
            @test VectorTuple(tuple...).starts == [1, 2, 5, 9, 13]
            @test VectorTuple(tuple...).ends == [1, 4, 8, 12, 15]
            expected = [nothing, CartesianIndices(b), CartesianIndices(c),
                        (CartesianIndices(c2), axes(c2)), [(1, 1), (2, 1), (2, 2)]]
            @test VectorTuple(tuple...).indices == expected
        end
        # Base.:(==)
        @testset "Base.:(==)" begin
            @test VectorTuple(a) == VectorTuple(a)
            @test VectorTuple(a, b, c, c2, d) == VectorTuple(a, b, c, c2, d)
            @test VectorTuple(a, b) != VectorTuple(b, c)
        end
        # VectorTuple constructor default
        @testset "Constructor (Default)" begin
            @test VectorTuple() == VectorTuple(Any[], Int64[], Int64[], Any[])
        end
    end

    @testset "Indexing" begin
        vt = VectorTuple(b, a, c2, d)
        values = copy(vt.values)
        # test tuple_length
        @testset "tuple_length" begin
            @test tuple_length(vt) == 4
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
            @test_throws ErrorException lastindex(vt, 2)
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
            vt2 = VectorTuple(b, d)
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
        vt = VectorTuple(a, b, c2)
        # test length
        @testset "Base.length" begin
            @test length(vt) == length(a) + length(b) + length(c2)
        end
        # test isempty
        @testset "Base.isempty" begin
            @test !isempty(vt)
            @test isempty(VectorTuple())
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
    end

    @testset "Modification" begin
        # test empty!
        @testset "Base.empty!" begin
            vt = VectorTuple(b, a, c2)
            @test empty!(vt) == VectorTuple()
            @test isempty(vt)
            @test isempty(vt.starts)
            @test isempty(vt.ends)
            @test isempty(vt.indices)
        end
        # test push! (single item)
        @testset "Base.push! (Single)" begin
            vt = VectorTuple(b, a, c2)
            @test push!(vt, 42) == VectorTuple(b, a, c2, 42)
        end
        # test push! (array)
        @testset "Base.push! (Array)" begin
            vt = VectorTuple(b, a, c2)
            @test push!(vt, c) == VectorTuple(b, a, c2, c)
        end
        # test push! (Multiple)
        @testset "Base.push! (Multiple)" begin
            vt = VectorTuple(b, a, c2)
            @test push!(vt, 42, d) == VectorTuple(b, a, c2, 42, d)
        end
        # test _update_indices (SparseAxisArray)
        @testset "_update_indices (SparseAxisArray)" begin
            inds = IC._get_indices(d)
            @test IC._update_indices(inds, 2) == [(1, 1), (2, 2)]
        end
        # test _update_indices (DenseAxisArray)
        @testset "_update_indices (DenseAxisArray)" begin
            inds = IC._get_indices(c2)
            @test IC._update_indices(inds, 2) == [(3, :q), (3, :z), (4, :z)]
        end
        # test _update_indices (Array)
        @testset "_update_indices (Array)" begin
            inds = IC._get_indices(c)
            @test IC._update_indices(inds, 2) == [(1, 1), (1, 2), (2, 2)]
        end
        # test delete!
        @testset "Base.delete!" begin
            # test with tuple index
            vt = VectorTuple(b, a, c2)
            @test deleteat!(vt, 2, tuple_index = true) == VectorTuple(b, c2)
            @test deleteat!(vt, 2, tuple_index = true) == VectorTuple(b)
            @test isempty(deleteat!(vt, 1, tuple_index = true))
            vt = VectorTuple(d, a, c2)
            @test deleteat!(vt, 3, tuple_index = true) == VectorTuple(d, a)
            @test deleteat!(vt, 1, tuple_index = true) == VectorTuple(a)
            @test isempty(deleteat!(vt, 1, tuple_index = true))
            # test with linear index and single element group
            vt = VectorTuple(b, a, c2)
            @test deleteat!(vt, 4) == VectorTuple(b, c2)
            vt = VectorTuple(a, b, c2)
            @test deleteat!(vt, 1) == VectorTuple(b, c2)
            vt = VectorTuple(b, c2, a)
            @test deleteat!(vt, 8) == VectorTuple(b, c2)
            # test with linear index and truncated Array
            vt = VectorTuple(b, a, c)
            inds = CartesianIndices(c)
            cs = JuMPC.SparseAxisArray(Dict(Tuple(inds[k]) => c[inds[k]] for k = [1, 3, 4]))
            @test deleteat!(vt, 6) == VectorTuple(b, a, cs)
            # test with linear index and truncated DenseAxisArray
            vt = VectorTuple(b, c2, a)
            inds = CartesianIndices(c2)
            new_inds = [k.I for k in JuMPC.DenseAxisArrayKeys(Base.Iterators.product(axes(c2)...))]
            cs = JuMPC.SparseAxisArray(Dict(new_inds[k] => c[inds[k]] for k = [1, 3, 4]))
            @test deleteat!(vt, 5) == VectorTuple(b, cs, a)
            # test with linear index and truncated DenseAxisArray
            vt = VectorTuple(d, c2, a)
            new_d = JuMPC.SparseAxisArray(delete!(copy(d.data), (1, 1)))
            @test deleteat!(vt, 1) == VectorTuple(new_d, c2, a)
        end
    end

    @testset "Tuple Creation" begin
        # test _make_array (Array)
        @testset "_make_array (Array)" begin
            vt = VectorTuple(b, a, c)
            @test IC._make_array(vt.values, 1, 3, vt.indices[1]) == b
            @test IC._make_array(vt.values, 5, 8, vt.indices[3]) == c
        end
        # test _make_array (DenseAxisArray)
        @testset "_make_array (DenseAxisArray)" begin
            vt = VectorTuple(c2, a)
            @test IC._make_array(vt.values, 1, 4, vt.indices[1]) == c2
        end
        # test _make_array (SparseAxisArray)
        @testset "_make_array (SparseAxisArray)" begin
            vt = VectorTuple(a, d)
            @test IC._make_array(vt.values, 2, 4, vt.indices[2]) == d
        end
        # test _make_array (Fallback)
        @testset "_make_array (Fallback)" begin
            vt = VectorTuple(a, d)
            @test IC._make_array(vt.values, 1, 1, vt.indices[1]) == a
            vt = VectorTuple(b, (1, 2))
            @test IC._make_array(vt.values, 4, 5, vt.indices[2]) == [1, 2]
        end
        # test Tuple
        @testset "Base.Tuple" begin
            # test reproducing original
            vt = VectorTuple(a, d, c2, b)
            @test Tuple(vt) == (a, d, c2, b)
            # test vector mode
            vt = VectorTuple(a, c)
            @test Tuple(vt, use_indices = false) == (a, [c...])
        end
    end

    @testset "Printing Methods" begin
        vt = VectorTuple(a, c)
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
