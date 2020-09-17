@testset "DualDict" begin
    # Test the constructors
    @testset "Definition" begin
        # test the default
        @testset "Default" begin
            @test IC.DualDict{String, Int, Symbol}() isa IC.DualDict{String, Int, Symbol}
        end
        # test the tuple constructor
        @testset "Tuple" begin
            @test IC.DualDict{String, Int, Symbol}(("ab" => :a, 3 => :w)) isa IC.DualDict
        end
        # test the splat constructor
        @testset "Pairs" begin
            @test IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w) isa IC.DualDict
        end
    end

    # test basic extensions
    @testset "Basic Extensions" begin
        dd1 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
        # test Base.:(==)
        @testset "Base.:(==)" begin
            @test dd1 == IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
        end
        # test Base.copy
        @testset "Base.copy" begin
            @test copy(dd1) == IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
        end
        # test Base.in
        @testset "Base.in" begin
            @test ("ab" => :a) in dd1
            @test (3 => :w) in dd1
            @test !(("ab" => :b) in dd1)
            @test !((3 => :k) in dd1)
        end
        # test Base.haskey
        @testset "Base.haskey" begin
            @test haskey(dd1, "ab")
            @test haskey(dd1, 3)
            @test !haskey(dd1, "bob")
            @test !haskey(dd1, 42)
        end
        # test Base.length
        @testset "Base.length" begin
            @test length(dd1) == 2
        end
        # test Base.isempty
        @testset "Base.isempty" begin
            @test !isempty(dd1)
            @test isempty(IC.DualDict{Int, Float64, Symbol}())
        end
    end

    # test indexing
    @testset "Indexing" begin
        dd1 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
        # test Base.getindex
        @testset "Base.getindex" begin
            @test dd1["ab"] == :a
            @test dd1[3] == :w
            @test_throws BoundsError dd1["bob"]
            @test_throws BoundsError dd1[54]
        end
        # test Base.setindex!
        @testset "Base.setindex!" begin
            @test (dd1["ab"] = :d) == :d
            @test (dd1[3] = :l) == :l
            @test (dd1["dog"] = :w) == :w
            @test (dd1[42] = :k) == :k
            @test length(dd1) == 4
        end
    end

    # test modification methods
    @testset "Modifcation" begin
        # test Base.empty!
        @testset "Base.empty!" begin
            dd1 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
            @test empty!(dd1) == IC.DualDict{String, Int, Symbol}()
        end
        # test Base.delete!
        @testset "Base.delete!" begin
            dd1 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
            @test delete!(dd1, "ab") == IC.DualDict{String, Int, Symbol}(3 => :w)
            @test delete!(dd1, 3) == IC.DualDict{String, Int, Symbol}()
        end
    end

    # test iteration methods
    @testset "Iteration" begin
        # setup
        dd1 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w)
        dd2 = IC.DualDict{String, Int, Symbol}(3 => :w)
        dd3 = IC.DualDict{String, Int, Symbol}("ab" => :a)
        dd4 = IC.DualDict{String, Int, Symbol}()
        dd5 = IC.DualDict{String, Int, Symbol}("ab" => :a, 3 => :w, "vf" => :b, 42 => :k)
        # test _process_iterate
        @testset "_process_iterate" begin
            @test IC._process_iterate(42, nothing) isa Nothing
            @test IC._process_iterate(42, (2 => 3, 17)) == (2=>3, (42, 17))
        end
        # test Base.iterate
        @testset "Base.iterate" begin
            @test [p for p in dd1] isa Vector
            @test [p for p in dd2] == [3 => :w]
            @test [p for p in dd3] == ["ab" => :a]
            @test [p for p in dd4] == []
            @test [p for p in dd5] isa Vector
         end
         # test Base.keys
         @testset "Base.keys" begin
             @test keys(dd5) isa Vector
         end
         # test Base.values
         @testset "Base.values" begin
             @test values(dd5) isa Vector{Symbol}
         end
    end

    # test Base.show
    @testset "Showing" begin
        dd1 = IC.DualDict{String, Int64, Symbol}()
        dd2 = IC.DualDict{String, Int64, Symbol}(Int64(3) => :w)
        # test show_string
        @testset "show_string" begin
            @test IC.show_string(dd1) == "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 0 entries:"
            @test IC.show_string(dd2) == "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 1 entry:\n  3 => :w"
        end
        # test Base.show (REPL)
        @testset "Base.show (REPL)" begin
            show_test(REPLMode, dd1, "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 0 entries:")
            show_test(REPLMode, dd2, "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 1 entry:\n  3 => :w")
        end
        # test Base.show (IJulia)
        @testset "Base.show (IJulia)" begin
            show_test(IJuliaMode, dd1, "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 0 entries:")
            show_test(IJuliaMode, dd2, "InfiniteOpt.Collections.DualDict{String,Int64,Symbol} with 1 entry:\n  3 => :w")
        end
    end
end
