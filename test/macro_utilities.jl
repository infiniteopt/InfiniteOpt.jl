# Test basic methods 
@testset "Basics" begin 
    # test _macro_error 
    @testset "_macro_error" begin 
        source = LineNumberNode(42, :bob)
        @test_throws ErrorException InfiniteOpt._macro_error(:test, (:x,), source, "hello")
    end
    # test _esc_non_constant
    @testset "_esc_non_constant" begin 
        @test InfiniteOpt._esc_non_constant(42) == 42
        @test InfiniteOpt._esc_non_constant(:(quote 42 end)) == :(quote 42 end)
        @test InfiniteOpt._esc_non_constant(:(x[i])) == esc(:(x[i]))
        @test InfiniteOpt._esc_non_constant(:x) == esc(:x)
    end
    # test _name_call 
    @testset "_name_call" begin 
        @test InfiniteOpt._name_call("bob", []) == "bob"
        @test InfiniteOpt._name_call("", []) == ""
        expected = :(string("x", "[", string($(Expr(:escape, "2"))), ",", string($(Expr(:escape, "1"))), "]"))
        @test InfiniteOpt._name_call("x", ["2", "1"]) == expected
    end
    # test _extract_kw_args
    @testset "_extract_kw_args" begin 
        args = (:x, :(d=42), :(container=Cat))
        @test InfiniteOpt._extract_kw_args(args)[1] == [:x]
        @test InfiniteOpt._extract_kw_args(args)[2] == [:(d=42)]
        @test InfiniteOpt._extract_kw_args(args)[3] == :(Cat)
    end
end
