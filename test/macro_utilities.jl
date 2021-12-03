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
    # test _get_name
    @testset "_get_name" begin 
        @test InfiniteOpt._get_name(:a) == :a
        @test InfiniteOpt._get_name(nothing) == ()
        @test InfiniteOpt._get_name("bob") == "bob"
        @test InfiniteOpt._get_name(:("a$(x)")) == :("a$(x)")
        @test InfiniteOpt._get_name(:(x[1:2])) == :x 
    end
    # test _name_call 
    @testset "_name_call" begin 
        @test InfiniteOpt._name_call("bob", []) == "bob"
        @test InfiniteOpt._name_call("", []) == ""
        expected = :(string("x", "[", string($(Expr(:escape, "2"))), ",", string($(Expr(:escape, "1"))), "]"))
        @test InfiniteOpt._name_call("x", ["2", "1"]) == expected
    end
    # test _extract_kwargs
    @testset "_extract_kwargs" begin 
        # test normal
        args = (:x, :(d=42), :(container=Cat))
        @test InfiniteOpt._extract_kwargs(args)[1] == [:x]
        @test InfiniteOpt._extract_kwargs(args)[2] == [:(d=42)]
        @test InfiniteOpt._extract_kwargs(args)[3] == :(Cat)
        @test InfiniteOpt._extract_kwargs(args)[4] === nothing
        # test with semi-colon keywords
        args = (Expr(:parameters, :(d=42), :(base_name = "cat")), :x, :(container=Cat))
        @test InfiniteOpt._extract_kwargs(args)[1] == [:x]
        @test InfiniteOpt._extract_kwargs(args)[2] == [:(d=42)]
        @test InfiniteOpt._extract_kwargs(args)[3] == :(Cat)
        @test InfiniteOpt._extract_kwargs(args)[4] == esc("cat")
    end
    # test _add_kw_args
    @testset "_add_kwargs" begin
        esc_kw = esc(Expr(:kw, :d, 2))
        call = :(f(a))
        @test InfiniteOpt._add_kwargs(call, [:(d = 2)]) isa Nothing
        @test call == :(f(a, $esc_kw))
    end
     # test _has_idxvars
     @testset "_has_idxvars" begin
        vars = Any[:i, :j]
        @test !InfiniteOpt._has_idxvars(:w, vars)
        @test InfiniteOpt._has_idxvars(:i, vars)
        @test InfiniteOpt._has_idxvars(:(a[i]), vars)
        @test !InfiniteOpt._has_idxvars(:(ones(a)), vars)
    end
    # test _valid_model
    @testset "_valid_model" begin
        m = InfiniteModel()
        m2 = 2
        @test InfiniteOpt._valid_model(error, m, :m) isa Nothing
        @test_throws ErrorException InfiniteOpt._valid_model(error, m2, :m2)
    end
    # test _error_if_cannot_register
    @testset "_error_if_cannot_register" begin
        m = InfiniteModel()
        @test InfiniteOpt._error_if_cannot_register(error, m, :x) isa Nothing
        object_dictionary(m)[:x] = 2
        @test_throws ErrorException InfiniteOpt._error_if_cannot_register(error, m, :x)
        @test_throws ErrorException InfiniteOpt._error_if_cannot_register(error, m, :((x, 2)))
    end
    # test _macro_assign_and_return
    @testset "_macro_assign_and_return" begin
        m = InfiniteModel()
        esc_m = esc(:m)
        code = :(GeneralVariableRef($esc_m, 1, FiniteVariableIndex))
        @test InfiniteOpt._macro_assign_and_return(error, code, :a, esc_m) isa Expr
    end
    # test _finalize_macro
    @testset "_finalize_macro" begin
        m = InfiniteModel()
        esc_m = esc(:m)
        code = :(GeneralVariableRef($esc_m, 1, FiniteVariableIndex))
        source = LineNumberNode(42, :bob)
        @test InfiniteOpt._finalize_macro(error, esc_m, code, source) isa Expr
    end
end
