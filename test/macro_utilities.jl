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

# Test container methods
@testset "Container Helpers" begin 
    # test _parse_idxset
    @testset "_parse_idxset" begin
        @test InfiniteOpt._parse_idxset(:(i = 1:2)) == (:i, esc(:(1:2)))
        @test InfiniteOpt._parse_idxset(:(i in 1:2)) == (:i, esc(:(1:2)))
        @test InfiniteOpt._parse_idxset(:(1:2))[1] isa Symbol
        @test InfiniteOpt._parse_idxset(:(1:2))[2] == esc(:(1:2))
        @test InfiniteOpt._parse_idxset(:a)[1] isa Symbol
        @test InfiniteOpt._parse_idxset(:a)[2] == esc(:a)
    end
    # test _explicit_oneto
    @testset "_explicit_oneto" begin
        @test InfiniteOpt._explicit_oneto(:(1:2)) == :(Base.OneTo(1:2))
        s = esc(:(1:2))
        @test InfiniteOpt._explicit_oneto(s) == :(Base.OneTo($s))
        @test InfiniteOpt._explicit_oneto(:(2:4)) == :(2:4)
    end
    # test _expr_is_splat
    @testset "_expr_is_splat" begin
        @test !InfiniteOpt._expr_is_splat(:(1:2))
        @test !InfiniteOpt._expr_is_splat(esc(:(1:2)))
        @test InfiniteOpt._expr_is_splat(:(a...))
        @test !InfiniteOpt._expr_is_splat(:a)
    end
    # test _parse_ref_sets
    @testset "_parse_ref_sets" begin
        # test errors
        bad_vcat = :([1:2; k; j])
        @test_throws ErrorException InfiniteOpt._parse_ref_sets(error, bad_vcat)
        bad_vect = :([1:2, 1:2; k; j])
        @test_throws ErrorException InfiniteOpt._parse_ref_sets(error, bad_vect)
        duplicate_ref = :(a[i in 1:2, i in 1:2])
        @test_throws ErrorException InfiniteOpt._parse_ref_sets(error, duplicate_ref)
        bad_type = :((i in 2:3))
        @test_throws ErrorException InfiniteOpt._parse_ref_sets(error, bad_type)
        # test vcat
        vcat = :([i = s; i <= 2])
        @test InfiniteOpt._parse_ref_sets(error, vcat)[1] == [:i]
        @test InfiniteOpt._parse_ref_sets(error, vcat)[2] == [esc(:s)]
        @test InfiniteOpt._parse_ref_sets(error, vcat)[3] == :(i <= 2)
        # test typed_vcat
        typed_vcat = :(a[i = 2:w; k])
        @test InfiniteOpt._parse_ref_sets(error, typed_vcat)[1] == [:i]
        @test InfiniteOpt._parse_ref_sets(error, typed_vcat)[2] == [esc(:(2:w))]
        @test InfiniteOpt._parse_ref_sets(error, typed_vcat)[3] == :k
        # test ref w/o condition
        ref = :(a[1:2, i = S])
        @test length(InfiniteOpt._parse_ref_sets(error, ref)[1]) == 2 
        @test InfiniteOpt._parse_ref_sets(error, ref)[2] == [esc(:(1:2)), esc(:S)]
        @test InfiniteOpt._parse_ref_sets(error, ref)[3] == :()
        # test ref w/ condition
        ref = :(a[i in 2:4, j = [:d, :e]; i >= 3])
        @test InfiniteOpt._parse_ref_sets(error, ref)[1] == [:i, :j]
        @test InfiniteOpt._parse_ref_sets(error, ref)[2] == [esc(:(2:4)), esc(:([:d, :e]))]
        @test InfiniteOpt._parse_ref_sets(error, ref)[3] == :(i >= 3)
        # test vect w/o condition
        vect = :([s in S, 4:5])
        @test InfiniteOpt._parse_ref_sets(error, vect)[1][1] == :s 
        @test InfiniteOpt._parse_ref_sets(error, vect)[1][2] isa Symbol
        @test InfiniteOpt._parse_ref_sets(error, vect)[2] == [esc(:S), esc(:(4:5))]
        @test InfiniteOpt._parse_ref_sets(error, vect)[3] == :()
        # test vect w/ condition
        vect = :([s in S, j = 4:5; s + j == 6])
        @test InfiniteOpt._parse_ref_sets(error, vect)[1] == [:s, :j]
        @test InfiniteOpt._parse_ref_sets(error, vect)[2] == [esc(:S), esc(:(4:5))]
        @test InfiniteOpt._parse_ref_sets(error, vect)[3] == :(s + j == 6)
        # test symbol
        @test InfiniteOpt._parse_ref_sets(error, :a) == (Any[], Any[], :())
    end
    # test _build_ref_sets
    @testset "_build_ref_sets" begin
        # test errors 
        @test_throws ErrorException InfiniteOpt._build_ref_sets(error, :([a...]))
        bad_vcat = :([1:2; k; j])
        @test_throws ErrorException InfiniteOpt._build_ref_sets(error, bad_vcat)
        # test dependent build
        cont = JuMPC.SparseAxisArray(Dict((1, 1) => 2, (2, 1) => 3, (2, 2) => 4))
        @test @gen_container([i = 1:2, j = 1:i], i + j) == cont 
        # test conditioned build 
        cont = JuMPC.SparseAxisArray(Dict((1, 2) => 3, (1, 3) => 4, (2, 2) => 4))
        @test @gen_container([i = 1:2, j = 2:3; i + j <= 4], i + j) == cont 
        # test easy build 
        cont = JuMPC.DenseAxisArray(ones(2, 2), Base.OneTo(2), 2:3)
        @test @gen_container([1:2, j = 2:3], 1.0) == cont 
    end
end
