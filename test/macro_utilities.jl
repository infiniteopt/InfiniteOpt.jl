# Test basic methods 
@testset "Basics" begin 
    # test _esc_non_constant
    @testset "_esc_non_constant" begin 
        @test InfiniteOpt._esc_non_constant(42) == 42
        @test InfiniteOpt._esc_non_constant(:(quote 42 end)) == :(quote 42 end)
        @test InfiniteOpt._esc_non_constant(:(x[i])) == esc(:(x[i]))
        @test InfiniteOpt._esc_non_constant(:x) == esc(:x)
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
    # test _finalize_macro
    @testset "_finalize_macro" begin
        m = InfiniteModel()
        esc_m = esc(:m)
        code = :(GeneralVariableRef($esc_m, 1, FiniteVariableIndex))
        source = LineNumberNode(42, :bob)
        @test InfiniteOpt._finalize_macro(error, esc_m, code, source, :test) isa Expr
    end
end
