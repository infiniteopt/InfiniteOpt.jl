# Test operator utilities
@testset "Operator Methods" begin
    # setup model data 
    m = InfiniteModel()
    @variable(m, y)
    # define functions for tests 
    f(a) = a^3
    g(a::Int) = 42
    h(a, b) = 42
    f1(a) = 32
    f2(a) = 10
    h1(a, b) = 13
    function hg(v, a, b)
        v[1] = 1
        v[2] = 2
        return 
    end
    function ∇²h(H, x...)
        H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
        H[2, 1] = -400 * x[1]
        H[2, 2] = 200.0
        return
    end
    # test JuMP.add_nonlinear_operator
    @testset "JuMP.add_nonlinear_operator" begin
        # test errors
        @test_throws ErrorException add_nonlinear_operator(m, 1, f, name = :max)
        m.op_lookup[:f] = (f, 1)
        @test_throws ErrorException add_nonlinear_operator(m, 1, f)
        empty!(m.op_lookup)
        @test_throws ErrorException add_nonlinear_operator(m, 2, f)
        # test normal 
        @test add_nonlinear_operator(m, 1, f).head == :f
        @test m.op_lookup[:f] == (f, 1)
        @test last(m.operators).f == f
        @test last(m.operators).dim == 1
        @test last(m.operators).name == :f
        @test add_nonlinear_operator(m, 2, h, hg).head == :h 
        @test m.op_lookup[:h] == (h, 2)
        @test last(m.operators).f == h
        @test last(m.operators).dim == 2
        @test last(m.operators).name == :h
        @test last(m.operators).∇f == hg
    end 
    # test name_to_operator
    @testset "name_to_operator" begin 
        @test name_to_operator(m, :f) == f
        @test name_to_operator(m, :h) == h
        @test name_to_operator(m, :bad) isa Nothing
    end
    # test all_nonlinear_operators
    @testset "all_nonlinear_operators" begin 
        @test all_nonlinear_operators(m) isa Vector{Symbol}
    end
    # test added_nonlinear_operators
    @testset "added_nonlinear_operators" begin 
        @test added_nonlinear_operators(m) isa Vector{NLPOperator}
    end
    empty!(m.operators)
    empty!(m.op_lookup)
    # test @operator
    @testset "@operator" begin 
        # test errors 
        @test @operator(m, f1, 1, f) isa NonlinearOperator
        @test @operator(m, f2, 1, f, f) isa NonlinearOperator
        @test @operator(m, f3, 1, f, f, f) isa NonlinearOperator
        @test @operator(m, h1, 2, h) isa NonlinearOperator
        @test @operator(m, h2, 2, h, hg) isa NonlinearOperator
        @test @operator(m, h3, 2, h, hg, ∇²h) isa NonlinearOperator
        # test operator scoping
        function scope_test()
            mt = InfiniteModel()
            @variable(mt, x)
            q(a) = 1
            @test @operator(mt, my_q, 1, q) isa NonlinearOperator
            q(x::GeneralVariableRef) = GenericNonlinearExpr{GeneralVariableRef}(:my_q, x)
            @test @expression(mt, q(x)) isa GenericNonlinearExpr
            return 
        end
        @test scope_test() isa Nothing 
        @test scope_test() isa Nothing 
    end
    # test @register
    @testset "@register" begin
        @test_macro_throws ErrorException @register(m, f(a))
    end
    # test add_operators_to_jump
    @testset "add_operators_to_jump" begin 
        # test normal 
        m1 = Model()
        @test add_operators_to_jump(m1, m) isa Nothing 
        attr_dict = backend(m1).model_cache.modattr
        @test length(attr_dict) == 6
        @test attr_dict[MOI.UserDefinedFunction(:f1, 1)] == (f,)
        @test attr_dict[MOI.UserDefinedFunction(:f2, 1)] == (f, f)
        @test attr_dict[MOI.UserDefinedFunction(:f3, 1)] == (f, f, f)
        @test attr_dict[MOI.UserDefinedFunction(:h1, 2)] == (h,)
        @test attr_dict[MOI.UserDefinedFunction(:h2, 2)] == (h, hg)
        @test attr_dict[MOI.UserDefinedFunction(:h3, 2)] == (h, hg, ∇²h)
    end
end 
