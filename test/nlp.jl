# Test registration utilities
@testset "Registration Methods" begin
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
    # test JuMP.add_user_defined_function
    @testset "JuMP.add_user_defined_function" begin
        # test errors
        @test_throws ErrorException add_user_defined_function(m, :f, 1)
        @test_throws ErrorException add_user_defined_function(m, :f, 1, f, 2)
        @test_throws ErrorException add_user_defined_function(m, :max, 1, f)
        m.func_lookup[:f] = (f, 1)
        @test_throws ErrorException add_user_defined_function(m, :f, 1, f)
        empty!(m.func_lookup)
        @test_throws ErrorException add_user_defined_function(m, :f, 2, f)
        # test normal 
        @test add_user_defined_function(m, :f, 1, f).head == :f
        @test m.func_lookup[:f] == (f, 1)
        @test last(m.registrations).f == f
        @test last(m.registrations).dim == 1
        @test last(m.registrations).op == :f
        @test add_user_defined_function(m, :h, 2, h, hg).head == :h 
        @test m.func_lookup[:h] == (h, 2)
        @test last(m.registrations).f == h
        @test last(m.registrations).dim == 2
        @test last(m.registrations).op == :h
        @test last(m.registrations).∇f == hg
    end 
    # test name_to_function
    @testset "name_to_function" begin 
        @test name_to_function(m, :f) == f
        @test name_to_function(m, :h) == h
        @test name_to_function(m, :bad) isa Nothing
    end
    # test all_registered_functions
    @testset "all_registered_functions" begin 
        @test all_registered_functions(m) isa Vector{Symbol}
    end
    # test user_registered_functions
    @testset "user_registered_functions" begin 
        @test user_registered_functions(m) isa Vector{RegisteredFunction}
    end
    empty!(m.registrations)
    empty!(m.func_lookup)
    # test @register
    @testset "@register" begin 
        # test errors 
        @test @register(m, f1, 1, f) isa UserDefinedFunction
        @test @register(m, f2, 1, f, f) isa UserDefinedFunction
        @test @register(m, f3, 1, f, f, f) isa UserDefinedFunction
        @test @register(m, h1, 2, h) isa UserDefinedFunction
        @test @register(m, h2, 2, h, hg) isa UserDefinedFunction
        @test @register(m, h3, 2, h, hg, ∇²h) isa UserDefinedFunction
        # test functional registration
        function registration_test()
            mt = InfiniteModel()
            @variable(mt, x)
            q(a) = 1
            @test @register(mt, q, 1, q) isa UserDefinedFunction
            @test @expression(mt, q(x)) isa GenericNonlinearExpr
            return 
        end
        @test registration_test() isa Nothing 
        @test registration_test() isa Nothing 
    end
    # test add_registered_to_jump
    @testset "add_registered_to_jump" begin 
        # test normal 
        m1 = Model()
        @test add_registered_to_jump(m1, m) isa Nothing 
        # TODO update checks below
        # r1 = m1.nlp_model.operators
        # @test length(r1.registered_univariate_operators) == 3
        # @test [r1.registered_univariate_operators[i].f for i in 1:3] == [f, f1, f2]
        # @test [r1.registered_univariate_operators[i].f′ for i in 2:3] == [f, f]
        # @test r1.registered_univariate_operators[3].f′′ == f1
        # @test length(r1.registered_multivariate_operators) == 2
    end
end 
