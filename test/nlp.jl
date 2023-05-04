# Test the basic expression generation via operators
@testset "Operator Definition" begin 
    # setup model data 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    # test the sum 
    @testset "sum" begin 
        # test empty sums
        @test sum(i for i in Int[]) == 0
        @test isequal(sum(NLPExpr[]), zero(NLPExpr))
        @test isequal(sum(GeneralVariableRef[]), 
                      zero(GenericAffExpr{Float64, GeneralVariableRef}))
        # test NLP sums 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        new = NLPExpr(InfiniteOpt._call_graph(:+, nlp, nlp))
        @test isequal(sum(i for i in [nlp, nlp]), new)
        @test isequal(sum([nlp, nlp]), new)
        @test_throws ErrorException sum(i for i in [nlp, nlp]; bad = 42)
        @test_throws ErrorException sum(i for i in [y, t]; bad = 42)
        # test other expressions 
        @test isequal(sum(i for i in [y, t]), y + t)
        @test isequal(sum([y, t]), y + t)
        # normal sums 
        @test sum(i for i in 1:3) == 6
    end
    # test the product
    @testset "prod" begin 
        # test empty sums
        @test prod(i for i in Int[]) == 1
        @test isequal(prod(NLPExpr[]), one(NLPExpr))
        # test NLP sums 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        new = NLPExpr(InfiniteOpt._call_graph(:*, nlp, nlp))
        @test isequal(prod(i for i in [nlp, nlp]), new)
        @test isequal(prod([nlp, nlp]), new)
        @test_throws ErrorException prod(i for i in [nlp, nlp]; bad = 42)
        @test_throws ErrorException prod(i for i in [y, t, z]; bad = 42)
        # test other expressions 
        new = NLPExpr(InfiniteOpt._call_graph(:*, y, t, z))
        @test isequal(prod(i for i in [y, t, z]), new)
        @test isequal(prod([y, t, z]), new)
        # test normal products 
        @test prod(i for i in 1:3) == 6
    end
    # prepare the test grid 
    aff = 2z - 2
    quad = y^2 + y
    p = Node(NodeData(:sin))
    addchild(p, NodeData(y))
    nlp = NLPExpr(p)
    # test the multiplication operator
    @testset "Multiplication Operator" begin
        for (i, iorder) in [(42, 0), (y, 1), (aff, 1), (quad, 2), (nlp, 3)]
            for (j, jorder) in [(42, 0), (y, 1), (aff, 1), (quad, 2), (nlp, 3)]
                if iorder + jorder >= 3
                    expected = NLPExpr(InfiniteOpt._call_graph(:*, i, j))
                    @test isequal(i * j, expected)
                end
            end
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:*, y, t, z,nlp))
        @test isequal(y * t * z * nlp, expected)
        @test isequal(*(nlp), nlp)
    end
    # test division operator 
    @testset "Division Operator" begin
        for i in [42, y, aff, quad, nlp]
            for j in [y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:/, i, j))
                @test isequal(i / j, expected)
            end
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:/, nlp, 42))
        @test isequal(nlp / 42, expected)
        @test isequal(nlp / 1, nlp)
        @test_throws ErrorException nlp / 0
    end
    # test the power operator 
    @testset "Power Operator" begin
        for i in [42, y, aff, quad, nlp]
            for j in [y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:^, i, j))
                @test isequal(i ^ j, expected)
            end
        end
        one_aff = one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i in [y, aff, quad, nlp]
            for f in [Float64, Int]
                @test isequal(i^zero(f), one_aff)
                @test isequal(i^one(f), i)
                if i isa NLPExpr
                    expected = NLPExpr(InfiniteOpt._call_graph(:^, i, f(2)))
                    @test isequal(i ^ f(2), expected)
                else
                    @test isequal(i^f(2), i * i)
                end
                expected = NLPExpr(InfiniteOpt._call_graph(:^, i, f(3)))
                @test isequal(i ^ f(3), expected)
            end
        end
        # extra tests 
        @test isequal(y^0, one_aff)
        @test isequal(y^1, y)
        @test isequal(y^2, y * y)
        @test isequal(y^0.0, one_aff)
        @test isequal(y^1.0, y)
        @test isequal(y^2.0, y * y)
    end
    # test the subtraction operator
    @testset "Subtraction Operator" begin
        for i in [42, y, aff, quad, nlp]
            expected = NLPExpr(InfiniteOpt._call_graph(:-, nlp, i))
            @test isequal(nlp - i, expected)
            expected = NLPExpr(InfiniteOpt._call_graph(:-, i, nlp))
            @test isequal(i - nlp, expected)
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:-, nlp))
        @test isequal(-nlp, expected)
        @test isequal(y - y, zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})) 
        @test (y - z).constant == 0 
        @test (y - z).terms[y] == 1
        @test (y - z).terms[z] == -1
    end
    # test the addition operator 
    @testset "Addition Operator" begin
        for i in [42, y, aff, quad, nlp]
            expected = NLPExpr(InfiniteOpt._call_graph(:+, nlp, i))
            @test isequal(nlp + i, expected)
            expected = NLPExpr(InfiniteOpt._call_graph(:+, i, nlp))
            @test isequal(i + nlp, expected)
        end
        @test isequal(+nlp, nlp)
    end
    # test the comparison operators 
    @testset "Comparison Operators" begin 
        for i in [42, y, aff, quad, nlp]
            for j in [42, y, aff, quad, nlp]
                for (n, f) in (:< => Base.:(<), :(==) => Base.:(==), 
                               :> => Base.:(>), :<= => Base.:(<=), 
                               :>= => Base.:(>=))
                    if !(isequal(i, 42) && isequal(j, 42)) && !(isequal(i, y) && isequal(j, y))
                        expected = NLPExpr(InfiniteOpt._call_graph(n, i, j))
                        @test isequal(f(i, j), expected)
                    end
                end
            end
        end
        for (n, f) in (:(==) => Base.:(==), :<= => Base.:(<=), :>= => Base.:(>=))
            expected = NLPExpr(InfiniteOpt._call_graph(n, y, z))
            @test isequal(f(y, z), expected)
            @test f(y, y)
        end
        for (n, f) in (:< => Base.:(<), :> => Base.:(>))
            expected = NLPExpr(InfiniteOpt._call_graph(n, y, z))
            @test isequal(f(y, z), expected)
            @test !f(y, y)
        end
        # extra tests
        @test (y == 0) isa NLPExpr
        @test (y <= 0) isa NLPExpr
        @test (y >= 0) isa NLPExpr
        @test (y > 0) isa NLPExpr
        @test (y < 0) isa NLPExpr
        @test (0 == y) isa NLPExpr
        @test (0 <= y) isa NLPExpr
        @test (0 >= y) isa NLPExpr
        @test (0 > y) isa NLPExpr
        @test (0 < y) isa NLPExpr
    end
    # test the logic operators 
    @testset "Logic Operators" begin 
        for i in [y, nlp]
            for j in [y, nlp]
                for (n, f) in (:&& => Base.:&, :|| => Base.:|)
                    expected = NLPExpr(InfiniteOpt._call_graph(n, i, j))
                    @test isequal(f(i, j), expected)
                end
            end
        end
        for i in [y, nlp]
            @test !(i & false)
            @test !(false & i)
            @test isequal(i & true, i)
            @test isequal(true & i, i)
            @test (i | true)
            @test (true | i)
            @test isequal(i | false, i)
            @test isequal(false | i, i)
        end
    end
    # test ifelse 
    @testset "ifelse" begin
        for i in [42, y, aff, quad, nlp]
            for j in [42, y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:ifelse, nlp, i, j))
                @test isequal(InfiniteOpt.ifelse(nlp, i, j), expected)
            end
        end
        @test isequal(InfiniteOpt.ifelse(true, y, z), y)
        @test isequal(InfiniteOpt.ifelse(false, y, z), z)
    end
    # test the default functions 
    @testset "Default Registered Functions" begin 
        for i in [y, aff, quad, nlp]
            for (n, f) in InfiniteOpt._Base1ArgFuncList
            expected = NLPExpr(InfiniteOpt._call_graph(n, i))
            @test isequal(f(i), expected)
            end
        end
        for i in [y, aff, quad, nlp]
            for (n, f) in InfiniteOpt._Special1ArgFuncList
            expected = NLPExpr(InfiniteOpt._call_graph(n, i))
            @test isequal(f(i), expected)
            end
        end
    end
end


# Test registration utilities
@testset "Registration Methods" begin
    # setup model data 
    m = InfiniteModel()
    @variable(m, y)
    # test name_to_function
    @testset "name_to_function" begin 
        @test name_to_function(m, :tan, 1) == tan
        @test name_to_function(m, :+, 10) == +
        @test name_to_function(m, :*, 2) == *
    end
    # test all_registered_functions
    @testset "all_registered_functions" begin 
        @test all_registered_functions(m) isa Vector{Function}
    end
    # test user_registered_functions
    @testset "user_registered_functions" begin 
        @test user_registered_functions(m) == RegisteredFunction[]
    end
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
    # test creation helper errors 
    @testset "Registration Helpers" begin 
        @test_throws ErrorException RegisteredFunction(:a, 1, f, g)
        @test_throws ErrorException RegisteredFunction(:a, 2, f, g)
        @test_throws ErrorException RegisteredFunction(:a, 1, f, g, f)
        @test_throws ErrorException RegisteredFunction(:a, 1, f, f, g)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :f, 1, 1)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :sin, 1, sin)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :g, 1, g)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :eta, 1, eta)
    end
    # test @register
    @testset "@register" begin 
        # test errors 
        @test_macro_throws ErrorException @register(Model(), f(a))
        @test_macro_throws ErrorException @register(m, f(a), bad = 42)
        @test_macro_throws ErrorException @register(m, f(a), g, h, 3)
        @test_macro_throws ErrorException @register(m, f[i](a))
        @test_macro_throws ErrorException @register(m, f(a::Int))
        @test_macro_throws ErrorException @register(m, f(a), 2)
        @test_macro_throws ErrorException @register(m, g(a))
        @test_macro_throws ErrorException @register(m, eta(a))
        @test_macro_throws ErrorException @register(m, f(a), g)
        @test_macro_throws ErrorException @register(m, f(a), g, g)
        @test_macro_throws ErrorException @register(m, f(a), f, g)
        @test_macro_throws ErrorException @register(m, h(a, b), h)
        @test_macro_throws ErrorException @register(m, sin(a))
        # test univariate function with no gradient or hessian
        @test @register(m, f(a)) isa Function
        @test f(y) isa NLPExpr
        @test f(y).tree_root.data.value == :f
        @test isequal(f(y).tree_root.child.data.value, y)
        @test f(2) == 8
        @test length(user_registered_functions(m)) == 1
        @test user_registered_functions(m)[1].name == :f
        @test user_registered_functions(m)[1].func == f
        @test user_registered_functions(m)[1].num_args == 1
        @test user_registered_functions(m)[1].gradient isa Nothing 
        @test user_registered_functions(m)[1].hessian isa Nothing 
        @test name_to_function(m, :f, 1) == f
        # test univariate function with gradient
        @test @register(m, f1(a), f) isa Function
        @test f1(y) isa NLPExpr
        @test f1(y).tree_root.data.value == :f1
        @test isequal(f1(y).tree_root.child.data.value, y)
        @test f1(2) == 32
        @test length(user_registered_functions(m)) == 2
        @test user_registered_functions(m)[2].name == :f1
        @test user_registered_functions(m)[2].func == f1
        @test user_registered_functions(m)[2].num_args == 1
        @test user_registered_functions(m)[2].gradient == f
        @test user_registered_functions(m)[2].hessian isa Nothing 
        @test name_to_function(m, :f1, 1) == f1
        # test univariate function with gradient and hessian
        @test @register(m, f2(a), f, f1) isa Function
        @test f2(y) isa NLPExpr
        @test f2(y).tree_root.data.value == :f2
        @test isequal(f2(y).tree_root.child.data.value, y)
        @test f2(2) == 10
        @test length(user_registered_functions(m)) == 3
        @test user_registered_functions(m)[3].name == :f2
        @test user_registered_functions(m)[3].func == f2
        @test user_registered_functions(m)[3].num_args == 1
        @test user_registered_functions(m)[3].gradient == f
        @test user_registered_functions(m)[3].hessian == f1
        @test name_to_function(m, :f2, 1) == f2
        # test multivariate function with no gradient
        @test @register(m, h(a, b)) isa Function
        @test h(2, y) isa NLPExpr
        @test h(2, y).tree_root.data.value == :h
        @test isequal(h(2, y).tree_root.child.data.value, 2)
        @test isequal(h(2, y).tree_root.child.sibling.data.value, y)
        @test h(y, y) isa NLPExpr
        @test h(y, 2) isa NLPExpr
        @test h(2, 2) == 42
        @test length(user_registered_functions(m)) == 4
        @test user_registered_functions(m)[4].name == :h
        @test user_registered_functions(m)[4].func == h
        @test user_registered_functions(m)[4].num_args == 2
        @test user_registered_functions(m)[4].gradient isa Nothing
        @test user_registered_functions(m)[4].hessian isa Nothing
        @test name_to_function(m, :h, 2) == h
        # test multivariate function with gradient
        @test @register(m, h1(a, b), hg) isa Function
        @test h1(2, y) isa NLPExpr
        @test h1(2, y).tree_root.data.value == :h1
        @test isequal(h1(2, y).tree_root.child.data.value, 2)
        @test isequal(h1(2, y).tree_root.child.sibling.data.value, y)
        @test h1(y, y) isa NLPExpr
        @test h1(y, 2) isa NLPExpr
        @test h1(2, 2) == 13
        @test length(user_registered_functions(m)) == 5
        @test user_registered_functions(m)[5].name == :h1
        @test user_registered_functions(m)[5].func == h1
        @test user_registered_functions(m)[5].num_args == 2
        @test user_registered_functions(m)[5].gradient == hg
        @test user_registered_functions(m)[5].hessian isa Nothing
        @test name_to_function(m, :h1, 2) == h1
        # test wrong model error 
        @test_throws ErrorException f(@variable(InfiniteModel()))
        # test functional registration
        function registration_test()
            mt = InfiniteModel()
            @variable(mt, x)
            q(a) = 1
            @test @register(mt, q(a)) isa Function 
            @test q(x) isa NLPExpr
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
        r1 = m1.nlp_model.operators
        @test length(r1.registered_univariate_operators) == 3
        @test [r1.registered_univariate_operators[i].f for i in 1:3] == [f, f1, f2]
        @test [r1.registered_univariate_operators[i].f′ for i in 2:3] == [f, f]
        @test r1.registered_univariate_operators[3].f′′ == f1
        @test length(r1.registered_multivariate_operators) == 2
        # test error 
        m2 = Model()
        h2(a, b) = 3
        @test @register(m, h2(a, b), hg, f) isa Function
        @test_throws ErrorException add_registered_to_jump(m2, m) isa Nothing 
        r2 = m2.nlp_model.operators
        @test length(r2.registered_univariate_operators) == 3
        @test [r2.registered_univariate_operators[i].f for i in 1:3] == [f, f1, f2]
        @test [r2.registered_univariate_operators[i].f′ for i in 2:3] == [f, f]
        @test r2.registered_univariate_operators[3].f′′ == f1
        @test length(r2.registered_multivariate_operators) == 2
    end
end 

# Test string methods
@testset "String Methods" begin
    # setup the model data
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    # test making strings
    @testset "String Creation" begin
        # test some simple ones
        @test string(sin(y) + 2) == "sin(y(t)) + 2"
        @test string((z*z) ^ 4) == "(z²)^4"
        @test string((cos(z) + sin(z)) / y) == "(cos(z) + sin(z)) / y(t)"
        @test string(-cos(z + y) * z^2.3) == "-cos(z + y(t)) * z^2.3"
        @test string((-InfiniteOpt.ifelse(z == 0, 0, z)) ^ 3) == "(-(ifelse(z == 0, 0, z))^3"
        # test AffExpr cases 
        aff0 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        @test string(aff0^4 + (2z - 3y + 42) ^ (-1z)) == "0^4 + (2 z - 3 y(t) + 42)^(-z)"
        @test string((1z) / (2.3 * y)) == "z / (2.3 y(t))"
        # test QuadExpr cases 
        quad0 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        @test string(quad0 / (z^2 + y * z)) == "0 / (z² + y(t)*z)"
        @test string((0z^2 + 42) * sin(y)) == "42 * sin(y(t))"
        @test string((z^2) / y) == "(z²) / y(t)"
    end
end
