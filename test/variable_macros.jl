# Test helper functions for infinite variable macro
@testset "Infinite Helpers" begin
    # _check_rhs
    @testset "_check_rhs" begin
        # test with reversed sides
        arg1 = :(data[i]); arg2 = :(x[i=1:2])
        @test InfiniteOpt._check_rhs(arg1, arg2) == (arg2, arg1)
        # test with normal case that shouldn't be swapped
        arg1 = :((x[i=1:2])(t)); arg2 = :(data[i])
        @test InfiniteOpt._check_rhs(arg1, arg2) == (arg1, arg2)
        # test reversed case that stays reversed because cannot be distinguished
        arg1 = :(data(i)); arg2 = :((x[i=1:2])(t))
        @test InfiniteOpt._check_rhs(arg1, arg2) == (arg1, arg2)
    end
    # _less_than_parse
    @testset "_less_than_parse" begin
        # test with reversed sides
        arg1 = :(data[i]); arg2 = :(x[i=1:2])
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:($(arg2) >= $(arg1)),
                                                           nothing)
        # test normal reference with parameter tuple
        arg1 = :(x[i=1:2](t)); arg2 = :(data[i])
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x[i=1:2] <= $(arg2)),
                                                           :((t,)))
        # test normal with parameter tuple
        arg1 = :(x(t)); arg2 = :(data[i])
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x <= $(arg2)),
                                                           :((t,)))
        # test normal without tuple
        arg1 = :(x); arg2 = :(data[i])
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x <= $(arg2)),
                                                           nothing)
        # test normal without tuple
        arg1 = :(x); arg2 = 1
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x <= $(arg2)),
                                                           nothing)
    end
    # _greater_than_parse
    @testset "_greater_than_parse" begin
        # test with reversed sides
        arg1 = :(data[i]); arg2 = :(x[i=1:2])
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:($(arg2) <= $(arg1)),
                                                              nothing)
        # test normal reference with parameter tuple
        arg1 = :(x[i=1:2](t)); arg2 = :(data[i])
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x[i=1:2] >= $(arg2)),
                                                              :((t,)))
        # test normal with parameter tuple
        arg1 = :(x(t)); arg2 = :(data[i])
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x >= $(arg2)),
                                                              :((t,)))
        # test normal without tuple
        arg1 = :(x); arg2 = :(data[i])
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x >= $(arg2)),
                                                              nothing)
        # test normal without tuple
        arg1 = :(x); arg2 = 1
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x >= $(arg2)),
                                                              nothing)
    end
    # _less_than_parse (number on lhs)
    @testset "_less_than_parse (reversed)" begin
        # test with reference
        arg1 = 1; arg2 = :(x[i=1:2])
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:($(arg2) >= $(arg1)),
                                                           nothing)
        # test with reference and parameter tuple
        arg1 = 1; arg2 = :(x[i=1:2](t))
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x[i=1:2] >= $(arg1)),
                                                           :((t,)))
        # test normal without tuple
        arg1 = 1; arg2 = :(x)
        @test InfiniteOpt._less_than_parse(arg1, arg2) == (:(x >= $(arg1)),
                                                           nothing)
    end
    # _greater_than_parse (number on lhs)
    @testset "_greater_than_parse (reversed)" begin
        # test with reference
        arg1 = 1; arg2 = :(x[i=1:2])
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:($(arg2) <= $(arg1)),
                                                              nothing)
        # test with reference and parameter tuple
        arg1 = 1; arg2 = :(x[i=1:2](t))
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x[i=1:2] <= $(arg1)),
                                                              :((t,)))
        # test normal without tuple
        arg1 = 1; arg2 = :(x)
        @test InfiniteOpt._greater_than_parse(arg1, arg2) == (:(x <= $(arg1)),
                                                              nothing)
    end
    # _equal_to_parse
    @testset "_equal_to_parse" begin
        # test with reversed sides
        arg1 = :(data[i]); arg2 = :(x[i=1:2])
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:($(arg2) == $(arg1)),
                                                          nothing)
        # test normal reference with parameter tuple
        arg1 = :(x[i=1:2](t)); arg2 = :(data[i])
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x[i=1:2] == $(arg2)),
                                                          :((t,)))
        # test normal with parameter tuple
        arg1 = :(x(t)); arg2 = :(data[i])
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x == $(arg2)),
                                                          :((t,)))
        # test normal without tuple
        arg1 = :(x); arg2 = :(data[i])
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x == $(arg2)),
                                                          nothing)
        # test normal without tuple
        arg1 = :(x); arg2 = 1
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x == $(arg2)),
                                                          nothing)
    end
    # _equal_to_parse (number on lhs)
    @testset "_equal_to_parse (reversed)" begin
        # test with reference
        arg1 = 1; arg2 = :(x[i=1:2])
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:($(arg2) == $(arg1)),
                                                          nothing)
        # test with reference and parameter tuple
        arg1 = 1; arg2 = :(x[i=1:2](t))
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x[i=1:2] == $(arg1)),
                                                          :((t,)))
        # test normal without tuple
        arg1 = 1; arg2 = :(x)
        @test InfiniteOpt._equal_to_parse(arg1, arg2) == (:(x == $(arg1)),
                                                          nothing)
    end
    # _parse_parameters (call)
    @testset "_parse_parameters (call)" begin
        # test less than parse
        expr = :(x(t, x) <= 1)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                            expr.args) == (:(x <= 1), :((t, x)))
        # test greater than parse
        expr = :(x[1:2](t) >= 1)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                         expr.args) == (:(x[1:2] >= 1), :((t,)))
        # test equal to parse
        expr = :(x(t) == d)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                            expr.args) == (:(x == d), :((t,)))
        # test only variable parse
        expr = :(x(t))
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                            expr.args) == (:(x), :((t,)))
        # test invalid use of operator
        expr = :(x(t) in 1)
        @test_throws ErrorException InfiniteOpt._parse_parameters(error,
                                                                  Val(expr.head),
                                                                  expr.args)
    end
    # _parse_parameters (comparison)
    @testset "_parse_parameters (compare)" begin
        # test with parameters
        expr = :(0 <= x(t, x) <= 1)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                       expr.args) == (:(0 <= x <= 1), :((t, x)))
        # test with parameters and references
        expr = :(0 <= x[1:2](t, x) <= 1)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                  expr.args) == (:(0 <= x[1:2] <= 1), :((t, x)))
        # test without parameters
        expr = :(0 <= x <= 1)
        @test InfiniteOpt._parse_parameters(error, Val(expr.head),
                                         expr.args) == (:(0 <= x <= 1), nothing)
    end
end

# Test the infinite variable macro
@testset "Infinite" begin
    # initialize model and infinite parameters
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= t <= 1)
    @infinite_parameter(m, -1 <= x[1:2] <= 1)
    # test single variable definition
    @testset "Single" begin
        # test basic defaults
        vref = InfiniteVariableRef(m, 1)
        @test @infinite_variable(m, parameter_refs = t) == vref
        @test name(vref) == "noname(t)"
        # test more tuple input and variable details
        vref = InfiniteVariableRef(m, 2)
        @test @infinite_variable(m, parameter_refs = (t, x), base_name = "test",
                                 binary = true) == vref
        @test name(vref) == "test(t, x)"
        @test is_binary(vref)
        # test nonanonymous with simple single arg
        vref = InfiniteVariableRef(m, 3)
        @test @infinite_variable(m, a(x)) == vref
        @test name(vref) == "a(x)"
        # test nonanonymous with complex single arg
        vref = InfiniteVariableRef(m, 4)
        @test @infinite_variable(m, 0 <= b(x) <= 1) == vref
        @test name(vref) == "b(x)"
        @test lower_bound(vref) == 0
        @test upper_bound(vref) == 1
        # test nonanonymous with reversed single arg
        vref = InfiniteVariableRef(m, 5)
        @test @infinite_variable(m, 0 <= c(t)) == vref
        @test name(vref) == "c(t)"
        @test lower_bound(vref) == 0
        # test multi-argument expr 1
        vref = InfiniteVariableRef(m, 6)
        @test @infinite_variable(m, d(t) == 0, Int, base_name = "test") == vref
        @test name(vref) == "test(t)"
        @test fix_value(vref) == 0
    end
    # test array variable definition
    @testset "Array" begin
        # test anonymous array
        vrefs = [InfiniteVariableRef(m, 7), InfiniteVariableRef(m, 8)]
        @test @infinite_variable(m, [1:2], parameter_refs = t) == vrefs
        @test name(vrefs[1]) == "noname(t)"
        # test basic param expression
        vrefs = [InfiniteVariableRef(m, 9), InfiniteVariableRef(m, 10)]
        @test @infinite_variable(m, e[1:2], parameter_refs = (t, x)) == vrefs
        @test name(vrefs[2]) == "e[2](t, x)"
        # test comparison without params
        vrefs = [InfiniteVariableRef(m, 11), InfiniteVariableRef(m, 12)]
        @test @infinite_variable(m, 0 <= f[1:2] <= 1,
                                 parameter_refs = (t, x)) == vrefs
        @test name(vrefs[2]) == "f[2](t, x)"
        @test lower_bound(vrefs[1]) == 0
        @test upper_bound(vrefs[2]) == 1
        # test comparison with call
        vrefs = [InfiniteVariableRef(m, 13), InfiniteVariableRef(m, 14)]
        @test @infinite_variable(m, 0 <= g[1:2](t) <= 1) == vrefs
        @test name(vrefs[1]) == "g[1](t)"
        @test lower_bound(vrefs[1]) == 0
        @test upper_bound(vrefs[2]) == 1
        # test fixed
        vrefs = [InfiniteVariableRef(m, 15), InfiniteVariableRef(m, 16)]
        @test @infinite_variable(m, h[i = 1:2](t) == ones(2)[i]) == vrefs
        @test name(vrefs[1]) == "h[1](t)"
        @test fix_value(vrefs[1]) == 1
        # test containers
        vrefs = [InfiniteVariableRef(m, 17), InfiniteVariableRef(m, 18)]
        vrefs = convert(JuMP.Containers.SparseAxisArray, vrefs)
        @test @infinite_variable(m, [1:2](t),
                                 container = SparseAxisArray) == vrefs
        @test name(vrefs[1]) == "noname(t)"
    end
    # test errors
    @testset "Errors" begin
        # test model assertion errors
        m2 = Model()
        @test_throws AssertionError @infinite_variable(m2, parameter_refs = t)
        @test_throws AssertionError @infinite_variable(m2, i(t))
        @test_throws AssertionError @infinite_variable(m2, i, Int)
        @test_throws AssertionError @infinite_variable(m2, i(t), Bin)
        # test double specification
        @test_macro_throws ErrorException @infinite_variable(m, i(t),
                                                             parameter_refs = x)
        # test undefined parameter error
        @test_macro_throws ErrorException @infinite_variable(m, i, Int)
        # test invalid keyword arguments
        @test_macro_throws ErrorException @infinite_variable(m, i(t),
                                                           parameter_values = 1)
        @test_macro_throws ErrorException @infinite_variable(m, i(t),
                              infinite_variable_ref = InfiniteVariableRef(m, 1))
        @test_macro_throws ErrorException @infinite_variable(m, i(t), bad = 42)
        # test name duplication
        @test_macro_throws ErrorException @infinite_variable(m, a(t), Int)
    end
end

# Test the point variable macro
@testset "Point" begin
    # initialize model, parameters, and infinite variables
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= t <= 1)
    @infinite_parameter(m, -1 <= x[1:2] <= 1)
    @infinite_variable(m, 0 <= z(t, x) <= 1, Int)
    @infinite_variable(m, z2[1:2](t) == 3)
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        vref = PointVariableRef(m, 4)
        @test @point_variable(m, infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_integer(vref)
        @test lower_bound(vref) == 0
        # test anon with changes to fixed
        vref = PointVariableRef(m, 5)
        @test @point_variable(m, infinite_variable_ref = z, lower_bound = -5,
                          parameter_values = (0, [0, 0]), binary = true) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test !is_integer(vref)
        @test is_binary(vref)
        @test lower_bound(vref) == -5
        # test regular with alias
        vref = PointVariableRef(m, 6)
        @test @point_variable(m, z(0, [0, 0]), z0, Bin) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
        # test regular with semi anon
        vref = PointVariableRef(m, 7)
        @test @point_variable(m, z(0, [0, 0]), base_name = "z0",
                              binary = true) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
    end
    # test array variable definition
    @testset "Array" begin
        # test anon array with one infvar
        vrefs = [PointVariableRef(m, 8), PointVariableRef(m, 9)]
        @test @point_variable(m, [1:2], infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_integer(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        # test anon array with different inf vars
        vrefs = [PointVariableRef(m, 10), PointVariableRef(m, 11)]
        @test @point_variable(m, [i = 1:2], infinite_variable_ref = z2[i],
                              parameter_values = 0) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test fix_value(vrefs[2]) == 3
        @test name(vrefs[1]) == "z2[1](0)"
        # test array with same infvar
        vrefs = [PointVariableRef(m, 12), PointVariableRef(m, 13)]
        @test @point_variable(m, z(0, [0, 0]), a[1:2], Bin) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_binary(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        @test name(vrefs[1]) == "a[1]"
        # test test array with differnt infvars
        vrefs = [PointVariableRef(m, 14), PointVariableRef(m, 15)]
        @test @point_variable(m, z2[i](0), b[i = 1:2] >= -5) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "b[1]"
        # test semi anon array
        vrefs = [PointVariableRef(m, 16), PointVariableRef(m, 17)]
        @test @point_variable(m, z2[i](0), [i = 1:2], lower_bound = -5) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "z2[1](0)"
    end
    # test errors
    @testset "Errors" begin
        # test model assertion errors
        m2 = Model()
        @test_throws AssertionError @point_variable(m2, infinite_variable_ref = z,
                                                 parameter_values = (0, [0, 0]))
        @test_throws AssertionError @point_variable(m2, z(0, [0, 0]), bob, Bin)
        @test_throws AssertionError @point_variable(m2, [1:2],
                      infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        # test double specification
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), bob,
                                                      infinite_variable_ref = z)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), bob,
                                                 parameter_values = (0, [0, 0]))
        # test the adding expressions to infvar
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) >= 0,
                                                          bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) == 0,
                                                          bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) in 0,
                                                          bob, Bin)
        # test other syntaxes
        @test_macro_throws ErrorException @point_variable(m,
                                               0 <= z(0, [0, 0]) <= 1, bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, [1:2], Bin,
                      infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        @test_macro_throws ErrorException @point_variable(m, bob,
                     infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        # test redefinition catch
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), z0,
                                                          Bin)
    end
end

# Test the hold variable macro if no bounds are provided
# The bounded case is tested after @BDconstraint in constraints.jl
@testset "Hold" begin
    # initialize model
    m = InfiniteModel()
    @testset "@hold_variable (Not Bounded)" begin
        # test regular
        vref = HoldVariableRef(m, 1)
        @test @hold_variable(m, x >= 1, Bin) == vref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test anan
        vref = HoldVariableRef(m, 2)
        @test @hold_variable(m, binary = true, lower_bound = 1,
                               base_name = "x") == vref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test array
        vrefs = [HoldVariableRef(m, 3), HoldVariableRef(m, 4)]
        @test @hold_variable(m, y[1:2] == 2, Int) == vrefs
        @test name(vrefs[1]) == "y[1]"
        @test fix_value(vrefs[2]) == 2
        @test is_integer(vrefs[1])
        # test errors
        @test_throws AssertionError @hold_variable(Model(), z >= 1, Bin)
        @test_macro_throws ErrorException @hold_variable(m, x >= 1, Bin)
    end
    # test _make_interval_set
    @testset "_make_interval_set" begin
        expected = :(IntervalSet(0, 1))
        @test InfiniteOpt._make_interval_set(error, :([0, 1])) == expected
        @test_throws ErrorException InfiniteOpt._make_interval_set(error, :([0]))
    end
    # test _process_parameter
    @testset "_process_parameter" begin
        @test InfiniteOpt._process_parameter(:x) == :(InfiniteOpt._make_vector(x))
    end
    # test _make_bound_pair (in)
    @testset "_make_bound_pair (in)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
    end
    # test _make_bound_pair (==)
    @testset "_make_bound_pair (==)" begin
        expr = :(t == 0)
        set = :(IntervalSet(0, 0))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:(==))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with <=)
    @testset "_make_bound_pair (<=)" begin
        expr = :(0 <= t <= 1)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with >=)
    @testset "_make_bound_pair (>=)" begin
        expr = :(1 >= t >= 0)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(>=)),
                                           Val(:(>=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison fallback)
    @testset "_make_bound_pair (Fallback)" begin
        expr = :(t <= 0)
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                                 Val(:(<=)))
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                         Val(:(<=)), Val(:(>=)))
    end
    # test _make_bound_pair (wrapper)
    @testset "_make_bound_pair (Wrapper)" begin
        # test :call expression
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
        # test :compariosn expression
        expr = :(0 <= t <= 1)
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
        # test error
        expr = :([2, 3])
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr)
    end
    # test _parse_parameter_bounds (Vector)
    @testset "_parse_parameter_bounds (Vector)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg), $(dict_arg)))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, args) == bounds
    end
    # test _parse_parameter_bounds (Expr)
    @testset "_parse_parameter_bounds (Expr)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg)))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, expr) == bounds
    end
    # test _extract_bounds (:call)
    @testset "_extract_bounds (:call)" begin
        # test single anonymous bound
        args = [:in, :t, :([0, 1])]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (nothing, bounds)
    end
    # test _extract_bounds (:tuple)
    @testset "_extract_bounds (:tuple)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg), $(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:tuple)) == (nothing, bounds)
    end
    # test _extract_bounds (:comparison)
    @testset "_extract_bounds (:comparison)" begin
        args = [0, :(<=), :t, :(<=), 1]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:comparison)) == (nothing, bounds)
    end
    # test @hold_variable with bounds
    @testset "@hold_variable (Bounded)" begin
        # prepare parameters
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        # test regular call
        vref = HoldVariableRef(m, 5)
        @test @hold_variable(m, xb >= 1, Bin,
                             parameter_bounds = (par == 0)) == vref
        @test name(vref) == "xb"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        @test m.has_hold_bounds
        @test m.vars[5].parameter_bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test supports(par) == [0]
        # test regular tuple
        vref = HoldVariableRef(m, 6)
        @test @hold_variable(m, yb, parameter_bounds = (pars == 0,
                                                       par in [0, 1])) == vref
        @test name(vref) == "yb"
        dict = Dict(pars[1] => IntervalSet(0, 0), pars[2] => IntervalSet(0, 0),
                    par => IntervalSet(0, 1))
        @test m.vars[6].parameter_bounds == ParameterBounds(dict)
        @test supports(pars[1]) == [0]
        # test comparison
        vref = HoldVariableRef(m, 7)
        @test @hold_variable(m, zb, parameter_bounds = (0 <= par <= 1),
                             base_name = "bob") == vref
        @test name(vref) == "bob"
        @test m.vars[7].parameter_bounds == ParameterBounds(Dict(par => IntervalSet(0, 1)))
        # test unrecognized format
        @test_macro_throws ErrorException @hold_variable(m, parameter_bounds = par)
        # test container specification
        vrefs = [HoldVariableRef(m, 8), HoldVariableRef(m, 9)]
        vrefs = convert(JuMPC.SparseAxisArray, vrefs)
        @test @hold_variable(m, [i = 1:2], parameter_bounds = (pars[i] == 0),
                             container = SparseAxisArray) == vrefs
        @test name(vrefs[1]) == ""
        @test m.vars[9].parameter_bounds == ParameterBounds(Dict(pars[2] => IntervalSet(0, 0)))
        # test wrong model type
        @test_macro_throws ErrorException @hold_variable(Model(),
                                                  parameter_bounds = (par == 0))
    end
end
