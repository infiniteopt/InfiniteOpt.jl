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

end
