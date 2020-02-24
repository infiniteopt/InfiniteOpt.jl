# Test _arg_promote
@testset "_arg_promote" begin
    # setup the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_variable(m, inf(par))
    @hold_variable(m, hold)
    # test GenericAffExpr
    @testset "AffExpr" begin
        aff = 2inf + 2
        @test isa(InfiniteOpt._arg_promote(GeneralVariableRef, aff), GenericAffExpr{Float64, GeneralVariableRef})
        @test InfiniteOpt._arg_promote(InfiniteVariableRef, aff) == aff
    end
    # test GenericQuadExpr
    @testset "QuadExpr" begin
        quad = inf^2 + 2inf
        @test isa(InfiniteOpt._arg_promote(GeneralVariableRef, quad), GenericQuadExpr{Float64, GeneralVariableRef})
        @test InfiniteOpt._arg_promote(InfiniteVariableRef, quad) == quad
    end
    # test Variable references
    @testset "Variables" begin
        @test InfiniteOpt._arg_promote(GeneralVariableRef, inf) == inf + 0
        @test InfiniteOpt._arg_promote(GeneralVariableRef, par) == par + 0
        @test InfiniteOpt._arg_promote(GeneralVariableRef, hold) == hold + 0
    end
end

# Test MutableArithmetics.mutable_operate!
@testset "MA.mutable_operate!" begin
    # setup the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_variable(m, inf(par))
    @hold_variable(m, hold)
    aff1 = 2inf + 3
    aff2 = hold - par + 2
    quad1 = hold^2 + 4hold
    quad2 = 2par^2 - par - 2
    # 2 argument adding
    @testset "+" begin
        # mixed aff + var
        expected = copy(aff1) + hold
        @test MA.mutable_operate!(+, copy(aff1), hold) == expected
        # mixed aff + aff
        expected = copy(aff1) + copy(aff2)
        @test MA.mutable_operate!(+, copy(aff1), copy(aff2)) == expected
        # same aff + var
        expected = copy(aff1) + inf
        @test MA.mutable_operate!(+, copy(aff1), inf) == expected
        # mixed quad + var
        expected = copy(quad1) + par
        @test MA.mutable_operate!(+, copy(quad1), par) == expected
        # mixed quad + aff
        expected = copy(quad2) + copy(aff2)
        @test MA.mutable_operate!(+, copy(quad2), copy(aff2)) == expected
        # test quad + constant (not this function)
        expected = copy(quad2) + 2
        @test MA.mutable_operate!(+, copy(quad2), 2) == expected
    end

    # 2 argument substracting
    @testset "-" begin
        # mixed aff + var
        expected = copy(aff1) - hold
        @test MA.mutable_operate!(-, copy(aff1), hold) == expected
        # mixed aff + aff
        expected = copy(aff1) - copy(aff2)
        @test MA.mutable_operate!(-, copy(aff1), copy(aff2)) == expected
        # same aff + var
        expected = copy(aff1) - inf
        @test MA.mutable_operate!(-, copy(aff1), inf) == expected
        # mixed quad + var
        expected = copy(quad1) - par
        @test MA.mutable_operate!(-, copy(quad1), par) == expected
        # mixed quad + aff
        expected = copy(quad2) - copy(aff2)
        @test MA.mutable_operate!(-, copy(quad2), copy(aff2)) == expected
        # test quad + constant (not this function)
        expected = copy(quad2) - 2
        @test MA.mutable_operate!(-, copy(quad2), 2) == expected
    end

    # 2 argument adding/multiplying
    @testset "MA.add_mul" begin
        # mixed aff + var
        expected = copy(aff1) + hold
        @test MA.mutable_operate!(MA.add_mul, copy(aff1), hold) == expected
        # mixed aff + aff
        expected = copy(aff1) + copy(aff2)
        @test MA.mutable_operate!(MA.add_mul, copy(aff1), copy(aff2)) == expected
        # same aff + var
        expected = copy(aff1) + inf
        @test MA.mutable_operate!(MA.add_mul, copy(aff1), inf) == expected
        # mixed quad + var
        expected = copy(quad1) + par
        @test MA.mutable_operate!(MA.add_mul, copy(quad1), par) == expected
        # mixed quad + aff
        expected = copy(quad2) + copy(aff2)
        @test MA.mutable_operate!(MA.add_mul, copy(quad2), copy(aff2)) == expected
        # test quad + constant (not this function)
        expected = copy(quad2) + 2
        @test MA.mutable_operate!(MA.add_mul, copy(quad2), 2) == expected
    end

    # 2 argument substracting/multiplying
    @testset "MA.sub_mul" begin
        # mixed aff + var
        expected = copy(aff1) - hold
        @test MA.mutable_operate!(MA.sub_mul, copy(aff1), hold) == expected
        # mixed aff + aff
        expected = copy(aff1) - copy(aff2)
        @test MA.mutable_operate!(MA.sub_mul, copy(aff1), copy(aff2)) == expected
        # same aff + var
        expected = copy(aff1) - inf
        @test MA.mutable_operate!(MA.sub_mul, copy(aff1), inf) == expected
        # mixed quad + var
        expected = copy(quad1) - par
        @test MA.mutable_operate!(MA.sub_mul, copy(quad1), par) == expected
        # mixed quad + aff
        expected = copy(quad2) - copy(aff2)
        @test MA.mutable_operate!(MA.sub_mul, copy(quad2), copy(aff2)) == expected
        # test quad + constant (not this function)
        expected = copy(quad2) - 2
        @test MA.mutable_operate!(MA.sub_mul, copy(quad2), 2) == expected
    end

    # 3 argument adding/multiplying (no constants)
    @testset "MA.add_mul (3 args no constants)" begin
        # mixed quad + var * var
        expected = copy(quad2) + hold * par
        @test MA.mutable_operate!(MA.add_mul, copy(quad2), hold, par) == expected
        # same quad + var * var
        expected = copy(quad1) + hold * hold
        @test MA.mutable_operate!(MA.add_mul, copy(quad1), hold, hold) == expected
    end

    # 3 argument adding/multiplying (terminal constant)
    @testset "MA.add_mul (3 args last constant)" begin
        # mixed quad + var * const
        expected = copy(quad2) + hold * 2
        @test MA.mutable_operate!(MA.add_mul, copy(quad2), hold, 2) == expected
        # same quad + var * var
        expected = copy(quad1) + hold * 3
        @test MA.mutable_operate!(MA.add_mul, copy(quad1), hold, 3) == expected
    end

    # 3 argument adding/multiplying (middle constant)
    @testset "MA.add_mul (3 args middle constant)" begin
        # mixed quad + 2 * var
        expected = copy(quad2) + 2 * par
        @test MA.mutable_operate!(MA.add_mul, copy(quad2), 2, par) == expected
        # same quad + 42 * var
        expected = copy(quad1) + 42 * hold
        @test MA.mutable_operate!(MA.add_mul, copy(quad1), 42, hold) == expected
        # all constants (different function)
        expected = copy(quad1) + 42 * 2
        @test MA.mutable_operate!(MA.add_mul, copy(quad1), 42, 2) == expected
    end

    # 3 argument substracting/multiplying (no constants)
    @testset "MA.sub_mul (3 args no constants)" begin
        # mixed quad + var * var
        expected = copy(quad2) - hold * par
        @test MA.mutable_operate!(MA.sub_mul, copy(quad2), hold, par) == expected
        # same quad + var * var
        expected = copy(quad1) - hold * hold
        @test MA.mutable_operate!(MA.sub_mul, copy(quad1), hold, hold) == expected
    end

    # 3 argument substracting/multiplying (terminal constant)
    @testset "MA.sub_mul (3 args last constant)" begin
        # mixed quad + var * const
        expected = copy(quad2) - hold * 2
        @test MA.mutable_operate!(MA.sub_mul, copy(quad2), hold, 2) == expected
        # same quad + var * var
        expected = copy(quad1) - hold * 3
        @test MA.mutable_operate!(MA.sub_mul, copy(quad1), hold, 3) == expected
    end

    # 3 argument substracting/multiplying (middle constant)
    @testset "MA.sub_mul (3 args middle constant)" begin
        # mixed quad + 2 * var
        expected = copy(quad2) - 2 * par
        @test MA.mutable_operate!(MA.sub_mul, copy(quad2), 2, par) == expected
        # same quad + 42 * var
        expected = copy(quad1) - 42 * hold
        @test MA.mutable_operate!(MA.sub_mul, copy(quad1), 42, hold) == expected
        # all constants (different function)
        expected = copy(quad1) - 42 * 2
        @test MA.mutable_operate!(MA.sub_mul, copy(quad1), 42, 2) == expected
    end
end
