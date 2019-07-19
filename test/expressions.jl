# Test _all_function_variables
@testset "_all_function_variables" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    # test for variable reference
    @testset "Variable" begin
        @test InfiniteOpt._all_function_variables(par) == [par]
        @test InfiniteOpt._all_function_variables(inf) == [inf]
        @test InfiniteOpt._all_function_variables(meas) == [meas]
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = meas + 2par + glob
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test InfiniteOpt._all_function_variables(aff1) == [meas, par, glob]
        @test InfiniteOpt._all_function_variables(aff2) == GeneralVariableRef[]
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt - meas + 2par + glob
        quad2 = pt^2 + inf * pt
        quad3 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        # test expressions
        @test InfiniteOpt._all_function_variables(quad1) == [meas, par, glob, pt, inf]
        @test InfiniteOpt._all_function_variables(quad2) == [pt, inf]
        @test InfiniteOpt._all_function_variables(quad3) == GeneralVariableRef[]
    end
end

# Test comparisons
@testset "Comparisons" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    red = InfiniteOpt._ReducedInfiniteRef(m, index(inf2), inf2, Dict(1 => 3))
    # test AffExpr comparison
    @testset "Base.:(==) AffExpr" begin
        @test par + par2 + inf - 2 == par + (par2 + inf) - 2
        @test 0.25red - inf == 0.25red - inf
        @test red + 3 - inf != par + red
        @test red + 3 - inf != red + 3 - glob
    end
    # test QuadExpr comparison
    @testset "Base.:(==) QuadExpr" begin
        @test par * inf + par2 + inf - 2 == par * inf + (par2 + inf) - 2
        @test red * inf - inf == red * inf - inf
        @test red * inf + 3 - inf != par * inf2 + red
        @test par * par2 + inf * inf2 == par * par2 + inf * inf2
        @test par * par2 + inf * inf2 != par * par2
        @test par * par2 + 2 * inf * inf2 !=  par * par2 + inf * inf2
        @test par * inf + par2 + inf - 2 != par * inf + (par2 + inf) - 3
    end
end

# Test _all_parameter_refs
@testset "_all_parameter_refs" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    red = InfiniteOpt._ReducedInfiniteRef(m, index(inf2), inf2, Dict(1 => 3))
    # test for finite variable reference
    @testset "FiniteVariable" begin
        @test InfiniteOpt._all_parameter_refs(pt) == ()
        @test InfiniteOpt._all_parameter_refs(glob) == ()
    end
    # test for infinite variable reference
    @testset "InfiniteVariable" begin
        @test InfiniteOpt._all_parameter_refs(inf) == (par, )
    end
    # test for parameter reference
    @testset "Parameter" begin
        @test InfiniteOpt._all_parameter_refs(par) == (par, )
    end
    # test for reduced infinite variable reference
    @testset "ReducedInfinite" begin
        @test InfiniteOpt._all_parameter_refs(red) == (par2, )
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = inf + inf2 + pt - 3
        aff2 = pt + glob - 2
        # test expressions
        @test InfiniteOpt._all_parameter_refs(aff1) == (par, par2)
        @test InfiniteOpt._all_parameter_refs(aff2) == ()
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = inf * inf2 + inf + inf2 + pt - 3 - par
        quad2 = pt * pt + pt + glob - 2
        # test expressions
        @test InfiniteOpt._all_parameter_refs(quad1) == (par, par2)
        @test InfiniteOpt._all_parameter_refs(quad2) == ()
    end
end

# Test _remove_variable
@testset "_remove_variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = meas + 2par + glob
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isa(InfiniteOpt._remove_variable(aff1, glob), Nothing)
        @test !haskey(aff1.terms, glob)
        @test isa(InfiniteOpt._remove_variable(aff2, inf), Nothing)
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt - meas + 2par + glob
        quad2 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        # test expressions
        quad = copy(quad1)
        @test isa(InfiniteOpt._remove_variable(quad, inf), Nothing)
        @test !haskey(quad.aff.terms, inf)
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(inf, pt))
        quad = copy(quad1)
        @test isa(InfiniteOpt._remove_variable(quad, pt), Nothing)
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(inf, pt))
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(pt, pt))
        @test isa(InfiniteOpt._remove_variable(quad2, inf), Nothing)
    end
end
