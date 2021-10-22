# Test operations between 2 different variables
@testset "Variable--Variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    dinf = @deriv(inf, par)
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    # test addition
    @testset "Base.:+" begin
        # test infinite + measure
        @test isa(inf + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (inf + meas).constant == 0.0
        @test (inf + meas).terms[meas] == 1
        @test (inf + meas).terms[inf] == 1
        # test point + finite
        @test isa(pt + finite, GenericAffExpr{Float64, GeneralVariableRef})
        @test (pt + finite).terms[finite] == 1
        # test other combos
        @test isa(pt + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + pt, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + meas, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (meas + meas).terms[meas] == 2
        @test isa(inf + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(pt + pt, GenericAffExpr{Float64, GeneralVariableRef})
        # test with derivative 
        @test isa(dinf + dinf, GenericAffExpr{Float64, GeneralVariableRef})
        @test (dinf + dinf).terms[dinf] == 2
        @test isa(dinf + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(dinf + pt, GenericAffExpr{Float64, GeneralVariableRef})
    end
    # test subtraction
    @testset "Base.:-" begin
        # test finite - measure
        @test isa(finite - meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (finite - meas).constant == 0.0
        @test (finite - meas).terms[meas] == -1
        @test (finite - meas).terms[finite] == 1
        # test point - finite
        @test isa(pt - finite, GenericAffExpr{Float64, GeneralVariableRef})
        @test (pt - finite).terms[finite] == -1
        # test other combos
        @test isa(pt - meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - pt, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - meas, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas - meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test isequal((meas - meas), zero(GenericAffExpr{Float64, GeneralVariableRef}))
        @test isa(inf - inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(pt - pt, GenericAffExpr{Float64, GeneralVariableRef})
    end
    # test multiplication
    @testset "Base.:*" begin
        # test point * measure
        @test isa(pt * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (pt * meas).aff.constant == 0.0
        pair = UnorderedPair{GeneralVariableRef}(pt, meas)
        @test (pt * meas).terms[pair] == 1
        # test point * finite
        @test isa(pt * finite, GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (pt * finite).terms[pair] == 1
        # test other combos
        @test isa(pt * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(meas, meas)
        @test (meas * meas).terms[pair] == 1
        @test isa(inf * inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(pt * pt, GenericQuadExpr{Float64, GeneralVariableRef})
    end
end

# Test operations between variables and GenericAffExpr
@testset "Variable--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    # test addition
    @testset "Base.:+" begin
        # test infinite + aff1
        @test isa(inf + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test (inf + aff1).constant == -2
        @test (inf + aff1).terms[pt] == 1
        @test (inf + aff1).terms[inf] == 2
        # test meas + aff2
        @test isa(meas + aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test (meas + aff2).terms[finite] == 1
        # test other combos
        @test isa(finite + aff3, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + aff3, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas + (meas + meas), GenericAffExpr{Float64, GeneralVariableRef})
        @test (meas + (meas + meas)).terms[meas] == 3
    end
    # test subtraction
    @testset "Base.:-" begin
        # test infinite - aff1
        @test isa(inf - aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test (inf - aff1).constant == 2
        @test (inf - aff1).terms[pt] == -1
        @test (inf - aff1).terms[inf] == 0
        # test meas + aff2
        @test isa(meas - aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test (meas - aff2).terms[finite] == -1
        # test other combos
        @test isa(finite - aff3, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - aff3, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas - (meas + meas), GenericAffExpr{Float64, GeneralVariableRef})
        @test (meas - (meas + meas)).terms[meas] == -1
    end
    # test multiplication
    @testset "Base.:*" begin
        # test infinite * aff1
        @test isa(inf * aff1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (inf * aff1).aff.constant == 0
        @test (inf * aff1).aff.terms[inf] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test (inf * aff1).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test (inf * aff1).terms[pair] == 1
        # test meas * aff2
        @test isa(meas * aff2, GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(meas, finite)
        @test (meas * aff2).terms[pair] == 1
        @test length((meas * aff2).aff.terms) == 0
        # test other combos
        @test isa(finite * aff3, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * aff2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * aff1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * aff3, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas * (meas + meas), GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(meas, meas)
        @test (meas * (meas + meas)).terms[pair] == 2
    end
end

# Test operations between GenericAffExpr and variable
@testset "AffExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    # test addition
    @testset "Base.:+" begin
        # test aff1 + inf
        @test isa(copy(aff1) + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff1) + inf).constant == -2
        @test (copy(aff1) + inf).terms[pt] == 1
        @test (copy(aff1) + inf).terms[inf] == 2
        # test aff2 + meas
        @test isa(copy(aff2) + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff2) + meas).terms[finite] == 1
        # test other combos
        @test isa(aff3 + finite, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff2 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff1 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff3 + par, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas) + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test ((meas + meas) + meas).terms[meas] == 3
    end
    # test subtraction
    @testset "Base.:-" begin
        # test aff1 - inf
        @test isa(copy(aff1) - inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff1) - inf).constant == -2
        @test (copy(aff1) - inf).terms[pt] == 1
        @test (copy(aff1) - inf).terms[inf] == 0
        # test aff2 - meas
        @test isa(copy(aff2) - meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff2) - meas).terms[finite] == 1
        # test other combos
        @test isa(aff3 - finite, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff2 - par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff1 - par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff3 - par, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas) - meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test ((meas + meas) - meas).terms[meas] == 1
    end
    # test multiplication
    @testset "Base.:*" begin
        # test aff1 * inf
        @test isa(aff1 * inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (aff1 * inf).aff.constant == 0
        @test (aff1 * inf).aff.terms[inf] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test (aff1 * inf).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test (aff1 * inf).terms[pair] == 1
        # test aff2 * meas
        @test isa(aff2 * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(meas, finite)
        @test (aff2 * meas).terms[pair] == 1
        @test length((aff2 * meas).aff.terms) == 0
        # test other combos
        @test isa(aff3 * finite, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff2 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff1 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff3 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas * meas), GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(meas, meas)
        @test ((meas + meas) * meas).terms[pair] == 2
    end
end

# Test operators on GenericAffExpr and GenericAffExpr
@testset "AffExpr--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff4 = (inf + inf) + 3
    aff5 = 4pt - 42
    aff_big = sum(test[i] for i = 1:length(test))
    # test addition
    @testset "Base.:+" begin
        # test large input
        @test isa(aff_big + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff1) + aff_big, GenericAffExpr{Float64, GeneralVariableRef})
        # test various mixed combination types
        @test isa(copy(aff1) + aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff2 + aff3, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) + aff5, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff4 + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test (copy(aff1) + aff2).constant == -2
        @test (copy(aff1) + aff2).terms[inf] == 1
        @test (copy(aff1) + aff2).terms[pt] == 0
        @test (aff4 + aff1).constant == 1
        @test (aff4 + aff1).terms[inf] == 3
        @test (aff4 + aff1).terms[pt] == 1
        # test results of measurefinite combinations
        @test (aff2 + aff3).constant == 1
        @test (aff2 + aff3).terms[meas] == -1
        @test (aff2 + aff3).terms[finite] == 2
        @test (copy(aff3) + aff5).constant == -41
        @test (copy(aff3) + aff5).terms[meas] == -1
        @test (copy(aff3) + aff5).terms[pt] == 4
        # test results of fininite combinations
        @test (copy(aff2) + aff5).constant == -42
        @test (copy(aff2) + aff5).terms[pt] == 3
        @test (copy(aff2) + aff5).terms[finite] == 1
        # test pure additions
        @test isa(copy(aff4) + aff4, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff4) + aff4).constant == 6
        @test (copy(aff4) + aff4).terms[inf] == 4
        @test isa(copy(aff5) + aff5, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff5) + aff5).constant == -84
        @test (copy(aff5) + aff5).terms[pt] == 8
    end
    # test subtraction
    @testset "Base.:-" begin
        # test various mixed combination types
        @test isa(copy(aff1) - aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff2 - aff3, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) - aff5, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff4 - aff1, GenericAffExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test (copy(aff1) - aff2).constant == -2
        @test (copy(aff1) - aff2).terms[inf] == 1
        @test (copy(aff1) - aff2).terms[pt] == 2
        @test (aff4 - aff1).constant == 5
        @test (aff4 - aff1).terms[inf] == 1
        @test (aff4 - aff1).terms[pt] == -1
        # test results of measurefinite combinations
        @test (aff2 - aff3).constant == -1
        @test (aff2 - aff3).terms[meas] == 1
        @test (aff2 - aff3).terms[finite] == 0
        @test (copy(aff3) - aff5).constant == 43
        @test (copy(aff3) - aff5).terms[meas] == -1
        @test (copy(aff3) - aff5).terms[pt] == -4
        # test results of fininite combinations
        @test (copy(aff2) - aff5).constant == 42
        @test (copy(aff2) - aff5).terms[pt] == -5
        @test (copy(aff2) - aff5).terms[finite] == 1
        # test pure additions
        @test isa(copy(aff4) - aff4, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff4) - aff4).constant == 0
        @test (copy(aff4) - aff4).terms[inf] == 0
        @test isa(copy(aff5) - aff5, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff5) - aff5).constant == 0
        @test (copy(aff5) - aff5).terms[pt] == 0
    end
    # test multiplication
    @testset "Base.:*" begin
        # test various mixed combination types
        @test isa(aff1 * aff2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff2 * aff3, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff2 * aff5, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff4 * aff1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test (aff1 * aff2).aff.constant == 0
        @test (aff1 * aff2).aff.terms[finite] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (aff1 * aff2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (aff1 * aff2).terms[pair] == -1
        @test (aff4 * aff1).aff.constant == -6
        @test (aff4 * aff1).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test (aff4 * aff1).terms[pair] == 2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test (aff4 * aff1).terms[pair] == 2
        # test results of measurefinite combinations
        @test (aff2 * aff3).aff.constant == 0
        @test (aff2 * aff3).aff.terms[pt] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, meas)
        @test (aff2 * aff3).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, finite)
        @test (aff2 * aff3).terms[pair] == 1
        @test (aff3 * aff5).aff.constant == -42
        @test (aff3 * aff5).aff.terms[meas] == 42
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (aff3 * aff5).terms[pair] == 4
        pair = UnorderedPair{GeneralVariableRef}(meas, pt)
        @test (aff3 * aff5).terms[pair] == -4
        # test results of fininite combinations
        @test (aff2 * aff5).aff.constant == 0
        @test (aff2 * aff5).aff.terms[finite] == -42
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (aff2 * aff5).terms[pair] == 4
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (aff2 * aff5).terms[pair] == -4
        # test pure additions
        @test isa(aff4 * aff4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (aff4 * aff4).aff.constant == 9
        @test (aff4 * aff4).aff.terms[inf] == 12
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test (aff4 * aff4).terms[pair] == 4
        @test isa(aff5 * aff5, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (aff5 * aff5).aff.constant == 42 * 42
        @test (aff5 * aff5).aff.terms[pt] == 2 * 4 * -42
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (aff5 * aff5).terms[pair] == 16
    end
end

# Test operations for GenericQuadExpr--GeneralVariableRef
@testset "QuadExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*finite - inf*pt + pt*finite - pt^2 - 2finite + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*finite - pt*meas + pt - 3)
    quad3 = finite * aff2       #(finite^2 - finite*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + meas, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) + pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad4) + par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) + inf).aff.constant == 0
        @test (copy(quad1) + inf).aff.terms[pt] == 2
        @test (copy(quad1) + inf).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) + inf).terms[pair] == 1
        @test (copy(quad4) + par).aff.constant == 0
        @test (copy(quad4) + par).aff.terms[pt] == -42
        @test (copy(quad4) + par).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + par).terms[pair] == 4
        # test measurefinite combos
        @test (copy(quad2) + meas).aff.constant == -3
        @test (copy(quad2) + meas).aff.terms[pt] == 1
        @test (copy(quad2) + meas).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) + meas).terms[pair] == 1
        @test (copy(quad2) + pt).aff.constant == -3
        @test (copy(quad2) + pt).aff.terms[pt] == 2
        # test finite combos
        @test (copy(quad3) + pt).aff.constant == 0
        @test (copy(quad3) + pt).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) + pt).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) + pt).aff.constant == 0
        @test (copy(quad4) + pt).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + pt).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - meas, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) - pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad4) - par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - inf).aff.constant == 0
        @test (copy(quad1) - inf).aff.terms[pt] == 2
        @test (copy(quad1) - inf).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) - inf).terms[pair] == 1
        @test (copy(quad4) - par).aff.constant == 0
        @test (copy(quad4) - par).aff.terms[pt] == -42
        @test (copy(quad4) - par).aff.terms[par] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - par).terms[pair] == 4
        # test measurefinite combos
        @test (copy(quad2) - meas).aff.constant == -3
        @test (copy(quad2) - meas).aff.terms[pt] == 1
        @test (copy(quad2) - meas).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) - meas).terms[pair] == 1
        @test (copy(quad2) - pt).aff.constant == -3
        @test (copy(quad2) - pt).aff.terms[pt] == 0
        # test finite combos
        @test (copy(quad3) - pt).aff.constant == 0
        @test (copy(quad3) - pt).aff.terms[pt] == -1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) - pt).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) - pt).aff.constant == 0
        @test (copy(quad4) - pt).aff.terms[pt] == -43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - pt).terms[pair] == 4
    end
end

# Test operations for GeneralVariableRef--GenericQuadExpr
@testset "Variable--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*finite - inf*pt + pt*finite - pt^2 - 2finite + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*finite - pt*meas + pt - 3)
    quad3 = finite * aff2       #(finite^2 - finite*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(inf + copy(quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(meas + copy(quad2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(pt + copy(quad3), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par + copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (inf + copy(quad1)).aff.constant == 0
        @test (inf + copy(quad1)).aff.terms[pt] == 2
        @test (inf + copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (inf + copy(quad1)).terms[pair] == 1
        @test (par + copy(quad4)).aff.constant == 0
        @test (par + copy(quad4)).aff.terms[pt] == -42
        @test (par + copy(quad4)).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (par + copy(quad4)).terms[pair] == 4
        # test measurefinite combos
        @test (meas + copy(quad2)).aff.constant == -3
        @test (meas + copy(quad2)).aff.terms[pt] == 1
        @test (meas + copy(quad2)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (meas + copy(quad2)).terms[pair] == 1
        @test (pt + copy(quad2)).aff.constant == -3
        @test (pt + copy(quad2)).aff.terms[pt] == 2
        # test finite combos
        @test (pt + copy(quad3)).aff.constant == 0
        @test (pt + copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (pt + copy(quad3)).terms[pair] == -1
        # test same types
        @test isa(pt + copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        @test (pt + copy(quad4)).aff.constant == 0
        @test (pt + copy(quad4)).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (pt + copy(quad4)).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(inf - copy(quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(meas - copy(quad2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(pt - copy(quad3), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par - copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (inf - copy(quad1)).aff.constant == 0
        @test (inf - copy(quad1)).aff.terms[pt] == -2
        @test (inf - copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (inf - copy(quad1)).terms[pair] == -1
        @test (par - copy(quad4)).aff.constant == 0
        @test (par - copy(quad4)).aff.terms[pt] == 42
        @test (par - copy(quad4)).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (par - copy(quad4)).terms[pair] == -4
        # test measurefinite combos
        @test (meas - copy(quad2)).aff.constant == 3
        @test (meas - copy(quad2)).aff.terms[pt] == -1
        @test (meas - copy(quad2)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (meas - copy(quad2)).terms[pair] == -1
        @test (pt - copy(quad2)).aff.constant == 3
        @test (pt - copy(quad2)).aff.terms[pt] == 0
        # test finite combos
        @test (pt - copy(quad3)).aff.constant == 0
        @test (pt - copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (pt - copy(quad3)).terms[pair] == 1
        # test same types
        @test isa(pt - copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        @test (pt - copy(quad4)).aff.constant == 0
        @test (pt - copy(quad4)).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (pt - copy(quad4)).terms[pair] == -4
    end
end

# Test operators on GenericAffExpr--GenericQuadExpr
@testset "AffExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*finite - inf*pt + pt*finite - pt^2 - 2finite + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*finite - pt*meas + pt - 3)
    quad3 = finite * aff2       #(finite^2 - finite*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(aff1) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) + quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff5) + quad3, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff3) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(aff1) + quad1).aff.constant == -2
        @test (copy(aff1) + quad1).aff.terms[pt] == 3
        @test (copy(aff1) + quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(aff1) + quad1).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff1) + quad1).terms[pair] == -1
        @test (copy(aff3) + quad1).aff.constant == 1
        @test (copy(aff3) + quad1).aff.terms[meas] == -1
        @test (copy(aff3) + quad1).aff.terms[finite] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(aff3) + quad1).terms[pair] == 1
        # test measurefinite combos
        @test (copy(aff3) + quad2).aff.constant == -2
        @test (copy(aff3) + quad2).aff.terms[pt] == 1
        @test (copy(aff3) + quad2).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(aff3) + quad2).terms[pair] == 1
        @test (copy(aff3) + quad4).aff.constant == 1
        @test (copy(aff3) + quad4).aff.terms[pt] == -42
        @test (copy(aff3) + quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff3) + quad4).terms[pair] == 4
        # test finite combos
        @test (copy(aff2) + quad3).aff.constant == 0
        @test (copy(aff2) + quad3).aff.terms[pt] == -1
        @test (copy(aff2) + quad3).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(aff2) + quad3).terms[pair] == -1
        # test same types
        @test isa(copy(aff5) + quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(aff5) + quad4).aff.constant == -42
        @test (copy(aff5) + quad4).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff5) + quad4).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(aff1) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) - quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff5) - quad3, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff3) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(aff1) - quad1).aff.constant == -2
        @test (copy(aff1) - quad1).aff.terms[pt] == -1
        @test (copy(aff1) - quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(aff1) - quad1).terms[pair] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff1) - quad1).terms[pair] == 1
        @test (copy(aff3) - quad1).aff.constant == 1
        @test (copy(aff3) - quad1).aff.terms[meas] == -1
        @test (copy(aff3) - quad1).aff.terms[finite] == 3
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(aff3) - quad1).terms[pair] == -1
        # test measurefinite combos
        @test (copy(aff3) - quad2).aff.constant == 4
        @test (copy(aff3) - quad2).aff.terms[pt] == -1
        @test (copy(aff3) - quad2).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(aff3) - quad2).terms[pair] == -1
        @test (copy(aff3) - quad4).aff.constant == 1
        @test (copy(aff3) - quad4).aff.terms[pt] == 42
        @test (copy(aff3) - quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff3) - quad4).terms[pair] == -4
        # test finite combos
        @test (copy(aff2) - quad3).aff.constant == 0
        @test (copy(aff2) - quad3).aff.terms[pt] == -1
        @test (copy(aff2) - quad3).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(aff2) - quad3).terms[pair] == 1
        # test same types
        @test isa(copy(aff5) - quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(aff5) - quad4).aff.constant == -42
        @test (copy(aff5) - quad4).aff.terms[pt] == 46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff5) - quad4).terms[pair] == -4
    end
end

# Test operators on GenericQuadExpr--GenericAffExpr
@testset "QuadExpr--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*finite - inf*pt + pt*finite - pt^2 - 2finite + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*finite - pt*meas + pt - 3)
    quad3 = finite * aff2       #(finite^2 - finite*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + copy(aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + copy(aff2) + quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) + copy(aff5), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad1) + copy(aff3), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1)+ copy(aff1)).aff.constant == -2
        @test (copy(quad1) + copy(aff1)).aff.terms[pt] == 3
        @test (copy(quad1) + copy(aff1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) + copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad1) + copy(aff3)).terms[pair] == -1
        @test (copy(quad1) + copy(aff3)).aff.constant == 1
        @test (copy(quad1) + copy(aff3)).aff.terms[meas] == -1
        @test (copy(quad1) + copy(aff3)).aff.terms[finite] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) + copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test (copy(quad2) + copy(aff3)).aff.constant == -2
        @test (copy(quad2) + copy(aff3)).aff.terms[pt] == 1
        @test (copy(quad2) + copy(aff3)).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) + copy(aff3)).terms[pair] == 1
        @test (copy(quad4) + copy(aff3)).aff.constant == 1
        @test (copy(quad4) + copy(aff3)).aff.terms[pt] == -42
        @test (copy(quad4) + copy(aff3)).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + copy(aff3)).terms[pair] == 4
        # test finite combos
        @test (copy(quad3) + copy(aff2)).aff.constant == 0
        @test (copy(quad3) + copy(aff2)).aff.terms[pt] == -1
        @test (copy(quad3) + copy(aff2)).aff.terms[finite] == 1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) + copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + copy(aff5), GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) + copy(aff5)).aff.constant == -42
        @test (copy(quad4) + copy(aff5)).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + copy(aff5)).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - copy(aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - copy(aff2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) - copy(aff5), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad1) - copy(aff3), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - copy(aff1)).aff.constant == 2
        @test (copy(quad1) - copy(aff1)).aff.terms[pt] == 1
        @test (copy(quad1) - copy(aff1)).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) - copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad1) - copy(aff1)).terms[pair] == -1
        @test (copy(quad1) - copy(aff3)).aff.constant == -1
        @test (copy(quad1) - copy(aff3)).aff.terms[meas] == 1
        @test (copy(quad1) - copy(aff3)).aff.terms[finite] == -3
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) - copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test (copy(quad2) - copy(aff3)).aff.constant == -4
        @test (copy(quad2) - copy(aff3)).aff.terms[pt] == 1
        @test (copy(quad2) - copy(aff3)).aff.terms[finite] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) - copy(aff3)).terms[pair] == 1
        @test (copy(quad4) - copy(aff3)).aff.constant == -1
        @test (copy(quad4) - copy(aff3)).aff.terms[pt] == -42
        @test (copy(quad4) - copy(aff3)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4)- copy(aff3)).terms[pair] == 4
        # test finite combos
        @test (copy(quad3) - copy(aff2)).aff.constant == 0
        @test (copy(quad3) - copy(aff2)).aff.terms[pt] == 1
        @test (copy(quad3) - copy(aff2)).aff.terms[finite] == -1
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) - copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - copy(aff5), GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) - copy(aff5)).aff.constant == 42
        @test (copy(quad4) - copy(aff5)).aff.terms[pt] == -46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - copy(aff5)).terms[pair] == 4
    end
end

# Test operators on GenericQuadExpr--GenericQuadExpr
@testset "QuadExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    @variable(m, test[1:51])
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    aff1 = (inf + pt) - 2
    aff2 = finite - pt
    aff3 = (finite - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*finite - inf*pt + pt*finite - pt^2 - 2finite + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*finite - pt*meas + pt - 3)
    quad3 = finite * aff2       #(finite^2 - finite*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) + quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad4) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) + quad2).aff.constant == -3
        @test (copy(quad1) + quad2).aff.terms[pt] == 3
        @test (copy(quad1) + quad2).aff.terms[finite] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) + quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad1) + quad2).terms[pair] == 2
        @test (copy(quad4) + quad1).aff.constant == 0
        @test (copy(quad4) + quad1).aff.terms[pt] == -40
        @test (copy(quad4) + quad1).aff.terms[finite] == -2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + quad1).terms[pair] == 3
        # test measurefinite combos
        @test (copy(quad2) + quad4).aff.constant == -3
        @test (copy(quad2) + quad4).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) + quad4).terms[pair] == 1
        # test finite combos
        @test (copy(quad3) + quad4).aff.constant == 0
        @test (copy(quad3) + quad4).aff.terms[pt] == -42
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) + quad4).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) + quad4).aff.constant == 0
        @test (copy(quad4) + quad4).aff.terms[pt] == -84
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + quad4).terms[pair] == 8
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad3) - quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad4) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - quad2).aff.constant == 3
        @test (copy(quad1) - quad2).aff.terms[pt] == 1
        @test (copy(quad1) - quad2).aff.terms[finite] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, finite)
        @test (copy(quad1) - quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad1) - quad2).terms[pair] == 0
        @test (copy(quad4) - quad1).aff.constant == 0
        @test (copy(quad4) - quad1).aff.terms[pt] == -44
        @test (copy(quad4) - quad1).aff.terms[finite] == 2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - quad1).terms[pair] == 5
        # test measurefinite combos
        @test (copy(quad2) - quad4).aff.constant == -3
        @test (copy(quad2) - quad4).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, finite)
        @test (copy(quad2) - quad4).terms[pair] == 1
        # test finite combos
        @test (copy(quad3) - quad4).aff.constant == 0
        @test (copy(quad3) - quad4).aff.terms[pt] == 42
        pair = UnorderedPair{GeneralVariableRef}(finite, pt)
        @test (copy(quad3) - quad4).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - quad4, GenericQuadExpr{Float64, GeneralVariableRef})
        @test (copy(quad4) - quad4).aff.constant == 0
        @test (copy(quad4) - quad4).aff.terms[pt] == 0
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - quad4).terms[pair] == 0
    end
end
