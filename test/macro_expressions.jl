# Test operations between 2 different variables
@testset "Variable--Variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    # test addition
    @testset "Add" begin
        # test infinite + measure
        @test isa(@expression(m, inf + meas), GenericAffExpr{Float64, GeneralVariableRef})
        @test @expression(m, inf + meas).constant == 0.0
        @test @expression(m, inf + meas).terms[meas] == 1
        @test @expression(m, inf + meas).terms[inf] == 1
        # test point + hold
        @test isa(@expression(m, pt + hold), GenericAffExpr{Float64, FiniteVariableRef})
        @test @expression(m, pt + hold).terms[hold] == 1
        # test other combos
        @test isa(@expression(m, pt + meas), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par + pt), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par + inf), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par + meas), GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas + meas), GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, meas + meas).terms[meas] == 2
        @test isa(@expression(m, inf + inf), GenericAffExpr{Float64, InfiniteVariableRef})
        @test isa(@expression(m, pt + pt), GenericAffExpr{Float64, PointVariableRef})
    end
    # test subtraction
    @testset "Subtract" begin
        # test hold - measure
        @test isa(@expression(m, hold - meas), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, hold - meas).constant == 0.0
        @test @expression(m, hold - meas).terms[meas] == -1
        @test @expression(m, hold - meas).terms[hold] == 1
        # test point - hold
        @test isa(@expression(m, pt - hold), GenericAffExpr{Float64, FiniteVariableRef})
        @test @expression(m, pt - hold).terms[hold] == -1
        # test other combos
        @test isa(@expression(m, pt - meas), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par - pt), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par - inf), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par - meas), GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas - meas), GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, meas - meas).terms[meas] == 0
        @test isa(@expression(m, inf - inf), GenericAffExpr{Float64, InfiniteVariableRef})
        @test isa(@expression(m, pt - pt), GenericAffExpr{Float64, PointVariableRef})
    end
    # test multiplication
    @testset "Multiply" begin
        # test point * measure
        @test isa(@expression(m, pt * meas), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, pt * meas).aff.constant == 0.0
        pair = UnorderedPair{MeasureFiniteVariableRef}(pt, meas)
        @test @expression(m, pt * meas).terms[pair] == 1
        # test point * hold
        @test isa(@expression(m, pt * hold), GenericQuadExpr{Float64, FiniteVariableRef})
        pair = UnorderedPair{FiniteVariableRef}(pt, hold)
        @test @expression(m, pt * hold).terms[pair] == 1
        # test other combos
        @test isa(@expression(m, pt * meas), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par * pt), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par * inf), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par * meas), GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas * meas), GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test @expression(m, meas * meas).terms[pair] == 1
        @test isa(@expression(m, inf * inf), GenericQuadExpr{Float64, InfiniteVariableRef})
        @test isa(@expression(m, pt * pt), GenericQuadExpr{Float64, PointVariableRef})
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, inf / par)
        @test_throws ErrorException @expression(m, inf / inf)
        @test_throws ErrorException @expression(m, pt / hold)
    end
end

# Test operations between variables and GenericAffExpr
@testset "Variable--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    # test addition
    @testset "Add" begin
        # test infinite + aff1
        @test isa(@expression(m, inf + aff1), GenericAffExpr{Float64, GeneralVariableRef})
        @test @expression(m, inf + aff1).constant == -2
        @test @expression(m, inf + aff1).terms[pt] == 1
        @test @expression(m, inf + aff1).terms[inf] == 2
        # test meas + aff2
        @test isa(@expression(m, meas + aff2), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, meas + aff2).terms[hold] == 1
        # test other combos
        @test isa(@expression(m, hold + aff3), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par + aff2), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par + aff1), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par + aff3), GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas + (meas + meas)), GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, meas + (meas + meas)).terms[meas] == 3
    end
    # test subtraction
    @testset "Subtract" begin
        # test infinite - aff1
        @test isa(@expression(m, inf - aff1), GenericAffExpr{Float64, GeneralVariableRef})
        @test @expression(m, inf - aff1).constant == 2
        @test @expression(m, inf - aff1).terms[pt] == -1
        @test @expression(m, inf - aff1).terms[inf] == 0
        # test meas + aff2
        @test isa(@expression(m, meas - aff2), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, meas - aff2).terms[hold] == -1
        # test other combos
        @test isa(@expression(m, hold - aff3), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par - aff2), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par - aff1), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par - aff3), GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas - (meas + meas)), GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, meas - (meas + meas)).terms[meas] == -1
    end
    # test multiplication
    @testset "Multiply" begin
        # test infinite * aff1
        @test isa(@expression(m, inf * aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test @expression(m, inf * aff1).aff.constant == 0
        @test @expression(m, inf * aff1).aff.terms[inf] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test @expression(m, inf * aff1).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test @expression(m, inf * aff1).terms[pair] == 1
        # test meas * aff2
        @test isa(@expression(m, meas * aff2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, hold)
        @test @expression(m, meas * aff2).terms[pair] == 1
        @test length(@expression(m, meas * aff2).aff.terms) == 0
        # test other combos
        @test isa(@expression(m, hold * aff3), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, par * aff2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par * aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, par * aff3), GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, meas * (meas + meas)), GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test @expression(m, meas * (meas + meas)).terms[pair] == 2
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, inf / aff1)
        @test_throws ErrorException @expression(m, inf / aff2)
        @test_throws ErrorException @expression(m, pt / aff3)
    end
end

# Test operations between GenericAffExpr and variable
@testset "AffExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    # test addition
    @testset "Add" begin
        # test aff1 + inf
        @test isa(copy(aff1) + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test @expression(m, copy(aff1) + inf).constant == -2
        @test @expression(m, copy(aff1) + inf).terms[pt] == 1
        @test @expression(m, copy(aff1) + inf).terms[inf] == 2
        # test aff2 + meas
        @test isa(copy(aff2) + meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, copy(aff2) + meas).terms[hold] == 1
        # test other combos
        @test isa(aff3 + hold, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(aff2 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff1 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff3 + par, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas) + meas, GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, (meas + meas) + meas).terms[meas] == 3
    end
    # test subtraction
    @testset "Subtract" begin
        # test aff1 - inf
        @test isa(@expression(m, copy(aff1) - inf), GenericAffExpr{Float64, GeneralVariableRef})
        @test @expression(m, copy(aff1) - inf).constant == -2
        @test @expression(m, copy(aff1) - inf).terms[pt] == 1
        @test @expression(m, copy(aff1) - inf).terms[inf] == 0
        # test aff2 - meas
        @test isa(@expression(m, copy(aff2) - meas), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test @expression(m, copy(aff2) - meas).terms[hold] == 1
        # test other combos
        @test isa(@expression(m, aff3 - hold), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, aff2 - par), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff1 - par), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff3 - par), GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, (meas + meas) - meas), GenericAffExpr{Float64, MeasureRef})
        @test @expression(m, (meas + meas) - meas).terms[meas] == 1
    end
    # test multiplication
    @testset "Multiply" begin
        # test aff1 * inf
        @test isa(@expression(m, aff1 * inf), GenericQuadExpr{Float64, GeneralVariableRef})
        @test @expression(m, aff1 * inf).aff.constant == 0
        @test @expression(m, aff1 * inf).aff.terms[inf] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test @expression(m, aff1 * inf).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test @expression(m, aff1 * inf).terms[pair] == 1
        # test aff2 * meas
        @test isa(@expression(m, aff2 * meas), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, hold)
        @test @expression(m, aff2 * meas).terms[pair] == 1
        @test length((aff2 * meas).aff.terms) == 0
        # test other combos
        @test isa(@expression(m, aff3 * hold), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, aff2 * par), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff1 * par), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff3 * par), GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(@expression(m, (meas + meas * meas)), GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test @expression(m, (meas + meas) * meas).terms[pair] == 2
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, aff1 / inf)
        @test_throws ErrorException @expression(m, aff2 / inf)
        @test_throws ErrorException @expression(m, aff3 / pt)
    end
end

# Test operators on GenericAffExpr and GenericAffExpr
@testset "AffExpr--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff4 = (inf + inf) + 3
    aff5 = 4pt - 42
    aff_big = sum(test[i] for i = 1:length(test))
    # test addition
    @testset "Add" begin
        # test large input
        @test isa(@expression(m, aff_big + aff1), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(aff1) + aff_big), GenericAffExpr{Float64, GeneralVariableRef})
        # test various mixed combination types
        @test isa(@expression(m, copy(aff1) + aff2), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff2 + aff3), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(aff2) + aff5), GenericAffExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, aff4 + aff1), GenericAffExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test @expression(m, copy(aff1) + aff2).constant == -2
        @test @expression(m, copy(aff1) + aff2).terms[inf] == 1
        @test @expression(m, copy(aff1) + aff2).terms[pt] == 0
        @test @expression(m, aff4 + aff1).constant == 1
        @test @expression(m, aff4 + aff1).terms[inf] == 3
        @test @expression(m, aff4 + aff1).terms[pt] == 1
        # test results of measurefinite combinations
        @test @expression(m, aff2 + aff3).constant == 1
        @test @expression(m, aff2 + aff3).terms[meas] == -1
        @test @expression(m, aff2 + aff3).terms[hold] == 2
        @test @expression(m, copy(aff3) + aff5).constant == -41
        @test @expression(m, copy(aff3) + aff5).terms[meas] == -1
        @test @expression(m, copy(aff3) + aff5).terms[pt] == 4
        # test results of fininite combinations
        @test @expression(m, copy(aff2) + aff5).constant == -42
        @test @expression(m, copy(aff2) + aff5).terms[pt] == 3
        @test @expression(m, copy(aff2) + aff5).terms[hold] == 1
        # test pure additions
        @test isa(@expression(m, copy(aff4) + aff4), GenericAffExpr{Float64, InfiniteVariableRef})
        @test @expression(m, copy(aff4) + aff4).constant == 6
        @test @expression(m, copy(aff4) + aff4).terms[inf] == 4
        @test isa(@expression(m, copy(aff5) + aff5), GenericAffExpr{Float64, PointVariableRef})
        @test @expression(m, copy(aff5) + aff5).constant == -84
        @test @expression(m, copy(aff5) + aff5).terms[pt] == 8
    end
    # test subtraction
    @testset "Subtract" begin
        # test various mixed combination types
        @test isa(@expression(m, copy(aff1) - aff2), GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff2 - aff3), GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(aff2) - aff5), GenericAffExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, aff4 - aff1), GenericAffExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test @expression(m, copy(aff1) - aff2).constant == -2
        @test @expression(m, copy(aff1) - aff2).terms[inf] == 1
        @test @expression(m, copy(aff1) - aff2).terms[pt] == 2
        @test @expression(m, aff4 - aff1).constant == 5
        @test @expression(m, aff4 - aff1).terms[inf] == 1
        @test @expression(m, aff4 - aff1).terms[pt] == -1
        # test results of measurefinite combinations
        @test @expression(m, aff2 - aff3).constant == -1
        @test @expression(m, aff2 - aff3).terms[meas] == 1
        @test @expression(m, aff2 - aff3).terms[hold] == 0
        @test @expression(m, copy(aff3) - aff5).constant == 43
        @test @expression(m, copy(aff3) - aff5).terms[meas] == -1
        @test @expression(m, copy(aff3) - aff5).terms[pt] == -4
        # test results of fininite combinations
        @test @expression(m, copy(aff2) - aff5).constant == 42
        @test @expression(m, copy(aff2) - aff5).terms[pt] == -5
        @test @expression(m, copy(aff2) - aff5).terms[hold] == 1
        # test pure additions
        @test isa(@expression(m, copy(aff4) - aff4), GenericAffExpr{Float64, InfiniteVariableRef})
        @test @expression(m, copy(aff4) - aff4).constant == 0
        @test @expression(m, copy(aff4) - aff4).terms[inf] == 0
        @test isa(@expression(m, copy(aff5) - aff5), GenericAffExpr{Float64, PointVariableRef})
        @test @expression(m, copy(aff5) - aff5).constant == 0
        @test @expression(m, copy(aff5) - aff5).terms[pt] == 0
    end
    # test multiplication
    @testset "Multiply" begin
        # test various mixed combination types
        @test isa(@expression(m, aff1 * aff2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, aff2 * aff3), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, aff2 * aff5), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, aff4 * aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test @expression(m, aff1 * aff2).aff.constant == 0
        @test @expression(m, aff1 * aff2).aff.terms[hold] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, aff1 * aff2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, aff1 * aff2).terms[pair] == -1
        @test @expression(m, aff4 * aff1).aff.constant == -6
        @test @expression(m, aff4 * aff1).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, inf)
        @test @expression(m, aff4 * aff1).terms[pair] == 2
        pair = UnorderedPair{GeneralVariableRef}(inf, pt)
        @test @expression(m, aff4 * aff1).terms[pair] == 2
        # test results of measurefinite combinations
        @test @expression(m, aff2 * aff3).aff.constant == 0
        @test @expression(m, aff2 * aff3).aff.terms[pt] == -1
        pair = UnorderedPair{MeasureFiniteVariableRef}(pt, meas)
        @test @expression(m, aff2 * aff3).terms[pair] == 1
        pair = UnorderedPair{MeasureFiniteVariableRef}(hold, hold)
        @test @expression(m, aff2 * aff3).terms[pair] == 1
        @test @expression(m, aff3 * aff5).aff.constant == -42
        @test @expression(m, aff3 * aff5).aff.terms[meas] == 42
        pair = UnorderedPair{MeasureFiniteVariableRef}(hold, pt)
        @test @expression(m, aff3 * aff5).terms[pair] == 4
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, pt)
        @test @expression(m, aff3 * aff5).terms[pair] == -4
        # test results of fininite combinations
        @test @expression(m, aff2 * aff5).aff.constant == 0
        @test @expression(m, aff2 * aff5).aff.terms[hold] == -42
        pair = UnorderedPair{FiniteVariableRef}(hold, pt)
        @test @expression(m, aff2 * aff5).terms[pair] == 4
        pair = UnorderedPair{FiniteVariableRef}(pt, pt)
        @test @expression(m, aff2 * aff5).terms[pair] == -4
        # test pure additions
        @test isa(@expression(m, aff4 * aff4), GenericQuadExpr{Float64, InfiniteVariableRef})
        @test @expression(m, aff4 * aff4).aff.constant == 9
        @test @expression(m, aff4 * aff4).aff.terms[inf] == 12
        pair = UnorderedPair{InfiniteVariableRef}(inf, inf)
        @test @expression(m, aff4 * aff4).terms[pair] == 4
        @test isa(@expression(m, aff5 * aff5), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, aff5 * aff5).aff.constant == 42 * 42
        @test @expression(m, aff5 * aff5).aff.terms[pt] == 2 * 4 * -42
        pair = UnorderedPair{PointVariableRef}(pt, pt)
        @test @expression(m, aff5 * aff5).terms[pair] == 16
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, aff1 / aff2)
        @test_throws ErrorException @expression(m, aff2 / aff4)
        @test_throws ErrorException @expression(m, aff3 / aff5)
        @test_throws ErrorException @expression(m, aff5 / aff5)
    end
end

# Test operations for GenericQuadExpr--GeneralVariableRef
@testset "QuadExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*hold - inf*pt + pt*hold - pt^2 - 2hold + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*hold - pt*meas + pt - 3)
    quad3 = hold * aff2       #(hold^2 - hold*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Add" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) + inf), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) + meas), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) + pt), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad4) + par), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) + inf).aff.constant == 0
        @test @expression(m, copy(quad1) + inf).aff.terms[pt] == 2
        @test @expression(m, copy(quad1) + inf).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) + inf).terms[pair] == 1
        @test @expression(m, copy(quad4) + par).aff.constant == 0
        @test @expression(m, copy(quad4) + par).aff.terms[pt] == -42
        @test @expression(m, copy(quad4) + par).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + par).terms[pair] == 4
        # test measurefinite combos
        @test @expression(m, copy(quad2) + meas).aff.constant == -3
        @test @expression(m, copy(quad2) + meas).aff.terms[pt] == 1
        @test @expression(m, copy(quad2) + meas).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) + meas).terms[pair] == 1
        @test @expression(m, copy(quad2) + pt).aff.constant == -3
        @test @expression(m, copy(quad2) + pt).aff.terms[pt] == 2
        # test finite combos
        @test @expression(m, copy(quad3) + pt).aff.constant == 0
        @test @expression(m, copy(quad3) + pt).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) + pt).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) + pt), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) + pt).aff.constant == 0
        @test @expression(m, copy(quad4) + pt).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + pt).terms[pair] == 4
    end
    # test subtraction
    @testset "Subtract" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) - inf), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) - meas), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) - pt), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad4) - par), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) - inf).aff.constant == 0
        @test @expression(m, copy(quad1) - inf).aff.terms[pt] == 2
        @test @expression(m, copy(quad1) - inf).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) - inf).terms[pair] == 1
        @test @expression(m, copy(quad4) - par).aff.constant == 0
        @test @expression(m, copy(quad4) - par).aff.terms[pt] == -42
        @test @expression(m, copy(quad4) - par).aff.terms[par] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) - par).terms[pair] == 4
        # test measurefinite combos
        @test @expression(m, copy(quad2) - meas).aff.constant == -3
        @test @expression(m, copy(quad2) - meas).aff.terms[pt] == 1
        @test @expression(m, copy(quad2) - meas).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) - meas).terms[pair] == 1
        @test @expression(m, copy(quad2) - pt).aff.constant == -3
        @test @expression(m, copy(quad2) - pt).aff.terms[pt] == 0
        # test finite combos
        @test @expression(m, copy(quad3) - pt).aff.constant == 0
        @test @expression(m, copy(quad3) - pt).aff.terms[pt] == -1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) - pt).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) - pt), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) - pt).aff.constant == 0
        @test @expression(m, copy(quad4) - pt).aff.terms[pt] == -43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) - pt).terms[pair] == 4
    end
    # test multiplication
    @testset "Multiply" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 * inf)
        @test_throws ErrorException @expression(m, quad2 * meas)
        @test_throws ErrorException @expression(m, quad3 * par)
        @test_throws ErrorException @expression(m, quad4 * pt)
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 / inf)
        @test_throws ErrorException @expression(m, quad2 / meas)
        @test_throws ErrorException @expression(m, quad3 / par)
        @test_throws ErrorException @expression(m, quad4 / pt)
    end
end

# Test operations for GeneralVariableRef--GenericQuadExpr
@testset "Variable--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*hold - inf*pt + pt*hold - pt^2 - 2hold + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*hold - pt*meas + pt - 3)
    quad3 = hold * aff2       #(hold^2 - hold*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Add" begin
        # test mixed combo types
        @test isa(@expression(m, inf + copy(quad1)), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, meas + copy(quad2)), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, pt + copy(quad3)), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, par + copy(quad4)), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, inf + copy(quad1)).aff.constant == 0
        @test @expression(m, inf + copy(quad1)).aff.terms[pt] == 2
        @test @expression(m, inf + copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, inf + copy(quad1)).terms[pair] == 1
        @test @expression(m, par + copy(quad4)).aff.constant == 0
        @test @expression(m, par + copy(quad4)).aff.terms[pt] == -42
        @test @expression(m, par + copy(quad4)).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, par + copy(quad4)).terms[pair] == 4
        # test measurefinite combos
        @test @expression(m, meas + copy(quad2)).aff.constant == -3
        @test @expression(m, meas + copy(quad2)).aff.terms[pt] == 1
        @test @expression(m, meas + copy(quad2)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, meas + copy(quad2)).terms[pair] == 1
        @test @expression(m, pt + copy(quad2)).aff.constant == -3
        @test @expression(m, pt + copy(quad2)).aff.terms[pt] == 2
        # test finite combos
        @test @expression(m, pt + copy(quad3)).aff.constant == 0
        @test @expression(m, pt + copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, pt + copy(quad3)).terms[pair] == -1
        # test same types
        @test isa(@expression(m, pt + copy(quad4)), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, pt + copy(quad4)).aff.constant == 0
        @test @expression(m, pt + copy(quad4)).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, pt + copy(quad4)).terms[pair] == 4
    end
    # test subtraction
    @testset "Subtract" begin
        # test mixed combo types
        @test isa(@expression(m, inf - copy(quad1)), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, meas - copy(quad2)), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, pt - copy(quad3)), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, par - copy(quad4)), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, inf - copy(quad1)).aff.constant == 0
        @test @expression(m, inf - copy(quad1)).aff.terms[pt] == -2
        @test @expression(m, inf - copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, inf - copy(quad1)).terms[pair] == -1
        @test @expression(m, par - copy(quad4)).aff.constant == 0
        @test @expression(m, par - copy(quad4)).aff.terms[pt] == 42
        @test @expression(m, par - copy(quad4)).aff.terms[par] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, par - copy(quad4)).terms[pair] == -4
        # test measurefinite combos
        @test @expression(m, meas - copy(quad2)).aff.constant == 3
        @test @expression(m, meas - copy(quad2)).aff.terms[pt] == -1
        @test @expression(m, meas - copy(quad2)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, meas - copy(quad2)).terms[pair] == -1
        @test @expression(m, pt - copy(quad2)).aff.constant == 3
        @test @expression(m, pt - copy(quad2)).aff.terms[pt] == 0
        # test finite combos
        @test @expression(m, pt - copy(quad3)).aff.constant == 0
        @test @expression(m, pt - copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, pt - copy(quad3)).terms[pair] == 1
        # test same types
        @test isa(@expression(m, pt - copy(quad4)), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, pt - copy(quad4)).aff.constant == 0
        @test @expression(m, pt - copy(quad4)).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, pt - copy(quad4)).terms[pair] == -4
    end
    # test multiplication
    @testset "Multiply" begin
        # test some combos
        @test_throws ErrorException @expression(m, inf * quad1)
        @test_throws ErrorException @expression(m, meas * quad2)
        @test_throws ErrorException @expression(m, par * quad3)
        @test_throws ErrorException @expression(m, pt * quad4)
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, inf / quad1)
        @test_throws ErrorException @expression(m, meas / quad2)
        @test_throws ErrorException @expression(m, par / quad3)
        @test_throws ErrorException @expression(m, pt / quad4)
    end
end

# Test operators on GenericAffExpr--GenericQuadExpr
@testset "AffExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*hold - inf*pt + pt*hold - pt^2 - 2hold + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*hold - pt*meas + pt - 3)
    quad3 = hold * aff2       #(hold^2 - hold*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Add" begin
        # test mixed combo types
        @test isa(@expression(m, copy(aff1) + quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(aff2) + quad2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(aff5) + quad3), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(aff3) + quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(aff1) + quad1).aff.constant == -2
        @test @expression(m, copy(aff1) + quad1).aff.terms[pt] == 3
        @test @expression(m, copy(aff1) + quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(aff1) + quad1).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff1) + quad1).terms[pair] == -1
        @test @expression(m, copy(aff3) + quad1).aff.constant == 1
        @test @expression(m, copy(aff3) + quad1).aff.terms[meas] == -1
        @test @expression(m, copy(aff3) + quad1).aff.terms[hold] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(aff3) + quad1).terms[pair] == 1
        # test measurefinite combos
        @test @expression(m, copy(aff3) + quad2).aff.constant == -2
        @test @expression(m, copy(aff3) + quad2).aff.terms[pt] == 1
        @test @expression(m, copy(aff3) + quad2).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(aff3) + quad2).terms[pair] == 1
        @test @expression(m, copy(aff3) + quad4).aff.constant == 1
        @test @expression(m, copy(aff3) + quad4).aff.terms[pt] == -42
        @test @expression(m, copy(aff3) + quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff3) + quad4).terms[pair] == 4
        # test finite combos
        @test @expression(m, copy(aff2) + quad3).aff.constant == 0
        @test @expression(m, copy(aff2) + quad3).aff.terms[pt] == -1
        @test @expression(m, copy(aff2) + quad3).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(aff2) + quad3).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(aff5) + quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(aff5) + quad4).aff.constant == -42
        @test @expression(m, copy(aff5) + quad4).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff5) + quad4).terms[pair] == 4
    end
    # test subtraction
    @testset "Subtract" begin
        # test mixed combo types
        @test isa(@expression(m, copy(aff1) - quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(aff2) - quad2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(aff5) - quad3), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(aff3) - quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(aff1) - quad1).aff.constant == -2
        @test @expression(m, copy(aff1) - quad1).aff.terms[pt] == -1
        @test @expression(m, copy(aff1) - quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(aff1) - quad1).terms[pair] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff1) - quad1).terms[pair] == 1
        @test @expression(m, copy(aff3) - quad1).aff.constant == 1
        @test @expression(m, copy(aff3) - quad1).aff.terms[meas] == -1
        @test @expression(m, copy(aff3) - quad1).aff.terms[hold] == 3
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(aff3) - quad1).terms[pair] == -1
        # test measurefinite combos
        @test @expression(m, copy(aff3) - quad2).aff.constant == 4
        @test @expression(m, copy(aff3) - quad2).aff.terms[pt] == -1
        @test @expression(m, copy(aff3) - quad2).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(aff3) - quad2).terms[pair] == -1
        @test @expression(m, copy(aff3) - quad4).aff.constant == 1
        @test @expression(m, copy(aff3) - quad4).aff.terms[pt] == 42
        @test @expression(m, copy(aff3) - quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff3) - quad4).terms[pair] == -4
        # test finite combos
        @test @expression(m, copy(aff2) - quad3).aff.constant == 0
        @test @expression(m, copy(aff2) - quad3).aff.terms[pt] == -1
        @test @expression(m, copy(aff2) - quad3).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(aff2) - quad3).terms[pair] == 1
        # test same types
        @test isa(@expression(m, copy(aff5) - quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(aff5) - quad4).aff.constant == -42
        @test @expression(m, copy(aff5) - quad4).aff.terms[pt] == 46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(aff5) - quad4).terms[pair] == -4
    end
    # test multiplication
    @testset "Multiply" begin
        # test some combos
        @test_throws ErrorException @expression(m, aff1 * quad1)
        @test_throws ErrorException @expression(m, aff2 * quad3)
        @test_throws ErrorException @expression(m, aff3 * quad4)
        @test_throws ErrorException @expression(m, aff5 * quad2)
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, aff1 / quad1)
        @test_throws ErrorException @expression(m, aff2 / quad3)
        @test_throws ErrorException @expression(m, aff3 / quad4)
        @test_throws ErrorException @expression(m, aff5 / quad2)
    end
end

# Test operators on GenericQuadExpr--GenericAffExpr
@testset "QuadExpr--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*hold - inf*pt + pt*hold - pt^2 - 2hold + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*hold - pt*meas + pt - 3)
    quad3 = hold * aff2       #(hold^2 - hold*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Add" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) + copy(aff1)), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) + copy(aff2) + quad2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) + copy(aff5)), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad1) + copy(aff3)), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) + copy(aff1)).aff.constant == -2
        @test @expression(m, copy(quad1) + copy(aff1)).aff.terms[pt] == 3
        @test @expression(m, copy(quad1) + copy(aff1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) + copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad1) + copy(aff3)).terms[pair] == -1
        @test @expression(m, copy(quad1) + copy(aff3)).aff.constant == 1
        @test @expression(m, copy(quad1) + copy(aff3)).aff.terms[meas] == -1
        @test @expression(m, copy(quad1) + copy(aff3)).aff.terms[hold] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) + copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test @expression(m, copy(quad2) + copy(aff3)).aff.constant == -2
        @test @expression(m, copy(quad2) + copy(aff3)).aff.terms[pt] == 1
        @test @expression(m, copy(quad2) + copy(aff3)).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) + copy(aff3)).terms[pair] == 1
        @test @expression(m, copy(quad4) + copy(aff3)).aff.constant == 1
        @test @expression(m, copy(quad4) + copy(aff3)).aff.terms[pt] == -42
        @test @expression(m, copy(quad4) + copy(aff3)).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + copy(aff3)).terms[pair] == 4
        # test finite combos
        @test @expression(m, copy(quad3) + copy(aff2)).aff.constant == 0
        @test @expression(m, copy(quad3) + copy(aff2)).aff.terms[pt] == -1
        @test @expression(m, copy(quad3) + copy(aff2)).aff.terms[hold] == 1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) + copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) + copy(aff5)), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) + copy(aff5)).aff.constant == -42
        @test @expression(m, copy(quad4) + copy(aff5)).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + copy(aff5)).terms[pair] == 4
    end
    # test subtraction
    @testset "Subtract" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) - copy(aff1)), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) - copy(aff2)), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) - copy(aff5)), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad1) - copy(aff3)), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) - copy(aff1)).aff.constant == 2
        @test @expression(m, copy(quad1) - copy(aff1)).aff.terms[pt] == 1
        @test @expression(m, copy(quad1) - copy(aff1)).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) - copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad1) - copy(aff1)).terms[pair] == -1
        @test @expression(m, copy(quad1) - copy(aff3)).aff.constant == -1
        @test @expression(m, copy(quad1) - copy(aff3)).aff.terms[meas] == 1
        @test @expression(m, copy(quad1) - copy(aff3)).aff.terms[hold] == -3
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) - copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test @expression(m, copy(quad2) - copy(aff3)).aff.constant == -4
        @test @expression(m, copy(quad2) - copy(aff3)).aff.terms[pt] == 1
        @test @expression(m, copy(quad2) - copy(aff3)).aff.terms[hold] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) - copy(aff3)).terms[pair] == 1
        @test @expression(m, copy(quad4) - copy(aff3)).aff.constant == -1
        @test @expression(m, copy(quad4) - copy(aff3)).aff.terms[pt] == -42
        @test @expression(m, copy(quad4) - copy(aff3)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4)- copy(aff3)).terms[pair] == 4
        # test finite combos
        @test @expression(m, copy(quad3) - copy(aff2)).aff.constant == 0
        @test @expression(m, copy(quad3) - copy(aff2)).aff.terms[pt] == 1
        @test @expression(m, copy(quad3) - copy(aff2)).aff.terms[hold] == -1
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) - copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) - copy(aff5)), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) - copy(aff5)).aff.constant == 42
        @test @expression(m, copy(quad4) - copy(aff5)).aff.terms[pt] == -46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) - copy(aff5)).terms[pair] == 4
    end
    # test multiplication
    @testset "Multiply" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 * aff1)
        @test_throws ErrorException @expression(m, quad3 * aff2)
        @test_throws ErrorException @expression(m, quad4 * aff3)
        @test_throws ErrorException @expression(m, quad2 * aff5)
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 / aff1)
        @test_throws ErrorException @expression(m, quad3 / aff2)
        @test_throws ErrorException @expression(m, quad4 / aff3)
        @test_throws ErrorException @expression(m, quad2 / aff5)
    end
end

# Test operators on GenericQuadExpr--GenericQuadExpr
@testset "QuadExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    @hold_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = hold - pt
    aff3 = (hold - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*hold - inf*pt + pt*hold - pt^2 - 2hold + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*hold - pt*meas + pt - 3)
    quad3 = hold * aff2       #(hold^2 - hold*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Add" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) + quad2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) + quad4), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) + quad4), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad4) + quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) + quad2).aff.constant == -3
        @test @expression(m, copy(quad1) + quad2).aff.terms[pt] == 3
        @test @expression(m, copy(quad1) + quad2).aff.terms[hold] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) + quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad1) + quad2).terms[pair] == 2
        @test @expression(m, copy(quad4) + quad1).aff.constant == 0
        @test @expression(m, copy(quad4) + quad1).aff.terms[pt] == -40
        @test @expression(m, copy(quad4) + quad1).aff.terms[hold] == -2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + quad1).terms[pair] == 3
        # test measurefinite combos
        @test @expression(m, copy(quad2) + quad4).aff.constant == -3
        @test @expression(m, copy(quad2) + quad4).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) + quad4).terms[pair] == 1
        # test finite combos
        @test @expression(m, copy(quad3) + quad4).aff.constant == 0
        @test @expression(m, copy(quad3) + quad4).aff.terms[pt] == -42
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) + quad4).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) + quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) + quad4).aff.constant == 0
        @test @expression(m, copy(quad4) + quad4).aff.terms[pt] == -84
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) + quad4).terms[pair] == 8
    end
    # test subtraction
    @testset "Subtract" begin
        # test mixed combo types
        @test isa(@expression(m, copy(quad1) - quad2), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(@expression(m, copy(quad2) - quad4), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(@expression(m, copy(quad3) - quad4), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(@expression(m, copy(quad4) - quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test @expression(m, copy(quad1) - quad2).aff.constant == 3
        @test @expression(m, copy(quad1) - quad2).aff.terms[pt] == 1
        @test @expression(m, copy(quad1) - quad2).aff.terms[hold] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, hold)
        @test @expression(m, copy(quad1) - quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad1) - quad2).terms[pair] == 0
        @test @expression(m, copy(quad4) - quad1).aff.constant == 0
        @test @expression(m, copy(quad4) - quad1).aff.terms[pt] == -44
        @test @expression(m, copy(quad4) - quad1).aff.terms[hold] == 2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) - quad1).terms[pair] == 5
        # test measurefinite combos
        @test @expression(m, copy(quad2) - quad4).aff.constant == -3
        @test @expression(m, copy(quad2) - quad4).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, hold)
        @test @expression(m, copy(quad2) - quad4).terms[pair] == 1
        # test finite combos
        @test @expression(m, copy(quad3) - quad4).aff.constant == 0
        @test @expression(m, copy(quad3) - quad4).aff.terms[pt] == 42
        pair = UnorderedPair{GeneralVariableRef}(hold, pt)
        @test @expression(m, copy(quad3) - quad4).terms[pair] == -1
        # test same types
        @test isa(@expression(m, copy(quad4) - quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test @expression(m, copy(quad4) - quad4).aff.constant == 0
        @test @expression(m, copy(quad4) - quad4).aff.terms[pt] == 0
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test @expression(m, copy(quad4) - quad4).terms[pair] == 0
    end
    # test multiplication
    @testset "Multiply" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 * quad2)
        @test_throws ErrorException @expression(m, quad3 * quad1)
        @test_throws ErrorException @expression(m, quad4 * quad4)
        @test_throws ErrorException @expression(m, quad2 * quad3)
    end
    # test division
    @testset "Divide" begin
        # test some combos
        @test_throws ErrorException @expression(m, quad1 / quad2)
        @test_throws ErrorException @expression(m, quad3 / quad1)
        @test_throws ErrorException @expression(m, quad4 / quad4)
        @test_throws ErrorException @expression(m, quad2 / quad3)
    end
end
