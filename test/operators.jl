# Test operations between 2 different variables
@testset "Variable--Variable" begin
    # initialize model and references
    m = InfiniteModel()
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    # test addition
    @testset "Base.:+" begin
        # test infinite + measure
        @test isa(inf + meas, GenericAffExpr{Float64, GeneralVariableRef})
        @test (inf + meas).constant == 0.0
        @test (inf + meas).terms[meas] == 1
        @test (inf + meas).terms[inf] == 1
        # test point + global
        @test isa(pt + glob, GenericAffExpr{Float64, FiniteVariableRef})
        @test (pt + glob).terms[glob] == 1
        # test other combos
        @test isa(pt + meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par + pt, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + meas, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas + meas, GenericAffExpr{Float64, MeasureRef})
        @test (meas + meas).terms[meas] == 2
        @test isa(inf + inf, GenericAffExpr{Float64, InfiniteVariableRef})
        @test isa(pt + pt, GenericAffExpr{Float64, PointVariableRef})
    end
    # test subtraction
    @testset "Base.:-" begin
        # test global - measure
        @test isa(glob - meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test (glob - meas).constant == 0.0
        @test (glob - meas).terms[meas] == -1
        @test (glob - meas).terms[glob] == 1
        # test point - global
        @test isa(pt - glob, GenericAffExpr{Float64, FiniteVariableRef})
        @test (pt - glob).terms[glob] == -1
        # test other combos
        @test isa(pt - meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par - pt, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - meas, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas - meas, GenericAffExpr{Float64, MeasureRef})
        @test (meas - meas).terms[meas] == 0
        @test isa(inf - inf, GenericAffExpr{Float64, InfiniteVariableRef})
        @test isa(pt - pt, GenericAffExpr{Float64, PointVariableRef})
    end
    # test multiplication
    @testset "Base.:*" begin
        # test point * measure
        @test isa(pt * meas, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test (pt * meas).aff.constant == 0.0
        pair = UnorderedPair{MeasureFiniteVariableRef}(pt, meas)
        @test (pt * meas).terms[pair] == 1
        # test point * global
        @test isa(pt * glob, GenericQuadExpr{Float64, FiniteVariableRef})
        pair = UnorderedPair{FiniteVariableRef}(pt, glob)
        @test (pt * glob).terms[pair] == 1
        # test other combos
        @test isa(pt * meas, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par * pt, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * meas, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas * meas, GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test (meas * meas).terms[pair] == 1
        @test isa(inf * inf, GenericQuadExpr{Float64, InfiniteVariableRef})
        @test isa(pt * pt, GenericQuadExpr{Float64, PointVariableRef})
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException inf / par
        @test_throws ErrorException inf / inf
        @test_throws ErrorException pt / glob
    end
end

# Test operations between variables and GenericAffExpr
@testset "Variable--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    # test addition
    @testset "Base.:+" begin
        # test infinite + aff1
        @test isa(inf + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test (inf + aff1).constant == -2
        @test (inf + aff1).terms[pt] == 1
        @test (inf + aff1).terms[inf] == 2
        # test meas + aff2
        @test isa(meas + aff2, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test (meas + aff2).terms[glob] == 1
        # test other combos
        @test isa(glob + aff3, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par + aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par + aff3, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas + (meas + meas), GenericAffExpr{Float64, MeasureRef})
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
        @test isa(meas - aff2, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test (meas - aff2).terms[glob] == -1
        # test other combos
        @test isa(glob - aff3, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par - aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - aff1, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(par - aff3, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas - (meas + meas), GenericAffExpr{Float64, MeasureRef})
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
        @test isa(meas * aff2, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, glob)
        @test (meas * aff2).terms[pair] == 1
        @test length((meas * aff2).aff.terms) == 0
        # test other combos
        @test isa(glob * aff3, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(par * aff2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * aff1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(par * aff3, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa(meas * (meas + meas), GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test (meas * (meas + meas)).terms[pair] == 2
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException inf / aff1
        @test_throws ErrorException inf / aff2
        @test_throws ErrorException pt / aff3
    end
end

# Test operations between GenericAffExpr and variable
@testset "AffExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    # test addition
    @testset "Base.:+" begin
        # test aff1 + inf
        @test isa(aff1 + inf, GenericAffExpr{Float64, GeneralVariableRef})
        result = ((inf + pt) - 2) + inf
        @test result.constant == -2
        @test result.terms[pt] == 1
        @test result.terms[inf] == 2
        # test aff2 + meas
        @test isa(aff2 + meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        result = (glob - pt) + meas
        @test result.terms[glob] == 1
        # test other combos
        @test isa(aff3 + glob, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(aff2 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff1 + par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff3 + par, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas) + meas, GenericAffExpr{Float64, MeasureRef})
        @test ((meas + meas) + meas).terms[meas] == 3
    end
    # test subtraction
    @testset "Base.:-" begin
        # test aff1 - inf
        @test isa(aff1 - inf, GenericAffExpr{Float64, GeneralVariableRef})
        result = ((inf + pt) - 2) - inf
        @test result.constant == -2
        @test result.terms[pt] == 1
        @test result.terms[inf] == 0
        # test aff2 - meas
        @test isa(aff2 - meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        result = (glob - pt) - meas
        @test result.terms[glob] == 1
        # test other combos
        @test isa(aff3 - glob, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(aff2 - par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff1 - par, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff3 - par, GenericAffExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas) - meas, GenericAffExpr{Float64, MeasureRef})
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
        @test isa(aff2 * meas, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, glob)
        @test (aff2 * meas).terms[pair] == 1
        @test length((aff2 * meas).aff.terms) == 0
        # test other combos
        @test isa(aff3 * glob, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(aff2 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff1 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff3 * par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test same
        @test isa((meas + meas * meas), GenericQuadExpr{Float64, MeasureRef})
        pair = UnorderedPair{MeasureRef}(meas, meas)
        @test ((meas + meas) * meas).terms[pair] == 2
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException aff1 / inf
        @test_throws ErrorException aff2 / inf
        @test_throws ErrorException aff3 / pt
    end
end

# Test operators on GenericAffExpr and GenericAffExpr
@testset "AffExpr--AffExpr" begin
    
end
