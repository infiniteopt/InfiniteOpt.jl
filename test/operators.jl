# Test operator methods
@testset "Internal Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    # test _var_type_parser
    @testset "_var_type_parser" begin
        # test same types
        @test InfiniteOpt._var_type_parser(MeasureRef, MeasureRef) == MeasureRef
        @test InfiniteOpt._var_type_parser(ParameterRef,
                                           ParameterRef) == ParameterRef
        @test InfiniteOpt._var_type_parser(InfiniteVariableRef,
                                     InfiniteVariableRef) == InfiniteVariableRef
        @test InfiniteOpt._var_type_parser(PointVariableRef,
                                           PointVariableRef) == PointVariableRef
        @test InfiniteOpt._var_type_parser(GlobalVariableRef,
                                         GlobalVariableRef) == GlobalVariableRef
        # test finite types
        @test InfiniteOpt._var_type_parser(PointVariableRef,
                                         GlobalVariableRef) == FiniteVariableRef
        @test InfiniteOpt._var_type_parser(GlobalVariableRef,
                                         FiniteVariableRef) == FiniteVariableRef
        # test measure finite types
        @test InfiniteOpt._var_type_parser(PointVariableRef,
                                         MeasureRef) == MeasureFiniteVariableRef
        @test InfiniteOpt._var_type_parser(GlobalVariableRef,
                                         MeasureRef) == MeasureFiniteVariableRef
        # test other combos
        @test InfiniteOpt._var_type_parser(InfiniteVariableRef,
                                               MeasureRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(ParameterRef,
                                           MeasureRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(InfiniteVariableRef,
                                           ParameterRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(InfiniteVariableRef,
                                         PointVariableRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(ParameterRef,
                                         PointVariableRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(InfiniteVariableRef,
                                        GlobalVariableRef) == GeneralVariableRef
        @test InfiniteOpt._var_type_parser(ParameterRef,
                                        GlobalVariableRef) == GeneralVariableRef
    end
    # test extension to add_to_expression!
    @testset "JuMP.add_to_expression!" begin
        # make expression
        quad = pt^2 + 3glob + 2
        # test with same types
        @test isa(add_to_expression!(copy(quad), 1., glob, pt),
                  GenericQuadExpr{Float64, FiniteVariableRef})
        pair = UnorderedPair{FiniteVariableRef}(glob, pt)
        @test add_to_expression!(copy(quad), 1., glob, pt).terms[pair] == 1
        # test with different types
        @test isa(add_to_expression!(copy(quad), 2., glob, inf),
                  GenericQuadExpr{Float64, GeneralVariableRef})
        pair = UnorderedPair{GeneralVariableRef}(glob, inf)
        @test add_to_expression!(copy(quad), 2., glob, inf).terms[pair] == 2
    end
end

# Test operations between 2 different variables
@testset "Variable--Variable" begin
    # initialize model and references
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
        @test isa(copy(aff1) + inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff1) + inf).constant == -2
        @test (copy(aff1) + inf).terms[pt] == 1
        @test (copy(aff1) + inf).terms[inf] == 2
        # test aff2 + meas
        @test isa(copy(aff2) + meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test (copy(aff2) + meas).terms[glob] == 1
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
        @test isa(copy(aff1) - inf, GenericAffExpr{Float64, GeneralVariableRef})
        @test (copy(aff1) - inf).constant == -2
        @test (copy(aff1) - inf).terms[pt] == 1
        @test (copy(aff1) - inf).terms[inf] == 0
        # test aff2 - meas
        @test isa(copy(aff2) - meas, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test (copy(aff2) - meas).terms[glob] == 1
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
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
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
        @test isa(aff2 + aff3, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(aff2) + aff5, GenericAffExpr{Float64, FiniteVariableRef})
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
        @test (aff2 + aff3).terms[glob] == 2
        @test (copy(aff3) + aff5).constant == -41
        @test (copy(aff3) + aff5).terms[meas] == -1
        @test (copy(aff3) + aff5).terms[pt] == 4
        # test results of fininite combinations
        @test (copy(aff2) + aff5).constant == -42
        @test (copy(aff2) + aff5).terms[pt] == 3
        @test (copy(aff2) + aff5).terms[glob] == 1
        # test pure additions
        @test isa(copy(aff4) + aff4, GenericAffExpr{Float64, InfiniteVariableRef})
        @test (copy(aff4) + aff4).constant == 6
        @test (copy(aff4) + aff4).terms[inf] == 4
        @test isa(copy(aff5) + aff5, GenericAffExpr{Float64, PointVariableRef})
        @test (copy(aff5) + aff5).constant == -84
        @test (copy(aff5) + aff5).terms[pt] == 8
    end
    # test subtraction
    @testset "Base.:-" begin
        # test various mixed combination types
        @test isa(copy(aff1) - aff2, GenericAffExpr{Float64, GeneralVariableRef})
        @test isa(aff2 - aff3, GenericAffExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(aff2) - aff5, GenericAffExpr{Float64, FiniteVariableRef})
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
        @test (aff2 - aff3).terms[glob] == 0
        @test (copy(aff3) - aff5).constant == 43
        @test (copy(aff3) - aff5).terms[meas] == -1
        @test (copy(aff3) - aff5).terms[pt] == -4
        # test results of fininite combinations
        @test (copy(aff2) - aff5).constant == 42
        @test (copy(aff2) - aff5).terms[pt] == -5
        @test (copy(aff2) - aff5).terms[glob] == 1
        # test pure additions
        @test isa(copy(aff4) - aff4, GenericAffExpr{Float64, InfiniteVariableRef})
        @test (copy(aff4) - aff4).constant == 0
        @test (copy(aff4) - aff4).terms[inf] == 0
        @test isa(copy(aff5) - aff5, GenericAffExpr{Float64, PointVariableRef})
        @test (copy(aff5) - aff5).constant == 0
        @test (copy(aff5) - aff5).terms[pt] == 0
    end
    # test multiplication
    @testset "Base.:*" begin
        # test various mixed combination types
        @test isa(aff1 * aff2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(aff2 * aff3, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(aff2 * aff5, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(aff4 * aff1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test results of general combinations
        @test (aff1 * aff2).aff.constant == 0
        @test (aff1 * aff2).aff.terms[glob] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
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
        pair = UnorderedPair{MeasureFiniteVariableRef}(pt, meas)
        @test (aff2 * aff3).terms[pair] == 1
        pair = UnorderedPair{MeasureFiniteVariableRef}(glob, glob)
        @test (aff2 * aff3).terms[pair] == 1
        @test (aff3 * aff5).aff.constant == -42
        @test (aff3 * aff5).aff.terms[meas] == 42
        pair = UnorderedPair{MeasureFiniteVariableRef}(glob, pt)
        @test (aff3 * aff5).terms[pair] == 4
        pair = UnorderedPair{MeasureFiniteVariableRef}(meas, pt)
        @test (aff3 * aff5).terms[pair] == -4
        # test results of fininite combinations
        @test (aff2 * aff5).aff.constant == 0
        @test (aff2 * aff5).aff.terms[glob] == -42
        pair = UnorderedPair{FiniteVariableRef}(glob, pt)
        @test (aff2 * aff5).terms[pair] == 4
        pair = UnorderedPair{FiniteVariableRef}(pt, pt)
        @test (aff2 * aff5).terms[pair] == -4
        # test pure additions
        @test isa(aff4 * aff4, GenericQuadExpr{Float64, InfiniteVariableRef})
        @test (aff4 * aff4).aff.constant == 9
        @test (aff4 * aff4).aff.terms[inf] == 12
        pair = UnorderedPair{InfiniteVariableRef}(inf, inf)
        @test (aff4 * aff4).terms[pair] == 4
        @test isa(aff5 * aff5, GenericQuadExpr{Float64, PointVariableRef})
        @test (aff5 * aff5).aff.constant == 42 * 42
        @test (aff5 * aff5).aff.terms[pt] == 2 * 4 * -42
        pair = UnorderedPair{PointVariableRef}(pt, pt)
        @test (aff5 * aff5).terms[pair] == 16
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException aff1 / aff2
        @test_throws ErrorException aff2 / aff4
        @test_throws ErrorException aff3 / aff5
        @test_throws ErrorException aff5 / aff5
    end
end

# Test operations for GenericQuadExpr--GeneralVariableRef
@testset "QuadExpr--Variable" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*glob - inf*pt + pt*glob - pt^2 - 2glob + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*glob - pt*meas + pt - 3)
    quad3 = glob * aff2       #(glob^2 - glob*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + meas, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) + pt, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad4) + par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) + inf).aff.constant == 0
        @test (copy(quad1) + inf).aff.terms[pt] == 2
        @test (copy(quad1) + inf).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
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
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) + meas).terms[pair] == 1
        @test (copy(quad2) + pt).aff.constant == -3
        @test (copy(quad2) + pt).aff.terms[pt] == 2
        # test finite combos
        @test (copy(quad3) + pt).aff.constant == 0
        @test (copy(quad3) + pt).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) + pt).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + pt, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) + pt).aff.constant == 0
        @test (copy(quad4) + pt).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + pt).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - inf, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - meas, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) - pt, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad4) - par, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - inf).aff.constant == 0
        @test (copy(quad1) - inf).aff.terms[pt] == 2
        @test (copy(quad1) - inf).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
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
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) - meas).terms[pair] == 1
        @test (copy(quad2) - pt).aff.constant == -3
        @test (copy(quad2) - pt).aff.terms[pt] == 0
        # test finite combos
        @test (copy(quad3) - pt).aff.constant == 0
        @test (copy(quad3) - pt).aff.terms[pt] == -1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) - pt).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - pt, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) - pt).aff.constant == 0
        @test (copy(quad4) - pt).aff.terms[pt] == -43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - pt).terms[pair] == 4
    end
    # test multiplication
    @testset "Base.:*" begin
        # test some combos
        @test_throws ErrorException quad1 * inf
        @test_throws ErrorException quad2 * meas
        @test_throws ErrorException quad3 * par
        @test_throws ErrorException quad4 * pt
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException quad1 / inf
        @test_throws ErrorException quad2 / meas
        @test_throws ErrorException quad3 / par
        @test_throws ErrorException quad4 / pt
    end
end

# Test operations for GeneralVariableRef--GenericQuadExpr
@testset "Variable--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*glob - inf*pt + pt*glob - pt^2 - 2glob + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*glob - pt*meas + pt - 3)
    quad3 = glob * aff2       #(glob^2 - glob*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(inf + copy(quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(meas + copy(quad2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(pt + copy(quad3), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(par + copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (inf + copy(quad1)).aff.constant == 0
        @test (inf + copy(quad1)).aff.terms[pt] == 2
        @test (inf + copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
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
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (meas + copy(quad2)).terms[pair] == 1
        @test (pt + copy(quad2)).aff.constant == -3
        @test (pt + copy(quad2)).aff.terms[pt] == 2
        # test finite combos
        @test (pt + copy(quad3)).aff.constant == 0
        @test (pt + copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (pt + copy(quad3)).terms[pair] == -1
        # test same types
        @test isa(pt + copy(quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test (pt + copy(quad4)).aff.constant == 0
        @test (pt + copy(quad4)).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (pt + copy(quad4)).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(inf - copy(quad1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(meas - copy(quad2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(pt - copy(quad3), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(par - copy(quad4), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (inf - copy(quad1)).aff.constant == 0
        @test (inf - copy(quad1)).aff.terms[pt] == -2
        @test (inf - copy(quad1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
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
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (meas - copy(quad2)).terms[pair] == -1
        @test (pt - copy(quad2)).aff.constant == 3
        @test (pt - copy(quad2)).aff.terms[pt] == 0
        # test finite combos
        @test (pt - copy(quad3)).aff.constant == 0
        @test (pt - copy(quad3)).aff.terms[pt] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (pt - copy(quad3)).terms[pair] == 1
        # test same types
        @test isa(pt - copy(quad4), GenericQuadExpr{Float64, PointVariableRef})
        @test (pt - copy(quad4)).aff.constant == 0
        @test (pt - copy(quad4)).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (pt - copy(quad4)).terms[pair] == -4
    end
    # test multiplication
    @testset "Base.:*" begin
        # test some combos
        @test_throws ErrorException inf * quad1
        @test_throws ErrorException meas * quad2
        @test_throws ErrorException par * quad3
        @test_throws ErrorException pt * quad4
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException inf / quad1
        @test_throws ErrorException meas / quad2
        @test_throws ErrorException par / quad3
        @test_throws ErrorException pt / quad4
    end
end

# Test operators on GenericAffExpr--GenericQuadExpr
@testset "AffExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*glob - inf*pt + pt*glob - pt^2 - 2glob + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*glob - pt*meas + pt - 3)
    quad3 = glob * aff2       #(glob^2 - glob*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(aff1) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) + quad2, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(aff5) + quad3, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(aff3) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(aff1) + quad1).aff.constant == -2
        @test (copy(aff1) + quad1).aff.terms[pt] == 3
        @test (copy(aff1) + quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(aff1) + quad1).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff1) + quad1).terms[pair] == -1
        @test (copy(aff3) + quad1).aff.constant == 1
        @test (copy(aff3) + quad1).aff.terms[meas] == -1
        @test (copy(aff3) + quad1).aff.terms[glob] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(aff3) + quad1).terms[pair] == 1
        # test measurefinite combos
        @test (copy(aff3) + quad2).aff.constant == -2
        @test (copy(aff3) + quad2).aff.terms[pt] == 1
        @test (copy(aff3) + quad2).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(aff3) + quad2).terms[pair] == 1
        @test (copy(aff3) + quad4).aff.constant == 1
        @test (copy(aff3) + quad4).aff.terms[pt] == -42
        @test (copy(aff3) + quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff3) + quad4).terms[pair] == 4
        # test finite combos
        @test (copy(aff2) + quad3).aff.constant == 0
        @test (copy(aff2) + quad3).aff.terms[pt] == -1
        @test (copy(aff2) + quad3).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(aff2) + quad3).terms[pair] == -1
        # test same types
        @test isa(copy(aff5) + quad4, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(aff5) + quad4).aff.constant == -42
        @test (copy(aff5) + quad4).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff5) + quad4).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(aff1) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(aff2) - quad2, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(aff5) - quad3, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(aff3) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(aff1) - quad1).aff.constant == -2
        @test (copy(aff1) - quad1).aff.terms[pt] == -1
        @test (copy(aff1) - quad1).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(aff1) - quad1).terms[pair] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff1) - quad1).terms[pair] == 1
        @test (copy(aff3) - quad1).aff.constant == 1
        @test (copy(aff3) - quad1).aff.terms[meas] == -1
        @test (copy(aff3) - quad1).aff.terms[glob] == 3
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(aff3) - quad1).terms[pair] == -1
        # test measurefinite combos
        @test (copy(aff3) - quad2).aff.constant == 4
        @test (copy(aff3) - quad2).aff.terms[pt] == -1
        @test (copy(aff3) - quad2).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(aff3) - quad2).terms[pair] == -1
        @test (copy(aff3) - quad4).aff.constant == 1
        @test (copy(aff3) - quad4).aff.terms[pt] == 42
        @test (copy(aff3) - quad4).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff3) - quad4).terms[pair] == -4
        # test finite combos
        @test (copy(aff2) - quad3).aff.constant == 0
        @test (copy(aff2) - quad3).aff.terms[pt] == -1
        @test (copy(aff2) - quad3).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(aff2) - quad3).terms[pair] == 1
        # test same types
        @test isa(copy(aff5) - quad4, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(aff5) - quad4).aff.constant == -42
        @test (copy(aff5) - quad4).aff.terms[pt] == 46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(aff5) - quad4).terms[pair] == -4
    end
    # test multiplication
    @testset "Base.:*" begin
        # test some combos
        @test_throws ErrorException aff1 * quad1
        @test_throws ErrorException aff2 * quad3
        @test_throws ErrorException aff3 * quad4
        @test_throws ErrorException aff5 * quad2
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException aff1 / quad1
        @test_throws ErrorException aff2 / quad3
        @test_throws ErrorException aff3 / quad4
        @test_throws ErrorException aff5 / quad2
    end
end

# Test operators on GenericQuadExpr--GenericAffExpr
@testset "QuadExpr--AffExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*glob - inf*pt + pt*glob - pt^2 - 2glob + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*glob - pt*meas + pt - 3)
    quad3 = glob * aff2       #(glob^2 - glob*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + copy(aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + copy(aff2) + quad2, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) + copy(aff5), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad1) + copy(aff3), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1)+ copy(aff1)).aff.constant == -2
        @test (copy(quad1) + copy(aff1)).aff.terms[pt] == 3
        @test (copy(quad1) + copy(aff1)).aff.terms[inf] == 1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) + copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad1) + copy(aff3)).terms[pair] == -1
        @test (copy(quad1) + copy(aff3)).aff.constant == 1
        @test (copy(quad1) + copy(aff3)).aff.terms[meas] == -1
        @test (copy(quad1) + copy(aff3)).aff.terms[glob] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) + copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test (copy(quad2) + copy(aff3)).aff.constant == -2
        @test (copy(quad2) + copy(aff3)).aff.terms[pt] == 1
        @test (copy(quad2) + copy(aff3)).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) + copy(aff3)).terms[pair] == 1
        @test (copy(quad4) + copy(aff3)).aff.constant == 1
        @test (copy(quad4) + copy(aff3)).aff.terms[pt] == -42
        @test (copy(quad4) + copy(aff3)).aff.terms[meas] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + copy(aff3)).terms[pair] == 4
        # test finite combos
        @test (copy(quad3) + copy(aff2)).aff.constant == 0
        @test (copy(quad3) + copy(aff2)).aff.terms[pt] == -1
        @test (copy(quad3) + copy(aff2)).aff.terms[glob] == 1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) + copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + copy(aff5), GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) + copy(aff5)).aff.constant == -42
        @test (copy(quad4) + copy(aff5)).aff.terms[pt] == -38
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + copy(aff5)).terms[pair] == 4
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - copy(aff1), GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - copy(aff2), GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) - copy(aff5), GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad1) - copy(aff3), GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - copy(aff1)).aff.constant == 2
        @test (copy(quad1) - copy(aff1)).aff.terms[pt] == 1
        @test (copy(quad1) - copy(aff1)).aff.terms[inf] == -1
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) - copy(aff1)).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad1) - copy(aff1)).terms[pair] == -1
        @test (copy(quad1) - copy(aff3)).aff.constant == -1
        @test (copy(quad1) - copy(aff3)).aff.terms[meas] == 1
        @test (copy(quad1) - copy(aff3)).aff.terms[glob] == -3
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) - copy(aff3)).terms[pair] == 1
        # test measurefinite combos
        @test (copy(quad2) - copy(aff3)).aff.constant == -4
        @test (copy(quad2) - copy(aff3)).aff.terms[pt] == 1
        @test (copy(quad2) - copy(aff3)).aff.terms[glob] == -1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) - copy(aff3)).terms[pair] == 1
        @test (copy(quad4) - copy(aff3)).aff.constant == -1
        @test (copy(quad4) - copy(aff3)).aff.terms[pt] == -42
        @test (copy(quad4) - copy(aff3)).aff.terms[meas] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4)- copy(aff3)).terms[pair] == 4
        # test finite combos
        @test (copy(quad3) - copy(aff2)).aff.constant == 0
        @test (copy(quad3) - copy(aff2)).aff.terms[pt] == 1
        @test (copy(quad3) - copy(aff2)).aff.terms[glob] == -1
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) - copy(aff2)).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - copy(aff5), GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) - copy(aff5)).aff.constant == 42
        @test (copy(quad4) - copy(aff5)).aff.terms[pt] == -46
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - copy(aff5)).terms[pair] == 4
    end
    # test multiplication
    @testset "Base.:*" begin
        # test some combos
        @test_throws ErrorException quad1 * aff1
        @test_throws ErrorException quad3 * aff2
        @test_throws ErrorException quad4 * aff3
        @test_throws ErrorException quad2 * aff5
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException quad1 / aff1
        @test_throws ErrorException quad3 / aff2
        @test_throws ErrorException quad4 / aff3
        @test_throws ErrorException quad2 / aff5
    end
end

# Test operators on GenericQuadExpr--GenericQuadExpr
@testset "QuadExpr--QuadExpr" begin
    # initialize model and references and expressions
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @global_variable(m, glob)
    @global_variable(m, test[1:51])
    meas = MeasureRef(m, 1)
    m.meas_to_name[1] = "meas"
    aff1 = (inf + pt) - 2
    aff2 = glob - pt
    aff3 = (glob - meas) + 1
    aff5 = 4pt - 42
    quad1 = copy(aff1) * aff2 #(inf*glob - inf*pt + pt*glob - pt^2 - 2glob + 2pt)
    quad2 = (pt * aff3) - 3   #(pt*glob - pt*meas + pt - 3)
    quad3 = glob * aff2       #(glob^2 - glob*pt)
    quad4 = pt * aff5         #(4pt^2 - 42pt)
    # test addition
    @testset "Base.:+" begin
        # test mixed combo types
        @test isa(copy(quad1) + quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) + quad4, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) + quad4, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad4) + quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) + quad2).aff.constant == -3
        @test (copy(quad1) + quad2).aff.terms[pt] == 3
        @test (copy(quad1) + quad2).aff.terms[glob] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) + quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad1) + quad2).terms[pair] == 2
        @test (copy(quad4) + quad1).aff.constant == 0
        @test (copy(quad4) + quad1).aff.terms[pt] == -40
        @test (copy(quad4) + quad1).aff.terms[glob] == -2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + quad1).terms[pair] == 3
        # test measurefinite combos
        @test (copy(quad2) + quad4).aff.constant == -3
        @test (copy(quad2) + quad4).aff.terms[pt] == -41
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) + quad4).terms[pair] == 1
        # test finite combos
        @test (copy(quad3) + quad4).aff.constant == 0
        @test (copy(quad3) + quad4).aff.terms[pt] == -42
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) + quad4).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) + quad4, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) + quad4).aff.constant == 0
        @test (copy(quad4) + quad4).aff.terms[pt] == -84
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) + quad4).terms[pair] == 8
    end
    # test subtraction
    @testset "Base.:-" begin
        # test mixed combo types
        @test isa(copy(quad1) - quad2, GenericQuadExpr{Float64, GeneralVariableRef})
        @test isa(copy(quad2) - quad4, GenericQuadExpr{Float64, MeasureFiniteVariableRef})
        @test isa(copy(quad3) - quad4, GenericQuadExpr{Float64, FiniteVariableRef})
        @test isa(copy(quad4) - quad1, GenericQuadExpr{Float64, GeneralVariableRef})
        # test general combos
        @test (copy(quad1) - quad2).aff.constant == 3
        @test (copy(quad1) - quad2).aff.terms[pt] == 1
        @test (copy(quad1) - quad2).aff.terms[glob] == -2
        pair = UnorderedPair{GeneralVariableRef}(inf, glob)
        @test (copy(quad1) - quad2).terms[pair] == 1
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad1) - quad2).terms[pair] == 0
        @test (copy(quad4) - quad1).aff.constant == 0
        @test (copy(quad4) - quad1).aff.terms[pt] == -44
        @test (copy(quad4) - quad1).aff.terms[glob] == 2
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - quad1).terms[pair] == 5
        # test measurefinite combos
        @test (copy(quad2) - quad4).aff.constant == -3
        @test (copy(quad2) - quad4).aff.terms[pt] == 43
        pair = UnorderedPair{GeneralVariableRef}(pt, glob)
        @test (copy(quad2) - quad4).terms[pair] == 1
        # test finite combos
        @test (copy(quad3) - quad4).aff.constant == 0
        @test (copy(quad3) - quad4).aff.terms[pt] == 42
        pair = UnorderedPair{GeneralVariableRef}(glob, pt)
        @test (copy(quad3) - quad4).terms[pair] == -1
        # test same types
        @test isa(copy(quad4) - quad4, GenericQuadExpr{Float64, PointVariableRef})
        @test (copy(quad4) - quad4).aff.constant == 0
        @test (copy(quad4) - quad4).aff.terms[pt] == 0
        pair = UnorderedPair{GeneralVariableRef}(pt, pt)
        @test (copy(quad4) - quad4).terms[pair] == 0
    end
    # test multiplication
    @testset "Base.:*" begin
        # test some combos
        @test_throws ErrorException quad1 * quad2
        @test_throws ErrorException quad3 * quad1
        @test_throws ErrorException quad4 * quad4
        @test_throws ErrorException quad2 * quad3
    end
    # test division
    @testset "Base.:/" begin
        # test some combos
        @test_throws ErrorException quad1 / quad2
        @test_throws ErrorException quad3 / quad1
        @test_throws ErrorException quad4 / quad4
        @test_throws ErrorException quad2 / quad3
    end
end
