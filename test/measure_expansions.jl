# Test variable reference makers
@testset "Variable Reference Makers/Deleters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [0, 2])
    @infinite_parameter(m, par2 in [0, 2])
    @variable(m, inf1, Infinite(par1))
    @variable(m, inf2, Infinite(par1, par2))
    d1 = deriv(inf1, par1)
    d2 = deriv(inf2, par1)
    # test make_point_variable_ref (InfiniteModel)
    @testset "make_point_variable_ref (InfiniteModel)" begin
        # test with inf1
        pvref = GeneralVariableRef(m, 1, PointVariableIndex)
        @test isequal(make_point_variable_ref(m, inf1, Float64[0]), pvref)
        @test parameter_values(pvref) == (0.,)
        @test isequal(make_point_variable_ref(m, inf1, Float64[0]), pvref)
        @test parameter_values(pvref) == (0.,)
        # test with inf2
        pvref = GeneralVariableRef(m, 2, PointVariableIndex)
        @test isequal(make_point_variable_ref(m, inf2, Float64[0, 0]), pvref)
        @test parameter_values(pvref) == (0., 0.)
        @test isequal(make_point_variable_ref(m, inf2, Float64[0, 0]), pvref)
        @test parameter_values(pvref) == (0., 0.)
        # test with d1
        pvref = GeneralVariableRef(m, 3, PointVariableIndex)
        @test isequal(make_point_variable_ref(m, d1, Float64[0]), pvref)
        @test parameter_values(pvref) == (0.,)
        @test isequal(make_point_variable_ref(m, d1, Float64[0]), pvref)
        @test parameter_values(pvref) == (0.,)
        # test with d2
        pvref = GeneralVariableRef(m, 4, PointVariableIndex)
        @test isequal(make_point_variable_ref(m, d2, Float64[0, 0]), pvref)
        @test parameter_values(pvref) == (0., 0.)
        @test isequal(make_point_variable_ref(m, d2, Float64[0, 0]), pvref)
        @test parameter_values(pvref) == (0., 0.)
    end
    # test make_point_variable_ref with infinite parameter functions 
    @testset "make_point_variable_ref (Parameter Function)" begin
        f = parameter_function(sin, par1)
        @test make_point_variable_ref(m, f, [0.]) == f(0)
    end
    # test add_point_variable
    @testset "add_point_variable" begin
        @test_throws ErrorException add_point_variable(TestBackend(), d1, [0.])
    end
    # test make_point_variable_ref (backend)
    @testset "make_point_variable_ref (backend)" begin
        @test_throws ErrorException make_point_variable_ref(TestBackend(), inf1, Float64[0])
    end
    # test make_semi_infinite_variable_ref (InfiniteModel)
    @testset "make_semi_infinite_variable_ref (InfiniteModel)" begin
        # test first addition
        rvref1 = GeneralVariableRef(m, 1, SemiInfiniteVariableIndex)
        @test isequal(make_semi_infinite_variable_ref(m, inf2, [1., NaN]), rvref1)
        @test isequal(infinite_variable_ref(rvref1), inf2)
        @test eval_support(rvref1)[1] == 1
        @test isequal(make_semi_infinite_variable_ref(m, inf2, [1., NaN]), rvref1)
        @test isequal(infinite_variable_ref(rvref1), inf2)
        @test isnan(eval_support(rvref1)[2])
        # test second addition
        rvref2 = GeneralVariableRef(m, 2, SemiInfiniteVariableIndex)
        @test isequal(make_semi_infinite_variable_ref(m, inf2, [0., NaN]), rvref2)
        @test isequal(infinite_variable_ref(rvref2), inf2)
        @test eval_support(rvref2)[1] == 0
        # test first addition with derivative
        rvref1 = GeneralVariableRef(m, 3, SemiInfiniteVariableIndex)
        @test isequal(make_semi_infinite_variable_ref(m, d2, [1, NaN]), rvref1)
        @test isequal(infinite_variable_ref(rvref1), d2)
        @test eval_support(rvref1)[1] == 1
        # test second addition with derivative
        rvref2 = GeneralVariableRef(m, 4, SemiInfiniteVariableIndex)
        @test isequal(make_semi_infinite_variable_ref(m, d2, [0, NaN]), rvref2)
        @test isequal(infinite_variable_ref(rvref2), d2)
        @test isnan(eval_support(rvref2)[2])
    end
    # test add_semi_infinite_variable
    @testset "add_semi_infinite_variable" begin
        @test_throws ErrorException add_semi_infinite_variable(TestBackend(), Bad())
    end
    # test make_semi_infinite_variable_ref (backend)
    @testset "make_semi_infinite_variable_ref (backend)" begin
        @test_throws ErrorException make_semi_infinite_variable_ref(TestBackend(), inf2, [1., NaN])
    end
    # test _process_aff_result
    @testset "_process_aff_result" begin
        @test InfiniteOpt._process_aff_result(42) == 42
        @test isequal(InfiniteOpt._process_aff_result(1par1), par1)
        @test isequal(InfiniteOpt._process_aff_result(2par1 + 3), 2par1 + 3)
    end
end

# Test measure expansion methods
@testset "expand_measure" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [1, 2])
    @infinite_parameter(m, par2 in [1, 2])
    @infinite_parameter(m, pars1[1:2] in [1, 2], independent = true)
    @infinite_parameter(m, pars2[1:2] in [1, 2])
    @variable(m, inf1, Infinite(par1))
    @variable(m, inf2, Infinite(par1, par2))
    @variable(m, inf3, Infinite(par2))
    @variable(m, inf4, Infinite(pars1...))
    @variable(m, inf5, Infinite(pars1..., pars2))
    @variable(m, inf6, Infinite(pars2))
    @variable(m, inf7, Infinite(par1, par2, pars1...))
    @variable(m, inf8, Infinite(pars1[1]))
    @variable(m, inf9, Infinite(pars1[1], par1))
    @variable(m, x)
    d1 = @deriv(inf1, par1)
    d2 = @deriv(inf2, par1)
    d3 = @deriv(inf3, par2)
    d4 = @deriv(inf4, pars1[1])
    d5 = @deriv(inf5, pars1[1])
    # d6 = @deriv(inf6, pars2[2]) # no longer supported
    d7 = @deriv(inf7, par1)
    d8 = @deriv(inf8, pars1[1])
    d9 = @deriv(inf9, pars1[1])
    f1 = parameter_function(sin, par1)
    f2 = parameter_function((a,b) -> 1, (par1, par2))
    f3 = parameter_function(cos, par2)
    f4 = parameter_function((a, b) -> -1, (pars1...,))
    f5 = parameter_function((a,b, c) -> 1, (pars1..., pars2))
    f6 = parameter_function(a -> -1, (pars2,))
    f7 = parameter_function((a,b,c, d) -> 1, (par1, par2, pars1...))
    f8 = parameter_function(cos, pars1[1])
    f9 = parameter_function((a,b) -> 1, (pars1[1], par1))
    # prepare measures
    w(t) = 3
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(par2, [0.5, 0.5], [1, 1], weight_function = w)
    data3 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data4 = DiscreteMeasureData(pars2, [2, 2], [[1, 1], [1, 1]])
    data5 = DiscreteMeasureData([pars1[2], pars1[1]], [1, 1], [[1, 1.2], [2, 2]])
    data6 = DiscreteMeasureData([pars1[2]], [1, 1], [[1], [2]])
    data7 = DiscreteMeasureData(pars1[1], [0.5, 0.5], [1, 2])
    meas1 = measure(inf1, data1)
    meas2 = measure(2 * inf3 * x - par1, data2)
    # test expand_measure (infinite variable) with DiscreteMeasureData
    @testset "Infinite Variable (1D DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        expected = 0.5 * (inf1(1) + inf1(2))
        @test isequal_canonical(InfiniteOpt.expand_measure(inf1, data1, m), expected)
        # test single param infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(inf1, data2, m), 3inf1)
        # test single param infinite var with measure param and others
        expected = 0.5 * (inf2(1, par2) + inf2(2, par2))
        @test isequal_canonical(InfiniteOpt.expand_measure(inf2, data1, m), expected)
        # test multi param infinite var with single element evaluation
        expected = 0.5 * (inf7(par1, par2, 1, pars1[2]) + inf7(par1, par2, 2, pars1[2]))
        @test isequal_canonical(InfiniteOpt.expand_measure(inf7, data7, m), expected)
    end
    # test expand_measure (derivatives) with DiscreteMeasureData
    @testset "Derivative (1D DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        expected = 0.5 * (d1(1) + d1(2))
        @test isequal_canonical(InfiniteOpt.expand_measure(d1, data1, m), expected)
        # test single param infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(d1, data2, m), 3d1)
        # test single param infinite var with measure param and others
        expected = 0.5 * (d2(1, par2) + d2(2, par2))
        @test isequal_canonical(InfiniteOpt.expand_measure(d2, data1, m), expected)
        # test multi param infinite var with single element evaluation
        expected = 0.5 * (d7(par1, par2, 1, pars1[2]) + d7(par1, par2, 2, pars1[2]))
        @test isequal_canonical(InfiniteOpt.expand_measure(d7, data7, m), expected)
    end
    # test expand_measure (infinite parameter functions) with DiscreteMeasureData
    @testset "Parameter Function (1D DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        @test InfiniteOpt.expand_measure(f1, data1, m) == 0.5 * (f1(1) + f1(2))
        # test single param infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(f1, data2, m), 3 * f1)
        # test single param infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, SemiInfiniteVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, SemiInfiniteVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, SemiInfiniteVariableIndex)
        @test isequal_canonical(InfiniteOpt.expand_measure(f2, data1, m), 0.5 * (rv1 + rv2))
        @test eval_support(rv2)[1] == 2
        # test multi param infinite var with single element evaluation
        rv1 = GeneralVariableRef(m, idx + 2, SemiInfiniteVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 3, SemiInfiniteVariableIndex)
        @test isequal_canonical(InfiniteOpt.expand_measure(f7, data7, m), 0.5 * (rv1 + rv2))
        @test eval_support(rv2)[3] == 2
    end
    # test expand_measure (infinite variable) with Multi DiscreteMeasureData
    @testset "Infinite Variable (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        expected = inf4(1, 1) + inf4(2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(inf4, data3, m), expected)
        # test infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(inf4, data4, m), 4inf4)
        # test infinite var with measure param and others
        expected = inf5(1, 1, pars2) + inf5(2, 2, pars2)
        @test isequal_canonical(InfiniteOpt.expand_measure(inf5, data3, m), expected)
        # test infinite var with measure param that is out of order
        expected = inf4(1.2, 1) + inf4(2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(inf4, data5, m), expected)
        # test infinite var with measure param that doesn't overlap
        @test isequal_canonical(InfiniteOpt.expand_measure(inf8, data6, m), 2inf8)
        # test infinite var with measure param that doesn't overlap and there are others
        @test isequal_canonical(InfiniteOpt.expand_measure(inf9, data6, m), 2inf9)
        # test infinite variable has subset of multi params
        expected = inf8(1) + inf8(2)
        @test isequal_canonical(InfiniteOpt.expand_measure(inf8, data3, m), expected)
        # test making semi-infinite vars with variable not containing all the measure prefs
        expected = inf9(1, par1) + inf9(2, par1)
        @test isequal_canonical(InfiniteOpt.expand_measure(inf9, data3, m), expected)
    end
    # test expand_measure (derivatives) with Multi DiscreteMeasureData
    @testset "Derivative (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        expected = d4(1, 1) + d4(2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(d4, data3, m), expected)
        # test infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(d4, data4, m), 4d4)
        # test infinite var with measure param and others
        expected = d5(1, 1, pars2) + d5(2, 2, pars2)
        @test isequal_canonical(InfiniteOpt.expand_measure(d5, data3, m), expected)
        # test infinite var with measure param that is out of order
        expected = d4(1.2, 1) + d4(2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(d4, data5, m), expected)
        # test infinite var with measure param that doesn't overlap
        @test isequal_canonical(InfiniteOpt.expand_measure(d8, data6, m), 2d8)
        # test infinite var with measure param that doesn't overlap and there are others
        @test isequal_canonical(InfiniteOpt.expand_measure(d9, data6, m), 2d9)
        # test infinite variable has subset of multi params
        expected = d8(1) + d8(2)
        @test isequal_canonical(InfiniteOpt.expand_measure(d8, data3, m), expected)
        # test making semi-infinite vars with variable not containing all the measure prefs
        expected = d9(1, par1) + d9(2, par1)
        @test isequal_canonical(InfiniteOpt.expand_measure(d9, data3, m), expected)
    end
    # test expand_measure (infinite parameter functions) with Multi DiscreteMeasureData
    @testset "Parameter Function (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        @test InfiniteOpt.expand_measure(f4, data3, m) == f4(1, 1) + f4(2, 2)
        # test infinite var without measure param
        @test isequal_canonical(InfiniteOpt.expand_measure(f4, data4, m), 4 * f4)
        # test infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, SemiInfiniteVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, SemiInfiniteVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, SemiInfiniteVariableIndex)
        @test isequal_canonical(InfiniteOpt.expand_measure(f5, data3, m), rv1 + rv2)
        @test eval_support(rv2)[[1, 2]] == [2, 2]
        # test infinite var with measure param that is out of order
        @test InfiniteOpt.expand_measure(f4, data5, m) == f4(1.2, 1) + f4(2, 2)
        # test infinite var with measure param that doesn't overlap
        @test isequal_canonical(InfiniteOpt.expand_measure(f8, data6, m), 2 * f8)
        # test infinite var with measure param that doesn't overlap and there are others
        @test isequal_canonical(InfiniteOpt.expand_measure(f9, data6, m), 2 * f9)
        # test infinite variable has subset of multi params
        @test InfiniteOpt.expand_measure(f8, data3, m) == f8(1) + f8(2)
        # test making semi-infinite vars with variable not containing all the measure prefs
        idx = length(InfiniteOpt._data_dictionary(m, SemiInfiniteVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, SemiInfiniteVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, SemiInfiniteVariableIndex)
        @test isequal_canonical(InfiniteOpt.expand_measure(f9, data3, m), rv1 + rv2)
        @test eval_support(rv2)[1] == 2
    end
    # test expand_measure (Semi-infinite variable)
    @testset "Semi-Infinite Variable (1D DiscreteMeasureData)" begin
        # test single param semi-infinite var without measure param
        v = inf2(1, par2) 
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data1, m), 0.5 * (v + v))
        # test single param semi-infinite var with measure param
        expected = 1.5 * (inf2(1, 1) + inf2(1, 1)) 
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data2, m), expected)
        # test single param semi-infinite var with measure param and others
        v = inf7(1, par2, pars1...)
        expected = 1.5 * (inf7(1, 1, pars1...) + inf7(1, 1, pars1...))
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data2, m), expected)
        # test single param semi-infinite var partially from array element
        v = inf7(1, par2, pars1...)
        expected = 0.5 * (inf7(1, par2, 1, pars1[2]) + inf7(1, par2, 2, pars1[2]))
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data7, m), expected)
    end
    # test expand_measure (semi-infinite variable)
    @testset "Semi-Infinite Variable (Multi DiscreteMeasureData)" begin
        # test array param semi-infinite var without measure param
        v = inf5(1, 1, pars2)
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data3, m), v + v)
        # test array param semi-infinite var with measure param
        expected = 2 * (inf5(1, 1, [1, 1]) + inf5(1, 1, [1, 1]))
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data4, m), expected)
        # test array param semi-infinite var with measure param and others
        v = inf7(1, par2, pars1...)
        expected = inf7(1, par2, 1, 1) + inf7(1, par2, 2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data3, m), expected)
        # test array param that is out of order and should make point variable
        v = inf5(pars1..., [1, 1])
        expected = inf5(1.2, 1, [1, 1]) + inf5(2, 2, [1, 1])
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data5, m), expected)
        # test array param that with no others that partially overlap
        v = inf5(pars1[1], 1, [1, 1])
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data6, m), 2v)
        # test array param with others that partially overlap
        v = inf5(pars1[1], 1, pars2)
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data6, m), 2v)
        # test partial overlap
        v = inf5(pars1[1], 1, pars2)
        expected = inf5(1, 1, pars2) + inf5(2, 1, pars2)
        @test isequal_canonical(InfiniteOpt.expand_measure(v, data3, m), expected)
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (1D DiscreteMeasureData)" begin
        # test with single parameter measure
        @test isequal_canonical(InfiniteOpt.expand_measure(x, data1, m), 0.5 * (x + x))
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (Multi DiscreteMeasureData)" begin
        # test with multi parameter measure
        @test isequal_canonical(InfiniteOpt.expand_measure(x, data3, m), x + x)
    end
    # test expand_measure (parameter)
    @testset "Infinite Parameter (1D DiscreteMeasureData)" begin
        # test with different parameter
        @test isequal_canonical(InfiniteOpt.expand_measure(par1, data2, m), 1.5 * (par1 + par1))
        # test with same parameter
        @test InfiniteOpt.expand_measure(par1, data1, m) == 0.5 * (1. + 2.)
    end
    # test expand_measure (multi measure data with parameter)
    @testset "Infinite Parameter (Multi DiscreteMeasureData)" begin
        # test with different parameter
        @test isequal_canonical(InfiniteOpt.expand_measure(par1, data3, m), par1 + par1)
        # test with same parameter
        @test InfiniteOpt.expand_measure(pars1[2], data3, m) == 3.
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (1D DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data2, m), expr * 3)
        # test QuadExpr
        expr = x^2 + 2x +3
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data2, m), expr * 3)
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (Multi DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expr * 2)
        # test QuadExpr
        expr = x^2 + 2x +3
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data4, m), expr * 4)
    end
#     # test expand_measure (AffExpr)
    @testset "AffExpr (1D DiscreteMeasureData)" begin
        # test single param AffExpr, no measures
        expr = 2inf1 + par1 - x + 3par2 - 3
        expected = inf1(1) + inf1(2) - x + 3par2 - 1.5
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        expected = 6 * inf3(1) * x - x - 3
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
    end
    # test expand_measure (AffExpr)
    @testset "AffExpr (Multi DiscreteMeasureData)" begin
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        expected = inf4(1, 1) + inf4(2, 2) + 2par1 - 2x + 7
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        expected = 12 * inf3(1) * x - 4par1 - 2x - inf4(1, 1) - inf4(2, 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (1D DiscreteMeasureData)" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        expected = 0.5 * (2 * inf1(1) * inf2(1, par2) - 2 * inf3 * inf4 + 
                          2 * inf1(2) * inf2(2, par2) + 2x + 4)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        expected = 0.5 * (3 * inf3 * inf1(1) + 3 * inf3 * inf1(2) + 2pars1[1] - 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        expected = 1.5 * (6 * inf2(par1, 1) * inf1 + 2pars1[1] + 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data2, m), expected)
        # test single parameter with first quadratic term becomes number
        expr = par1 * inf1 + 2
        expected = 0.5 * (inf1(1) + 2 * inf1(2) + 4)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected) 
        # test single parameter with first quadratic term becomes number, 2nd constant
        expr = par1 * x + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
        # test single parameter with both quadratic terms become numbers
        expr = par1 * par1 + 2
        expected = 0.5 * (1 + 2 * 2 + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == Float64(expected)
        # test single parameter with 2nd quadratic term becomes number
        expr = inf1 * par1 + 2
        expected = 0.5 * (inf1(1) + 2 * inf1(2) + 4)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
        # test single parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * par1 + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data1, m), expected)
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (Multi DiscreteMeasureData)" begin
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        expected = 2 * inf4(1, 1) * inf5(1, 1, pars2) - 2 * inf1 * inf2 + 
                   2 * inf4(2, 2) * inf5(2, 2, pars2) + 2x + 4
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        expected = 3 * inf1 * inf4(1, 1) + 3 * inf1 * inf4(2, 2) + 2par1 - 2
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        expected = 2 * (6 * inf5(pars1..., [1, 1]) * inf1 + 2pars1[1] + 2)
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data4, m), expected)
        # test array parameter with first quadratic term becomes number
        expr = pars1[1] * inf4 + 2
        expected = inf4(1, 1) + 2 * inf4(2, 2) + 4
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array parameter with first quadratic term becomes number, 2nd constant
        expr = pars1[1] * x + 2
        expected = x + 2 * x + 4
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array parameter with both quadratic terms become numbers
        expr = pars1[1] * pars1[2] + 2
        expected = 1 + 2 * 2 + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == Float64(expected)
        # test array parameter with 2nd quadratic term becomes number
        expr = inf4 * pars1[1] + 2
        expected = inf4(1, 1) + 2 * inf4(2, 2) + 4
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
        # test array parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * pars1[1] + 2
        expected = x + 2 * x + 4
        @test isequal_canonical(InfiniteOpt.expand_measure(expr, data3, m), expected)
    end
    # test expand_measure (NonlinearExpr univariate)
    @testset "GenericNonlinearExpr (1D DiscreteMeasureData)" begin
        # test simple
        expr = sin(inf1)
        expected = 0.5 * sin(inf1(1)) + 0.5 * sin(inf1(2))
        @test isequal(expand_measure(expr, data1, m), expected)
        # test with parameter 
        expr = sin(inf1) + par1
        expected = 0.5 * (sin(inf1(1)) + 1) + 0.5 * (sin(inf1(2)) + 2)
        @test isequal(expand_measure(expr, data1, m), expected)
    end
    # test expand_measure (NonlinearExpr multivariate)
    @testset "GenericNonlinearExpr (Multi DiscreteMeasureData)" begin
        # test simple
        expr = sin(inf5)
        expected = 1 * sin(inf5(1, 1, pars2)) + 1 * sin(inf5(2, 2, pars2))
        @test isequal(expand_measure(expr, data3, m), expected)
    end
    # test expand_measure (measure)
    @testset "Measure" begin
        # test single param measure
        expected = 0.5 * (inf1(1) + inf1(2))
        @test isequal_canonical(InfiniteOpt.expand_measure(meas1, data1, m), expected)
        # test another single param
        expected = 0.5 * (12 * inf3(1) * x - 3 * (1 + 2))
        @test isequal_canonical(InfiniteOpt.expand_measure(meas2, data1, m), expected)
        # test array param measure
        expected = inf1(1) + inf1(2)
        @test isequal_canonical(InfiniteOpt.expand_measure(meas1, data3, m), expected)
    end
    # test expand_measure (FunctionalDiscreteMeasureData)
    @testset "FunctionalDiscreteMeasureData" begin
        # test without generative supports 
        coefa(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(par1, coefa, 2, All, 
                                             is_expect = true)
        @test isequal_canonical(InfiniteOpt.expand_measure(x, data, m), 2x)
        # test with generative supports 
        coefb(a) = ones(length(a) + 1)
        info = UniformGenerativeInfo([0.5], InternalLabel)
        data = FunctionalDiscreteMeasureData(par1, coefb, 2, All, 
                                             generative_support_info = info)
        InfiniteOpt._set_generative_support_info(dispatch_variable_ref(par1), info)
        @test isequal_canonical(InfiniteOpt.expand_measure(x, data, m), 3x)
    end
    # test _prep_generative_supps
    @testset "_prep_generative_supps" begin 
        @test InfiniteOpt._prep_generative_supps(par1, NoGenerativeSupports) isa Nothing 
        @test InfiniteOpt._prep_generative_supps(par1, UniformGenerativeInfo) isa Nothing 
    end
    # test expand_measure (other)
    @testset "Other" begin
        # prepare test
        @variable(Model(), y)
        # test it
        @test_throws ErrorException InfiniteOpt.expand_measure(y, BadData(), m)
    end
end

# Test analytic_expansion
@testset "analytic_expansion" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [1, 2])
    @infinite_parameter(m, par2 in [1, 2])
    @infinite_parameter(m, pars1[1:2] in [1, 2], num_supports = 2)
    @variable(m, inf3, Infinite(par2))
    @variable(m, x)
    # test with 1D data
    @testset "1D Discrete Data" begin
        # test with bounds
        data = DiscreteMeasureData(par1, [1], [1], lower_bound = 1,
                                   upper_bound = 1.5)
        @test isequal_canonical(InfiniteOpt.analytic_expansion(2x + inf3, data, m), 0.5 * (2x + inf3))
        # test as expectation
        coef1(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(par2, coef1, 0, WeightedSample,
                                             is_expect = true)
        @test isequal_canonical(InfiniteOpt.analytic_expansion(inf3, data, m), inf3)
        # test other
        data = DiscreteMeasureData(par1, [1, 0.5], [1, 1.5])
        @test isequal_canonical(InfiniteOpt.analytic_expansion(x, data, m), 1.5 * x)
    end
    # test with multi dimensional data
    @testset "Multi Discrete Data" begin
        # test with bounds
        coef2(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(pars1, coef2, 0, All,
                    lower_bounds = [1, 1], upper_bounds = [1.5, 2])
        @test isequal_canonical(InfiniteOpt.analytic_expansion(x, data, m), 0.5x)
        # test as expectation
        data = DiscreteMeasureData(pars1, [1], [[1, 1]], is_expect = true)
        @test isequal_canonical(InfiniteOpt.analytic_expansion(inf3, data, m), inf3)
        # test other
        coef3(a) = ones(size(a, 2))
        data = FunctionalDiscreteMeasureData(pars1, coef3, 2, All)
        @test isequal_canonical(InfiniteOpt.analytic_expansion(x, data, m), 2x)
    end
    # test with fallback
    @testset "Fallback" begin
        data = TestData(par1, 1, 2)
        @test_throws ErrorException InfiniteOpt.analytic_expansion(x, data, m)
    end
end

# Test expand_measures
@testset "expand_measures" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [1, 2])
    @infinite_parameter(m, par2 in [1, 2])
    @infinite_parameter(m, pars1[1:2] in [1, 2])
    @variable(m, inf1 >= 1, Infinite(par1))
    @variable(m, inf2, Infinite(par1, par2))
    @variable(m, inf3, Infinite(par2))
    @variable(m, inf4, Infinite(pars1))
    @variable(m, x)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2], is_expect = true)
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand_measures (measure)
    @testset "MeasureRef" begin
        # test the first measure
        expected = 0.5 * (inf1(1) + inf1(2) + 6x - inf2(1, par2) - 
                          inf2(2, par2) + 2inf3 - 4)
        @test isequal_canonical(InfiniteOpt.expand_measures(meas1, m), expected)
        # test the second measure
        expected = 2 * inf4([1, 1]) * x + 2 * inf4([2, 2]) * x - 3 + 2inf2
        @test isequal_canonical(InfiniteOpt.expand_measures(meas2, m), expected)
        # test analytic expansion
        meas3 = measure(x, data1)
        @test isequal_canonical(InfiniteOpt.expand_measures(meas3, m), x)
    end
    # test expand_measures (GeneralVariableRef)
    @testset "GeneralVariableRef" begin
        @test isequal_canonical(InfiniteOpt.expand_measures(x, m), x)
        @test isequal_canonical(InfiniteOpt.expand_measures(par1, m), par1)
        @test isequal_canonical(InfiniteOpt.expand_measures(inf1, m), inf1)
    end
    # test expand_measures (Finite Expression)
    @testset "Expression (Finite)" begin
        expr = 2x + 2
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), expr)
        expr = x^2 + 2x + 2
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), expr)
    end
    # test expand_measures (GenericAffExpr)
    @testset "AffExpr (General)" begin
        # test returning AffExpr
        expr = 2inf4 + x + 2meas1
        expected = 2inf4 + 7x + inf1(1) + inf1(2) - inf2(1, par2) - 
                   inf2(2, par2) + 2inf3 - 4
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), expected)
        # test returning QuadExpr
        expr = x + 3par1 + meas2
        expected = x + 3par1 + 2 * inf4([1, 1]) * x + 2 * inf4([2, 2]) * x - 3 + 2inf2
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), expected)
    end
    # test expand_measures (GenericQuadExpr)
    @testset "QuadExpr (General)" begin
        # prepare first measure expansion
        expr = 2 * meas1 * x - x + 1
        m1 = 0.5 * (inf1(1) + inf1(2) + 6x - inf2(1, par2) - inf2(2, par2) + 
                    2inf3 - 4)
        # test with this thorough expression
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), 2 * m1 * x - x + 1)
    end
     # test expand_measures (NLPExpr)
     @testset "NLPExpr (General)" begin
        # prepare first measure expansion
        expr = sin(meas1)
        m1 = 0.5 * (inf1(1) + inf1(2) + 6x - inf2(1, par2) - inf2(2, par2) + 
                    2inf3 - 4)
        # test with this thorough expression
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), sin(m1))
    end
    # test expand_measures (AbstractArray)
    @testset "AbstractArray" begin
        # prepare first measure expansion
        expr = [2x + 1, x]
        @test isequal_canonical(InfiniteOpt.expand_measures(expr, m), expr)
    end
    # test expand_measures (Fallback)
    @testset "Fallback" begin
        @test_throws ErrorException InfiniteOpt.expand_measures(1, m)
    end
end

# Test user expansion methods
@testset "User Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [1, 2])
    @infinite_parameter(m, par2 in [1, 2])
    @infinite_parameter(m, pars1[1:2] in [1, 2])
    @variable(m, inf1 >= 1, Infinite(par1))
    @variable(m, inf2 >= 0, Infinite(par1, par2))
    @variable(m, inf3, Infinite(par2))
    @variable(m, inf4, Infinite(pars1))
    @variable(m, x >= 0)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2], lower_bound = 1, 
                                upper_bound = 2)
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data3 = DiscreteMeasureData(par2, [2, 2], [1, 2])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand
    @testset "expand" begin
        # test the first measure
        expected = 0.5 * (inf1(1) + inf1(2) + 6x - inf2(1, par2) - 
                   inf2(2, par2) + 2inf3 - 4)
        @test isequal_canonical(expand(meas1), expected)
        # test the second measure
        expected = 2 * inf4([1, 1]) * x + 2 * inf4([2, 2]) * x - 3 + 2inf2
        @test isequal_canonical(expand(meas2), expected)
        # test analytic
        meas = measure(x, data1)
        @test is_analytic(meas)
        @test isequal_canonical(expand(meas), 1*x)
    end
    # test expand_all_measures!
    @testset "expand_all_measures!" begin
        # prepare the model for testing
        @variable(m, y)
        set_objective_function(m, x + measure(inf1 * par1 + 3, data1))
        @constraint(m, c1, inf1 + x >= 42.)
        @constraint(m, c2, 2x - measure(measure(inf1 * x + par1 + inf2, data1),
                                        data3) == 0.)
        @constraint(m, c3, measure(inf2, data1) + inf1 >= 0., 
                    DomainRestrictions(par2 => [1, 1.5]))
        @constraint(m, c4, measure(inf4 + y, data2) <= 3., 
                    DomainRestrictions(par2 => 1))
        @constraint(m, c5, [measure(inf2, data1), x] in MOI.Zeros(2))
        # prepare comparison expressions
        obj = x + 0.5inf1(1) + inf1(2) + 3
        c2_expected = 2x - (2 * inf1(1) * x + 2 * inf1(2) * x + inf2(1, 1) + 
                            inf2(1, 2) + inf2(2, 1) + inf2(2, 2))
        c3_expected = 0.5inf2(1, par2) + 0.5inf2(2, par2) + inf1
        c4_expected = inf4([1, 1]) + inf4([2, 2]) + 2y
        c5_expected = AbstractJuMPScalar[0.5inf2(1, par2) + 0.5inf2(2, par2), x]
        # test the expansion
        @test isa(expand_all_measures!(m), Nothing)
        @test isequal_canonical(objective_function(m), obj)
        @test is_valid(m, c2)
        @test is_valid(m, c3)
        @test is_valid(m, c4)
        @test is_valid(m, c5)
        @test name(c1) == "c1"
        @test name(c2) == "c2"
        @test name(c3) == "c3"
        @test name(c4) == "c4"
        @test name(c5) == "c5"
        @test isequal_canonical(constraint_object(c1).func, inf1 + x)
        @test isequal_canonical(constraint_object(c2).func, c2_expected)
        @test isequal_canonical(constraint_object(c3).func, c3_expected)
        @test isequal_canonical(constraint_object(c4).func, c4_expected)
        @test isequal_canonical(constraint_object(c5).func, c5_expected)
        @test constraint_object(c1).set == MOI.GreaterThan(42.)
        @test constraint_object(c2).set == MOI.EqualTo(6.)
        @test constraint_object(c3).set == MOI.GreaterThan(0.)
        @test constraint_object(c4).set == MOI.LessThan(3.)
        @test constraint_object(c5).set == MOI.Zeros(2)
        @test domain_restrictions(c1) == DomainRestrictions()
        @test domain_restrictions(c2) == DomainRestrictions()
        @test domain_restrictions(c3) == DomainRestrictions(par2 => [1, 1.5])
        @test domain_restrictions(c4) == DomainRestrictions(par2 => 1)
        @test domain_restrictions(c5) == DomainRestrictions()
    end
end
