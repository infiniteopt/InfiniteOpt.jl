# Test methods for reduced infinite variables
@testset "Reduced Variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_variable(m, 0 <= inf1(par1, par2) <= 1, Int)
    @infinite_variable(m, inf2(par1, par2) == 1, Bin, start = 0)
    index = m.next_var_index + 1
    rvref1 = ReducedInfiniteVariableRef(m, index)
    m.reduced_info[index] = ReducedInfiniteInfo(inf1, Dict(1 => 1))
    rvref2 = ReducedInfiniteVariableRef(m, index + 1)
    m.reduced_info[index + 1] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
    # test _reduced_info
    @testset "_reduced_info" begin
        @test InfiniteOpt._reduced_info(rvref1).infinite_variable_ref == inf1
        @test InfiniteOpt._reduced_info(rvref2).infinite_variable_ref == inf2
        @test InfiniteOpt._reduced_info(rvref1).eval_supports == Dict(1 => 1)
        @test InfiniteOpt._reduced_info(rvref2).eval_supports == Dict(2 => 1)
    end
    # test infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(rvref1) == inf1
        @test infinite_variable_ref(rvref2) == inf2
    end
    # test eval_supports
    @testset "eval_supports" begin
        @test eval_supports(rvref1) == Dict(1 => 1)
        @test eval_supports(rvref2) == Dict(2 => 1)
    end
    # test name
    @testset "JuMP.name" begin
        @test name(rvref1) == "inf1(1, par2)"
        @test name(rvref2) == "inf2(par1, 1)"
    end
    # parameter_refs (reduced infinite variable)
    @testset "parameter_refs (reduced infinite)" begin
        # test the references
        @test parameter_refs(rvref1) == (par2, )
        @test parameter_refs(rvref2) == (par1, )
    end
    # has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(rvref1)
        @test !has_lower_bound(rvref2)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(rvref1) == 0
        @test_throws ErrorException lower_bound(rvref2)
    end
    # lower_bound_index
    @testset "JuMP.lower_bound_index" begin
        @test lower_bound_index(rvref1) == lower_bound_index(inf1)
        @test_throws ErrorException lower_bound_index(rvref2)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        @test LowerBoundRef(rvref1) == LowerBoundRef(inf1)
        @test_throws ErrorException LowerBoundRef(rvref2)
    end
    # has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(rvref1)
        @test !has_upper_bound(rvref2)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(rvref1) == 1
        @test_throws ErrorException upper_bound(rvref2)
    end
    # upper_bound_index
    @testset "JuMP.upper_bound_index" begin
        @test upper_bound_index(rvref1) == upper_bound_index(inf1)
        @test_throws ErrorException upper_bound_index(rvref2)
    end
    # UpperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        @test UpperBoundRef(rvref1) == UpperBoundRef(inf1)
        @test_throws ErrorException UpperBoundRef(rvref2)
    end
    # is_fixed
    @testset "JuMP.is_fixed" begin
        @test is_fixed(rvref2)
        @test !is_fixed(rvref1)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(rvref2) == 1
        @test_throws ErrorException fix_value(rvref1)
    end
    # fix_index
    @testset "JuMP.fix_index" begin
        @test fix_index(rvref2) == fix_index(inf2)
        @test_throws ErrorException fix_index(rvref1)
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        @test FixRef(rvref2) == FixRef(inf2)
        @test_throws ErrorException FixRef(rvref1)
    end
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(rvref2) == 0
    end
    # is_binary
    @testset "JuMP.is_binary" begin
        @test is_binary(rvref2)
        @test !is_binary(rvref1)
    end
    # binary_index
    @testset "JuMP.binary_index" begin
        @test binary_index(rvref2) == binary_index(inf2)
        @test_throws ErrorException binary_index(rvref1)
    end
    # BinaryRef
    @testset "JuMP.BinaryRef" begin
        @test BinaryRef(rvref2) == BinaryRef(inf2)
        @test_throws ErrorException BinaryRef(rvref1)
    end
    # is_integer
    @testset "JuMP.is_integer" begin
        @test is_integer(rvref1)
        @test !is_integer(rvref2)
    end
    # integer_index
    @testset "JuMP.integer_index" begin
        @test integer_index(rvref1) == integer_index(inf1)
        @test_throws ErrorException integer_index(rvref2)
    end
    # IntegerRef
    @testset "JuMP.IntegerRef" begin
        @test IntegerRef(rvref1) == IntegerRef(inf1)
        @test_throws ErrorException IntegerRef(rvref2)
    end
    # used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(rvref1)
        # prepare use case
        m.reduced_to_constrs[JuMP.index(rvref1)] = [1]
        # test used
        @test used_by_constraint(rvref1)
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(rvref1)
        # prepare use case
        m.reduced_to_meas[JuMP.index(rvref1)] = [1]
        # test used
        @test used_by_measure(rvref1)
    end
    # test is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, rvref1)
        @test is_valid(m, rvref2)
        @test !is_valid(m, ReducedInfiniteVariableRef(m, 100))
    end
end

# Test variable reference makers
@testset "Variable Reference Makers" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    # test _make_point_variable
    @testset "_make_point_variable" begin
        index = m.next_var_index + 1
        @test InfiniteOpt._make_point_variable(inf1) == PointVariableRef(m,
                                                                         index)
    end
    # test _make_reduced_variable (from ivref)
    @testset "_make_reduced_variable (from ivref)" begin
        index = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, index)
        @test InfiniteOpt._make_reduced_variable(inf2, 1, 1) == rvref
        @test m.reduced_info[index].infinite_variable_ref == inf2
        @test m.reduced_info[index].eval_supports == Dict(1 => 1)
    end
    # test _make_reduced_variable (from rvref)
    @testset "_make_reduced_variable (from rvref)" begin
        index = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, index)
        @test InfiniteOpt._make_reduced_variable(inf2, Dict(1 => 1)) == rvref
        @test m.reduced_info[index].infinite_variable_ref == inf2
        @test m.reduced_info[index].eval_supports == Dict(1 => 1)
    end
end

# Test _get_param_value_list
@testset "_get_param_value_list" begin
    # initialize model and data
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 0 <= pars1[1:2] <= 2)
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 0.5], [2, 1.5]])
    # test with DiscreteMeasureData
    @testset "DiscreteMeasureData" begin
        @test InfiniteOpt._get_param_value_list(par1, data1) == [1, 2]
    end
    # test with MultiDiscreteMeasureData
    @testset "MultiDiscreteMeasureData" begin
        @test InfiniteOpt._get_param_value_list(pars1[1], data2) == [1, 2]
        @test InfiniteOpt._get_param_value_list(pars1[2], data2) == [0.5, 1.5]
    end
end

# Test measure expansion methods
@testset "_expand_measure" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_parameter(m, 1 <= pars2[1:2] <= 2)
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @infinite_variable(m, inf5(pars1, pars2))
    @infinite_variable(m, inf6(pars2))
    @infinite_variable(m, inf7(par1, par2, pars1))
    @global_variable(m, x)
    # prepare measures
    w(t) = 3
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(par2, [0.5, 0.5], [1, 1], weight_function = w)
    data3 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data4 = DiscreteMeasureData(pars2, [2, 2], [[1, 1], [1, 1]])
    meas1 = measure(inf1, data1)
    meas2 = measure(2 * inf3 * x - par1, data2)
    # prepare point_mapper function and info
    function mapper(m2, pvref, ivref, support)
        set_name(pvref, string(InfiniteOpt._root_name(ivref), support))
        return
    end
    trans_model = Model()
    map_args = (trans_model, mapper)
    # test _expand_measure (infinite variable)
    @testset "Infinite Variable" begin
        # test single param infinite var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(inf1, data1, map_args...) == expected
        # test single param infinite var without measure param
        expected = 3inf1
        @test InfiniteOpt._expand_measure(inf1, data2, map_args...) == expected
        # test single param infinite var with measure param and others
        index = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (rv1 + rv2)
        @test InfiniteOpt._expand_measure(inf2, data1, map_args...) == expected
        # test array param infinite var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(inf4, data3, map_args...) == expected
        # test array param infinite var without measure param
        expected = 4inf4
        @test InfiniteOpt._expand_measure(inf4, data4, map_args...) == expected
        # test array param infinite var with measure param and others
        index = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf5, Dict(1 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf5, Dict(1 => supp2))
        expected = rv1 + rv2
        @test InfiniteOpt._expand_measure(inf5, data3, map_args...) == expected
    end
    # test _expand_measure (reduced infinite variable)
    @testset "Reduced Variable" begin
        # test single param reduced var without measure param
        rv = ReducedInfiniteVariableRef(m, -1)
        m.reduced_info[-1] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        expected = 0.5 * (rv + rv)
        @test InfiniteOpt._expand_measure(rv, data1, map_args...) == expected
        # test single param reduced var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 1.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data2, map_args...) == expected
        # test single param reduced var with measure param and others
        rv = ReducedInfiniteVariableRef(m, -2)
        m.reduced_info[-2] = ReducedInfiniteInfo(inf7, Dict(1 => 1))
        index = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf7, Dict(1 => 1, 2 => 1))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf7,
                                                         Dict(1 => 2, 2 => 1))
        expected = 1.5 * (rv1 + rv2)
        @test InfiniteOpt._expand_measure(rv, data2, map_args...) == expected
        # test array param reduced var without measure param
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        rv = ReducedInfiniteVariableRef(m, -1)
        m.reduced_info[-1] = ReducedInfiniteInfo(inf5, Dict(1 => supp1))
        expected = rv + rv
        @test InfiniteOpt._expand_measure(rv, data3, map_args...) == expected
        # test array param reduced var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 2 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data4, map_args...) == expected
        # test array param reduced var with measure param and others
        rv = ReducedInfiniteVariableRef(m, -2)
        m.reduced_info[-2] = ReducedInfiniteInfo(inf7, Dict(1 => 1))
        index = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf7,
                                                     Dict(1 => 1, 3 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf7,
                                                       Dict(1 => 2, 3 => supp2))
        expected = rv1 + rv2
        @test InfiniteOpt._expand_measure(rv, data3, map_args...) == expected
    end
    # test _expand_measure (finite variable)
    @testset "Finite Variable" begin
        # test with single parameter measure
        expected = 0.5 * (x + x)
        @test InfiniteOpt._expand_measure(x, data1, map_args...) == expected
        # test with multi parameter measure
        expected = x + x
        @test InfiniteOpt._expand_measure(x, data3, map_args...) == expected
    end
    # test _expand_measure (parameter)
    @testset "Parameter" begin
        # test with different parameter
        expected = 1.5 * (par1 + par1)
        @test InfiniteOpt._expand_measure(par1, data2, map_args...) == expected
        # test with same parameter
        expected = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        expected.constant = 0.5 * (1 + 2)
        @test InfiniteOpt._expand_measure(par1, data1, map_args...) == expected
    end
    # test _expand_measure (multi measure data with parameter)
    @testset "MultiMeasureData Parameter" begin
        # test with different parameter
        expected = par1 + par1
        @test InfiniteOpt._expand_measure(par1, data3, map_args...) == expected
        # test with same parameter
        expected = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        expected.constant = 1 + 2
        @test InfiniteOpt._expand_measure(pars1[2], data3, map_args...) == expected
    end
    # test _expand_measure (AffExpr)
    @testset "AffExpr" begin
        # test single param AffExpr, no measures
        expr = 2inf1 + par1 - x + 3par2 - 3
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (2pts[1] + 2pts[2] + 1 + 2 - x - x + 3par2 + 3par2 - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (6 * (pts[1] * x + pts[2] * x) - 2x - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2] + 2par1 - 2x + 9 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1),
               PointVariableRef(m, index + 2), PointVariableRef(m, index + 3)]
        expected = 6 * (pts[1] * x + pts[2] * x) - 4par1 - 2x - pts[3] - pts[4]
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
    end
    # test _expand_measure (QuadExpr)
    @testset "QuadExpr" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        rv1 = ReducedInfiniteVariableRef(m, index + 2)
        m.reduced_info[index + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, index + 3)
        m.reduced_info[index + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf3 *
                          inf4 + 2x + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (3 * inf3 * pts[1] + 3 * inf3 * pts[2] + 2pars1[1] - 2)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        index = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
        expected = 1.5 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data2, map_args...) == expected
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, index + 2)
        m.reduced_info[index + 2] = ReducedInfiniteInfo(inf5, Dict(1 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, index + 3)
        m.reduced_info[index + 3] = ReducedInfiniteInfo(inf5, Dict(1 => supp2))
        expected = 2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf1 * inf2 + 2x + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 3 * inf1 * pts[1] + 3 * inf1 * pts[2] + 2par1 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        index = m.next_var_index + 1
        supp = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        rv1 = ReducedInfiniteVariableRef(m, index)
        m.reduced_info[index] = ReducedInfiniteInfo(inf5, Dict(2 => supp))
        rv2 = ReducedInfiniteVariableRef(m, index + 1)
        m.reduced_info[index + 1] = ReducedInfiniteInfo(inf5, Dict(2 => supp))
        expected = 2 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data4, map_args...) == expected
    end
    # test _expand_measure (QuadExpr) with parameter quadratics that become numbers
    @testset "QuadExpr with Number Quadratics" begin
        # test single parameter with first quadratic term becomes number
        expr = par1 * inf1 + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...)== expected
        # test single parameter with first quadratic term becomes number, 2nd constant
        expr = par1 * x + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 0.5 * (x + 2 * x + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single parameter with both quadratic terms become numbers
        expr = par1 * par1 + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 0.5 * (1 + 2 * 2 + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single parameter with 2nd quadratic term becomes number
        expr = inf1 * par1 + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * par1 + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 0.5 * (x + 2 * x + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test array parameter with first quadratic term becomes number
        expr = pars1[1] * inf4 + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = pts[1] + 2 * pts[2] + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array parameter with first quadratic term becomes number, 2nd constant
        expr = pars1[1] * x + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = x + 2 * x + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array parameter with both quadratic terms become numbers
        expr = pars1[1] * pars1[2] + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = 1 + 2 * 2 + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array parameter with 2nd quadratic term becomes number
        expr = inf4 * pars1[1] + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = pts[1] + 2 * pts[2] + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * pars1[1] + 2
        expected = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        expected.aff = x + 2 * x + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
    end
    # test _expand_measure (measure)
    @testset "Measure" begin
        # test single param measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(meas1, data1, map_args...) == expected
        # test another single param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (6 * pts[1] * x + 6 * pts[2] * x - 3 * (1 + 2))
        @test InfiniteOpt._expand_measure(meas2, data1, map_args...) == expected
        # test array param measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(meas1, data3, map_args...) == expected
    end
    # test _expand_measure (other)
    @testset "Other" begin
        # prepare test
        @variable(trans_model, y)
        struct new <: AbstractMeasureData end
        # test it
        @test_throws ErrorException InfiniteOpt._expand_measure(y, new(),
                                                                map_args...)
    end
end

# Test user expansion methods
@testset "User Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_variable(m, inf1(par1) >= 1)
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @global_variable(m, x)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test _add_mapped_point_variable
    @testset "_add_mapped_point_variable" begin
        # prepare the partially made point variable reference
        pref = InfiniteOpt._make_point_variable(inf1)
        # try mapping/adding it
        @test isa(InfiniteOpt._add_mapped_point_variable(m, pref, inf1, (1,)),
                  Nothing)
        @test name(pref) == "inf1(1)"
        @test haskey(m.vars, JuMP.index(pref))
        @test lower_bound(pref) == 1
    end
    # test expand
    @testset "expand" begin
        # test the first measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        rv1 = ReducedInfiniteVariableRef(m, index + 2)
        m.reduced_info[index + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, index + 3)
        m.reduced_info[index + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand(meas1) == expected
        # test the second measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand(meas2) == expected
        # test error
        m2 = Model()
        @variable(m2, y)
        struct new <: AbstractMeasureData end
        m.measures[1] = Measure(y, new())
        @test_throws ErrorException expand(meas1)
    end
    # TODO test _expand_measures
    # TODO test expand_all_measures!
end
