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
        idx = m.next_var_index + 1
        @test InfiniteOpt._make_point_variable(inf1) == PointVariableRef(m,
                                                                         idx)
    end
    # test _make_reduced_variable (from ivref)
    @testset "_make_reduced_variable (from ivref)" begin
        # test first addition
        idx = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, idx)
        @test InfiniteOpt._make_reduced_variable(inf2, 1, 1) == rvref
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => 1)
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx]
        # test second addition
        idx = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, idx)
        @test InfiniteOpt._make_reduced_variable(inf2, 1, 1) == rvref
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => 1)
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx - 1, idx]
        # undo changes
        delete!(m.infinite_to_reduced, JuMP.index(inf2))
    end
    # test _make_reduced_variable (from rvref)
    @testset "_make_reduced_variable (from rvref)" begin
        # test first addition
        idx = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, idx)
        @test InfiniteOpt._make_reduced_variable(inf2, Dict(1 => 1)) == rvref
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => 1)
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx]
        # test second addition
        idx = m.next_var_index + 1
        rvref = ReducedInfiniteVariableRef(m, idx)
        @test InfiniteOpt._make_reduced_variable(inf2, Dict(1 => 1)) == rvref
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => 1)
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx - 1, idx]
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(inf1, data1, map_args...) == expected
        # test single param infinite var without measure param
        expected = 3inf1
        @test InfiniteOpt._expand_measure(inf1, data2, map_args...) == expected
        # test single param infinite var with measure param and others
        idx = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (rv1 + rv2)
        @test InfiniteOpt._expand_measure(inf2, data1, map_args...) == expected
        # test array param infinite var with measure param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(inf4, data3, map_args...) == expected
        # test array param infinite var without measure param
        expected = 4inf4
        @test InfiniteOpt._expand_measure(inf4, data4, map_args...) == expected
        # test array param infinite var with measure param and others
        idx = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf5, Dict(1 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf5, Dict(1 => supp2))
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 1.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data2, map_args...) == expected
        # test single param reduced var with measure param and others
        rv = ReducedInfiniteVariableRef(m, -2)
        m.reduced_info[-2] = ReducedInfiniteInfo(inf7, Dict(1 => 1))
        idx = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf7, Dict(1 => 1, 2 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf7,
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data4, map_args...) == expected
        # test array param reduced var with measure param and others
        rv = ReducedInfiniteVariableRef(m, -2)
        m.reduced_info[-2] = ReducedInfiniteInfo(inf7, Dict(1 => 1))
        idx = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf7,
                                                     Dict(1 => 1, 3 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf7,
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (2pts[1] + 2pts[2] + 1 + 2 - x - x + 3par2 + 3par2 - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (6 * (pts[1] * x + pts[2] * x) - 2x - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2] + 2par1 - 2x + 9 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1),
               PointVariableRef(m, idx + 2), PointVariableRef(m, idx + 3)]
        expected = 6 * (pts[1] * x + pts[2] * x) - 4par1 - 2x - pts[3] - pts[4]
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
    end
    # test _expand_measure (QuadExpr)
    @testset "QuadExpr" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf3 *
                          inf4 + 2x + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (3 * inf3 * pts[1] + 3 * inf3 * pts[2] + 2pars1[1] - 2)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        idx = m.next_var_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
        expected = 1.5 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data2, map_args...) == expected
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf5, Dict(1 => supp1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf5, Dict(1 => supp2))
        expected = 2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf1 * inf2 + 2x + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 3 * inf1 * pts[1] + 3 * inf1 * pts[2] + 2par1 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        idx = m.next_var_index + 1
        supp = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        rv1 = ReducedInfiniteVariableRef(m, idx)
        m.reduced_info[idx] = ReducedInfiniteInfo(inf5, Dict(2 => supp))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf5, Dict(2 => supp))
        expected = 2 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data4, map_args...) == expected
    end
    # test _expand_measure (QuadExpr) with parameter quadratics that become numbers
    @testset "QuadExpr with Number Quadratics" begin
        # test single parameter with first quadratic term becomes number
        expr = par1 * inf1 + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(meas1, data1, map_args...) == expected
        # test another single param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (6 * pts[1] * x + 6 * pts[2] * x - 3 * (1 + 2))
        @test InfiniteOpt._expand_measure(meas2, data1, map_args...) == expected
        # test array param measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(meas1, data3, map_args...) == expected
    end
    # test _expand_measure (other)
    @testset "Other" begin
        # prepare test
        @variable(trans_model, y)
        # test it
        @test_throws ErrorException InfiniteOpt._expand_measure(y, BadData(),
                                                                map_args...)
    end
end

# Test _expand_measures
@testset "_expand_measures" begin
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
    # prepare point_mapper function and info
    function mapper(m2, pvref, ivref, support)
        set_name(pvref, string(InfiniteOpt._root_name(ivref), support))
        return
    end
    trans_model = Model()
    map_args = (trans_model, mapper)
    # test _expand_measures (measure)
    @testset "MeasureRef" begin
        # test the first measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test InfiniteOpt._expand_measures(meas1, map_args...) == expected
        # test the second measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test InfiniteOpt._expand_measures(meas2, map_args...) == expected
        # test error
        m2 = Model()
        @variable(m2, y)
        m.measures[1] = Measure(y, BadData())
        @test_throws ErrorException InfiniteOpt._expand_measures(meas1,
                                                                 map_args...)
        # Undo changes
        m.measures[1] = Measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    end
    # test _expand_measures (GenericAffExpr)
    @testset "AffExpr" begin
        # test returning AffExpr
        expr = 2inf4 + x + 2meas1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 2inf4 + 7x + pts[1] + pts[2] - rv1 - rv2 + 2inf3 - 4
        @test InfiniteOpt._expand_measures(expr, map_args...) == expected
        # test returning QuadExpr
        expr = x + 3par1 + meas2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = x + 3par1 + 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test InfiniteOpt._expand_measures(expr, map_args...) == expected
    end
    # test _expand_measures (GenericQuadExpr)
    @testset "QuadExpr" begin
        # prepare first measure expansion
        expr = 2 * meas1 * x + inf1 * meas1 - x + meas1 - 3
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        m1 = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        # prepare second measure expansion
        pts = [PointVariableRef(m, idx + 4), PointVariableRef(m, idx + 5)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 6)
        m.reduced_info[idx + 6] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 7)
        m.reduced_info[idx + 7] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        m2 = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        # prepare third measure expansion
        pts = [PointVariableRef(m, idx + 8), PointVariableRef(m, idx + 9)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 10)
        m.reduced_info[idx + 10] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 11)
        m.reduced_info[idx + 11] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        m3 = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        # test with this thorough expression
        quad = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        quad.aff = - x + m1 - 3
        expected = quad + 2 * m2 * x + inf1 * m3
        @test InfiniteOpt._expand_measures(expr, map_args...) == expected
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
    @infinite_variable(m, inf2(par1, par2) >= 0)
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @global_variable(m, x >= 0)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data3 = DiscreteMeasureData(par2, [2, 2], [1, 2])
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
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        rv1 = ReducedInfiniteVariableRef(m, idx + 2)
        m.reduced_info[idx + 2] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 3)
        m.reduced_info[idx + 3] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand(meas1) == expected
        # test the second measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand(meas2) == expected
        # test error
        m2 = Model()
        @variable(m2, y)
        m.measures[1] = Measure(y, BadData())
        @test_throws ErrorException expand(meas1)
        # Undo changes
        m.measures[1] = Measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    end
    # test expand_all_measures!
    @testset "expand_all_measures!" begin
        # test error
        @objective(m, Min, x + meas1)
        @test_throws ErrorException expand_all_measures!(m)
        # prepare the model for testing
        set_objective_function(m, x + measure(inf1 * par1 + 3, data1))
        @constraint(m, c1, inf1 + x >= 42.)
        @constraint(m, c2, 2x - measure(measure(inf1 * x + par1 + inf2, data1),
                                        data3) == 0.)
        @constraint(m, c3, measure(inf2, data1) + inf1 >= 0.,
                    parameter_bounds = Dict(par2 => IntervalSet(1, 1.5)))
        @constraint(m, c4, measure(inf4, data2) <= 3.)
        # prepare comparison expressions
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        obj = x + 0.5pts[1] + pts[2] + 3
        pts = [PointVariableRef(m, idx + 4), PointVariableRef(m, idx + 5),
               PointVariableRef(m, idx + 6), PointVariableRef(m, idx + 7),
               PointVariableRef(m, idx + 8), PointVariableRef(m, idx + 9)]
        c2_expected = 2x - (2 * pts[1] * x + 2 * pts[2] * x + pts[3] + pts[4] +
                      pts[5] + pts[6])
        rv1 = ReducedInfiniteVariableRef(m, idx + 10)
        m.reduced_info[idx + 10] = ReducedInfiniteInfo(inf2, Dict(1 => 1))
        rv2 = ReducedInfiniteVariableRef(m, idx + 11)
        m.reduced_info[idx + 11] = ReducedInfiniteInfo(inf2, Dict(1 => 2))
        c3_expected = 0.5rv1 + 0.5rv2 + inf1
        pts = [PointVariableRef(m, idx + 12), PointVariableRef(m, idx + 13)]
        c4_expected = pts[1] + pts[2]
        # test the expansion
        @test isa(expand_all_measures!(m), Nothing)
        @test objective_function(m) == obj
        @test is_valid(m, c2)
        @test is_valid(m, c3)
        @test is_valid(m, c4)
        @test name(c1) == "c1"
        @test name(c2) == "c2"
        @test name(c3) == "c3"
        @test name(c4) == "c4"
        @test m.constrs[JuMP.index(c1)].func == inf1 + x
        @test m.constrs[JuMP.index(c2)].func == c2_expected
        @test m.constrs[JuMP.index(c3)].func == c3_expected
        @test m.constrs[JuMP.index(c4)].func == c4_expected
        @test m.constrs[JuMP.index(c1)].set == MOI.GreaterThan(42.)
        @test m.constrs[JuMP.index(c2)].set == MOI.EqualTo(6.)
        @test m.constrs[JuMP.index(c3)].set == MOI.GreaterThan(0.)
        @test m.constrs[JuMP.index(c4)].set == MOI.LessThan(3.)
        @test isa(m.constrs[JuMP.index(c1)], ScalarConstraint)
        @test isa(m.constrs[JuMP.index(c2)], ScalarConstraint)
        @test isa(m.constrs[JuMP.index(c3)], BoundedScalarConstraint)
        @test isa(m.constrs[JuMP.index(c4)], ScalarConstraint)
        @test m.constrs[JuMP.index(c3)].bounds == Dict(par2 => IntervalSet(1, 1.5))
    end
end
