# Test variable reference makers
@testset "Variable Reference Makers/Deleters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [0, 2])
    @infinite_parameter(m, par2 in [0, 2])
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    # test make_point_variable_ref (InfiniteModel)
    @testset "make_point_variable_ref (InfiniteModel)" begin
        # test with inf1
        idx = m.next_var_index + 1
        supp = VectorTuple{Float64}(0)
        pref = PointVariableRef(m, idx)
        @test make_point_variable_ref(m, inf1, supp) == pref
        @test name(pref) == "inf1(0)"
        @test delete(m, pref) isa Nothing
        # test with inf2
        idx += 1
        supp = VectorTuple{Float64}(0, 0)
        pref = PointVariableRef(m, idx)
        @test make_point_variable_ref(m, inf2, supp) == pref
        @test name(pref) == "inf2(0, 0)"
        @test delete(m, pref) isa Nothing
    end
    # test add_measure_variable
    @testset "add_measure_variable" begin
        @test_throws ErrorException add_measure_variable(Model(), Bad(), Val(:some_key))
    end
    # test make_point_variable_ref (optmizer_model)
    @testset "make_point_variable_ref (optimizer_model)" begin
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        supp = VectorTuple{Float64}(0)
        @test_throws ErrorException make_point_variable_ref(opt_m, inf1, supp)
    end
    # test make_point_variable_ref (Tuple input)
    @testset "make_point_variable_ref (Tuple input)" begin
        # test with inf1
        idx = m.next_var_index + 1
        pref = PointVariableRef(m, idx)
        @test make_point_variable_ref(m, inf1, (0,)) == pref
        @test name(pref) == "inf1(0)"
        @test delete(m, pref) isa Nothing
        # test with inf2
        idx += 1
        pref = PointVariableRef(m, idx)
        @test make_point_variable_ref(m, inf2, (0, 0)) == pref
        @test name(pref) == "inf2(0, 0)"
        @test delete(m, pref) isa Nothing
    end
    # test make_reduced_variable_ref (InfiniteModel)
    @testset "make_reduced_variable_ref (InfiniteModel)" begin
        # test first addition
        idx = m.next_reduced_index + 1
        rvref1 = ReducedInfiniteVariableRef(m, idx)
        @test make_reduced_variable_ref(m, inf2, [1], Float64[1]) == rvref1
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => Float64(1))
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx]
        # test second addition
        idx = m.next_reduced_index + 1
        rvref2 = ReducedInfiniteVariableRef(m, idx)
        @test make_reduced_variable_ref(m, inf2, [1], Float64[0]) == rvref2
        @test m.reduced_info[idx].infinite_variable_ref == inf2
        @test m.reduced_info[idx].eval_supports == Dict(1 => Float64(0))
        @test m.infinite_to_reduced[JuMP.index(inf2)] == [idx - 1, idx]
    end
    # test make_reduced_variable_ref (optimizer_model)
    @testset "make_reduced_variable_ref (optimizer_model)" begin
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        @test_throws ErrorException make_reduced_variable_ref(opt_m, inf2, [1], Float64[1])
    end
    # test delete_internal_reduced_variable (InfiniteModel)
    @testset "delete_internal_reduced_variable (InfiniteModel)" begin
        idx = m.next_reduced_index - 1
        rvref1 = ReducedInfiniteVariableRef(m, idx)
        rvref2 = ReducedInfiniteVariableRef(m, idx + 1)
        # test cannot delete
        m.reduced_to_constrs[idx] = [1]
        @test delete_internal_reduced_variable(m, rvref1) isa Nothing
        @test is_valid(m, rvref1)
        delete!(m.reduced_to_constrs, idx)
        # test can delete
        @test delete_internal_reduced_variable(m, rvref1) isa Nothing
        @test delete_internal_reduced_variable(m, rvref2) isa Nothing
        @test !is_valid(m, rvref1)
        @test !is_valid(m, rvref2)
    end
    # test delete_reduced_variable
    @testset "delete_reduced_variable" begin
        rvref = ReducedInfiniteVariableRef(m, m.next_reduced_index - 1)
        warn = "'delete_reduced_variable' not extended for reduced variable type " *
              "`ReducedInfiniteVariableRef` and optimizer model with key `bad`."
        @test_logs (:warn, warn) delete_reduced_variable(Model(), rvref, Val(:bad))
    end
    # test delete_internal_reduced_variable (optimizer_model)
    @testset "delete_internal_reduced_variable (optimizer_model)" begin
        rvref = ReducedInfiniteVariableRef(m, m.next_reduced_index - 1)
        @test !is_valid(m, rvref)
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        warn = "'delete_reduced_variable' not extended for reduced variable type " *
              "`ReducedInfiniteVariableRef` and optimizer model with key `my_key`."
        @test_logs (:warn, warn) delete_internal_reduced_variable(opt_m, rvref)
    end
end

# Test measure expansion methods
@testset "expand_measure" begin
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
    @infinite_variable(m, inf8(pars1[1]))
    @infinite_variable(m, inf9(pars1[1], par1))
    @hold_variable(m, x)
    # prepare measures
    w(t) = 3
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(par2, [0.5, 0.5], [1, 1], weight_function = w)
    data3 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data4 = DiscreteMeasureData(pars2, [2, 2], [[1, 1], [1, 1]])
    data5 = DiscreteMeasureData([pars1[2], pars1[1]], [1, 1], [[1, 1], [2, 2]])
    data6 = DiscreteMeasureData([pars1[2]], [1, 1], [[1], [2]])
    data7 = DiscreteMeasureData(pars1[1], [0.5, 0.5], [1, 2])
    meas1 = measure(inf1, data1)
    meas2 = measure(2 * inf3 * x - par1, data2)
    # test expand_measure (infinite variable) with DiscreteMeasureData
    @testset "Infinite Variable (DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test expand_measure(inf1, data1, m) == expected
        # test single param infinite var without measure param
        expected = 3inf1
        @test expand_measure(inf1, data2, m) == expected
        # test single param infinite var with measure param and others
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        # m.reduced_info[idx] = ReducedInfiniteInfo(inf2, Dict(1 => Float64(1)))
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        # m.reduced_info[idx + 1] = ReducedInfiniteInfo(inf2, Dict(1 => Float64(2)))
        expected = 0.5 * (rv1 + rv2)
        @test expand_measure(inf2, data1, m) == expected
        # test multi param infinite var with single element evaluation
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 0.5 * (rv1 + rv2)
        @test expand_measure(inf7, data7, m) == expected
    end
    # test expand_measure (infinite variable) with MultiDiscreteMeasureData
    @testset "Infinite Variable (MultiDiscreteMeasureData)" begin
        # test infinite var with measure param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2]
        @test expand_measure(inf4, data3, m) == expected
        # test infinite var without measure param
        expected = 4inf4
        @test expand_measure(inf4, data4, m) == expected
        # test infinite var with measure param and others
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = rv1 + rv2
        @test expand_measure(inf5, data3, m) == expected
        # test infinite var with measure param that is out of order
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2]
        @test expand_measure(inf4, data5, m) == expected
        # test infinite var with measure param that doesn't overlap
        @test_throws ErrorException expand_measure(inf8, data6, m)
        # test infinite var with measure param tat doesn't overlap and there are others
        @test_throws ErrorException expand_measure(inf9, data6, m)
        # test dimension mismatch
        @test_throws ErrorException expand_measure(inf8, data3, m)
    end
    @testset "_make_point_support (DiscreteMeasureData)" begin
        orig_prefs = raw_parameter_refs(inf2)
        support_dict = Dict{Int, Float64}(2 => 42)
        expected = VectorTuple{Float64}(23, 42)
        @test InfiniteOpt._make_point_support(orig_prefs, support_dict, 1, Float64(23)) == expected
    end
    # test expand_measure (reduced infinite variable)
    @testset "Reduced Variable (DiscreteMeasureData)" begin
        # test single param reduced var without measure param
        rv = make_reduced_variable_ref(m, inf2, [1], Float64[1])
        expected = 0.5 * (rv + rv)
        @test expand_measure(rv, data1, m) == expected
        # test single param reduced var with measure param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 1.5 * (pts[1] + pts[2])
        @test expand_measure(rv, data2, m) == expected
        @test !is_valid(m, rv)
        # test single param reduced var with measure param and others
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 1.5 * (rv1 + rv2)
        @test expand_measure(rv, data2, m) == expected
        @test !is_valid(m, rv)
        # test single param reduced var partially from array element
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 0.5 * (rv1 + rv2)
        @test expand_measure(rv, data7, m) == expected
        @test !is_valid(m, rv)
    end
    @testset "_make_point_support (MultiDiscreteMeasureData)" begin
        orig_prefs = raw_parameter_refs(inf7)
        support_dict = Dict{Int, Float64}(1 => 23, 2 => 42)
        expected = VectorTuple{Float64}((23, 42, [3, 2]))
        index = 3
        values = Float64[3, 2]
        @test InfiniteOpt._make_point_support(orig_prefs, support_dict, index, values) == expected
    end
    # test expand_measure (reduced infinite variable)
    @testset "Reduced Variable (MultiDiscreteMeasureData)" begin
        # test array param reduced var without measure param
        rv = make_reduced_variable_ref(m, inf5, [1, 2], Float64[1, 1])
        expected = rv + rv
        @test expand_measure(rv, data3, m) == expected
        # test array param reduced var with measure param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * (pts[1] + pts[2])
        @test expand_measure(rv, data4, m) == expected
        @test !is_valid(m, rv)
        # test array param reduced var with measure param and others
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = rv1 + rv2
        @test expand_measure(rv, data3, m) == expected
        @test !is_valid(m, rv)
        # test array param that is out of order and should make point variable
        rv = make_reduced_variable_ref(m, inf5, [3, 4], Float64[1, 1])
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = (pts[1] + pts[2])
        @test expand_measure(rv, data5, m) == expected
        @test !is_valid(m, rv)
        # test array param that with no others that partially overlap
        rv = make_reduced_variable_ref(m, inf5, [2, 3, 4], Float64[1, 1, 1])
        @test_throws ErrorException expand_measure(rv, data6, m)
        # test array param with others that partially overlap
        rv = make_reduced_variable_ref(m, inf5, [2], Float64[1])
        @test_throws ErrorException expand_measure(rv, data6, m)
        # test dimension mismatch
        rv = make_reduced_variable_ref(m, inf5, [2], Float64[1])
        @test_throws ErrorException expand_measure(rv, data3, m)
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (DiscreteMeasureData)" begin
        # test with single parameter measure
        expected = 0.5 * (x + x)
        @test expand_measure(x, data1, m) == expected
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (MultiDiscreteMeasureData)" begin
        # test with multi parameter measure
        expected = x + x
        @test expand_measure(x, data3, m) == expected
    end
    # test expand_measure (parameter)
    @testset "Parameter (DiscreteMeasureData)" begin
        # test with different parameter
        expected = 1.5 * (par1 + par1)
        @test expand_measure(par1, data2, m) == expected
        # test with same parameter
        @test expand_measure(par1, data1, m) == Float64(0.5 * (1 + 2))
    end
    # test expand_measure (multi measure data with parameter)
    @testset "Parameter (MultiDiscreteMeasureData)" begin
        # test with different parameter
        expected = par1 + par1
        @test expand_measure(par1, data3, m) == expected
        # test with same parameter
        @test expand_measure(pars1[2], data3, m) == Float64(3)
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test expand_measure(expr, data2, m) == expr * 3
        # test QuadExpr
        expr = x^2 + 2x +3
        @test expand_measure(expr, data2, m) == expr * 3
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test expand_measure(expr, data3, m) == expr * 2
        # test QuadExpr
        expr = x^2 + 2x +3
        @test expand_measure(expr, data4, m) == expr * 4
    end
    # test expand_measure (AffExpr)
    @testset "AffExpr (DiscreteMeasureData)" begin
        # test single param AffExpr, no measures
        expr = 2inf1 + par1 - x + 3par2 - 3
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (2pts[1] + 2pts[2] + 1 + 2 - x - x + 3par2 + 3par2 - 6)
        @test expand_measure(expr, data1, m) == expected
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (6 * (pts[1] * x + pts[2] * x) - 2x - 6)
        @test expand_measure(expr, data1, m) == expected
    end
    # test expand_measure (AffExpr)
    @testset "AffExpr (MultiDiscreteMeasureData)" begin
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2] + 2par1 - 2x + 9 - 2
        @test expand_measure(expr, data3, m) == expected
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1),
               PointVariableRef(m, idx + 2), PointVariableRef(m, idx + 3)]
        expected = 6 * (pts[1] * x + pts[2] * x) - 4par1 - 2x - pts[3] - pts[4]
        @test expand_measure(expr, data3, m) == expected
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (DiscreteMeasureData)" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 0.5 * (2 * pts[1] * rv1 + - 2 * inf3 * inf4 + 2 * pts[2] *
                          rv2 + 2x + 4)
        @test expand_measure(expr, data1, m) == expected
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (3 * inf3 * pts[1] + 3 * inf3 * pts[2] + 2pars1[1] - 2)
        @test expand_measure(expr, data1, m) == expected
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 1.5 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test expand_measure(expr, data2, m) == expected
        # test single parameter with first quadratic term becomes number
        expr = par1 * inf1 + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test expand_measure(expr, data1, m) == expected
        # test single parameter with first quadratic term becomes number, 2nd constant
        expr = par1 * x + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test expand_measure(expr, data1, m) == expected
        # test single parameter with both quadratic terms become numbers
        expr = par1 * par1 + 2
        expected = 0.5 * (1 + 2 * 2 + 4)
        @test expand_measure(expr, data1, m) == Float64(expected)
        # test single parameter with 2nd quadratic term becomes number
        expr = inf1 * par1 + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test expand_measure(expr, data1, m) == expected
        # test single parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * par1 + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test expand_measure(expr, data1, m) == expected
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (MultiDiscreteMeasureData)" begin
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 2 * pts[1] * rv1 - 2 * inf1 * inf2 + 2 * pts[2] * rv2 + 2x + 4
        @test expand_measure(expr, data3, m) == expected
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 3 * inf1 * pts[1] + 3 * inf1 * pts[2] + 2par1 - 2
        @test expand_measure(expr, data3, m) == expected
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 2 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test expand_measure(expr, data4, m) == expected
        # test array parameter with first quadratic term becomes number
        expr = pars1[1] * inf4 + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + 2 * pts[2] + 4
        @test expand_measure(expr, data3, m) == expected
        # test array parameter with first quadratic term becomes number, 2nd constant
        expr = pars1[1] * x + 2
        expected = x + 2 * x + 4
        @test expand_measure(expr, data3, m) == expected
        # test array parameter with both quadratic terms become numbers
        expr = pars1[1] * pars1[2] + 2
        expected = 1 + 2 * 2 + 4
        @test expand_measure(expr, data3, m) == Float64(expected)
        # test array parameter with 2nd quadratic term becomes number
        expr = inf4 * pars1[1] + 2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + 2 * pts[2] + 4
        @test expand_measure(expr, data3, m) == expected
        # test array parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * pars1[1] + 2
        expected = x + 2 * x + 4
        @test expand_measure(expr, data3, m) == expected
    end
    # test expand_measure (measure)
    @testset "Measure" begin
        # test single param measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test expand_measure(meas1, data1, m) == expected
        # test another single param
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 0.5 * (6 * pts[1] * x + 6 * pts[2] * x - 3 * (1 + 2))
        @test expand_measure(meas2, data1, m) == expected
        # test array param measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = pts[1] + pts[2]
        @test expand_measure(meas1, data3, m) == expected
    end
    # test expand_measure (other)
    @testset "Other" begin
        # prepare test
        @variable(Model(), y)
        # test it
        @test_throws ErrorException expand_measure(y, BadData(), m)
    end
end

# Test expand_measures
@testset "expand_measures" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_variable(m, inf1(par1) >= 1)
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @hold_variable(m, x)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand_measures (measure)
    @testset "MeasureRef" begin
        # test the first measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand_measures(meas1, m) == expected
        # test the second measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand_measures(meas2, m) == expected
        # test error
        m2 = Model()
        @variable(m2, y)
        m.measures[1] = Measure(y, BadData())
        @test_throws ErrorException expand_measures(meas1, m)
        # Undo changes
        m.measures[1] = Measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    end
    # test expand_measures (GeneralVariableRef)
    @testset "GeneralVariableRef" begin
        @test expand_measures(x, m) == x
        @test expand_measures(par1, m) == par1
        @test expand_measures(inf1, m) == inf1
    end
    # test expand_measures (Finite Expression)
    @testset "Expression (Finite)" begin
        expr = 2x + 2
        @test expand_measures(expr, m) == expr
        expr = x^2 + 2x + 2
        @test expand_measures(expr, m) == expr
    end
    # test expand_measures (GenericAffExpr)
    @testset "AffExpr (General)" begin
        # test returning AffExpr
        expr = 2inf4 + x + 2meas1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 2inf4 + 7x + pts[1] + pts[2] - rv1 - rv2 + 2inf3 - 4
        @test expand_measures(expr, m) == expected
        # test returning QuadExpr
        expr = x + 3par1 + meas2
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = x + 3par1 + 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand_measures(expr, m) == expected
    end
    # test expand_measures (GenericQuadExpr)
    @testset "QuadExpr (General)" begin
        # prepare first measure expansion
        expr = 2 * meas1 * x - x + 1
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx2 = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx2)
        rv2 = ReducedInfiniteVariableRef(m, idx2 + 1)
        m1 = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        # test with this thorough expression
        @test expand_measures(expr, m) == 2 * m1 * x - x + 1
    end
    # test expand_measures (Fallback)
    @testset "Fallback" begin
        @test_throws ErrorException expand_measures(1, m)
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
    @hold_variable(m, x >= 0)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data3 = DiscreteMeasureData(par2, [2, 2], [1, 2])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand
    @testset "expand" begin
        # test the first measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        idx = m.next_reduced_index + 1
        rv1 = ReducedInfiniteVariableRef(m, idx)
        rv2 = ReducedInfiniteVariableRef(m, idx + 1)
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand(meas1) == expected
        # test the second measure
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand(meas2) == expected
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
        @BDconstraint(m, c3(par2 in [1, 1.5]), measure(inf2, data1) + inf1 >= 0.)
        @constraint(m, c4, measure(inf4, data2) <= 3.)
        # prepare comparison expressions
        idx = m.next_var_index + 1
        pts = [PointVariableRef(m, idx), PointVariableRef(m, idx + 1)]
        obj = x + 0.5pts[1] + pts[2] + 3
        pts = [PointVariableRef(m, idx + 2), PointVariableRef(m, idx + 3),
               PointVariableRef(m, idx + 4), PointVariableRef(m, idx + 5),
               PointVariableRef(m, idx + 6), PointVariableRef(m, idx + 7)]
        c2_expected = 2x - (2 * pts[1] * x + 2 * pts[2] * x + pts[3] + pts[4] +
                      pts[5] + pts[6])
        idx2 = m.next_reduced_index + 3
        rv1 = ReducedInfiniteVariableRef(m, idx2)
        rv2 = ReducedInfiniteVariableRef(m, idx2 + 1)
        c3_expected = 0.5rv1 + 0.5rv2 + inf1
        pts = [PointVariableRef(m, idx + 8), PointVariableRef(m, idx + 9)]
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
        # @test m.constrs[JuMP.index(c2)].func == c2_expected --> remove for 32 bit testing
        # @test m.constrs[JuMP.index(c3)].func == c3_expected
        # @test m.constrs[JuMP.index(c4)].func == c4_expected
        @test m.constrs[JuMP.index(c2)].func isa JuMP.GenericQuadExpr{Float64, FiniteVariableRef}
        @test m.constrs[JuMP.index(c3)].func isa JuMP.GenericAffExpr{Float64, GeneralVariableRef}
        @test m.constrs[JuMP.index(c4)].func isa JuMP.GenericAffExpr{Float64, PointVariableRef}
        @test m.constrs[JuMP.index(c1)].set == MOI.GreaterThan(42.)
        @test m.constrs[JuMP.index(c2)].set == MOI.EqualTo(6.)
        @test m.constrs[JuMP.index(c3)].set == MOI.GreaterThan(0.)
        @test m.constrs[JuMP.index(c4)].set == MOI.LessThan(3.)
        @test isa(m.constrs[JuMP.index(c1)], ScalarConstraint)
        @test isa(m.constrs[JuMP.index(c2)], ScalarConstraint)
        @test isa(m.constrs[JuMP.index(c3)], BoundedScalarConstraint)
        @test isa(m.constrs[JuMP.index(c4)], ScalarConstraint)
        @test m.constrs[JuMP.index(c3)].bounds.intervals == Dict(par2 => IntervalSet(1, 1.5))
    end
end
