# Test variable reference makers
@testset "Variable Reference Makers/Deleters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [0, 2])
    @infinite_parameter(m, par2 in [0, 2])
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    d1 = @deriv(inf1, par1)
    d2 = @deriv(inf2, par1)
    # test make_point_variable_ref (InfiniteModel)
    @testset "make_point_variable_ref (InfiniteModel)" begin
        # test with inf1
        pvref = GeneralVariableRef(m, 1, PointVariableIndex)
        @test make_point_variable_ref(m, inf1, Float64[0]) == pvref
        @test parameter_values(pvref) == (0.,)
        # test with inf2
        pvref = GeneralVariableRef(m, 2, PointVariableIndex)
        @test make_point_variable_ref(m, inf2, Float64[0, 0]) == pvref
        @test parameter_values(pvref) == (0., 0.)
        # test with d1
        pvref = GeneralVariableRef(m, 3, PointVariableIndex)
        @test make_point_variable_ref(m, d1, Float64[0]) == pvref
        @test parameter_values(pvref) == (0.,)
        # test with d2
        pvref = GeneralVariableRef(m, 4, PointVariableIndex)
        @test make_point_variable_ref(m, d2, Float64[0, 0]) == pvref
        @test parameter_values(pvref) == (0., 0.)
    end
    # test make_point_variable_ref with infinite parameter functions 
    @testset "make_point_variable_ref (Parameter Function)" begin
        f = parameter_function(sin, par1)
        @test make_point_variable_ref(m, f, [0.]) == 0 
        @test make_point_variable_ref(Model(), f, [0.]) == 0 
    end
    # test add_measure_variable
    @testset "add_measure_variable" begin
        @test_throws ErrorException add_measure_variable(Model(), Bad(), Val(:some_key))
    end
    # test make_point_variable_ref (optmizer_model)
    @testset "make_point_variable_ref (optimizer_model)" begin
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        @test_throws ErrorException make_point_variable_ref(opt_m, inf1, Float64[0])
    end
    # test make_reduced_variable_ref (InfiniteModel)
    @testset "make_reduced_variable_ref (InfiniteModel)" begin
        # test first addition
        rvref1 = GeneralVariableRef(m, 1, ReducedVariableIndex)
        @test make_reduced_variable_ref(m, inf2, [1], Float64[1]) == rvref1
        @test infinite_variable_ref(rvref1) == inf2
        @test eval_supports(rvref1) == Dict(1 => Float64(1))
        # test second addition
        rvref2 = GeneralVariableRef(m, 2, ReducedVariableIndex)
        @test make_reduced_variable_ref(m, inf2, [1], Float64[0]) == rvref2
        @test infinite_variable_ref(rvref2) == inf2
        @test eval_supports(rvref2) == Dict(1 => Float64(0))
        # test first addition with derivative
        rvref1 = GeneralVariableRef(m, 3, ReducedVariableIndex)
        @test make_reduced_variable_ref(m, d2, [1], Float64[1]) == rvref1
        @test infinite_variable_ref(rvref1) == d2
        @test eval_supports(rvref1) == Dict(1 => Float64(1))
        # test second addition with derivative
        rvref2 = GeneralVariableRef(m, 4, ReducedVariableIndex)
        @test make_reduced_variable_ref(m, d2, [1], Float64[0]) == rvref2
        @test infinite_variable_ref(rvref2) == d2
        @test eval_supports(rvref2) == Dict(1 => Float64(0))
    end
    # test make_reduced_variable_ref (optimizer_model)
    @testset "make_reduced_variable_ref (optimizer_model)" begin
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        @test_throws ErrorException make_reduced_variable_ref(opt_m, inf2, [1], Float64[1])
    end
    # test delete_internal_reduced_variable (InfiniteModel)
    @testset "delete_internal_reduced_variable (InfiniteModel)" begin
        rvref1 = ReducedVariableRef(m, ReducedVariableIndex(1))
        rvref2 = ReducedVariableRef(m, ReducedVariableIndex(2))
        rvref3 = ReducedVariableRef(m, ReducedVariableIndex(3))
        rvref4 = ReducedVariableRef(m, ReducedVariableIndex(4))
        # test cannot delete
        push!(InfiniteOpt._constraint_dependencies(rvref1), ConstraintIndex(1))
        @test delete_internal_reduced_variable(m, rvref1) isa Nothing
        @test is_valid(m, rvref1)
        empty!(InfiniteOpt._constraint_dependencies(rvref1))
        # test can delete
        @test delete_internal_reduced_variable(m, rvref1) isa Nothing
        @test delete_internal_reduced_variable(m, rvref2) isa Nothing
        @test delete_internal_reduced_variable(m, rvref3) isa Nothing
        @test delete_internal_reduced_variable(m, rvref4) isa Nothing
        @test !is_valid(m, rvref1)
        @test !is_valid(m, rvref2)
        @test !is_valid(m, rvref3)
        @test !is_valid(m, rvref4)
    end
    # test delete_reduced_variable
    @testset "delete_reduced_variable" begin
        rvref = ReducedVariableRef(m, ReducedVariableIndex(1))
        warn = "'delete_reduced_variable' not extended for reduced variable type " *
              "`ReducedVariableRef` and optimizer model with key `bad`."
        @test_logs (:warn, warn) delete_reduced_variable(Model(), rvref, Val(:bad))
    end
    # test delete_internal_reduced_variable (optimizer_model)
    @testset "delete_internal_reduced_variable (optimizer_model)" begin
        rvref = ReducedVariableRef(m, ReducedVariableIndex(1))
        @test !is_valid(m, rvref)
        opt_m = Model()
        opt_m.ext[:my_key] = 42
        warn = "'delete_reduced_variable' not extended for reduced variable type " *
              "`ReducedVariableRef` and optimizer model with key `my_key`."
        @test_logs (:warn, warn) delete_internal_reduced_variable(opt_m, rvref)
    end
end

# Test measure expansion methods
@testset "expand_measure" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2, independent = true)
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
    d1 = @deriv(inf1, par1)
    d2 = @deriv(inf2, par1)
    d3 = @deriv(inf3, par2)
    d4 = @deriv(inf4, pars1[1])
    d5 = @deriv(inf5, pars1[1])
    d6 = @deriv(inf6, pars2[1])
    d7 = @deriv(inf7, par1)
    d8 = @deriv(inf8, pars1[1])
    d9 = @deriv(inf9, pars1[1])
    f1 = parameter_function(sin, par1)
    f2 = parameter_function((a,b) -> 1, (par1, par2))
    f3 = parameter_function(cos, par2)
    f4 = parameter_function(a -> -1, (pars1,))
    f5 = parameter_function((a,b) -> 1, (pars1, pars2))
    f6 = parameter_function(a -> -1, (pars2,))
    f7 = parameter_function((a,b,c) -> 1, (par1, par2, pars1))
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
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(inf1, data1, m) == 0.5 * (pts[1] + pts[2])
        @test parameter_values.(pts) == [(1.,), (2.,)]
        # test single param infinite var without measure param
        @test InfiniteOpt.expand_measure(inf1, data2, m) == 3inf1
        # test single param infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(inf2, data1, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2.)
        # test multi param infinite var with single element evaluation
        rv1 = GeneralVariableRef(m, idx + 2, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 3, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(inf7, data7, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(3 => 2.)
    end
    # test expand_measure (derivatives) with DiscreteMeasureData
    @testset "Derivative (1D DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(d1, data1, m) == 0.5 * (pts[1] + pts[2])
        @test parameter_values.(pts) == [(1.,), (2.,)]
        # test single param infinite var without measure param
        @test InfiniteOpt.expand_measure(d1, data2, m) == 3d1
        # test single param infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(d2, data1, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2.)
        # test multi param infinite var with single element evaluation
        rv1 = GeneralVariableRef(m, idx + 2, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 3, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(d7, data7, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(3 => 2.)
    end
    # test expand_measure (infinite parameter functions) with DiscreteMeasureData
    @testset "Parameter Function (1D DiscreteMeasureData)" begin
        # test single param infinite var with measure param
        @test InfiniteOpt.expand_measure(f1, data1, m) == 0.5 * (sin(1) + sin(2))
        # test single param infinite var without measure param
        @test InfiniteOpt.expand_measure(f1, data2, m) == 3 * f1
        # test single param infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(f2, data1, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2.)
        # test multi param infinite var with single element evaluation
        rv1 = GeneralVariableRef(m, idx + 2, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 3, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(f7, data7, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(3 => 2.)
    end
    # test expand_measure (infinite variable) with Multi DiscreteMeasureData
    @testset "Infinite Variable (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(inf4, data3, m) == pts[1] + pts[2]
        @test parameter_values(pts[1]) == (Float64[1., 1.],)
        # test infinite var without measure param
        @test InfiniteOpt.expand_measure(inf4, data4, m) == 4inf4
        # test infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(inf5, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2, 2 => 2)
        # test infinite var with measure param that is out of order
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(inf4, data5, m) == pts[1] + pts[2]
        @test parameter_values(pts[1]) == (Float64[1.2, 1],)
        # test infinite var with measure param that doesn't overlap
        @test InfiniteOpt.expand_measure(inf8, data6, m) == 2inf8
        # test infinite var with measure param that doesn't overlap and there are others
        @test InfiniteOpt.expand_measure(inf9, data6, m) == 2inf9
        # test infinite variable has subset of multi params
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(inf8, data3, m) == pts[1] + pts[2]
        @test parameter_values(pts[2]) == (2.,)
        # test making reduced vars with variable not containing all the measure prefs
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(inf9, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2)
    end
    # test expand_measure (derivatives) with Multi DiscreteMeasureData
    @testset "Derivative (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(d4, data3, m) == pts[1] + pts[2]
        @test parameter_values(pts[1]) == (Float64[1., 1.],)
        # test infinite var without measure param
        @test InfiniteOpt.expand_measure(d4, data4, m) == 4d4
        # test infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(d5, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2, 2 => 2)
        # test infinite var with measure param that is out of order
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(d4, data5, m) == pts[1] + pts[2]
        @test parameter_values(pts[1]) == (Float64[1.2, 1],)
        # test infinite var with measure param that doesn't overlap
        @test InfiniteOpt.expand_measure(d8, data6, m) == 2d8
        # test infinite var with measure param that doesn't overlap and there are others
        @test InfiniteOpt.expand_measure(d9, data6, m) == 2d9
        # test infinite variable has subset of multi params
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(d8, data3, m) == pts[1] + pts[2]
        @test parameter_values(pts[2]) == (2.,)
        # test making reduced vars with variable not containing all the measure prefs
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(d9, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2)
    end
    # test expand_measure (infinite parameter functions) with Multi DiscreteMeasureData
    @testset "Parameter Function (Multi DiscreteMeasureData)" begin
        # test infinite var with measure param
        @test InfiniteOpt.expand_measure(f4, data3, m) == -2
        # test infinite var without measure param
        @test InfiniteOpt.expand_measure(f4, data4, m) == 4 * f4
        # test infinite var with measure param and others
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(f5, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2, 2 => 2)
        # test infinite var with measure param that is out of order
        @test InfiniteOpt.expand_measure(f4, data5, m) == -2
        # test infinite var with measure param that doesn't overlap
        @test InfiniteOpt.expand_measure(f8, data6, m) == 2 * f8
        # test infinite var with measure param that doesn't overlap and there are others
        @test InfiniteOpt.expand_measure(f9, data6, m) == 2 * f9
        # test infinite variable has subset of multi params
        @test InfiniteOpt.expand_measure(f8, data3, m) == cos(1) + cos(2)
        # test making reduced vars with variable not containing all the measure prefs
        idx = length(InfiniteOpt._data_dictionary(m, ReducedVariable)) + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(f9, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 2)
    end
    @testset "_make_point_support (1D DiscreteMeasureData)" begin
        orig_prefs = parameter_list(inf2)
        support_dict = Dict{Int, Float64}(2 => 42)
        expected = Float64[23, 42]
        @test InfiniteOpt._make_point_support(orig_prefs, support_dict, 1, Float64(23)) == expected
    end
    # test expand_measure (reduced infinite variable)
    @testset "Reduced Variable (1D DiscreteMeasureData)" begin
        # test single param reduced var without measure param
        rv = make_reduced_variable_ref(m, inf2, [1], Float64[1])
        @test InfiniteOpt.expand_measure(rv, data1, m) == 0.5 * (rv + rv)
        # test single param reduced var with measure param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(rv, data2, m) == 1.5 * (pts[1] + pts[2])
        @test parameter_values(pts[1]) == (1., 1.)
        @test !is_valid(m, rv)
        # test single param reduced var with measure param and others
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = index(rv).value + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(rv, data2, m) == 1.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 1, 2 => 1)
        @test !is_valid(m, rv)
        # test single param reduced var partially from array element
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = index(rv).value + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(rv, data7, m) == 0.5 * (rv1 + rv2)
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 1, 3 => 2)
        @test !is_valid(m, rv)
    end
    @testset "_make_point_support (Multi DiscreteMeasureData)" begin
        orig_prefs = parameter_list(inf7)
        support_dict = Dict{Int, Float64}(1 => 23, 2 => 42)
        expected = Float64[23, 42, 3, 2]
        vals = Float64[3, 2]
        @test InfiniteOpt._make_point_support(orig_prefs, support_dict, [3, 4], vals) == expected
    end
    # test expand_measure (reduced infinite variable)
    @testset "Reduced Variable (Multi DiscreteMeasureData)" begin
        # test array param reduced var without measure param
        rv = make_reduced_variable_ref(m, inf5, [1, 2], Float64[1, 1])
        @test InfiniteOpt.expand_measure(rv, data3, m) == rv + rv
        # test array param reduced var with measure param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(rv, data4, m) == 2 * (pts[1] + pts[2])
        @test parameter_values(pts[2]) == (Float64[1, 1], Float64[1, 1])
        @test !is_valid(m, rv)
        # test array param reduced var with measure param and others
        rv = make_reduced_variable_ref(m, inf7, [1], Float64[1])
        idx = index(rv).value + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(rv, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(1 => 1, 3 => 2, 4 => 2)
        @test !is_valid(m, rv)
        # test array param that is out of order and should make point variable
        rv = make_reduced_variable_ref(m, inf5, [3, 4], Float64[1, 1])
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        @test InfiniteOpt.expand_measure(rv, data5, m) == pts[1] + pts[2]
        @test parameter_values(pts[1]) == (Float64[1.2, 1], Float64[1, 1])
        @test !is_valid(m, rv)
        # test array param that with no others that partially overlap
        rv = make_reduced_variable_ref(m, inf5, [2, 3, 4], Float64[1, 1, 1])
        @test InfiniteOpt.expand_measure(rv, data6, m) == 2rv
        # test array param with others that partially overlap
        rv = make_reduced_variable_ref(m, inf5, [2], Float64[1])
        @test InfiniteOpt.expand_measure(rv, data6, m) == 2rv
        # test partial overlap
        rv = make_reduced_variable_ref(m, inf5, [2], Float64[1])
        idx = index(rv).value + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        @test InfiniteOpt.expand_measure(rv, data3, m) == rv1 + rv2
        @test eval_supports(rv2) == Dict{Int, Float64}(2 => 1, 1 => 2)
        @test !is_valid(m, rv)
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (1D DiscreteMeasureData)" begin
        # test with single parameter measure
        @test InfiniteOpt.expand_measure(x, data1, m) == 0.5 * (x + x)
    end
    # test expand_measure (finite variable)
    @testset "Finite Variable (Multi DiscreteMeasureData)" begin
        # test with multi parameter measure
        @test InfiniteOpt.expand_measure(x, data3, m) == x + x
    end
    # test expand_measure (parameter)
    @testset "Infinite Parameter (1D DiscreteMeasureData)" begin
        # test with different parameter
        @test InfiniteOpt.expand_measure(par1, data2, m) == 1.5 * (par1 + par1)
        # test with same parameter
        @test InfiniteOpt.expand_measure(par1, data1, m) == 0.5 * (1. + 2.)
    end
    # test expand_measure (multi measure data with parameter)
    @testset "Infinite Parameter (Multi DiscreteMeasureData)" begin
        # test with different parameter
        @test InfiniteOpt.expand_measure(par1, data3, m) == par1 + par1
        # test with same parameter
        @test InfiniteOpt.expand_measure(pars1[2], data3, m) == 3.
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (1D DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test InfiniteOpt.expand_measure(expr, data2, m) == expr * 3
        # test QuadExpr
        expr = x^2 + 2x +3
        @test InfiniteOpt.expand_measure(expr, data2, m) == expr * 3
    end
    # test expand_measure (Finite Expr) with DiscreteMeasureData
    @testset "Finite Expression (Multi DiscreteMeasureData)"  begin
        # test AffExpr
        expr = 2x +3
        @test InfiniteOpt.expand_measure(expr, data3, m) == expr * 2
        # test QuadExpr
        expr = x^2 + 2x +3
        @test InfiniteOpt.expand_measure(expr, data4, m) == expr * 4
    end
    # test expand_measure (AffExpr)
    @testset "AffExpr (1D DiscreteMeasureData)" begin
        # test single param AffExpr, no measures
        expr = 2inf1 + par1 - x + 3par2 - 3
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (2pts[1] + 2pts[2] + 1 + 2 - x - x + 3par2 + 3par2 - 6)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (6 * (pts[1] * x + pts[2] * x) - 2x - 6)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
    end
    # test expand_measure (AffExpr)
    @testset "AffExpr (Multi DiscreteMeasureData)" begin
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = pts[1] + pts[2] + 2par1 - 2x + 9 - 2
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex),
               GeneralVariableRef(m, idx + 2, PointVariableIndex),
               GeneralVariableRef(m, idx + 3, PointVariableIndex)]
        expected = 6 * (pts[1] * x + pts[2] * x) - 4par1 - 2x - pts[3] - pts[4]
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (1D DiscreteMeasureData)" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 0.5 * (2 * pts[1] * rv1 + - 2 * inf3 * inf4 + 2 * pts[2] *
                          rv2 + 2x + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (3 * inf3 * pts[1] + 3 * inf3 * pts[2] + 2pars1[1] - 2)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 1.5 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt.expand_measure(expr, data2, m) == expected
        # test single parameter with first quadratic term becomes number
        expr = par1 * inf1 + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected 
        # test single parameter with first quadratic term becomes number, 2nd constant
        expr = par1 * x + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
        # test single parameter with both quadratic terms become numbers
        expr = par1 * par1 + 2
        expected = 0.5 * (1 + 2 * 2 + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == Float64(expected)
        # test single parameter with 2nd quadratic term becomes number
        expr = inf1 * par1 + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (pts[1] + 2 * pts[2] + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
        # test single parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * par1 + 2
        expected = 0.5 * (x + 2 * x + 4)
        @test InfiniteOpt.expand_measure(expr, data1, m) == expected
    end
    # test expand_measure (QuadExpr)
    @testset "QuadExpr (Multi DiscreteMeasureData)" begin
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
               idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
               rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
               rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 2 * pts[1] * rv1 - 2 * inf1 * inf2 + 2 * pts[2] * rv2 + 2x + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 3 * inf1 * pts[1] + 3 * inf1 * pts[2] + 2par1 - 2
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 2 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt.expand_measure(expr, data4, m) == expected
        # test array parameter with first quadratic term becomes number
        expr = pars1[1] * inf4 + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = pts[1] + 2 * pts[2] + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array parameter with first quadratic term becomes number, 2nd constant
        expr = pars1[1] * x + 2
        expected = x + 2 * x + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array parameter with both quadratic terms become numbers
        expr = pars1[1] * pars1[2] + 2
        expected = 1 + 2 * 2 + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == Float64(expected)
        # test array parameter with 2nd quadratic term becomes number
        expr = inf4 * pars1[1] + 2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = pts[1] + 2 * pts[2] + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
        # test array parameter with 2nd quadratic term becomes number, 1st constant
        expr = x * pars1[1] + 2
        expected = x + 2 * x + 4
        @test InfiniteOpt.expand_measure(expr, data3, m) == expected
    end
    # test expand_measure (measure)
    @testset "Measure" begin
        # test single param measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt.expand_measure(meas1, data1, m) == expected
        # test another single param
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 0.5 * (6 * pts[1] * x + 6 * pts[2] * x - 3 * (1 + 2))
        @test InfiniteOpt.expand_measure(meas2, data1, m) == expected
        # test array param measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt.expand_measure(meas1, data3, m) == expected
    end
    # test expand_measure (FunctionalDiscreteMeasureData)
    @testset "FunctionalDiscreteMeasureData" begin
        coef(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(par1, coef, 2, All, is_expect = true)
        @test InfiniteOpt.expand_measure(x, data, m) == 2x
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
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2, num_supports = 2)
    @infinite_variable(m, inf3(par2))
    @hold_variable(m, x)
    # test with 1D data
    @testset "1D Discrete Data" begin
        # test with bounds
        data = DiscreteMeasureData(par1, [1], [1], lower_bound = 1,
                                   upper_bound = 1.5)
        @test InfiniteOpt.analytic_expansion(2x + inf3, data, m) == 0.5 * (2x + inf3)
        # test as expectation
        coef(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(par2, coef, 0, WeightedSample,
                                             is_expect = true)
        @test InfiniteOpt.analytic_expansion(inf3, data, m) == inf3
        # test other
        data = DiscreteMeasureData(par1, [1, 0.5], [1, 1.5])
        @test InfiniteOpt.analytic_expansion(x, data, m) == 1.5 * x
    end
    # test with multi dimensional data
    @testset "Multi Discrete Data" begin
        # test with bounds
        coef(a) = ones(length(a))
        data = FunctionalDiscreteMeasureData(pars1, coef, 0, All,
                    lower_bounds = [1, 1], upper_bounds = [1.5, 2])
        @test InfiniteOpt.analytic_expansion(x, data, m) == 0.5x
        # test as expectation
        data = DiscreteMeasureData(pars1, [1], [[1, 1]], is_expect = true)
        @test InfiniteOpt.analytic_expansion(inf3, data, m) == inf3
        # test other
        coef(a) = ones(size(a, 2))
        data = FunctionalDiscreteMeasureData(pars1, coef, 2, All)
        @test InfiniteOpt.analytic_expansion(x, data, m) == 2x
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
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_variable(m, inf1(par1) >= 1)
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @hold_variable(m, x)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2], is_expect = true)
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand_measures (measure)
    @testset "MeasureRef" begin
        # test the first measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test InfiniteOpt.expand_measures(meas1, m) == expected
        # test the second measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test InfiniteOpt.expand_measures(meas2, m) == expected
        # test analytic expansion
        meas3 = measure(x, data1)
        @test InfiniteOpt.expand_measures(meas3, m) == x
    end
    # test expand_measures (GeneralVariableRef)
    @testset "GeneralVariableRef" begin
        @test InfiniteOpt.expand_measures(x, m) == x
        @test InfiniteOpt.expand_measures(par1, m) == par1
        @test InfiniteOpt.expand_measures(inf1, m) == inf1
    end
    # test expand_measures (Finite Expression)
    @testset "Expression (Finite)" begin
        expr = 2x + 2
        @test InfiniteOpt.expand_measures(expr, m) == expr
        expr = x^2 + 2x + 2
        @test InfiniteOpt.expand_measures(expr, m) == expr
    end
    # test expand_measures (GenericAffExpr)
    @testset "AffExpr (General)" begin
        # test returning AffExpr
        expr = 2inf4 + x + 2meas1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 2inf4 + 7x + pts[1] + pts[2] - rv1 - rv2 + 2inf3 - 4
        @test InfiniteOpt.expand_measures(expr, m) == expected
        # test returning QuadExpr
        expr = x + 3par1 + meas2
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = x + 3par1 + 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test InfiniteOpt.expand_measures(expr, m) == expected
    end
    # test expand_measures (GenericQuadExpr)
    @testset "QuadExpr (General)" begin
        # prepare first measure expansion
        expr = 2 * meas1 * x - x + 1
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        m1 = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        # test with this thorough expression
        @test InfiniteOpt.expand_measures(expr, m) == 2 * m1 * x - x + 1
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
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_variable(m, inf1(par1) >= 1)
    @infinite_variable(m, inf2(par1, par2) >= 0)
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @hold_variable(m, x >= 0)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2], lower_bound = 1, upper_bound = 2)
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data3 = DiscreteMeasureData(par2, [2, 2], [1, 2])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test expand
    @testset "expand" begin
        # test the first measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        idx = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 1
        rv1 = GeneralVariableRef(m, idx, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx + 1, ReducedVariableIndex)
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand(meas1) == expected
        @test parameter_values(pts[2]) == (2.,)
        @test eval_supports(rv1) == Dict{Int, Float64}(1 => 1)
        # test the second measure
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand(meas2) == expected
        @test parameter_values(pts[2]) == (Float64[2, 2],)
        # test analytic
        meas = measure(x, data1)
        @test is_analytic(meas)
        @test expand(meas) == 1*x
        # test integral with only 1 support 
        delete(m, meas1)
        delete(m, meas)
        delete_supports(par1)
        add_supports(par1, 1)
        @test_throws ErrorException expand(integral(inf1, par1))
    end

    # test expand_all_measures!
    @testset "expand_all_measures!" begin
        # prepare the model for testing
        @hold_variable(m, y, parameter_bounds = (par2 == 1))
        set_objective_function(m, x + measure(inf1 * par1 + 3, data1))
        @constraint(m, c1, inf1 + x >= 42.)
        @constraint(m, c2, 2x - measure(measure(inf1 * x + par1 + inf2, data1),
                                        data3) == 0.)
        @BDconstraint(m, c3(par2 in [1, 1.5]), measure(inf2, data1) + inf1 >= 0.)
        @constraint(m, c4, measure(inf4 + y, data2) <= 3.)
        # prepare comparison expressions
        idx = length(InfiniteOpt._data_dictionary(m, PointVariable)) + 1
        pts = [GeneralVariableRef(m, idx, PointVariableIndex),
               GeneralVariableRef(m, idx + 1, PointVariableIndex)]
        obj = x + 0.5pts[1] + pts[2] + 3
        pts = [GeneralVariableRef(m, idx + 2, PointVariableIndex),
               GeneralVariableRef(m, idx + 3, PointVariableIndex),
               GeneralVariableRef(m, idx + 4, PointVariableIndex),
               GeneralVariableRef(m, idx + 5, PointVariableIndex),
               GeneralVariableRef(m, idx + 6, PointVariableIndex),
               GeneralVariableRef(m, idx + 7, PointVariableIndex)]
        c2_expected = 2x - (2 * pts[1] * x + 2 * pts[2] * x + pts[3] + pts[4] +
                      pts[5] + pts[6])
        idx2 = InfiniteOpt._data_dictionary(m, ReducedVariable).last_index + 3
        rv1 = GeneralVariableRef(m, idx2, ReducedVariableIndex)
        rv2 = GeneralVariableRef(m, idx2 + 1, ReducedVariableIndex)
        c3_expected = 0.5rv1 + 0.5rv2 + inf1
        pts = [GeneralVariableRef(m, idx + 8, PointVariableIndex),
               GeneralVariableRef(m, idx + 9, PointVariableIndex)]
        c4_expected = pts[1] + pts[2] + 2y
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
        @test constraint_object(c1).func == inf1 + x
        @test constraint_object(c2).func == c2_expected
        @test constraint_object(c3).func == c3_expected
        @test constraint_object(c4).func == c4_expected
        @test constraint_object(c2).func isa JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
        @test constraint_object(c3).func isa JuMP.GenericAffExpr{Float64, GeneralVariableRef}
        @test constraint_object(c4).func isa JuMP.GenericAffExpr{Float64, GeneralVariableRef}
        @test constraint_object(c1).set == MOI.GreaterThan(42.)
        @test constraint_object(c2).set == MOI.EqualTo(6.)
        @test constraint_object(c3).set == MOI.GreaterThan(0.)
        @test constraint_object(c4).set == MOI.LessThan(3.)
        @test isa(constraint_object(c1), ScalarConstraint)
        @test isa(constraint_object(c2), ScalarConstraint)
        @test isa(constraint_object(c3), BoundedScalarConstraint)
        @test isa(constraint_object(c4), BoundedScalarConstraint)
        @test constraint_object(c3).bounds.intervals == Dict(par2 => IntervalSet(1, 1.5))
        @test constraint_object(c4).bounds.intervals == Dict(par2 => IntervalSet(1, 1))
    end
end
