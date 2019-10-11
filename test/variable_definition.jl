# Test name methods
@testset "Infinite Variable Name" begin
    # initialize model and infinite variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test1")
    pref2 = add_parameter(m, param, "test2")
    var = InfiniteVariable(info, (pref, pref2))
    m.vars[1] = var
    m.var_to_name[1] = "var"
    vref = InfiniteVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == (pref, pref2)
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # make extra infinite variable
        vref2 = InfiniteVariableRef(m, 2)
        m.vars[2] = var
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new(test1, test2)"
        # test default
        @test isa(set_name(vref2, ""), Nothing)
        @test name(vref2) == "noname(test1, test2)"
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, 1) == vref
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new(test1, test2)") == vref
        @test isa(variable_by_name(m, "test(test1, test2)"), Nothing)
        # prepare variable with same name
        m.vars[2] = var
        m.var_to_name[2] = "new(test1, test2)"
        m.name_to_var = nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "new(test1, test2)")
    end
    # _root_name
    @testset "_root_name" begin
        @test InfiniteOpt._root_name(vref) == "new"
    end
end

# Test variable definition methods
@testset "Infinite Variable Definition" begin
    # initialize model and infinite variable info
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    pref2 = add_parameter(m, param, "θ")
    prefs = @infinite_parameter(m, x[1:2], set = IntervalSet(0, 1))
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    info2 = VariableInfo(true, 0, true, 0, true, 0, true, 0, true, false)
    info3 = VariableInfo(true, 0, true, 0, true, 0, true, 0, false, true)
    # _check_parameter_tuple
    @testset "_check_parameter_tuple" begin
        @test isa(InfiniteOpt._check_parameter_tuple(error, (pref, prefs)),
                  Nothing)
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error,
                                                               (pref, prefs, 2))
    end
    # _make_formatted_tuple
    @testset "_make_formatted_tuple" begin
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs)), Tuple)
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs))[2],
                  JuMP.Containers.SparseAxisArray)
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs))[1],
                  ParameterRef)
    end
    # _check_tuple_groups
    @testset "_check_tuple_groups" begin
        # prepare param tuple
        tuple = InfiniteOpt._make_formatted_tuple((pref, prefs))
        # test normal
        @test isa(InfiniteOpt._check_tuple_groups(error, tuple), Nothing)
        # prepare bad param tuple
        tuple = InfiniteOpt._make_formatted_tuple(([pref; pref2], prefs))
        # test for errors
        @test_throws ErrorException InfiniteOpt._check_tuple_groups(error,
                                                                    tuple)
        @test_throws ErrorException InfiniteOpt._check_tuple_groups(error,
                                                                   (pref, pref))
    end
    # _make_variable
    @testset "_make_variable" begin
        # test for each error message
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Infinite),
                                                   bob = 42)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, :bad)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Infinite))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Infinite),
                                                   parameter_refs = (pref, 2))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Infinite),
                                                  parameter_refs = (pref, pref))
        # defined expected output
        expected = InfiniteVariable(info, (pref,))
        # test for expected output
        @test InfiniteOpt._make_variable(error, info, Val(Infinite),
                             parameter_refs = pref).info == expected.info
        @test InfiniteOpt._make_variable(error, info, Val(Infinite),
                parameter_refs = pref).parameter_refs == expected.parameter_refs
        # test various types of param tuples
        @test InfiniteOpt._make_variable(error, info, Val(Infinite),
                 parameter_refs = (pref, pref2)).parameter_refs == (pref, pref2)
        tuple = InfiniteOpt._make_formatted_tuple((pref, prefs))
        @test InfiniteOpt._make_variable(error, info, Val(Infinite),
                         parameter_refs = (pref, prefs)).parameter_refs == tuple
        tuple = InfiniteOpt._make_formatted_tuple((prefs,))
        @test InfiniteOpt._make_variable(error, info, Val(Infinite),
                             parameter_refs = prefs).parameter_refs == tuple
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for each error message
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   bob = 42)
        @test_throws ErrorException build_variable(error, info, :bad)
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_refs = pref)
        @test_throws ErrorException build_variable(error, info, Infinite)
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   parameter_refs = (pref, 2))
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                  parameter_refs = (pref, pref),
                                                  error = error)
        # defined expected output
        expected = InfiniteVariable(info, (pref,))
        # test for expected output
        @test build_variable(error, info, Infinite,
                             parameter_refs = pref).info == expected.info
        @test build_variable(error, info, Infinite,
                parameter_refs = pref).parameter_refs == expected.parameter_refs
        # test various types of param tuples
        @test build_variable(error, info, Infinite,
                 parameter_refs = (pref, pref2)).parameter_refs == (pref, pref2)
        tuple = InfiniteOpt._make_formatted_tuple((pref, prefs))
        @test build_variable(error, info, Infinite,
                         parameter_refs = (pref, prefs)).parameter_refs == tuple
        tuple = InfiniteOpt._make_formatted_tuple((prefs,))
        @test build_variable(error, info, Infinite,
                             parameter_refs = prefs).parameter_refs == tuple
    end
    # _update_param_var_mapping
    @testset "_update_param_var_mapping" begin
        # initialize secondary model and infinite variable
        m2 = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref3 = add_parameter(m2, param, "test")
        prefs2 = @infinite_parameter(m2, x[1:2], set = IntervalSet(0, 1))
        ivref = InfiniteVariableRef(m2, 1)
        ivref2 = InfiniteVariableRef(m2, 2)
        # prepare tuple
        tuple = (pref3, prefs2)
        tuple = InfiniteOpt._make_formatted_tuple(tuple)
        # test normal
        @test isa(InfiniteOpt._update_param_var_mapping(ivref, tuple), Nothing)
        @test m2.param_to_vars[1] == [1]
        @test m2.param_to_vars[2] == [1]
        @test m2.param_to_vars[3] == [1]
        @test isa(InfiniteOpt._update_param_var_mapping(ivref2, tuple), Nothing)
        @test m2.param_to_vars[1] == [1, 2]
        @test m2.param_to_vars[2] == [1, 2]
        @test m2.param_to_vars[3] == [1, 2]
    end
    # _check_parameters_valid
    @testset "_check_parameters_valid" begin
        # prepare param tuple
        tuple = (pref, prefs, copy(pref2, InfiniteModel()))
        tuple = InfiniteOpt._make_formatted_tuple(tuple)
        # test that catches error
        @test_throws ErrorException InfiniteOpt._check_parameters_valid(m, tuple)
        # test normal
        @test isa(InfiniteOpt._check_parameters_valid(m, (pref, pref2)), Nothing)
    end
    # _check_make_variable_ref
    @testset "_check_make_variable_ref" begin
        # prepare secondary model and parameter and variable
        m2 = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref3 = add_parameter(m2, param, "test")
        v = build_variable(error, info, Infinite,
                           parameter_refs = pref3)
        # test for error of invalid variable
        @test_throws ErrorException InfiniteOpt._check_make_variable_ref(m, v)
        # prepare normal variable
        v = build_variable(error, info, Infinite, parameter_refs = pref)
        # test normal
        @test InfiniteOpt._check_make_variable_ref(m, v) == InfiniteVariableRef(m, 0)
        @test m.param_to_vars[1] == [0]
        delete!(m.param_to_vars, 1)
        # test with other variable object
        @test_throws ErrorException InfiniteOpt._check_make_variable_ref(m, :bad)
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        # prepare secondary model and parameter and variable
        m2 = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref3 = add_parameter(m2, param, "test")
        v = build_variable(error, info, Infinite,
                           parameter_refs = pref3)
        # test for error of invalid variable
        @test_throws ErrorException add_variable(m, v)
        # prepare normal variable
        v = build_variable(error, info, Infinite, parameter_refs = pref)
        # test normal
        @test add_variable(m, v, "name") == InfiniteVariableRef(m, 2)
        @test haskey(m.vars, 2)
        @test m.param_to_vars[1] == [2]
        @test m.var_to_name[2] == "name(test)"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Infinite, parameter_refs = pref)
        # test info addition functions
        vref = InfiniteVariableRef(m, 3)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        # lower bound
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == 1
        @test isa(m.constrs[1], ScalarConstraint{InfiniteVariableRef,
                                                 MOI.GreaterThan{Float64}})
        @test m.constr_in_var_info[1]
        # upper bound
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == 2
        @test isa(m.constrs[2], ScalarConstraint{InfiniteVariableRef,
                                                 MOI.LessThan{Float64}})
        @test m.constr_in_var_info[2]
        # fix
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == 3
        @test isa(m.constrs[3], ScalarConstraint{InfiniteVariableRef,
                                                 MOI.EqualTo{Float64}})
        @test m.constr_in_var_info[3]
        # binary
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == 4
        @test isa(m.constrs[4], ScalarConstraint{InfiniteVariableRef,
                                                 MOI.ZeroOne})
        @test m.constr_in_var_info[4]
        @test m.var_to_constrs[3] == [1, 2, 3, 4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Infinite, parameter_refs = pref)
        # test integer addition functions
        vref = InfiniteVariableRef(m, 4)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == 8
        @test isa(m.constrs[8], ScalarConstraint{InfiniteVariableRef,
                                                 MOI.Integer})
        @test m.constr_in_var_info[8]
        @test m.var_to_constrs[4] == [5, 6, 7, 8]
    end
end

# Test name methods
@testset "Point Variable Name" begin
    # initialize model and point variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test1")
    pref2 = add_parameter(m, param, "test2")
    ivar = InfiniteVariable(info, (pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    var = PointVariable(info, ivref, (0.5, 0.5))
    m.vars[2] = var
    m.var_to_name[2] = "var"
    vref = PointVariableRef(m, 2)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(vref) == ivref
    end
    # _make_str_value (Number)
    @testset "_make_str_value (Number)" begin
        @test InfiniteOpt._make_str_value(1.0) == "1"
    end
    # _make_str_value (Array)
    @testset "_make_str_value (Array)" begin
        # test short array
        values = convert(JuMPC.SparseAxisArray, [1., 2., 3.])
        @test InfiniteOpt._make_str_value(values) == "[1, 2, 3]"
        # test long array
        values = convert(JuMPC.SparseAxisArray, [1., 2., 3., 4., 5., 6.])
        @test InfiniteOpt._make_str_value(values) == "[1, ..., 6]"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # prepare a secondary point variable
        vref2 = PointVariableRef(m, 3)
        m.vars[3] = var
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        # test default
        @test isa(set_name(vref2, ""), Nothing)
        @test name(vref2) == "ivar(0.5, 0.5)"
        # test other default
        m.next_var_index = 3
        ivref2 = add_variable(m, InfiniteVariable(info, (pref,)), "ivar2")
        m.vars[5] = PointVariable(info, ivref2, (0.5,))
        m.var_to_name[5] = "var42"
        vref3 = PointVariableRef(m, 5)
        @test isa(set_name(vref3, ""), Nothing)
        @test name(vref3) == "ivar2(0.5)"
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, 2) == vref
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new") == vref
        @test isa(variable_by_name(m, "test"), Nothing)
        # make variable with duplicate name
        m.vars[3] = var
        m.var_to_name[3] = "new"
        m.name_to_var = nothing
        # test for multiple name error
        @test_throws ErrorException variable_by_name(m, "new")
    end
end

# Test variable definition methods
@testset "Point Variable Definition" begin
    # initialize model and infinite variables
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    pref2 = add_parameter(m, param, "θ")
    prefs = @infinite_parameter(m, x[1:2], set = IntervalSet(0, 1))
    info = VariableInfo(false, 0., false, 0., false, 0., false, NaN, false, false)
    info2 = VariableInfo(true, 0, true, 0, true, 0, true, 0, true, false)
    info3 = VariableInfo(true, 0, true, 0, true, 0, true, 0, false, true)
    ivar = InfiniteVariable(info, (pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    ivar2 = build_variable(error, info, Infinite, parameter_refs = (pref, prefs))
    ivref2 = add_variable(m, ivar2, "ivar2")
    # _check_tuple_shape
    @testset "_check_tuple_shape" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref, (0.5, 0.5)),
                  Nothing)
        # prepare param value tuple
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0.5, 0.5]))
        # test normal with array
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref2, tuple), Nothing)
        # test for errors in shape
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                                   (0.5,))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                                   (0.5, [0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                                   (0.5, 0.5))
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0.5, 0.5, 0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                                   tuple)
    end
    # _check_tuple_values
    @testset "_check_tuple_values" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_values(error, ivref, (0.5, 0.5)),
                  Nothing)
        # prepare array test
        tuple = InfiniteOpt._make_formatted_tuple((0, [0.5, 1]))
        # test normal with array
        @test isa(InfiniteOpt._check_tuple_values(error, ivref2, tuple), Nothing)
        # test for out of bound errors
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, ivref,
                                                                    (0, 2))
        tuple = InfiniteOpt._make_formatted_tuple((0, [2, 1]))
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, ivref2,
                                                                    tuple)
    end
    # _update_point_info
    @testset "_update_point_info" begin
        # prepare info for test
        new_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                false)
        InfiniteOpt._update_variable_info(ivref, new_info)
        # test with current info
        @test InfiniteOpt._update_point_info(info, ivref) == new_info
        # prepare info for test
        new_info = VariableInfo(false, 0., false, 0., true, 0., true, 0., false,
                                true)
        InfiniteOpt._update_variable_info(ivref, new_info)
        # test with current info
        @test InfiniteOpt._update_point_info(info, ivref) == new_info
        # prepare info for test
        curr_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                 false)
        # test with current info
        @test InfiniteOpt._update_point_info(curr_info, ivref) == curr_info
        # undo info changes
        InfiniteOpt._update_variable_info(ivref, info)
    end
    # test _make_variable
    @testset "_make_variable" begin
        # test for all errors
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                              Val(Point), parameter_refs = pref)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                                               Val(Point))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                      Val(Point), infinite_variable_ref = ivref)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                               Val(Point), parameter_values = 3)
        # test a variety of builds
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).infinite_variable_ref == ivref
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).parameter_values == (0.5, 0.5)
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Point),
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 0]))
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == tuple
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Point),
                                            infinite_variable_ref = ivref2,
                                            parameter_values = (0.5, [0, 0, 0]))
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for all errors
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                  infinite_variable_ref = ivref)
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   parameter_values = 3)
        @test_throws ErrorException build_variable(error, info, Point)
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref)
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_values = 3)
        # test a variety of builds
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).infinite_variable_ref == ivref
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).parameter_values == (0.5, 0.5)
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 0]))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == tuple
        @test_throws ErrorException build_variable(error, info, Point,
                                            infinite_variable_ref = ivref2,
                                            parameter_values = (0.5, [0, 0, 0]))
    end
    # _update_param_supports
    @testset "_update_param_supports" begin
        # test normal
        @test isa(InfiniteOpt._update_param_supports(ivref, (0.5, 1)), Nothing)
        @test supports(pref) == [0.5]
        @test supports(pref2) == [1]
        # prepare array tuple
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 1]))
        # test normal with array
        @test isa(InfiniteOpt._update_param_supports(ivref2, tuple), Nothing)
        @test supports(pref) == [0.5]
        @test supports(prefs[1]) == [0]
        @test supports(prefs[2]) == [1]
    end
    # _update_infinite_point_mapping
    @testset "_update_infinite_point_mapping" begin
        # test first addition
        pvref = PointVariableRef(m, 12)
        @test isa(InfiniteOpt._update_infinite_point_mapping(pvref, ivref),
                  Nothing)
        @test m.infinite_to_points[JuMP.index(ivref)] == [12]
        # test second addition
        pvref = PointVariableRef(m, 42)
        @test isa(InfiniteOpt._update_infinite_point_mapping(pvref, ivref),
                  Nothing)
        @test m.infinite_to_points[JuMP.index(ivref)] == [12, 42]
        # undo changes
        delete!( m.infinite_to_points, JuMP.index(ivref))
    end
    # _check_make_variable_ref
    @testset "_check_make_variable_ref" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        pref3 = add_parameter(m2, param, "test")
        ivar3 = InfiniteVariable(info, (pref3,))
        ivref3 = add_variable(m2, ivar3, "ivar")
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        # test for invalid variable error
        @test_throws ErrorException InfiniteOpt._check_make_variable_ref(m, v)
        # test normal
        v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        @test InfiniteOpt._check_make_variable_ref(m, v) == PointVariableRef(m, 2)
        @test supports(pref) == [0.5, 0]
        @test supports(pref2) == [1]
        @test m.infinite_to_points[JuMP.index(ivref)] == [2]
        delete!(m.infinite_to_points, JuMP.index(ivref))
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        pref3 = add_parameter(m2, param, "test")
        ivar3 = InfiniteVariable(info, (pref3,))
        ivref3 = add_variable(m2, ivar3, "ivar")
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        # test for invalid variable error
        @test_throws ErrorException add_variable(m, v)
        # test normal
        v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        @test add_variable(m, v, "name") == PointVariableRef(m, 4)
        @test haskey(m.vars, 4)
        @test supports(pref) == [0.5, 0]
        @test supports(pref2) == [1]
        @test m.var_to_name[4] == "name"
        @test m.infinite_to_points[JuMP.index(ivref)] == [4]
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        # test info addition functions
        vref = PointVariableRef(m, 5)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        # lower bound
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == 1
        @test isa(m.constrs[1], ScalarConstraint{PointVariableRef,
                                                 MOI.GreaterThan{Float64}})
        @test m.constr_in_var_info[1]
        # upper bound
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == 2
        @test isa(m.constrs[2], ScalarConstraint{PointVariableRef,
                                                 MOI.LessThan{Float64}})
        @test m.constr_in_var_info[2]
        # fix
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == 3
        @test isa(m.constrs[3], ScalarConstraint{PointVariableRef,
                                                 MOI.EqualTo{Float64}})
        @test m.constr_in_var_info[3]
        # binary
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == 4
        @test isa(m.constrs[4], ScalarConstraint{PointVariableRef,
                                                 MOI.ZeroOne})
        @test m.constr_in_var_info[4]
        @test m.var_to_constrs[5] == [1, 2, 3, 4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        # test integer addition functions
        vref = PointVariableRef(m, 6)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == 8
        @test isa(m.constrs[8], ScalarConstraint{PointVariableRef,
                                                 MOI.Integer})
        @test m.constr_in_var_info[8]
        @test m.var_to_constrs[6] == [5, 6, 7, 8]
    end
end

# Test name methods
@testset "Hold Variable Name" begin
    # initialize model and variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    dict = Dict{ParameterRef, IntervalSet}()
    var = HoldVariable(info, dict)
    m.vars[1] = var
    m.var_to_name[1] = "test"
    vref = HoldVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, 1) == vref
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new") == vref
        @test isa(variable_by_name(m, "test2"), Nothing)
        # prepare variable with duplicate name
        m.vars[2] = var
        m.var_to_name[2] = "new"
        m.name_to_var = nothing
        # test for duplciate name error
        @test_throws ErrorException variable_by_name(m, "new")
    end
end

# Test variable definition methods
@testset "Hold Variable Definition" begin
    # initialize model and info
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    info2 = VariableInfo(true, 0, true, 0, true, 0, true, 0, true, false)
    info3 = VariableInfo(true, 0, true, 0, true, 0, true, 0, false, true)
    dict = Dict{ParameterRef, IntervalSet}()
    @infinite_parameter(m, 0 <= par <= 10)
    @infinite_parameter(m, 0 <= pars[1:2] <= 10)
    # test _check_bounds
    @testset "_check_bounds" begin
        # test normal
        @test isa(InfiniteOpt._check_bounds(Dict(par => IntervalSet(0, 1))),
                                                                        Nothing)
        # test errors
        @test_throws ErrorException InfiniteOpt._check_bounds(
                                                Dict(par => IntervalSet(-1, 1)))
        @test_throws ErrorException InfiniteOpt._check_bounds(
                                                Dict(par => IntervalSet(0, 11)))
    end
    # test _expand_parameter_dict(Dict{ParameterRef,IntervalSet}))
    @testset "_expand_parameter_dict (acceptable Form)" begin
        d = Dict(par => IntervalSet(0, 1))
        @test InfiniteOpt._expand_parameter_dict(d) == d
    end
    # test _expand_parameter_dict(Dict{Any,IntervalSet}))
    @testset "_expand_parameter_dict (Array Form)" begin
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(InfiniteOpt._expand_parameter_dict(d),
                  Dict{ParameterRef, IntervalSet})
    end
    # test _expand_parameter_dict(Dict))
    @testset "_expand_parameter_dict (Fallback)" begin
        d = Dict(pars => 1, par => 2)
        @test_throws ErrorException InfiniteOpt._expand_parameter_dict(d)
    end
    # _make_variable
    @testset "_make_variable" begin
        # test normal
        expected = HoldVariable(info, dict)
        @test InfiniteOpt._make_variable(error, info, Val(Hold)).info == expected.info
        # test errors
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                                Val(Hold), parameter_values = 3)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test normal
        expected = HoldVariable(info, dict)
        @test build_variable(error, info, Hold).info == expected.info
        # test errors
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_bounds = dict)
    end
    # _validate_bounds
    @testset "_validate_bounds" begin
        # test normal
        @test isa(InfiniteOpt._validate_bounds(m, Dict(par => IntervalSet(0, 1))),
                  Nothing)
        # test error
        par2 = ParameterRef(InfiniteModel(), 1)
        @test_throws ErrorException InfiniteOpt._validate_bounds(m,
                                                Dict(par2 => IntervalSet(0, 1)))
        # test support addition
        @test isa(InfiniteOpt._validate_bounds(m, Dict(par => IntervalSet(0, 0))),
                  Nothing)
        @test supports(par) == [0]
    end
    # _check_make_variable_ref
    @testset "_check_make_variable_ref" begin
        # test normal
        v = build_variable(error, info, Hold)
        @test InfiniteOpt._check_make_variable_ref(m, v) == HoldVariableRef(m, 0)
        # test bad bounds
        @infinite_parameter(InfiniteModel(), par2 in [0, 2])
        v = build_variable(error, info, Hold,
                           parameter_bounds = Dict(par2 => IntervalSet(0, 1)))
        @test_throws ErrorException InfiniteOpt._check_make_variable_ref(m, v)
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        v = build_variable(error, info, Hold)
        @test add_variable(m, v, "name") == HoldVariableRef(m, 1)
        @test haskey(m.vars, 1)
        @test m.var_to_name[1] == "name"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Hold)
        # test info addition functions
        vref = HoldVariableRef(m, 2)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        # lower bound
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == 1
        @test isa(m.constrs[1], ScalarConstraint{HoldVariableRef,
                                                 MOI.GreaterThan{Float64}})
        @test m.constr_in_var_info[1]
        # upper bound
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == 2
        @test isa(m.constrs[2], ScalarConstraint{HoldVariableRef,
                                                 MOI.LessThan{Float64}})
        @test m.constr_in_var_info[2]
        # fix
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == 3
        @test isa(m.constrs[3], ScalarConstraint{HoldVariableRef,
                                                 MOI.EqualTo{Float64}})
        @test m.constr_in_var_info[3]
        # binary
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == 4
        @test isa(m.constrs[4], ScalarConstraint{HoldVariableRef,
                                                 MOI.ZeroOne})
        @test m.constr_in_var_info[4]
        @test m.var_to_constrs[2] == [1, 2, 3, 4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Hold)
        # test integer addition functions
        vref = HoldVariableRef(m, 3)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == 8
        @test isa(m.constrs[8], ScalarConstraint{HoldVariableRef,
                                                 MOI.Integer})
        @test m.constr_in_var_info[8]
        @test m.var_to_constrs[3] == [5, 6, 7, 8]
    end
end
