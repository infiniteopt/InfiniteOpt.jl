# Test basic methods
@testset "Basics" begin
    # initialize model and point variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])

    ivar = InfiniteVariable(info, VectorTuple(pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    var = PointVariable(info, ivref, VectorTuple{Float64}(0.5, 0.5))
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
        values = [1., 2., 3.]
        @test InfiniteOpt._make_str_value(values) == "[1, 2, 3]"
        # test long array
        values = [1., 2., 3., 4., 5., 6.]
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
        ivref2 = add_variable(m, InfiniteVariable(info, VectorTuple(pref)), "ivar2")
        m.vars[5] = PointVariable(info, ivref2, VectorTuple{Float64}(0.5))
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
    ivar = InfiniteVariable(info, VectorTuple(pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    ivar2 = build_variable(error, info, Infinite, parameter_refs = (pref, prefs))
    ivref2 = add_variable(m, ivar2, "ivar2")
    # _check_tuple_shape
    @testset "_check_tuple_shape" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref, VectorTuple{Float64}(0.5, 0.5)),
                  Nothing)
        # prepare param value tuple
        tuple = VectorTuple{Float64}(0.5, [0.5, 0.5])
        # test normal with array
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref2, tuple), Nothing)
        # test for errors in shape
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                              VectorTuple{Float64}(0.5,))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                        VectorTuple{Float64}(0.5, [0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                          VectorTuple{Float64}(0.5, 0.5))
        tuple = VectorTuple{Float64}(0.5, [0.5, 0.5, 0.5])
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                                   tuple)
    end
    # _check_tuple_values
    @testset "_check_tuple_values" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_values(error, ivref, VectorTuple{Float64}(0.5, 0.5)),
                  Nothing)
        # test normal with array
        tuple = VectorTuple{Float64}(0, [0.5, 1])
        @test isa(InfiniteOpt._check_tuple_values(error, ivref2, tuple), Nothing)
        # test for out of bound errors
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, ivref,
                                                              VectorTuple{Float64}(0, 2))
        tuple = VectorTuple{Float64}(0, [2, 1])
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
                   parameter_values = (0.5, 0.5)).parameter_values == VectorTuple{Float64}(0.5, 0.5)
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Point),
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = VectorTuple{Float64}(0.5, [0, 0])
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
                   parameter_values = (0.5, 0.5)).parameter_values == VectorTuple{Float64}(0.5, 0.5)
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = VectorTuple{Float64}(0.5, [0, 0])
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == tuple
        @test_throws ErrorException build_variable(error, info, Point,
                                            infinite_variable_ref = ivref2,
                                            parameter_values = (0.5, [0, 0, 0]))
    end
    # _update_param_supports
    @testset "_update_param_supports" begin
        # test normal
        @test isa(InfiniteOpt._update_param_supports(ivref, VectorTuple{Float64}(0.5, 1)),
                  Nothing)
        @test supports(pref) == Float64[0.5]
        @test supports(pref2) == Float64[1]
        # prepare array tuple
        tuple = VectorTuple{Float64}(0.5, [0, 1])
        # test normal with array
        @test isa(InfiniteOpt._update_param_supports(ivref2, tuple), Nothing)
        @test supports(pref) == Float64[0.5]
        @test supports(prefs[1]) == Float64[0]
        @test supports(prefs[2]) == Float64[1]
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
        delete!(m.infinite_to_points, JuMP.index(ivref))
    end
    # _check_and_make_variable_ref
    @testset "_check_and_make_variable_ref" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        pref3 = add_parameter(m2, param, "test")
        ivar3 = InfiniteVariable(info, VectorTuple(pref3))
        ivref3 = add_variable(m2, ivar3, "ivar")
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        # test for invalid variable error
        @test_throws ErrorException InfiniteOpt._check_and_make_variable_ref(m, v)
        # test normal
        v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        @test InfiniteOpt._check_and_make_variable_ref(m, v) == PointVariableRef(m, 2)
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
        ivar3 = InfiniteVariable(info, VectorTuple(pref3))
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

# Test the point variable macro
@testset "Point" begin
    # initialize model, parameters, and infinite variables
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= t <= 1)
    @infinite_parameter(m, -1 <= x[1:2] <= 1)
    @infinite_variable(m, 0 <= z(t, x) <= 1, Int)
    @infinite_variable(m, z2[1:2](t) == 3)
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        vref = PointVariableRef(m, 4)
        @test @point_variable(m, infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_integer(vref)
        @test lower_bound(vref) == 0
        # test anon with changes to fixed
        vref = PointVariableRef(m, 5)
        @test @point_variable(m, infinite_variable_ref = z, lower_bound = -5,
                          parameter_values = (0, [0, 0]), binary = true) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test !is_integer(vref)
        @test is_binary(vref)
        @test lower_bound(vref) == -5
        # test regular with alias
        vref = PointVariableRef(m, 6)
        @test @point_variable(m, z(0, [0, 0]), z0, Bin) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
        # test regular with semi anon
        vref = PointVariableRef(m, 7)
        @test @point_variable(m, z(0, [0, 0]), base_name = "z0",
                              binary = true) == vref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
    end
    # test array variable definition
    @testset "Array" begin
        # test anon array with one infvar
        vrefs = [PointVariableRef(m, 8), PointVariableRef(m, 9)]
        @test @point_variable(m, [1:2], infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_integer(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        # test anon array with different inf vars
        vrefs = [PointVariableRef(m, 10), PointVariableRef(m, 11)]
        @test @point_variable(m, [i = 1:2], infinite_variable_ref = z2[i],
                              parameter_values = 0) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test fix_value(vrefs[2]) == 3
        @test name(vrefs[1]) == "z2[1](0)"
        # test array with same infvar
        vrefs = [PointVariableRef(m, 12), PointVariableRef(m, 13)]
        @test @point_variable(m, z(0, [0, 0]), a[1:2], Bin) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_binary(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        @test name(vrefs[1]) == "a[1]"
        # test test array with differnt infvars
        vrefs = [PointVariableRef(m, 14), PointVariableRef(m, 15)]
        @test @point_variable(m, z2[i](0), b[i = 1:2] >= -5) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "b[1]"
        # test semi anon array
        vrefs = [PointVariableRef(m, 16), PointVariableRef(m, 17)]
        @test @point_variable(m, z2[i](0), [i = 1:2], lower_bound = -5) == vrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "z2[1](0)"
    end
    # test errors
    @testset "Errors" begin
        # test model assertion errors
        m2 = Model()
        @test_throws AssertionError @point_variable(m2, infinite_variable_ref = z,
                                                 parameter_values = (0, [0, 0]))
        @test_throws AssertionError @point_variable(m2, z(0, [0, 0]), bob, Bin)
        @test_throws AssertionError @point_variable(m2, [1:2],
                      infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        # test double specification
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), bob,
                                                      infinite_variable_ref = z)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), bob,
                                                 parameter_values = (0, [0, 0]))
        # test the adding expressions to infvar
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) >= 0,
                                                          bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) == 0,
                                                          bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]) in 0,
                                                          bob, Bin)
        # test other syntaxes
        @test_macro_throws ErrorException @point_variable(m,
                                               0 <= z(0, [0, 0]) <= 1, bob, Bin)
        @test_macro_throws ErrorException @point_variable(m, [1:2], Bin,
                      infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        @test_macro_throws ErrorException @point_variable(m, bob,
                     infinite_variable_ref = z, parameter_values = (0, [0, 0]))
        # test redefinition catch
        @test_macro_throws ErrorException @point_variable(m, z(0, [0, 0]), z0,
                                                          Bin)
    end
end
