# Test basic methods
@testset "Basics" begin
    # initialize model and point variable
    m = InfiniteModel()
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])
    @infinite_variable(m, ivref(a, b, c))
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, 0., true, false)
    var = PointVariable(info, ivref, [0., 0., 1., 1., 0.])
    object = VariableData(var, "var")
    idx = PointVariableIndex(1)
    vref = PointVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, PointVariableIndex)
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(vref) === m
        @test owner_model(gvref) === m
    end
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(vref) == idx
        @test index(gvref) == idx
    end
    # dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == vref
        @test dispatch_variable_ref(gvref) == vref
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(m, PointVariable) === m.point_vars
        @test InfiniteOpt._data_dictionary(vref) === m.point_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.point_vars
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, vref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(vref) === object
        @test InfiniteOpt._data_object(gvref) === object
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) === var
        @test InfiniteOpt._core_variable_object(gvref) === var
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
    end
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(vref) == info
    end
    # _update_variable_info
    @testset "_update_variable_info" begin
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(vref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(vref) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == ConstraintIndex[]
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
        @test name(gvref) == "var"
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(vref) == ivref
        @test infinite_variable_ref(gvref) == ivref
    end
    # raw_parameter_values
    @testset "raw_parameter_values" begin
        @test raw_parameter_values(vref) == Float64[0., 0., 1., 1., 0.]
        @test raw_parameter_values(gvref) == Float64[0., 0., 1., 1., 0.]
    end
    # parameter_values
    @testset "parameter_values" begin
        @test parameter_values(vref) == (Float64(0), Float64[0, 1], Float64[1, 0])
        @test parameter_values(gvref) == (Float64(0), Float64[0, 1], Float64[1, 0])
    end
    # _update_variable_param_values
    @testset "_update_variable_param_values" begin
        @test isa(InfiniteOpt._update_variable_param_values(vref, ones(5)),
                  Nothing)
        @test raw_parameter_values(vref) == ones(5)
        @test isa(InfiniteOpt._update_variable_param_values(vref, [0., 0., 1., 1., 0.]),
                  Nothing)
    end
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == ()
        @test parameter_refs(gvref) == ()
    end
    # _make_str_value (Number)
    @testset "_make_str_value (Number)" begin
        @test InfiniteOpt._make_str_value(1.0) == "1"
    end
    # _make_str_value (Array)
    @testset "_make_str_value (Array)" begin
        # test single
        @test InfiniteOpt._make_str_value([1.1]) == "1.1"
        # test short array
        values = [1., 2., 3.]
        @test InfiniteOpt._make_str_value(values) == "[1, 2, 3]"
        # test long array
        values = [1., 2., 3., 4., 5., 6.]
        @test InfiniteOpt._make_str_value(values) == "[1, ..., 6]"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        # test default
        @test isa(set_name(gvref, ""), Nothing)
        @test name(vref) == "ivref(0, [0, 1], [1, 0])"
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, idx) == gvref
    end
    # _var_name_dict
    @testset "_var_name_dict" begin
        @test InfiniteOpt._var_name_dict(m) isa Nothing
    end
    # _update_var_name_dict
    @testset "_update_var_name_dict" begin
        m.name_to_var = Dict{String, ObjectIndex}()
        dict = InfiniteOpt._data_dictionary(vref)
        @test InfiniteOpt._update_var_name_dict(m, dict) isa Nothing
        @test InfiniteOpt._var_name_dict(m) == Dict(name(vref) => idx)
        m.name_to_var = nothing
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "ivref(0, [0, 1], [1, 0])") == gvref
        @test variable_by_name(m, "test(0, [0, 1], [1, 0])") isa Nothing
        # prepare variable with same name
        idx2 = PointVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = PointVariableRef(m, idx2)
        @test set_name(vref2, "") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "ivref(0, [0, 1], [1, 0])")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(vref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(vref)) == 1
        @test !is_valid(m, vref)
    end
end

# Test variable definition methods
@testset "Definition" begin
    # initialize model and infinite variables
    m = InfiniteModel()
    @independent_parameter(m, pref in [0, 1])
    @independent_parameter(m, pref2 in [0, 1])
    @dependent_parameters(m, prefs[1:2] in [0, 1])
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, num, false, true)
    @infinite_variable(m, ivref(pref, pref2))
    @infinite_variable(m, ivref2(pref, prefs))
    divref = dispatch_variable_ref(ivref)
    divref2 = dispatch_variable_ref(ivref2)
    # _check_tuple_shape
    @testset "_check_tuple_shape" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_shape(error, divref, IC.VectorTuple{Float64}(0.5, 0.5)),
                  Nothing)
        # prepare param value tuple
        tuple = IC.VectorTuple{Float64}(0.5, [0.5, 0.5])
        # test normal with array
        @test isa(InfiniteOpt._check_tuple_shape(error, divref2, tuple), Nothing)
        # test for errors in shape
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref,
                                                              IC.VectorTuple{Float64}(0.5,))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref,
                                                        IC.VectorTuple{Float64}(0.5, [0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref2,
                                                          IC.VectorTuple{Float64}(0.5, 0.5))
        tuple = IC.VectorTuple{Float64}(0.5, [0.5, 0.5, 0.5])
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref2,
                                                                   tuple)
    end
    # _check_element_support
    @testset "_check_element_support (IndependentParameterRef)" begin
        # test normal
        ps = [dispatch_variable_ref(pref)]
        vals = Float64[0, 2]
        @test InfiniteOpt._check_element_support(error, ps, vals, 1) == 2
        # test error
        ps = [dispatch_variable_ref(pref), dispatch_variable_ref(pref2)]
        @test_throws ErrorException InfiniteOpt._check_element_support(error, ps, vals, 1)
    end
    # _check_element_support
    @testset "_check_element_support (DependentParameterRef)" begin
        # test normal
        ps = dispatch_variable_ref.(prefs)
        vals = Float64[0, 1, 0, 1]
        @test InfiniteOpt._check_element_support(error, ps, vals, 3) == 5
        # test error
        vals = Float64[0, 0, 2, 1]
        @test_throws ErrorException InfiniteOpt._check_element_support(error, ps, vals, 3)
    end
    # _check_tuple_values
    @testset "_check_tuple_values" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_values(error, divref, Float64[0.5, 0.5]),
                  Nothing)
        # test normal with array
        vals = Float64[0, 0.5, 1]
        @test isa(InfiniteOpt._check_tuple_values(error, divref2, vals), Nothing)
        # test for out of bound errors
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, divref,
                                                                    Float64[0, 2])
        vals = Float64[0, 2, 1]
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, divref2,
                                                                    vals)
    end
    # _update_point_info
    @testset "_update_point_info" begin
        # prepare info for test
        new_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                false)
        InfiniteOpt._update_variable_info(divref, new_info)
        # test with current info
        @test InfiniteOpt._update_point_info(info, divref) == new_info
        # prepare info for test
        new_info = VariableInfo(false, 0., false, 0., true, 0., true, 0., false,
                                true)
        InfiniteOpt._update_variable_info(divref, new_info)
        # test with current info
        @test InfiniteOpt._update_point_info(info, divref) == new_info
        # prepare info for test
        curr_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                 false)
        # test with current info
        @test InfiniteOpt._update_point_info(curr_info, divref) == curr_info
        # undo info changes
        info = VariableInfo{Float64, Float64, Float64, Float64}(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
        InfiniteOpt._update_variable_info(divref, info)
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
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                               Val(Point), parameter_values = 3,
                                               infinite_variable_ref = pref)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                               Val(Point), parameter_values = (1, 1, 1),
                                               infinite_variable_ref = ivref)
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                               Val(Point), parameter_values = (1, 2),
                                               infinite_variable_ref = ivref)
        # test a variety of builds
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).infinite_variable_ref == ivref
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).parameter_values == [0.5, 0.5]
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException InfiniteOpt._make_variable(error, info, Val(Point),
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        @test InfiniteOpt._make_variable(error, info, Val(Point), infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == [0.5, 0, 0]
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
                   parameter_values = (0.5, 0.5)).parameter_values == [0.5, 0.5]
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                             parameter_values = (0.5, 0.5)).info == info
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == [0.5, 0, 0]
        @test_throws ErrorException build_variable(error, info, Point,
                                            infinite_variable_ref = ivref2,
                                            parameter_values = (0.5, [0, 0, 0]))
    end
    # _add_point_support
    @testset "_add_point_support (IndependentParameterRef)" begin
        # test normal
        ps = [dispatch_variable_ref(pref)]
        vals = Float64[0, 2]
        @test InfiniteOpt._add_point_support(ps, vals, 1) == 2
        @test supports(pref, label = UserDefined) == [0]
        @test delete_supports(pref) isa Nothing
        # test other
        ps = [dispatch_variable_ref(pref), dispatch_variable_ref(pref2)]
        vals = Float64[0, 1]
        @test InfiniteOpt._add_point_support(ps, vals, 1) == 3
        @test supports(pref, label = UserDefined) == [0]
        @test supports(pref2, label = UserDefined) == [1]
        @test delete_supports(pref) isa Nothing
        @test delete_supports(pref2) isa Nothing
    end
    # _add_point_support
    @testset "_add_point_support (DependentParameterRef)" begin
        # test normal
        ps = dispatch_variable_ref.(prefs)
        vals = Float64[0, 0]
        @test InfiniteOpt._add_point_support(ps, vals, 1) == 3
        @test supports(prefs, label = UserDefined) == zeros(2, 1)
        @test delete_supports(prefs) isa Nothing
    end
    # _update_param_supports
    @testset "_update_param_supports" begin
        # test normal
        @test isa(InfiniteOpt._update_param_supports(divref, Float64[0.5, 1]),
                  Nothing)
        @test supports(pref) == [0.5]
        @test supports(pref2) == [1]
        # prepare array tuple
        tuple = IC.VectorTuple{Float64}(0.5, [0, 1])
        # test normal with array
        @test isa(InfiniteOpt._update_param_supports(divref2, Float64[0.5, 0, 1]),
                  Nothing)
        @test supports(pref) == [0.5]
        @test supports(prefs[1]) == [0]
        @test supports(prefs[2]) == [1]
        @test delete_supports(pref) isa Nothing
        @test delete_supports(prefs) isa Nothing
    end
    # _update_infinite_point_mapping
    @testset "_update_infinite_point_mapping" begin
        # test first addition
        idx1 = PointVariableIndex(12)
        pvref = PointVariableRef(m, idx1)
        @test isa(InfiniteOpt._update_infinite_point_mapping(pvref, divref),
                  Nothing)
        @test InfiniteOpt._point_variable_dependencies(ivref) == [idx1]
        # test second addition
        idx2 = PointVariableIndex(42)
        pvref = PointVariableRef(m, idx2)
        @test isa(InfiniteOpt._update_infinite_point_mapping(pvref, divref),
                  Nothing)
        @test InfiniteOpt._point_variable_dependencies(ivref) == [idx1, idx2]
        # undo changes
        empty!(InfiniteOpt._point_variable_dependencies(ivref))
    end
    # _check_and_make_variable_ref
    @testset "_check_and_make_variable_ref" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        @independent_parameter(m2, pref3 in [0, 1])
        @infinite_variable(m2, ivref3(pref3))
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        # test for invalid variable error
        @test_throws VariableNotOwned{InfiniteVariableRef} InfiniteOpt._check_and_make_variable_ref(m, v)
        # test normal
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0)
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m2, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m2, v) == vref
        @test supports(pref3) == [0]
        @test InfiniteOpt._point_variable_dependencies(ivref3) == [idx]
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        @independent_parameter(m2, pref3 in [0, 1])
        @infinite_variable(m2, ivref3(pref3))
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        # test for invalid variable error
        @test_throws VariableNotOwned{InfiniteVariableRef} add_variable(m, v)
        # test normal
        v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0]
        @test supports(pref2) == [1]
        @test name(vref) == "name"
        @test InfiniteOpt._point_variable_dependencies(ivref) == [idx]
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        # test info addition functions
        idx = PointVariableIndex(2)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        # lower bound
        cindex = ConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex, ScalarShape())
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = ConstraintIndex(2)
        cref = InfOptConstraintRef(m, cindex, ScalarShape())
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = ConstraintIndex(3)
        cref = InfOptConstraintRef(m, cindex, ScalarShape())
        @test has_upper_bound(vref)
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = ConstraintIndex(4)
        cref = InfOptConstraintRef(m, cindex, ScalarShape())
        @test has_upper_bound(vref)
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 1:4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        # test integer addition functions
        idx = PointVariableIndex(3)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        cindex = ConstraintIndex(8)
        cref = InfOptConstraintRef(m, cindex, ScalarShape())
        @test has_upper_bound(vref)
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.Integer}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 5:8]
    end
end

# Test the point variable macro
@testset "Macro Definition" begin
    # initialize model, parameters, and infinite variables
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @infinite_variable(m, 0 <= z(t, x) <= 1, Int)
    @infinite_variable(m, z2[1:2](t) == 3)
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @point_variable(m, infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == gvref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_integer(vref)
        @test lower_bound(vref) == 0
        # test anon with changes to fixed
        idx = PointVariableIndex(2)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @point_variable(m, infinite_variable_ref = z, lower_bound = -5,
                          parameter_values = (0, [0, 0]), binary = true) == gvref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test !is_integer(vref)
        @test is_binary(vref)
        @test lower_bound(vref) == -5
        # test regular with alias
        idx = PointVariableIndex(3)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @point_variable(m, z(0, [0, 0]), z0, Bin) == gvref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
        # test regular with semi anon
        idx = PointVariableIndex(4)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @point_variable(m, z(0, [0, 0]), base_name = "z0",
                              binary = true) == gvref
        @test infinite_variable_ref(vref) == z
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
    end
    # test array variable definition
    @testset "Array" begin
        # test anon array with one infvar
        idxs = [PointVariableIndex(5), PointVariableIndex(6)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @point_variable(m, [1:2], infinite_variable_ref = z,
                              parameter_values = (0, [0, 0])) == gvrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_integer(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        # test anon array with different inf vars
        idxs = [PointVariableIndex(7), PointVariableIndex(8)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @point_variable(m, [i = 1:2], infinite_variable_ref = z2[i],
                              parameter_values = 0) == gvrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test fix_value(vrefs[2]) == 3
        @test name(vrefs[1]) == "z2[1](0)"
        # test array with same infvar
        idxs = [PointVariableIndex(9), PointVariableIndex(10)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @point_variable(m, z(0, [0, 0]), a[1:2], Bin) == gvrefs
        @test infinite_variable_ref(vrefs[1]) == z
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_binary(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        @test name(vrefs[1]) == "a[1]"
        # test test array with differnt infvars
        idxs = [PointVariableIndex(11), PointVariableIndex(12)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @point_variable(m, z2[i](0), b[i = 1:2] >= -5) == gvrefs
        @test infinite_variable_ref(vrefs[1]) == z2[1]
        @test infinite_variable_ref(vrefs[2]) == z2[2]
        @test parameter_values(vrefs[2]) == (0,)
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "b[1]"
        # test semi anon array
        idxs = [PointVariableIndex(13), PointVariableIndex(14)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @point_variable(m, z2[i](0), [i = 1:2], lower_bound = -5) == gvrefs
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

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @independent_parameter(m, t in [0, 1])
    @dependent_parameters(m, x[1:2] in [-1, 1])
    @infinite_variable(m, y(t, x))
    @point_variable(m, y(0, [0, 0]), y0)
    vref = dispatch_variable_ref(y0)
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(y0)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test used_by_constraint(y0)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(y0)
        @test !used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = true
        @test used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = false
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test is_used(y0)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by objective
        InfiniteOpt._data_object(vref).in_objective = true
        @test is_used(vref)
    end
end
