# Test basic methods
@testset "Basics" begin
    # initialize model and point variable
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1])
    @infinite_parameter(m, b[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, c[1:2] in [0, 1])
    @variable(m, ivref, Infinite(a, b, c))
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
        @test isequal(dispatch_variable_ref(m, idx), vref)
        @test isequal(dispatch_variable_ref(gvref), vref)
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
        @test core_object(vref) === var
        @test core_object(gvref) === var
    end
    # _set_core_object
    @testset "_set_core_object" begin
        @test InfiniteOpt._set_core_object(vref, var) isa Nothing
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
        @test InfiniteOpt._constraint_dependencies(vref) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == InfOptConstraintIndex[]
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
        @test name(gvref) == "var"
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test isequal(infinite_variable_ref(vref), ivref)
        @test isequal(infinite_variable_ref(gvref), ivref)
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
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        # test default
        @test isa(set_name(gvref, "a"), Nothing)
        @test name(vref) == "a"
    end
    # GeneralVariableRef
    @testset "GeneralVariableRef" begin
        @test isequal(InfiniteOpt.GeneralVariableRef(m, idx), gvref)
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
        @test isequal(variable_by_name(m, "a"), gvref)
        @test variable_by_name(m, "test") isa Nothing
        # prepare variable with same name
        idx2 = PointVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = PointVariableRef(m, idx2)
        @test set_name(vref2, "a") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "a")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(vref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(vref)) == 1
        @test !is_valid(m, vref)
        @test_throws ErrorException InfiniteOpt._data_object(vref)
    end
end

# Test variable definition methods
@testset "Definition" begin
    # initialize model and infinite variables
    m = InfiniteModel()
    @infinite_parameter(m, pref in [0, 1])
    @infinite_parameter(m, pref2 in [0, 1])
    @infinite_parameter(m, prefs[1:2] in [0, 1])
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, num, false, true)
    info4 = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    num = Float64(1)
    info5 = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    @variable(m, ivref, Infinite(pref, pref2))
    @variable(m, ivref2, Infinite(pref, prefs))
    divref = dispatch_variable_ref(ivref)
    divref2 = dispatch_variable_ref(ivref2)
    # test Point
    @testset "Point{V, T}" begin 
        @test isequal(Point(ivref, [0, 0]).infinite_variable_ref, ivref)
        @test Point(ivref, [0, 0]).parameter_values.values == [0, 0]
    end
    # _check_tuple_shape
    @testset "_check_tuple_shape" begin
        # test normal
        @test isa(InfiniteOpt._check_tuple_shape(error, divref, IC.VectorTuple(0.5, 0.5)),
                  Nothing)
        # prepare param value tuple
        tuple = IC.VectorTuple(0.5, [0.5, 0.5])
        # test normal with array
        @test isa(InfiniteOpt._check_tuple_shape(error, divref2, tuple), Nothing)
        # test for errors in shape
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref,
                                                              IC.VectorTuple(0.5,))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref,
                                                        IC.VectorTuple(0.5, [0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, divref2,
                                                          IC.VectorTuple(0.5, 0.5))
        tuple = IC.VectorTuple(0.5, [0.5, 0.5, 0.5])
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
        basic_func = (a::Vector) -> 1
        # prepare info for test
        new_info = VariableInfo(true, 0., true, 0., false, 0., false, basic_func,
                                 true, false)
        InfiniteOpt._update_variable_info(divref, new_info)
        expected = VariableInfo(true, 0., true, 0., false, 0., false, 0., true, false)
        # test with current info
        @test InfiniteOpt._update_point_info(info, divref, Float64[0, 0]) == expected
        # prepare info for test
        new_info = VariableInfo(false, 0., false, 0., true, 0., true, basic_func,
                                false, true)
        InfiniteOpt._update_variable_info(divref, new_info)
        expected = VariableInfo(false, 0., false, 0., true, 0., true, 1, false, true)
        # test with current info
        @test InfiniteOpt._update_point_info(info, divref, Float64[0, 0]) == expected
        # prepare info for test
        curr_info = VariableInfo(true, 0., true, 0., false, 0., true, 0.,
                                 true, false)
        # test with current info
        @test InfiniteOpt._update_point_info(curr_info, divref, Float64[0, 0]) == curr_info
        # undo info changes
        basic_func = (a::Vector) -> NaN
        old_info = VariableInfo(false, NaN, false, NaN, false, NaN, false, basic_func,
                            false, false)
        InfiniteOpt._update_variable_info(divref, old_info)
        # test with user defined start function
        @test set_start_value_function(divref, (a, b) -> a + b) isa Nothing
        expected = VariableInfo(false, 0., false, 0., false, 0., true, 3, false, false)
        @test InfiniteOpt._update_point_info(info, divref, Float64[1, 2]) == expected
        @test reset_start_value_function(divref) isa Nothing
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for all errors
        @test_throws ErrorException build_variable(error, info, Point(2))
        @test_throws ErrorException build_variable(error, info, Point(pref, 0))
        @test_throws ErrorException build_variable(error, info, Point(ivref))
        @test_throws ErrorException build_variable(error, info, Point(ivref, "d"))
        @test_throws ErrorException build_variable(error, info, Point(ivref, 0.5, 0.5), 
                                                   bad = 42)
        # test a variety of builds
        @test isequal(build_variable(error, info, Point(ivref, 0.5, 0.5)).infinite_variable_ref, ivref)
        @test build_variable(error, info, Point(ivref, 0.5, 0.5)).parameter_values == [0.5, 0.5]
        @test build_variable(error, info, Point(ivref, 0.5, 0.5)).info == info
        @test_throws ErrorException build_variable(error, info, Point(ivref, 0.5, 2))
        @test isequal(build_variable(error, info, Point(ivref2, 0.5, [0, 0])).infinite_variable_ref, ivref2)
        @test build_variable(error, info, Point(ivref2, 0.5, [0, 0])).parameter_values == [0.5, 0, 0]
        @test_throws ErrorException build_variable(error, info, Point(ivref2, 0.5, [0, 0, 0]))
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
        tuple = IC.VectorTuple(0.5, [0, 1])
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
    # add_variable
    @testset "JuMP.add_variable" begin
        # prepare secondary model and infinite variable
        m2 = InfiniteModel()
        @infinite_parameter(m2, pref3 in [0, 1])
        @variable(m2, ivref3, Infinite(pref3))
        v = build_variable(error, info, Point(ivref3, 0.5))
        # test for invalid variable error
        @test_throws VariableNotOwned{InfiniteVariableRef} add_variable(m, v)
        # test normal
        v = build_variable(error, info, Point(ivref, 0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0]
        @test supports(pref2) == [1]
        @test name(vref) == "name"
        @test InfiniteOpt._point_variable_dependencies(ivref) == [idx]
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Point(ivref, 0, 0))
        # test info addition functions
        idx = PointVariableIndex(2)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test !transformation_backend_ready(m)
        # lower bound
        cindex = InfOptConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex)
        @test has_lower_bound(vref)
        @test InfiniteOpt._lower_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = InfOptConstraintIndex(2)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test InfiniteOpt._upper_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = InfOptConstraintIndex(3)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test is_fixed(vref)
        @test InfiniteOpt._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = InfOptConstraintIndex(4)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test is_binary(vref)
        @test InfiniteOpt._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [InfOptConstraintIndex(i)
                                                             for i = 1:4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Point(ivref, 1, 1))
        # test integer addition functions
        idx = PointVariableIndex(3)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test !transformation_backend_ready(m)
        cindex = InfOptConstraintIndex(8)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test is_integer(vref)
        @test InfiniteOpt._integer_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.Integer}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [InfOptConstraintIndex(i)
                                                             for i = 5:8]
        # test redundant add with same info 
        v = build_variable(error, info, Point(ivref, 0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name2"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0, 1]
        @test supports(pref2) == [0, 1]
        @test name(vref) == "name2"
        @test !has_upper_bound(vref)
        # test redundant add with all new info 
        v = build_variable(error, info4, Point(ivref, 0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0, 1]
        @test supports(pref2) == [0, 1]
        @test name(vref) == "name2"
        @test lower_bound(vref) == 0
        @test upper_bound(vref) == 0
        @test fix_value(vref) == 0
        @test start_value(vref) == 0
        @test is_binary(vref)
        @test is_integer(vref)
        @test LowerBoundRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(9))
        @test UpperBoundRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(10))
        @test FixRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(11))
        @test BinaryRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(12))
        @test IntegerRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(13))
        # test redundant add with all new values info 
        v = build_variable(error, info5, Point(ivref, 0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name4"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0, 1]
        @test supports(pref2) == [0, 1]
        @test name(vref) == "name4"
        @test lower_bound(vref) == 1
        @test upper_bound(vref) == 1
        @test fix_value(vref) == 1
        @test start_value(vref) == 1
        @test is_binary(vref)
        @test is_integer(vref)
        @test LowerBoundRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(9))
        @test UpperBoundRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(10))
        @test FixRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(11))
        @test BinaryRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(12))
        @test IntegerRef(vref) == InfOptConstraintRef(m, InfOptConstraintIndex(13))
        # test redundant add with new info that deletes 
        v = build_variable(error, info, Point(ivref, 0, 1))
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name5"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0, 1]
        @test supports(pref2) == [0, 1]
        @test name(vref) == "name5"
        @test !has_lower_bound(vref)
        @test !has_upper_bound(vref)
        @test !is_fixed(vref)
        @test start_value(vref) === nothing
        @test !is_binary(vref)
        @test !is_integer(vref)
        # test add with dependent parameters 
        v = build_variable(error, info, Point(ivref2, 0, [1, 1]))
        idx = PointVariableIndex(4)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test supports(pref) == [0, 1]
        @test supports(prefs) == ones(2, 1)
        @test name(vref) == "name"
    end
end

# Test the point variable macro
@testset "Macro Definition" begin
    # initialize model, parameters, and infinite variables
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, 0 <= z <= 1, Infinite(t, x), Int, start = (a, x) -> a + sum(x))
    @variable(m, z2[1:2] == 3, Infinite(t))
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Point(z, 0, [0, 0])), gvref)
        @test isequal(infinite_variable_ref(vref), z)
        @test parameter_values(vref) == (0, [0, 0])
        @test is_integer(vref)
        @test lower_bound(vref) == 0
        @test start_value(vref) == 0
        # test anon with changes to fixed
        idx = PointVariableIndex(2)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Point(z, 0, [1, 1]), 
                        lower_bound = -5, binary = true), gvref)
        @test isequal(infinite_variable_ref(vref), z)
        @test parameter_values(vref) == (0, [1, 1])
        @test !is_integer(vref)
        @test is_binary(vref)
        @test lower_bound(vref) == -5
        # test regular with alias
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, z0, Point(z, 0, [0, 0]), Bin), gvref)
        @test isequal(infinite_variable_ref(vref), z)
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
        # test regular with semi anon
        idx = PointVariableIndex(1)
        vref = PointVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Point(z, 0, [0, 0]), base_name = "z0",
                        binary = true), gvref)
        @test isequal(infinite_variable_ref(vref), z)
        @test parameter_values(vref) == (0, [0, 0])
        @test is_binary(vref)
        @test lower_bound(vref) == 0
        @test name(vref) == "z0"
    end
    # test array variable definition
    @testset "Array" begin
        # test anon array with one infvar
        idxs = [PointVariableIndex(1), PointVariableIndex(1)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [1:2], Point(z, 0, [0, 0])), gvrefs)
        @test isequal(infinite_variable_ref(vrefs[1]), z)
        @test parameter_values(vrefs[2]) == (0, [0, 0])
        @test is_integer(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        # test anon array with different inf vars
        idxs = [PointVariableIndex(3), PointVariableIndex(4)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [i = 1:2], Point(z2[i], 0)), gvrefs)
        @test isequal(infinite_variable_ref(vrefs[1]), z2[1])
        @test isequal(infinite_variable_ref(vrefs[2]), z2[2])
        @test parameter_values(vrefs[2]) == (0,)
        @test fix_value(vrefs[2]) == 3
        @test name(vrefs[1]) == ""
        # test array with same infvar
        idxs = [PointVariableIndex(1), PointVariableIndex(5)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, a[i = 1:2], Point(z, -1 + i, [0, 0]), Bin), gvrefs)
        @test isequal(infinite_variable_ref(vrefs[1]), z)
        @test parameter_values(vrefs[2]) == (1, [0, 0])
        @test is_binary(vrefs[1])
        @test lower_bound(vrefs[2]) == 0
        @test name(vrefs[1]) == "a[1]"
        # test test array with differnt infvars
        idxs = [PointVariableIndex(3), PointVariableIndex(4)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, b[i = 1:2] >= -5, Point(z2[i], 0)), gvrefs)
        @test isequal(infinite_variable_ref(vrefs[1]), z2[1])
        @test isequal(infinite_variable_ref(vrefs[2]), z2[2])
        @test parameter_values(vrefs[2]) == (0,)
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == "b[1]"
        # test semi anon array
        idxs = [PointVariableIndex(6), PointVariableIndex(7)]
        vrefs = [PointVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [i = 1:2], Point(z2[i], 1), lower_bound = -5), gvrefs)
        @test isequal(infinite_variable_ref(vrefs[1]), z2[1])
        @test isequal(infinite_variable_ref(vrefs[2]), z2[2])
        @test lower_bound(vrefs[2]) == -5
        @test name(vrefs[1]) == ""
        @test parameter_values.(vrefs) == [(1,), (1,)]
    end
    # test errors
    @testset "Errors" begin
        # test other syntaxes
        @test_macro_throws ErrorException @variable(m, w, Point(z2[1], 0), 
                                                    bad = 1)
        # test redefinition catch
        @test_macro_throws ErrorException @variable(m, z0, Point(z, 0, [0, 0]))
    end
end

# test restriciton definition 
@testset "Restriction Definition" begin 
    # setup the info 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [0, 1])
    @variable(m, y, Infinite(t, x))
    # test errors 
    @test_throws ErrorException restrict(y, 0, [2, 2])
    @test_throws ErrorException restrict(y, 0, 2)
    @test_throws ErrorException y(0, 2)
    # test normal wth restrict
    vref = GeneralVariableRef(m, 1, PointVariableIndex)
    @test isequal(restrict(y, 0, [0, 0]), vref)
    @test parameter_values(vref) == (0, [0, 0])
    # test normal functionally
    vref = GeneralVariableRef(m, 2, PointVariableIndex)
    @test isequal(y(0.5, [0, 0]), vref)
    @test parameter_values(vref) == (0.5, [0, 0])
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, y, Infinite(t, x))
    @variable(m, y0, Point(y, 0, [0, 0]))
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
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
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
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test is_used(y0)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by objective
        InfiniteOpt._data_object(vref).in_objective = true
        @test is_used(vref)
    end
end
