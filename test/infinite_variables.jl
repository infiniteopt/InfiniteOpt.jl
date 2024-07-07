# Test name methods
@testset "Basics" begin
    # setup data
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1])
    @infinite_parameter(m, b[1:3] in [0, 1], independent = true)
    @infinite_parameter(m, c[1:2] in [0, 1])
    @infinite_parameter(m, d in [0, 1])
    idx = InfiniteVariableIndex(1)
    num = Float64(0)
    func = (x) -> NaN
    func2 = (x...) -> 1 
    info = VariableInfo(false, num, false, num, false, num, false, func, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, func2, true, false)
    new_info2 = VariableInfo(true, 0, true, 0, true, 0, true, func2, true, false)
    var = InfiniteVariable(info, IC.VectorTuple(a, b[1:2], c, [b[3], d]),
                           [1:7...], [1:6...], true)
    var2 = InfiniteVariable(new_info, IC.VectorTuple(a, b[1:2], c, [b[3], d]),
                            [1:7...], [1:6...], true)
    object = VariableData(var, "var")
    vref = InfiniteVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
    bad_vref = InfiniteVariableRef(m, InfiniteVariableIndex(-1))
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
        @test InfiniteOpt._data_dictionary(m, InfiniteVariable) === m.infinite_vars
        @test InfiniteOpt._data_dictionary(vref) === m.infinite_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.infinite_vars
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
        @test_throws ErrorException InfiniteOpt._data_object(bad_vref)
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) === var
        @test InfiniteOpt._core_variable_object(gvref) === var
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
        @test InfiniteOpt._set_core_variable_object(vref, var2) isa Nothing
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
    end
    # parameter_group_int_indices
    @testset "parameter_group_int_indices" begin
        @test InfiniteOpt.parameter_group_int_indices(vref) == [1:6...]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(vref) == [1:7...]
    end
    # test _variable_info
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(vref) == info
    end
    # test _is_vector_start
    @testset "_is_vector_start" begin
        @test InfiniteOpt._is_vector_start(vref)
    end
    # _update_variable_info
    @testset "_update_variable_info" begin
        @test isa(InfiniteOpt._update_variable_info(vref, new_info2), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info2
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
    # _semi_infinite_variable_dependencies
    @testset "_semi_infinite_variable_dependencies" begin
        @test InfiniteOpt._semi_infinite_variable_dependencies(vref) == SemiInfiniteVariableIndex[]
        @test InfiniteOpt._semi_infinite_variable_dependencies(gvref) == SemiInfiniteVariableIndex[]
    end
    # _point_variable_dependencies
    @testset "_point_variable_dependencies" begin
        @test InfiniteOpt._point_variable_dependencies(vref) == PointVariableIndex[]
        @test InfiniteOpt._point_variable_dependencies(gvref) == PointVariableIndex[]
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(vref) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(gvref) == DerivativeIndex[]
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
        @test name(gvref) == "var"
        @test name(bad_vref) == ""
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        @test isequal(raw_parameter_refs(vref), IC.VectorTuple(a, b[1:2], c, [b[3], d]))
        @test isequal(raw_parameter_refs(gvref), IC.VectorTuple(a, b[1:2], c, [b[3], d]))
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test isequal(parameter_refs(vref), (a, b[1:2], c, [b[3], d]))
        @test isequal(parameter_refs(gvref), (a, b[1:2], c, [b[3], d]))
    end
    # parameter_list
    @testset "parameter_list" begin
        @test isequal(parameter_list(vref), [a; b[1:2]; c; [b[3], d]])
        @test isequal(parameter_list(gvref), [a; b[1:2]; c; [b[3], d]])
        @test isequal(parameter_list(raw_parameter_refs(vref)), [a; b[1:2]; c; [b[3], d]])
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test default
        @test isa(set_name(vref, ""), Nothing)
        @test name(vref) == ""
        @test_throws ErrorException set_name(bad_vref, "")
        # test normal
        @test isa(set_name(gvref, "new"), Nothing)
        @test name(vref) == "new"
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
        @test isequal(variable_by_name(m, "new"), gvref)
        @test variable_by_name(m, "test") isa Nothing
        # prepare variable with same name
        idx2 = InfiniteVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = InfiniteVariableRef(m, idx2)
        @test set_name(vref2, "new") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "new")
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
    # initialize model and infinite variable info
    m = InfiniteModel()
    @infinite_parameter(m, pref in [0, 1])
    @infinite_parameter(m, pref2 in [0, 1])
    @infinite_parameter(m, prefs[1:2] in [0, 1])
    @finite_parameter(m, fin == 42)
    num = Float64(0)
    func1 = (a) -> 2
    info = VariableInfo(false, num, false, num, false, num, false, NaN, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, func1, true, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, 42, false, true)
    # test Infinite 
    @testset "Infinite{VT}" begin 
        @test Infinite(pref, prefs).parameter_refs isa IC.VectorTuple
        @test isequal(Infinite(pref).parameter_refs.values, [pref])
        @test isequal(Infinite(pref, prefs, pref2).parameter_refs.values, [pref, prefs..., pref2])
    end
    # _check_tuple_element (IndependentParameterRefs)
    @testset "_check_tuple_element (IndependentParameterRefs)" begin
        iprefs = dispatch_variable_ref.([pref, pref2])
        @test InfiniteOpt._check_tuple_element(error, iprefs) isa Nothing
    end
    # _check_tuple_element (DependentParameterRefs)
    @testset "_check_tuple_element (DependentParameterRefs)" begin
        dprefs = dispatch_variable_ref.(prefs)
        @test InfiniteOpt._check_tuple_element(error, dprefs) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_tuple_element(error, dprefs[1:1])
    end
    # _check_tuple_element (Fallback)
    @testset "_check_tuple_element (Fallback)" begin
        refs = dispatch_variable_ref.([fin, pref])
        @test_throws ErrorException InfiniteOpt._check_tuple_element(error, refs)
    end
    # _check_parameter_tuple
    @testset "_check_parameter_tuple" begin
        # test normal
        tuple = IC.VectorTuple((pref, prefs, pref2))
        @test InfiniteOpt._check_parameter_tuple(error, tuple) isa Nothing
        tuple = IC.VectorTuple(([pref, pref2], prefs))
        @test InfiniteOpt._check_parameter_tuple(error, tuple) isa Nothing
        # test bad
        tuple = IC.VectorTuple((pref, pref))
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error, tuple)
        tuple = IC.VectorTuple((pref, [prefs[1], pref2]))
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error, tuple)
        tuple = IC.VectorTuple(fin)
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error, tuple)
        tuple = IC.VectorTuple(2)
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error, tuple)
    end
    # _check_and_format_infinite_info (Real)
    @testset "_check_and_format_infinite_info (Real)" begin
        tuple = IC.VectorTuple((pref, prefs, pref2))
        @test !InfiniteOpt._check_and_format_infinite_info(error, info, tuple)[1].has_start
        @test InfiniteOpt._check_and_format_infinite_info(error, info, tuple)[1].start isa Function
        @test InfiniteOpt._check_and_format_infinite_info(error, info, tuple)[2]
        @test InfiniteOpt._check_and_format_infinite_info(error, info3, tuple)[1].has_start
        @test InfiniteOpt._check_and_format_infinite_info(error, info3, tuple)[1].start isa Function
        @test InfiniteOpt._check_and_format_infinite_info(error, info3, tuple)[2]
    end
    # _check_and_format_infinite_info (Function)
    @testset "_check_and_format_infinite_info (Function)" begin
        tuple = IC.VectorTuple((pref,))
        # test normal
        @test InfiniteOpt._check_and_format_infinite_info(error, info2, tuple)[1].has_start
        @test InfiniteOpt._check_and_format_infinite_info(error, info2, tuple)[1].start == func1
        @test !InfiniteOpt._check_and_format_infinite_info(error, info2, tuple)[2]
        # test multiple methods
        mult_func(a,b,c) = 2
        mult_func(a) = 2
        mult_info = VariableInfo(true, num, true, num, true, num, true, mult_func, true, false)
        @test InfiniteOpt._check_and_format_infinite_info(error, mult_info, tuple)[1] isa VariableInfo 
        # test bad method
        bad_func(a, b) = 2
        bad_info = VariableInfo(true, num, true, num, true, num, true, bad_func, true, false)
        @test_throws ErrorException InfiniteOpt._check_and_format_infinite_info(error, bad_info, tuple)
        # test bad args
        func = (a, b) -> a + sum(b)
        bad_info = VariableInfo(true, num, true, num, true, num, true, func, true, false)
        @test_throws ErrorException InfiniteOpt._check_and_format_infinite_info(error, bad_info, tuple)
    end
    # _check_and_format_infinite_info (Fallback)
    @testset "_check_and_format_infinite_info (Fallback)" begin
        tuple = IC.VectorTuple((pref, prefs, pref2))
        bad_info = VariableInfo(true, num, true, num, true, num, true, Complex(1), true, false)
        @test_throws ErrorException InfiniteOpt._check_and_format_infinite_info(error, bad_info, tuple)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for each error message
        @test_throws ErrorException build_variable(error, info, Infinite(pref),
                                                   bob = 42)
        @test_throws ErrorException build_variable(error, info, Infinite())
        @test_throws ErrorException build_variable(error, info, Infinite(pref, fin))
        @test_throws ErrorException build_variable(error, info, Infinite(pref, pref))
        # test for expected output
        @test build_variable(error, info, Infinite(pref)).info isa JuMP.VariableInfo{Float64, Float64, Float64, <:Function}
        @test isequal(build_variable(error, info, Infinite(pref)).parameter_refs, IC.VectorTuple(pref))
        # test various types of param tuples
        tuple = IC.VectorTuple(pref, prefs)
        @test isequal(build_variable(error, info, Infinite(pref, prefs)).parameter_refs, tuple)
        tuple = IC.VectorTuple(prefs)
        @test isequal(build_variable(error, info, Infinite(prefs)).parameter_refs, tuple)
        @test isequal(build_variable(error, info, Infinite(prefs)).parameter_nums, [3, 4])
        @test isequal(build_variable(error, info, Infinite(prefs)).object_nums, [3])
    end
    # _check_parameters_valid
    @testset "_check_parameters_valid" begin
        # prepare param tuple
        @infinite_parameter(InfiniteModel(), pref3 in [0, 1])
        tuple = IC.VectorTuple((pref, prefs, pref3))
        # test that catches error
        @test_throws VariableNotOwned{GeneralVariableRef} InfiniteOpt._check_parameters_valid(m, tuple)
        # test normal
        tuple = IC.VectorTuple(pref, pref2)
        @test InfiniteOpt._check_parameters_valid(m, tuple) isa Nothing
    end
    # _update_param_var_mapping
    @testset "_update_param_var_mapping" begin
        # prepare tuple
        tuple = IC.VectorTuple(pref, prefs)
        idx = InfiniteVariableIndex(1)
        ivref = InfiniteVariableRef(m, idx)
        # test normal
        @test InfiniteOpt._update_param_var_mapping(ivref, tuple) isa Nothing
        @test InfiniteOpt._infinite_variable_dependencies(pref) == [idx]
        @test InfiniteOpt._infinite_variable_dependencies(prefs[1]) == [idx]
        @test InfiniteOpt._infinite_variable_dependencies(prefs[2]) == [idx]
        # undo changes
        empty!(InfiniteOpt._infinite_variable_dependencies(pref))
        empty!(InfiniteOpt._infinite_variable_dependencies(prefs[1]))
        empty!(InfiniteOpt._infinite_variable_dependencies(prefs[2]))
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        # prepare secondary model and parameter and variable
        m2 = InfiniteModel()
        @infinite_parameter(m2, pref3 in [0, 1])
        v = build_variable(error, info, Infinite(pref3))
        # test for error of invalid variable
        @test_throws VariableNotOwned{GeneralVariableRef} add_variable(m, v)
        # prepare normal variable
        v = build_variable(error, info, Infinite(pref))
        # test normal
        idx = InfiniteVariableIndex(1)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test InfiniteOpt._core_variable_object(vref) == v
        @test InfiniteOpt._infinite_variable_dependencies(pref) == [idx]
        @test name(vref) == "name"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Infinite(pref))
        # test info addition functions
        idx = InfiniteVariableIndex(2)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test !transformation_backend_ready(m)
        @test !InfiniteOpt._is_vector_start(vref)
        @test InfiniteOpt._variable_info(vref).start == func1
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
        @test is_fixed(vref)
        @test InfiniteOpt._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = InfOptConstraintIndex(4)
        cref = InfOptConstraintRef(m, cindex)
        @test is_binary(vref)
        @test InfiniteOpt._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [InfOptConstraintIndex(i)
                                                             for i = 1:4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Infinite(pref))
        # test integer addition functions
        idx = InfiniteVariableIndex(3)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, v, "name"), gvref)
        @test !transformation_backend_ready(m)
        @test InfiniteOpt._is_vector_start(vref)
        @test InfiniteOpt._variable_info(vref).start isa Function
        cindex = InfOptConstraintIndex(8)
        cref = InfOptConstraintRef(m, cindex)
        @test is_integer(vref)
        @test InfiniteOpt._integer_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.Integer}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [InfOptConstraintIndex(i)
                                                             for i = 5:8]
    end
end

# Test the infinite variable macro
@testset "Macro Definition" begin
    # initialize model and infinite parameters
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    # test single variable definition
    @testset "Single" begin
        # test basic defaults
        idx = InfiniteVariableIndex(1)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Infinite(t)), gvref)
        @test name(vref) == ""
        # test more tuple input and variable details
        idx = InfiniteVariableIndex(2)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Infinite(t, x), base_name = "test",
                        binary = true), gvref)
        @test name(vref) == "test"
        @test is_binary(vref)
        # test nonanonymous with simple single arg
        idx = InfiniteVariableIndex(3)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, a, Infinite(x)), gvref)
        @test name(vref) == "a"
        # test nonanonymous with complex single arg
        idx = InfiniteVariableIndex(4)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, 0 <= b <= 1, Infinite(x)), gvref)
        @test name(vref) == "b"
        @test lower_bound(vref) == 0
        @test upper_bound(vref) == 1
        # test nonanonymous with reversed single arg
        idx = InfiniteVariableIndex(5)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, 0 <= c, Infinite(t)), gvref)
        @test name(vref) == "c"
        @test lower_bound(vref) == 0
        # test multi-argument expr 1
        idx = InfiniteVariableIndex(6)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, d == 0, Infinite(t), Int, base_name = "test"), gvref)
        @test name(vref) == "test"
        @test fix_value(vref) == 0
        # test single start value
        idx = InfiniteVariableIndex(7)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        @test isequal(@variable(m, aa, Infinite(t), Int, start = 42), gvref)
        @test name(vref) == "aa"
        @test start_value_function(vref) isa Function
        @test start_value_function(vref)([0]) == 42
        @test_throws MethodError start_value_function(vref)(0)
        # test function start value
        idx = InfiniteVariableIndex(8)
        vref = InfiniteVariableRef(m, idx)
        gvref = InfiniteOpt.GeneralVariableRef(m, idx)
        func = (a) -> a
        @test isequal(@variable(m, ab, Infinite(t), Int, start = func), gvref)
        @test name(vref) == "ab"
        @test start_value_function(vref) isa Function
        @test start_value_function(vref)(0) == 0
    end
    # test array variable definition
    @testset "Array" begin
        # test anonymous array
        idxs = [InfiniteVariableIndex(9), InfiniteVariableIndex(10)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [1:2], variable_type = Infinite(t)), gvrefs)
        @test name(vrefs[1]) == ""
        # test basic param expression
        idxs = [InfiniteVariableIndex(11), InfiniteVariableIndex(12)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, e[1:2], variable_type = Infinite(t, x)), gvrefs)
        @test name(vrefs[2]) == "e[2]"
        # test comparison without params
        idxs = [InfiniteVariableIndex(13), InfiniteVariableIndex(14)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, 0 <= f[1:2] <= 1, Infinite(t, x)), gvrefs)
        @test name(vrefs[2]) == "f[2]"
        @test lower_bound(vrefs[1]) == 0
        @test upper_bound(vrefs[2]) == 1
        # test comparison with call
        idxs = [InfiniteVariableIndex(15), InfiniteVariableIndex(16)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test @variable(m, 0 <= g[1:2] <= 1, Infinite(t)) == gvrefs
        @test name(vrefs[1]) == "g[1]"
        @test lower_bound(vrefs[1]) == 0
        @test upper_bound(vrefs[2]) == 1
        # test fixed
        idxs = [InfiniteVariableIndex(17), InfiniteVariableIndex(18)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, h[i = 1:2] == ones(2)[i], Infinite(t)), gvrefs)
        @test name(vrefs[1]) == "h[1]"
        @test fix_value(vrefs[1]) == 1
        # test containers
        idxs = [InfiniteVariableIndex(19), InfiniteVariableIndex(20)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        svrefs = convert(JuMP.Containers.SparseAxisArray, gvrefs)
        @test @variable(m, [1:2], Infinite(t),
                        container = SparseAxisArray) isa JuMPC.SparseAxisArray
        @test name(vrefs[1]) == ""
        # test scalar starts
        idxs = [InfiniteVariableIndex(21), InfiniteVariableIndex(22)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        starts = [2, 5]
        @test isequal(@variable(m, w[i = 1:2], Infinite(t), start = starts[i]), gvrefs)
        @test start_value_function(vrefs[1])([0]) == 2
        @test start_value_function(vrefs[2])([0]) == 5
        # test mixed starts
        idxs = [InfiniteVariableIndex(23), InfiniteVariableIndex(24)]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        starts = [(a) -> 2, 5]
        @test isequal(@variable(m, cat[i = 1:2], Infinite(t), 
                        start = starts[i]), gvrefs)
        @test start_value_function(vrefs[1])(0) == 2
        @test start_value_function(vrefs[2])([0]) == 5
        # test with set 
        idxs = [InfiniteVariableIndex(i) for i in 25:27]
        vrefs = [InfiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt.GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(gvrefs, @variable(m, [1:3], Infinite(t), set = SecondOrderCone()))
        @test num_constraints(m, Vector{GeneralVariableRef}, MOI.SecondOrderCone) == 1
    end
    # test errors
    @testset "Errors" begin
        # test invalid keyword arguments
        @test_macro_throws ErrorException @variable(m, i, Infinite(t),
                                                    parameter_values = 1)
        # test name duplication
        @test_macro_throws ErrorException @variable(m, a, Infinite(t), Int)
        # test bad start function
        @test_macro_throws ErrorException @variable(m, zz, Infinite(t), 
                                                    start = (a, b) -> 2)
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
    @test_throws ErrorException restrict(y, t)
    @test_throws ErrorException restrict(t, x)
    @test_throws ErrorException y(t)
    # test warning 
    warn = "Unnecessary use of functional infinite variable restriction syntax " *
           "that will cause performance degredations. This was probably caused " *
           "by using syntax like `y(t, x)` inside expressions. Instead just " *
           "use the infinite variable reference (e.g. `y`)."
    @test_logs (:warn, warn) restrict(y, t, x)
    @test_logs (:warn, warn) y(t, x)
    # Note: the rest of the cases are tested with point variables and semi-infinite variables
end

# test collocation restrictions
@testset "constant_over_collocation" begin 
    # setup the info 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, xi[1:2] in [0, 1])
    @variable(m, y, Infinite(t, x, xi))
    @variable(m, q, Infinite(t))
    # test errors 
    @test_throws ErrorException constant_over_collocation(y, xi[1])
    @test_throws ErrorException constant_over_collocation(q, x[1])
    # test normal 
    @test constant_over_collocation(y, t) isa Nothing
    @test constant_over_collocation(q, t) isa Nothing
    @test m.piecewise_vars[index(t)] == Set([index(y), index(q)])
    @test constant_over_collocation(y, x[2]) isa Nothing
    @test m.piecewise_vars[index(x[2])] == Set(index(y))
end

# Test variable(s) constrained on creation 
@testset "Creation Constraints" begin 
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
    v1 = build_variable(error, info, Infinite(t))
    v2 = build_variable(error, info, Infinite(t, x))
    # test add_variable with VariableConstrainedOnCreation
    @testset "Adding VariableConstrainedOnCreation" begin 
        v = VariableConstrainedOnCreation(v1, MOI.Integer())
        vref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
        @test isequal(add_variable(m, v, "test"), vref)
        @test name(vref) == "test"
        @test num_constraints(m, MOI.Integer) == 1
    end
    # test add_variable with VariablesConstrainedOnCreation
    @testset "Adding VariablesConstrainedOnCreation" begin 
        v = VariablesConstrainedOnCreation([v1, v2, v1], MOI.SecondOrderCone(3), 
                                           VectorShape())
        vrefs = [GeneralVariableRef(m, i, InfiniteVariableIndex) for i in 2:4]
        names = ["test[$i]" for i in 1:3]
        @test isequal(add_variable(m, v, names), vrefs)
        @test name.(vrefs) == names
        @test num_constraints(m, MOI.SecondOrderCone) == 1
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, y, Infinite(t, x))
    vref = dispatch_variable_ref(y)
    # test used_by_semi_infinite_variable
    @testset "used_by_semi_infinite_variable" begin
        @test !used_by_semi_infinite_variable(vref)
        push!(InfiniteOpt._semi_infinite_variable_dependencies(vref),
              SemiInfiniteVariableIndex(1))
        @test used_by_semi_infinite_variable(y)
        @test used_by_semi_infinite_variable(vref)
        empty!(InfiniteOpt._semi_infinite_variable_dependencies(vref))
    end
    # test used_by_point_variable
    @testset "used_by_point_variable" begin
        @test !used_by_point_variable(vref)
        push!(InfiniteOpt._point_variable_dependencies(vref), PointVariableIndex(1))
        @test used_by_point_variable(y)
        @test used_by_point_variable(vref)
        empty!(InfiniteOpt._point_variable_dependencies(vref))
    end
    # test used_by_derivative
    @testset "used_by_derivative" begin
        @test !used_by_derivative(vref)
        push!(InfiniteOpt._derivative_dependencies(vref), DerivativeIndex(1))
        @test used_by_derivative(y)
        @test used_by_derivative(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
    end
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(y)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test used_by_constraint(y)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(y)
        @test !used_by_objective(vref)
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test is_used(y)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by point variable
        num = Float64(0)
        info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
        var = PointVariable(info, y, [0., 0., 0.])
        object = VariableData(var, "var")
        idx = PointVariableIndex(1)
        pvref = PointVariableRef(m, idx)
        @test InfiniteOpt._add_data_object(m, object) == idx
        push!(InfiniteOpt._point_variable_dependencies(vref), idx)
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(pvref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._point_variable_dependencies(vref))
        # test used by semi-infinite variable
        eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1)
        var = SemiInfiniteVariable(y, eval_supps, [2], [2])
        object = VariableData(var, "var")
        idx = SemiInfiniteVariableIndex(1)
        rvref = SemiInfiniteVariableRef(m, idx)
        @test InfiniteOpt._add_data_object(m, object) == idx
        push!(InfiniteOpt._semi_infinite_variable_dependencies(vref), idx)
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(rvref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._semi_infinite_variable_dependencies(vref))
        # test used by derivative
        func = (x) -> NaN
        num = 0.
        info = VariableInfo(true, num, true, num, true, num, false, func, true, true)
        deriv = Derivative(info, true, y, t, 1)
        object = VariableData(deriv)
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        @test InfiniteOpt._add_data_object(m, object) == idx
        push!(InfiniteOpt._derivative_dependencies(vref), idx)
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(dref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
    end
end
