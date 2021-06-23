# Test the basics 
@testset "Basics" begin
    # setup data
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1], derivative_method = TestMethod())
    @variable(m, x, Infinite(a))
    idx = DerivativeIndex(1)
    num = Float64(0)
    func = (x) -> NaN
    func2 = (a...) -> 1
    info = VariableInfo(false, num, false, num, false, num, false, func, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, func2, false, false)
    deriv = Derivative(info, true, x, a)
    deriv2 = Derivative(new_info, true, x, a)
    object = VariableData(deriv, "var")
    dref = DerivativeRef(m, idx)
    gvref = GeneralVariableRef(m, 1, DerivativeIndex)
    bad_dref = DerivativeRef(m, DerivativeIndex(-1))
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(dref) === m
        @test owner_model(gvref) === m
    end
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(dref) == idx
        @test index(gvref) == idx
    end
    # dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == dref
        @test dispatch_variable_ref(gvref) == dref
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(m, Derivative) === m.derivatives
        @test InfiniteOpt._data_dictionary(dref) === m.derivatives
        @test InfiniteOpt._data_dictionary(gvref) === m.derivatives
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, dref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(dref) === object
        @test InfiniteOpt._data_object(gvref) === object
        @test_throws ErrorException InfiniteOpt._data_object(bad_dref)
    end
    # _existing_derivative_index
    @testset "_existing_derivative_index" begin
        @test InfiniteOpt._existing_derivative_index(x, a) === nothing
        m.deriv_lookup[(x, a)] = idx
        @test InfiniteOpt._existing_derivative_index(x, a) === idx
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(dref) === deriv
        @test InfiniteOpt._core_variable_object(gvref) === deriv
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(dref, deriv) isa Nothing
        @test InfiniteOpt._set_core_variable_object(dref, deriv2) isa Nothing
        @test InfiniteOpt._set_core_variable_object(dref, deriv) isa Nothing
    end
    # _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(dref) == [1]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(dref) == [1]
    end
    # test _variable_info
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(dref) == info
    end
    # test _is_vector_start
    @testset "_is_vector_start" begin
        @test InfiniteOpt._is_vector_start(dref)
    end
    # _update_variable_info
    @testset "_update_variable_info" begin
        @test isa(InfiniteOpt._update_variable_info(dref, new_info), Nothing)
        @test InfiniteOpt._variable_info(dref) == new_info
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(dref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(dref) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == InfOptConstraintIndex[]
    end
    # _semi_infinite_variable_dependencies
    @testset "_semi_infinite_variable_dependencies" begin
        @test InfiniteOpt._semi_infinite_variable_dependencies(dref) == SemiInfiniteVariableIndex[]
        @test InfiniteOpt._semi_infinite_variable_dependencies(gvref) == SemiInfiniteVariableIndex[]
    end
    # _point_variable_dependencies
    @testset "_point_variable_dependencies" begin
        @test InfiniteOpt._point_variable_dependencies(dref) == PointVariableIndex[]
        @test InfiniteOpt._point_variable_dependencies(gvref) == PointVariableIndex[]
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(dref) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(gvref) == DerivativeIndex[]
    end
    # _derivative_constraint_dependencies
    @testset "_derivative_constraint_dependencies" begin
        @test InfiniteOpt._derivative_constraint_dependencies(dref) == InfOptConstraintIndex[]
        @test InfiniteOpt._derivative_constraint_dependencies(gvref) == InfOptConstraintIndex[]
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(dref) == "var"
        @test name(gvref) == "var"
        @test name(bad_dref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test default
        @test isa(set_name(dref, ""), Nothing)
        @test name(dref) == ""
        @test_throws ErrorException set_name(bad_dref, "")
        # test normal
        @test isa(set_name(gvref, "new"), Nothing)
        @test name(dref) == "new"
    end
    # derivative_argument
    @testset "derivative_argument" begin
        @test derivative_argument(dref) == x
        @test derivative_argument(gvref) == x
    end
    # operator_parameter
    @testset "operator_parameter" begin
        @test operator_parameter(dref) == a
        @test operator_parameter(gvref) == a
    end
    # derivative_method
    @testset "derivative_method" begin
        @test derivative_method(dref) == TestMethod()
        @test derivative_method(gvref) == TestMethod()
    end
    # set_derivative_method
    @testset "set_derivative_method" begin
        @test_throws ErrorException set_derivative_method(dref, TestGenMethod())
        @test_throws ErrorException set_derivative_method(gvref, TestGenMethod())
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        @test raw_parameter_refs(dref) == IC.VectorTuple(a)
        @test raw_parameter_refs(gvref) == IC.VectorTuple(a)
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(dref) == (a, )
        @test parameter_refs(gvref) == (a, )
    end
    # parameter_list
    @testset "parameter_list" begin
        @test parameter_list(dref) == [a]
        @test parameter_list(gvref) == [a]
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, idx) == gvref
    end
    # _var_name_dict
    @testset "_var_name_dict" begin
        @test InfiniteOpt._var_name_dict(m) isa Nothing
    end
    # test has_derivative_constraints
    @testset "has_derivative_constraints" begin 
        @test !has_derivative_constraints(dref)
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(dref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(dref)) == 0
        @test !is_valid(m, dref)
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
    @variable(m, x, Infinite(pref, prefs))
    @variable(m, y, Infinite(pref2))
    num = Float64(0)
    func1 = (a, b) -> 2
    info = VariableInfo(false, num, false, num, false, num, false, NaN, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, func1, false, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, 42, false, true)
    info4 = VariableInfo(true, num, false, num, false, num, false, s -> NaN, false, false)
    # build_derivative
    @testset "build_derivative" begin 
        # test errors 
        @test_throws ErrorException build_derivative(error, info, y, fin)
        @test_throws ErrorException build_derivative(error, info, y, pref)
        @test_throws ErrorException build_derivative(error, info, pref, pref)
        func = (a, b, c) -> a + sum(c)
        bad_info = VariableInfo(true, num, true, num, true, num, true, func, false, false)
        @test_throws ErrorException build_derivative(error, bad_info, x, pref)
        # test for expected output
        @test build_derivative(error, info, x, pref).info isa VariableInfo
        @test build_derivative(error, info, x, pref).is_vector_start
        @test build_derivative(error, info, x, pref).variable_ref == x 
        @test build_derivative(error, info, x, pref).parameter_ref == pref
    end
    # Deriv{V, P}
    @testset "Deriv{V, P}" begin
        @test Deriv(x, pref).argument == x
        @test Deriv(x, pref).operator_parameter == pref
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for each error message
        @test_throws ErrorException build_variable(error, info, Deriv(x, prefs[1]),
                                                   bob = 42)
        @test_throws ErrorException build_variable(error, info, Deriv(0, prefs[1]))
        @test_throws ErrorException build_variable(error, info, Deriv(x, 0))
        # test for expected output
        @test build_variable(error, info, Deriv(x, prefs[1])) isa Derivative
    end
    # add_derivative
    @testset "add_derivative" begin
        # prepare secondary model and parameter and variable
        m2 = InfiniteModel()
        @infinite_parameter(m2, pref3 in [0, 1])
        @variable(m2, z, Infinite(pref3))
        d = build_derivative(error, info, z, pref3)
        # test for error of invalid variable
        @test_throws VariableNotOwned{InfiniteVariableRef} add_derivative(m, d)
        # prepare normal variable
        d = build_derivative(error, info, x, pref)
        # test normal
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_derivative(m, d, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test InfiniteOpt._core_variable_object(dref) == d
        @test InfiniteOpt._derivative_dependencies(pref) == [idx]
        @test InfiniteOpt._derivative_dependencies(x) == [idx]
        @test name(dref) == "name"
        @test InfiniteOpt._existing_derivative_index(x, pref) == idx
        # prepare infinite variable with all the possible info additions
        d = build_derivative(error, info2, x, prefs[1])
        # test info addition functions
        idx = DerivativeIndex(2)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_derivative(m, d, "name") == gvref
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._derivative_dependencies(prefs[1]) == [idx]
        @test InfiniteOpt._derivative_dependencies(x) == [DerivativeIndex(i) for i = 1:2]
        @test !InfiniteOpt._is_vector_start(dref)
        @test InfiniteOpt._variable_info(dref).start == func1
        @test InfiniteOpt._existing_derivative_index(x, prefs[1]) == idx
        # lower bound
        cindex = InfOptConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex)
        @test has_lower_bound(dref)
        @test InfiniteOpt._lower_bound_index(dref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = InfOptConstraintIndex(2)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(dref)
        @test InfiniteOpt._upper_bound_index(dref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = InfOptConstraintIndex(3)
        cref = InfOptConstraintRef(m, cindex)
        @test is_fixed(dref)
        @test InfiniteOpt._fix_index(dref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(dref) == [InfOptConstraintIndex(i)
                                                             for i = 1:3]
        # test redundant build 
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        d = Derivative(info4, true, x, pref)
        @test add_derivative(m, d, "name2") == gvref
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test InfiniteOpt._core_variable_object(dref) == d
        @test name(dref) == "name2"
        @test lower_bound(dref) == 0
    end
    # test JuMP.add_variable 
    @testset "JuMP.add_variable" begin 
        # prepare normal variable
        d = build_derivative(error, info, x, prefs[2])
        # test normal
        idx = DerivativeIndex(3)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, d, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test InfiniteOpt._core_variable_object(dref) == d
        @test InfiniteOpt._derivative_dependencies(prefs[2]) == [idx]
        @test InfiniteOpt._derivative_dependencies(x) == [DerivativeIndex(i) for i = 1:3]
        @test name(dref) == "name"
        @test InfiniteOpt._existing_derivative_index(x, prefs[2]) == idx
    end
    # test "_build_deriv_expr (GeneralVariableRef)"
    @testset "_build_deriv_expr (GeneralVariableRef)" begin
        # test with same parameters 
        @test InfiniteOpt._build_deriv_expr(pref, pref) == 1
        @test InfiniteOpt._build_deriv_expr(prefs[2], prefs[2]) == 1
        # test with new derivative 
        gvref = GeneralVariableRef(m, 4, DerivativeIndex)
        @test InfiniteOpt._build_deriv_expr(y, pref2) == gvref 
        # test returning existing derivative 
        @test InfiniteOpt._build_deriv_expr(y, pref2) == gvref 
        # test others 
        @test InfiniteOpt._build_deriv_expr(pref, pref2) == 0
        @test InfiniteOpt._build_deriv_expr(fin, pref2) == 0
        @test InfiniteOpt._build_deriv_expr(x, pref2) == 0
    end
    # test "_build_deriv_expr (AffExpr)"
    @testset "_build_deriv_expr (AffExpr)" begin
        @test InfiniteOpt._build_deriv_expr(2x+42, pref2) == 0
        @test InfiniteOpt._build_deriv_expr(2pref2 - 23, pref2) == 2
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test InfiniteOpt._build_deriv_expr(-4x, pref) == -4gvref
    end
    # test "_build_deriv_expr (QuadExpr)"
    @testset "_build_deriv_expr (QuadExpr)" begin
        @test InfiniteOpt._build_deriv_expr(2x^2 + 42, pref2) == @expression(m, 0*x)
        @test InfiniteOpt._build_deriv_expr(2pref2^2 -pref2 - 23, pref2) == 4pref2 - 1
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test InfiniteOpt._build_deriv_expr(-4x^2, pref) == -8 * gvref * x
    end
    # test "_build_deriv_expr (Number)"
    @testset "_build_deriv_expr (Number)" begin
        @test InfiniteOpt._build_deriv_expr(42, pref2) == 0
    end
    # test "_build_deriv_expr (Fallback)"
    @testset "_build_deriv_expr (Fallback)" begin
        @test_throws ErrorException InfiniteOpt._build_deriv_expr("bob", pref2)
    end
    # test "_recursive_deriv_build"
    @testset "_recursive_deriv_build" begin
        @test InfiniteOpt._recursive_deriv_build(2pref^2 + pref, [pref, pref]) == 4
        @test InfiniteOpt._recursive_deriv_build(2pref^2 + pref, [pref, pref2]) == 0.
        gvref = GeneralVariableRef(m, 5, DerivativeIndex)
        @test InfiniteOpt._recursive_deriv_build(x + pref, [pref, pref]) == gvref + 0
    end
    # test "deriv"
    @testset "deriv" begin
        # test errors
        @test_throws ErrorException deriv(x, pref, fin)
        @test_throws ErrorException deriv(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test deriv(x, pref) == gvref 
        @test deriv(pref, pref) == 1 
        @test deriv(x^2 + pref, pref) == 2 * gvref * x + 1
        @test deriv(pref^2, pref, pref) == 2
        @test deriv(pref^2, pref, pref, pref) == 0
        @test deriv(42, pref, pref) == 0
    end
    # test "∂"
    @testset "∂" begin
        # test errors
        @test_throws ErrorException ∂(x, pref, fin)
        @test_throws ErrorException ∂(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test ∂(x, pref) == gvref 
        @test ∂(pref, pref) == 1 
        @test ∂(x^2 + pref, pref) == 2 * gvref * x + 1
        @test ∂(pref^2, pref, pref) == 2
        @test ∂(pref^2, pref, pref, pref) == 0
        @test ∂(42, pref, pref) == 0
    end
    # test "@deriv"
    @testset "@deriv" begin
        # test errors
        @test_throws ErrorException @deriv(x, pref, fin)
        @test_throws ErrorException @deriv(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        gvref2 = GeneralVariableRef(m, 5, DerivativeIndex)
        @test @deriv(x, pref) == gvref 
        @test @deriv(pref, pref) == 1 
        @test @deriv(x^2, pref^2) == 2 * gvref2 * x + 2gvref^2
        @test @deriv(pref^2, pref^2) == 2
        @test @deriv(pref^2, pref^3) == 0
        @test @deriv(42, pref^1) == 0
        @test @deriv(x, pref^2, prefs[2]) isa GeneralVariableRef
        # test with semi_infinite variables 
        rv = build_variable(error, gvref, Dict{Int, Float64}(1 => 0.))
        rvref = add_variable(m, rv)
        @test @deriv(rvref, prefs[1]) isa GeneralVariableRef    
        # test with measures 
        meas = Measure(y, TestData(pref, 0, 1), [3], [3, 4], false)
        data = MeasureData(meas, "meas")
        @test InfiniteOpt._add_data_object(m, data) isa MeasureIndex
        mref = GeneralVariableRef(m, 1, MeasureIndex)
        @test @deriv(mref, prefs[2]) isa GeneralVariableRef
    end
    # test "@∂"
    @testset "@∂" begin
        # test errors
        @test_throws ErrorException @∂(x, pref, fin)
        @test_throws ErrorException @∂(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        gvref2 = GeneralVariableRef(m, 5, DerivativeIndex)
        @test @∂(x, pref) == gvref 
        @test @∂(pref, pref) == 1 
        @test @∂(x^2, pref^2) == 2 * gvref2 * x + 2gvref^2
        @test @∂(pref^2, pref^2) == 2
        @test @∂(pref^2, pref^3) == 0
        @test @∂(42, pref^1) == 0
        @test @∂(x, pref^2, prefs[2]) isa GeneralVariableRef
    end
end

# Test the info methods 
@testset "Info Methods" begin 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, x, Infinite(t))
    @variable(m, dx >= 0, Deriv(x, t))
    # test lower bound stuff 
    @testset "Lower Bound Methods" begin 
        @test has_lower_bound(dx)
        @test lower_bound(dx) == 0 
        @test LowerBoundRef(dx) isa InfOptConstraintRef 
        @test set_lower_bound(dx, 2) isa Nothing 
        @test lower_bound(dx) == 2
    end
    # test upper bound stuff 
    @testset "Upper Bound Methods" begin 
        @test !has_upper_bound(dx)
        @test set_upper_bound(dx, 10) isa Nothing 
        @test upper_bound(dx) == 10
        @test UpperBoundRef(dx) isa InfOptConstraintRef 
    end
    # test fix stuff 
    @testset "Fix Methods" begin 
        @test !is_fixed(dx)
        @test fix(dx, 42, force = true) isa Nothing 
        @test fix_value(dx) == 42
        @test FixRef(dx) isa InfOptConstraintRef 
    end
end

# Test Start Methods 
@testset "Start Value" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, inf, Infinite(t))
    @variable(m, d1, Deriv(inf, t), start = 0)
    @variable(m, d2, Deriv(d1, t))
    dref = dispatch_variable_ref(d1)
    # start_value
    @testset "JuMP.start_value" begin
        @test_throws ErrorException start_value(d1)
        @test_throws ErrorException start_value(dref)
    end
    # set_start_value
    @testset "JuMP.set_start_value" begin
        @test_throws ErrorException set_start_value(d1, 0)
        @test_throws ErrorException set_start_value(dref, 0)
    end
    # start_value_function
    @testset "start_value_function" begin
        @test start_value_function(d1)([1]) == 0
        @test start_value_function(dref)([1]) == 0
        @test InfiniteOpt._is_vector_start(dref)
        @test start_value_function(d2) isa Nothing
    end
    # set_start_value_function
    @testset "set_start_value_function" begin
        @test set_start_value_function(d1, 1.5) isa Nothing
        @test start_value_function(d1)([0]) == 1.5
        @test InfiniteOpt._is_vector_start(dref)
        @test !optimizer_model_ready(m)
        func = (a) -> 1
        @test set_start_value_function(d1, func) isa Nothing
        @test start_value_function(d1)(0) == 1
        @test !InfiniteOpt._is_vector_start(dref)
    end
    # reset_start_value_function
    @testset "reset_start_value_function" begin
        @test reset_start_value_function(d1) isa Nothing
        @test start_value_function(d1) isa Nothing
        @test InfiniteOpt._is_vector_start(dref)
        @test !optimizer_model_ready(m)
    end
end

# Test Macro Definition 
@testset "Macro Definition" begin 
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, y, Infinite(t, x))
    @variable(m, q[1:3], Infinite(t, x))
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @variable(m, variable_type = Deriv(y, t)) == gvref
        @test derivative_argument(dref) == y
        @test operator_parameter(dref) == t
        # test anon with changes to fixed
        idx = DerivativeIndex(2)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @variable(m, variable_type = Deriv(y, x[1]), lower_bound = -5, 
                        upper_bound = 10, start = 2) == gvref
        @test derivative_argument(dref) == y
        @test operator_parameter(dref) == x[1]
        @test lower_bound(dref) == -5
        @test upper_bound(dref) == 10
        @test start_value_function(dref) isa Function 
        # test regular with alias
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @variable(m, dx <= 0, Deriv(y, t)) == gvref
        @test derivative_argument(dref) == y
        @test operator_parameter(dref) == t
        @test used_by_derivative(y)
        @test used_by_derivative(t)
        @test upper_bound(dref) == 0
        @test name(dref) == "dx"
    end
    # test array variable definition
    @testset "Array" begin
        # test anon array with one parameter
        idxs = [DerivativeIndex(3), DerivativeIndex(4)]
        drefs = [DerivativeRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @variable(m, [i = 1:2], Deriv(q[i], t), lower_bound = 0) == gvrefs
        @test derivative_argument(drefs[1]) == q[1]
        @test operator_parameter(drefs[2]) == t
        @test lower_bound(drefs[2]) == 0
        # test explicit array 
        idxs = [DerivativeIndex(5), DerivativeIndex(6)]
        drefs = [DerivativeRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @variable(m, q4[i = 1:2] == 3, Deriv(q[i+1], x[i])) == gvrefs
        @test derivative_argument(drefs[1]) == q[2]
        @test operator_parameter(drefs[2]) == x[2]
        @test fix_value(drefs[2]) == 3
        @test name(drefs[1]) == "q4[1]"
        # test semi anon array
        idxs = [DerivativeIndex(7), DerivativeIndex(8)]
        drefs = [DerivativeRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @variable(m, [i = 1:2], Deriv(q[i], x[i]), lower_bound = -5) == gvrefs
        @test derivative_argument(drefs[2]) == q[2]
        @test operator_parameter(drefs[1]) == x[1]
        @test lower_bound(drefs[2]) == -5
        @test name(drefs[1]) == ""
    end
    # test errors
    @testset "Errors" begin
        # test same name error
        @test_macro_throws ErrorException @variable(m, y, Deriv(y, t))
    end
    # test the deprecations 
    @testset "@derivative_variable" begin 
        @test_macro_throws ErrorException @derivative_variable(m, d(y)/d(t))
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, y, Infinite(t, x))
    d = @deriv(y, t)
    vref = dispatch_variable_ref(d)
    # test used_by_semi_infinite_variable
    @testset "used_by_semi_infinite_variable" begin
        @test !used_by_semi_infinite_variable(vref)
        push!(InfiniteOpt._semi_infinite_variable_dependencies(vref),
              SemiInfiniteVariableIndex(1))
        @test used_by_semi_infinite_variable(d)
        @test used_by_semi_infinite_variable(vref)
        empty!(InfiniteOpt._semi_infinite_variable_dependencies(vref))
    end
    # test used_by_point_variable
    @testset "used_by_point_variable" begin
        @test !used_by_point_variable(vref)
        push!(InfiniteOpt._point_variable_dependencies(vref), PointVariableIndex(1))
        @test used_by_point_variable(d)
        @test used_by_point_variable(vref)
        empty!(InfiniteOpt._point_variable_dependencies(vref))
    end
    # test used_by_derivative
    @testset "used_by_derivative" begin
        @test !used_by_derivative(vref)
        push!(InfiniteOpt._derivative_dependencies(vref), DerivativeIndex(1))
        @test used_by_derivative(d)
        @test used_by_derivative(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
    end
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(d)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test used_by_constraint(d)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(d)
        @test !used_by_objective(vref)
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test is_used(d)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by point variable
        num = Float64(0)
        info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
        var = PointVariable(info, d, [0., 0., 0.])
        object = VariableData(var, "var")
        idx = PointVariableIndex(1)
        pvref = PointVariableRef(m, idx)
        @test InfiniteOpt._add_data_object(m, object) == idx
        push!(InfiniteOpt._point_variable_dependencies(vref), idx)
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(pvref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._point_variable_dependencies(vref))
        # test used by semi_infinite variable
        eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1)
        var = SemiInfiniteVariable(d, eval_supps, [2], [2])
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
        info = VariableInfo{Float64, Float64, Float64, Function}(true, num, true,
                                         num, true, num, false, func, true, true)
        d2 = @deriv(d, t)
        dref = DerivativeRef(m, index(d2))
        push!(InfiniteOpt._derivative_dependencies(vref), index(d2))
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(dref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
    end
end

# Test General Queries 
@testset "Other Queries" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [0, 1])
    @variable(m, y, Infinite(t, x))
    d1 = @deriv(y, t)
    d2 = @deriv(y, x[1])
    d3 = @deriv(d2, t)
    cref = @constraint(m, d1 == 0)
    push!(InfiniteOpt._derivative_constraint_dependencies(d1), index(cref))
    InfiniteOpt._set_has_derivative_constraints(t, true)
    # test num_derivatives 
    @testset "num_derivative" begin 
        @test num_derivatives(m) == 3
    end
    # test all_derivatives 
    @testset "all_derivative" begin 
        @test all_derivatives(m) == [d1, d2, d3]
    end
    # test derivative_constraints
    @testset "derivative_constraints" begin 
        @test derivative_constraints(d1) == [cref]
    end
    # test delete_derivative_constraints
    @testset "delete_derivative_constraints" begin 
        @test delete_derivative_constraints(d1) isa Nothing
        @test !has_derivative_constraints(d1)
        @test !has_derivative_constraints(t)
        @test !is_valid(m, cref)
    end
end

# Test the definition of semi-infinite and point derivatives 
@testset "Semi-Infinite/Point Derivatives" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [0, 1])
    @variable(m, y, Infinite(t, x))
    d = @deriv(y, t)
    # test semi-infinite derivatives
    @testset "Semi-infinite Derivatives" begin 
        dref = GeneralVariableRef(m, 1, SemiInfiniteVariableIndex)
        @test @variable(m, variable_type = SemiInfinite(d, 1, x)) == dref
        @test parameter_refs(dref) == (x,)
        @test infinite_variable_ref(dref) == d
        @test eval_supports(dref)[1] == 1
    end
    # test point derivatives
    @testset "Point Derivatives" begin 
        dref = GeneralVariableRef(m, 1, PointVariableIndex)
        @test @variable(m, variable_type = Point(d, 0, [0, 0])) == dref
        @test parameter_refs(dref) == ()
        @test infinite_variable_ref(dref) == d
        @test parameter_values(dref) == (0, [0, 0])
    end
    # test using restrict
    @testset "Restriction Definition" begin 
        # test semi-infinite
        dref = GeneralVariableRef(m, 2, SemiInfiniteVariableIndex)
        @test restrict(d, 0, x) == dref
        @test parameter_refs(dref) == (x,)
        @test infinite_variable_ref(dref) == d
        @test eval_supports(dref)[1] == 0
        # test point 
        dref = GeneralVariableRef(m, 2, PointVariableIndex)
        @test d(1, [0, 0]) == dref
        @test parameter_refs(dref) == ()
        @test infinite_variable_ref(dref) == d
        @test parameter_values(dref) == (1, [0, 0])
    end
end
