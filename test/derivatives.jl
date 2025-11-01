# Test the basics 
@testset "Basics" begin
    # setup data
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1], derivative_method = TestMethod())
    @variable(m, x, Infinite(a))
    idx = DerivativeIndex(1)
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, num, false, false)
    deriv = Derivative(info, x, a, 1)
    deriv2 = Derivative(new_info, x, a, 2)
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
        @test isequal(dispatch_variable_ref(m, idx), dref)
        @test isequal(dispatch_variable_ref(gvref), dref)
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
        @test InfiniteOpt._existing_derivative_index(x, a, 1) === nothing
        m.deriv_lookup[(x, a, 1)] = idx
        @test InfiniteOpt._existing_derivative_index(x, a, 1) === idx
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test core_object(dref) === deriv
        @test core_object(gvref) === deriv
    end
    # _set_core_object
    @testset "_set_core_object" begin
        @test InfiniteOpt._set_core_object(dref, deriv) isa Nothing
        @test InfiniteOpt._set_core_object(dref, deriv2) isa Nothing
        @test InfiniteOpt._set_core_object(dref, deriv) isa Nothing
    end
    # parameter_group_int_indices
    @testset "parameter_group_int_indices" begin
        @test InfiniteOpt.parameter_group_int_indices(dref) == [1]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(dref) == [1]
    end
    # test _variable_info
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(dref) == info
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
        @test isequal(derivative_argument(dref), x)
        @test isequal(derivative_argument(gvref), x)
    end
    # operator_parameter
    @testset "operator_parameter" begin
        @test isequal(operator_parameter(dref), a)
        @test isequal(operator_parameter(gvref), a)
    end
    # derivative_order
    @testset "derivative_order" begin
        @test isequal(derivative_order(dref), 1)
        @test isequal(derivative_order(gvref), 1)
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
        @test isequal(raw_parameter_refs(dref), IC.VectorTuple(a))
        @test isequal(raw_parameter_refs(gvref), IC.VectorTuple(a))
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test isequal(parameter_refs(dref), (a, ))
        @test isequal(parameter_refs(gvref), (a, ))
    end
    # parameter_list
    @testset "parameter_list" begin
        @test isequal(parameter_list(dref), [a])
        @test isequal(parameter_list(gvref), [a])
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
    @variable(m, w, Infinite(pref, pref2))
    f = parameter_function(sin, pref)
    num = Float64(0)
    func1 = (a, b) -> 2
    info = VariableInfo(false, num, false, num, false, num, false, NaN, false, false)
    info2 = VariableInfo(true, func1, true, num, true, func1, true, func1, false, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, 42, false, true)
    info4 = VariableInfo(true, num, false, num, false, num, true, num, false, false)
    info5 = VariableInfo(false, num, false, num, false, num, false, NaN, false, false)
    # build_derivative
    @testset "build_derivative" begin 
        # test errors 
        @test_throws ErrorException build_derivative(error, info, y, fin)
        @test_throws ErrorException build_derivative(error, info, y, pref)
        @test_throws ErrorException build_derivative(error, info, pref, pref)
        @test_throws ErrorException build_derivative(error, info, y, pref2, -1)
        @test_throws ErrorException build_derivative(error, info, x, prefs[1])
        @test_throws ErrorException build_derivative(error, info, f, pref)
        @test_throws ErrorException build_derivative(error, info3, x, pref)
        func = (a, b, c) -> a + sum(c)
        bad_info = VariableInfo(true, num, true, num, true, num, true, func, false, false)
        @test_throws ErrorException build_derivative(error, bad_info, x, pref)
        # test for expected output
        @test build_derivative(error, info, x, pref).info isa VariableInfo
        @test build_derivative(error, info, x, pref).order == 1
        @test isequal(build_derivative(error, info, x, pref).variable_ref, x)
        @test isequal(build_derivative(error, info, x, pref).parameter_ref, pref)
        @test build_derivative(error, info2, x, pref, 2).order == 2
    end
    # Deriv{V, P}
    @testset "Deriv{V, P}" begin
        @test isequal(Deriv(x, pref).argument, x)
        @test isequal(Deriv(x, pref).operator_parameter, pref)
        @test isequal(Deriv(x, pref).order, 1)
        @test isequal(Deriv(x, pref, 3).order, 3)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test for each error message
        @test_throws ErrorException build_variable(error, info, Deriv(x, pref),
                                                   bob = 42)
        @test_throws ErrorException build_variable(error, info, Deriv(0, pref))
        @test_throws ErrorException build_variable(error, info, Deriv(x, 0))
        @test_throws ErrorException build_variable(error, info, Deriv(x, pref, 0))
        # test for expected output
        @test build_variable(error, info, Deriv(x, pref)) isa Derivative
        @test build_variable(error, info, Deriv(x, pref, 3)) isa Derivative
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
        gvref = GeneralVariableRef(m, idx)
        @test isequal(add_derivative(m, d, "name"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test core_object(dref) == d
        @test InfiniteOpt._derivative_dependencies(pref) == [idx]
        @test InfiniteOpt._derivative_dependencies(x) == [idx]
        @test name(dref) == "name"
        @test InfiniteOpt._existing_derivative_index(x, pref, 1) == idx
        # prepare infinite variable with all the possible info additions
        d = build_derivative(error, info2, x, pref, 2)
        # test info addition functions
        idx2 = DerivativeIndex(2)
        dref = DerivativeRef(m, idx2)
        gvref = GeneralVariableRef(m, idx2)
        @test isequal(add_derivative(m, d, "name"), gvref)
        @test !transformation_backend_ready(m)
        @test InfiniteOpt._derivative_dependencies(pref) == [idx, idx2]
        @test InfiniteOpt._derivative_dependencies(x) == [DerivativeIndex(i) for i = 1:2]
        @test InfiniteOpt._variable_info(dref).start.func == func1
        @test InfiniteOpt._existing_derivative_index(x, pref, 2) == idx2
        # lower bound
        cindex = InfOptConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex)
        @test has_lower_bound(dref)
        @test InfiniteOpt._lower_bound_index(dref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                           <:MOI.GreaterThan{<:ParameterFunctionRealWrapper}}
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
                                                           <:MOI.EqualTo{<:ParameterFunctionRealWrapper}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(dref) == [InfOptConstraintIndex(i)
                                                             for i = 1:3]
        # test redundant build 
        idx = DerivativeIndex(2)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        d = Derivative(InfiniteOpt._DefaultDerivativeInfo, x, pref, 2)
        @test isequal(add_derivative(m, d, "name2"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test name(dref) == "name2"
        @test lower_bound(dref) isa ParameterFunction
        # test redundant build with info changes
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        d = Derivative(info4, x, pref, 1)
        @test add_derivative(m, d, "name3") == gvref
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test core_object(dref) == d
        @test name(dref) == "name3"
        @test lower_bound(dref) == 0
        @test !has_upper_bound(dref)
        @test start_value(dref) == 0
    end
    # test JuMP.add_variable 
    @testset "JuMP.add_variable" begin 
        # prepare normal variable
        d = build_derivative(error, info, x, pref, 4)
        # test normal
        idx = DerivativeIndex(3)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        @test isequal(add_variable(m, d, "name"), gvref)
        @test haskey(InfiniteOpt._data_dictionary(dref), idx)
        @test core_object(dref) == d
        @test InfiniteOpt._derivative_dependencies(pref) == [DerivativeIndex(i) for i in 1:3]
        @test InfiniteOpt._derivative_dependencies(x) == [DerivativeIndex(i) for i = 1:3]
        @test name(dref) == "name"
        @test InfiniteOpt._existing_derivative_index(x, pref, 4) == idx
    end
    # test "_build_deriv_expr (GeneralVariableRef)"
    @testset "_build_deriv_expr (GeneralVariableRef)" begin
        # test zero order 
        @test InfiniteOpt._build_deriv_expr(pref, pref, 0) == pref
        # test with same parameters 
        @test InfiniteOpt._build_deriv_expr(pref, pref, 1) == 1
        @test InfiniteOpt._build_deriv_expr(pref, pref, 2) == 0
        # test with new derivative 
        gvref = GeneralVariableRef(m, 4, DerivativeIndex)
        @test isequal(InfiniteOpt._build_deriv_expr(y, pref2, 1), gvref)
        # test returning existing derivative 
        @test isequal(InfiniteOpt._build_deriv_expr(y, pref2, 1), gvref)
        # test others 
        @test InfiniteOpt._build_deriv_expr(pref, pref2, 1) == 0
        @test InfiniteOpt._build_deriv_expr(fin, pref2, 1) == 0
        @test InfiniteOpt._build_deriv_expr(x, pref2, 1) == 0
    end
    # test "_build_deriv_expr (AffExpr)"
    @testset "_build_deriv_expr (AffExpr)" begin
        @test InfiniteOpt._build_deriv_expr(2x+42, pref2, 1) == 0
        @test InfiniteOpt._build_deriv_expr(2pref2 - 23, pref2, 1) == 2
        @test InfiniteOpt._build_deriv_expr(2pref2 - 23, pref2, 2) == 0
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test isequal(InfiniteOpt._build_deriv_expr(-4x, pref, 1), -4gvref)
    end
    # test "_build_deriv_expr (QuadExpr)"
    @testset "_build_deriv_expr (QuadExpr)" begin
        @test isequal(InfiniteOpt._build_deriv_expr(2x^2 + 42, pref2, 1), @expression(m, 0*x))
        @test isequal(InfiniteOpt._build_deriv_expr(2pref2^2 -pref2 - 23, pref2, 1), 4pref2 - 1)
        @test isequal(InfiniteOpt._build_deriv_expr(2pref2^2 -pref2 - 23, pref2, 2), 4)
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test isequal(InfiniteOpt._build_deriv_expr(-4x^2, pref, 1).terms, (-8 * x * gvref).terms)
        gvref2 = GeneralVariableRef(m, 2, DerivativeIndex)
        @test isequal(InfiniteOpt._build_deriv_expr(-4x^2, pref, 2).terms, (-8 * gvref2 * x - 8 * gvref^2).terms)
    end
    # test "_build_deriv_expr (Number)"
    @testset "_build_deriv_expr (Number)" begin
        @test InfiniteOpt._build_deriv_expr(42, pref2, 1) == 0
    end
    # test "_build_deriv_expr (Fallback)"
    @testset "_build_deriv_expr (Fallback)" begin
        @test_throws ErrorException InfiniteOpt._build_deriv_expr("bob", pref2, 2)
    end
    # test "_recursive_deriv_build"
    @testset "_recursive_deriv_build" begin
        @test InfiniteOpt._recursive_deriv_build(2pref^2 + pref, [pref], [2]) == 4
        @test InfiniteOpt._recursive_deriv_build(2pref^2 + pref, [pref, pref2], [1, 1]) == 0.
        gvref = GeneralVariableRef(m, 2, DerivativeIndex)
        @test isequal(InfiniteOpt._recursive_deriv_build(x + pref, [pref], [2]), gvref + 0)
    end
    # test "deriv"
    @testset "deriv" begin
        # test errors
        @test_throws ErrorException deriv(x, pref, fin)
        @test_throws ErrorException deriv(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test isequal(deriv(x, pref), gvref)
        @test deriv(pref, pref) == 1 
        @test isequal_canonical(deriv(x^2 + pref, pref), 2 * gvref * x + 1)
        @test deriv(pref^2, pref, pref) == 2
        @test deriv(pref^2, pref, pref, pref2) == 0
        @test deriv(pref^2, pref, pref2) == 0
        @test deriv(pref^2, pref, pref, pref) == 0
        @test deriv(42, pref, pref) == 0
        @test deriv(deriv(x, pref), pref) == deriv(x, pref, pref) # test unnesting
    end
    # test "∂"
    @testset "∂" begin
        # test errors
        @test_throws ErrorException ∂(x, pref, fin)
        @test_throws ErrorException ∂(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        @test isequal(∂(x, pref), gvref)
        @test ∂(pref, pref) == 1 
        @test isequal_canonical(∂(x^2 + pref, pref), 2 * gvref * x + 1)
        @test ∂(pref^2, pref, pref) == 2
        @test ∂(pref^2, pref, pref, pref) == 0
        @test ∂(42, pref, pref) == 0
    end
    # test "@deriv"
    @testset "@deriv" begin
        # test errors
        @test_throws ErrorException @deriv(x, pref, fin)
        @test_throws ErrorException @deriv(x)
        @test_macro_throws ErrorException @deriv(x, pref, bad = 1)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        gvref2 = GeneralVariableRef(m, 2, DerivativeIndex)
        @test isequal(@deriv(x, pref), gvref)
        @test @deriv(pref, pref) == 1 
        @test isequal_canonical(@deriv(x^2, pref^2), 2 * gvref2 * x + 2gvref^2)
        @test @deriv(pref^2, pref^2) == 2
        @test @deriv(pref^2, pref^3) == 0
        @test @deriv(42, pref^1) == 0
        @test @deriv(w, pref^2, pref2) isa GeneralVariableRef
        # test with semi_infinite variables 
        rvref = gvref(pref, [1, 0])
        @test @deriv(rvref, pref) isa GeneralVariableRef    
        # test with measures 
        meas = Measure(w, TestData(pref, 0, 1), [2], [2], false)
        data = MeasureData(meas, "meas")
        @test InfiniteOpt._add_data_object(m, data) isa MeasureIndex
        mref = GeneralVariableRef(m, 1, MeasureIndex)
        @test @deriv(mref, pref2) isa GeneralVariableRef
    end
    # test "@∂"
    @testset "@∂" begin
        # test errors
        @test_throws ErrorException @∂(x, pref, fin)
        @test_throws ErrorException @∂(x)
        # test normal 
        gvref = GeneralVariableRef(m, 1, DerivativeIndex)
        gvref2 = GeneralVariableRef(m, 2, DerivativeIndex)
        @test isequal(@∂(x, pref), gvref)
        @test @∂(pref, pref) == 1 
        @test isequal_canonical(@∂(x^2, pref^2), 2 * gvref2 * x + 2gvref^2)
        @test @∂(pref^2, pref^2) == 2
        @test @∂(pref^2, pref^3) == 0
        @test @∂(42, pref^1) == 0
        @test @∂(w, pref^2, pref2) isa GeneralVariableRef
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
        @test set_lower_bound(dx, sin) isa Nothing 
        @test lower_bound(dx)(0.5) == sin(0.5)
    end
    # test upper bound stuff 
    @testset "Upper Bound Methods" begin 
        @test !has_upper_bound(dx)
        @test set_upper_bound(dx, cos) isa Nothing 
        @test upper_bound(dx)(0.5) == cos(0.5)
        @test UpperBoundRef(dx) isa InfOptConstraintRef 
    end
    # test fix stuff 
    @testset "Fix Methods" begin 
        @test !is_fixed(dx)
        @test fix(dx, 42, force = true) isa Nothing 
        @test fix_value(dx) == 42
        @test FixRef(dx) isa InfOptConstraintRef 
    end
    # test integality errors
    @testset "Integrality Errors" begin 
        @test_throws ErrorException set_binary(dx)
        @test_throws ErrorException set_integer(dx)
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
        @test start_value(d1) == 0
        @test start_value(dref) == 0
        @test start_value(d2) isa Nothing
        @test_deprecated start_value_function(d1) == 0
    end
    # set_start_value
    @testset "JuMP.set_start_value" begin
        @test set_start_value(d1, sin)  isa Nothing
        @test start_value(d1)(0.5) == sin(0.5)
        @test set_start_value(d2, 42)  isa Nothing
        @test start_value(d2) == 42
        @test_deprecated set_start_value_function(d1, cos) isa Nothing
    end
end

# Test Macro Definition 
@testset "Macro Definition" begin 
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1], independent = true)
    @variable(m, y, Infinite(t, x...))
    @variable(m, q[1:3], Infinite(t, x...))
    # test single variable definition
    @testset "Single" begin
        # test simple anon case
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Deriv(y, t)), gvref)
        @test isequal(derivative_argument(dref), y)
        @test isequal(operator_parameter(dref), t)
        @test isequal(derivative_order(dref), 1)
        # test anon with changes to fixed
        idx = DerivativeIndex(2)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        @test isequal(@variable(m, variable_type = Deriv(y, x[1], 3), lower_bound = -5, 
                        upper_bound = 10, start = 2), gvref)
        @test isequal(derivative_argument(dref), y)
        @test isequal(operator_parameter(dref), x[1])
        @test isequal(derivative_order(dref), 3)
        @test lower_bound(dref) == -5
        @test upper_bound(dref) == 10
        @test start_value(dref) == 2
        # test regular with alias
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        gvref = GeneralVariableRef(m, idx)
        @test isequal(@variable(m, dx <= 0, Deriv(y, t)), gvref)
        @test isequal(derivative_argument(dref), y)
        @test isequal(operator_parameter(dref), t)
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
        gvrefs = [GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [i = 1:2], Deriv(q[i], t), lower_bound = 0), gvrefs)
        @test isequal(derivative_argument(drefs[1]), q[1])
        @test isequal(operator_parameter(drefs[2]), t)
        @test lower_bound(drefs[2]) == 0
        # test explicit array 
        idxs = [DerivativeIndex(5), DerivativeIndex(6)]
        drefs = [DerivativeRef(m, idx) for idx in idxs]
        gvrefs = [GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, q4[i = 1:2] == 3, Deriv(q[i+1], x[i])), gvrefs)
        @test isequal(derivative_argument(drefs[1]), q[2])
        @test isequal(operator_parameter(drefs[2]), x[2])
        @test fix_value(drefs[2]) == 3
        @test name(drefs[1]) == "q4[1]"
        # test semi anon array
        idxs = [DerivativeIndex(7), DerivativeIndex(8)]
        drefs = [DerivativeRef(m, idx) for idx in idxs]
        gvrefs = [GeneralVariableRef(m, idx) for idx in idxs]
        @test isequal(@variable(m, [i = 1:2], Deriv(q[i], x[i]), lower_bound = (x...) -> 3), gvrefs)
        @test isequal(derivative_argument(drefs[2]), q[2])
        @test isequal(operator_parameter(drefs[1]), x[1])
        @test lower_bound(drefs[2])(0, [0, 0]) == 3
        @test name(drefs[1]) == ""
    end
    # test errors
    @testset "Errors" begin
        # test same name error
        @test_macro_throws ErrorException @variable(m, y, Deriv(y, t))
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1], independent = true)
    @variable(m, y, Infinite(t, x...))
    d = deriv(y, t)
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
        pvref = dispatch_variable_ref(d(0, 0, 0))
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(pvref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._point_variable_dependencies(vref))
        # test used by semi_infinite variable
        rvref = dispatch_variable_ref(d(0.5, x[1], 1.0))
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(rvref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._semi_infinite_variable_dependencies(vref))
        # test used by derivative
        d2 = deriv(d, x[1])
        dref = DerivativeRef(m, index(d2))
        @test !isempty(InfiniteOpt._derivative_dependencies(vref))
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
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    @variable(m, y, Infinite(t, x...))
    d1 = deriv(y, t)
    d2 = deriv(y, x[1])
    d3 = deriv(d2, t)
    cref = @constraint(m, d1 == 0)
    push!(InfiniteOpt._derivative_constraint_dependencies(d1), index(cref))
    InfiniteOpt._set_has_derivative_constraints(t, true)
    # test num_derivatives 
    @testset "num_derivative" begin 
        @test num_derivatives(m) == 3
    end
    # test all_derivatives 
    @testset "all_derivative" begin 
        @test isequal(all_derivatives(m), [d1, d2, d3])
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
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    @variable(m, y, Infinite(t, x...))
    d = deriv(y, t)
    # test semi-infinite derivatives
    @testset "Semi-infinite Derivatives" begin 
        dref = GeneralVariableRef(m, 1, SemiInfiniteVariableIndex)
        @test isequal(@variable(m, variable_type = SemiInfinite(d, 1, x...)), dref)
        @test isequal(parameter_refs(dref), (x...,))
        @test isequal(infinite_variable_ref(dref), d)
        @test eval_support(dref)[1] == 1
    end
    # test point derivatives
    @testset "Point Derivatives" begin 
        dref = GeneralVariableRef(m, 1, PointVariableIndex)
        @test isequal(@variable(m, variable_type = Point(d, 0, 0, 0)), dref)
        @test parameter_refs(dref) == ()
        @test isequal(infinite_variable_ref(dref), d)
        @test parameter_values(dref) == (0, 0, 0)
    end
    # test using restrict
    @testset "Restriction Definition" begin 
        # test semi-infinite
        dref = GeneralVariableRef(m, 2, SemiInfiniteVariableIndex)
        @test isequal(restrict(d, 0, x...), dref)
        @test isequal(parameter_refs(dref), (x...,))
        @test isequal(infinite_variable_ref(dref), d)
        @test eval_support(dref)[1] == 0
        # test point 
        dref = GeneralVariableRef(m, 2, PointVariableIndex)
        @test isequal(d(1, 0, 0), dref)
        @test parameter_refs(dref) == ()
        @test isequal(infinite_variable_ref(dref), d)
        @test parameter_values(dref) == (1, 0, 0)
    end
end
