# Test internal info methods
@testset "Internal Info Get/Set" begin
    # initialize model and infinite variable
    m = InfiniteModel()
    info = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, 0., true, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param)
    var = InfiniteVariable(info, (pref,))
    m.vars[1] = var
    vref = InfiniteVariableRef(m, 1)
    bounds = ParameterBounds()
    # _variable_info
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(vref) == info
    end
    # _update_variable_info (Infinite)
    @testset "_update_variable_info (Infinite)" begin
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
    # _update_variable_info (Point)
    @testset "_update_variable_info (Point)" begin
        # initialize point variable
        var = PointVariable(info, vref, (0,))
        m.vars[2] = var
        vref = PointVariableRef(m, 2)
        # test function
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
    # _update_variable_info (Hold)
    @testset "_update_variable_info (Hold)" begin
        # initialize hold variable
        var = HoldVariable(info, bounds)
        m.vars[3] = var
        vref = HoldVariableRef(m, 3)
        # test function
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
end

# Test the lower bound methods
@testset "Lower Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    info2 = VariableInfo(true, 0., false, 0., false, 0., false, 0., false, false)
    info3 = VariableInfo(false, 0., false, 0., true, 0., false, 0., false, false)
    bounds = ParameterBounds()
    var1 = HoldVariable(info1, bounds)
    var2 = HoldVariable(info2, bounds)
    var3 = HoldVariable(info3, bounds)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_lower_bound[2] = 1
    vref1 = HoldVariableRef(m, 1)
    vref2 = HoldVariableRef(m, 2)
    vref3 = HoldVariableRef(m, 3)
    # _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # test normal with no previous constrs
        @test isa(InfiniteOpt._update_var_constr_mapping([vref2], 1), Nothing)
        @test m.var_to_constrs[JuMP.index(vref2)] == [1]
        # test normal with previous constrs
        @test isa(InfiniteOpt._update_var_constr_mapping([vref2], 2), Nothing)
        @test m.var_to_constrs[JuMP.index(vref2)] == [1, 2]
        # undo changes
        delete!(m.var_to_constrs, JuMP.index(vref2))
    end
    # set_name
    @testset "JuMP.set_name" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.GreaterThan(0.))
        cref = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test isa(set_name(cref, "bob"), Nothing)
        @test m.constr_to_name[1] == "bob"
        # undo changes
        delete!(m.constr_to_name, 1)
    end
    # index
    @testset "JuMP.index" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.GreaterThan(0.))
        cref = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test JuMP.index(cref) == 1
    end
    # add_constraint
    @testset "JuMP.add_constraint" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.GreaterThan(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test add_constraint(m, constr, "bob") == expected
        @test m.constrs[1] == constr
        @test m.constr_to_name[1] == "bob"
        @test !m.constr_in_var_info[1]
        m.constr_in_var_info[1] = true
    end
    # is_fixed
    @testset "JuMP.is_fixed" begin
        @test is_fixed(vref3)
        @test !is_fixed(vref2)
    end
    # has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test !has_lower_bound(vref1)
        @test has_lower_bound(vref2)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(vref2) == 0
        @test_throws ErrorException lower_bound(vref1)
    end
    # _lower_bound_index
    @testset "JuMP._lower_bound_index" begin
        @test JuMP._lower_bound_index(vref2) == 1
        @test_throws ErrorException JuMP._lower_bound_index(vref1)
    end
    # set_lower_bound_index
    @testset "JuMP.set_lower_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_lower_bound_index(vref2, 2), Nothing)
        @test m.var_to_lower_bound[2] == 2
        # undo changes
        m.var_to_lower_bound[2] = 1
    end
    # set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        # test adding new lower bound
        @test isa(set_lower_bound(vref1, 1), Nothing)
        @test has_lower_bound(vref1)
        @test lower_bound(vref1) == 1
        @test m.constr_in_var_info[2]
        @test !optimizer_model_ready(m)
        @test m.constrs[2] == ScalarConstraint(vref1, MOI.GreaterThan(1.))
        # test changing lower bound
        @test isa(set_lower_bound(vref2, 1.5), Nothing)
        @test has_lower_bound(vref2)
        @test lower_bound(vref2) == 1.5
        @test m.constr_in_var_info[1]
        @test !optimizer_model_ready(m)
        @test m.constrs[1] == ScalarConstraint(vref2, MOI.GreaterThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_lower_bound(vref3, 0)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        # prepare constraint reference
        m.constr_to_name[1] = ""
        constr = ScalarConstraint(vref2, MOI.GreaterThan(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test for finite variable
        @test LowerBoundRef(vref2) == expected
        # prepare infinite variable
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref = add_parameter(m, param)
        var = InfiniteVariable(info1, (pref,))
        m.vars[4] = var
        m.var_to_name[4] = "var"
        vref = InfiniteVariableRef(m, 4)
        # prepare lower bound constraint
        set_lower_bound(vref, 0)
        constr = ScalarConstraint(vref, MOI.GreaterThan(0.))
        expected = InfiniteConstraintRef(m, 3, JuMP.shape(constr))
        # test for infinite variable
        @test LowerBoundRef(vref) == expected
    end
end

# Test the upper bound methods
@testset "Upper Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    info2 = VariableInfo(false, 0., true, 0., false, 0., false, 0., false, false)
    info3 = VariableInfo(false, 0., false, 0., true, 0., false, 0., false, false)
    bounds = ParameterBounds()
    var1 = HoldVariable(info1, bounds)
    var2 = HoldVariable(info2, bounds)
    var3 = HoldVariable(info3, bounds)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_upper_bound[2] = 1
    vref1 = HoldVariableRef(m, 1)
    vref2 = HoldVariableRef(m, 2)
    vref3 = HoldVariableRef(m, 3)
    # add_constraint
    @testset "JuMP.add_constraint" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.LessThan(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test add_constraint(m, constr, "bob") == expected
        @test m.constrs[1] == constr
        @test m.constr_to_name[1] == "bob"
        @test !m.constr_in_var_info[1]
        m.constr_in_var_info[1] = true
    end
    # has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test !has_upper_bound(vref1)
        @test has_upper_bound(vref2)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(vref2) == 0
        @test_throws ErrorException upper_bound(vref1)
    end
    # _upper_bound_index
    @testset "JuMP._upper_bound_index" begin
        @test JuMP._upper_bound_index(vref2) == 1
        @test_throws ErrorException JuMP._upper_bound_index(vref1)
    end
    # _set_upper_bound_index
    @testset "_set_upper_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_upper_bound_index(vref2, 2), Nothing)
        @test m.var_to_upper_bound[2] == 2
        # undo changes
        m.var_to_upper_bound[2] = 1
    end
    # set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        # test adding new upper bound
        @test isa(set_upper_bound(vref1, 1), Nothing)
        @test has_upper_bound(vref1)
        @test upper_bound(vref1) == 1
        @test m.constr_in_var_info[2]
        @test !optimizer_model_ready(m)
        @test m.constrs[2] == ScalarConstraint(vref1, MOI.LessThan(1.))
        # test changing upper bound
        @test isa(set_upper_bound(vref2, 1.5), Nothing)
        @test has_upper_bound(vref2)
        @test upper_bound(vref2) == 1.5
        @test m.constr_in_var_info[1]
        @test !optimizer_model_ready(m)
        @test m.constrs[1] == ScalarConstraint(vref2, MOI.LessThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_upper_bound(vref3, 0)
    end
    # upperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        # prepare constraint reference
        m.constr_to_name[1] = ""
        constr = ScalarConstraint(vref2, MOI.LessThan(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test for finite variable
        @test UpperBoundRef(vref2) == expected
        # prepare infinite variable
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref = add_parameter(m, param)
        var = InfiniteVariable(info1, (pref,))
        m.vars[4] = var
        m.var_to_name[4] = "var"
        vref = InfiniteVariableRef(m, 4)
        # prepare Upper bound constraint
        set_upper_bound(vref, 0)
        constr = ScalarConstraint(vref, MOI.LessThan(0.))
        expected = InfiniteConstraintRef(m, 3, JuMP.shape(constr))
        # test for infinite variable
        @test UpperBoundRef(vref) == expected
    end
end

# Test the fix methods
@testset "Fix" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    info2 = VariableInfo(false, 0., false, 0., true, 0., false, 0., false, false)
    bounds = ParameterBounds()
    var1 = HoldVariable(info1, bounds)
    var2 = HoldVariable(info2, bounds)
    var3 = HoldVariable(info1, bounds)
    var4 = HoldVariable(info1, bounds)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.vars[4] = var4
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_name[4] = "var4"
    m.var_to_fix[2] = 1
    vref1 = HoldVariableRef(m, 1)
    vref2 = HoldVariableRef(m, 2)
    vref3 = HoldVariableRef(m, 3)
    vref4 = HoldVariableRef(m, 4)
    # add_constraint
    @testset "JuMP.add_constraint" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.EqualTo(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test add_constraint(m, constr, "bob") == expected
        @test m.constrs[1] == constr
        @test m.constr_to_name[1] == "bob"
        @test !m.constr_in_var_info[1]
        m.constr_in_var_info[1] = true
    end
    # is_fixed
    @testset "is_fixed" begin
        @test !is_fixed(vref1)
        @test is_fixed(vref2)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(vref2) == 0
        @test_throws ErrorException fix_value(vref1)
    end
    # _fix_index
    @testset "JuMP._fix_index" begin
        @test JuMP._fix_index(vref2) == 1
        @test_throws ErrorException JuMP._fix_index(vref1)
    end
    # _set_fix_index
    @testset "_set_fix_index" begin
        # test function
        @test isa(InfiniteOpt._set_fix_index(vref2, 2), Nothing)
        @test m.var_to_fix[2] == 2
        # undo changes
        m.var_to_fix[2] = 1
    end
    # fix
    @testset "JuMP.fix" begin
        # test adding new fix
        @test isa(fix(vref1, 1), Nothing)
        @test is_fixed(vref1)
        @test fix_value(vref1) == 1
        @test m.constr_in_var_info[2]
        @test !optimizer_model_ready(m)
        @test m.constrs[2] == ScalarConstraint(vref1, MOI.EqualTo(1.))
        # test changing fix
        @test isa(fix(vref2, 1.5), Nothing)
        @test is_fixed(vref2)
        @test fix_value(vref2) == 1.5
        @test m.constr_in_var_info[1]
        @test !optimizer_model_ready(m)
        @test m.constrs[1] == ScalarConstraint(vref2, MOI.EqualTo(1.5))
        # add lower and upper bounds to vars 3 and 4
        set_lower_bound(vref3, 0.0)
        set_upper_bound(vref4, 0.0)
        # test lower/upper bound  error
        @test_throws ErrorException fix(vref3, 0)
        @test_throws ErrorException fix(vref4, 0)
        # test forcing with lower
        @test isa(fix(vref3, 1, force = true), Nothing)
        @test is_fixed(vref3)
        @test fix_value(vref3) == 1
        @test m.constr_in_var_info[5]
        @test !optimizer_model_ready(m)
        @test m.constrs[5] == ScalarConstraint(vref3, MOI.EqualTo(1.))
        # test forcing with upper
        @test isa(fix(vref4, 1.5, force = true), Nothing)
        @test is_fixed(vref4)
        @test fix_value(vref4) == 1.5
        @test m.constr_in_var_info[6]
        @test !optimizer_model_ready(m)
        @test m.constrs[6] == ScalarConstraint(vref4, MOI.EqualTo(1.5))
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        # prepare constraint reference
        m.constr_to_name[1] = ""
        constr = ScalarConstraint(vref2, MOI.EqualTo(0.))
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test for finite variable
        @test FixRef(vref2) == expected
        # prepare infinite variable
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref = add_parameter(m, param)
        var = InfiniteVariable(info1, (pref,))
        m.vars[5] = var
        m.var_to_name[5] = "var"
        vref = InfiniteVariableRef(m, 5)
        # prepare Upper bound constraint
        fix(vref, 0)
        constr = ScalarConstraint(vref, MOI.EqualTo(0.))
        expected = InfiniteConstraintRef(m, 7, JuMP.shape(constr))
        # test for infinite variable
        @test FixRef(vref) == expected
    end
end

# Test the start value methods
@testset "Start Value" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., true, 0., false, false)
    bounds = ParameterBounds()
    var = HoldVariable(info1, bounds)
    m.vars[1] = var
    vref = HoldVariableRef(m, 1)
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(vref) == 0
    end
    # set_start_value
    @testset "JuMP.set_start_value" begin
        @test isa(set_start_value(vref, 1.5), Nothing)
        @test start_value(vref) == 1.5
        @test !optimizer_model_ready(m)
    end
end

# Test the binary methods
@testset "Binary" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    info2 = VariableInfo(false, 0., false, 0., false, 0., false, 0., true, false)
    info3 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, true)
    bounds = ParameterBounds()
    var1 = HoldVariable(info1, bounds)
    var2 = HoldVariable(info2, bounds)
    var3 = HoldVariable(info3, bounds)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_zero_one[2] = 1
    vref1 = HoldVariableRef(m, 1)
    vref2 = HoldVariableRef(m, 2)
    vref3 = HoldVariableRef(m, 3)
    # add_constraint
    @testset "JuMP.add_constraint" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.ZeroOne())
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test add_constraint(m, constr, "bob") == expected
        @test m.constrs[1] == constr
        @test m.constr_to_name[1] == "bob"
        @test !m.constr_in_var_info[1]
        m.constr_in_var_info[1] = true
    end
    # is_binary
    @testset "JuMP.is_binary" begin
        @test !is_binary(vref1)
        @test is_binary(vref2)
    end
    # _binary_index
    @testset "JuMP._binary_index" begin
        @test JuMP._binary_index(vref2) == 1
        @test_throws ErrorException JuMP._binary_index(vref1)
    end
    # _set_binary_index
    @testset "_set_binary_index" begin
        # test function
        @test isa(InfiniteOpt._set_binary_index(vref2, 2), Nothing)
        @test m.var_to_zero_one[2] == 2
        # undo changes
        m.var_to_zero_one[2] = 1
    end
    # set_binary
    @testset "JuMP.set_binary" begin
        # test adding binary
        @test isa(set_binary(vref1), Nothing)
        @test is_binary(vref1)
        @test m.constr_in_var_info[2]
        @test !optimizer_model_ready(m)
        @test m.constrs[2] == ScalarConstraint(vref1, MOI.ZeroOne())
        # test setting binary again
        @test isa(set_binary(vref2), Nothing)
        @test is_binary(vref2)
        @test m.constr_in_var_info[1]
        @test !optimizer_model_ready(m)
        @test m.constrs[1] == ScalarConstraint(vref2, MOI.ZeroOne())
        # test integer variable error
        @test_throws ErrorException set_binary(vref3)
    end
    # BinaryRef
    @testset "BinaryRef" begin
        # prepare constraint reference
        m.constr_to_name[1] = ""
        constr = ScalarConstraint(vref2, MOI.ZeroOne())
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test for finite variable
        @test BinaryRef(vref2) == expected
        # prepare infinite variable
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref = add_parameter(m, param)
        var = InfiniteVariable(info1, (pref,))
        m.vars[4] = var
        m.var_to_name[4] = "var"
        vref = InfiniteVariableRef(m, 4)
        # prepare Upper bound constraint
        set_binary(vref)
        constr = ScalarConstraint(vref, MOI.ZeroOne())
        expected = InfiniteConstraintRef(m, 3, JuMP.shape(constr))
        # test for infinite variable
        @test BinaryRef(vref) == expected
    end
end

# Test the integer methods
@testset "Integer" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    info1 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    info2 = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, true)
    info3 = VariableInfo(false, 0., false, 0., false, 0., false, 0., true, false)
    bounds = ParameterBounds()
    var1 = HoldVariable(info1, bounds)
    var2 = HoldVariable(info2, bounds)
    var3 = HoldVariable(info3, bounds)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_integrality[2] = 1
    vref1 = HoldVariableRef(m, 1)
    vref2 = HoldVariableRef(m, 2)
    vref3 = HoldVariableRef(m, 3)
    # add_constraint
    @testset "JuMP.add_constraint" begin
        # prepare constraint reference
        constr = ScalarConstraint(vref2, MOI.Integer())
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test function
        @test add_constraint(m, constr, "bob") == expected
        @test m.constrs[1] == constr
        @test m.constr_to_name[1] == "bob"
        @test !m.constr_in_var_info[1]
        m.constr_in_var_info[1] = true
    end
    # is_integer
    @testset "JuMP.is_integer" begin
        @test !is_integer(vref1)
        @test is_integer(vref2)
    end
    # _integer_index
    @testset "JuMP._integer_index" begin
        @test JuMP._integer_index(vref2) == 1
        @test_throws ErrorException JuMP._integer_index(vref1)
    end
    # _set_integer_index
    @testset "_set_integer_index" begin
        # test function
        @test isa(InfiniteOpt._set_integer_index(vref2, 2), Nothing)
        @test m.var_to_integrality[2] == 2
        # undo changes
        m.var_to_integrality[2] = 1
    end
    # set_integer
    @testset "JuMP.set_integer" begin
        # test adding integer
        @test isa(set_integer(vref1), Nothing)
        @test is_integer(vref1)
        @test m.constr_in_var_info[2]
        @test !optimizer_model_ready(m)
        @test m.constrs[2] == ScalarConstraint(vref1, MOI.Integer())
        # test setting integer again
        @test isa(set_integer(vref2), Nothing)
        @test is_integer(vref2)
        @test m.constr_in_var_info[1]
        @test !optimizer_model_ready(m)
        @test m.constrs[1] == ScalarConstraint(vref2, MOI.Integer())
        # test integer variable error
        @test_throws ErrorException set_integer(vref3)
    end
    # IntegerRef
    @testset "IntegerRef" begin
        # prepare constraint reference
        m.constr_to_name[1] = ""
        constr = ScalarConstraint(vref2, MOI.Integer())
        expected = FiniteConstraintRef(m, 1, JuMP.shape(constr))
        # test for finite variable
        @test IntegerRef(vref2) == expected
        # prepare infinite variable
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref = add_parameter(m, param)
        var = InfiniteVariable(info1, (pref,))
        m.vars[4] = var
        m.var_to_name[4] = "var"
        vref = InfiniteVariableRef(m, 4)
        # prepare Upper bound constraint
        set_integer(vref)
        constr = ScalarConstraint(vref, MOI.Integer())
        expected = InfiniteConstraintRef(m, 3, JuMP.shape(constr))
        # test for infinite variable
        @test IntegerRef(vref) == expected
    end
end
