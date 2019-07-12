# Test extensions to basic Base methods
@testset "Base Extensions" begin
    # initialize models and references
    m = InfiniteModel()
    m2 = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    pvref = PointVariableRef(m, 2)
    gvref = GlobalVariableRef(m, 3)
    pref = ParameterRef(m, 1)
    # variable compare
    @testset "(==)" begin
        @test ivref == ivref
        @test pvref == pvref
        @test gvref == gvref
        @test ivref == InfiniteVariableRef(m, 1)
        @test pvref == PointVariableRef(m, 2)
        @test gvref == GlobalVariableRef(m, 3)
        @test !(ivref == InfiniteVariableRef(m, 2))
        @test !(ivref == InfiniteVariableRef(m2, 1))
        @test !(ivref != InfiniteVariableRef(m, 1))
        @test !(pref == ivref)
    end
    # copy(v)
    @testset "copy(v)" begin
        @test copy(ivref) == ivref
        @test copy(pvref) == pvref
        @test copy(gvref) == gvref
    end
    # copy(v, m)
    @testset "copy(v, m)" begin
        @test copy(ivref, m2) == InfiniteVariableRef(m2, 1)
        @test copy(pvref, m2) == PointVariableRef(m2, 2)
        @test copy(gvref, m2) == GlobalVariableRef(m2, 3)
    end
    # broadcastable
    @testset "broadcastable" begin
        @test isa(Base.broadcastable(ivref), Base.RefValue{InfiniteVariableRef})
        @test isa(Base.broadcastable(pvref), Base.RefValue{PointVariableRef})
        @test isa(Base.broadcastable(gvref), Base.RefValue{GlobalVariableRef})
    end
end

# Test core JuMP methods
@testset "Core JuMP Extensions" begin
    # initialize models and references
    m = InfiniteModel()
    m2 = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    pvref = PointVariableRef(m, 2)
    gvref = GlobalVariableRef(m, 3)
    pref = ParameterRef(m, 1)
    # isequal_canonical
    @testset "JuMP.isequal_canonical" begin
        @test isequal_canonical(ivref, ivref)
        @test isequal_canonical(pvref, pvref)
        @test isequal_canonical(gvref, gvref)
        @test !isequal_canonical(ivref, InfiniteVariableRef(m2, 1))
        @test !isequal_canonical(ivref, InfiniteVariableRef(m, 2))
    end
    # variable_type(m)
    @testset "JuMP.variable_type(m)" begin
        @test variable_type(m) == GeneralVariableRef
    end
    # variable_type(m, t)
    @testset "JuMP.variable_type(m, t)" begin
        @test variable_type(m, Infinite) == InfiniteVariableRef
        @test variable_type(m, Point) == PointVariableRef
        @test variable_type(m, Global) == GlobalVariableRef
        @test variable_type(m, Parameter) == ParameterRef
        @test_throws ErrorException variable_type(m, :bad)
    end
end

# Test precursor functions needed for add_parameter
@testset "Basic Reference Queries" begin
    # initialize model and infinite variable
    m = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    m.vars[1] = InfiniteVariable(info, (pref, ))
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(ivref) == 1
    end
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(ivref) == m
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, ivref)
        @test !is_valid(InfiniteModel(), ivref)
        @test !is_valid(m, InfiniteVariableRef(m, 5))
    end
end

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
    # _update_variable_info (Global)
    @testset "_update_variable_info (Global)" begin
        # initialize global variable
        var = GlobalVariable(info)
        m.vars[3] = var
        vref = GlobalVariableRef(m, 3)
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
    var1 = GlobalVariable(info1)
    var2 = GlobalVariable(info2)
    var3 = GlobalVariable(info3)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_lower_bound[2] = 1
    vref1 = GlobalVariableRef(m, 1)
    vref2 = GlobalVariableRef(m, 2)
    vref3 = GlobalVariableRef(m, 3)
    # _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # test normal with no previous constrs
        @test isa(InfiniteOpt._update_var_constr_mapping([vref2], 1), Nothing)
        @test m.var_to_constrs[index(vref2)] == [1]
        # test normal with previous constrs
        @test isa(InfiniteOpt._update_var_constr_mapping([vref2], 2), Nothing)
        @test m.var_to_constrs[index(vref2)] == [1, 2]
        # undo changes
        delete!(m.var_to_constrs, index(vref2))
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
        @test index(cref) == 1
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
    # lower_bound_index
    @testset "JuMP.lower_bound_index" begin
        @test lower_bound_index(vref2) == 1
        @test_throws ErrorException lower_bound_index(vref1)
    end
    # set_lower_bound_index
    @testset "JuMP.set_lower_bound_index" begin
        # test function
        @test isa(set_lower_bound_index(vref2, 2), Nothing)
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
    var1 = GlobalVariable(info1)
    var2 = GlobalVariable(info2)
    var3 = GlobalVariable(info3)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_upper_bound[2] = 1
    vref1 = GlobalVariableRef(m, 1)
    vref2 = GlobalVariableRef(m, 2)
    vref3 = GlobalVariableRef(m, 3)
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
    # upper_bound_index
    @testset "JuMP.upper_bound_index" begin
        @test upper_bound_index(vref2) == 1
        @test_throws ErrorException upper_bound_index(vref1)
    end
    # set_upper_bound_index
    @testset "JuMP.set_upper_bound_index" begin
        # test function
        @test isa(set_upper_bound_index(vref2, 2), Nothing)
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
    var1 = GlobalVariable(info1)
    var2 = GlobalVariable(info2)
    var3 = GlobalVariable(info1)
    var4 = GlobalVariable(info1)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.vars[4] = var4
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_name[4] = "var4"
    m.var_to_fix[2] = 1
    vref1 = GlobalVariableRef(m, 1)
    vref2 = GlobalVariableRef(m, 2)
    vref3 = GlobalVariableRef(m, 3)
    vref4 = GlobalVariableRef(m, 4)
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
    # fix_index
    @testset "JuMP.fix_index" begin
        @test fix_index(vref2) == 1
        @test_throws ErrorException fix_index(vref1)
    end
    # set_fix_index
    @testset "JuMP.set_fix_index" begin
        # test function
        @test isa(set_fix_index(vref2, 2), Nothing)
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
    var = GlobalVariable(info1)
    m.vars[1] = var
    vref = GlobalVariableRef(m, 1)
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
    var1 = GlobalVariable(info1)
    var2 = GlobalVariable(info2)
    var3 = GlobalVariable(info3)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_zero_one[2] = 1
    vref1 = GlobalVariableRef(m, 1)
    vref2 = GlobalVariableRef(m, 2)
    vref3 = GlobalVariableRef(m, 3)
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
    # binary_index
    @testset "JuMP.binary_index" begin
        @test binary_index(vref2) == 1
        @test_throws ErrorException binary_index(vref1)
    end
    # set_binary_index
    @testset "JuMP.set_binary_index" begin
        # test function
        @test isa(set_binary_index(vref2, 2), Nothing)
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
    var1 = GlobalVariable(info1)
    var2 = GlobalVariable(info2)
    var3 = GlobalVariable(info3)
    m.vars[1] = var1
    m.vars[2] = var2
    m.vars[3] = var3
    m.var_to_name[1] = "var1"
    m.var_to_name[2] = "var2"
    m.var_to_name[3] = "var3"
    m.var_to_integrality[2] = 1
    vref1 = GlobalVariableRef(m, 1)
    vref2 = GlobalVariableRef(m, 2)
    vref3 = GlobalVariableRef(m, 3)
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
    # integer_index
    @testset "JuMP.integer_index" begin
        @test integer_index(vref2) == 1
        @test_throws ErrorException integer_index(vref1)
    end
    # set_integer_index
    @testset "JuMP.set_integer_index" begin
        # test function
        @test isa(set_integer_index(vref2, 2), Nothing)
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
    # add_variable
    # TODO test info updates
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
    info = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
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
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 0]))
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
    # add_variable
    # TODO test info updates
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
        @test supports(pref) == [0, 0.5]
        @test supports(pref2) == [1]
        @test m.var_to_name[4] == "name"
    end
end

# Test name methods
@testset "Global Variable Name" begin
    # initialize model and variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    var = GlobalVariable(info)
    m.vars[1] = var
    m.var_to_name[1] = "test"
    vref = GlobalVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
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
@testset "Global Variable Definition" begin
    # initialize model and info
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    # build_variable
    @testset "JuMP.build_variable" begin
        expected = GlobalVariable(info)
        @test build_variable(error, info, Global) == expected
    end
    # add_variable
    # TODO test info updates
    @testset "JuMP.add_variable" begin
        v = build_variable(error, info, Global)
        @test add_variable(m, v, "name") == GlobalVariableRef(m, 1)
        @test haskey(m.vars, 1)
        @test m.var_to_name[1] == "name"
    end
end
