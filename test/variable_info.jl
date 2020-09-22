# Test the lower bound methods
@testset "Lower Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @hold_variable(m, var1)
    gvref2 = @hold_variable(m, var2 >= 0)
    gvref3 = @hold_variable(m, var3 == 0)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    # _temp_constraint_ref
    @testset "_temp_constraint_ref" begin
        cindex = ConstraintIndex(1)
        expected = InfOptConstraintRef(m, cindex, ScalarShape())
        @test InfiniteOpt._temp_constraint_ref(m, cindex) == expected
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
        @test has_lower_bound(gvref2)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(vref2) == 0
        @test lower_bound(gvref2) == 0
        @test_throws ErrorException lower_bound(vref1)
    end
    # _lower_bound_index
    @testset "InfiniteOpt._lower_bound_index" begin
        @test InfiniteOpt._lower_bound_index(vref2) == ConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._lower_bound_index(vref1)
    end
    # _set_lower_bound_index
    @testset "_set_lower_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_lower_bound_index(vref2, ConstraintIndex(2)),
                  Nothing)
        @test InfiniteOpt._lower_bound_index(vref2) == ConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_lower_bound_index(vref2, ConstraintIndex(1)),
                  Nothing)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, ConstraintIndex(1), JuMP.ScalarShape())
        # test for finite variable
        @test LowerBoundRef(vref2) == expected
        @test LowerBoundRef(gvref2) == expected
        # prepare infinite variable
        @infinite_parameter(m, t in [0, 1])
        @infinite_variable(m, ivref(t) >= 0)
        expected = InfOptConstraintRef(m, ConstraintIndex(3), JuMP.ScalarShape())
        @test LowerBoundRef(ivref) == expected
    end
    # set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        # test adding new lower bound
        @test isa(set_lower_bound(vref1, 1), Nothing)
        @test has_lower_bound(vref1)
        @test lower_bound(vref1) == 1
        cref = LowerBoundRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref1, MOI.GreaterThan(1.))
        # test changing lower bound
        @test isa(set_lower_bound(vref2, 1.5), Nothing)
        @test has_lower_bound(vref2)
        @test lower_bound(gvref2) == 1.5
        cref = LowerBoundRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref2, MOI.GreaterThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_lower_bound(vref3, 0)
    end
end

# Test the upper bound methods
@testset "Upper Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @hold_variable(m, var1)
    gvref2 = @hold_variable(m, var2 <= 0)
    gvref3 = @hold_variable(m, var3 == 0)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    # is_fixed
    @testset "JuMP.is_fixed" begin
        @test is_fixed(vref3)
        @test !is_fixed(vref2)
    end
    # has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test !has_upper_bound(vref1)
        @test has_upper_bound(vref2)
        @test has_upper_bound(gvref2)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(vref2) == 0
        @test upper_bound(gvref2) == 0
        @test_throws ErrorException upper_bound(vref1)
    end
    # _upper_bound_index
    @testset "InfiniteOpt._upper_bound_index" begin
        @test InfiniteOpt._upper_bound_index(vref2) == ConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._upper_bound_index(vref1)
    end
    # _set_upper_bound_index
    @testset "_set_upper_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_upper_bound_index(vref2, ConstraintIndex(2)),
                  Nothing)
        @test InfiniteOpt._upper_bound_index(vref2) == ConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_upper_bound_index(vref2, ConstraintIndex(1)),
                  Nothing)
    end
    # UpperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, ConstraintIndex(1), JuMP.ScalarShape())
        # test for finite variable
        @test UpperBoundRef(vref2) == expected
        @test UpperBoundRef(gvref2) == expected
        # prepare infinite variable
        @infinite_parameter(m, t in [0, 1])
        @infinite_variable(m, ivref(t) <= 0)
        expected = InfOptConstraintRef(m, ConstraintIndex(3), JuMP.ScalarShape())
        @test UpperBoundRef(ivref) == expected
    end
    # set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        # test adding new upper bound
        @test isa(set_upper_bound(vref1, 1), Nothing)
        @test has_upper_bound(vref1)
        @test upper_bound(vref1) == 1
        cref = UpperBoundRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref1, MOI.LessThan(1.))
        # test changing upper bound
        @test isa(set_upper_bound(vref2, 1.5), Nothing)
        @test has_upper_bound(vref2)
        @test upper_bound(gvref2) == 1.5
        cref = UpperBoundRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref2, MOI.LessThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_upper_bound(vref3, 0)
    end
end

# Test the fix methods
@testset "Fix" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    gvref1 = @hold_variable(m, var1)
    gvref2 = @hold_variable(m, var2 == 0)
    gvref3 = @hold_variable(m, var3)
    gvref4 = @hold_variable(m, var4)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    vref4 = dispatch_variable_ref(gvref4)
    # is_fixed
    @testset "is_fixed" begin
        @test !is_fixed(vref1)
        @test is_fixed(vref2)
        @test is_fixed(gvref2)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(vref2) == 0
        @test fix_value(gvref2) == 0
        @test_throws ErrorException fix_value(vref1)
    end
    # _fix_index
    @testset "InfiniteOpt._fix_index" begin
        @test InfiniteOpt._fix_index(vref2) == ConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._fix_index(vref1)
    end
    # _set_fix_index
    @testset "_set_fix_index" begin
        # test function
        @test isa(InfiniteOpt._set_fix_index(vref2, ConstraintIndex(2)), Nothing)
        @test InfiniteOpt._fix_index(vref2) == ConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_fix_index(vref2, ConstraintIndex(1)), Nothing)
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, ConstraintIndex(1), JuMP.ScalarShape())
        # test for finite variable
        @test FixRef(vref2) == expected
        # prepare infinite variable
        @infinite_parameter(m, t in [0, 1])
        @infinite_variable(m, ivref(t) == 0)
        expected = InfOptConstraintRef(m, ConstraintIndex(2), JuMP.ScalarShape())
        # test for infinite variable
        @test FixRef(ivref) == expected
    end
    # fix
    @testset "JuMP.fix" begin
        # test adding new fix
        @test isa(fix(vref1, 1), Nothing)
        @test is_fixed(vref1)
        @test fix_value(vref1) == 1
        cref = FixRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref1, MOI.EqualTo(1.))
        # test changing fix
        @test isa(fix(gvref2, 1.5), Nothing)
        @test is_fixed(vref2)
        @test fix_value(vref2) == 1.5
        cref = FixRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref2, MOI.EqualTo(1.5))
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
        cref = FixRef(vref3)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref3, MOI.EqualTo(1.))
        # test forcing with upper
        @test isa(fix(gvref4, 1.5, force = true), Nothing)
        @test is_fixed(vref4)
        @test fix_value(vref4) == 1.5
        cref = FixRef(vref4)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref4, MOI.EqualTo(1.5))
    end
end

# Test the start value methods
@testset "Start Value" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    gvref = @hold_variable(m, var1, start = 0)
    vref = dispatch_variable_ref(gvref)
    gvref2 = @hold_variable(m, var2)
    vref2 = dispatch_variable_ref(gvref2)
    @infinite_parameter(m, t in [0, 1])
    @infinite_variable(m, inf(t), start = 0)
    dinf = dispatch_variable_ref(inf)
    @infinite_variable(m, inf2(t))
    dinf2 = dispatch_variable_ref(inf2)
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(vref) == 0
        @test start_value(gvref) == 0
        @test start_value(vref2) isa Nothing
        @test start_value(gvref2) isa Nothing
        @test_throws ErrorException start_value(inf)
        @test_throws ErrorException start_value(dinf)
    end
    # set_start_value
    @testset "JuMP.set_start_value" begin
        @test isa(set_start_value(vref, 1.5), Nothing)
        @test start_value(vref) == 1.5
        @test !optimizer_model_ready(m)
        @test isa(set_start_value(gvref, 1), Nothing)
        @test_throws ErrorException set_start_value(inf, 0)
        @test_throws ErrorException set_start_value(dinf, 0)
    end
    # start_value_function
    @testset "start_value_function" begin
        @test start_value_function(inf)([1]) == 0
        @test start_value_function(dinf)([1]) == 0
        @test InfiniteOpt._is_vector_start(dinf)
        @test start_value_function(inf2) isa Nothing
        @test start_value_function(dinf2) isa Nothing
    end
    # set_start_value_function
    @testset "set_start_value_function" begin
        @test set_start_value_function(inf, 1.5) isa Nothing
        @test start_value_function(inf)([0]) == 1.5
        @test InfiniteOpt._is_vector_start(dinf)
        @test !optimizer_model_ready(m)
        func = (a) -> 1
        @test set_start_value_function(inf, func) isa Nothing
        @test start_value_function(inf)(0) == 1
        @test !InfiniteOpt._is_vector_start(dinf)
    end
    # reset_start_value_function
    @testset "reset_start_value_function" begin
        @test reset_start_value_function(inf) isa Nothing
        @test start_value_function(inf) isa Nothing
        @test InfiniteOpt._is_vector_start(dinf)
        @test !optimizer_model_ready(m)
    end
end

# Test the binary methods
@testset "Binary" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @hold_variable(m, var1)
    gvref2 = @hold_variable(m, var2, Bin)
    gvref3 = @hold_variable(m, var3, Int)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    # is_binary
    @testset "JuMP.is_binary" begin
        @test !is_binary(vref1)
        @test is_binary(vref2)
        @test is_binary(gvref2)
    end
    # _binary_index
    @testset "InfiniteOpt._binary_index" begin
        @test InfiniteOpt._binary_index(vref2) == ConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._binary_index(vref1)
    end
    # _set_binary_index
    @testset "_set_binary_index" begin
        # test function
        @test isa(InfiniteOpt._set_binary_index(vref2, ConstraintIndex(2)), Nothing)
        @test InfiniteOpt._binary_index(vref2) == ConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_binary_index(vref2, ConstraintIndex(1)), Nothing)
    end
    # BinaryRef
    @testset "BinaryRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, ConstraintIndex(1), JuMP.ScalarShape())
        # test for finite variable
        @test BinaryRef(vref2) == expected
        # prepare infinite variable
        @infinite_parameter(m, t in [0, 1])
        @infinite_variable(m, ivref(t), Bin)
        expected = InfOptConstraintRef(m, ConstraintIndex(3), JuMP.ScalarShape())
        # test for infinite variable
        @test BinaryRef(ivref) == expected
    end
    # set_binary
    @testset "JuMP.set_binary" begin
        # test adding binary
        @test isa(set_binary(vref1), Nothing)
        @test is_binary(vref1)
        cref = BinaryRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref1, MOI.ZeroOne())
        # test setting binary again
        @test isa(set_binary(gvref2), Nothing)
        @test is_binary(vref2)
        cref = BinaryRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref2, MOI.ZeroOne())
        # test integer variable error
        @test_throws ErrorException set_binary(vref3)
    end
end

# Test the integer methods
@testset "Integer" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @hold_variable(m, var1)
    gvref2 = @hold_variable(m, var2, Int)
    gvref3 = @hold_variable(m, var3, Bin)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    # is_integer
    @testset "JuMP.is_integer" begin
        @test !is_integer(vref1)
        @test is_integer(vref2)
        @test is_integer(gvref2)
    end
    # _integer_index
    @testset "InfiniteOpt._integer_index" begin
        @test InfiniteOpt._integer_index(vref2) == ConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._integer_index(vref1)
    end
    # _set_integer_index
    @testset "_set_integer_index" begin
        # test function
        @test isa(InfiniteOpt._set_integer_index(vref2, ConstraintIndex(2)), Nothing)
        @test InfiniteOpt._integer_index(vref2) == ConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_integer_index(vref2, ConstraintIndex(1)), Nothing)
    end
    # IntegerRef
    @testset "IntegerRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, ConstraintIndex(1), JuMP.ScalarShape())
        # test for finite variable
        @test IntegerRef(vref2) == expected
        # prepare infinite variable
        @infinite_parameter(m, t in [0, 1])
        @infinite_variable(m, ivref(t), Int)
        expected = InfOptConstraintRef(m, ConstraintIndex(3), JuMP.ScalarShape())
        # test for infinite variable
        @test IntegerRef(ivref) == expected
    end
    # set_integer
    @testset "JuMP.set_integer" begin
        # test adding integer
        @test isa(set_integer(vref1), Nothing)
        @test is_integer(vref1)
        cref = IntegerRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref1, MOI.Integer())
        # test setting integer again
        @test isa(set_integer(gvref2), Nothing)
        @test is_integer(vref2)
        cref = IntegerRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._core_constraint_object(cref) == ScalarConstraint(gvref2, MOI.Integer())
        # test integer variable error
        @test_throws ErrorException set_integer(vref3)
    end
end

# Test relaxation 
@testset "JuMP.relax_integrality" begin 
    # setup the model 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_variable(m, y(t), Int)
    @hold_variable(m, 0.25 <= x <= 1.25, Bin)
    @hold_variable(m, z, Bin)
    # test relaxing it 
    result_store = []
    @test push!(result_store, relax_integrality(m)) isa Vector
    @test !is_integer(y)
    @test !is_binary(x)
    @test !is_binary(z)
    @test lower_bound(x) == 0.25
    @test upper_bound(x) == 1
    @test lower_bound(z) == 0
    @test upper_bound(z) == 1
    # test unrelaxing it 
    @test result_store[1]() isa Nothing 
    @test is_integer(y)
    @test is_binary(x)
    @test is_binary(z)
    @test lower_bound(x) == 0.25
    @test upper_bound(x) == 1.25
    @test !has_lower_bound(z)
    @test !has_upper_bound(z)
end
