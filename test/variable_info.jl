# Test the lower bound methods
@testset "Lower Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @variable(m, var1)
    gvref2 = @variable(m, var2 >= 0)
    gvref3 = @variable(m, var3 == 0)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    f(a, b) = 42
    @variable(m, y >= f, Infinite(t, x))
    y0 = y(0, x)
    yb = y(0, -1)
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
        @test has_lower_bound(y)
        @test has_lower_bound(y0)
        @test has_lower_bound(yb)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(vref2) == 0
        @test lower_bound(gvref2) == 0
        @test_throws ErrorException lower_bound(vref1)
        @test lower_bound(y)(0, 0) == 42
        @test lower_bound(y0)(0) == 42
        @test lower_bound(yb) == 42
    end
    # _lower_bound_index
    @testset "InfiniteOpt._lower_bound_index" begin
        @test InfiniteOpt._lower_bound_index(vref2) == InfOptConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._lower_bound_index(vref1)
    end
    # _set_lower_bound_index
    @testset "_set_lower_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_lower_bound_index(vref2, InfOptConstraintIndex(2)),
                  Nothing)
        @test InfiniteOpt._lower_bound_index(vref2) == InfOptConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_lower_bound_index(vref2, InfOptConstraintIndex(1)),
                  Nothing)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(1))
        # test
        @test LowerBoundRef(vref2) == expected
        @test LowerBoundRef(gvref2) == expected
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(3))
        @test LowerBoundRef(y) == expected
        @test LowerBoundRef(y0) == expected
        @test LowerBoundRef(yb) == expected
    end
    # set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        # test adding new lower bound
        @test isa(set_lower_bound(vref1, 1), Nothing)
        @test has_lower_bound(vref1)
        @test lower_bound(vref1) == 1
        cref = LowerBoundRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref1, MOI.GreaterThan(1.))
        # test changing lower bound
        @test isa(set_lower_bound(vref2, 1.5), Nothing)
        @test has_lower_bound(vref2)
        @test lower_bound(gvref2) == 1.5
        cref = LowerBoundRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref2, MOI.GreaterThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_lower_bound(vref3, 0)
        # test infinite variable lower bound
        @test isa(set_lower_bound(y, 1), Nothing)
        @test lower_bound(y) == 1
        @test lower_bound(y0) == 1
        @test lower_bound(yb) == 1
        @test isa(set_lower_bound(y, f), Nothing)
        @test lower_bound(y)(0, 0) == 42
        # test restricted variables
        @test set_lower_bound(y0, 2) isa Nothing
        @test lower_bound(y0) == 2
        @test set_lower_bound(y0, 3) isa Nothing
        @test lower_bound(y0) == 3
        @test fix(yb, 4, force = true) isa Nothing
        @test_throws AssertionError set_lower_bound(yb, 3)
        @test delete_lower_bound(y0) isa Nothing
        @test !has_lower_bound(y0)
        @test_throws ErrorException lower_bound(y0)
    end
end

# Test the upper bound methods
@testset "Upper Bound" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @variable(m, var1)
    gvref2 = @variable(m, var2 <= 0)
    gvref3 = @variable(m, var3 == 0)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    f(a, b) = 42
    @variable(m, y <= f, Infinite(t, x))
    y0 = y(0, x)
    yb = y(0, -1)
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
        @test has_upper_bound(y)
        @test has_upper_bound(y0)
        @test has_upper_bound(yb)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(vref2) == 0
        @test upper_bound(gvref2) == 0
        @test_throws ErrorException upper_bound(vref1)
        @test upper_bound(y)(0, 0) == 42
        @test upper_bound(y0)(0) == 42
        @test upper_bound(yb) == 42
    end
    # _upper_bound_index
    @testset "InfiniteOpt._upper_bound_index" begin
        @test InfiniteOpt._upper_bound_index(vref2) == InfOptConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._upper_bound_index(vref1)
    end
    # _set_upper_bound_index
    @testset "_set_upper_bound_index" begin
        # test function
        @test isa(InfiniteOpt._set_upper_bound_index(vref2, InfOptConstraintIndex(2)),
                  Nothing)
        @test InfiniteOpt._upper_bound_index(vref2) == InfOptConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_upper_bound_index(vref2, InfOptConstraintIndex(1)),
                  Nothing)
    end
    # UpperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(1))
        # test
        @test UpperBoundRef(vref2) == expected
        @test UpperBoundRef(gvref2) == expected
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(3))
        @test UpperBoundRef(y) == expected
        @test UpperBoundRef(y0) == expected
        @test UpperBoundRef(yb) == expected
    end
    # set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        # test adding new upper bound
        @test isa(set_upper_bound(vref1, 1), Nothing)
        @test has_upper_bound(vref1)
        @test upper_bound(vref1) == 1
        cref = UpperBoundRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref1, MOI.LessThan(1.))
        @test_throws ErrorException set_upper_bound(vref1, sin)
        # test changing upper bound
        @test isa(set_upper_bound(vref2, 1.5), Nothing)
        @test has_upper_bound(vref2)
        @test upper_bound(gvref2) == 1.5
        cref = UpperBoundRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref2, MOI.LessThan(1.5))
        # test fixed variable error
        @test_throws AssertionError set_upper_bound(vref3, 0)
        # test infinite variable lower bound
        @test isa(set_upper_bound(y, 1), Nothing)
        @test upper_bound(y) == 1
        @test upper_bound(y0) == 1
        @test upper_bound(yb) == 1
        @test isa(set_upper_bound(y, f), Nothing)
        @test upper_bound(y)(0, 0) == 42
        # test restricted variables
        @test set_upper_bound(y0, 2) isa Nothing
        @test upper_bound(y0) == 2
        @test set_upper_bound(y0, 3) isa Nothing
        @test upper_bound(y0) == 3
        @test fix(yb, 4, force = true) isa Nothing
        @test_throws AssertionError set_upper_bound(yb, 3)
        @test delete_upper_bound(y0) isa Nothing
        @test !has_upper_bound(y0)
        @test_throws ErrorException upper_bound(y0)
    end
end

# Test the fix methods
@testset "Fix" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    gvref1 = @variable(m, var1)
    gvref2 = @variable(m, var2 == 0)
    gvref3 = @variable(m, var3)
    gvref4 = @variable(m, var4)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    vref4 = dispatch_variable_ref(gvref4)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    f(a, b) = 42
    @variable(m, y == f, Infinite(t, x))
    y0 = y(0, x)
    yb = y(0, -1)
    # is_fixed
    @testset "is_fixed" begin
        @test !is_fixed(vref1)
        @test is_fixed(vref2)
        @test is_fixed(gvref2)
        @test is_fixed(y)
        @test is_fixed(y0)
        @test is_fixed(yb)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(vref2) == 0
        @test fix_value(gvref2) == 0
        @test_throws ErrorException fix_value(vref1)
        @test fix_value(y)(0, 0) == 42
        @test fix_value(y0)(0) == 42
        @test fix_value(yb) == 42
    end
    # _fix_index
    @testset "InfiniteOpt._fix_index" begin
        @test InfiniteOpt._fix_index(vref2) == InfOptConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._fix_index(vref1)
    end
    # _set_fix_index
    @testset "_set_fix_index" begin
        # test function
        @test isa(InfiniteOpt._set_fix_index(vref2, InfOptConstraintIndex(2)), Nothing)
        @test InfiniteOpt._fix_index(vref2) == InfOptConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_fix_index(vref2, InfOptConstraintIndex(1)), Nothing)
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(1))
        # test
        @test FixRef(vref2) == expected
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(2))
        @test FixRef(y) == expected
        @test FixRef(y0) == expected
        @test FixRef(yb) == expected
    end
    # fix
    @testset "JuMP.fix" begin
        # test adding new fix
        @test isa(fix(vref1, 1), Nothing)
        @test is_fixed(vref1)
        @test fix_value(vref1) == 1
        cref = FixRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref1, MOI.EqualTo(1.))
        # test changing fix
        @test isa(fix(gvref2, 1.5), Nothing)
        @test is_fixed(vref2)
        @test fix_value(vref2) == 1.5
        cref = FixRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref2, MOI.EqualTo(1.5))
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
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref3, MOI.EqualTo(1.))
        # test forcing with upper
        @test isa(fix(gvref4, 1.5, force = true), Nothing)
        @test is_fixed(vref4)
        @test fix_value(vref4) == 1.5
        cref = FixRef(vref4)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref4, MOI.EqualTo(1.5))
        # test restricted and infinite variables
        @test unfix(y) isa Nothing
        @test !is_fixed(y)
        @test set_lower_bound(y, 2) isa Nothing
        @test set_upper_bound(y, 4) isa Nothing
        @test_throws ErrorException fix(y, 3)
        @test_throws ErrorException fix(y0, 5)
        @test fix(y, f, force = true) isa Nothing
        @test fix_value(y)(0, 0) == 42
        @test fix(y0, 6) isa Nothing
        @test fix_value(y0) == 6
        @test fix(y0, 3) isa Nothing
        @test fix_value(y0) == 3
        @test fix(yb, 2) isa Nothing
        @test fix_value(yb) == 2
        @test unfix(yb) isa Nothing
        @test !is_fixed(yb)
        @test_throws ErrorException fix_value(yb)
    end
end

# Test the start value methods
@testset "Start Value" begin
    # initialize model and 4 test variables
    m = InfiniteModel()
    gvref = @variable(m, var1, start = 0)
    vref = dispatch_variable_ref(gvref)
    gvref2 = @variable(m, var2)
    vref2 = dispatch_variable_ref(gvref2)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    f(a, b) = 42
    @variable(m, y, Infinite(t, x), start = f)
    y0 = y(0, x)
    yb = y(0, -1)
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(vref) == 0
        @test start_value(gvref) == 0
        @test start_value(vref2) isa Nothing
        @test start_value(gvref2) isa Nothing
        @test start_value(y)(0, 0) == 42
        @test start_value(y0)(0) == 42
        @test start_value(yb) == 42
    end
    # set_start_value
    @testset "JuMP.set_start_value" begin
        @test isa(set_start_value(vref, 1.5), Nothing)
        @test start_value(vref) == 1.5
        @test !transformation_backend_ready(m)
        @test isa(set_start_value(gvref, 1), Nothing)
        @test set_start_value(y, 2) isa Nothing
        @test start_value(y) == 2
        @test start_value(y0) == 2
        @test start_value(yb) == 2
        @test isa(set_start_value(y, f), Nothing)
        @test start_value(y)(0, 0) == 42
        @test isa(set_start_value(y0, 3), Nothing)
        @test start_value(y0) == 3
        @test isa(set_start_value(yb, 4), Nothing)
        @test start_value(yb) == 4
    end
end

# Test the binary methods
@testset "Binary" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @variable(m, var1)
    gvref2 = @variable(m, var2, Bin)
    gvref3 = @variable(m, var3, Int)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    @variable(m, y, Infinite(t, x), Bin)
    y0 = y(0, x)
    yb = y(0, -1)
    # is_binary
    @testset "JuMP.is_binary" begin
        @test !is_binary(vref1)
        @test is_binary(vref2)
        @test is_binary(gvref2)
        @test is_binary(y)
        @test is_binary(y0)
        @test is_binary(yb)
    end
    # _binary_index
    @testset "InfiniteOpt._binary_index" begin
        @test InfiniteOpt._binary_index(vref2) == InfOptConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._binary_index(vref1)
    end
    # _set_binary_index
    @testset "_set_binary_index" begin
        # test function
        @test isa(InfiniteOpt._set_binary_index(vref2, InfOptConstraintIndex(2)), Nothing)
        @test InfiniteOpt._binary_index(vref2) == InfOptConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_binary_index(vref2, InfOptConstraintIndex(1)), Nothing)
    end
    # BinaryRef
    @testset "BinaryRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(1))
        # test
        @test BinaryRef(vref2) == expected
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(3))
        @test BinaryRef(y) == expected
        @test BinaryRef(y0) == expected
        @test BinaryRef(yb) == expected
    end
    # set_binary
    @testset "JuMP.set_binary" begin
        # test adding binary
        @test isa(set_binary(vref1), Nothing)
        @test is_binary(vref1)
        cref = BinaryRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref1, MOI.ZeroOne())
        # test setting binary again
        @test isa(set_binary(gvref2), Nothing)
        @test is_binary(vref2)
        cref = BinaryRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref2, MOI.ZeroOne())
        # test integer variable error
        @test_throws ErrorException set_binary(vref3)
        @test_throws ErrorException set_binary(y0)
        @test_throws ErrorException set_binary(yb)
    end
end

# Test the integer methods
@testset "Integer" begin
    # initialize model and 3 test variables
    m = InfiniteModel()
    gvref1 = @variable(m, var1)
    gvref2 = @variable(m, var2, Int)
    gvref3 = @variable(m, var3, Bin)
    vref1 = dispatch_variable_ref(gvref1)
    vref2 = dispatch_variable_ref(gvref2)
    vref3 = dispatch_variable_ref(gvref3)
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    @variable(m, y, Infinite(t, x), Int)
    y0 = y(0, x)
    yb = y(0, -1)
    # is_integer
    @testset "JuMP.is_integer" begin
        @test !is_integer(vref1)
        @test is_integer(vref2)
        @test is_integer(gvref2)
        @test is_integer(y)
        @test is_integer(y0)
        @test is_integer(yb)
    end
    # _integer_index
    @testset "InfiniteOpt._integer_index" begin
        @test InfiniteOpt._integer_index(vref2) == InfOptConstraintIndex(1)
        @test_throws ErrorException InfiniteOpt._integer_index(vref1)
    end
    # _set_integer_index
    @testset "_set_integer_index" begin
        # test function
        @test isa(InfiniteOpt._set_integer_index(vref2, InfOptConstraintIndex(2)), Nothing)
        @test InfiniteOpt._integer_index(vref2) == InfOptConstraintIndex(2)
        # undo changes
        @test isa(InfiniteOpt._set_integer_index(vref2, InfOptConstraintIndex(1)), Nothing)
    end
    # IntegerRef
    @testset "IntegerRef" begin
        # prepare constraint reference
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(1))
        # test
        @test IntegerRef(vref2) == expected
        expected = InfOptConstraintRef(m, InfOptConstraintIndex(3))
        @test IntegerRef(y) == expected
        @test IntegerRef(y0) == expected
        @test IntegerRef(yb) == expected
    end
    # set_integer
    @testset "JuMP.set_integer" begin
        # test adding integer
        @test isa(set_integer(vref1), Nothing)
        @test is_integer(vref1)
        cref = IntegerRef(vref1)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref1, MOI.Integer())
        # test setting integer again
        @test isa(set_integer(gvref2), Nothing)
        @test is_integer(vref2)
        cref = IntegerRef(vref2)
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test !transformation_backend_ready(m)
        @test constraint_object(cref) == ScalarConstraint(gvref2, MOI.Integer())
        # test integer variable error
        @test_throws ErrorException set_integer(vref3)
        @test_throws ErrorException set_integer(y0)
        @test_throws ErrorException set_integer(yb)
    end
end

# Test relaxation 
@testset "JuMP.relax_integrality" begin 
    # setup the model 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t), Int)
    @variable(m, 0.25 <= x <= 1.25, Bin)
    @variable(m, z, Bin)
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
