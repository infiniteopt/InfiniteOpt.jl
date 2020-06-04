# Test constraint deletion
@testset "JuMP.delete (Constraints)" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x)
    var = build_variable(error, inf, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var, define_name = false)
    data = TestData(par, 0, 1)
    meas = measure(inf + par - x, data)
    @constraint(m, cref, par - inf + pt + 2x - rv + meas <= 1)
    # test normal deletion
    @test isa(delete(m, cref), Nothing)
    @test !is_valid(m, cref)
    @test !used_by_constraint(par)
    @test !used_by_constraint(inf)
    @test !used_by_constraint(pt)
    @test !used_by_constraint(x)
    @test !used_by_constraint(rv)
    @test !used_by_constraint(meas)
    @test isempty(InfiniteOpt._data_dictionary(cref))
    # test assertion error
    @test_throws AssertionError delete(m, cref)
end

# Test deleting variable information
@testset "Variable Information" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, 0 <= inf1(par) <= 1, Bin)
    @infinite_variable(m, inf2(par) == 1, Int)
    @point_variable(m, inf1(0.5), pt1)
    @point_variable(m, inf2(0.5), pt2)
    @hold_variable(m, 0 <= gb1 <= 1, Bin)
    @hold_variable(m, gb2 == 1, Int)
    # test delete_lower_bound
    @testset "JuMP.delete_lower_bound" begin
        # test with infinite variable
        @test has_lower_bound(inf1)
        @test isa(delete_lower_bound(inf1), Nothing)
        @test !has_lower_bound(inf1)
        @test_throws ErrorException delete_lower_bound(inf2)
        # test with point variable
        @test has_lower_bound(pt1)
        @test isa(delete_lower_bound(pt1), Nothing)
        @test !has_lower_bound(pt1)
        @test_throws ErrorException delete_lower_bound(pt2)
        # test with hold variable
        @test has_lower_bound(gb1)
        @test isa(delete_lower_bound(gb1), Nothing)
        @test !has_lower_bound(gb1)
        @test_throws ErrorException delete_lower_bound(gb2)
    end
    # test delete_upper_bound
    @testset "JuMP.delete_upper_bound" begin
        # test with infinite variable
        @test has_upper_bound(inf1)
        @test isa(delete_upper_bound(inf1), Nothing)
        @test !has_upper_bound(inf1)
        @test_throws ErrorException delete_upper_bound(inf2)
        # test with point variable
        @test has_upper_bound(pt1)
        @test isa(delete_upper_bound(pt1), Nothing)
        @test !has_upper_bound(pt1)
        @test_throws ErrorException delete_upper_bound(pt2)
        # test with hold variable
        @test has_upper_bound(gb1)
        @test isa(delete_upper_bound(gb1), Nothing)
        @test !has_upper_bound(gb1)
        @test_throws ErrorException delete_upper_bound(gb2)
    end
    # test unfix
    @testset "JuMP.unfix" begin
        # test with infinite variable
        @test is_fixed(inf2)
        @test isa(unfix(inf2), Nothing)
        @test !is_fixed(inf2)
        @test_throws ErrorException unfix(inf1)
        # test with point variable
        @test is_fixed(pt2)
        @test isa(unfix(pt2), Nothing)
        @test !is_fixed(pt2)
        @test_throws ErrorException unfix(pt1)
        # test with hold variable
        @test is_fixed(gb2)
        @test isa(unfix(gb2), Nothing)
        @test !is_fixed(gb2)
        @test_throws ErrorException unfix(gb1)
    end
    # test unset_binary
    @testset "JuMP.unset_binary" begin
        # test with infinite variable
        @test is_binary(inf1)
        @test isa(unset_binary(inf1), Nothing)
        @test !is_binary(inf1)
        @test_throws ErrorException unset_binary(inf2)
        # test with point variable
        @test is_binary(pt1)
        @test isa(unset_binary(pt1), Nothing)
        @test !is_binary(pt1)
        @test_throws ErrorException unset_binary(pt2)
        # test with hold variable
        @test is_binary(gb1)
        @test isa(unset_binary(gb1), Nothing)
        @test !is_binary(gb1)
        @test_throws ErrorException unset_binary(gb2)
    end
    # test unset_integer
    @testset "JuMP.unset_integer" begin
        # test with infinite variable
        @test is_integer(inf2)
        @test isa(unset_integer(inf2), Nothing)
        @test !is_integer(inf2)
        @test_throws ErrorException unset_integer(inf1)
        # test with point variable
        @test is_integer(pt2)
        @test isa(unset_integer(pt2), Nothing)
        @test !is_integer(pt2)
        @test_throws ErrorException unset_integer(pt1)
        # test with hold variable
        @test is_integer(gb2)
        @test isa(unset_integer(gb2), Nothing)
        @test !is_integer(gb2)
        @test_throws ErrorException unset_integer(gb1)
    end
end

# Test scalar parameter deletion
@testset "Scalar Parameters" begin
    # initialize model
    m = InfiniteModel()
    # setup the parameters
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_parameter(m, 0 <= par3 <= 1)
    @finite_parameter(m, fin, 42)
    @finite_parameter(m, fin2, 42)
    dpar = dispatch_variable_ref(par)
    dpars = dispatch_variable_ref.(pars)
    dpar2 = dispatch_variable_ref(par2)
    dpar3 = dispatch_variable_ref(par3)
    dfin = dispatch_variable_ref(fin)
    dfin2 = dispatch_variable_ref(fin2)
    # setup the variables
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par, pars))
    @infinite_variable(m, inf4(par, par2, pars))
    @point_variable(m, inf(0.5), pt)
    pt2 = @point_variable(m, inf2(0.5, 0.5))
    @point_variable(m, inf3(0, [0, 0]), pt3)
    @hold_variable(m, x)
    var = build_variable(error, inf4, Dict{Int, Float64}(2 => 0.5), check = false)
    rv = add_variable(m, var, define_name = false)
    dinf = dispatch_variable_ref(inf)
    dinf2 = dispatch_variable_ref(inf2)
    dinf3 = dispatch_variable_ref(inf3)
    dinf4 = dispatch_variable_ref(inf4)
    dpt = dispatch_variable_ref(pt)
    dpt2 = dispatch_variable_ref(pt2)
    dpt3 = dispatch_variable_ref(pt3)
    dx = dispatch_variable_ref(x)
    drv = dispatch_variable_ref(rv)
    # setup the measure
    data = TestData(par, 0, 1)
    mref = measure(inf + par - x + rv + par2 + par3 + fin, data)
    mref2 = measure(par2, data)
    dmref = dispatch_variable_ref(mref)
    dmref2 = dispatch_variable_ref(mref2)
    # setup the constraints
    @constraint(m, con, inf2 + inf4 - par2 + par3 + fin <= 0)
    constr = ScalarConstraint(par2, MOI.GreaterThan(0.))
    con2 = add_constraint(m, constr)
    # setup objective
    set_objective(m, MOI.MIN_SENSE, fin)
    # test _check_param_in_data
    @testset "_check_param_in_data" begin
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, data)
        @test isa(InfiniteOpt._check_param_in_data(par2, data), Nothing)
        @test isa(InfiniteOpt._check_param_in_data(pars[1], data), Nothing)
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, BadData())
    end
    # test _update_reduced_variable
    @testset "_update_reduced_variable" begin
        # test removing single parameter that is not reduced
        _update_variable_param_refs(dinf4, IC.VectorTuple(par2, pars))
        @test isa(InfiniteOpt._update_reduced_variable(drv, 1:1), Nothing)
        @test infinite_variable_ref(drv) == inf4
        @test eval_supports(drv) == Dict(1 => 0.5)
        @test set_name(drv, "") isa Nothing
        @test name(drv) == "inf4(0.5, [pars[1], pars[2]])"
        # Undo changes
        InfiniteOpt._set_core_variable_object(drv, var)
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, par2, pars))
        # test removing single parameter that is reduced
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, pars))
        @test eval_supports(drv) == Dict(2 => 0.5)
        @test isa(InfiniteOpt._update_reduced_variable(drv, 2:2), Nothing)
        @test infinite_variable_ref(drv) == inf4
        @test eval_supports(drv) == Dict{Int, Float64}()
        @test set_name(drv, "") isa Nothing
        @test name(drv) == "inf4(par, [pars[1], pars[2]])"
        # Undo changes
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, par2, pars))
        # test removing a different single parameter
        eval_supports(drv)[1] = 0.5
        eval_supports(drv)[2] = 0.5
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, pars))
        @test isa(InfiniteOpt._update_reduced_variable(drv, 2:2), Nothing)
        @test infinite_variable_ref(drv) == inf4
        @test eval_supports(drv) == Dict(1 => 0.5)
        @test set_name(drv, "") isa Nothing
        @test name(drv) == "inf4(0.5, [pars[1], pars[2]])"
        # Undo changes
        InfiniteOpt._set_core_variable_object(drv, var)
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, par2, pars))
        # test removing array element
        eval_supports(drv)[3] = 0.2
        eval_supports(drv)[4] = 0.1
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, par2))
        @test isa(InfiniteOpt._update_reduced_variable(drv, 3:4), Nothing)
        @test infinite_variable_ref(drv) == inf4
        @test eval_supports(drv) == Dict(2 => 0.5)
        @test set_name(drv, "") isa Nothing
        @test name(drv) == "inf4(par, 0.5)"
        # Undo changes
        var = build_variable(error, inf4, Dict{Int, Float64}(2 => 0.5), check = false)
        InfiniteOpt._set_core_variable_object(drv, var)
        _update_variable_param_refs(dinf4, IC.VectorTuple(par, par2, pars))
    end
    # test _update_measures
    @testset "_update_measures" begin
        @test InfiniteOpt._update_measures(m, par2) isa Nothing
        @test measure_function(dmref) == inf + par - x + rv + par3 + fin
        @test measure_function(dmref2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        @test InfiniteOpt._object_numbers(dmref) == [2, 3, 4]
        @test InfiniteOpt._object_numbers(dmref2) == []
        @test InfiniteOpt._parameter_numbers(dmref) == [2, 3, 4, 5]
        @test InfiniteOpt._parameter_numbers(dmref2) == []
        # undo changes
        meas = Measure(inf + par - x + rv + par2 + par3 + fin, data, [2, 4], [2, 5], false)
        InfiniteOpt._set_core_variable_object(dmref, meas)
        meas = Measure(par2, data, [2], [2], false)
        InfiniteOpt._set_core_variable_object(dmref2, meas)
    end
    # test _update_constraints
    @testset "_update_constraints" begin
        @test InfiniteOpt._update_constraints(m, par2) isa Nothing
        @test jump_function(constraint_object(con)) == inf2 + inf4 + par3 + fin
        @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2, 3, 4]))
        @test InfiniteOpt._object_numbers(con2) == []
        # undo changes
        constr = ScalarConstraint(inf2 + inf4 - par2 + par3 + fin, MOI.LessThan(0.))
        InfiniteOpt._set_core_constraint_object(con, constr)
        constr = ScalarConstraint(par2, MOI.GreaterThan(0.))
        InfiniteOpt._set_core_constraint_object(con2, constr)
    end
    # test _update_number_list
    @testset "_update_number_list" begin
        list = [1, 2, 3, 4, 5]
        @test InfiniteOpt._update_number_list([2, 3], list) isa Nothing
        @test list == [1, 2, 3]
        list = [1, 3, 4, 5]
        @test InfiniteOpt._update_number_list([2], list) isa Nothing
        @test list == [1, 2, 3, 4]
    end
    # test JuMP.delete for IndependentParameters
    @testset "JuMP.delete (IndependentParameterRef)" begin
        # test normal usage
        @test isa(delete(m, par2), Nothing)
        @test !is_valid(m, par2)
        @test measure_function(dmref) == inf + par - x + rv + par3 + fin
        @test measure_function(dmref2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        @test parameter_refs(inf2) == (par,)
        @test parameter_refs(inf4) == (par, pars)
        @test parameter_values(pt2) == (0.5,)
        @test parameter_refs(rv) == (par, pars)
        @test name(inf2) == "inf2(par)"
        @test name(inf4) == "inf4(par, pars)"
        @test name(pt2) == "inf2(0.5)"
        @test set_name(rv, "") isa Nothing
        @test name(rv) == "inf4(par, [pars[1], pars[2]])"
        @test jump_function(constraint_object(con)) == inf2 + inf4 + par3 + fin
        @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        expected = [IndependentParameterIndex(1), DependentParametersIndex(1),
                    IndependentParameterIndex(3)]
        @test InfiniteOpt._param_object_indices(m) == expected
        @test InfiniteOpt._last_param_num(m) == 4
        @test InfiniteOpt._object_number(dpar) == 1
        @test InfiniteOpt._object_number(dpars[1]) == 2
        @test InfiniteOpt._object_number(dpar3) == 3
        @test InfiniteOpt._parameter_number(dpar3) == 4
        @test InfiniteOpt._parameter_number(dpars[2]) == 3
        @test InfiniteOpt._object_numbers(dinf2) == [1]
        @test InfiniteOpt._parameter_numbers(dinf4) == [1, 2, 3]
        @test isempty(setdiff(InfiniteOpt._object_numbers(inf3), [1, 2]))
        @test InfiniteOpt._parameter_numbers(dinf) == [1]
        @test isempty(setdiff(InfiniteOpt._object_numbers(drv), [1, 2]))
        @test InfiniteOpt._object_numbers(dmref) == [3]
        @test InfiniteOpt._parameter_numbers(dmref) == [4]
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2, 3]))
        @test isempty(InfiniteOpt._object_numbers(con2))
        # test special measure case with single parameter (possible through
        # multiple deletions) and with single parameter in constraint
        @test isa(delete(m, par3), Nothing)
        @test measure_function(dmref) == inf + par - x + rv + fin
        @test jump_function(constraint_object(con)) == inf2 + inf4 + fin
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2]))
        expected = [IndependentParameterIndex(1), DependentParametersIndex(1)]
        @test InfiniteOpt._param_object_indices(m) == expected
        @test InfiniteOpt._last_param_num(m) == 3
        @test InfiniteOpt._parameter_number(dpars[2]) == 3
        @test InfiniteOpt._object_numbers(dinf2) == [1]
        @test InfiniteOpt._parameter_numbers(dinf4) == [1, 2, 3]
        # test invalid parameter
        @test_throws AssertionError delete(m, par2)
        @test_throws AssertionError delete(m, par3)
        # test measure data check
        @test_throws ErrorException delete(m, par)
    end
    # test delete for finite parameters
    @testset "JuMP.delete (FiniteParameterRef)" begin
        # test first case
        @test delete(m, fin) isa Nothing
        @test measure_function(dmref) == inf + par - x + rv
        @test jump_function(constraint_object(con)) == inf2 + inf4
        @test objective_function(m) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        @test objective_sense(m) == MOI.FEASIBILITY_SENSE
        # test error
        @test_throws AssertionError delete(m, fin)
        # test different objective
        @objective(m, Min, fin2 + x + 1)
        @test delete(m, fin2) isa Nothing
        @test objective_function(m) == x + 1
        @test objective_sense(m) == MOI.MIN_SENSE
    end
 end

 # test deletion of dependent parameters
 @testset "Dependent Parameters" begin
     # initialize model
     m = InfiniteModel()
     # setup the parameters
     @infinite_parameter(m, 0 <= par <= 1)
     @infinite_parameter(m, 0 <= pars[1:2] <= 1)
     @infinite_parameter(m, 0 <= par2 <= 1)
     dpar = dispatch_variable_ref(par)
     dpars = dispatch_variable_ref.(pars)
     dpar2 = dispatch_variable_ref(par2)
     # setup the variables
     @infinite_variable(m, inf(par))
     @infinite_variable(m, inf2(par, par2))
     @infinite_variable(m, inf3(par, pars))
     @infinite_variable(m, inf4(par, par2, pars))
     @point_variable(m, inf(0.5), pt)
     pt2 = @point_variable(m, inf2(0.5, 0.5))
     @point_variable(m, inf3(0, [0, 0]), pt3)
     @hold_variable(m, x)
     var = build_variable(error, inf4, Dict{Int, Float64}(2 => 0.5, 3 => 0, 4 => 4),
                          check = false)
     rv = add_variable(m, var, define_name = false)
     dinf = dispatch_variable_ref(inf)
     dinf2 = dispatch_variable_ref(inf2)
     dinf3 = dispatch_variable_ref(inf3)
     dinf4 = dispatch_variable_ref(inf4)
     dpt = dispatch_variable_ref(pt)
     dpt2 = dispatch_variable_ref(pt2)
     dpt3 = dispatch_variable_ref(pt3)
     dx = dispatch_variable_ref(x)
     drv = dispatch_variable_ref(rv)
     # setup the measure
     data = TestData(par, 0, 1)
     mref = measure(inf + par - x + rv + pars[1], data)
     data = TestData(pars[2], 0, 1)
     mref2 = measure(inf4, data)
     dmref = dispatch_variable_ref(mref)
     dmref2 = dispatch_variable_ref(mref2)
     # setup the constraints
     @constraint(m, con, inf2 + inf4 - par2 + pars[1] + pars[2] <= 0)
     # test delete for dependent parameters
     @testset "JuMP.delete (DependentParameterRefs)" begin
         # test measure error
         @test_throws ErrorException delete(m, pars)
         mindex = index(dmref2)
         filter!(e -> e != mindex, InfiniteOpt._measure_dependencies(inf4))
         filter!(e -> e != mindex, InfiniteOpt._measure_dependencies(pars[2]))
         InfiniteOpt._delete_data_object(dmref2)
         # test regular
         @test delete(m, pars) isa Nothing
         @test parameter_refs(dinf) == (par,)
         @test parameter_refs(dinf3) == (par,)
         @test parameter_refs(dinf4) == (par, par2)
         @test name(dinf) == "inf(par)"
         @test name(dinf4) == "inf4(par, par2)"
         @test set_name(drv, "") isa Nothing
         @test name(drv) == "inf4(par, 0.5)"
         @test parameter_values(dpt3) == (0,)
         @test parameter_values(dpt2) == (0.5, 0.5)
         @test measure_function(dmref) == inf + par - x + rv
         @test jump_function(constraint_object(con)) == inf2 + inf4 - par2
         expected = [IndependentParameterIndex(1), IndependentParameterIndex(2)]
         @test InfiniteOpt._param_object_indices(m) == expected
         @test InfiniteOpt._last_param_num(m) == 2
         @test InfiniteOpt._object_number(dpar) == 1
         @test InfiniteOpt._object_number(dpar2) == 2
         @test InfiniteOpt._parameter_number(dpar2) == 2
         @test isempty(setdiff(InfiniteOpt._object_numbers(dinf2), [1, 2]))
         @test InfiniteOpt._parameter_numbers(dinf4) == [1, 2]
         @test InfiniteOpt._object_numbers(inf3) == [1]
         @test InfiniteOpt._parameter_numbers(dinf) == [1]
         @test InfiniteOpt._object_numbers(drv) == [1]
         @test InfiniteOpt._object_numbers(dmref) == []
         @test InfiniteOpt._parameter_numbers(dmref) == []
         @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2]))
         # test assertion error
         @test_throws AssertionError delete(m, pars)
     end
  end

 # Test reduced variable deletion
 @testset "JuMP.delete (Reduced Variables)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @infinite_parameter(m, 0 <= par2 <= 1)
     @infinite_variable(m, inf(par, par2))
     @point_variable(m, inf(0.5, 0.5), pt)
     @hold_variable(m, x)
     var = build_variable(error, inf, Dict{Int, Float64}(2 => 0.5), check = false)
     rv = add_variable(m, var, define_name = false)
     rv2 = add_variable(m, var, define_name = false)
     data = TestData(par, 0, 1)
     meas = measure(inf + par - x + rv, data)
     meas2 = measure(rv2, data)
     @constraint(m, con, x + rv <= 0)
     constr = ScalarConstraint(rv2, MOI.GreaterThan(0.))
     con2 = add_constraint(m, constr)
     # test normal deletion
     @test isa(delete(m, rv), Nothing)
     @test measure_function(meas) == inf + par - x
     @test InfiniteOpt._object_numbers(meas) == [2]
     @test jump_function(constraint_object(con)) == x + 0
     @test InfiniteOpt._object_numbers(con) == []
     @test InfiniteOpt._reduced_variable_dependencies(inf) == [JuMP.index(rv2)]
     @test !haskey(InfiniteOpt._data_dictionary(m, ReducedVariable), JuMP.index(rv))
     # test deletion of special cases
     @test isa(delete(m, rv2), Nothing)
     @test measure_function(meas2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(meas2) == []
     @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(con2) == []
     @test InfiniteOpt._reduced_variable_dependencies(inf) == []
     @test !haskey(InfiniteOpt._data_dictionary(m, ReducedVariable), JuMP.index(rv2))
     # test error
     @test_throws AssertionError delete(m, rv)
     @test_throws AssertionError delete(m, rv2)
 end

 # Test variable deletion
 @testset "JuMP.delete (Hold Variables)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @hold_variable(m, 0 <= x <= 1, Bin, parameter_bounds = (par == 0))
     @hold_variable(m, y == 1, Int)
     data = TestData(par, 0, 0)
     meas1 = measure(x + y + par, data)
     meas2 = measure(y, data)
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     @objective(m, Min, x + y)
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test num_constraints(m) == 4
     @test measure_function(meas1) == y + par
     @test jump_function(constraint_object(con1)) == y + par
     @test objective_function(m) == y + 0
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(x))
     # test deletion of y
     set_objective_function(m, y)
     @test isa(delete(m, y), Nothing)
     @test num_constraints(m) == 2
     @test measure_function(meas1) == par + 0
     @test measure_function(meas2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test jump_function(constraint_object(con1)) == par + 0
     @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test objective_function(m) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(y))
     # test errors
     @test_throws AssertionError delete(m, x)
     @test_throws AssertionError delete(m, y)
 end

 # Test variable deletion
 @testset "JuMP.delete (Point Variables)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @infinite_variable(m, inf(par))
     @point_variable(m, inf(0), 0 <= x <= 1, Bin)
     @point_variable(m, inf(1), y == 1, Int)
     data = TestData(par, 0, 1)
     meas1 = measure(x + y + par, data)
     meas2 = measure(y, data)
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     @objective(m, Min, x + y)
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test num_constraints(m) == 4
     @test measure_function(meas1) == y + par
     @test InfiniteOpt._object_numbers(meas1) == []
     @test jump_function(constraint_object(con1)) == y + par
     @test InfiniteOpt._object_numbers(con1) == [1]
     @test objective_function(m) == y + 0
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(x))
     # test deletion of y
     set_objective_function(m, y)
     @test isa(delete(m, y), Nothing)
     @test num_constraints(m) == 2
     @test measure_function(meas1) == par + 0
     @test measure_function(meas2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test jump_function(constraint_object(con1)) == par + 0
     @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test objective_function(m) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(y))
     # test errors
     @test_throws AssertionError delete(m, x)
     @test_throws AssertionError delete(m, y)
 end

 # Test variable deletion
 @testset "JuMP.delete (Infinite Variables)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @infinite_variable(m, 0 <= x(par) <= 1, Bin)
     @infinite_variable(m, y(par) == 1, Int)
     @point_variable(m, x(0), x0)
     var = build_variable(error, x, Dict{Int, Float64}(1 => 0.5), check = false)
     rv = add_variable(m, var, define_name = false)
     data = TestData(par, 0, 1)
     meas1 = measure(x + y + par, data)
     meas2 = measure(y, data)
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test num_constraints(m) == 4
     @test measure_function(meas1) == y + par
     @test InfiniteOpt._object_numbers(meas1) == []
     @test jump_function(constraint_object(con1)) == y + par
     @test InfiniteOpt._object_numbers(con1) == [1]
     @test InfiniteOpt._infinite_variable_dependencies(par) == [index(y)]
     @test !is_valid(m, rv)
     @test !is_valid(m, x0)
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(x))
     # test deletion of y
     @test isa(delete(m, y), Nothing)
     @test num_constraints(m) == 2
     @test measure_function(meas1) == par + 0
     @test measure_function(meas2) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(meas1) == []
     @test jump_function(constraint_object(con1)) == par + 0
     @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(con1) == [1]
     @test InfiniteOpt._object_numbers(con2) == []
     @test InfiniteOpt._infinite_variable_dependencies(par) == []
     @test !haskey(InfiniteOpt._data_dictionary(m, HoldVariable), JuMP.index(y))
     # test errors
     @test_throws AssertionError delete(m, x)
     @test_throws AssertionError delete(m, y)
 end

 # Test variable deletion
 @testset "JuMP.delete (Measures)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @infinite_parameter(m, 0 <= par2 <= 1)
     @infinite_parameter(m, 0 <= pars[1:2] <= 1)
     @infinite_variable(m, w(pars))
     @infinite_variable(m, x(par))
     @infinite_variable(m, y(par2))
     @point_variable(m, x(0), x0)
     var = build_variable(error, x, Dict{Int, Float64}(1 => 0.5), check = false)
     rv = add_variable(m, var, define_name = false)
     data = TestData(par, 0, 1)
     data2 = TestData(pars, [0, 0], [1, 1])
     meas = measure(x, data)
     meas1 = measure(x + x0 + rv + par + meas, data)
     meas2 = measure(x + x0 + rv + par, data)
     meas3 = measure(meas1 + x0, data)
     meas4 = measure(meas2, data)
     meas5 = measure(w, data2)
     @constraint(m, con1, x0 + meas1 <= 0)
     con2 = add_constraint(m, ScalarConstraint(meas2, MOI.LessThan(0.)))
     @objective(m, Min, meas1 + x0)
     # test deletion of meas1
     @test isa(delete(m, meas1), Nothing)
     @test measure_function(meas3) == x0 + 0
     @test InfiniteOpt._object_numbers(meas3) == []
     @test jump_function(constraint_object(con1)) == x0 + 0
     @test InfiniteOpt._object_numbers(con1) == []
     @test objective_function(m) == x0 + 0
     @test InfiniteOpt._measure_dependencies(x) == [JuMP.index(meas), JuMP.index(meas2)]
     @test InfiniteOpt._measure_dependencies(y) == []
     @test length(InfiniteOpt._measure_dependencies(par)) == 4
     @test InfiniteOpt._measure_dependencies(rv) == [JuMP.index(meas2)]
     @test !haskey(InfiniteOpt._data_dictionary(m, Measure), JuMP.index(meas1))
     # test deletion of meas2
     set_objective_function(m, meas2)
     @test isa(delete(m, meas2), Nothing)
     @test measure_function(meas4) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(meas4) == []
     @test jump_function(constraint_object(con2)) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._object_numbers(con2) == []
     @test objective_function(m) == zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
     @test InfiniteOpt._measure_dependencies(x) == [JuMP.index(meas)]
     @test InfiniteOpt._measure_dependencies(y) == []
     @test length(InfiniteOpt._measure_dependencies(par)) == 3
     @test InfiniteOpt._measure_dependencies(rv) == []
     @test !haskey(InfiniteOpt._data_dictionary(m, Measure), JuMP.index(meas2))
     # test deletion of meas5
     @test isa(delete(m, meas5), Nothing)
     @test InfiniteOpt._measure_dependencies(pars[1]) == []
     @test InfiniteOpt._measure_dependencies(pars[2]) == []
     # test errors
     @test_throws AssertionError delete(m, meas1)
     @test_throws AssertionError delete(m, meas2)
 end
