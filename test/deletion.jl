# Test constraint deletion
@testset "JuMP.delete (Constraints)" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, rv, SemiInfinite(inf, 0.5))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, x)
    data = TestData(par, 0, 1)
    meas = measure(inf + par - x, data)
    d1 = @deriv(inf, par)
    @constraint(m, c1, par - inf + pt + 2x - rv + meas - d1 <= 1)
    @constraint(m, c2, inf <= 9, DomainRestrictions(par => 0))
    @constraint(m, c3, [inf, x] in MOI.Zeros(2))
    # test normal deletion
    @test isa(delete(m, c1), Nothing)
    @test !is_valid(m, c1)
    @test !used_by_constraint(par)
    @test !(index(c1) in InfiniteOpt._constraint_dependencies(inf))
    @test !used_by_constraint(pt)
    @test !(index(c1) in InfiniteOpt._constraint_dependencies(x))
    @test !used_by_constraint(rv)
    @test !used_by_constraint(meas)
    @test !used_by_constraint(d1)
    @test length(InfiniteOpt._data_dictionary(c1)) == 2
    # test assertion error
    @test_throws AssertionError delete(m, c1)
    # test restricted deletion
    @test isa(delete(m, c2), Nothing)
    @test !is_valid(m, c2)
    @test !(index(c2) in InfiniteOpt._constraint_dependencies(inf))
    @test length(InfiniteOpt._data_dictionary(c2)) == 1
    @test isempty(m.constraint_restrictions)
    # test vector constraint deletion 
    @test isa(delete(m, c3), Nothing)
    @test !is_valid(m, c3)
    @test !used_by_constraint(inf)
    @test !used_by_constraint(x)
    @test isempty(InfiniteOpt._data_dictionary(c3))
end

# Test deleting variable information
@testset "Variable Information" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, 0 <= inf1 <= 1, Infinite(par), Bin)
    @variable(m, inf2 == 1, Infinite(par), Int)
    @variable(m, pt1, Point(inf1, 0.5))
    @variable(m, pt2, Point(inf2, 0.5))
    @variable(m, 0 <= gb1 <= 1, Bin)
    @variable(m, gb2 == 1, Int)
    @variable(m, 0 <= d1 <= 1, Deriv(inf1, par))
    @variable(m, d2 == 1, Deriv(inf2, par))
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
        # test with finite variable
        @test has_lower_bound(gb1)
        @test isa(delete_lower_bound(gb1), Nothing)
        @test !has_lower_bound(gb1)
        @test_throws ErrorException delete_lower_bound(gb2)
        # test with derivative
        @test has_lower_bound(d1)
        @test isa(delete_lower_bound(d1), Nothing)
        @test !has_lower_bound(d1)
        @test_throws ErrorException delete_lower_bound(d2)
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
        # test with finite variable
        @test has_upper_bound(gb1)
        @test isa(delete_upper_bound(gb1), Nothing)
        @test !has_upper_bound(gb1)
        @test_throws ErrorException delete_upper_bound(gb2)
        # test with derivative
        @test has_upper_bound(d1)
        @test isa(delete_upper_bound(d1), Nothing)
        @test !has_upper_bound(d1)
        @test_throws ErrorException delete_upper_bound(d2)
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
        # test with finite variable
        @test is_fixed(gb2)
        @test isa(unfix(gb2), Nothing)
        @test !is_fixed(gb2)
        @test_throws ErrorException unfix(gb1)
        # test with derivative
        @test is_fixed(d2)
        @test isa(unfix(d2), Nothing)
        @test !is_fixed(d2)
        @test_throws ErrorException unfix(d1)
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
        # test with finite variable
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
        # test with finite variable
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
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, par3 in [0, 1])
    @finite_parameter(m, fin == 42)
    @finite_parameter(m, fin2 == 42)
    dpar = dispatch_variable_ref(par)
    dpars = dispatch_variable_ref.(pars)
    dpar2 = dispatch_variable_ref(par2)
    dpar3 = dispatch_variable_ref(par3)
    dfin = dispatch_variable_ref(fin)
    dfin2 = dispatch_variable_ref(fin2)
    # setup the variables
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, par2))
    @variable(m, inf3, Infinite(par, pars))
    @variable(m, inf4, Infinite(par, pars))
    @variable(m, pt, Point(inf, 0.5))
    pt2 = @variable(m, variable_type = Point(inf2, 0.5, 0.5))
    @variable(m, pt3, Point(inf3, 0, [0, 0]))
    @variable(m, x)
    var = build_variable(error, inf4, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var)
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
    # setup a derivative 
    d1 = @deriv(mref, par2)
    # setup the constraints
    @constraint(m, con, inf2 + inf4 - par2 + par3 + fin <= 0)
    constr = ScalarConstraint(par2, MOI.GreaterThan(0.))
    con2 = add_constraint(m, constr)
    @constraint(m, con3, [par2, par2] in MOI.Zeros(2))
    # setup objective
    set_objective(m, MOI.MIN_SENSE, fin)
    # test _check_param_in_data
    @testset "_check_param_in_data" begin
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, data)
        @test isa(InfiniteOpt._check_param_in_data(par2, data), Nothing)
        @test isa(InfiniteOpt._check_param_in_data(pars[1], data), Nothing)
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, BadData())
    end
    # test _update_measures
    @testset "_update_measures" begin
        @test InfiniteOpt._update_measures(m, par2) isa Nothing
        @test isequal_canonical(measure_function(dmref), inf + par - x + rv + par3 + fin)
        @test isequal(measure_function(dmref2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
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
        @test isequal_canonical(jump_function(constraint_object(con)), inf2 + inf4 + par3 + fin)
        @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2, 3, 4]))
        @test InfiniteOpt._object_numbers(con2) == []
        @test !is_valid(m, con3)
        # undo changes
        constr = ScalarConstraint(inf2 + inf4 - par2 + par3 + fin, MOI.LessThan(0.))
        InfiniteOpt._set_core_constraint_object(con, constr)
        constr = ScalarConstraint(par2, MOI.GreaterThan(0.))
        InfiniteOpt._set_core_constraint_object(con2, constr)
        unregister(m, :con3)
        @constraint(m, con3, [par2, par2] in MOI.Zeros(2))
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
        # add parameter function for object number updating 
        pfunc = parameter_function(sin, par)
        # test used by infinite variable
        @test_throws ErrorException delete(m, par2)
        @test delete(m, inf2) isa Nothing
        @test !is_valid(m, inf2)
        @test !is_valid(m, pt2)
        # test parameter function dependency 
        push!(InfiniteOpt._parameter_function_dependencies(par2), ParameterFunctionIndex(1))
        @test_throws ErrorException delete(m, par2)
        popfirst!(InfiniteOpt._parameter_function_dependencies(par2)) 
        # test normal usage
        @test isa(delete(m, par2), Nothing)
        @test !is_valid(m, par2)
        @test !is_valid(m, d1)
        @test isequal_canonical(measure_function(dmref), inf + par - x + rv + par3 + fin)
        @test isequal_canonical(measure_function(dmref2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
        @test isequal(parameter_refs(inf4), (par, pars))
        @test isequal(parameter_refs(rv), (pars,))
        @test isequal_canonical(jump_function(constraint_object(con)), inf4 + par3 + fin)
        @test isequal(jump_function(constraint_object(con2)),zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
        expected = [IndependentParameterIndex(1), DependentParametersIndex(1),
                    IndependentParameterIndex(3)]
        @test InfiniteOpt._param_object_indices(m) == expected
        @test InfiniteOpt._last_param_num(m) == 4
        @test InfiniteOpt._object_number(dpar) == 1
        @test InfiniteOpt._object_number(dpars[1]) == 2
        @test InfiniteOpt._object_number(dpar3) == 3
        @test InfiniteOpt._parameter_number(dpar3) == 4
        @test InfiniteOpt._parameter_number(dpars[2]) == 3
        @test InfiniteOpt._parameter_numbers(dinf4) == [1, 2, 3]
        @test isempty(setdiff(InfiniteOpt._object_numbers(inf3), [1, 2]))
        @test InfiniteOpt._parameter_numbers(dinf) == [1]
        @test isempty(setdiff(InfiniteOpt._object_numbers(drv), [1, 2]))
        @test InfiniteOpt._object_numbers(dmref) == [3]
        @test InfiniteOpt._parameter_numbers(dmref) == [4]
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2, 3]))
        @test isempty(InfiniteOpt._object_numbers(con2))
        @test !is_valid(m, con3)
        # test special measure case with single parameter (possible through
        # multiple deletions) and with single parameter in constraint
        @test isa(delete(m, par3), Nothing)
        @test isequal_canonical(measure_function(dmref), inf + par - x + rv + fin)
        @test isequal_canonical(jump_function(constraint_object(con)), inf4 + fin)
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2]))
        expected = [IndependentParameterIndex(1), DependentParametersIndex(1)]
        @test InfiniteOpt._param_object_indices(m) == expected
        @test InfiniteOpt._last_param_num(m) == 3
        @test InfiniteOpt._parameter_number(dpars[2]) == 3
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
        @test isequal_canonical(measure_function(dmref), inf + par - x + rv)
        @test isequal_canonical(jump_function(constraint_object(con)), inf4 + 0.)
        @test isequal(objective_function(m), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
        @test objective_sense(m) == MOI.FEASIBILITY_SENSE
        # test error
        @test_throws AssertionError delete(m, fin)
        # test different objective
        @objective(m, Min, fin2 + x + 1)
        @test delete(m, fin2) isa Nothing
        @test isequal_canonical(objective_function(m), x + 1)
        @test objective_sense(m) == MOI.MIN_SENSE
    end
end

# test deletion of dependent parameters
@testset "Dependent Parameters" begin
    # initialize model
    m = InfiniteModel()
    # setup the parameters
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    dpar = dispatch_variable_ref(par)
    dpars = dispatch_variable_ref.(pars)
    dpar2 = dispatch_variable_ref(par2)
    # setup the variables
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, par2))
    @variable(m, inf3, Infinite(par, pars))
    @variable(m, pt, Point(inf, 0.5))
    pt2 = @variable(m, variable_type = Point(inf2, 0.5, 0.5))
    @variable(m, pt3, Point(inf3, 0, [0, 0]))
    @variable(m, x)
    dinf = dispatch_variable_ref(inf)
    dinf2 = dispatch_variable_ref(inf2)
    dinf3 = dispatch_variable_ref(inf3)
    dpt = dispatch_variable_ref(pt)
    dpt2 = dispatch_variable_ref(pt2)
    dpt3 = dispatch_variable_ref(pt3)
    dx = dispatch_variable_ref(x)
    # setup the measure
    data = TestData(par, 0, 1)
    mref = measure(inf + par - x + pars[1], data)
    data = TestData(pars[2], 0, 1)
    mref2 = measure(inf2, data)
    dmref = dispatch_variable_ref(mref)
    dmref2 = dispatch_variable_ref(mref2)
    # setup the derivatives 
    d1 = @deriv(mref, pars[1])
    # setup the constraints
    @constraint(m, con, inf2 - par2 + pars[1] + pars[2] <= 0)
    @constraint(m, con2, pars in MOI.Zeros(2))
    # test delete for dependent parameters
    @testset "JuMP.delete (DependentParameterRefs)" begin
        # test measure error
        @test_throws ErrorException delete(m, pars)
        mindex = index(dmref2)
        filter!(e -> e != mindex, InfiniteOpt._measure_dependencies(inf2))
        filter!(e -> e != mindex, InfiniteOpt._measure_dependencies(pars[2]))
        InfiniteOpt._delete_data_object(dmref2)
        # test infinite variable error
        @test_throws ErrorException delete(m, pars)
        @test delete(m, inf3) isa Nothing 
        @test !is_valid(m, inf3)
        @test !is_valid(m, pt3)
        # test parameter function dependency 
        data = InfiniteOpt._data_object(first(pars))
        push!(data.parameter_func_indices, ParameterFunctionIndex(1))
        @test_throws ErrorException delete(m, pars)
        empty!(data.parameter_func_indices)
        # test regular
        @test delete(m, pars) isa Nothing
        @test !is_valid(m, d1)
        @test isequal(parameter_refs(dinf), (par,))
        @test string(dinf) == "inf(par)"
        @test parameter_values(dpt2) == (0.5, 0.5)
        @test isequal_canonical(measure_function(dmref), inf + par - x)
        @test isequal_canonical(jump_function(constraint_object(con)), inf2 - par2)
        expected = [IndependentParameterIndex(1), IndependentParameterIndex(2)]
        @test InfiniteOpt._param_object_indices(m) == expected
        @test InfiniteOpt._last_param_num(m) == 2
        @test InfiniteOpt._object_number(dpar) == 1
        @test InfiniteOpt._object_number(dpar2) == 2
        @test InfiniteOpt._parameter_number(dpar2) == 2
        @test isempty(setdiff(InfiniteOpt._object_numbers(dinf2), [1, 2]))
        @test InfiniteOpt._parameter_numbers(dinf) == [1]
        @test InfiniteOpt._object_numbers(dmref) == []
        @test InfiniteOpt._parameter_numbers(dmref) == []
        @test isempty(setdiff(InfiniteOpt._object_numbers(con), [1, 2]))
        @test !is_valid(m, con2)
        # test assertion error
        @test_throws AssertionError delete(m, pars)
    end
end

# Test semi_infinite variable deletion
@testset "JuMP.delete (SemiInfinite Variables)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @variable(m, inf, Infinite(par, par2))
    @variable(m, pt, Point(inf, 0.5, 0.5))
    @variable(m, x)
    var = build_variable(error, inf, Dict{Int, Float64}(2 => 0.5), check = false)
    rv = add_variable(m, var)
    var = build_variable(error, inf, Dict{Int, Float64}(2 => 0), check = false)
    rv2 = add_variable(m, var)
    data = TestData(par, 0, 1)
    meas = measure(inf + par - x + rv, data)
    meas2 = measure(rv2, data)
    d1 = @deriv(rv, par)
    @constraint(m, con, x + rv <= 0)
    constr = ScalarConstraint(rv2, MOI.GreaterThan(0.))
    con2 = add_constraint(m, constr)
    @constraint(m, con3, [rv, rv] in MOI.Zeros(2))
    # test normal deletion
    @test isa(delete(m, rv), Nothing)
    @test !is_valid(m, d1)
    @test isequal_canonical(measure_function(meas), inf + par - x)
    @test InfiniteOpt._object_numbers(meas) == [2]
    @test isequal_canonical(jump_function(constraint_object(con)), x + 0)
    @test InfiniteOpt._object_numbers(con) == []
    @test InfiniteOpt._semi_infinite_variable_dependencies(inf) == [JuMP.index(rv2)]
    @test !haskey(InfiniteOpt._data_dictionary(m, SemiInfiniteVariable), JuMP.index(rv))
    @test !is_valid(m, rv)
    # test deletion of special cases
    @test isa(delete(m, rv2), Nothing)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(meas2) == []
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(con2) == []
    @test InfiniteOpt._semi_infinite_variable_dependencies(inf) == []
    @test !haskey(InfiniteOpt._data_dictionary(m, SemiInfiniteVariable), JuMP.index(rv2))
    @test isempty(m.semi_lookup)
    # test error
    @test_throws AssertionError delete(m, rv)
    @test_throws AssertionError delete(m, rv2)
end

# Test variable deletion
@testset "JuMP.delete (Finite Variables)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, 0 <= x <= 1, Bin)
    @variable(m, y == 1, Int)
    data = TestData(par, 0, 0)
    meas1 = measure(x + y + par, data)
    meas2 = measure(y, data)
    @constraint(m, con1, x + y + par <= 0)
    con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
    @constraint(m, con3, [x, y] in MOI.Zeros(2))
    @constraint(m, con4, sin(x) <= 0)
    @objective(m, Min, x + y)
    # test deletion of x
    @test isa(delete(m, x), Nothing)
    @test num_constraints(m) == 5
    @test isequal_canonical(measure_function(meas1), y + par)
    @test isequal_canonical(jump_function(constraint_object(con1)), y + par)
    @test isequal_canonical(objective_function(m), y + 0)
    @test isequal_canonical(jump_function(constraint_object(con4)), GenericNonlinearExpr{GeneralVariableRef}(:-, Any[GenericNonlinearExpr{GeneralVariableRef}(:sin, Any[0.0]), 0.0]))
    @test !is_valid(m, con3)
    @test !haskey(InfiniteOpt._data_dictionary(m, FiniteVariable), JuMP.index(x))
    # test deletion of y
    set_objective_function(m, y)
    @test isa(delete(m, y), Nothing)
    @test num_constraints(m) == 3
    @test isequal_canonical(measure_function(meas1), par + 0)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test isequal_canonical(jump_function(constraint_object(con1)), par + 0)
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test isequal(objective_function(m), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test !haskey(InfiniteOpt._data_dictionary(m, FiniteVariable), JuMP.index(y))
    # test errors
    @test_throws AssertionError delete(m, x)
    @test_throws AssertionError delete(m, y)
end

# Test variable deletion
@testset "JuMP.delete (Point Variables)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, 0 <= x <= 1, Point(inf, 0), Bin)
    @variable(m, y == 1, Point(inf, 1), Int)
    data = TestData(par, 0, 1)
    meas1 = measure(x + y + par, data)
    meas2 = measure(y, data)
    @constraint(m, con1, x + y + par <= 0)
    con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
    @constraint(m, con3, [x, y] in MOI.Zeros(2))
    @objective(m, Min, x + y)
    # test deletion of x
    @test isa(delete(m, x), Nothing)
    @test num_constraints(m) == 4
    @test isequal_canonical(measure_function(meas1), y + par)
    @test InfiniteOpt._object_numbers(meas1) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), y + par)
    @test InfiniteOpt._object_numbers(con1) == [1]
    @test isequal_canonical(objective_function(m), y + 0)
    @test !haskey(InfiniteOpt._data_dictionary(m, PointVariable), JuMP.index(x))
    @test !is_valid(m, con3)
    # test deletion of y
    set_objective_function(m, y)
    @test isa(delete(m, y), Nothing)
    @test num_constraints(m) == 2
    @test isequal_canonical(measure_function(meas1), par + 0)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test isequal_canonical(jump_function(constraint_object(con1)), par + 0)
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test isequal(objective_function(m), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test !haskey(InfiniteOpt._data_dictionary(m, PointVariable), JuMP.index(y))
    @test isempty(m.point_lookup)
    # test errors
    @test_throws AssertionError delete(m, x)
    @test_throws AssertionError delete(m, y)
end

# Test variable deletion
@testset "JuMP.delete (Infinite Variables)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, 0 <= x <= 1, Infinite(par), Bin)
    @variable(m, y == 1, Infinite(par), Int)
    @variable(m, x0, Point(x, 0))
    var = build_variable(error, x, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var)
    data = TestData(par, 0, 1)
    meas1 = measure(x + y + par, data)
    meas2 = measure(y, data)
    d1 = @deriv(x, par)
    @constraint(m, con1, x + y + par <= 0)
    con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
    @constraint(m, con3, [x, y] in MOI.Zeros(2))
    constant_over_collocation(x, par)
    constant_over_collocation(y, par)
    # test deletion of x
    @test isa(delete(m, x), Nothing)
    @test num_constraints(m) == 4
    @test isequal_canonical(measure_function(meas1), y + par)
    @test InfiniteOpt._object_numbers(meas1) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), y + par)
    @test InfiniteOpt._object_numbers(con1) == [1]
    @test InfiniteOpt._infinite_variable_dependencies(par) == [index(y)]
    @test !is_valid(m, rv)
    @test !is_valid(m, x0)
    @test !is_valid(m, d1)
    @test !haskey(InfiniteOpt._data_dictionary(m, InfiniteVariable), JuMP.index(x))
    @test !is_valid(m, con3)
    @test m.piecewise_vars[index(par)] == Set(index(y))
    # test deletion of y
    @test isa(delete(m, y), Nothing)
    @test num_constraints(m) == 2
    @test isequal_canonical(measure_function(meas1), par + 0)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(meas1) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), par + 0)
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(con1) == [1]
    @test InfiniteOpt._object_numbers(con2) == []
    @test InfiniteOpt._infinite_variable_dependencies(par) == []
    @test !haskey(InfiniteOpt._data_dictionary(m, InfiniteVariable), JuMP.index(y))
    @test isempty(m.piecewise_vars)
    # test errors
    @test_throws AssertionError delete(m, x)
    @test_throws AssertionError delete(m, y)
end

# Test derivative deletion 
@testset "JuMP.delete (Derivatives)" begin 
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, 0 <= x <= 1, Infinite(par), Bin)
    @variable(m, y == 1, Infinite(par), Int)
    d1 = @deriv(x, par)
    d2 = @deriv(y, par^2)
    d3 = @deriv(d1, par)
    @variable(m, dx0, Point(d1, 0))
    var = build_variable(error, d1, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var)
    data = TestData(par, 0, 1)
    meas1 = measure(d1 + y + par, data)
    meas2 = measure(d2, data)
    @constraint(m, con1, d1 + d2 + par <= 0)
    con2 = add_constraint(m, ScalarConstraint(d2, MOI.LessThan(0.)))
    cref = @constraint(m, d1 == 0)
    push!(InfiniteOpt._derivative_constraint_dependencies(d1), index(cref))
    InfiniteOpt._set_has_derivative_constraints(par, true)
    @constraint(m, con3, [d1, d2] in MOI.Zeros(2))
    # test deletion of d1
    @test isa(delete(m, d1), Nothing)
    @test num_constraints(m) == 7
    @test isequal_canonical(measure_function(meas1), y + par)
    @test InfiniteOpt._object_numbers(meas1) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), d2 + par)
    @test InfiniteOpt._object_numbers(con1) == [1]
    @test InfiniteOpt._derivative_dependencies(par) == [index(d2)]
    @test !is_valid(m, rv)
    @test !is_valid(m, dx0)
    @test !is_valid(m, d3)
    @test !is_valid(m, cref)
    @test !is_valid(m, con3)
    @test !haskey(InfiniteOpt._data_dictionary(m, Derivative), JuMP.index(d1))
    @test !has_derivative_constraints(par)
    @test !haskey(m.deriv_lookup, (x, par, 1))
    # test deletion of d2
    @test isa(delete(m, d2), Nothing)
    @test num_constraints(m) == 7
    @test isequal_canonical(measure_function(meas1), y + par)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(meas1) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), par + 0)
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(con1) == [1]
    @test InfiniteOpt._object_numbers(con2) == []
    @test InfiniteOpt._derivative_dependencies(par) == []
    @test !haskey(InfiniteOpt._data_dictionary(m, Derivative), JuMP.index(d2))
    @test !haskey(m.deriv_lookup, (y, par, 2))
    # test errors
    @test_throws AssertionError delete(m, d1)
    @test_throws AssertionError delete(m, d2)
end

# Test infinite parameter function deletion
@testset "JuMP.delete (Parameter Function)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    f = parameter_function(sin, par)
    d = deriv(f, par)
    d = deriv(f, par)
    rv = add_variable(m, build_variable(error, f, Dict(1 => 0.)))
    data = TestData(par, 0, 0)
    meas1 = measure(f + par, data)
    meas2 = measure(f, data)
    @constraint(m, con1, f + par <= 0)
    con2 = add_constraint(m, ScalarConstraint(f, MOI.LessThan(0.)))
    @constraint(m, con3, [f, f] in MOI.Zeros(2))
    # test deletion of x
    @test isa(delete(m, f), Nothing)
    @test num_constraints(m) == 2
    @test isequal_canonical(measure_function(meas1), par + 0)
    @test isequal(measure_function(meas2), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test isequal_canonical(jump_function(constraint_object(con1)), par + 0)
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test !haskey(InfiniteOpt._data_dictionary(f), JuMP.index(f))
    @test !is_valid(m, f)
    @test !is_valid(m, d)
    @test !is_valid(m, rv)
    @test !is_valid(m, con3)
    # test errors
    @test_throws AssertionError delete(m, f)
end

# Test variable deletion
@testset "JuMP.delete (Measures)" begin
    # intialize the model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, w, Infinite(pars))
    @variable(m, x, Infinite(par))
    @variable(m, y, Infinite(par2))
    @variable(m, x0, Point(x, 0))
    var = build_variable(error, x, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var)
    data = TestData(par, 0, 1)
    data2 = TestData(pars, [0, 0], [1, 1])
    meas = measure(x, data)
    meas1 = measure(x + x0 + rv + par + meas, data)
    meas2 = measure(x + x0 + rv + par, data)
    meas3 = measure(meas1 + x0, data)
    meas4 = measure(meas2, data)
    meas5 = measure(w, data2)
    meas6 = measure(x + y, data)
    d1 = @deriv(meas6, par2)
    @constraint(m, con1, x0 + meas1 <= 0)
    con2 = add_constraint(m, ScalarConstraint(meas2, MOI.LessThan(0.)))
    @constraint(m, con3, [meas1, meas1] in MOI.Zeros(2))
    @objective(m, Min, meas1 + x0)
    # test deletion of meas1
    @test isa(delete(m, meas1), Nothing)
    @test isequal_canonical(measure_function(meas3), x0 + 0)
    @test InfiniteOpt._object_numbers(meas3) == []
    @test isequal_canonical(jump_function(constraint_object(con1)), x0 + 0)
    @test InfiniteOpt._object_numbers(con1) == []
    @test isequal_canonical(objective_function(m), x0 + 0)
    @test InfiniteOpt._measure_dependencies(x) == [JuMP.index(meas), JuMP.index(meas2), JuMP.index(meas6)]
    @test InfiniteOpt._measure_dependencies(y) == [JuMP.index(meas6)]
    @test length(InfiniteOpt._measure_dependencies(par)) == 5
    @test InfiniteOpt._measure_dependencies(rv) == [JuMP.index(meas2)]
    @test !haskey(InfiniteOpt._data_dictionary(m, Measure), JuMP.index(meas1))
    @test !is_valid(m, con3)
    # test deletion of meas2
    set_objective_function(m, meas2)
    @test isa(delete(m, meas2), Nothing)
    @test isequal(measure_function(meas4), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(meas4) == []
    @test isequal(jump_function(constraint_object(con2)), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._object_numbers(con2) == []
    @test isequal(objective_function(m), zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @test InfiniteOpt._measure_dependencies(x) == [JuMP.index(meas), JuMP.index(meas6)]
    @test InfiniteOpt._measure_dependencies(y) == [JuMP.index(meas6)]
    @test length(InfiniteOpt._measure_dependencies(par)) == 4
    @test InfiniteOpt._measure_dependencies(rv) == []
    @test !haskey(InfiniteOpt._data_dictionary(m, Measure), JuMP.index(meas2))
    # test deletion of meas5
    @test isa(delete(m, meas5), Nothing)
    @test InfiniteOpt._measure_dependencies(pars[1]) == []
    @test InfiniteOpt._measure_dependencies(pars[2]) == []
    # test deletion of meas6
    @test isa(delete(m, meas6), Nothing)
    @test !is_valid(m, d1)
    # test errors
    @test_throws AssertionError delete(m, meas1)
    @test_throws AssertionError delete(m, meas2)
end
