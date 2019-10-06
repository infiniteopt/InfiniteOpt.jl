# Test constraint deletion
@testset "JuMP.delete (Constraints)" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
    @constraint(m, cref, par - inf + pt + 2x - rv + meas <= 1)
    # test normal deletion
    @test isa(delete(m, cref), Nothing)
    @test !is_valid(m, cref)
    @test !haskey(m.param_to_constrs, JuMP.index(par))
    @test !haskey(m.var_to_constrs, JuMP.index(inf))
    @test !haskey(m.var_to_constrs, JuMP.index(pt))
    @test !haskey(m.var_to_constrs, JuMP.index(x))
    @test !haskey(m.reduced_to_constrs, JuMP.index(rv))
    @test !haskey(m.meas_to_constrs, JuMP.index(meas))
    @test !haskey(m.constrs, JuMP.index(cref))
    @test !haskey(m.constr_to_name, JuMP.index(cref))
    @test !haskey(m.constr_in_var_info, JuMP.index(cref))
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
    @global_variable(m, 0 <= gb1 <= 1, Bin)
    @global_variable(m, gb2 == 1, Int)
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
        # test with global variable
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
        # test with global variable
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
        # test with global variable
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
        # test with global variable
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
        # test with global variable
        @test is_integer(gb2)
        @test isa(unset_integer(gb2), Nothing)
        @test !is_integer(gb2)
        @test_throws ErrorException unset_integer(gb1)
    end
end

# Test parameter deletion
@testset "Parameters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= par3 <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par, pars))
    @infinite_variable(m, inf4(par, par2, pars))
    @point_variable(m, inf(0.5), pt)
    pt2 = @point_variable(m, inf2(0.5, 0.5))
    @point_variable(m, inf3(0, [0, 0]), pt3)
    @global_variable(m, x)
    m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
    m.infinite_to_reduced[JuMP.index(inf4)] = [-1]
    rv = ReducedInfiniteVariableRef(m, -1)
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    meas = measure(inf + par - x + rv, data)
    meas2 = measure(inf3 - x + pt2, data2)
    meas3 = measure(inf2 + par2, data)
    meas4 = add_measure(m, Measure(par3, data))
    meas5 = measure(pars[2] + par, data)
    @constraint(m, con, inf2 + inf4 - par2 <= 0)
    m.constrs[-1] = ScalarConstraint(par3, MOI.GreaterThan(0.))
    m.param_to_constrs[JuMP.index(par3)] = [-1]
    con2 = InfiniteConstraintRef(m, Int64(-1), JuMP.shape(m.constrs[-1]))
    set_name(con2, "")

    # test _check_param_in_data (scalar)
    @testset "_check_param_in_data (Scalar)" begin
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, data)
        @test isa(InfiniteOpt._check_param_in_data(par2, data), Nothing)
    end
    # test _check_param_in_data (array)
    @testset "_check_param_in_data (Array)" begin
        @test_throws ErrorException InfiniteOpt._check_param_in_data(pars[1],
                                                                     data2)
        @test isa(InfiniteOpt._check_param_in_data(par, data2), Nothing)
    end
    # test _check_param_in_data (fallback)
    @testset "_check_param_in_data (Fallback)" begin
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, Bad())
    end
    # test _contains_pref (scalar)
    @testset "_contains_pref (Scalar)" begin
        @test InfiniteOpt._contains_pref(par, par)
        @test !InfiniteOpt._contains_pref(par, par2)
    end
    # test _contains_pref (array)
    @testset "_contains_pref (Array)" begin
        @test InfiniteOpt._contains_pref(pars, pars[1])
        @test !InfiniteOpt._contains_pref(pars, par)
    end
    # test _remove_parameter
    @testset "_remove_parameter" begin
        # test removing single element
        @test InfiniteOpt._remove_parameter((par,), par) == ((), (1,))
        @test InfiniteOpt._remove_parameter((par2, par), par) == ((par2,), (2,))
        # test removing element from array
        new_pars = JuMP.Containers.SparseAxisArray(Dict((2,) => pars[2]))
        @test InfiniteOpt._remove_parameter((pars,),
                                            pars[1]) == ((new_pars,), (1, (1,)))
        @test InfiniteOpt._remove_parameter((par, pars),
                                        pars[1]) == ((par, new_pars), (2, (1,)))
        # test removing from array that becomes empty
        @test InfiniteOpt._remove_parameter((new_pars,), pars[2]) == ((), (1,))
        @test InfiniteOpt._remove_parameter((par, new_pars),
                                            pars[2]) == ((par,), (2,))
    end
    # test _update_infinite_variable
    @testset "_update_infinite_variable" begin
        # test updating when used by measure
        @test isa(InfiniteOpt._update_infinite_variable(inf3, (pars,)), Nothing)
        @test parameter_refs(inf3) == (pars,)
        @test name(inf3) == "inf3(pars)"
        @test name(meas2) == "measure(inf3(pars) - x + inf2(0.5, 0.5))"
        # Undo the changes
        @test isa(InfiniteOpt._update_infinite_variable(inf3, (par, pars)),
                  Nothing)
        @test parameter_refs(inf3) == (par, pars)
        @test name(inf3) == "inf3(par, pars)"
        @test name(meas2) == "measure(inf3(par, pars) - x + inf2(0.5, 0.5))"
    end
    # test _remove_parameter_values
    @testset "_remove_parameter_values" begin
        # test remove single element
        @test InfiniteOpt._remove_parameter_values((0, 1), (1,)) == (1,)
        @test InfiniteOpt._remove_parameter_values((0, 1), (2,)) == (0,)
        @test InfiniteOpt._remove_parameter_values(([0, 0],), (1,)) == ()
        # test remove array element
        arr = convert(JuMP.Containers.SparseAxisArray, [0, 1, 0])
        arr2 = JuMP.Containers.SparseAxisArray(Dict((1,) => 0, (2,) => 1,
                                                    (3,) => 0))
        filter!(x -> x.first != (2,), arr2.data)
        @test InfiniteOpt._remove_parameter_values((arr, 3), (1, (2,))) == (arr2, 3)
        arr3 = JuMP.Containers.SparseAxisArray(filter(x -> x.first != (1,),
                                                      arr2.data))
        @test InfiniteOpt._remove_parameter_values((arr2,), (1, (1,))) == (arr3,)
    end
    # test _update_point_variable
    @testset "_update_point_variable" begin
        # test point variable with single parameter
        @test isa(InfiniteOpt._update_point_variable(pt, ()), Nothing)
        @test parameter_values(pt) == ()
        @test name(pt) == "pt"
        # Undo changes
        @test isa(InfiniteOpt._update_point_variable(pt, (0.5,)), Nothing)
        @test parameter_values(pt) == (0.5,)
        @test name(pt) == "pt"
        # test point variable with multiple parameters
        @test isa(InfiniteOpt._update_point_variable(pt2, (0.5,)), Nothing)
        @test parameter_values(pt2) == (0.5,)
        @test name(pt2) == "inf2(0.5)"
        @test name(meas2) == "measure(inf3(par, pars) - x + inf2(0.5))"
        # Undo changes
        @test isa(InfiniteOpt._update_point_variable(pt2, (0.5, 0.5)), Nothing)
        @test parameter_values(pt2) == (0.5, 0.5)
        @test name(pt2) == "inf2(0.5, 0.5)"
        @test name(meas2) == "measure(inf3(par, pars) - x + inf2(0.5, 0.5))"
    end
    # test _update_reduced_variable
    @testset "_update_reduced_variable" begin
        # test removing single parameter that is not reduced
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par2, pars)),
                                                        Nothing)
        @test isa(InfiniteOpt._update_reduced_variable(rv, (1, )), Nothing)
        @test m.reduced_info[JuMP.index(rv)].infinite_variable_ref == inf4
        @test m.reduced_info[JuMP.index(rv)].eval_supports == Dict(1 => 0.5)
        @test name(rv) == "inf4(0.5, pars)"
        @test name(meas) == "measure(inf(par) + par - x + inf4(0.5, pars))"
        # Undo changes
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, par2, pars)),
                                                        Nothing)
        # test removing single parameter that is reduced
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, pars)),
                                                        Nothing)
        @test isa(InfiniteOpt._update_reduced_variable(rv, (2, )), Nothing)
        @test m.reduced_info[JuMP.index(rv)].infinite_variable_ref == inf4
        @test m.reduced_info[JuMP.index(rv)].eval_supports == Dict{Int,
                               Union{Number, JuMP.Containers.SparseAxisArray}}()
        @test name(rv) == "inf4(par, pars)"
        @test name(meas) == "measure(inf(par) + par - x + inf4(par, pars))"
        # Undo changes
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, par2, pars)),
                                                        Nothing)
        # test removing a different single parameter
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(1 => 0.5, 2 => 0.5))
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, pars)),
                                                        Nothing)
        @test isa(InfiniteOpt._update_reduced_variable(rv, (2, )), Nothing)
        @test m.reduced_info[JuMP.index(rv)].infinite_variable_ref == inf4
        @test m.reduced_info[JuMP.index(rv)].eval_supports == Dict(1 => 0.5)
        @test name(rv) == "inf4(0.5, pars)"
        @test name(meas) == "measure(inf(par) + par - x + inf4(0.5, pars))"
        # Undo changes
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, par2, pars)),
                                                        Nothing)
        # prepare for removing array element
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0])
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5, 3 => supp))
        new_supp = JuMP.Containers.SparseAxisArray(Dict((1,) => 0, (2,) => 0))
        filter!(x -> x.first != (1,), new_supp.data)
        # test removing array element
        @test isa(InfiniteOpt._update_reduced_variable(rv, (3, (1,))), Nothing)
        @test m.reduced_info[JuMP.index(rv)].infinite_variable_ref == inf4
        @test m.reduced_info[JuMP.index(rv)].eval_supports == Dict(2 => 0.5, 3 => new_supp)
        @test name(rv) == "inf4(par, 0.5, " * string(new_supp) * ")"
        @test name(meas) == "measure(inf(par) + par - x + " * name(rv) * ")"
        # Undo changes
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, par2, pars)),
                                                        Nothing)
    end
    # test JuMP.delete for parameters
    @testset "JuMP.delete (Parameters)" begin
        # test normal usage
        idx = JuMP.index(par2)
        @test isa(delete(m, par2), Nothing)
        @test name(meas3) == "measure(inf2(par))"
        @test !haskey(m.param_to_meas, idx)
        @test parameter_refs(inf2) == (par,)
        @test parameter_refs(inf4) == (par, pars)
        @test parameter_values(pt2) == (0.5,)
        @test parameter_refs(rv) == (par, pars)
        @test name(inf2) == "inf2(par)"
        @test name(inf4) == "inf4(par, pars)"
        @test name(pt2) == "inf2(0.5)"
        @test name(rv) == "inf4(par, pars)"
        @test !haskey(m.param_to_vars, idx)
        @test string(m.constrs[JuMP.index(con)].func) == "inf2(par) + " *
                                                         "inf4(par, pars)"
        @test !haskey(m.param_to_constrs, idx)
        @test !haskey(m.params, idx)
        @test !haskey(m.param_to_name, idx)
        @test !haskey(m.param_to_group_id, idx)
        # test special measure case with single parameter (possible through
        # multiple deletions) and with single parameter in constraint
        idx = JuMP.index(par3)
        @test isa(delete(m, par3), Nothing)
        @test name(meas4) == "measure(0)"
        @test !haskey(m.param_to_meas, idx)
        @test string(m.constrs[JuMP.index(con2)].func) == "0"
        @test !haskey(m.param_to_constrs, idx)
        @test !haskey(m.params, idx)
        @test !haskey(m.param_to_name, idx)
        @test !haskey(m.param_to_group_id, idx)
        # test invalid parameter
        @test_throws AssertionError delete(m, par2)
        @test_throws AssertionError delete(m, par3)
        # test measure data check
        @test_throws ErrorException delete(m, par)
        @test_throws ErrorException delete(m, pars[1])
        # test array element deletion
        delete!(m.measures, JuMP.index(meas2))
        m.param_to_meas[JuMP.index(pars[2])] = [JuMP.index(meas5)]
        delete!(m.var_to_meas, JuMP.index(inf3))
        idx = JuMP.index(pars[2])
        @test isa(delete(m, pars[2]), Nothing)
        @test name(meas5) == "measure(par)"
        @test !haskey(m.param_to_meas, idx)
        new_pars = convert(JuMP.Containers.SparseAxisArray, [pars[1]])
        @test parameter_refs(inf3) == (par, new_pars)
        @test parameter_refs(inf4) == (par, new_pars)
        new_vals = convert(JuMP.Containers.SparseAxisArray, [0])
        @test parameter_values(pt3) == (0, new_vals)
        @test parameter_refs(rv) == (par, new_pars)
        @test name(inf3) == "inf3(par, pars)"
        @test name(inf4) == "inf4(par, pars)"
        @test name(pt3) == "pt3"
        @test name(rv) == "inf4(par, pars)"
        @test !haskey(m.param_to_vars, idx)
        @test string(m.constrs[JuMP.index(con)].func) == "inf2(par) + " *
                                                         "inf4(par, pars)"
        @test !haskey(m.param_to_constrs, idx)
        @test !haskey(m.params, idx)
        @test !haskey(m.param_to_name, idx)
        @test !haskey(m.param_to_group_id, idx)
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
     @global_variable(m, x)
     m.reduced_info[-1] = ReducedInfiniteInfo(inf, Dict(2 => 0.5))
     m.infinite_to_reduced[JuMP.index(inf)] = [-1]
     rv = ReducedInfiniteVariableRef(m, -1)
     m.reduced_info[-2] = ReducedInfiniteInfo(inf, Dict(2 => 0.5))
     m.infinite_to_reduced[JuMP.index(inf)] = [-2]
     rv2 = ReducedInfiniteVariableRef(m, -2)
     data = DiscreteMeasureData(par, [1], [1])
     meas = measure(inf + par - x + rv, data)
     meas2 = measure(rv2, data)
     @constraint(m, con, x + rv <= 0)
     m.constrs[-1] = ScalarConstraint(rv2, MOI.GreaterThan(0.))
     m.reduced_to_constrs[JuMP.index(rv2)] = [-1]
     con2 = InfiniteConstraintRef(m, Int64(-1), JuMP.shape(m.constrs[-1]))
     set_name(con2, "")
     # test normal deletion
     @test isa(delete(m, rv), Nothing)
     @test name(meas) == "measure(inf(par, par2) + par - x)"
     @test !haskey(m.reduced_to_meas, JuMP.index(rv))
     @test string(m.constrs[JuMP.index(con)].func) == "x"
     @test !haskey(m.reduced_to_constrs, JuMP.index(rv))
     @test m.infinite_to_reduced[JuMP.index(inf)] == [JuMP.index(rv2)]
     @test !haskey(m.reduced_info, JuMP.index(rv))
     # test deletion of special cases
     @test isa(delete(m, rv2), Nothing)
     @test name(meas2) == "measure(0)"
     @test !haskey(m.reduced_to_meas, JuMP.index(rv2))
     @test string(m.constrs[JuMP.index(con2)].func) == "0"
     @test !haskey(m.reduced_to_constrs, JuMP.index(rv2))
     @test !haskey(m.infinite_to_reduced, JuMP.index(inf))
     @test !haskey(m.reduced_info, JuMP.index(rv2))
     # test error
     @test_throws AssertionError delete(m, rv)
     @test_throws AssertionError delete(m, rv2)
 end

 # Test variable deletion
 @testset "JuMP.delete (Global Variables)" begin
     # intialize the model
     m = InfiniteModel()
     @infinite_parameter(m, 0 <= par <= 1)
     @global_variable(m, 0 <= x <= 1, Bin)
     @global_variable(m, y == 1, Int)
     data = DiscreteMeasureData(par, [1], [1])
     meas1 = measure(x + y + par, data)
     meas2 = add_measure(m, Measure(y, data))
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     @objective(m, Min, x + y)
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test !haskey(m.var_to_lower_bound, JuMP.index(x))
     @test !haskey(m.var_to_upper_bound, JuMP.index(x))
     @test !haskey(m.var_to_zero_one, JuMP.index(x))
     @test name(meas1) == "measure(y + par)"
     @test !haskey(m.var_to_meas, JuMP.index(x))
     @test string(m.constrs[JuMP.index(con1)].func) == "y + par"
     @test !haskey(m.var_to_constrs, JuMP.index(x))
     @test string(objective_function(m)) == "y"
     @test !haskey(m.var_in_objective, JuMP.index(x))
     @test !haskey(m.vars, JuMP.index(x))
     @test !haskey(m.var_to_name, JuMP.index(x))
     # test deletion of y
     set_objective_function(m, y)
     @test isa(delete(m, y), Nothing)
     @test !haskey(m.var_to_fix, JuMP.index(y))
     @test !haskey(m.var_to_integrality, JuMP.index(y))
     @test name(meas1) == "measure(par)"
     @test name(meas2) == "measure(0)"
     @test !haskey(m.var_to_meas, JuMP.index(y))
     @test string(m.constrs[JuMP.index(con1)].func) == "par"
     @test string(m.constrs[JuMP.index(con2)].func) == "0"
     @test !haskey(m.var_to_constrs, JuMP.index(y))
     @test string(objective_function(m)) == "0"
     @test !haskey(m.var_in_objective, JuMP.index(y))
     @test !haskey(m.vars, JuMP.index(y))
     @test !haskey(m.var_to_name, JuMP.index(y))
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
     data = DiscreteMeasureData(par, [1], [1])
     meas1 = measure(x + y + par, data)
     meas2 = add_measure(m, Measure(y, data))
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     @objective(m, Min, x + y)
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test !haskey(m.var_to_lower_bound, JuMP.index(x))
     @test !haskey(m.var_to_upper_bound, JuMP.index(x))
     @test !haskey(m.var_to_zero_one, JuMP.index(x))
     @test name(meas1) == "measure(y + par)"
     @test !haskey(m.var_to_meas, JuMP.index(x))
     @test string(m.constrs[JuMP.index(con1)].func) == "y + par"
     @test !haskey(m.var_to_constrs, JuMP.index(x))
     @test string(objective_function(m)) == "y"
     @test m.infinite_to_points[JuMP.index(inf)] == [JuMP.index(y)]
     @test !haskey(m.var_in_objective, JuMP.index(x))
     @test !haskey(m.vars, JuMP.index(x))
     @test !haskey(m.var_to_name, JuMP.index(x))
     # test deletion of y
     set_objective_function(m, y)
     @test isa(delete(m, y), Nothing)
     @test !haskey(m.var_to_fix, JuMP.index(y))
     @test !haskey(m.var_to_integrality, JuMP.index(y))
     @test name(meas1) == "measure(par)"
     @test name(meas2) == "measure(0)"
     @test !haskey(m.var_to_meas, JuMP.index(y))
     @test string(m.constrs[JuMP.index(con1)].func) == "par"
     @test string(m.constrs[JuMP.index(con2)].func) == "0"
     @test !haskey(m.var_to_constrs, JuMP.index(y))
     @test string(objective_function(m)) == "0"
     @test !haskey(m.infinite_to_points, JuMP.index(inf))
     @test !haskey(m.var_in_objective, JuMP.index(y))
     @test !haskey(m.vars, JuMP.index(y))
     @test !haskey(m.var_to_name, JuMP.index(y))
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
     m.reduced_info[-1] = ReducedInfiniteInfo(x, Dict(1 => 0.5))
     m.infinite_to_reduced[JuMP.index(x)] = [-1]
     rv = ReducedInfiniteVariableRef(m, -1)
     data = DiscreteMeasureData(par, [1], [1])
     meas1 = measure(x + y + par, data)
     meas2 = add_measure(m, Measure(y, data))
     @constraint(m, con1, x + y + par <= 0)
     con2 = add_constraint(m, ScalarConstraint(y, MOI.LessThan(0.)))
     # test deletion of x
     @test isa(delete(m, x), Nothing)
     @test !haskey(m.var_to_lower_bound, JuMP.index(x))
     @test !haskey(m.var_to_upper_bound, JuMP.index(x))
     @test !haskey(m.var_to_zero_one, JuMP.index(x))
     @test name(meas1) == "measure(y(par) + par)"
     @test !haskey(m.var_to_meas, JuMP.index(x))
     @test string(m.constrs[JuMP.index(con1)].func) == "y(par) + par"
     @test !haskey(m.var_to_constrs, JuMP.index(x))
     @test m.param_to_vars[JuMP.index(par)] == [JuMP.index(y)]
     @test !is_valid(m, x0)
     @test !haskey(m.infinite_to_points, JuMP.index(x))
     @test !is_valid(m, rv)
     @test !haskey(m.infinite_to_reduced, JuMP.index(x))
     @test !haskey(m.var_in_objective, JuMP.index(x))
     @test !haskey(m.vars, JuMP.index(x))
     @test !haskey(m.var_to_name, JuMP.index(x))
     # test deletion of y
     @test isa(delete(m, y), Nothing)
     @test !haskey(m.var_to_fix, JuMP.index(y))
     @test !haskey(m.var_to_integrality, JuMP.index(y))
     @test name(meas1) == "measure(par)"
     @test name(meas2) == "measure(0)"
     @test !haskey(m.var_to_meas, JuMP.index(y))
     @test string(m.constrs[JuMP.index(con1)].func) == "par"
     @test string(m.constrs[JuMP.index(con2)].func) == "0"
     @test !haskey(m.var_to_constrs, JuMP.index(y))
     @test !haskey(m.param_to_vars, JuMP.index(par))
     @test !haskey(m.var_in_objective, JuMP.index(y))
     @test !haskey(m.vars, JuMP.index(y))
     @test !haskey(m.var_to_name, JuMP.index(y))
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
     @infinite_variable(m, x(par))
     @infinite_variable(m, y(par2))
     @point_variable(m, x(0), x0)
     m.reduced_info[-1] = ReducedInfiniteInfo(x, Dict(1 => 0.5))
     m.infinite_to_reduced[JuMP.index(x)] = [-1]
     rv = ReducedInfiniteVariableRef(m, -1)
     data = DiscreteMeasureData(par, [1], [1])
     meas = measure(x, data)
     meas1 = measure(x + x0 + rv + par + meas + y + par2, data)
     meas2 = measure(x + x0 + rv + par + y + par2, data)
     meas3 = measure(meas1 + x0, data)
     meas4 = measure(meas2, data)
     @constraint(m, con1, x0 + meas1 <= 0)
     con2 = add_constraint(m, ScalarConstraint(meas2, MOI.LessThan(0.)))
     @objective(m, Min, meas1 + x0)
     # test deletion of meas1
     @test isa(delete(m, meas1), Nothing)
     @test name(meas3) == "measure(x0)"
     @test !haskey(m.meas_to_meas, JuMP.index(meas1))
     @test string(m.constrs[JuMP.index(con1)].func) == "x0"
     @test !haskey(m.meas_to_constrs, JuMP.index(meas1))
     @test string(objective_function(m)) == "x0"
     @test m.var_to_meas[JuMP.index(x)] == [JuMP.index(meas), JuMP.index(meas2)]
     @test m.var_to_meas[JuMP.index(y)] == [JuMP.index(meas2)]
     @test m.param_to_meas[JuMP.index(par)] == [1, 3, 4, 5]
     @test m.param_to_meas[JuMP.index(par2)] == [JuMP.index(meas2)]
     @test !haskey(m.meas_to_meas, JuMP.index(meas))
     @test m.reduced_to_meas[JuMP.index(rv)] == [JuMP.index(meas2)]
     @test !haskey(m.meas_in_objective, JuMP.index(meas1))
     @test !haskey(m.measures, JuMP.index(meas1))
     @test !haskey(m.meas_to_name, JuMP.index(meas1))
     # test deletion of meas2
     set_objective_function(m, meas2)
     @test isa(delete(m, meas2), Nothing)
     @test name(meas4) == "measure(0)"
     @test !haskey(m.meas_to_meas, JuMP.index(meas1))
     @test string(m.constrs[JuMP.index(con2)].func) == "0"
     @test !haskey(m.meas_to_constrs, JuMP.index(meas2))
     @test string(objective_function(m)) == "0"
     @test m.var_to_meas[JuMP.index(x)] == [JuMP.index(meas)]
     @test !haskey(m.var_to_meas, JuMP.index(y))
     @test m.param_to_meas[JuMP.index(par)] == [1, 4, 5]
     @test !haskey(m.param_to_meas, JuMP.index(par2))
     @test !haskey(m.reduced_to_meas, JuMP.index(rv))
     @test !haskey(m.meas_in_objective, JuMP.index(meas2))
     @test !haskey(m.measures, JuMP.index(meas2))
     @test !haskey(m.meas_to_name, JuMP.index(meas2))
     # test errors
     @test_throws AssertionError delete(m, meas1)
     @test_throws AssertionError delete(m, meas2)
 end
