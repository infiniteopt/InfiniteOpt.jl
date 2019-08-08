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

# Test parameter deletion
@testset "Parameters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par, pars))
    @infinite_variable(m, inf4(par, par2, pars))
    @point_variable(m, inf(0.5), pt)
    pt2 = @point_variable(m, inf2(0.5, 0.5))
    @global_variable(m, x)
    rv = ReducedInfiniteVariableRef(m, -1)
    m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    meas = measure(inf + par - x + rv, data)
    meas2 = measure(inf3 - x + pt2, data2)

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
        struct bob end
        @test_throws ErrorException InfiniteOpt._check_param_in_data(par, bob())
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
        new_pars = JuMP.Containers.SparseAxisArray(copy(pars.data))
        filter!(x -> x.second != new_pars[1], new_pars.data)
        pars_copy = JuMP.Containers.SparseAxisArray(copy(pars.data))
        @test InfiniteOpt._remove_parameter((pars_copy,),
                                            pars[1]) == ((new_pars,), (1, (1,)))
        pars_copy = JuMP.Containers.SparseAxisArray(copy(pars.data))
        @test InfiniteOpt._remove_parameter((par, pars_copy),
                                        pars[1]) == ((par, new_pars), (2, (1,)))
        # test removing from array that becomes empty
        @test InfiniteOpt._remove_parameter((new_pars,), pars[2]) == ((), (1,))
        new_pars = JuMP.Containers.SparseAxisArray(copy(pars.data))
        filter!(x -> x.second != new_pars[1], new_pars.data)
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
        arr2 = JuMP.Containers.SparseAxisArray(copy(arr.data))
        filter!(x -> x.first != (2,), arr2.data)
        @test InfiniteOpt._remove_parameter_values((arr, 3), (1, (2,))) == (arr2, 3)
        filter!(x -> x.first != (1,), arr2.data)
        @test InfiniteOpt._remove_parameter_values((arr,), (1, (1,))) == (arr2,)
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
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5))
        @test isa(InfiniteOpt._update_infinite_variable(inf4, (par, par2, pars)),
                                                        Nothing)
        # prepare for removing array element
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0])
        m.reduced_info[-1] = ReducedInfiniteInfo(inf4, Dict(2 => 0.5, 3 => supp))
        new_supp = JuMP.Containers.SparseAxisArray(copy(supp.data))
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
        # TODO test
    end
end
