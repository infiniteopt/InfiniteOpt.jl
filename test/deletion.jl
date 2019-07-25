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
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
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
end
