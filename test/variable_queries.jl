# Test if used
@testset "Used" begin
    # initialize model, parameter, and variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    v = build_variable(error, info, Global)
    vref = add_variable(m, v, "name")
    # used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(vref)
        # prepare use case
        m.var_to_constrs[index(vref)] = [1]
        # test used
        @test used_by_constraint(vref)
        # undo changes
        delete!(m.var_to_constrs, index(vref))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(vref)
        # prepare use case
        m.var_to_meas[index(vref)] = [1]
        # test used
        @test used_by_measure(vref)
        # undo changes
        delete!(m.var_to_meas, index(vref))
    end
    # used_by_objective
    @testset "used_by_objective" begin
        # test not used
        @test !used_by_objective(vref)
        # prepare use case
        m.var_in_objective[index(vref)] = true
        # test used
        @test used_by_objective(vref)
        # undo changes
        m.var_in_objective[index(vref)] = false
    end
    # is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # prepare use case
        m.var_to_constrs[index(vref)] = [1]
        # test used
        @test is_used(vref)
        # undo changes
        delete!(m.var_to_constrs, index(vref))
        # prepare use case
        m.var_to_meas[index(vref)] = [1]
        # test used
        @test is_used(vref)
        # undo changes
        delete!(m.var_to_meas, index(vref))
        # prepare use case
        m.var_in_objective[index(vref)] = true
        # test used
        @test is_used(vref)
        # undo changes
        m.var_in_objective[index(vref)] = false
    end
    # initialize infinite variable to test infinite methods
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    v = build_variable(error, info, Infinite, parameter_refs = (pref, ))
    vref = add_variable(m, v, "name")
    # used_by_point_variable
    @testset "used_by_point_variable" begin
        # test not used
        @test !used_by_point_variable(vref)
        # build point variable
        v = build_variable(error, info, Point, infinite_variable_ref = vref,
                           parameter_values = 0.5)
        pvref = add_variable(m, v, "name2")
        m.var_in_objective[index(pvref)] = true
        # test used
        @test used_by_point_variable(vref)
        # undo changes
        delete!(m.vars, index(pvref))
        m.var_in_objective[index(pvref)] = false
    end
    # is_used (infinite variable)
    @testset "is_used (infinite var)" begin
        # test not used
        @test !is_used(vref)
        # prepare use case
        m.var_to_constrs[index(vref)] = [1]
        # test used
        @test is_used(vref)
        # undo changes
        delete!(m.var_to_constrs, index(vref))
        # prepare use case
        m.var_to_meas[index(vref)] = [1]
        # test used
        @test is_used(vref)
        # undo changes
        delete!(m.var_to_meas, index(vref))
        # prepare use case
        v = build_variable(error, info, Point, infinite_variable_ref = vref,
                           parameter_values = 0.5)
        pvref = add_variable(m, v, "name2")
        m.var_in_objective[index(pvref)] = true
        # test used
        @test is_used(vref)
    end
end

# Test variable queries for models
@testset "Model Queries" begin
    # initialize model, parameter, and variables
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    v = build_variable(error, info, Infinite, parameter_refs = (pref, ))
    ivref = add_variable(m, v, "name")
    v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                       parameter_values = 0.5)
    pvref = add_variable(m, v, "name2")
    v = build_variable(error, info, Global)
    gvref = add_variable(m, v, "name")
    # num_variables
    @testset "JuMP.num_variables" begin
        @test num_variables(m) == 3
    end
    # all_variables
    @testset "JuMP.all_variables" begin
        @test all_variables(m) == [ivref, pvref, gvref]
    end
end

# Test other variable queries
@testset "Other" begin
    # initialize model, parameter, and variables
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    pref2 = add_parameter(m, param, "θ")
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    v = build_variable(error, info, Infinite, parameter_refs = (pref, ))
    ivref = add_variable(m, v, "name")
    v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                       parameter_values = 0.5)
    pvref = add_variable(m, v, "name2")
    # parameter_values
    @testset "parameter_values" begin
        @test parameter_values(pvref) == (0.5, )
    end
    # _update_variable_param_refs
    @testset "_update_variable_param_refs" begin
        @test isa(InfiniteOpt._update_variable_param_refs(ivref, (pref2, )),
                  Nothing)
        @test parameter_refs(ivref) == (pref2, )
    end
    # set_parameter_refs
    @testset "set_parameter_refs" begin
        # test normal with 1 param
        @test isa(set_parameter_refs(ivref, (pref2, )), Nothing)
        @test parameter_refs(ivref) == (pref2, )
        @test name(ivref) == "name(θ)"
        # test double specify
        @test_throws ErrorException set_parameter_refs(ivref, (pref2, pref2))
        # test with used variable
        m.var_to_meas[index(ivref)] = [1]
        @test isa(set_parameter_refs(ivref, (pref, )), Nothing)
        @test parameter_refs(ivref) == (pref, )
        @test name(ivref) == "name(test)"
        delete!(m.var_to_meas, index(ivref))
    end
    # add_parameter_ref
    @testset "add_parameter_ref" begin
        # test normal use
        @test isa(add_parameter_ref(ivref, pref2), Nothing)
        @test parameter_refs(ivref) == (pref, pref2)
        @test name(ivref) == "name(test, θ)"
        # test duplication error
        @test_throws ErrorException add_parameter_ref(ivref, pref2)
    end
end
