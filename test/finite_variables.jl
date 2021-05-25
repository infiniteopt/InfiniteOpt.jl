# Test name methods
@testset "Basics" begin
    # initialize model and variable
    m = InfiniteModel()
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, 0., true, false)
    var = ScalarVariable(info)
    object = VariableData(var, "test")
    idx = FiniteVariableIndex(1)
    vref = FiniteVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, FiniteVariableIndex)
    bad_vref = FiniteVariableRef(m, FiniteVariableIndex(-1))
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(vref) === m
        @test owner_model(gvref) === m
    end
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(vref) == idx
        @test index(gvref) == idx
    end
    # dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == vref
        @test dispatch_variable_ref(gvref) == vref
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(m, FiniteVariable) === m.finite_vars
        @test InfiniteOpt._data_dictionary(vref) === m.finite_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.finite_vars
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, vref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(vref) === object
        @test InfiniteOpt._data_object(gvref) === object
        @test_throws ErrorException InfiniteOpt._data_object(bad_vref)
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) === var
        @test InfiniteOpt._core_variable_object(gvref) === var
    end
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(vref) == info
    end
    # _update_variable_info
    @testset "_update_variable_info" begin
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
        @test name(gvref) == "test"
        @test name(bad_vref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        @test isa(set_name(gvref, "new2"), Nothing)
        @test name(vref) == "new2"
        @test_throws ErrorException set_name(bad_vref, "")
    end
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == ()
        @test parameter_refs(gvref) == ()
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, idx) == gvref
    end
    # _var_name_dict
    @testset "_var_name_dict" begin
        @test InfiniteOpt._var_name_dict(m) isa Nothing
    end
    # _update_var_name_dict
    @testset "_update_var_name_dict" begin
        m.name_to_var = Dict{String, ObjectIndex}()
        dict = InfiniteOpt._data_dictionary(vref)
        @test InfiniteOpt._update_var_name_dict(m, dict) isa Nothing
        @test InfiniteOpt._var_name_dict(m) == Dict(name(vref) => idx)
        m.name_to_var = nothing
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new2") == gvref
        @test variable_by_name(m, "bob") isa Nothing
        # prepare variable with same name
        idx2 = FiniteVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = FiniteVariableRef(m, idx2)
        @test set_name(vref2, "new2") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "new2")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(vref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(vref)) == 1
        @test !is_valid(m, vref)
    end
end

# Test variable definition methods
@testset "Definition" begin
    # initialize model and info
    m = InfiniteModel()
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, num, false, true)
    # build_variable
    @testset "JuMP.build_variable" begin
        # test normal
        @test build_variable(error, info).info == info
    end
    # _check_and_make_variable_ref
    @testset "_check_and_make_variable_ref" begin
        # test normal
        v = build_variable(error, info)
        idx = FiniteVariableIndex(1)
        vref = FiniteVariableRef(m, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m, v, "") == vref
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        idx = FiniteVariableIndex(2)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        v = build_variable(error, info)
        @test add_variable(m, v, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test name(vref) == "name"
        # prepare variable with all the possible info additions
        v = build_variable(error, info2)
        # test info addition functions
        idx = FiniteVariableIndex(3)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        # lower bound
        cindex = InfOptConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex)
        @test has_lower_bound(vref)
        @test InfiniteOpt._lower_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = InfOptConstraintIndex(2)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test InfiniteOpt._upper_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = InfOptConstraintIndex(3)
        cref = InfOptConstraintRef(m, cindex)
        @test is_fixed(vref)
        @test InfiniteOpt._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = InfOptConstraintIndex(4)
        cref = InfOptConstraintRef(m, cindex)
        @test is_binary(vref)
        @test InfiniteOpt._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [InfOptConstraintIndex(i)
                                                             for i = 1:4]
    end
end

# Test the finite variable macro if no bounds are provided
@testset "Macro" begin
    # initialize model
    m = InfiniteModel()
    # test the deprecations 
    @testset "@hold_variable" begin 
        @test_macro_throws ErrorException @hold_variable(m, w)
    end
    @testset "@variable" begin
        # test regular
        idx = FiniteVariableIndex(1)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @variable(m, x >= 1, Bin) == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test anan
        idx = FiniteVariableIndex(2)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @variable(m, binary = true, lower_bound = 1,
                        base_name = "x") == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test array
        idxs = [FiniteVariableIndex(3), FiniteVariableIndex(4)]
        vrefs = [FiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @variable(m, y[1:2] == 2, Int) == gvrefs
        @test name(vrefs[1]) == "y[1]"
        @test fix_value(vrefs[2]) == 2
        @test is_integer(vrefs[1])
        # test errors
        @test_macro_throws ErrorException @variable(m, x >= 1, Bin)
    end
end

# Test variable(s) constrained on creation 
@testset "Creation Constraints" begin 
    # initialize model and stuff
    m = InfiniteModel()
    info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
    v1 = build_variable(error, info)
    # test add_variable with VariableConstrainedOnCreation
    @testset "Adding VariableConstrainedOnCreation" begin 
        v = VariableConstrainedOnCreation(v1, MOI.Integer())
        vref = GeneralVariableRef(m, 1, FiniteVariableIndex)
        @test add_variable(m, v, "test") == vref
        @test name(vref) == "test"
        @test num_constraints(m, MOI.Integer) == 1
    end
    # test add_variable with VariablesConstrainedOnCreation
    @testset "Adding VariablesConstrainedOnCreation" begin 
        v = VariablesConstrainedOnCreation([v1, v1, v1], MOI.SecondOrderCone(3), 
                                           VectorShape())
        vrefs = [GeneralVariableRef(m, i, FiniteVariableIndex) for i in 2:4]
        names = ["test[$i]" for i in 1:3]
        @test add_variable(m, v, names) == vrefs
        @test name.(vrefs) == names
        @test num_constraints(m, MOI.SecondOrderCone) == 1
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @variable(m, y0)
    vref = dispatch_variable_ref(y0)
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(y0)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test used_by_constraint(y0)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(y0)
        @test !used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = true
        @test used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = false
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test is_used(y0)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by objective
        InfiniteOpt._data_object(vref).in_objective = true
        @test is_used(vref)
    end
end

# test deprecations
@testset "Parameter Bound Deprecations" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, z)
    rs = DomainRestrictions(t => 0)
    @test_throws ErrorException has_parameter_bounds(z)
    @test_throws ErrorException parameter_bounds(z)
    @test_throws ErrorException set_parameter_bounds(z, rs)
    @test_throws ErrorException add_parameter_bounds(z, rs)
    @test_throws ErrorException delete_parameter_bounds(z)
end
