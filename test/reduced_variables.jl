# Test basic methods
@testset "Basics" begin
    # initialize model and reduced variable
    m = InfiniteModel()
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])
    @infinite_variable(m, ivref(a, b, c))
    eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1, 4 => 0, 5 => 0)
    var = ReducedVariable(ivref, eval_supps, [2])
    object = VariableData(var, "var")
    idx = ReducedVariableIndex(1)
    vref = ReducedVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, ReducedVariableIndex)
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
        @test InfiniteOpt._data_dictionary(m, ReducedVariable) === m.reduced_vars
        @test InfiniteOpt._data_dictionary(vref) === m.reduced_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.reduced_vars
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
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) === var
        @test InfiniteOpt._core_variable_object(gvref) === var
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(vref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(vref) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == ConstraintIndex[]
    end
    # _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(vref) == [2]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(vref) == [2]
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(vref) == ivref
        @test infinite_variable_ref(gvref) == ivref
    end
    # eval_supports
    @testset "eval_supports" begin
        @test eval_supports(vref) === eval_supps
        @test eval_supports(gvref) === eval_supps
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        @test raw_parameter_refs(vref) == VectorTuple(b[1])
        @test raw_parameter_refs(gvref) == VectorTuple(b[1])
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == (b[1],)
        @test parameter_refs(gvref) == (b[1],)
    end
    # parameter_list
    @testset "parameter_list" begin
        @test parameter_list(vref) == [b[1]]
        @test parameter_list(gvref) == [b[1]]
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        # test default
        @test isa(set_name(gvref, ""), Nothing)
        @test name(vref) == "ivref(0.5, [b[1], 1], [0, 0])"
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "ivref(0.5, [b[1], 1], [0, 0])"
        InfiniteOpt._data_object(vref).name = ""
        @test name(gvref) == "ivref(0.5, [b[1], 1], [0, 0])"
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
        @test variable_by_name(m, "ivref(0.5, [b[1], 1], [0, 0])") == gvref
        @test variable_by_name(m, "test(0.5, [b[1], 1], [0, 0])") isa Nothing
        # prepare variable with same name
        idx2 = ReducedVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = ReducedVariableRef(m, idx2)
        @test set_name(vref2, "") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "ivref(0.5, [b[1], 1], [0, 0])")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(vref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(vref)) == 1
        @test !is_valid(m, vref)
    end
end

# Test definition methods
@testset "Definition" begin
    # initialize model and reduced variable
    m = InfiniteModel()
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])
    @infinite_variable(m, ivref(a, b, c))
    eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1, 4 => 0, 5 => 0)
    # test JuMP.build_variable
    @testset "JuMP.build_variable" begin
        # test errors
        dvref = InfiniteOpt._make_variable_ref(m, HoldVariableIndex(1))
        @test_throws ErrorException build_variable(error, dvref, eval_supps)
        eval_supps[6] = 1
        @test_throws ErrorException build_variable(error, ivref, eval_supps)
        delete!(eval_supps, 6)
        eval_supps[4] = 2
        @test_throws ErrorException build_variable(error, ivref, eval_supps)
        eval_supps[4] = 0
        # test normal
        @test build_variable(error, ivref, eval_supps).infinite_variable_ref == ivref
        @test build_variable(error, ivref, eval_supps).eval_supports === eval_supps
        @test build_variable(error, ivref, eval_supps).object_nums == [2]
        @test build_variable(error, ivref, eval_supps, check = false).object_nums == [2]
    end
    # test JuMP.add_variable
    @testset "JuMP.add_variable" begin
        # test errors
        m2 = InfiniteModel()
        var = build_variable(error, ivref, eval_supps, check = false)
        @test_throws VariableNotOwned{InfiniteVariableRef} add_variable(m2, var)
        # test normal
        idx = ReducedVariableIndex(1)
        vref = ReducedVariableRef(m, idx)
        gvref = GeneralVariableRef(m, 1, ReducedVariableIndex)
        @test add_variable(m, var) == gvref
        @test name(vref) == "ivref(0.5, [b[1], 1], [0, 0])"
        @test eval_supports(vref) === eval_supps
        @test InfiniteOpt._object_numbers(vref) == [2]
        @test InfiniteOpt._reduced_variable_dependencies(ivref) == [idx]
        # test with set name
        idx = ReducedVariableIndex(2)
        vref = ReducedVariableRef(m, idx)
        gvref = GeneralVariableRef(m, 2, ReducedVariableIndex)
        @test add_variable(m, var, "cat") == gvref
        @test name(vref) == "cat"
        # test with no name
        idx = ReducedVariableIndex(3)
        vref = ReducedVariableRef(m, idx)
        gvref = GeneralVariableRef(m, 3, ReducedVariableIndex)
        @test add_variable(m, var, define_name = false) == gvref
        @test InfiniteOpt._data_object(vref).name == ""
    end
end

# Test the info methods
@testset "Info" begin
    # initialize model and reduced variable
    m = InfiniteModel()
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])
    @infinite_variable(m,  0 <= ivref1(a, b, c) <= 1, Int)
    @infinite_variable(m,  ivref2(a, b, c) == 1, Bin, start = 0)
    eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1, 4 => 0, 5 => 0)
    var1 = build_variable(error, ivref1, eval_supps, check = false)
    gvref1 = add_variable(m, var1, define_name = false)
    rvref1 = dispatch_variable_ref(gvref1)
    var2 = build_variable(error, ivref2, eval_supps, check = false)
    gvref2 = add_variable(m, var2, define_name = false)
    rvref2 = dispatch_variable_ref(gvref2)
    dvref1 = dispatch_variable_ref(ivref1)
    dvref2 = dispatch_variable_ref(ivref2)
    # has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(rvref1)
        @test has_lower_bound(gvref1)
        @test !has_lower_bound(rvref2)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(rvref1) == 0
        @test lower_bound(gvref1) == 0
        @test_throws ErrorException lower_bound(rvref2)
    end
    # _lower_bound_index
    @testset "JuMP._lower_bound_index" begin
        @test JuMP._lower_bound_index(rvref1) == JuMP._lower_bound_index(dvref1)
        @test_throws ErrorException JuMP._lower_bound_index(rvref2)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        @test LowerBoundRef(rvref1) == LowerBoundRef(dvref1)
        @test_throws ErrorException LowerBoundRef(rvref2)
        @test_throws ErrorException LowerBoundRef(gvref2)
    end
    # has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(rvref1)
        @test has_upper_bound(gvref1)
        @test !has_upper_bound(rvref2)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(rvref1) == 1
        @test upper_bound(gvref1) == 1
        @test_throws ErrorException upper_bound(rvref2)
    end
    # _upper_bound_index
    @testset "JuMP._upper_bound_index" begin
        @test JuMP._upper_bound_index(rvref1) == JuMP._upper_bound_index(dvref1)
        @test_throws ErrorException JuMP._upper_bound_index(rvref2)
    end
    # UpperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        @test UpperBoundRef(rvref1) == UpperBoundRef(dvref1)
        @test_throws ErrorException UpperBoundRef(rvref2)
        @test_throws ErrorException UpperBoundRef(gvref2)
    end
    # is_fixed
    @testset "JuMP.is_fixed" begin
        @test is_fixed(rvref2)
        @test !is_fixed(rvref1)
        @test !is_fixed(gvref1)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(rvref2) == 1
        @test_throws ErrorException fix_value(rvref1)
        @test_throws ErrorException fix_value(gvref1)
    end
    # _fix_index
    @testset "JuMP._fix_index" begin
        @test JuMP._fix_index(rvref2) == JuMP._fix_index(dvref2)
        @test_throws ErrorException JuMP._fix_index(rvref1)
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        @test FixRef(rvref2) == FixRef(dvref2)
        @test FixRef(rvref2) == FixRef(dvref2)
        @test_throws ErrorException FixRef(rvref1)
    end
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(rvref2) == 0
        @test start_value(gvref2) == 0
    end
    # is_binary
    @testset "JuMP.is_binary" begin
        @test is_binary(rvref2)
        @test is_binary(gvref2)
        @test !is_binary(rvref1)
    end
    # _binary_index
    @testset "JuMP._binary_index" begin
        @test JuMP._binary_index(rvref2) == JuMP._binary_index(dvref2)
        @test_throws ErrorException JuMP._binary_index(rvref1)
    end
    # BinaryRef
    @testset "JuMP.BinaryRef" begin
        @test BinaryRef(rvref2) == BinaryRef(dvref2)
        @test BinaryRef(gvref2) == BinaryRef(dvref2)
        @test_throws ErrorException BinaryRef(rvref1)
    end
    # is_integer
    @testset "JuMP.is_integer" begin
        @test is_integer(rvref1)
        @test !is_integer(rvref2)
        @test !is_integer(gvref2)
    end
    # _integer_index
    @testset "JuMP._integer_index" begin
        @test JuMP._integer_index(rvref1) == JuMP._integer_index(dvref1)
        @test_throws ErrorException JuMP._integer_index(rvref2)
    end
    # IntegerRef
    @testset "JuMP.IntegerRef" begin
        @test IntegerRef(rvref1) == IntegerRef(dvref1)
        @test IntegerRef(gvref1) == IntegerRef(dvref1)
        @test_throws ErrorException IntegerRef(rvref2)
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @independent_parameter(m, t in [0, 1])
    @dependent_parameters(m, x[1:2] in [-1, 1])
    @infinite_variable(m, y(t, x))
    eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1)
    var = build_variable(error, y, eval_supps, check = false)
    gvref = add_variable(m, var, define_name = false)
    vref = dispatch_variable_ref(gvref)
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(gvref)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test used_by_constraint(gvref)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(gvref)
        @test !used_by_objective(vref)
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(gvref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test is_used(y)
    end
end

# Test variable queries for models
@testset "Model Queries" begin
    # initialize model, parameter, and variables
    m = InfiniteModel()
    @independent_parameter(m, a in [0, 1])
    @independent_parameter(m, b[1:2] in [0, 1])
    @dependent_parameters(m, c[1:2] in [0, 1])
    @infinite_variable(m, ivref(a, b, c))
    @point_variable(m, ivref(0, [0, 0], [0, 0]), pvref)
    @hold_variable(m, hvref)
    eval_supps = Dict{Int, Float64}(1 => 0.5, 3 => 1, 4 => 0, 5 => 0)
    var = build_variable(error, ivref, eval_supps, check = false)
    rvref = add_variable(m, var, define_name = false)
    # num_variables
    @testset "JuMP.num_variables" begin
        @test num_variables(m, InfiniteVariable) == 1
        @test num_variables(m, PointVariable) == 1
        @test num_variables(m, HoldVariable) == 1
        @test num_variables(m, ReducedVariable) == 1
        @test num_variables(m) == 4
        @test num_variables(m, InfOptVariable) == 4
    end
    # all_variables
    @testset "JuMP.all_variables" begin
        @test all_variables(m, InfiniteVariable) == [ivref]
        @test all_variables(m, PointVariable) == [pvref]
        @test all_variables(m, HoldVariable) == [hvref]
        @test all_variables(m, ReducedVariable) == [rvref]
        @test all_variables(m) == [ivref, rvref, pvref, hvref]
        @test all_variables(m, InfOptVariable) == [ivref, rvref, pvref, hvref]
    end
end
