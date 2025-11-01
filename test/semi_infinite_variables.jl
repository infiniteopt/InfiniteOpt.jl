# Test basic methods
@testset "Basics" begin
    # initialize model and semi_infinite variable
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1])
    @infinite_parameter(m, b[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, c[1:2] in [0, 1])
    @variable(m, ivref, Infinite(a, b..., c))
    eval_supp = [0.5, NaN, 1., 0., 0.]
    var = SemiInfiniteVariable(RestrictedDomainInfo(), ivref, eval_supp, [2], [2])
    object = VariableData(var, "var")
    idx = SemiInfiniteVariableIndex(1)
    vref = SemiInfiniteVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, SemiInfiniteVariableIndex)
    bad_vref = SemiInfiniteVariableRef(m, SemiInfiniteVariableIndex(-1))
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
        @test isequal(dispatch_variable_ref(m, idx), vref)
        @test isequal(dispatch_variable_ref(gvref), vref)
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(m, SemiInfiniteVariable) === m.semi_infinite_vars
        @test InfiniteOpt._data_dictionary(vref) === m.semi_infinite_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.semi_infinite_vars
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
        @test core_object(vref) === var
        @test core_object(gvref) === var
    end
    # _set_core_object
    @testset "_set_core_object" begin
        @test InfiniteOpt._set_core_object(vref, var) isa Nothing
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(vref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(vref) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == InfOptConstraintIndex[]
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(vref) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(gvref) == DerivativeIndex[]
    end
    # parameter_group_int_indices
    @testset "parameter_group_int_indices" begin
        @test InfiniteOpt.parameter_group_int_indices(vref) == [2]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(vref) == [2]
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test isequal(infinite_variable_ref(vref), ivref)
        @test isequal(infinite_variable_ref(gvref), ivref)
    end
    # eval_support
    @testset "eval_support" begin
        @test eval_support(vref) === eval_supp
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        @test isequal(raw_parameter_refs(vref), IC.VectorTuple(b[1]))
        @test isequal(raw_parameter_refs(gvref), IC.VectorTuple(b[1]))
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test isequal(parameter_refs(vref), (b[1],))
        @test isequal(parameter_refs(gvref), (b[1],))
    end
    # parameter_list
    @testset "parameter_list" begin
        @test isequal(parameter_list(vref), [b[1]])
        @test isequal(parameter_list(gvref), [b[1]])
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        # test normal
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        @test_throws ErrorException set_name(bad_vref, "")
        # test default
        @test isa(set_name(gvref, "a"), Nothing)
        @test name(vref) == "a"
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "a"
        @test name(bad_vref) == ""
    end
    # GeneralVariableRef
    @testset "GeneralVariableRef" begin
        @test isequal(InfiniteOpt.GeneralVariableRef(m, idx), gvref)
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
        @test isequal(variable_by_name(m, "a"), gvref)
        @test variable_by_name(m, "test") isa Nothing
        # prepare variable with same name
        idx2 = SemiInfiniteVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = SemiInfiniteVariableRef(m, idx2)
        @test set_name(vref2, "a") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "a")
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
    # initialize model and semi_infinite variable
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1])
    @infinite_parameter(m, b[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, c[1:2] in [0, 1])
    @variable(m, ivref, Infinite(a, b..., c))
    eval_supp = [0.5, NaN, 1.0, 0.0, 0.0]
    info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
    info2 = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, true)
    info3 = VariableInfo(true, 0.0, true, 0., true, 0., true, 0., false, false)
    # test _process_value 
    @testset "_process_value" begin 
        @test isequal(InfiniteOpt._process_value(0), 0)
        @test isequal(InfiniteOpt._process_value(0a + 1), 1)
        @test isequal(InfiniteOpt._process_value(1a), a)
        @test isequal(InfiniteOpt._process_value(a + 2), a + 2)
    end
    # test SemiInfinite{V, T}
    @testset "SemiInfinite" begin 
        @test isequal(SemiInfinite(ivref, 0, b..., c).infinite_variable_ref, ivref)
        @test SemiInfinite(ivref, 0, b[1], 0, c).parameter_values isa IC.VectorTuple
    end
    # test JuMP.build_variable for macros
    @testset "JuMP.build_variable macro wrapper" begin 
        # test errors
        @test_throws ErrorException build_variable(error, info, 
                                                   SemiInfinite(ivref, 0, b..., c), 
                                                   bad = 42)
        @test_throws ErrorException build_variable(error, info2, 
                                                   SemiInfinite(ivref, 0, b..., a))
        # @test_throws ErrorException build_variable(error, info, SemiInfinite(0))
        @test_throws ErrorException build_variable(error, info, 
                                                   SemiInfinite(ivref, 0, b..., 
                                                   dispatch_variable_ref.(c)))
        @test_throws ErrorException build_variable(error, info, 
                                                   SemiInfinite(ivref, 0, b...))
        @test_throws ErrorException build_variable(error, info, 
                                                   SemiInfinite(ivref, 0, b..., [0, c[2]]))
        # test normal 
        @test isequal(build_variable(error, info, 
                             SemiInfinite(ivref, 0, b..., c)).infinite_variable_ref, ivref)
        @test build_variable(error, info, 
                             SemiInfinite(ivref, 0, b..., c)).eval_support[1] == 0
        @test build_variable(error, info3, SemiInfinite(ivref, 0, b..., c)).info.lower_bound == 0
    end
    # test JuMP.build_variable
    @testset "JuMP.build_variable" begin
        # test errors
        dvref = InfiniteOpt.GeneralVariableRef(m, FiniteVariableIndex(1))
        @test_throws ErrorException build_variable(error, dvref, eval_supp)
        @test_throws ErrorException build_variable(error, ivref, [eval_supp; NaN])
        eval_supp[2] = 3.
        @test_throws ErrorException build_variable(error, ivref, eval_supp)
        eval_supp[2] = NaN
        # test normal
        @test isequal(build_variable(error, ivref, eval_supp).infinite_variable_ref, ivref)
        @test build_variable(error, ivref, eval_supp).eval_support === eval_supp
        new_info = RestrictedDomainInfo(true, 1., false, NaN, false, NaN, true, 42.)
        @test build_variable(error, ivref, eval_supp, new_info).info.lower_bound == 1
        @test build_variable(error, ivref, eval_supp, check = false).group_int_idxs == [2]
        @test build_variable(error, ivref, eval_supp, check = false).parameter_nums == [2]
    end
    # test JuMP.add_variable
    @testset "JuMP.add_variable" begin
        # test errors
        m2 = InfiniteModel()
        var = build_variable(error, ivref, eval_supp, check = false)
        @test_throws VariableNotOwned{InfiniteVariableRef} add_variable(m2, var)
        # test normal
        idx = SemiInfiniteVariableIndex(1)
        vref = SemiInfiniteVariableRef(m, idx)
        gvref = GeneralVariableRef(m, 1, SemiInfiniteVariableIndex)
        @test isequal(add_variable(m, var), gvref)
        @test name(vref) == ""
        @test eval_support(vref) === eval_supp
        @test supports(a) == [0.5]
        @test supports(b[2]) == [1]
        @test supports(c) == zeros(2, 1)
        @test supports(b[1]) == []
        @test InfiniteOpt.parameter_group_int_indices(vref) == [2]
        @test InfiniteOpt._semi_infinite_variable_dependencies(ivref) == [idx]
        # test with set name (redundant add)
        @test isequal(add_variable(m, var, "cat"), gvref)
        @test name(vref) == "cat"
        # test adding info constraints
        var = build_variable(error, info3, SemiInfinite(ivref, 0, b..., c))
        idx = SemiInfiniteVariableIndex(2)
        vref = SemiInfiniteVariableRef(m, idx)
        gvref = GeneralVariableRef(vref)
        @test isequal(add_variable(m, var), gvref)
        @test isequal(eval_support(vref), [0; fill(NaN, 4)])
        @test lower_bound(vref) == 0
        @test upper_bound(vref) == 0
        @test fix_value(vref) == 0
        @test start_value(vref) == 0
        @test LowerBoundRef(vref) isa InfOptConstraintRef
        @test UpperBoundRef(vref) isa InfOptConstraintRef
        @test FixRef(vref) isa InfOptConstraintRef
        # test redundant add that removes constraints
        var = build_variable(error, info, SemiInfinite(ivref, 0, b..., c))
        @test isequal(add_variable(m, var), gvref)
        @test isequal(eval_support(vref), [0; fill(NaN, 4)])
        @test !has_lower_bound(vref)
        @test !has_upper_bound(vref)
        @test !is_fixed(vref)
        @test start_value(vref) isa Nothing
    end
    # test macro definition
    @testset "Macro Definition" begin 
        # anonymous definition
        vref = GeneralVariableRef(m, 3, SemiInfiniteVariableIndex)
        @test isequal(@variable(m, variable_type = SemiInfinite(ivref, 0, b[1], 0, c)), vref)
        @test isequal(parameter_refs(vref), (b[1], c))
        @test eval_support(vref)[[1, 3]] == [0.0, 0.0]
        @test supports(a) == [0, 0.5]
        @test supports(b[2]) == [0, 1]
        @test supports(b[1]) == []
        # explicit definition (redundant)
        vref = GeneralVariableRef(m, 3, SemiInfiniteVariableIndex)
        @test isequal(@variable(m, test, SemiInfinite(ivref, 0, b[1], 0, c)), vref)
        @test isequal(parameter_refs(vref), (b[1], c))
        @test count(isnan.(eval_support(vref))) == 3
        @test name(vref) == "test"
        # array definition
        vrefs = [GeneralVariableRef(m, i, SemiInfiniteVariableIndex) for i in [2, 4]]
        @test isequal(@variable(m, [i = 1:2], SemiInfinite(ivref, i - 1, b..., c)), vrefs)
        @test isequal(parameter_refs(vrefs[1]), (b..., c))
        @test first.(eval_support.(vrefs)) == [0.0, 1.0]
    end
    # test restriciton definition 
    @testset "Restriction" begin 
        # test errors 
        @test_throws ErrorException restrict(ivref, 0, b...)
        @test_throws ErrorException ivref(0, b...)
        # test normal wth restrict
        vref = GeneralVariableRef(m, 5, SemiInfiniteVariableIndex)
        @test isequal(restrict(ivref, 0.5, b[1], 0, c), vref)
        @test isequal(parameter_refs(vref), (b[1], c))
        @test eval_support(vref)[1] == 0.5
        # test normal functionally
        vref = GeneralVariableRef(m, 5, SemiInfiniteVariableIndex)
        @test isequal(ivref(0.5, b[1], 0, c), vref)
        @test isequal(parameter_refs(vref), (b[1], c))
        @test isnan(eval_support(vref)[2])
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    @variable(m, y, Infinite(t, x))
    eval_supp = [0.5, NaN, NaN]
    var = build_variable(error, y, eval_supp, check = false)
    gvref = add_variable(m, var)
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
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test used_by_constraint(gvref)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_derivative
    @testset "used_by_derivative" begin
        @test !used_by_derivative(vref)
        push!(InfiniteOpt._derivative_dependencies(vref), DerivativeIndex(1))
        @test used_by_derivative(gvref)
        @test used_by_derivative(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
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
        push!(InfiniteOpt._constraint_dependencies(vref), InfOptConstraintIndex(1))
        @test is_used(y)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by derivative
        num = 0.
        info = VariableInfo(true, num, true, num, true, num, false, num, true, true)
        deriv = Derivative(info, y, t, 1)
        object = VariableData(deriv)
        idx = DerivativeIndex(1)
        dref = DerivativeRef(m, idx)
        @test InfiniteOpt._add_data_object(m, object) == idx
        push!(InfiniteOpt._derivative_dependencies(vref), idx)
        @test !is_used(vref)
        push!(InfiniteOpt._constraint_dependencies(dref), InfOptConstraintIndex(2))
        @test is_used(vref)
        empty!(InfiniteOpt._derivative_dependencies(vref))
    end
end

# Test variable queries for models
@testset "Model Queries" begin
    # initialize model, parameter, and variables
    m = InfiniteModel()
    @infinite_parameter(m, a in [0, 1])
    @infinite_parameter(m, b[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, c[1:2] in [0, 1])
    @variable(m, ivref, Infinite(a, b..., c))
    @variable(m, pvref, Point(ivref, 0, 0, 0, [0, 0]))
    @variable(m, hvref)
    eval_supp = [0.5, NaN, 1, 0, 0]
    var = build_variable(error, ivref, eval_supp, check = false)
    rvref = add_variable(m, var)
    # num_variables
    @testset "JuMP.num_variables" begin
        @test num_variables(m, InfiniteVariable) == 1
        @test num_variables(m, PointVariable) == 1
        @test num_variables(m, FiniteVariable) == 1
        @test num_variables(m, SemiInfiniteVariable) == 1
        @test num_variables(m) == 4
    end
    # all_variables
    @testset "JuMP.all_variables" begin
        @test all_variables(m, InfiniteVariable) == [ivref]
        @test all_variables(m, PointVariable) == [pvref]
        @test all_variables(m, FiniteVariable) == [hvref]
        @test all_variables(m, SemiInfiniteVariable) == [rvref]
        @test all_variables(m) == [ivref, rvref, pvref, hvref]
    end
end
