# Test Core data accessors
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    ind_idx = IndependentParameterIndex(1)
    fin_idx = FiniteParameterIndex(1)
    set = IntervalSet(0, 1)
    supps_dict = SortedDict{Float64, Set{Symbol}}(0. => Set{Symbol}([UserDefined]))
    ind_param = IndependentParameter(set, supps_dict)
    fin_param = FiniteParameter(42)
    ind_object = ScalarParameterData(ind_param, 1, 1, "ind")
    fin_object = ScalarParameterData(fin_param, -1, -1, "fin")
    ind_pref = IndependentParameterRef(m, ind_idx)
    fin_pref = FiniteParameterRef(m, fin_idx)
    ind_gvref = GeneralVariableRef(m, 1, IndependentParameterIndex)
    fin_gvref = GeneralVariableRef(m, 1, FiniteParameterIndex)
    # test dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, ind_idx) == ind_pref
        @test dispatch_variable_ref(ind_gvref) == ind_pref
        @test dispatch_variable_ref(m, fin_idx) == fin_pref
        @test dispatch_variable_ref(fin_gvref) == fin_pref
    end
    # test _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, ind_object) == ind_idx
        @test InfiniteOpt._add_data_object(m, fin_object) == fin_idx
        @test InfiniteOpt._param_object_indices(m)[end] == ind_idx
    end
    # test _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(ind_pref) === m.independent_params
        @test InfiniteOpt._data_dictionary(ind_gvref) === m.independent_params
        @test InfiniteOpt._data_dictionary(fin_pref) === m.finite_params
        @test InfiniteOpt._data_dictionary(fin_gvref) === m.finite_params
        @test InfiniteOpt._data_dictionary(m, IndependentParameter) === m.independent_params
        @test InfiniteOpt._data_dictionary(m, FiniteParameter) === m.finite_params
    end
    # test _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(ind_pref) == ind_object
        @test InfiniteOpt._data_object(ind_gvref) == ind_object
        @test InfiniteOpt._data_object(fin_pref) == fin_object
        @test InfiniteOpt._data_object(fin_gvref) == fin_object
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(ind_pref) == ind_param
        @test InfiniteOpt._core_variable_object(ind_gvref) == ind_param
        @test InfiniteOpt._core_variable_object(fin_pref) == fin_param
        @test InfiniteOpt._core_variable_object(fin_gvref) == fin_param
    end
    # test _parameter_number
    @testset "_parameter_number" begin
        @test InfiniteOpt._parameter_number(ind_pref) == 1
        @test InfiniteOpt._parameter_number(ind_gvref) == 1
    end
    # test _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(ind_pref) == [1]
        @test InfiniteOpt._parameter_numbers(ind_gvref) == [1]
    end
    # test _object_number
    @testset "_object_number" begin
        @test InfiniteOpt._object_number(ind_pref) == 1
        @test InfiniteOpt._object_number(ind_gvref) == 1
    end
    # test _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(ind_pref) == [1]
        @test InfiniteOpt._object_numbers(ind_gvref) == [1]
    end
    # test _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(ind_pref, ind_param) isa Nothing
        @test InfiniteOpt._set_core_variable_object(fin_pref, fin_param) isa Nothing
    end
end

# Test macro methods
@testset "Macro Helpers" begin
    @testset "Symbol Methods" begin
        @test InfiniteOpt._is_set_keyword(:(lower_bound = 0))
    end
    # test ParameterInfoExpr datatype
    @testset "_ParameterInfoExpr" begin
        @test InfiniteOpt._ParameterInfoExpr isa DataType
        @test InfiniteOpt._ParameterInfoExpr(ones(Bool, 8)...).has_lb
        @test !InfiniteOpt._ParameterInfoExpr().has_lb
    end
    # test JuMP._set_lower_bound_or_error
    @testset "JuMP._set_lower_bound_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        # test normal operation
        @test isa(JuMP._set_lower_bound_or_error(error, info, 0), Nothing)
        @test info.has_lb && info.lower_bound == 0
        # test double/lack of input errors
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
        info.has_lb = false; info.has_dist = true
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
        info.has_dist = false; info.has_set = true
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
    end
    # test JuMP._set_upper_bound_or_error
    @testset "JuMP._set_upper_bound_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        # test normal operation
        @test isa(JuMP._set_upper_bound_or_error(error, info, 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        # test double/lack of input errors
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
        info.has_ub = false; info.has_dist = true
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
        info.has_dist = false; info.has_set = true
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
    end
    # _dist_or_error
    @testset "_dist_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        # test normal operation
        @test isa(InfiniteOpt._dist_or_error(error, info, 0), Nothing)
        @test info.has_dist && info.distribution == 0
        # test double/lack of input errors
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
        info.has_dist = false; info.has_lb = true
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
        info.has_lb = false; info.has_set = true
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
    end
    # _set_or_error
    @testset "_set_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        # test normal operation
        @test isa(InfiniteOpt._set_or_error(error, info, 0), Nothing)
        @test info.has_set && info.set == 0
        # test double/lack of input errors
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
        info.has_set = false; info.has_lb = true
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
        info.has_lb = false; info.has_dist = true
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
    end
    # _constructor_set
    @testset "InfiniteOpt._constructor_set" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test_throws ErrorException InfiniteOpt._constructor_set(error, info)
        info.has_lb = true; info.lower_bound = 0
        @test_throws ErrorException InfiniteOpt._constructor_set(error, info)
        info.has_ub = true; info.upper_bound = 1
        check = :(isa($(info.lower_bound), Real))
        expected = :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a number."))
        @test InfiniteOpt._constructor_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(distribution = Normal())
        check = :(isa($(info.distribution), Distributions.UnivariateDistribution))
        expected = :($(check) ? UniDistributionSet($(info.distribution)) : error("Distribution must be a Distributions.UnivariateDistribution."))
        @test InfiniteOpt._constructor_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(set = IntervalSet(0, 1))
        check1 = :(isa($(info.set), InfiniteScalarSet))
        check2 = :(isa($(info.set), Distributions.UnivariateDistribution))
        expected = :($(check1) ? $(info.set) : ($(check2) ? UniDistributionSet($(info.set)) : error("Set must be a subtype of InfiniteScalarSet.")))
        @test InfiniteOpt._constructor_set(error, info) == expected
    end

    # _parse_one_operator_parameter
    @testset "_parse_one_operator_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:<=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:>=), 0), Nothing)
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:in), esc(0)), Nothing)
        @test info.has_set && info.set == esc(0)
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:in), esc(:([0, 1]))), Nothing)
        @test info.has_lb && info.has_ub
        @test_throws ErrorException InfiniteOpt._parse_one_operator_parameter(error, info,
                                                                              Val(:d), 0)
    end
    # _parse_ternary_parameter
    @testset "_parse_ternary_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_ternary_parameter(error, info, Val(:<=), 0,
                                                       Val(:<=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_ternary_parameter(error, info, Val(:>=), 0,
                                                       Val(:>=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
        @test_throws ErrorException InfiniteOpt._parse_ternary_parameter(error,
                                                 info, Val(:<=), 0, Val(:>=), 0)
    end
    # _parse_parameter
    @testset "_parse_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test InfiniteOpt._parse_parameter(error, info, :<=, :x, 0) == :x
        @test info.has_ub && info.upper_bound == 0
        @test InfiniteOpt._parse_parameter(error, info, :<=, 0, :x) == :x
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test InfiniteOpt._parse_parameter(error, info, 0, :<=, :x, :<=, 0) == :x
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
    end
end

# Test parameter definition methods
@testset "Definition" begin
    # _check_supports_in_bounds
    @testset "_check_supports_in_bounds" begin
        set = IntervalSet(0, 1)
        @test isa(InfiniteOpt._check_supports_in_bounds(error, 0, set), Nothing)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                        -1, set)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                         2, set)
        set = UniDistributionSet(Uniform())
        @test isa(InfiniteOpt._check_supports_in_bounds(error, 0, set), Nothing)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                        -1, set)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                         2, set)
    end
    # build_independent_parameter
    @testset "build_parameter (IndependentParameter)" begin
        set = IntervalSet(0, 1)
        supps = 0.
        supps_dict = SortedDict{Float64, Set{Symbol}}(0. => Set([UserDefined]))
        @test build_parameter(error, set, supports = supps).set == set
        @test build_parameter(error, set, supports = supps).supports == supps_dict
        @test_throws ErrorException build_parameter(error, set, bob = 42)
        warn = "Ignoring num_supports since supports is not empty."
        @test_logs (:warn, warn) build_parameter(error, set,
                                            supports = [0, 1], num_supports = 2)
        repeated_supps = [1, 1]
        expected = IndependentParameter(set, SortedDict{Float64, Set{Symbol}}(1. => Set{Symbol}()))
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) build_parameter(error, set, supports = repeated_supps) == expected
        set = UniDistributionSet(Normal())
        @test length(build_parameter(error, set, num_supports = 5).supports) == 5
    end
    # build_finite_parameter
    @testset "build_parameter (FiniteParameter)" begin
        @test_throws ErrorException build_parameter(error, 1, bob = 42)
        expected = FiniteParameter(1)
        @test build_parameter(error, 1) == expected
    end

    # add_parameter
    @testset "add_parameter" begin
        m = InfiniteModel()
        param = IndependentParameter(IntervalSet(0, 1),
                                             SortedDict{Float64, Set{Symbol}}())
        expected = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
        @test add_parameter(m, param) == expected
        @test InfiniteOpt._core_variable_object(expected) == param
        @test InfiniteOpt._param_object_indices(m)[InfiniteOpt._object_number(expected)] == index(expected)
        param = FiniteParameter(1.5)
        expected = GeneralVariableRef(m, 1, FiniteParameterIndex, -1)
        @test add_parameter(m, param) == expected
        @test InfiniteOpt._core_variable_object(expected) == param
    end
end

# Test Reference Queries
@testset "Basic Reference Queries" begin
    m = InfiniteModel()
    p = build_parameter(error, IntervalSet(0, 1))
    pref = add_parameter(m, p)
    dpref = dispatch_variable_ref(pref)
    # JuMP.index
    @testset "JuMP.index" begin
        @test JuMP.index(pref) == IndependentParameterIndex(1)
        @test JuMP.index(dpref) == IndependentParameterIndex(1)
    end
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(pref) === m
        @test owner_model(dpref) === m
    end
end

# Test name methods
@testset "Name" begin
    m = InfiniteModel()
    param = build_parameter(error, IntervalSet(0, 1))
    pref = add_parameter(m, param, "test")
    dpref = dispatch_variable_ref(pref)
    bad = TestVariableRef(m, TestIndex(-1))
    # JuMP.name
    @testset "JuMP.name" begin
        @test_throws ArgumentError name(bad)
        @test name(pref) == "test"
        @test name(dpref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test_throws ArgumentError set_name(bad, "test")
        @test isa(set_name(pref, "new"), Nothing)
        @test name(pref) == "new"
        @test isa(set_name(dpref, "test"), Nothing)
        @test name(dpref) == "test"
    end
    # _make_parameter_ref
    @testset "_make_parameter_ref" begin
        @test InfiniteOpt._make_parameter_ref(m, IndependentParameterIndex(1)) == pref
    end
    # _param_name_dict
    @testset "_param_name_dict" begin
        @test isa(InfiniteOpt._param_name_dict(m), Nothing)
    end
    # _update_param_name_dict
    @testset "_update_param_name_dict" begin
        m.name_to_param = Dict{String, AbstractInfOptIndex}()
        @test InfiniteOpt._update_param_name_dict(m, m.independent_params) isa Nothing
        @test InfiniteOpt._update_param_name_dict(m, m.dependent_params) isa Nothing
        @test m.name_to_param["test"] == IndependentParameterIndex(1)
        m.name_to_param = nothing
    end
    # parameter_by_name
    @testset "parameter_by_name" begin
        @test parameter_by_name(m, "test") == pref
        @test isa(parameter_by_name(m, "test2"), Nothing)
        pref = add_parameter(m, param, "test")
        m.name_to_param = nothing
        @test_throws ErrorException parameter_by_name(m, "test")
    end
end

# Test the parameter macro
@testset "Macro" begin
    m = InfiniteModel()
    # single parameter
    @testset "Single" begin
        pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
        @test @independent_parameter(m, 0 <= a <= 1) == pref
        @test InfiniteOpt._core_variable_object(pref).set == IntervalSet(0, 1)
        @test name(pref) == "a"
        pref = GeneralVariableRef(m, 2, IndependentParameterIndex, -1)
        @test @independent_parameter(m, b in Normal(), supports = [1; 2]) == pref
        @test InfiniteOpt._core_variable_object(pref).set == UniDistributionSet(Normal())
        @test InfiniteOpt._core_variable_object(pref).supports == SortedDict(i => Set{Symbol}([UserDefined]) for i in [1,2])
        pref = GeneralVariableRef(m, 3, IndependentParameterIndex, -1)
        @test @independent_parameter(m, c in IntervalSet(0, 1)) == pref
        @test InfiniteOpt._core_variable_object(pref).set == IntervalSet(0, 1)
        pref = GeneralVariableRef(m, 4, IndependentParameterIndex, -1)
        @test @independent_parameter(m, set = IntervalSet(0, 1),
                                  base_name = "d") == pref
        @test name(pref) == "d"
        pref = GeneralVariableRef(m, 5, IndependentParameterIndex, -1)
        @test @independent_parameter(m, set = IntervalSet(0, 1)) == pref
        @test name(pref) == ""
        pref = GeneralVariableRef(m, 6, IndependentParameterIndex, -1)
        @test @independent_parameter(m, z in [0, 1]) == pref
        @test InfiniteOpt._core_variable_object(pref).set == IntervalSet(0, 1)
        @test name(pref) == "z"
        @test InfiniteOpt._param_object_indices(m)[InfiniteOpt._object_number(pref)] == index(pref)
    end
    # multiple parameters
    @testset "Array" begin
        prefs = [GeneralVariableRef(m, 7, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 8, IndependentParameterIndex, -1)]
        @test @independent_parameter(m, 0 <= e[1:2] <= 1) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(0, 1)
        prefs = [GeneralVariableRef(m, 9, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 10, IndependentParameterIndex, -1)]
        @test @independent_parameter(m, [1:2], set = IntervalSet(0, 1)) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(0, 1)
        prefs = [GeneralVariableRef(m, 11, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 12, IndependentParameterIndex, -1)]
        sets = [IntervalSet(0, 1), IntervalSet(-1, 2)]
        @test @independent_parameter(m, f[i = 1:2], set = sets[i]) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(-1, 2)
        prefs = [GeneralVariableRef(m, 13, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 14, IndependentParameterIndex, -1)]
        @test @independent_parameter(m, [i = 1:2], set = sets[i]) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(-1, 2)
        prefs = [GeneralVariableRef(m, 15, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 16, IndependentParameterIndex, -1)]
        @test @independent_parameter(m, [0, -1][i] <= g[i = 1:2] <= [1, 2][i]) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(-1, 2)
        prefs = [GeneralVariableRef(m, 17, IndependentParameterIndex, -1),
                 GeneralVariableRef(m, 18, IndependentParameterIndex, -1)]
        prefs = convert(JuMP.Containers.SparseAxisArray, prefs)
        @test @independent_parameter(m, 0 <= i[1:2] <= 1,
                                  container = SparseAxisArray) == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).set == IntervalSet(0, 1)
        @test InfiniteOpt._core_variable_object(prefs[2]).set == IntervalSet(0, 1)
        @test InfiniteOpt._param_object_indices(m)[InfiniteOpt._object_number(prefs[2])] == index(prefs[2])
    end
    # test for errors
    @testset "Errors" begin
        @test_macro_throws ErrorException @independent_parameter(m, 0 <= [1:2] <= 1)
        @test_macro_throws ErrorException @independent_parameter(m, 0 <= "bob" <= 1)
        @test_macro_throws ErrorException @independent_parameter(m, 0 <= a <= 1)
        @test_macro_throws ErrorException @independent_parameter(m, 0 <= j)
        @test_macro_throws ErrorException @independent_parameter(m, j)
        @test_macro_throws ErrorException @independent_parameter(m, j, foo = 42)
        @test_macro_throws ErrorException @independent_parameter(m, j in Multinomial(3, [1/3, 1/3]))
        @test_macro_throws ErrorException @independent_parameter(m, 0 <= k <= 1, Int)
    end
end

# Test if used
@testset "Used" begin
    m = InfiniteModel()
    p1 = IndependentParameter(IntervalSet(0, 1), SortedDict{Float64, Set{Symbol}}())
    pref1 = add_parameter(m, p1, "p1")
    dpref1 = dispatch_variable_ref(pref1)
    p2 = FiniteParameter(1)
    pref2 = add_parameter(m, p2, "p2")
    dpref2 = dispatch_variable_ref(pref2)
    # _infinite_variable_dependencies
    @testset "_infinite_variable_dependencies" begin
        @test InfiniteOpt._infinite_variable_dependencies(pref1) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(dpref1) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(pref2) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(dpref2) == InfiniteVariableIndex[]
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(pref1) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(dpref1) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(pref2) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(dpref2) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(pref1) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(dpref1) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(pref2) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(dpref2) == ConstraintIndex[]
    end
    # used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(pref1)
        @test !used_by_constraint(pref2)
        @test !used_by_constraint(dpref1)
        @test !used_by_constraint(dpref2)
        push!(InfiniteOpt._constraint_dependencies(dpref1), ConstraintIndex(1))
        push!(InfiniteOpt._constraint_dependencies(dpref2), ConstraintIndex(1))
        @test used_by_constraint(pref1)
        @test used_by_constraint(pref2)
        @test used_by_constraint(dpref1)
        @test used_by_constraint(dpref2)
        popfirst!(InfiniteOpt._constraint_dependencies(dpref1))
        popfirst!(InfiniteOpt._constraint_dependencies(dpref2))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(pref1)
        @test !used_by_measure(pref2)
        @test !used_by_measure(dpref1)
        @test !used_by_measure(dpref2)
        push!(InfiniteOpt._measure_dependencies(dpref1), MeasureIndex(1))
        push!(InfiniteOpt._measure_dependencies(dpref2), MeasureIndex(1))
        @test used_by_measure(pref1)
        @test used_by_measure(pref2)
        @test used_by_measure(dpref1)
        @test used_by_measure(dpref2)
        popfirst!(InfiniteOpt._measure_dependencies(dpref1))
        popfirst!(InfiniteOpt._measure_dependencies(dpref2))
    end
    # used_by_variable
    @testset "used_by_infinite_variable" begin
        @test !used_by_infinite_variable(pref1)
        @test !used_by_infinite_variable(pref2)
        @test !used_by_infinite_variable(dpref1)
        @test !used_by_infinite_variable(dpref2)
        push!(InfiniteOpt._infinite_variable_dependencies(dpref1), InfiniteVariableIndex(1))
        @test used_by_infinite_variable(pref1)
        @test used_by_infinite_variable(dpref1)
        popfirst!(InfiniteOpt._infinite_variable_dependencies(dpref1))
    end
    # used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(pref1)
        @test !used_by_objective(dpref1)
        @test !used_by_objective(pref2)
        @test !used_by_objective(dpref2)
        InfiniteOpt._data_object(pref2).in_objective = true
        @test used_by_objective(pref2)
        @test used_by_objective(dpref2)
        InfiniteOpt._data_object(pref2).in_objective = false
    end
    # is_used
    @testset "is_used" begin
        @test !is_used(pref1)
        @test !is_used(pref2)
        @test !is_used(dpref1)
        @test !is_used(dpref2)
        push!(InfiniteOpt._constraint_dependencies(dpref1), ConstraintIndex(1))
        push!(InfiniteOpt._constraint_dependencies(dpref2), ConstraintIndex(1))
        @test is_used(pref1)
        @test is_used(pref2)
        @test is_used(dpref1)
        @test is_used(dpref2)
        popfirst!(InfiniteOpt._constraint_dependencies(dpref1))
        popfirst!(InfiniteOpt._constraint_dependencies(dpref2))
        push!(InfiniteOpt._measure_dependencies(dpref1), MeasureIndex(1))
        push!(InfiniteOpt._measure_dependencies(dpref2), MeasureIndex(1))
        @test is_used(pref1)
        @test is_used(pref2)
        @test is_used(dpref1)
        @test is_used(dpref2)
        popfirst!(InfiniteOpt._measure_dependencies(dpref1))
        popfirst!(InfiniteOpt._measure_dependencies(dpref2))
        push!(InfiniteOpt._infinite_variable_dependencies(dpref1), InfiniteVariableIndex(1))
        @test is_used(pref1)
        @test is_used(dpref1)
        popfirst!(InfiniteOpt._infinite_variable_dependencies(dpref1))
    end
end

# Test parameter set methods
@testset "Infinite Set" begin
    m = InfiniteModel()
    param = IndependentParameter(IntervalSet(0, 1), SortedDict{Float64, Set{Symbol}}())
    pref_gen = add_parameter(m, param, "test")
    pref_disp = dispatch_variable_ref(pref_gen)
    bad = Bad()
    # _parameter_set
    @testset "_parameter_set" begin
        @test InfiniteOpt._parameter_set(pref_disp) == IntervalSet(0, 1)
    end
    # _update_parameter_set
    @testset "_update_parameter_set " begin
        @test isa(InfiniteOpt._update_parameter_set(pref_disp,
                                                    IntervalSet(1, 2)), Nothing)
        @test InfiniteOpt._parameter_set(pref_disp) == IntervalSet(1, 2)
    end
    # infinite_set
    @testset "infinite_set" begin
        @test_throws ArgumentError infinite_set(bad)
        @test infinite_set(pref_disp) == IntervalSet(1, 2)
        @test infinite_set(pref_gen) == IntervalSet(1, 2)
    end
    # set_infinite_set
    @testset "set_infinite_set" begin
        @test_throws ArgumentError set_infinite_set(bad, IntervalSet(0, 1))
        @test isa(set_infinite_set(pref_disp, IntervalSet(2, 3)), Nothing)
        @test infinite_set(pref_disp) == IntervalSet(2, 3)
        @test isa(set_infinite_set(pref_gen, IntervalSet(1, 3)), Nothing)
        @test infinite_set(pref_gen) == IntervalSet(1, 3)
        push!(InfiniteOpt._data_object(pref_gen).measure_indices, MeasureIndex(1))
        @test_throws ErrorException set_infinite_set(pref_gen, IntervalSet(1, 3))
        @test_throws ErrorException set_infinite_set(pref_gen, UniDistributionSet(Normal()))
    end
end

# Test parameter support methods
@testset "Supports" begin
    m = InfiniteModel()
    param = IndependentParameter(IntervalSet(0, 1), SortedDict{Float64, Set{Symbol}}())
    pref = add_parameter(m, param, "test")
    pref_disp = dispatch_variable_ref(pref)
    bad = Bad()
    # _parameter_supports
    @testset "_parameter_supports" begin
        @test InfiniteOpt._parameter_supports(pref_disp) == SortedDict{Float64, Set{Symbol}}()
    end
    @testset "_parameter_support_values" begin
        @test InfiniteOpt._parameter_support_values(pref_disp) == Float64[]
    end
    # _update_parameter_supports
    @testset "_update_parameter_supports " begin
        dict = SortedDict{Float64, Set{Symbol}}(1. => Set{Symbol}())
        @test isa(InfiniteOpt._update_parameter_supports(pref_disp, dict), Nothing)
        @test InfiniteOpt._parameter_support_values(pref_disp) == [1.]
    end
    # num_supports
    @testset "num_supports" begin
        @test_throws ArgumentError num_supports(bad)
        @test num_supports(pref_disp) == 1
        @test num_supports(pref_disp, label = UserDefined) == 0
        @test num_supports(pref) == 1
        @test num_supports(pref, label = UserDefined) == 0
    end
    # has_supports
    @testset "has_supports" begin
        @test_throws ArgumentError has_supports(bad)
        @test has_supports(pref_disp)
        @test has_supports(pref)
        InfiniteOpt._update_parameter_supports(pref_disp, SortedDict{Float64, Set{Symbol}}())
        @test !has_supports(pref_disp)
        @test !has_supports(pref)
    end
    # supports
    @testset "supports" begin
        @test_throws ArgumentError supports(bad)
        @test_throws ErrorException supports(pref_disp)
        @test_throws ErrorException supports(pref)
        dict = SortedDict{Float64, Set{Symbol}}(1. => Set{Symbol}())
        InfiniteOpt._update_parameter_supports(pref_disp, dict)
        @test supports(pref_disp) == [1.]
        @test supports(pref) == [1.]
    end
    # set_supports
    @testset "set_supports" begin
        @test_throws ArgumentError set_supports(bad, [0, 1])
        @test isa(set_supports(pref_disp, [0, 1], force = true), Nothing)
        @test isa(set_supports(pref, [0, 1], force = true), Nothing)
        @test supports(pref) == [0., 1.]
        @test_throws ErrorException set_supports(pref, [2, 3])
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) set_supports(pref, [1, 1], force = true)
        @test_throws ErrorException set_supports(pref, [0.5])
    end
    # add_supports
    @testset "add_supports" begin
        @test_throws ArgumentError add_supports(bad, 0.5)
        @test isa(add_supports(pref_disp, 0.25), Nothing)
        @test isa(add_supports(pref, 0.5), Nothing)
        @test supports(pref) == [0.25, 0.5, 1.]
        @test isa(add_supports(pref, [0, 0.25, 1]), Nothing)
        @test supports(pref) == [0, 0.25, 0.5, 1.]
    end
    # delete_supports
    @testset "delete_supports" begin
        @test_throws ArgumentError delete_supports(bad)
        @test isa(delete_supports(pref_disp), Nothing)
        add_supports(pref, 0.1)
        @test isa(delete_supports(pref), Nothing)
        @test_throws ErrorException supports(pref)
        push!(InfiniteOpt._data_object(pref).measure_indices, MeasureIndex(1))
        @test_throws ErrorException delete_supports(pref)
    end
end

# Test lower bound functions
@testset "Lower Bound" begin
    m = InfiniteModel()
    p1 = IndependentParameter(IntervalSet(0, 1), SortedDict{Float64, Set{Symbol}}())
    pref1 = add_parameter(m, p1)
    p2 = IndependentParameter(UniDistributionSet(Normal()), SortedDict{Float64, Set{Symbol}}())
    pref2 = add_parameter(m, p2)
    p3 = IndependentParameter(BadScalarSet(), SortedDict{Float64, Set{Symbol}}())
    pref3 = add_parameter(m, p3)
    bad = TestVariableRef(m, TestIndex(-1))
    # JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test_throws ArgumentError has_lower_bound(bad)
        @test has_lower_bound(dispatch_variable_ref(pref1))
        @test has_lower_bound(pref1)
        @test has_lower_bound(dispatch_variable_ref(pref2))
        @test has_lower_bound(pref2)
    end
    # JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        @test_throws ArgumentError lower_bound(bad)
        @test lower_bound(dispatch_variable_ref(pref1)) == 0
        @test lower_bound(dispatch_variable_ref(pref2)) == -Inf
        @test lower_bound(pref1) == 0
        @test lower_bound(pref2) == -Inf
        @test_throws ErrorException lower_bound(pref3)
    end
    # JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test_throws ArgumentError set_lower_bound(bad, 2)
        @test_throws ErrorException set_lower_bound(dispatch_variable_ref(pref1), 2)
        @test_throws ErrorException set_lower_bound(dispatch_variable_ref(pref3), 2)
        @test_throws ErrorException set_lower_bound(pref1, 2)
        @test_throws ErrorException set_lower_bound(pref3, 2)
        @test isa(set_lower_bound(dispatch_variable_ref(pref1), -1), Nothing)
        @test lower_bound(pref1) == -1
        @test isa(set_lower_bound(pref1, -2), Nothing)
        @test lower_bound(pref1) == -2
    end
end

# Test upper bound functions
@testset "Upper Bound" begin
    m = InfiniteModel()
    p1 = IndependentParameter(IntervalSet(0, 1), SortedDict{Float64, Set{Symbol}}())
    pref1 = add_parameter(m, p1)
    p2 = IndependentParameter(UniDistributionSet(Normal()), SortedDict{Float64, Set{Symbol}}())
    pref2 = add_parameter(m, p2)
    p3 = IndependentParameter(BadScalarSet(), SortedDict{Float64, Set{Symbol}}())
    pref3 = add_parameter(m, p3)
    bad = TestVariableRef(m, TestIndex(-1))
    # JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test_throws ArgumentError has_upper_bound(bad)
        @test has_upper_bound(dispatch_variable_ref(pref1))
        @test has_upper_bound(dispatch_variable_ref(pref2))
        @test has_upper_bound(pref1)
        @test has_upper_bound(pref2)
    end
    # JuMP.lower_bound
    @testset "JuMP.upper_bound" begin
        @test_throws ArgumentError upper_bound(bad)
        @test upper_bound(dispatch_variable_ref(pref1)) == 1
        @test upper_bound(dispatch_variable_ref(pref2)) == Inf
        @test upper_bound(pref1) == 1
        @test upper_bound(pref2) == Inf
        @test_throws ErrorException upper_bound(pref3)
    end
    # JuMP.set_lower_bound
    @testset "JuMP.set_upper_bound" begin
        @test_throws ArgumentError set_upper_bound(bad, -1)
        @test_throws ErrorException set_upper_bound(dispatch_variable_ref(pref1), -1)
        @test_throws ErrorException set_upper_bound(dispatch_variable_ref(pref3), -1)
        @test_throws ErrorException set_upper_bound(pref1, -1)
        @test_throws ErrorException set_upper_bound(pref3, -1)
        @test isa(set_upper_bound(dispatch_variable_ref(pref1), 2), Nothing)
        @test upper_bound(pref1) == 2
        @test isa(set_upper_bound(pref1, 3), Nothing)
        @test upper_bound(pref1) == 3
    end
end

# Test methods for finite parameters
@testset "Finite Parameters" begin
    # initialize the model
    m = InfiniteModel()
    bad = TestVariableRef(m, TestIndex(-1))
    # test @finite_parameter
    @testset "@finite_parameter" begin
        m2 = Model()
        # test errors
        @test_macro_throws ErrorException @finite_parameter(m)
        @test_macro_throws ErrorException @finite_parameter(m, a, 2, 3)
        @test_macro_throws ErrorException @finite_parameter(m, (2, 3, 4), 2)
        @test_macro_throws ErrorException @finite_parameter(m2, 2)
        @test_macro_throws ErrorException @finite_parameter(m, "bob")
        @test_macro_throws ErrorException @finite_parameter(m2, test, 2)
        @test_macro_throws ErrorException @finite_parameter(m, test, 2, bob = 2)
        # test anonymous definition
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        @test @finite_parameter(m, 42) == pref
        @test InfiniteOpt._core_variable_object(pref).value == 42
        # test vector anonymous definition
        prefs = [GeneralVariableRef(m, 2, FiniteParameterIndex),
                 GeneralVariableRef(m, 3, FiniteParameterIndex)]
        @test @finite_parameter(m, [1:2], 42, base_name = "a") == prefs
        @test InfiniteOpt._core_variable_object(prefs[1]).value == 42
        @test name.(prefs) == ["a[1]", "a[2]"]
        # test named definition
        pref = GeneralVariableRef(m, 4, FiniteParameterIndex)
        @test @finite_parameter(m, b, 42) == pref
        @test InfiniteOpt._core_variable_object(pref).value == 42
        @test name(pref) == "b"
        # test named vector definition
        prefs = [GeneralVariableRef(m, 5, FiniteParameterIndex),
                 GeneralVariableRef(m, 6, FiniteParameterIndex)]
        prefs = convert(JuMPC.SparseAxisArray, prefs)
        @test @finite_parameter(m, c[i = 1:2], [3, 7][i],
                                container = SparseAxisArray) == prefs
        @test InfiniteOpt._core_variable_object(prefs[2]).value == 7
        @test name(prefs[2]) == "c[2]"
    end
    # initialize the model
    m = InfiniteModel()
    # test JuMP.value
    @testset "JuMP.value" begin
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        dpref = dispatch_variable_ref(pref)
        @test_throws ArgumentError value(bad)
        @test @finite_parameter(m, g, 1) == pref
        @test value(dpref) == 1
        @test value(pref) == 1
    end
    # test JuMP.set_value
    @testset "JuMP.set_value" begin
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        dpref = dispatch_variable_ref(pref)
        @test_throws ArgumentError set_value(bad, 42)
        @test isa(set_value(dpref, 42), Nothing)
        @test value(pref) == 42
        @test isa(set_value(pref, 41), Nothing)
        @test value(pref) == 41
    end
end

# Test support flll-in and geneartion functions
@testset "Support Fill-in and Generation" begin
    @testset "generate_and_add_supports! (AbstractInfiniteSet)" begin
        m = InfiniteModel()
        gvref1 = @independent_parameter(m, 0 <= a <= 1)
        pref1 = dispatch_variable_ref(gvref1)
        set1 = infinite_set(pref1)
        dist = Normal(0., 1.)
        gvref2 = @independent_parameter(m, c in dist)
        pref2 = dispatch_variable_ref(gvref2)
        set2 = infinite_set(pref2)
        @test generate_and_add_supports!(pref1, set1, num_supports = 10) isa Nothing
        @test generate_and_add_supports!(pref2, set2, num_supports = 10) isa Nothing
        @test length(supports(pref1)) == 10
        @test length(supports(pref2)) == 10
    end
    # fill_in_supports! (ParameterRef)
    @testset "fill_in_supports! (ParameterRef)" begin
        m = InfiniteModel()
        pref1 = @independent_parameter(m, 0 <= a <= 1)
        pref2 = @independent_parameter(m, 0 <= b[1:2] <= 1)
        dist = Normal(0., 1.)
        pref3 = @independent_parameter(m, c in dist, supports = [-0.5, 0.5])
        @test fill_in_supports!(pref1, num_supports = 10, sig_figs = 4) isa Nothing
        @test fill_in_supports!.(pref2, num_supports = 10, sig_figs = 4) isa Array{Nothing}
        @test fill_in_supports!(pref3, num_supports = 10, sig_figs = 4) isa Nothing
        @test length(supports(pref1)) == 10
        @test length(supports(pref2[1])) == 10
        @test length(supports(pref2[2])) == 10
        @test length(supports(pref3)) == 10
        @test -0.5 in supports(pref3)
        @test 0.5 in supports(pref3)
        @test fill_in_supports!(pref1, num_supports = 20) isa Nothing
        @test length(supports(pref1)) == 20
    end
end
