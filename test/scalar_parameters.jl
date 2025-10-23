# Test Core data accessors
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    ind_idx = IndependentParameterIndex(1)
    fin_idx = FiniteParameterIndex(1)
    domain = IntervalDomain(0, 1)
    supps_dict = SortedDict{Float64, Set{DataType}}(0. => Set{DataType}([UserDefined]))
    method = InfiniteOpt.DefaultDerivativeMethod
    info = NoGenerativeSupports()
    ind_param = IndependentParameter(domain, supps_dict, 5, method, info)
    fin_param = FiniteParameter(42)
    ind_object = ScalarParameterData(ind_param, 1, 1, "ind")
    fin_object = ScalarParameterData(fin_param, -1, -1, "fin")
    ind_pref = IndependentParameterRef(m, ind_idx)
    fin_pref = FiniteParameterRef(m, fin_idx)
    ind_gvref = GeneralVariableRef(m, 1, IndependentParameterIndex)
    fin_gvref = GeneralVariableRef(m, 1, FiniteParameterIndex)
    bad_ind_pref = IndependentParameterRef(m, IndependentParameterIndex(-1))
    bad_fin_pref = FiniteParameterRef(m, FiniteParameterIndex(-1))
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
        @test InfiniteOpt.parameter_group_indices(m)[end] == ind_idx
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
        @test_throws ErrorException InfiniteOpt._data_object(bad_ind_pref)
        @test_throws ErrorException InfiniteOpt._data_object(bad_fin_pref)
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test core_object(ind_pref) == ind_param
        @test core_object(ind_gvref) == ind_param
        @test core_object(fin_pref) == fin_param
        @test core_object(fin_gvref) == fin_param
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
    # test parameter_group_int_index
    @testset "parameter_group_int_index" begin
        @test InfiniteOpt.parameter_group_int_index(ind_pref) == 1
        @test InfiniteOpt.parameter_group_int_index(ind_gvref) == 1
    end
    # test parameter_group_int_indices
    @testset "parameter_group_int_indices" begin
        @test InfiniteOpt.parameter_group_int_indices(ind_pref) == [1]
        @test InfiniteOpt.parameter_group_int_indices(ind_gvref) == [1]
    end
    # test _adaptive_data_update
    @testset "_adaptive_data_update" begin
        # test change of same type
        p = IndependentParameter(IntervalDomain(0, 2), supps_dict, 12, method, info)
        data = InfiniteOpt._data_object(ind_pref)
        @test InfiniteOpt._adaptive_data_update(ind_pref, p, data) isa Nothing
        @test core_object(ind_pref) == p 
        # test change of different types 
        p = IndependentParameter(IntervalDomain(0, 2), supps_dict, 12, TestMethod(), info)
        data = InfiniteOpt._data_object(ind_pref)
        @test InfiniteOpt._adaptive_data_update(ind_pref, p, data) isa Nothing
        @test core_object(ind_pref) == p 
    end
    # test _set_core_object
    @testset "_set_core_object" begin
        @test InfiniteOpt._set_core_object(ind_pref, ind_param) isa Nothing
        @test InfiniteOpt._set_core_object(fin_pref, fin_param) isa Nothing
    end
    # test _delete_data_object
    @testset "_delete_data_object" begin
        @test is_valid(m, ind_pref)
        @test is_valid(m, fin_pref)
        @test InfiniteOpt._delete_data_object(ind_pref) isa Nothing
        @test InfiniteOpt._delete_data_object(fin_pref) isa Nothing
        @test !is_valid(m, ind_pref)
        @test !is_valid(m, fin_pref)
    end
end

# Test macro methods
@testset "Macro Helpers" begin
    # test _distribution_or_error
    @testset "_distribution_or_error" begin
        # univariate
        dom = UniDistributionDomain(Normal())
        @test InfiniteOpt._distribution_or_error(error, Normal()) == dom 
        # multivariate
        dist = MvNormal(ones(2), [1 0; 0 1])
        dom = MultiDistributionDomain(dist)
        @test InfiniteOpt._distribution_or_error(error, dist) == dom 
        # test not distribution 
        @test_throws ErrorException InfiniteOpt._distribution_or_error(error, 2)
    end
    # test _domain_or_error
    @testset "_domain_or_error" begin
        # domains
        dom = IntervalDomain(0, 1)
        @test InfiniteOpt._domain_or_error(error, dom) == dom 
        # vector 
        @test InfiniteOpt._domain_or_error(error, [0, 1]) == IntervalDomain(0, 1)
        @test_throws ErrorException InfiniteOpt._domain_or_error(error, [0, 1, 2])
        # distribution 
        @test_throws ErrorException InfiniteOpt._domain_or_error(error, Normal())
        # something else 
        @test_throws ErrorException InfiniteOpt._domain_or_error(error, 2)
    end
    # test _make_distribution_call
    @testset "_make_distribution_call" begin
        expected = :( _distribution_or_error($error, $(esc(:d))) )
        @test InfiniteOpt._make_distribution_call(error, :d) == expected
    end
    # test _make_domain_call
    @testset "_make_domain_call" begin
        expected = :( _domain_or_error($error, $(esc(:d))) )
        @test InfiniteOpt._make_domain_call(error, :d) == expected
    end
    # test _parse_parameter
    @testset "_parse_parameter" begin
        # domain
        expr = :( _domain_or_error($error, $(esc(:d))) )
        @test InfiniteOpt._parse_parameter(error, Val(:in), :p, :d) == (:p, expr)
        # distribution
        expr = :( _distribution_or_error($error, $(esc(:d))) )
        @test InfiniteOpt._parse_parameter(error, Val(:~), :p, :d) == (:p, expr)
        # fallback 
        @test_throws ErrorException InfiniteOpt._parse_parameter(error, Val(:(==)), :p, :d)
    end
end

# Test parameter definition methods
@testset "Definition" begin
    # _check_supports_in_bounds
    @testset "_process_scalar_supports" begin
        domain = IntervalDomain(0, 1)
        @test InfiniteOpt._process_scalar_supports(error, 0, domain, 8) == 0
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, -1, domain, 8)
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, 2, domain, 8)
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, "2", domain, 8)
        @test InfiniteOpt._process_scalar_supports(error, [0.5, 1], domain, 8) == [0.5, 1]
        @test InfiniteOpt._process_scalar_supports(error, 0:1, domain, 8) == [0, 1]
        @test InfiniteOpt._process_scalar_supports(error, 0:0.5:1, domain, 8) == [0, 0.5, 1]
        @test InfiniteOpt._process_scalar_supports(error, (0, 0.5), domain, 8) == [0.0, 0.5]
        @test InfiniteOpt._process_scalar_supports(error, (i for i in (0, 1)), domain, 8) == [0, 1]
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, (0, 2), domain, 8)
        domain = UniDistributionDomain(Uniform())
        @test InfiniteOpt._process_scalar_supports(error, 0, domain, 8) == 0
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, -1, domain, 8)
        @test_throws ErrorException InfiniteOpt._process_scalar_supports(error, 2, domain, 8)
    end
    # build_independent_parameter
    @testset "build_parameter (IndependentParameter)" begin
        domain = IntervalDomain(0, 1)
        supps = 0.
        supps_dict = SortedDict{Float64, Set{DataType}}(0. => Set([UserDefined]))
        method = TestMethod()
        geninfo = NoGenerativeSupports()
        @test build_parameter(error, domain, supports = supps).domain == domain
        @test build_parameter(error, domain, supports = supps).supports == supps_dict
        @test_throws ErrorException build_parameter(error, domain, bob = 42)
        warn = "Ignoring `num_supports` since `supports` is not empty."
        @test_logs (:warn, warn) build_parameter(error, domain,
                                            supports = [0, 1], num_supports = 2)
        repeated_supps = [1, 1]
        expected = IndependentParameter(domain, SortedDict{Float64, Set{DataType}}(1. => Set{DataType}()), 5, method, geninfo)
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) build_parameter(error, domain, supports = repeated_supps, 
                                                 derivative_method = method) == expected
        domain = UniDistributionDomain(Normal())
        @test length(build_parameter(error, domain, num_supports = 5).supports) == 5
        @test build_parameter(error, domain, derivative_method = method).derivative_method == method
        @test collect(keys(build_parameter(error, domain, supports = -1:0.5:1).supports)) == [-1, -0.5, 0, 0.5, 1]
    end
    # build_finite_parameter
    @testset "build_parameter (FiniteParameter)" begin
        @test_throws ErrorException build_parameter(error, 1, bob = 42)
        expected = FiniteParameter(1)
        @test build_parameter(error, 1) == expected
    end
    # build_parameter fallbacks
    @testset "build_parameter (Fallbacks)" begin
        @test_throws ErrorException build_parameter(error, CollectionDomain([IntervalDomain(0, 1)]))
        @test_throws ErrorException build_parameter(error, :bad)
    end
    # add_parameter
    @testset "add_parameter" begin
        m = InfiniteModel()
        method = InfiniteOpt.DefaultDerivativeMethod
        geninfo = NoGenerativeSupports()
        param = IndependentParameter(IntervalDomain(0, 1),
                                    SortedDict{Float64, Set{DataType}}(), 5, 
                                    method, geninfo)
        expected = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
        @test isequal(add_parameter(m, param), expected)
        @test core_object(expected) == param
        @test InfiniteOpt.parameter_group_indices(m)[InfiniteOpt.parameter_group_int_index(expected)] == index(expected)
        param = FiniteParameter(1.5)
        expected = GeneralVariableRef(m, 1, FiniteParameterIndex, -1)
        @test isequal(add_parameter(m, param), expected)
        @test core_object(expected) == param
    end
end

# Test Reference Queries
@testset "Basic Reference Queries" begin
    m = InfiniteModel()
    p = build_parameter(error, IntervalDomain(0, 1), derivative_method = TestMethod())
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
    # has_derivative_supports
    @testset "has_generative_supports" begin
        @test !has_generative_supports(pref)
        InfiniteOpt._data_object(pref).has_generative_supports = true
        @test has_generative_supports(dpref)
    end
    # _set_has_generative_supports
    @testset "_set_has_generative_supports" begin
        @test InfiniteOpt._set_has_generative_supports(pref, true) isa Nothing
        @test has_generative_supports(pref)
        @test InfiniteOpt._set_has_generative_supports(dpref, false) isa Nothing
        @test !has_generative_supports(pref)
    end
    # test generative_support_info
    @testset "generative_support_info" begin 
        @test generative_support_info(pref) == NoGenerativeSupports()
    end
    # derivative_method
    @testset "derivative_method" begin
        @test derivative_method(pref) isa TestMethod
        @test derivative_method(dpref) isa TestMethod
    end
    # has_internal_supports
    @testset "has_internal_supports" begin
        @test !has_internal_supports(pref)
        InfiniteOpt._data_object(pref).has_internal_supports = true
        @test has_internal_supports(dpref)
    end
    # _set_has_internal_supports
    @testset "_set_has_internal_supports" begin
        @test InfiniteOpt._set_has_internal_supports(pref, true) isa Nothing
        @test has_internal_supports(pref)
        @test InfiniteOpt._set_has_internal_supports(dpref, false) isa Nothing
        @test !has_internal_supports(pref)
    end
    # has_derivative_constraints
    @testset "has_derivative_constraints" begin
        @test !has_derivative_constraints(pref)
        InfiniteOpt._data_object(pref).has_deriv_constrs = true
        @test has_derivative_constraints(dpref)
    end
    # _set_has_derivative_constraints
    @testset "_set_has_derivative_constraints" begin
        @test InfiniteOpt._set_has_derivative_constraints(pref, true) isa Nothing
        @test has_derivative_constraints(pref)
        @test InfiniteOpt._set_has_derivative_constraints(dpref, false) isa Nothing
        @test !has_derivative_constraints(pref)
    end
end

# Test name methods
@testset "Name" begin
    m = InfiniteModel()
    param = build_parameter(error, IntervalDomain(0, 1))
    pref = add_parameter(m, param, "test")
    dpref = dispatch_variable_ref(pref)
    bad_pref = FiniteParameterRef(m, FiniteParameterIndex(-1))
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(pref) == "test"
        @test name(dpref) == "test"
        @test name(bad_pref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(pref, "new"), Nothing)
        @test name(pref) == "new"
        @test isa(set_name(dpref, "test"), Nothing)
        @test name(dpref) == "test"
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
        @test isequal(parameter_by_name(m, "test"), pref)
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
        pref = GeneralVariableRef(m, 1, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m, a in [0, 1]), pref)
        @test core_object(pref).domain == IntervalDomain(0, 1)
        @test name(pref) == "a"

        pref = GeneralVariableRef(m, 2, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m, b ~ Normal(), supports = [1; 2]), pref)
        @test core_object(pref).domain == UniDistributionDomain(Normal())
        @test core_object(pref).supports == SortedDict(i => Set{DataType}([UserDefined]) for i in [1,2])

        pref = GeneralVariableRef(m, 3, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m, c ∈ IntervalDomain(0, 1)), pref)
        @test core_object(pref).domain == IntervalDomain(0, 1)

        pref = GeneralVariableRef(m, 4, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m; domain = IntervalDomain(0, 1),
                      base_name = "d"), pref)
        @test name(pref) == "d"

        pref = GeneralVariableRef(m, 5, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m, distribution = Normal()), pref)
        @test name(pref) == ""
        @test core_object(pref).domain == UniDistributionDomain(Normal())

        pref = GeneralVariableRef(m, 6, IndependentParameterIndex)
        @test isequal(@infinite_parameter(m, z in [0, 1], derivative_method = TestMethod()), pref)
        @test core_object(pref).domain == IntervalDomain(0, 1)
        @test name(pref) == "z"
        @test derivative_method(pref) isa TestMethod
        @test InfiniteOpt.parameter_group_indices(m)[InfiniteOpt.parameter_group_int_index(pref)] == index(pref)
    end
    # multiple parameters
    @testset "Array" begin
        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 7:8]
        @test isequal(@infinite_parameter(m, e[1:2] in [0, 1], independent = true), prefs)
        @test core_object(prefs[1]).domain == IntervalDomain(0, 1)
        @test core_object(prefs[2]).domain == IntervalDomain(0, 1)

        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 9:10]
        @test isequal(@infinite_parameter(m, [1:2]; domain = IntervalDomain(0, 1), 
                                  independent = true), prefs)
        @test core_object(prefs[1]).domain == IntervalDomain(0, 1)
        @test core_object(prefs[2]).domain == IntervalDomain(0, 1)

        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 11:12]
        domains = [IntervalDomain(0, 1), IntervalDomain(-1, 2)]
        @test isequal(@infinite_parameter(m, f[i = 1:2], domain = domains[i], 
                                  independent = true), prefs)
        @test core_object(prefs[1]).domain == IntervalDomain(0, 1)
        @test core_object(prefs[2]).domain == IntervalDomain(-1, 2)

        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 13:14]
        @test isequal(@infinite_parameter(m, [i = 1:2], domain = domains[i], 
                                  independent = true), prefs)
        @test core_object(prefs[1]).domain == IntervalDomain(0, 1)
        @test core_object(prefs[2]).domain == IntervalDomain(-1, 2)

        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 15:16]
        @test isequal(@infinite_parameter(m, g[i = 1:2] in [1 - i, i], 
                                  independent = true), prefs)
        @test core_object(prefs[1]).domain == IntervalDomain(0, 1)
        @test core_object(prefs[2]).domain == IntervalDomain(-1, 2)

        prefs = [GeneralVariableRef(m, i, IndependentParameterIndex) for i in 17:18]
        @test all(isequal.([@infinite_parameter(m, i[1:2] ~ Normal(); independent = true,
                                  container = SparseAxisArray)...], prefs))
        @test core_object(prefs[1]).domain == UniDistributionDomain(Normal())
        @test core_object(prefs[2]).domain == UniDistributionDomain(Normal())
        @test InfiniteOpt.parameter_group_indices(m)[InfiniteOpt.parameter_group_int_index(prefs[2])] == index(prefs[2])
    end
    # test for errors
    @testset "Errors" begin
        @test_macro_throws ErrorException @infinite_parameter(m, z in [0, 1], 3)
        @test_macro_throws ErrorException @infinite_parameter(m, z in [0, 1], 
                                                   domain = IntervalDomain(0, 1))
        @test_macro_throws ErrorException @infinite_parameter(m, z in [0, 1], 
                                                         distribution = Normal())
        @test_macro_throws ErrorException @infinite_parameter(m, 
                                                         distribution = Normal(), 
                                                   domain = IntervalDomain(0, 1))
        @test_macro_throws ErrorException @infinite_parameter(m)
        @test_macro_throws ErrorException @infinite_parameter()
        @test_macro_throws ErrorException @infinite_parameter(m, 0 <= z <= 1)
        @test_macro_throws ErrorException @infinite_parameter(m, [1:2] in [0, 1], 
                                                              independent = a)
        @test_macro_throws ErrorException @infinite_parameter(m, "bob" in [0, 1])
        @test_macro_throws ErrorException @infinite_parameter(m, [m = 1:2] in [0, 1])
        @test_macro_throws ErrorException @infinite_parameter(Model(), z in [0, 1])
        @test_macro_throws ErrorException @infinite_parameter(m, z <= 0)
        @test_macro_throws ErrorException @infinite_parameter(m, z in 23)
        @test_macro_throws ErrorException @infinite_parameter(m, z ~ 42)
        @test_macro_throws ErrorException @infinite_parameter(m, z[a...] in [0, 1])
        @test_macro_throws ErrorException @infinite_parameter(m, z[i = s, i in d] in [0, 1])
        @test_macro_throws ErrorException @infinite_parameter(m, z in [0, 1], bad = 2)
        @test_macro_throws ErrorException @infinite_parameter(m, z ~ MvNormal([0, 0], [1 0; 0 1]))
    end
end

# Test if used
@testset "Used" begin
    m = InfiniteModel()
    @infinite_parameter(m, pref1 in [0, 1])
    @finite_parameter(m, pref2 == 1)
    dpref1 = dispatch_variable_ref(pref1)
    dpref2 = dispatch_variable_ref(pref2)
    bad_pref = IndependentParameterRef(m, IndependentParameterIndex(-1))
    # _infinite_variable_dependencies
    @testset "_infinite_variable_dependencies" begin
        @test InfiniteOpt._infinite_variable_dependencies(pref1) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(dpref1) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(pref2) == InfiniteVariableIndex[]
        @test InfiniteOpt._infinite_variable_dependencies(dpref2) == InfiniteVariableIndex[]
        @test_throws ErrorException InfiniteOpt._infinite_variable_dependencies(bad_pref)
    end
    # _parameter_function_dependencies
    @testset "_parameter_function_dependencies" begin
        @test InfiniteOpt._parameter_function_dependencies(pref1) == ParameterFunctionIndex[]
        @test InfiniteOpt._parameter_function_dependencies(dpref1) == ParameterFunctionIndex[]
        @test InfiniteOpt._parameter_function_dependencies(pref2) == ParameterFunctionIndex[]
        @test InfiniteOpt._parameter_function_dependencies(dpref2) == ParameterFunctionIndex[]
        @test_throws ErrorException InfiniteOpt._parameter_function_dependencies(bad_pref)
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(pref1) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(dpref1) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(pref2) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(dpref2) == DerivativeIndex[]
        @test_throws ErrorException InfiniteOpt._derivative_dependencies(bad_pref)
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(pref1) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(dpref1) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(pref2) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(dpref2) == MeasureIndex[]
        @test_throws ErrorException InfiniteOpt._measure_dependencies(bad_pref)
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(pref1) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(dpref1) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(pref2) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(dpref2) == InfOptConstraintIndex[]
        @test_throws ErrorException InfiniteOpt._constraint_dependencies(bad_pref)
    end
    # _generative_measures
    @testset "_generative_measures" begin
        @test InfiniteOpt._generative_measures(pref1) == MeasureIndex[]
        @test InfiniteOpt._generative_measures(dpref1) == MeasureIndex[]
        @test InfiniteOpt._generative_measures(pref2) == MeasureIndex[]
        @test InfiniteOpt._generative_measures(dpref2) == MeasureIndex[]
        @test_throws ErrorException InfiniteOpt._generative_measures(bad_pref)
    end
    # used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(pref1)
        @test !used_by_constraint(pref2)
        @test !used_by_constraint(dpref1)
        @test !used_by_constraint(dpref2)
        push!(InfiniteOpt._constraint_dependencies(dpref1), InfOptConstraintIndex(1))
        push!(InfiniteOpt._constraint_dependencies(dpref2), InfOptConstraintIndex(1))
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
    # used_by_infinite_variable
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
    # used_by_parameter_function
    @testset "used_by_parameter_function" begin
        @test !used_by_parameter_function(pref1)
        @test !used_by_parameter_function(pref2)
        @test !used_by_parameter_function(dpref1)
        @test !used_by_parameter_function(dpref2)
        push!(InfiniteOpt._parameter_function_dependencies(dpref1), ParameterFunctionIndex(1))
        @test used_by_parameter_function(pref1)
        @test used_by_parameter_function(dpref1)
        popfirst!(InfiniteOpt._parameter_function_dependencies(dpref1))
    end
    # used_by_derivative
    @testset "used_by_derivative" begin
        @test !used_by_derivative(pref1)
        @test !used_by_derivative(pref2)
        @test !used_by_derivative(dpref1)
        @test !used_by_derivative(dpref2)
        push!(InfiniteOpt._derivative_dependencies(dpref1), DerivativeIndex(1))
        @test used_by_derivative(pref1)
        @test used_by_derivative(dpref1)
        popfirst!(InfiniteOpt._derivative_dependencies(dpref1))
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
        push!(InfiniteOpt._constraint_dependencies(dpref1), InfOptConstraintIndex(1))
        push!(InfiniteOpt._constraint_dependencies(dpref2), InfOptConstraintIndex(1))
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
        push!(InfiniteOpt._parameter_function_dependencies(dpref1), ParameterFunctionIndex(1))
        @test is_used(pref1)
        @test is_used(dpref1)
        popfirst!(InfiniteOpt._parameter_function_dependencies(dpref1))
    end
end

# Test generative support info methods 
@testset "Generative Support Info" begin 
    # setup info 
    m = InfiniteModel()
    @infinite_parameter(m, pref in [0, 5])
    dpref = dispatch_variable_ref(pref)
    info = UniformGenerativeInfo([0.5], InternalLabel)
    method = FiniteDifference()
    supps = SortedDict{Float64, Set{DataType}}(0 => Set([All]), 
                                               2.5 => Set([InternalLabel]), 
                                               5 => Set([All]))
    new_param = IndependentParameter(IntervalDomain(0, 5), supps, 6, method, info)
    InfiniteOpt._set_core_object(dpref, new_param)
    push!(InfiniteOpt._measure_dependencies(dpref), MeasureIndex(-1))
    # test Base.copy 
    @testset "Base.copy" begin 
        @test copy(info) == info 
        @test copy(info) !== info 
        @test copy(NoGenerativeSupports()) == NoGenerativeSupports()
    end
    # test support label 
    @testset "support_label" begin 
        @test_throws ErrorException support_label(TestGenInfo())
        @test support_label(info) == InternalLabel
        @test support_label(NoGenerativeSupports()) == InfiniteOpt._NoLabel
    end
    # test _reset_generative_supports
    @testset "_reset_generative_supports" begin 
        @test InfiniteOpt._set_has_generative_supports(dpref, true) isa Nothing
        @test has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 2.5, 5]
        @test InfiniteOpt._reset_generative_supports(dpref) isa Nothing
        @test !has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 5]
    end
    # test _set_generative_support_info
    @testset "_set_generative_support_info" begin 
        # test 1
        @test InfiniteOpt._set_generative_support_info(dpref, NoGenerativeSupports()) isa Nothing
        @test generative_support_info(pref) == NoGenerativeSupports()
        # test 2
        @test InfiniteOpt._set_generative_support_info(dpref, info) isa Nothing
        @test generative_support_info(pref) == info
    end
    # test make_generative_supports
    @testset "make_generative_supports" begin 
        # test errors
        @test_throws ErrorException make_generative_supports(TestGenInfo(), dpref, [0, 1])
        @test_throws ErrorException make_generative_supports(info, dpref, [0.0])
        # test normal
        @test make_generative_supports(info, dpref, [0.0, 2.0, 5.0]) == [1, 3.5]
        @test make_generative_supports(info, dpref, [0.0, 5.0]) == [2.5]
    end
    # test _add_generative_supports 
    @testset "_add_generative_supports " begin 
        # test with NoGenerativeSupports
        @test InfiniteOpt._add_generative_supports(dpref, NoGenerativeSupports()) isa Nothing 
        @test !has_generative_supports(dpref)
        # test UniformGenerativeInfo
        @test InfiniteOpt._add_generative_supports(dpref, info) isa Nothing 
        @test has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 2.5, 5]
        # test subsequent run
        @test InfiniteOpt._add_generative_supports(dpref, info) isa Nothing 
        @test has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 2.5, 5]
        # undo changes 
        @test delete_supports(dpref, label = InternalLabel) isa Nothing 
        @test !has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 5]
    end
    # test add_generative_supports
    @testset "add_generative_supports" begin 
        # test normal 
        @test add_generative_supports(dpref) isa Nothing 
        @test has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 2.5, 5]
        # test subsequent run
        @test add_generative_supports(dpref) isa Nothing 
        @test has_generative_supports(dpref)
        @test supports(dpref, label = All) == [0, 2.5, 5]
    end
end

# Test derivative methods 
@testset "Derivative Methods" begin 
    m = InfiniteModel()
    @infinite_parameter(m, pref in [0, 1])
    dpref = dispatch_variable_ref(pref)
    func = (x) -> NaN
    num = 0.
    info = VariableInfo(true, num, true, num, true, num, false, func, false, false)
    d = Derivative(info, true, pref, pref, 1) # this is wrong but that is ok
    object = VariableData(d)
    idx = InfiniteOpt._add_data_object(m, object)
    push!(InfiniteOpt._derivative_dependencies(pref), idx)
    dref = DerivativeRef(m, idx)
    gdref = GeneralVariableRef(m, idx.value, DerivativeIndex)
    cref = @constraint(m, gdref == 0)
    info = UniformGenerativeInfo([0.5], InternalGaussLobatto)
    # test _reset_derivative_constraints
    @testset "_reset_derivative_constraints" begin 
        # test empty 
        @test InfiniteOpt._reset_derivative_constraints(dpref) isa Nothing
        # test warning 
        InfiniteOpt._set_has_derivative_constraints(pref, true) 
        @test push!(InfiniteOpt._derivative_constraint_dependencies(dref), index(cref)) isa Vector
        warn = "Support/method changes will invalidate existing derivative evaluation " *
               "constraints that have been added to the InfiniteModel. Thus, " *
               "these are being deleted."
        @test_logs (:warn, warn) InfiniteOpt._reset_derivative_constraints(dpref) isa Nothing
        @test !has_derivative_constraints(pref)
        @test !is_valid(m, cref)
    end
    # test set_derivative_method
    @testset "set_derivative_method" begin 
        # test NonGenerativeDerivativeMethod with no generative measures
        push!(InfiniteOpt._constraint_dependencies(dpref), InfOptConstraintIndex(1))
        @test set_derivative_method(pref, TestMethod()) isa Nothing
        @test derivative_method(dpref) isa TestMethod
        # test NonGenerativeDerivativeMethod with generative measures
        push!(InfiniteOpt._generative_measures(dpref), MeasureIndex(1))
        @test set_derivative_method(pref, TestMethod()) isa Nothing
        @test derivative_method(dpref) isa TestMethod
        # test GenerativeDerivativeMethod with generative measures
        @test InfiniteOpt._set_generative_support_info(dpref, info) isa Nothing 
        @test_throws ErrorException set_derivative_method(pref, OrthogonalCollocation(4))
        @test set_derivative_method(pref, OrthogonalCollocation(3)) isa Nothing
        @test derivative_method(dpref) isa OrthogonalCollocation
        # test GenerativeDerivativeMethod without generative measures
        empty!(InfiniteOpt._generative_measures(dpref))
        @test set_derivative_method(pref, OrthogonalCollocation(4)) isa Nothing
        @test derivative_method(dpref) isa OrthogonalCollocation
    end
    @testset "set_all_derivative_methods" begin
        @test set_all_derivative_methods(m, TestMethod()) isa Nothing
        @test derivative_method(dpref) isa TestMethod
    end
end

# Test parameter domain methods
@testset "Infinite Domain" begin
    m = InfiniteModel()
    @infinite_parameter(m, pref_gen in [0, 1])
    pref_disp = dispatch_variable_ref(pref_gen)
    bad = Bad()
    bad_pref = IndependentParameterRef(m, IndependentParameterIndex(-1))
    # _parameter_domain
    @testset "_parameter_domain" begin
        @test InfiniteOpt._parameter_domain(pref_disp) == IntervalDomain(0, 1)
        @test_throws ErrorException InfiniteOpt._parameter_domain(bad_pref)
    end
    # _update_parameter_domain
    @testset "_update_parameter_domain " begin
        push!(InfiniteOpt._constraint_dependencies(pref_disp), InfOptConstraintIndex(1))
        @test isa(InfiniteOpt._update_parameter_domain(pref_disp,
                                                    IntervalDomain(1, 2)), Nothing)
        @test InfiniteOpt._parameter_domain(pref_disp) == IntervalDomain(1, 2)
    end
    # infinite_domain
    @testset "infinite_domain" begin
        @test_throws ArgumentError infinite_domain(bad)
        @test infinite_domain(pref_disp) == IntervalDomain(1, 2)
        @test infinite_domain(pref_gen) == IntervalDomain(1, 2)
        @test_throws ErrorException infinite_domain(bad_pref)
    end
    # set_infinite_domain
    @testset "set_infinite_domain" begin
        @test_throws ArgumentError set_infinite_domain(bad, IntervalDomain(0, 1))
        @test isa(set_infinite_domain(pref_disp, IntervalDomain(2, 3)), Nothing)
        @test infinite_domain(pref_disp) == IntervalDomain(2, 3)
        @test isa(set_infinite_domain(pref_gen, IntervalDomain(1, 3)), Nothing)
        @test infinite_domain(pref_gen) == IntervalDomain(1, 3)
        @test set_infinite_domain(pref_gen, UniDistributionDomain(Normal())) isa Nothing
        @test infinite_domain(pref_disp) isa UniDistributionDomain 
        push!(InfiniteOpt._data_object(pref_gen).measure_indices, MeasureIndex(1))
        @test_throws ErrorException set_infinite_domain(pref_gen, IntervalDomain(1, 3))
    end
end

# Test parameter support methods
@testset "Supports" begin
    m = InfiniteModel()
    @infinite_parameter(m, pref in [0, 1], sig_digits = 5)
    @infinite_parameter(m, pref2 in [0, 1])
    pref_disp = dispatch_variable_ref(pref)
    bad = Bad()
    push!(InfiniteOpt._data_object(pref).constraint_indices, InfOptConstraintIndex(1))
    # _parameter_supports
    @testset "_parameter_supports" begin
        @test InfiniteOpt._parameter_supports(pref_disp) == SortedDict{Float64, Set{DataType}}()
    end
    @testset "_parameter_support_values" begin
        @test InfiniteOpt._parameter_support_values(pref_disp) == Float64[]
    end
    # _update_parameter_supports
    @testset "_update_parameter_supports " begin
        dict = SortedDict{Float64, Set{DataType}}(1. => Set{DataType}([MCSample]))
        @test isa(InfiniteOpt._update_parameter_supports(pref_disp, dict), Nothing)
        @test InfiniteOpt._parameter_support_values(pref_disp) == [1.]
    end
    # significant_digits
    @testset "significant_digits" begin
        @test significant_digits(pref) == 5
    end
    # num_supports
    @testset "num_supports" begin
        @test_throws ArgumentError num_supports(bad)
        @test num_supports(pref_disp) == 1
        @test num_supports(pref_disp, label = UserDefined) == 0
        @test num_supports(pref) == 1
        @test num_supports(pref, label = UserDefined) == 0
        @test num_supports(pref, label = SampleLabel) == 1
        @test num_supports(pref, label = All) == 1
    end
    # has_supports
    @testset "has_supports" begin
        @test_throws ArgumentError has_supports(bad)
        @test has_supports(pref_disp)
        @test has_supports(pref)
        InfiniteOpt._update_parameter_supports(pref_disp, SortedDict{Float64, Set{DataType}}())
        @test !has_supports(pref_disp)
        @test !has_supports(pref)
    end
    # supports
    @testset "supports" begin
        @test supports(pref_disp) == []
        @test supports(pref) == []
        dict = SortedDict{Float64, Set{DataType}}(1. => Set{DataType}([MCSample]))
        InfiniteOpt._update_parameter_supports(pref_disp, dict)
        @test supports(pref_disp) == [1.]
        @test supports(pref) == [1.]
        @test supports(pref, label = MCSample) == [1.]
        @test supports(pref, label = All) == [1.]
    end
    # supports (vector)
    @testset "supports (vector)" begin
        # test simple case
        @test supports([pref_disp], label = SampleLabel) == ones(1, 1)
        # test typical combinatorial case
        supps = [[-1, 0, 1], [-1, 1]] 
        @infinite_parameter(m, x[i = 1:2] in [-1, 1], supports = supps[i], independent = true)
        @test supports(x) == Float64[-1 0 1 -1 0 1; -1 -1 -1 1 1 1]
        # test non-combinatorial case 
        @test_throws ErrorException supports(x, use_combinatorics = false)
        @test set_supports(x[2], supps[1], force = true) isa Nothing 
        @test supports(x, use_combinatorics = false) == Float64[-1 0 1; -1 0 1]
    end
    # set_supports
    @testset "set_supports" begin
        @test_throws ArgumentError set_supports(bad, [0, 1])
        @test isa(set_supports(pref_disp, [0, 1], force = true), Nothing)
        @test isa(set_supports(pref, (0, 1), force = true), Nothing)
        @test supports(pref) == [0., 1.]
        @test_throws ErrorException set_supports(pref, [2, 3])
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) set_supports(pref, [1, 1], force = true)
        @test_throws ErrorException set_supports(pref, [0.5])
        @test !has_internal_supports(pref)
    end
    # add_supports
    @testset "add_supports" begin
        @test_throws ArgumentError add_supports(bad, 0.5)
        @test isa(add_supports(pref_disp, 0.25), Nothing)
        @test isa(add_supports(pref, (i for i in [0.5])), Nothing)
        @test supports(pref) == [0.25, 0.5, 1.]
        @test isa(add_supports(pref, [0, 0.25, 1], check = false), Nothing)
        @test supports(pref) == [0, 0.25, 0.5, 1.]
        @test add_supports(pref, 0.2, label = InternalGaussLobatto) isa Nothing
        @test supports(pref) == [0, 0.25, 0.5, 1.]
        @test supports(pref, label = All) == [0, 0.2, 0.25, 0.5, 1.]
        @test has_internal_supports(pref)
    end
    # delete_supports
    @testset "delete_supports" begin
        # test bad parameter 
        @test_throws ArgumentError delete_supports(bad)
        # test label deletion
        set_derivative_method(pref, OrthogonalCollocation(3))
        add_supports(pref, 0.1, label = UniformGrid)
        InfiniteOpt._set_has_generative_supports(pref, true)
        @test delete_supports(pref_disp, label = UniformGrid) isa Nothing
        @test !has_generative_supports(pref)
        @test !has_internal_supports(pref)
        @test supports(pref, label = All) == [0, 0.25, 0.5, 1.]
        # test total deletion
        @test isa(delete_supports(pref), Nothing)
        @test !has_generative_supports(pref)
        @test !has_internal_supports(pref)
        @test supports(pref) == []
        # test array input 
        @test isa(delete_supports([pref]), Nothing)
        # prepare to test derivative constraints 
        func = (x) -> NaN
        num = 0.
        info = VariableInfo(true, num, true, num, true, num, false, func, false, false)
        deriv = Derivative(info, true, pref, pref, 1)
        object = VariableData(deriv)
        idx = InfiniteOpt._add_data_object(m, object)
        push!(InfiniteOpt._derivative_dependencies(pref), idx)
        dref = DerivativeRef(m, idx)
        gdref = GeneralVariableRef(m, 1, DerivativeIndex)
        cref = @constraint(m, gdref == 0)
        InfiniteOpt._set_has_derivative_constraints(pref, true) 
        @test push!(InfiniteOpt._derivative_constraint_dependencies(dref), index(cref)) isa Vector
        # test derivative constraint warning 
        warn = "Deleting supports invalidated derivative evaluations. Thus, these " * 
               "are being deleted as well."
        @test_logs (:warn, warn) delete_supports(pref) isa Nothing
        @test !has_derivative_constraints(pref)
        @test !is_valid(m, cref)
        # test measure error
        push!(InfiniteOpt._data_object(pref).measure_indices, MeasureIndex(1))
        @test_throws ErrorException delete_supports(pref)
    end
end

# Test lower bound functions
@testset "Lower Bound" begin
    m = InfiniteModel()
    @infinite_parameter(m, pref1 in [0, 1])
    @infinite_parameter(m, pref2 ~ Normal())
    @infinite_parameter(m, pref3 in BadScalarDomain())
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
    @infinite_parameter(m, pref1 in [0, 1])
    @infinite_parameter(m, pref2 ~ Normal())
    @infinite_parameter(m, pref3 in BadScalarDomain())
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
        # test errors
        @test_macro_throws ErrorException @finite_parameter(m)
        @test_macro_throws ErrorException @finite_parameter()
        @test_macro_throws ErrorException @finite_parameter(m, a, 2)
        @test_macro_throws ErrorException @finite_parameter(m, a ~ 42)
        @test_macro_throws ErrorException @finite_parameter(Model(), 2)
        @test_macro_throws ErrorException @finite_parameter(m, "bob")
        @test_macro_throws ErrorException @finite_parameter(m, "bob" == 42)
        @test_macro_throws ErrorException @finite_parameter(m, test == 2, bob = 2)
        @test_macro_throws ErrorException @finite_parameter(m, test[m = s] == 2)
        @test_macro_throws ErrorException @finite_parameter(m, [i = s; k; j] == 2)
        # test anonymous definition
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        @test @finite_parameter(m, 42) == pref
        @test core_object(pref).value == 42
        @test name(pref) == ""
        # test vector anonymous definition
        prefs = [GeneralVariableRef(m, i, FiniteParameterIndex) for i in 2:3]
        @test @finite_parameter(m, [1:2] == 42; base_name = "a") == prefs
        @test core_object(prefs[1]).value == 42
        @test name.(prefs) == ["a[1]", "a[2]"]
        # test named definition
        pref = GeneralVariableRef(m, 4, FiniteParameterIndex)
        @test @finite_parameter(m, b == 42) == pref
        @test core_object(pref).value == 42
        @test name(pref) == "b"
        @test m[:b] == pref
        # test named vector definition
        prefs = [GeneralVariableRef(m, i, FiniteParameterIndex) for i in 5:6]
        @test @finite_parameter(m, c[i = 1:2] == [3, 7][i],
                                container = SparseAxisArray)[2] == prefs[2]
        @test core_object(prefs[2]).value == 7
        @test name(prefs[2]) == "c[2]"
    end
    # initialize the model
    m = InfiniteModel()
    # test parameter_value
    @testset "JuMP.parameter_value" begin
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        dpref = dispatch_variable_ref(pref)
        @test_throws ArgumentError parameter_value(bad)
        @test @finite_parameter(m, g == 1) == pref
        @test parameter_value(dpref) == 1
        @test parameter_value(pref) == 1
    end
    # test JuMP.set_parameter_value
    @testset "JuMP.set_parameter_value" begin
        pref = GeneralVariableRef(m, 1, FiniteParameterIndex)
        dpref = dispatch_variable_ref(pref)
        push!(InfiniteOpt._constraint_dependencies(dpref), InfOptConstraintIndex(1))
        @test_throws ArgumentError set_parameter_value(bad, 42)
        @test isa(set_parameter_value(dpref, 42), Nothing)
        @test parameter_value(pref) == 42
        @test isa(set_parameter_value(pref, 41), Nothing)
        @test parameter_value(pref) == 41
        @test_deprecated isa(set_value(pref, 40), Nothing)
        # test resetting update
        set_transformation_backend_ready(m, true)
        @test isa(set_parameter_value(pref, 39), Nothing)
        @test parameter_value(pref) == 39
        @test !transformation_backend_ready(m)
    end
end

# Test support flll-in and geneartion functions
@testset "Support Fill-in and Generation" begin
    @testset "generate_and_add_supports! (AbstractInfiniteDomain)" begin
        m = InfiniteModel()
        gvref1 = @infinite_parameter(m, a in [0, 1])
        pref1 = dispatch_variable_ref(gvref1)
        domain1 = infinite_domain(pref1)
        dist = Normal(0., 1.)
        gvref2 = @infinite_parameter(m, c ~ dist)
        pref2 = dispatch_variable_ref(gvref2)
        domain2 = infinite_domain(pref2)
        @test generate_and_add_supports!(pref1, domain1, num_supports = 4) isa Nothing
        @test generate_and_add_supports!(pref2, domain2, num_supports = 4) isa Nothing
        @test length(supports(pref1)) == 4
        @test length(supports(pref2)) == 4
    end
    # fill_in_supports! (ParameterRef)
    @testset "fill_in_supports! (ParameterRef)" begin
        m = InfiniteModel()
        pref1 = @infinite_parameter(m, a in [0, 1])
        pref2 = @infinite_parameter(m, b[1:2] in [0, 1], independent = true)
        dist = Normal(0., 1.)
        pref3 = @infinite_parameter(m, c ~ dist, supports = [-0.5, 0.5])
        @test fill_in_supports!(pref1, num_supports = 11) isa Nothing
        @test fill_in_supports!.(pref2, num_supports = 11) isa Array{Nothing}
        @test fill_in_supports!(pref3, num_supports = 11) isa Nothing
        @test length(supports(pref1)) == 11
        @test length(supports(pref2[1])) == 11
        @test length(supports(pref2[2])) == 11
        @test length(supports(pref3)) == 11
        @test -0.5 in supports(pref3)
        @test 0.5 in supports(pref3)
        @test fill_in_supports!(pref1, num_supports = 20) isa Nothing
        @test length(supports(pref1)) == 20
    end
    # test sigfig changes on infinite domain
    @testset "Sigfig support addition test" begin
        m = InfiniteModel()
        supps = [0.8236475079774124, 0.9103565379264364]
        @infinite_parameter(m, p in [supps[1], supps[2]])
        @test fill_in_supports!(p, num_supports = 8) isa Nothing
        @test length(supports(p)) == 8
    end
end
