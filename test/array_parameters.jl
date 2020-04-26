# Test core DispatchVariableRef extensions
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    obj_idx = DependentParametersIndex(1)
    idx = DependentParameterIndex(obj_idx, 1)
    set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 1)])
    params = DependentParameters(set, zeros(Float64, 2, 0), Set{Symbol}[])
    object = MultiParameterData(params, 1, 1:2, ["p1", "p2"])
    pref = DependentParameterRef(m, idx)
    gvref = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
    # test dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == pref
        @test dispatch_variable_ref(gvref) == pref
    end
    # test _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == obj_idx
    end
    # test _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(pref) == m.dependent_params
        @test InfiniteOpt._data_dictionary(gvref) == m.dependent_params
    end
    # test _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(pref) == object
        @test InfiniteOpt._data_object(gvref) == object
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(pref) == params
        @test InfiniteOpt._core_variable_object(gvref) == params
    end
    # test _num_parameters
    @testset "_num_parameters" begin
        @test InfiniteOpt._num_parameters(pref) == 2
    end
end

# Test Definition
@testset "Definition" begin
    # Setup data
    m = InfiniteModel();
    sset1 = IntervalSet(0, 1)
    sset2 = UniDistributionSet(Uniform())
    set1 = CollectionSet([sset1, sset1])
    set2 = MultiDistributionSet(MvNormal(ones(2)))
    set3 = MultiDistributionSet(MatrixBeta(2, 2, 2))
    set4 = CollectionSet([sset1, sset2])
    # test _DependentParameter
    @testset "_DependentParameter" begin
        @test InfiniteOpt._DependentParameter isa UnionAll
        @test isa(InfiniteOpt._DependentParameter(sset1, Float64[], "p"),
                  InfiniteOpt._DependentParameter{IntervalSet})
        @test isa(InfiniteOpt._DependentParameter(set1, 1, "p"),
                  InfiniteOpt._DependentParameter{<:CollectionSet})
        @test isa(InfiniteOpt._DependentParameter(set3, [1, 3], "p"),
                  InfiniteOpt._DependentParameter{<:MultiDistributionSet})
    end
    # prepare containers of _DependentParameters
    raw_params1 = [InfiniteOpt._DependentParameter(sset1, Float64[], "p") for i = 1:2]
    raw_params2 = [InfiniteOpt._DependentParameter(set1, [0, 1], "p") for i = 1:2]
    raw_params3 = [InfiniteOpt._DependentParameter(set2, 1, "p") for i = 1:2]
    raw_params4 = [InfiniteOpt._DependentParameter(set3, Float64[], "p") for i = 1:2]
    raw_params5 = [InfiniteOpt._DependentParameter(set3, Float64[], "p")
                   for i = CartesianIndices((1:2, 1:2))]
    raw_params6 = [InfiniteOpt._DependentParameter(set4, Float64[], "p") for i = 1:2]
    raw_params7 = [InfiniteOpt._DependentParameter(set2, Float64[], "p"),
                   InfiniteOpt._DependentParameter(set3, Float64[], "p")]
    raw_params8 = [InfiniteOpt._DependentParameter(set1, [0, 1], "p") for i = 1:3]
    raw_params9 = [InfiniteOpt._DependentParameter(set1, Float64[], "p"),
                   InfiniteOpt._DependentParameter(set4, Float64[], "p")]
    raw_params10 = [InfiniteOpt._DependentParameter(set1, Float64[], "p"),
                    InfiniteOpt._DependentParameter(set2, Float64[], "p")]
    raw_params11 = [InfiniteOpt._DependentParameter(BadArraySet(), 1, "p") for i = 1:2]
    raw_params12 = [InfiniteOpt._DependentParameter(sset1, Float64[], "p"),
                    InfiniteOpt._DependentParameter(sset2, Float64[], "p")]
    # test _check_param_sets (InfiniteScalarSet)
    @testset "_check_param_sets (InfiniteScalarSet)" begin
        @test InfiniteOpt._check_param_sets(error, raw_params1) isa Nothing
    end
    # test _check_param_sets (General MultiDistributionSet)
    @testset "_check_param_sets (General MultiDistributionSet)" begin
        @test InfiniteOpt._check_param_sets(error, raw_params5) isa Nothing
        @test InfiniteOpt._check_param_sets(error, raw_params3) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, raw_params4)
    end
    # test _check_param_sets (SparseAxisArray MultiDistributionSet)
    @testset "_check_param_sets (SparseAxisArray MultiDistributionSet)" begin
        new_params = convert(JuMPC.SparseAxisArray, raw_params5)
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, new_params)
    end
    # test _check_param_sets (CollectionSet)
    @testset "_check_param_sets (CollectionSet)" begin
        @test InfiniteOpt._check_param_sets(error, raw_params2) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, raw_params8)
        new_params = convert(JuMPC.SparseAxisArray, raw_params2)
        warn = "CollectionSet order may not match the given `SparseAxisArray` " *
               "of specified dependent infinite parameters, consider instead " *
               "specifying the `InfiniteScalarSet` for each parameter using " *
               "the `set` keyword and the appropriate indices."
        @test_logs (:warn, warn) InfiniteOpt._check_param_sets(error, new_params)
    end
    # test _check_param_sets (InfiniteArraySet)
    @testset "_check_param_sets (InfiniteArraySet)" begin
        @test InfiniteOpt._check_param_sets(error, raw_params11) isa Nothing
    end
    # test _check_param_sets (Mixed)
    @testset "_check_param_sets (Mixed)" begin
        @test InfiniteOpt._check_param_sets(error, raw_params12) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, raw_params10)
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, raw_params7)
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, raw_params9)
    end
    # test _check_param_sets (Fallback)
    @testset "_check_param_sets (Fallback)" begin
        @test_throws ErrorException InfiniteOpt._check_param_sets(error, ones(2))
    end
    # test _make_array_set (InfiniteArraySet)
    @testset "_make_array_set (InfiniteArraySet)" begin
        @test InfiniteOpt._make_array_set(raw_params3) == set2
        @test InfiniteOpt._make_array_set(raw_params2) == set1
    end
    # test _make_array_set (InfiniteScalarSet)
    @testset "_make_array_set (InfiniteScalarSet)" begin
        @test InfiniteOpt._make_array_set(raw_params1) isa CollectionSet{IntervalSet}
    end
    # test _build_parameters
    @testset "_build_parameters" begin
        # test errors
        @test_throws ErrorException InfiniteOpt._build_parameters(error, raw_params1, bob = 2)
        @test_throws ErrorException InfiniteOpt._build_parameters(error, raw_params4)
        @test_throws ErrorException InfiniteOpt._build_parameters(error, raw_params9)
        raw_params12 = [InfiniteOpt._DependentParameter(set1, 0, "p"),
                        InfiniteOpt._DependentParameter(set1, Float64[], "p")]
        @test_throws ErrorException InfiniteOpt._build_parameters(error, raw_params12)
        raw_params13 = [InfiniteOpt._DependentParameter(set1, [0, 2], "p") for i = 1:2]
        @test_throws ErrorException InfiniteOpt._build_parameters(error, raw_params13)
        # test has supports
        @test InfiniteOpt._build_parameters(error,
                            raw_params2)[1] isa DependentParameters
        @test InfiniteOpt._build_parameters(error, raw_params2)[1].set == set1
        @test InfiniteOpt._build_parameters(error,
                            raw_params2)[1].supports == Float64[0 1; 0 1]
        @test InfiniteOpt._build_parameters(error,
                            raw_params2)[1].labels == [Set([UserDefined]) for i = 1:2]
        @test InfiniteOpt._build_parameters(error, raw_params2)[2] == ["p", "p"]
        @test InfiniteOpt._build_parameters(error,
                            raw_params2)[3] == CartesianIndices(1:2)
        # test support generation
        @test InfiniteOpt._build_parameters(error, raw_params6,
                            num_supports = 4)[1] isa DependentParameters
        @test InfiniteOpt._build_parameters(error, raw_params6,
                            num_supports = 4)[1].set == set4
        @test InfiniteOpt._build_parameters(error, raw_params6,
                            num_supports = 4)[1].supports isa Array{Float64, 2}
        @test InfiniteOpt._build_parameters(error, raw_params6,
                            num_supports = 4)[1].labels == [Set([Mixture]) for i = 1:4]
        @test InfiniteOpt._build_parameters(error, raw_params2,
                            num_supports = 4)[2] == ["p", "p"]
        @test InfiniteOpt._build_parameters(error, raw_params2,
                            num_supports = 4)[3] == CartesianIndices(1:2)
        # test with no supports
        @test InfiniteOpt._build_parameters(error,
                            raw_params5)[1] isa DependentParameters
        @test InfiniteOpt._build_parameters(error, raw_params5)[1].set == set3
        @test InfiniteOpt._build_parameters(error,
                            raw_params5)[1].supports == zeros(4, 0)
        @test InfiniteOpt._build_parameters(error,
                            raw_params5)[1].labels == Set{Symbol}[]
        @test InfiniteOpt._build_parameters(error,
                            raw_params5)[2] == ["p", "p", "p", "p"]
        @test InfiniteOpt._build_parameters(error,
                            raw_params5)[3] == CartesianIndices((1:2, 1:2))
    end
    # test add_parameters
    @testset "add_parameters" begin
        # test default
        pref1 = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
        params = InfiniteOpt._build_parameters(error, raw_params1)[1]
        @test add_parameters(m, params) == [pref1, pref2]
        @test InfiniteOpt._last_object_num(m) == 1
        @test InfiniteOpt._last_param_num(m) == 2
        @test name(pref1) == "noname"
        # test vector build
        pref1 = GeneralVariableRef(m, 2, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 2, DependentParameterIndex, 2)
        inputs = InfiniteOpt._build_parameters(error, raw_params2)
        @test add_parameters(m, inputs...) == [pref1, pref2]
        @test InfiniteOpt._last_object_num(m) == 2
        @test InfiniteOpt._last_param_num(m) == 4
        @test name(pref1) == "p"
        # test array build
        pref1 = GeneralVariableRef(m, 3, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 3, DependentParameterIndex, 2)
        pref3 = GeneralVariableRef(m, 3, DependentParameterIndex, 3)
        pref4 = GeneralVariableRef(m, 3, DependentParameterIndex, 4)
        inputs = InfiniteOpt._build_parameters(error, raw_params5)
        @test add_parameters(m, inputs...) == [pref1 pref3; pref2 pref4]
        @test InfiniteOpt._last_object_num(m) == 3
        @test InfiniteOpt._last_param_num(m) == 8
        @test name(pref1) == "p"
        # test DenseAxisArray build
        pref1 = GeneralVariableRef(m, 4, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 4, DependentParameterIndex, 2)
        pref3 = GeneralVariableRef(m, 4, DependentParameterIndex, 3)
        pref4 = GeneralVariableRef(m, 4, DependentParameterIndex, 4)
        dense_params = JuMPC.DenseAxisArray(raw_params5, axes(raw_params5)...)
        inputs = InfiniteOpt._build_parameters(error, dense_params)
        expected = JuMPC.DenseAxisArray([pref1 pref3; pref2 pref4], axes(raw_params5)...)
        @test add_parameters(m, inputs...) == expected
        @test InfiniteOpt._last_object_num(m) == 4
        @test InfiniteOpt._last_param_num(m) == 12
        @test name(pref1) == "p"
        # test SparseAxisArray
        pref1 = GeneralVariableRef(m, 5, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 5, DependentParameterIndex, 2)
        sparse_params = convert(JuMPC.SparseAxisArray, raw_params1)
        inputs = InfiniteOpt._build_parameters(error, sparse_params)
        expected = convert(JuMPC.SparseAxisArray, [pref1, pref2])
        @test add_parameters(m, inputs...) == expected
        @test InfiniteOpt._last_object_num(m) == 5
        @test InfiniteOpt._last_param_num(m) == 14
        @test name(pref1) == "p"
    end
end

# test macro definition
@testset "Macro Definition" begin
    # Setup data
    m = InfiniteModel();
    dist1 = Uniform()
    dist2 = MvNormal(ones(2))
    dist3 = MatrixBeta(2, 2, 2)
    sset1 = IntervalSet(0, 1)
    sset2 = UniDistributionSet(dist1)
    set1 = CollectionSet([sset1, sset1])
    set2 = MultiDistributionSet(dist2)
    set3 = MultiDistributionSet(dist3)
    set4 = CollectionSet([sset1, sset2])
    set5 = CollectionSet([sset2, sset2])
    # test _construct_array_set
    @testset "_construct_array_set" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test_throws ErrorException InfiniteOpt._construct_array_set(error, info)
        info.has_lb = true; info.lower_bound = 0
        @test_throws ErrorException InfiniteOpt._construct_array_set(error, info)
        info.has_ub = true; info.upper_bound = 1
        check = :(isa($(info.lower_bound), Real))
        expected = :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a real number."))
        @test InfiniteOpt._construct_array_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(distribution = Normal())
        check = :(isa($(info.distribution), Distributions.UnivariateDistribution))
        expected = :($(check) ? UniDistributionSet($(info.distribution)) : MultiDistributionSet($(info.distribution)))
        @test InfiniteOpt._construct_array_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(set = IntervalSet(0, 1))
        check1 = :(isa($(info.set), AbstractInfiniteSet))
        check2 = :(isa($(info.set), Distributions.UnivariateDistribution))
        expected = :($(check1) ? $(info.set) : ($(check2) ? UniDistributionSet($(info.set)) : MultiDistributionSet($(info.set))))
        @test InfiniteOpt._construct_array_set(error, info) == expected
    end
    # test @dependent_parameters
    @testset "@dependent_parameters" begin
        # test macro errors
        @test_macro_throws ErrorException @dependent_parameters(m)
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2], 42)
        @test_macro_throws ErrorException @dependent_parameters(m, param)
        @test_macro_throws ErrorException @dependent_parameters(m, [1:2] in dist1)
        @test_macro_throws ErrorException @dependent_parameters(m, "bob"[1:2])
        # test set errors
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2],
                                            lower_bound = 3)
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2],
                                            upper_bound = 3)
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2],
                                            lower_bound = Complex(1),
                                            upper_bound = 2)
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2])
        # test build errors
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2] in dist1, bob = 42)
        @test_macro_throws ErrorException @dependent_parameters(m, p[i = 1:2] in [set1, set2][i])
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2] in dist3)
        @test_macro_throws ErrorException @dependent_parameters(m, p[i = 1:3] in set1)
        @test_macro_throws ErrorException @dependent_parameters(m, p[1:2] in set1, supports = 4)
        @test_macro_throws ErrorException @dependent_parameters(m, p[i = 1:2] in set1, supports = [[1, 0], 1][i])
        # test simple explict build
        pref1 = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
        @test @dependent_parameters(m, a[1:2] in dist1,
                num_supports = 10) == [pref1, pref2]
        @test InfiniteOpt._data_object(pref1).names == ["a[1]", "a[2]"]
        @test size(InfiniteOpt._core_variable_object(pref1).supports) == (2, 10)
        @test InfiniteOpt._core_variable_object(pref1).labels[1] == Set([McSample])
        @test InfiniteOpt._core_variable_object(pref1).set.sets == set5.sets
        # test another explicit build
        pref1 = GeneralVariableRef(m, 2, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 2, DependentParameterIndex, 2)
        expected = JuMPC.DenseAxisArray([pref1, pref2], 3:4)
        @test @dependent_parameters(m, b[3:4] in set2,
                supports = 0) == expected
        @test InfiniteOpt._data_object(pref1).names == ["b[3]", "b[4]"]
        @test InfiniteOpt._core_variable_object(pref1).supports == zeros(2, 1)
        # test explicit build with some args
        pref1 = GeneralVariableRef(m, 3, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 3, DependentParameterIndex, 2)
        expected = convert(JuMPC.SparseAxisArray, [pref1, pref2])
        @test @dependent_parameters(m, c[1:2] in sset1,
                supports = 0, base_name = "z",
                container = SparseAxisArray) == expected
        @test InfiniteOpt._data_object(pref1).names == ["z[1]", "z[2]"]
        # test explicit comparison
        pref1 = GeneralVariableRef(m, 4, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 4, DependentParameterIndex, 2)
        @test @dependent_parameters(m, 0 <= d[1:2] <= 1,
                num_supports = 10) == [pref1, pref2]
        @test InfiniteOpt._data_object(pref1).names == ["d[1]", "d[2]"]
        @test size(InfiniteOpt._core_variable_object(pref1).supports) == (2, 10)
        @test InfiniteOpt._core_variable_object(pref1).labels[1] == Set([UniformGrid])
        @test InfiniteOpt._core_variable_object(pref1).set.sets == set1.sets
        # test test anonymous
        pref1 = GeneralVariableRef(m, 5, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 5, DependentParameterIndex, 2)
        pref3 = GeneralVariableRef(m, 5, DependentParameterIndex, 3)
        pref4 = GeneralVariableRef(m, 5, DependentParameterIndex, 4)
        @test @dependent_parameters(m, [1:2, 1:2],
                 distribution = dist3) == [pref1 pref3; pref2 pref4]
        @test InfiniteOpt._data_object(pref1).names == ["noname[1,1]", "noname[2,1]", "noname[1,2]", "noname[2,2]"]
        @test InfiniteOpt._core_variable_object(pref1).labels == Set{Symbol}[]
        # test anonymous with set keyword
        pref1 = GeneralVariableRef(m, 6, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 6, DependentParameterIndex, 2)
        @test @dependent_parameters(m, [1:2], set = sset2) == [pref1, pref2]
        # test anonymous with dist keyword and base_name
        pref1 = GeneralVariableRef(m, 7, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 7, DependentParameterIndex, 2)
        @test @dependent_parameters(m, [1:2],
                 distribution = dist1, base_name = "zz") == [pref1, pref2]
        @test InfiniteOpt._data_object(pref1).names == ["zz[1]", "zz[2]"]
        # test anonymous with bounds keywords
        pref1 = GeneralVariableRef(m, 8, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 8, DependentParameterIndex, 2)
        @test @dependent_parameters(m, [i = 1:2], upper_bound = [1, 2][i],
                  lower_bound = 0) == [pref1, pref2]
    end
    # test @infinite_parameter
    @testset "@infinite_parameter" begin
        m = InfiniteModel();
        # test anonymous single
        pref = GeneralVariableRef(m, 1, IndependentParameterIndex)
        @test @infinite_parameter(m, distribution = dist1) == pref
        # test explicit single
        pref = GeneralVariableRef(m, 2, IndependentParameterIndex)
        @test @infinite_parameter(m, a in [0, 1]) == pref
        # test independent multi
        pref1 = GeneralVariableRef(m, 3, IndependentParameterIndex)
        pref2 = GeneralVariableRef(m, 4, IndependentParameterIndex)
        @test @infinite_parameter(m, b[1:2] in [0, 1], num_supports = 2,
                                  independent = true) == [pref1, pref2]
        # test dependent multi
        pref1 = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
        @test @infinite_parameter(m, c[1:2] in dist2,
                                  base_name = "bob") == [pref1, pref2]
        # test anonymous multi dependent
        pref1 = GeneralVariableRef(m, 2, DependentParameterIndex, 1)
        pref2 = GeneralVariableRef(m, 2, DependentParameterIndex, 2)
        @test @infinite_parameter(m, [1:2], set = set2,
                                  base_name = "bob") == [pref1, pref2]
        # test multi with variable independent
        pref1 = GeneralVariableRef(m, 5, IndependentParameterIndex)
        pref2 = GeneralVariableRef(m, 6, IndependentParameterIndex)
        indep = true
        @test @infinite_parameter(m, d[1:2] in [0, 1], num_supports = 2,
                                  independent = indep) == [pref1, pref2]
    end
end

# test naming stuff
@testset "Naming" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @dependent_parameters(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    # test _param_index
    @testset "_param_index" begin
        @test InfiniteOpt._param_index(prefs[1]) == 1
        @test InfiniteOpt._param_index(prefs[2]) == 2
    end
    # test name
    @testset "JuMP.name" begin
        @test name(prefs[1]) == "a[1]"
        @test name(prefs[2]) == "a[2]"
        @test name(gvrefs[1]) == "a[1]"
        @test name(gvrefs[2]) == "a[2]"
    end
    # test set_name
    @testset "JuMP.set_name" begin
        @test set_name(prefs[1], "joe") isa Nothing
        @test name(prefs[1]) == "joe"
        @test set_name(prefs[2], "joe") isa Nothing
        @test name(prefs[2]) == "joe"
        @test set_name(gvrefs[2], "joe2") isa Nothing
        @test name(gvrefs[2]) == "joe2"
    end
    # test is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, prefs[1])
        @test !is_valid(InfiniteModel(), prefs[1])
        @test is_valid(m, gvrefs[1])
    end
end

# test dependency functions
@testset "Parameter Dependencies" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @dependent_parameters(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    data = InfiniteOpt._data_object(first(prefs))
    # test _infinite_variable_dependencies
    @testset "_infinite_variable_dependencies" begin
        @test InfiniteOpt._infinite_variable_dependencies(prefs[1]) == data.infinite_var_indices[1]
        @test InfiniteOpt._infinite_variable_dependencies(prefs[2]) == data.infinite_var_indices[2]
        @test InfiniteOpt._infinite_variable_dependencies(gvrefs[1]) == data.infinite_var_indices[1]
    end
    # test _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(prefs[1]) == data.measure_indices[1]
        @test InfiniteOpt._measure_dependencies(prefs[2]) == data.measure_indices[2]
        @test InfiniteOpt._measure_dependencies(gvrefs[1]) == data.measure_indices[1]
    end
    # test _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(prefs[1]) == data.constraint_indices[1]
        @test InfiniteOpt._constraint_dependencies(prefs[2]) == data.constraint_indices[2]
        @test InfiniteOpt._constraint_dependencies(gvrefs[1]) == data.constraint_indices[1]
    end
    # test used_by_infinite_variable
    @testset "used_by_infinite_variable" begin
        # test not used
        @test !used_by_infinite_variable(prefs[1])
        @test !used_by_infinite_variable(prefs[2])
        @test !used_by_infinite_variable(gvrefs[1])
        # test used
        push!(data.infinite_var_indices[1], InfiniteVariableIndex(1))
        @test used_by_infinite_variable(prefs[1])
        # undo changes
        empty!(data.infinite_var_indices[1])
    end
    # test used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(prefs[1])
        @test !used_by_measure(prefs[2])
        @test !used_by_measure(gvrefs[1])
        # test used
        push!(data.measure_indices[1], MeasureIndex(1))
        @test used_by_measure(prefs[1])
        # undo changes
        empty!(data.measure_indices[1])
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(prefs[1])
        @test !used_by_constraint(prefs[2])
        @test !used_by_constraint(gvrefs[1])
        # test used
        push!(data.constraint_indices[1], ConstraintIndex(1))
        @test used_by_constraint(prefs[1])
        # undo changes
        empty!(data.constraint_indices[1])
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(prefs[1])
        @test !is_used(prefs[2])
        @test !is_used(gvrefs[1])
        # test used
        push!(data.infinite_var_indices[1], InfiniteVariableIndex(1))
        @test is_used(prefs[1])
        # undo changes
        empty!(data.infinite_var_indices[1])
    end
end

# test Parameter Object Methods
@testset "Object Methods" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @dependent_parameters(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    data = InfiniteOpt._data_object(first(prefs))
    set = CollectionSet([IntervalSet(0, 2), IntervalSet(0, 2)])
    params = DependentParameters(set, zeros(Float64, 2, 0), [Set{Symbol}()])
    # test _parameter_number
    @testset "_parameter_number" begin
        @test InfiniteOpt._parameter_number(prefs[1]) == 1
        @test InfiniteOpt._parameter_number(prefs[2]) == 2
        @test InfiniteOpt._parameter_number(gvrefs[1]) == 1
    end
    # test _object_number
    @testset "_object_number" begin
        @test InfiniteOpt._object_number(prefs[1]) == 1
        @test InfiniteOpt._object_number(prefs[2]) == 1
        @test InfiniteOpt._object_number(gvrefs[1]) == 1
    end
    # test _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(prefs[1]) == [1]
        @test InfiniteOpt._object_numbers(prefs[2]) == [1]
        @test InfiniteOpt._object_numbers(gvrefs[1]) == [1]
    end
    # test _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(prefs[1], params) isa Nothing
        @test InfiniteOpt._set_core_variable_object(prefs[2], params) isa Nothing
    end
end

# test Infinite Set Methods
@testset "Infinite Set Methods" begin
    # setup data
    m = InfiniteModel();
    gvrefs1 = @dependent_parameters(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @dependent_parameters(m, b[1:2, 1:2] in MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    push!(InfiniteOpt._constraint_dependencies(prefs1[1]), ConstraintIndex(1))
    # test _parameter_set (raw set)
    @testset "_parameter_set (Raw Set)" begin
        @test InfiniteOpt._parameter_set(prefs1[1]) isa CollectionSet
        @test InfiniteOpt._parameter_set(prefs2[1, 1]) isa MultiDistributionSet
    end
    # test _parameter_set (CollectionSet)
    @testset "_parameter_set (CollectionSet)" begin
        set = InfiniteOpt._parameter_set(prefs1[1])
        @test InfiniteOpt._parameter_set(set, prefs1[1]) == IntervalSet(0, 1)
    end
    # test _parameter_set (Fallback)
    @testset "_parameter_set (Fallback)" begin
        set = InfiniteOpt._parameter_set(prefs2[1])
        @test_throws ErrorException InfiniteOpt._parameter_set(set, prefs2[1])
    end
    # test infinite_set (single)
    @testset "infinite_set (Single)" begin
        @test infinite_set(prefs1[1]) == IntervalSet(0, 1)
        @test_throws ErrorException infinite_set(prefs2[1])
        @test infinite_set(gvrefs1[1]) == IntervalSet(0, 1)
    end
    # test _check_complete_param_array
    @testset "_check_complete_param_array" begin
        @test InfiniteOpt._check_complete_param_array(prefs1) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_complete_param_array(prefs2[:, 1])
    end
    # test infinite_set (Array)
    @testset "infinite_set (Array)" begin
        @test infinite_set(prefs1) isa CollectionSet
        @test_throws ErrorException infinite_set(prefs2[:, 1])
        @test infinite_set(prefs2) isa MultiDistributionSet
        @test infinite_set(gvrefs1) isa CollectionSet
    end
    # test _update_parameter_set
    @testset "_update_parameter_set" begin
        old_set = infinite_set(prefs1)
        new_set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 2)])
        @test InfiniteOpt._update_parameter_set(prefs1[1], new_set) isa Nothing
        @test num_supports(prefs1) == 0
        @test !optimizer_model_ready(m)
        @test infinite_set(prefs1) == new_set
        @test InfiniteOpt._update_parameter_set(prefs1[1], old_set) isa Nothing
    end
    # test set_infinite_set (Single)
    @testset "set_infinite_set (Single)" begin
        # test normal
        @test set_infinite_set(prefs1[2], IntervalSet(-1, 1)) isa Nothing
        @test infinite_set(prefs1[2]) == IntervalSet(-1, 1)
        @test set_infinite_set(gvrefs1[2], IntervalSet(0, 1)) isa Nothing
        # test errors
        @test_throws ErrorException set_infinite_set(prefs2[1], IntervalSet(0, 2))
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException set_infinite_set(prefs1[1], IntervalSet(0, 2))
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
        @test_throws ErrorException set_infinite_set(prefs1[1], UniDistributionSet(Normal()))
    end
    # test set_infinite_set (Array)
    @testset "set_infinite_set (Array)" begin
        # test normal
        old_set = infinite_set(prefs1)
        new_set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 2)])
        @test set_infinite_set(prefs1, new_set) isa Nothing
        @test infinite_set(prefs1) == new_set
        @test set_infinite_set(gvrefs1, old_set) isa Nothing
        # test errors
        @test_throws ErrorException set_infinite_set(prefs2[:, 1], new_set)
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException set_infinite_set(prefs1, new_set)
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
        new_set = CollectionSet([IntervalSet(0, 1), UniDistributionSet(Normal())])
        @test_throws ErrorException set_infinite_set(prefs1, new_set)
    end
    # test JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(prefs1[1])
        @test has_lower_bound(gvrefs1[2])
        @test !has_lower_bound(prefs2[2])
    end
    # test JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(prefs1[1]) == 0
        @test lower_bound(gvrefs1[2]) == 0
        @test_throws ErrorException lower_bound(prefs2[2])
    end
    # test JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test set_lower_bound(prefs1[1], -1) isa Nothing
        @test lower_bound(prefs1[1]) == -1
        @test set_lower_bound(gvrefs1[2], 0) isa Nothing
        @test lower_bound(prefs1[2]) == 0
        @test_throws ErrorException set_lower_bound(prefs2[2], -3)
    end
    # test JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(prefs1[1])
        @test has_upper_bound(gvrefs1[2])
        @test !has_upper_bound(prefs2[2])
    end
    # test JuMP.upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(prefs1[1]) == 1
        @test upper_bound(gvrefs1[2]) == 1
        @test_throws ErrorException upper_bound(prefs2[2])
    end
    # test JuMP.set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        @test set_upper_bound(prefs1[1], 42) isa Nothing
        @test upper_bound(prefs1[1]) == 42
        @test set_upper_bound(gvrefs1[2], 0) isa Nothing
        @test upper_bound(prefs1[2]) == 0
        @test_throws ErrorException set_upper_bound(prefs2[2], 42)
    end
end

# test Support Methods
@testset "Support Methods" begin
    # Setup data
    m = InfiniteModel();
    gvrefs1 = @dependent_parameters(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @dependent_parameters(m, b[1:2, 1:2] in MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    push!(InfiniteOpt._constraint_dependencies(prefs1[1]), ConstraintIndex(1))
    # test _parameter_supports
    @testset "_parameter_supports" begin
        @test InfiniteOpt._parameter_supports(prefs1[1]) == Float64[0 1; 0 1]
        @test InfiniteOpt._parameter_supports(prefs2[2]) == zeros(Float64, 4, 0)
    end
    # test _parameter_support_labels
    @testset "_parameter_support_labels" begin
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == [Set([UniformGrid]) for i in 1:2]
        @test InfiniteOpt._parameter_support_labels(prefs2[2]) == Set{Symbol}[]
    end
    # test num_supports (Single)
    @testset "num_supports (Single)" begin
        # test default
        @test num_supports(prefs1[1]) == 2
        @test num_supports(gvrefs1[2]) == 2
        @test num_supports(prefs2[2]) == 0
        # test label
        @test num_supports(prefs1[1], label = McSample) == 0
        @test num_supports(prefs1[1], label = UniformGrid) == 2
        @test num_supports(gvrefs1[2], label = McSample) == 0
    end
    # test num_supports (Array)
    @testset "num_supports (Array)" begin
        # test default
        @test num_supports(prefs1) == 2
        @test num_supports(gvrefs1) == 2
        @test num_supports(prefs2) == 0
        # test label
        @test num_supports(prefs1, label = McSample) == 0
        @test num_supports(prefs1, label = UniformGrid) == 2
        @test num_supports(gvrefs1, label = McSample) == 0
        # test error
        @test_throws ErrorException num_supports(prefs2[:, 1])
    end
    # test has_supports (Single)
    @testset "has_supports (Single)" begin
        @test has_supports(prefs1[1])
        @test has_supports(gvrefs1[2])
        @test !has_supports(prefs2[2])
    end
    # test has_supports (Array)
    @testset "has_supports (Array)" begin
        @test has_supports(prefs1)
        @test has_supports(gvrefs1)
        @test !has_supports(prefs2)
        # test error
        @test_throws ErrorException has_supports(prefs2[:, 1])
    end
    # test supports (Single)
    @testset "supports (Single)" begin
        # test default
        @test supports(prefs1[1]) == Float64[0, 1]
        @test supports(gvrefs1[2]) == Float64[0, 1]
        @test supports(prefs2[2]) == Float64[]
        # test label
        @test supports(prefs1[1], label = McSample) == Float64[]
        @test supports(prefs1[1], label = UniformGrid) == Float64[0, 1]
        @test supports(gvrefs1[2], label = McSample) == Float64[]
    end
    # test supports (AbstractArray)
    @testset "supports (AbstractArray)" begin
        # test default
        new_prefs1 = convert(JuMPC.SparseAxisArray, prefs1)
        new_gvrefs1 = convert(JuMPC.SparseAxisArray, gvrefs1)
        expected = convert(JuMPC.SparseAxisArray, [Float64[0, 1] for i = 1:2])
        expected2 = convert(JuMPC.SparseAxisArray, [Float64[] for i = 1:2])
        @test supports(new_prefs1) == expected
        @test supports(new_gvrefs1) == expected
        @test supports(prefs2) == [Float64[] for i in CartesianIndices(prefs2)]
        # test label
        @test supports(new_prefs1, label = McSample) == expected2
        @test supports(new_prefs1, label = UniformGrid) == expected
        @test supports(new_gvrefs1, label = McSample) == expected2
    end
    # test supports (Vector)
    @testset "supports (Vector)" begin
        # test default
        @test supports(prefs1) == Float64[0 1; 0 1]
        @test supports(gvrefs1) == Float64[0 1; 0 1]
        # test label
        @test supports(prefs1, label = McSample) == zeros(Float64, 2, 0)
        @test supports(prefs1, label = UniformGrid) == Float64[0 1; 0 1]
        @test supports(gvrefs1, label = McSample) == zeros(Float64, 2, 0)
        # test error
        @test_throws ErrorException supports(prefs2[:, 1])
    end
    # test _update_parameter_supports
    @testset "_update_parameter_supports" begin
        old_supports = supports(prefs1)
        @test InfiniteOpt._update_parameter_supports(prefs1, ones(Int, 2, 3),
                                  [Set([UserDefined]) for i in 1:3]) isa Nothing
        @test supports(prefs1) == ones(Float64, 2, 3)
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._update_parameter_supports(prefs1, old_supports,
                                  [Set([UniformGrid]) for i in 1:2]) isa Nothing
    end
    # test _make_support_matrix
    @testset "_make_support_matrix" begin
        # test errors
        supps = JuMPC.DenseAxisArray([[0] for i in CartesianIndices(prefs2)], 3:4, 1:2)
        @test_throws ErrorException InfiniteOpt._make_support_matrix(prefs2, supps)
        supps = [i[1] != 1 ? [0] : [0, 1] for i in CartesianIndices(prefs2)]
        @test_throws ErrorException InfiniteOpt._make_support_matrix(prefs2, supps)
        # test normal
        supps = [[1, 1, 1] for i in CartesianIndices(prefs2)]
        @test InfiniteOpt._make_support_matrix(prefs2, supps) == ones(Float64, 4, 3)
    end
    # test set_supports (Single)
    @testset "set_supports (Single)" begin
        @test_throws ErrorException set_supports(prefs1[1], [0, 1])
        @test_throws ErrorException set_supports(gvrefs1[1], [0, 1], label = Mixture)
    end
    # test set_supports (Vector)
    @testset "set_supports (Vector)" begin
        # test errors
        @test_throws ErrorException set_supports(prefs2[:, 1], ones(2, 2))
        @test_throws ErrorException set_supports(prefs1, ones(2, 2))
        @test_throws ErrorException set_supports(prefs1, ones(2, 2) * 2, force = true)
        # test default
        old_supports = supports(prefs1)
        @test set_supports(prefs1, ones(2, 3), force = true) isa Nothing
        @test supports(prefs1) == ones(Float64, 2, 3)
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == [Set([UserDefined]) for i = 1:3]
        # test keywords
        @test set_supports(gvrefs1, old_supports, force = true,
                           label = UniformGrid) isa Nothing
        @test supports(prefs1) == old_supports
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == [Set([UniformGrid]) for i = 1:2]
    end
    # test set_supports (AbstractArray)
    @testset "set_supports (AbstractArray)" begin
        # default
        supps = [[1, 1, 1] for i in CartesianIndices(prefs2)]
        @test set_supports(prefs2, supps) isa Nothing
        @test supports(prefs2) == [ones(Float64, 3) for i in CartesianIndices(prefs2)]
        @test InfiniteOpt._parameter_support_labels(prefs2[1]) == [Set([UserDefined]) for i = 1:3]
        # test keywords
        supps = [Float64[] for i in CartesianIndices(prefs2)]
        @test set_supports(gvrefs2, supps, force = true, label = Mixture) isa Nothing
        @test supports(prefs2) == [ones(Float64, 0) for i in CartesianIndices(prefs2)]
        @test InfiniteOpt._parameter_support_labels(prefs2[1]) == Set{Symbol}[]
    end
    # test add_supports (Single)
    @testset "add_supports (Single)" begin
        @test_throws ErrorException add_supports(prefs1[1], [0, 1])
        @test_throws ErrorException add_supports(gvrefs1[1], [0, 1], label = Mixture)
    end
    # test add_supports (Vector)
    @testset "add_supports (Vector)" begin
        # test errors
        @test_throws ErrorException add_supports(prefs2[:, 1], ones(2, 2))
        @test_throws ErrorException add_supports(prefs1, ones(2, 2) * 2)
        # test default
        @test add_supports(prefs1, ones(2, 1)) isa Nothing
        @test supports(prefs1) == Float64[0 1 1; 0 1 1]
        expected = [Set([UniformGrid]), Set([UniformGrid]), Set([UserDefined])]
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == expected
        # test keywords
        @test add_supports(gvrefs1, zeros(2, 1), check = false,
                           label = UniformGrid) isa Nothing
        @test supports(prefs1) == Float64[0 1 1 0; 0 1 1 0]
        push!(expected, Set([UniformGrid]))
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == expected
    end
    # test add_supports (AbstractArray)
    @testset "add_supports (AbstractArray)" begin
        # default
        supps = [[1, 1, 1] for i in CartesianIndices(prefs2)]
        @test add_supports(prefs2, supps) isa Nothing
        @test supports(prefs2) == [ones(Float64, 3) for i in CartesianIndices(prefs2)]
        expected = [Set([UserDefined]) for i = 1:3]
        @test InfiniteOpt._parameter_support_labels(prefs2[1]) == expected
        # test keywords
        supps = [[1] for i in CartesianIndices(prefs2)]
        @test add_supports(gvrefs2, supps, check = false, label = McSample) isa Nothing
        @test supports(prefs2) == [ones(Float64, 4) for i in CartesianIndices(prefs2)]
        push!(expected, Set([McSample]))
        @test InfiniteOpt._parameter_support_labels(prefs2[1]) == expected
    end
    # test delete_supports (Single)
    @testset "delete_supports (Single)" begin
        @test_throws ErrorException delete_supports(prefs1[1])
        @test_throws ErrorException delete_supports(gvrefs1[1])
    end
    # test delete_supports (AbstractArray)
    @testset "delete_supports (AbstractArray)" begin
        # test errors
        @test_throws ErrorException delete_supports(prefs2[:, 1])
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException delete_supports(prefs1)
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
        # normal
        @test delete_supports(prefs1) isa Nothing
        @test supports(prefs1) == zeros(Float64, 2, 0)
        @test delete_supports(gvrefs2) isa Nothing
        @test supports(prefs2) == [Float64[] for i in CartesianIndices(prefs2)]
    end
end

# test Support Filling
@testset "Support Filling" begin
    # Setup data
    m = InfiniteModel();
    gvrefs1 = @dependent_parameters(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @dependent_parameters(m, b[1:2, 1:2] in MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    pref = @independent_parameter(m, c in [0, 1])
    # test generate_and_add_supports!
    @testset "generate_and_add_supports!" begin
        old_supports = supports(prefs1)
        set = infinite_set(prefs1)
        @test generate_and_add_supports!(prefs1, set, num_supports = 2) isa Nothing
        @test supports(prefs1) == Float64[0 1 0 1; 0 1 0 1]
        expected = [Set([UniformGrid]) for i = 1:4]
        @test InfiniteOpt._parameter_support_labels(prefs1[1]) == expected
        @test set_supports(prefs1, old_supports, force = true) isa Nothing
    end
    # test fill_in_supports! (Single)
    @testset "fill_in_supports! (Single)" begin
        @test_throws ErrorException fill_in_supports!(prefs1[1])
        @test_throws ErrorException fill_in_supports!(gvrefs1[1])
    end
    # test fill_in_supports! (AbstractArray)
    @testset "fill_in_supports! (AbstractArray)" begin
        # test error
        @test_throws ErrorException fill_in_supports!(prefs2[1, :])
        # test default
        @test fill_in_supports!(prefs1) isa Nothing
        @test num_supports(prefs1) == 10
        @test fill_in_supports!(gvrefs2) isa Nothing
        @test num_supports(prefs2) == 10
        # test keywords
        @test fill_in_supports!(prefs1, num_supports = 2, sig_figs = 3) isa Nothing
        @test num_supports(prefs1) == 10
        @test fill_in_supports!(gvrefs2, modify = false, num_supports = 20,
                                sig_figs = 3) isa Nothing
        @test num_supports(prefs2) == 10
        # delete additions
        @test delete_supports(prefs1) isa Nothing
        @test delete_supports(prefs2) isa Nothing
    end
    # test fill_in_supports! (InfiniteModel)
    @testset "fill_in_supports! (InfiniteModel)" begin
        # TODO implement when independent parameter fill_in_supports! is done
    end
end
