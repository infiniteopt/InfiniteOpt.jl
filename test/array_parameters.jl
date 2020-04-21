# Test core DispatchVariableRef extensions
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    obj_idx = DependentParametersIndex(1)
    idx = DependentParameterIndex(obj_idx, 1)
    set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 1)])
    params = DependentParameters(set, zeros(2, 0), Set{Symbol}[])
    object = MultiParameterData(params, 1, 1:2, ["p1", "p2"])
    pref = DependentParameterRef(m, idx)
    # test dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == pref
    end
    # test _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == obj_idx
    end
    # test _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(pref) == m.dependent_params
    end
    # test _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(pref) == object
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(pref) == params
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
end

# test naming stuff
@testset "Naming" begin
    # setup data
    m = InfiniteModel();
    prefs = @dependent_parameters(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(prefs)
    # test _param_index
    @testset "_param_index" begin
        @test InfiniteOpt._param_index(prefs[1]) == 1
        @test InfiniteOpt._param_index(prefs[2]) == 2
    end
    # test name
    @testset "JuMP.name" begin
        @test name(prefs[1]) == "a[1]"
        @test name(prefs[2]) == "a[2]"
    end
    # test set_name
    @testset "JuMP.set_name" begin
        @test set_name(prefs[1], "joe") isa Nothing
        @test name(prefs[1]) == "joe"
        @test set_name(prefs[2], "joe") isa Nothing
        @test name(prefs[2]) == "joe"
    end
end
