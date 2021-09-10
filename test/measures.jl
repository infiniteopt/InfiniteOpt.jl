# Test Core data accessors
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    w = t -> 1
    meas_data = DiscreteMeasureData(t, ones(2), ones(2), All, w, NaN, NaN, false)
    meas = Measure(zero(AffExpr),meas_data, [1], [1], false)
    object = MeasureData(meas)
    idx = MeasureIndex(1)
    mref = MeasureRef(m, idx)
    gvref = GeneralVariableRef(m, 1, MeasureIndex)
    # test dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test isequal(dispatch_variable_ref(gvref), mref)
    end
    # test _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # test _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(mref) === m.measures
        @test InfiniteOpt._data_dictionary(m, Measure) === m.measures
    end
    # test _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(mref) == object
        @test InfiniteOpt._data_object(gvref) == object
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(mref) == meas
        @test InfiniteOpt._core_variable_object(gvref) == meas
    end
    # test _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(mref) == [1]
        @test InfiniteOpt._parameter_numbers(gvref) == [1]
    end
    # test _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(mref) == [1]
        @test InfiniteOpt._object_numbers(gvref) == [1]
    end
    # test _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(mref, meas) isa Nothing
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(mref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(mref) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == InfOptConstraintIndex[]
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(mref) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(gvref) == DerivativeIndex[]
    end
    # test _delete_data_object
    @testset "_delete_data_object" begin
        @test is_valid(m, mref)
        @test InfiniteOpt._delete_data_object(mref) isa Nothing
        @test !is_valid(m, mref)
    end
end

# Test measure data constructors
@testset "Measure Data Constructors" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1]) # Vector
    @infinite_parameter(m, pars2[1:2, 1:3] in [0, 1]) # Matrix
    @infinite_parameter(m, pars3[2:3] in [0, 1]) # DenseAxisArray
    @infinite_parameter(m, pars4[i = 1:2, j = 1:2; i <= j] in [0, 1]) # SparseAxisArray
    @finite_parameter(m, fpar == 42)
    coeff_func(x) = 1
    # test default weight function
    @testset "default_weight" begin
        @test default_weight(42) == 1
        @test default_weight(zeros(10)) == 1
    end
    # test DiscreteMeasureData constructor for scalar parameters
    @testset "DiscreteMeasureData (scalar)" begin
        # test normal usage
        @test isa(DiscreteMeasureData(par, [1], [1]), DiscreteMeasureData)
        @test DiscreteMeasureData(par, [1], [1]).weight_function(42) == 1
        # test errors
        @test_throws ErrorException DiscreteMeasureData(par, [1], Int[])
        @test_throws ErrorException DiscreteMeasureData(par, [1], [2])
        @test_throws ErrorException DiscreteMeasureData(fpar, [1], [1])
    end
    # test scalar FunctionalDiscreteMeasureData constructor
    @testset "FunctionalDiscreteMeasureData (scalar)" begin
        @test isa(FunctionalDiscreteMeasureData(par, coeff_func, 10, UniqueMeasure{:a}), FunctionalDiscreteMeasureData)
        @test_throws ErrorException FunctionalDiscreteMeasureData(pars[1], coeff_func, 10, UniqueMeasure{:a})
        @test_throws ErrorException FunctionalDiscreteMeasureData(par, coeff_func, -1, UniqueMeasure{:a})
        @test_throws ErrorException FunctionalDiscreteMeasureData(par, coeff_func, 10, MCSample, lower_bound = -1, upper_bound = 2)
        @test isa(FunctionalDiscreteMeasureData(par, coeff_func, 10, UniqueMeasure{:a}, 
                                                generative_support_info = TestGenInfo()), FunctionalDiscreteMeasureData)
    end
    # test _check_multidim_params
    @testset "_check_params" begin
        @test_throws ErrorException InfiniteOpt._check_params(fpar)
        @test_throws ErrorException InfiniteOpt._check_params([par, pars[1]])
        @test_throws ErrorException InfiniteOpt._check_params([pars[1], pars2[1]])
        @test_throws ErrorException InfiniteOpt._check_params([pars4[1, 1], pars4[1, 2]])
        @test InfiniteOpt._check_params(pars) isa Nothing
    end
    # test _prepare_supports
    @testset "_prepare_supports" begin
        # vector
        supp = [[0, 0], [1, 1]]
        @test InfiniteOpt._prepare_supports(supp, IC.indices(pars)) == [0 1; 0 1]
        # array
        supp = [zeros(2, 3), ones(2, 3)]
        expected = [zeros(6, 1) ones(6, 1)]
        @test InfiniteOpt._prepare_supports(supp, IC.indices(pars2)) == expected
        # DenseAxisArray
        supp = [JuMPC.DenseAxisArray([0, 0], axes(pars3)),
                JuMPC.DenseAxisArray([1, 1], axes(pars3))]
        @test InfiniteOpt._prepare_supports(supp, IC.indices(pars3)) == [0 1; 0 1]
        # SparseAxisArrays
        supp = [JuMPC.SparseAxisArray(Dict(k => 0 for k in keys(pars4.data))),
                JuMPC.SparseAxisArray(Dict(k => 1 for k in keys(pars4.data)))]
        @test isequal(InfiniteOpt._prepare_supports(supp, IC.indices(pars4)), [0 1; 0 1; 0 1])
    end
    # test _check_supports_in_bounds (IndependentParameterRefs)
    @testset "_check_supports_in_bounds (IndependentParameterRefs)" begin
        dpar = dispatch_variable_ref(par)
        @test InfiniteOpt._check_supports_in_bounds([dpar, dpar], [0.3 0.2; 0.7 0.6]) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds([dpar, dpar], [-1 0.5; 0.7 2])
    end
    # test _check_supports_in_bounds (DependentParameterRefs)
    @testset "_check_supports_in_bounds (DependentParameterRefs)" begin
        dpars = dispatch_variable_ref.(pars)
        @test InfiniteOpt._check_supports_in_bounds(dpars, [0.3 0.2; 0.7 0.6]) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(dpars, [-1 0.5; 0.7 2])
    end
    # test multidim DiscreteMeasureData constructor
    @testset "DiscreteMeasureData (multi-dim)" begin
        # test normal usage
        @test isa(DiscreteMeasureData(pars, [1], [[1,1]]),
                  DiscreteMeasureData)
        @test DiscreteMeasureData(pars, [1], [[1, 1]]).weight_function([0, 1]) == 1
        # test errors
        @test_throws ErrorException DiscreteMeasureData(pars, [1], [[1]])
        @test_throws ErrorException DiscreteMeasureData(pars, [1, 2], [[0, 0]])
        @test_throws ErrorException DiscreteMeasureData(pars2, [1], [[0, 0, 0]])
        @test_throws ErrorException DiscreteMeasureData([fpar, fpar], [1],
                                                        [[1, 1]])
        @test_throws ErrorException DiscreteMeasureData([par, pars[1]], [1], [[1, 1]])
        @test_throws ErrorException DiscreteMeasureData(pars3, [1], [[1, 2]])
        @test_throws ErrorException DiscreteMeasureData(pars3, [1], [JuMPC.DenseAxisArray([1,2],2:3)], lower_bounds = [NaN, NaN])
    end
    # test _check_bounds_in_domain (IndependentParameter)
    @testset "_check_bounds_in_domain (independent parameters)" begin
        disp_par = dispatch_variable_ref(par)
        @test_throws ErrorException InfiniteOpt._check_bounds_in_domain([disp_par, disp_par], [-1, -1], [1, 1])
        @test_throws ErrorException InfiniteOpt._check_bounds_in_domain([disp_par, disp_par], [0, 0], [2, 2])
        @test InfiniteOpt._check_bounds_in_domain([disp_par, disp_par], [0.3, 0.2], [0.7, 0.8]) isa Nothing
    end
    # test _check_bounds_in_domain (DependentParameter)
    @testset "_check_bounds_in_domain (dependent parameters)" begin
        disp_vars = dispatch_variable_ref.(pars)
        @test_throws ErrorException InfiniteOpt._check_bounds_in_domain(disp_vars, [-1, -1], [1, 1])
        @test_throws ErrorException InfiniteOpt._check_bounds_in_domain(disp_vars, [0, 0], [2, 2])
        @test InfiniteOpt._check_bounds_in_domain(disp_vars, [0.3, 0.2], [0.7, 0.8]) isa Nothing
    end
    # test multidim FunctionalDiscreteMeasureData constructor
    @testset "FunctionalDiscreteMeasureData (multidim)" begin
        @test isa(FunctionalDiscreteMeasureData(pars, coeff_func, 10, All, lower_bounds = [0.3, 0.3], upper_bounds = [0.7, 0.7]), FunctionalDiscreteMeasureData)
        @test_throws ErrorException FunctionalDiscreteMeasureData(pars, coeff_func, -1, All)
        @test_throws ErrorException FunctionalDiscreteMeasureData(pars3, coeff_func, 5, MCSample, lower_bounds = [NaN, NaN])

    end
    # test copy methods
    @testset "Base.copy" begin
        data1 = DiscreteMeasureData(par, [1], [1])
        data2 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample)
        @test isequal(data1.parameter_refs, copy(data1).parameter_refs)
        @test data1.coefficients == copy(data1).coefficients
        @test data1.label == copy(data1).label
        @test data1.is_expect == copy(data1).is_expect
        @test data1.supports == copy(data1).supports
        @test data1.weight_function == copy(data1).weight_function
        @test isnan(copy(data1).lower_bounds)
        @test isnan(copy(data1).upper_bounds)
        @test isequal(data2.parameter_refs, copy(data2).parameter_refs)
        @test data2.coeff_function == copy(data2).coeff_function
        @test data2.label == copy(data2).label
        @test data2.is_expect == copy(data2).is_expect
        @test data2.min_num_supports == copy(data2).min_num_supports
        @test data2.weight_function == copy(data2).weight_function
        @test isnan(copy(data2).lower_bounds)
        @test isnan(copy(data2).upper_bounds)
    end
end
# Test data access methods
@testset "Data Queries" begin
    # initialize model and references
    m = InfiniteModel()
    f(x) = [1]
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    dpar = dispatch_variable_ref(par)
    info = UniformGenerativeInfo([0.5], InternalGaussLobatto)
    InfiniteOpt._set_generative_support_info(dpar, info)
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data3 = FunctionalDiscreteMeasureData(par, f, 10, MCSample, 
                                          generative_support_info = info)
    data4 = FunctionalDiscreteMeasureData(pars, f, 10, MCSample)
    data5 = FunctionalDiscreteMeasureData(par, f, 10, All, lower_bound = 0, 
                                          upper_bound = 1)
    data6 = FunctionalDiscreteMeasureData(pars, f, 10, MCSample, lower_bounds = [0, 0], 
                                          upper_bounds = [1, 1])
    data7 = FunctionalDiscreteMeasureData(par, f, 10, All, lower_bound = 0, 
                                          upper_bound = 0.8)
    data8 = FunctionalDiscreteMeasureData(par, f, 10, All)
    # parameter_refs (DiscreteMeasureData)
    @testset "parameter_refs (Single)" begin
        @test isequal(parameter_refs(data), par)
    end
    # parameter_refs (MultiDiscreteMeasureData)
    @testset "parameter_refs (Multi)" begin
        @test isequal(parameter_refs(data2), pars)
    end
    # parameter_refs (Fallback)
    @testset "parameter_refs (Fallback)" begin
        @test_throws ErrorException parameter_refs(BadData())
    end
    # support_label
    @testset "support_label" begin
        @test support_label(data) <:AbstractSupportLabel
        @test_throws ErrorException support_label(BadData())
    end
    # JuMP.lower_bound for AbstractMeasureData
    @testset "JuMP.lower_bound" begin
        @test isnan(JuMP.lower_bound(data))
        @test isnan(JuMP.lower_bound(BadData()))
    end
    # JuMP.upper_bound for AbstractMeasureData
    @testset "JuMP.upper_bound" begin
        @test isnan(JuMP.upper_bound(data))
        @test isnan(JuMP.upper_bound(BadData()))
    end
    # generative_support_info
    @testset "generative_support_info" begin
        @test generative_support_info(data) == NoGenerativeSupports()
        @test generative_support_info(data3) == info
    end
    # _is_expect
    @testset "_is_expect" begin
        @test !InfiniteOpt._is_expect(data)
        @test !InfiniteOpt._is_expect(BadData())
    end
    # supports (1D DiscreteMeasureData)
    @testset "supports (Single)" begin
        @test supports(data) == Float64[1]
    end
    # supports (DiscreteMeasureData)
    @testset "supports (Multi)" begin
        @test supports(data2) == ones(Float64, 2, 1)
    end
    # test _get_supports 
    @testset "_get_supports" begin 
        # setup the supports 
        add_supports(dpar, 0, label = MCSample)
        add_supports(dpar, 1, label = InternalLabel)
        add_generative_supports(dpar)
        @test supports(dpar, label = All) == [0, 0.5, 1]
        # test include generative and uses All
        @test InfiniteOpt._get_supports(dpar, Val(true), info, All) == [0, 0.5, 1]
        # test include generative and doesn't use All
        @test InfiniteOpt._get_supports(dpar, Val(true), info, MCSample) == [0, 0.5]
        # test exclude generative and uses All 
        @test InfiniteOpt._get_supports(dpar, Val(false), info, All) == [0, 1]
        # test exclude generative and doesn't use All 
        @test InfiniteOpt._get_supports(dpar, Val(false), info, MCSample) == [0]
        # delete the supports 
        delete_supports(dpar, label = All)
        @test !has_generative_supports(dpar)
    end
    # supports (1D FunctionalDiscreteMeasureData)
    @testset "supports (Single Functional)" begin
        # add supports 
        add_supports(dpar, [0, 1], label = MCSample)
        add_generative_supports(dpar)
        # test with non generative
        @test supports(data5) == [0, 0.5, 1]
        @test supports(data7) == [0, 0.5]
        @test supports(data8) == [0, 0.5, 1]
        # test generative
        @test supports(data3, include_generative = false) == [0, 1]
        @test supports(data3) == [0, 0.5, 1]
        # remove the supports 
        delete_supports(dpar, label = All)
        @test !has_generative_supports(dpar)
    end
    # supports (FunctionalDiscreteMeasureData)
    @testset "supports (Multi Functional)" begin
        @test supports(data4) == zeros(2, 0)
        @test supports(data6) == zeros(2, 0)
    end
    # supports (Fallback)
    @testset "supports (Fallback)" begin
        @test supports(BadData()) == Float64[]
    end
    # num_supports
    @testset "num_supports" begin
        @test num_supports(data) == 1
        @test num_supports(data3) == 0
        @test num_supports(BadData()) == 0
    end
    # min_num_supports (FunctionalDiscreteMeasureData)
    @testset "min_num_supports" begin
        @test min_num_supports(data3) == 10
        @test min_num_supports(data) == 1
    end
    # coefficient_function (FunctionalDiscreteMeasureData)
    @testset "coefficient_function" begin
        @test coefficient_function(data3) == f
        @test_throws ErrorException coefficient_function(BadData())
    end
    # coefficients (DiscreteMeasureData)
    @testset "coefficients (Single)" begin
        @test coefficients(data) == Float64[1]
    end
    # coefficients (MultiDiscreteMeasureData)
    @testset "coefficients (Multi)" begin
        @test coefficients(data2) == Float64[1]
    end
    # coefficients (FunctionalDiscreteMeasureData)
    @testset "coefficients (Functional)" begin
        @test coefficients(data3) == [1]
    end
    # coefficients (Fallback)
    @testset "coefficients (Fallback)" begin
        @test coefficients(BadData()) == Float64[]
    end
    # weight_function (DiscreteMeasureData)
    @testset "weight_function (Single)" begin
        @test weight_function(data) == default_weight
    end
    # weight_function (MultiDiscreteMeasureData)
    @testset "weight_function (Multi)" begin
        @test weight_function(data2) == default_weight
    end
    # weight_function (Fallback)
    # TODO: this might need to be changed
    @testset "weight_function (Fallback)" begin
        @test_throws ErrorException weight_function(BadData())
    end
end

# Test measure construction methods
@testset "Measure Construction" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, par2))
    @variable(m, inf3, Infinite(par2))
    @variable(m, inf4, Infinite(pars))
    @variable(m, x)
    coeff_func(x) = 1
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(par2, [1], [1], lower_bound = 0, upper_bound = 1)
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data4 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample)
    data5 = FunctionalDiscreteMeasureData(pars, coeff_func, 5, MCSample)
    # test build_measure
    @testset "build_measure" begin
        @test isa(build_measure(inf, data), Measure)
        @test isa(build_measure(inf, data2), Measure)
        @test isa(build_measure(inf4, data3), Measure)
        @test isa(build_measure(inf2, data4), Measure)
        @test isa(build_measure(inf4, data5), Measure)
        @test !build_measure(inf, data).constant_func
        @test build_measure(inf, data2).constant_func
        @test !build_measure(inf4, data3).constant_func
        @test !build_measure(inf2, data4).constant_func
        @test !build_measure(inf4, data5).constant_func
    end
end

# Test add_measure and helper functions
@testset "Measure Addition" begin
    m = InfiniteModel()
    coeff_func(x) = 1
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @infinite_parameter(m, par3 in [0, 1])
    @infinite_parameter(m, pars2[1:2] in [0, 1])
    @infinite_parameter(m, pars3[1:2] in [0, 1], independent = true)
    @variable(m, inf, Infinite(par))
    @variable(m, x)
    dpar = dispatch_variable_ref(par)
    dpars = dispatch_variable_ref.(pars)
    info = UniformGenerativeInfo([0.5], InternalLabel)
    data1 = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data3 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample)
    data4 = FunctionalDiscreteMeasureData(pars, coeff_func, 5, MCSample)
    data5 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample, lower_bound = 0.3, upper_bound = 0.7)
    data6 = FunctionalDiscreteMeasureData(pars[1], coeff_func, 10, MCSample, NoGenerativeSupports(),
                                          default_weight, NaN, NaN, false)
    data7 = FunctionalDiscreteMeasureData(par, coeff_func, 0, All, lower_bound = 0.3, upper_bound = 0.7)
    data8 = FunctionalDiscreteMeasureData(pars, coeff_func, 0, All, lower_bounds = [.25, 0], upper_bounds = [.5,1])
    data9 = FunctionalDiscreteMeasureData(pars3, coeff_func, 0, All, lower_bounds = [.25, 0], upper_bounds = [.5,1])
    data10 = FunctionalDiscreteMeasureData(pars, coeff_func, 0, All, lower_bounds = [0, 0], upper_bounds = [1,NaN])
    data11 = FunctionalDiscreteMeasureData(pars3, coeff_func, 0, All, lower_bounds = [0, 0], upper_bounds = [1,NaN])
#    meas = Measure(par + 2inf - x, data, [1], [1], false)
    # test _check_and_set_generative_info 
    @testset "_check_and_set_generative_info" begin 
        # test DependentParameterRef and NoGenerativeSupports
        @test InfiniteOpt._check_and_set_generative_info(dpars[1], NoGenerativeSupports()) isa Nothing
        # test DependentParameterRef and not NoGenerativeSupports
        @test_throws ErrorException InfiniteOpt._check_and_set_generative_info(dpars[1], TestGenInfo())
        # test IndependentParameterRef and NoGenerativeSupports
        @test InfiniteOpt._check_and_set_generative_info(dpar, NoGenerativeSupports()) isa Nothing
        # test IndependentParameterRef w/ NoGenerativeSupports
        @test InfiniteOpt._check_and_set_generative_info(dpar, NoGenerativeSupports(), TestGenInfo()) isa Nothing
        # test IndependentParameterRef w/ info already
        @test_throws ErrorException InfiniteOpt._check_and_set_generative_info(dpar, TestGenInfo(), info)
        @test InfiniteOpt._check_and_set_generative_info(dpar, info, info) isa Nothing
        # test IndependentParameterRef and not NoGenerativeSupports
        @test InfiniteOpt._check_and_set_generative_info(dpar, TestGenInfo()) isa Nothing
        # undo changes 
        @test InfiniteOpt._set_generative_support_info(dpar, NoGenerativeSupports()) isa Nothing
    end
    # test _add_supports_to_multiple_parameters (independent)
    @testset "_add_supports_to_multiple_parameters (independent)" begin
        prefs = dispatch_variable_ref.([par, par])
        @test isa(InfiniteOpt._add_supports_to_multiple_parameters(prefs, [0.1 0.2;0.3 0.4], MCSample), Nothing)
        delete_supports(par)
    end
    # test _add_supports_to_multiple_parameters (dependent)
    @testset "_add_supports_to_multiple_parameters (dependent)" begin
        prefs = dispatch_variable_ref.(pars)
        @test isa(InfiniteOpt._add_supports_to_multiple_parameters(prefs, [0. 1.; 0. 1.], MCSample), Nothing)
        delete_supports(pars)
    end
    # test add_supports_to_parameters (scalar DiscreteMeasureData)
    @testset "add_supports_to_parameters (scalar DiscreteMeasureData)" begin
        # test functionality
        @test isa(add_supports_to_parameters(data1), Nothing)
        @test supports(par) == [1.]
        # clear supports
        delete_supports(par)
    end
    # test add_supports_to_parameters (multi DiscreteMeasureData)
    @testset "add_supports_to_parameters (multi DiscreteMeasureData)" begin
        # test functionality
        @test isa(add_supports_to_parameters(data2), Nothing)
        @test supports(pars[1]) == [1.]
        @test supports(pars[2]) == [1.]
        # clear supports
        delete_supports(pars)
    end
    # test add_supports_to_parameters (scalar FunctionalDiscreteMeasureData)
    @testset "add_supports_to_parameters (scalar FunctionalDiscreteMeasureData)" begin
        # test functionality
        @test isa(add_supports_to_parameters(data3), Nothing)
        @test num_supports(par) == 5
        # clear supports
        delete_supports(par)
        # test functionality
        @test isa(add_supports_to_parameters(data5), Nothing)
        @test num_supports(par) == 5
        # test error 
        @test_throws ErrorException add_supports_to_parameters(data6)
        # clear supports
        delete_supports(par)
        add_supports_to_parameters(data7)
        @test supports(par, label = MeasureBound) == [0.3, 0.7]
        delete_supports(par)
    end
    # test _generate_multiple_functional_supports (DependentParameterRefs)
    @testset "_generate_multiple_functional_supports (DependentParameterRefs)" begin
        dpars = dispatch_variable_ref.(pars2)
        @test InfiniteOpt._generate_multiple_functional_supports(dpars, 5, MCSample, [NaN], [NaN]) isa Nothing
        @test InfiniteOpt._generate_multiple_functional_supports(dpars, 5, MCSample, [0.3, 0.3], [0.7, 0.7]) isa Nothing
        @test num_supports(pars2[1]) == 10
    end
    # test _generate_multiple_functional_supports (IndependentParameterRefs)
    @testset "_generate_multiple_functional_supports (IndependentParameterRefs)" begin
        dpars = dispatch_variable_ref.([par2, par3])
        @test InfiniteOpt._generate_multiple_functional_supports(dpars, 5, MCSample, [NaN, NaN], [NaN, NaN]) isa Nothing
        @test InfiniteOpt._generate_multiple_functional_supports(dpars, 5, MCSample, [0.3, 0.3], [0.7, 0.7]) isa Nothing
        @test num_supports(par2) == 10
    end
    # test add_supports_to_parameters (multi FunctionalDiscreteMeasureData)
    @testset "add_supports_to_parameters (multi FunctionalDiscreteMeasureData)" begin
        # test functionality
        @test isa(add_supports_to_parameters(data4), Nothing)
        @test num_supports(pars) == 5 # wrong with CollectionDomain support generation
        # clear supports
        delete_supports(pars)
        add_supports_to_parameters(data8)
        @test sortcols(supports(pars, label = MeasureBound)) ==  [0.25 0.5; 0.0 1.0]
        delete_supports(pars)
        add_supports_to_parameters(data10)
        @test supports(pars, label = MeasureBound) == zeros(2,0)
        delete_supports(pars)
        add_supports_to_parameters(data9)
        @test sortcols(supports(pars3, label = MeasureBound)) == [[.25,0] [.25, 1] [.5,0] [.5, 1]]
        delete_supports(pars3)
        add_supports_to_parameters(data11)
        @test supports(pars3, label = MeasureBound) == zeros(2,0)
        delete_supports(pars3)
    end
    # test add_supports_to_parameters (fallbacks)
    @testset "add_supports_to_parameters (fallbacks)" begin
        @test_throws ErrorException add_supports_to_parameters(BadData())
    end
    # test _update_generative_measures
    @testset "_update_generative_measures" begin
        @test InfiniteOpt._update_generative_measures(dpars, NoGenerativeSupports(), MeasureIndex(1)) isa Nothing 
        @test InfiniteOpt._update_generative_measures(dpar, info, MeasureIndex(1)) isa Nothing 
        @test InfiniteOpt._generative_measures(dpar) == [MeasureIndex(1)]
        empty!(InfiniteOpt._generative_measures(dpar))
    end
    # test add_measure
    @testset "add_measure" begin
        meas1 = Measure(par + 2inf - x, data1, [1], [1], false)
        mref1 = MeasureRef(m, MeasureIndex(1))
        @test isequal(dispatch_variable_ref(add_measure(m, meas1, "measure1")), mref1)
        @test supports(par) == [1.0]
        @test InfiniteOpt._measure_dependencies(par) == [MeasureIndex(1)]
        @test InfiniteOpt._data_object(mref1).name == "measure1"
        @test !InfiniteOpt._data_object(mref1).in_objective
    end
end

# Test queries and used functions
@testset "Queries and Used" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = Measure(par + 2inf - x, data, [1], [1], false)
    mref = add_measure(m, meas)
    # test measure_function
    @testset "measure_function" begin
        @test isequal_canonical(measure_function(mref), par + 2inf - x)
    end
    # test measure_data
    @testset "measure_data" begin
        @test measure_data(mref) == data
    end
    # test is_analytic
    @testset "is_analytic" begin
        @test !is_analytic(mref)
    end
    # test _make_param_tuple_element (IndependentParameterIndex)
    @testset "_make_param_tuple_element (IndependentParameterIndex)" begin
        @test isequal(InfiniteOpt._make_param_tuple_element(m, IndependentParameterIndex(1), [1]), par)
    end
    # test _make_param_tuple_element (DependentParameterIndex)
    @testset "_make_param_tuple_element (DependentParameterIndex)" begin
        @test isequal(InfiniteOpt._make_param_tuple_element(m, DependentParametersIndex(1), [2, 3]), pars)
        @test isequal(InfiniteOpt._make_param_tuple_element(m, DependentParametersIndex(1), [2]), pars[1])
    end
    # test parameter_refs for measure
    @testset "parameter_refs (measure)" begin
        @test isequal(parameter_refs(mref), (par,))
    end
    # test raw_parameter_refs for measure
    @testset "raw_parameter_refs (measure)" begin
        @test isequal(raw_parameter_refs(mref), InfiniteOpt.Collections.VectorTuple(par))
    end
    # test parameter_list for measure
    @testset "parameter_list (measure)" begin
        @test isequal(parameter_list(mref), [par])
    end
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(mref)
        # prepare use case
        push!(InfiniteOpt._constraint_dependencies(mref), InfOptConstraintIndex(1))
        # test used
        @test used_by_constraint(mref)
        # undo changes
        popfirst!(InfiniteOpt._constraint_dependencies(mref))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(mref)
        # prepare use case
        push!(InfiniteOpt._measure_dependencies(mref), MeasureIndex(1))
        # test used
        @test used_by_measure(mref)
        # undo changes
        popfirst!(InfiniteOpt._measure_dependencies(mref))
    end
    # used_by_derivative
    @testset "used_by_derivative" begin
        # test not used
        @test !used_by_derivative(mref)
        # prepare use case
        push!(InfiniteOpt._derivative_dependencies(mref), DerivativeIndex(1))
        # test used
        @test used_by_derivative(mref)
        # undo changes
        popfirst!(InfiniteOpt._derivative_dependencies(mref))
    end
    # used_by_objective
    @testset "used_by_objective" begin
        # test not used
        @test !used_by_objective(mref)
        # prepare use case
        InfiniteOpt._data_object(mref).in_objective = true
        # test used
        @test used_by_objective(mref)
        # undo changes
        InfiniteOpt._data_object(mref).in_objective = false
    end
    # is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(mref)
        # prepare use case
        push!(InfiniteOpt._data_object(mref).constraint_indices, InfOptConstraintIndex(1))
        # test used
        @test is_used(mref)
        # undo changes
        InfiniteOpt._data_object(mref).constraint_indices = InfOptConstraintIndex[]
        # prepare use case
        push!(InfiniteOpt._data_object(mref).measure_indices, MeasureIndex(1))
        # test used
        @test is_used(mref)
        # undo changes
        InfiniteOpt._data_object(mref).measure_indices = MeasureIndex[]
        # prepare use case
        InfiniteOpt._data_object(mref).in_objective = true
        # test used
        @test is_used(mref)
        # undo changes
        InfiniteOpt._data_object(mref).in_objective = false
    end
end

# Test user definition methods
@testset "User Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, pars2[1:2] ~ MvNormal(ones(2), Float64[1 0; 0 1]))
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, par2))
    @variable(m, inf3, Infinite(par2))
    @variable(m, inf4, Infinite(pars))
    @variable(m, inf5, Infinite(pars2))
    @variable(m, x)
    coeff_func(x) = 1
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1], is_expect = true)
    data2 = DiscreteMeasureData(par2, [1], [1])
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data4 = FunctionalDiscreteMeasureData(par, coeff_func, 10, UniformGrid)
    data5 = FunctionalDiscreteMeasureData(pars2, coeff_func, 10, WeightedSample)
    data6 = DiscreteMeasureData(par, [1], [1], lower_bound = 0, upper_bound = 1)
    # test measure
    @testset "measure" begin
        # test scalar IndependentParameter
        mref = GeneralVariableRef(m, 1, MeasureIndex)
        @test isequal(measure(inf + x, data), mref)
        @test isequal_canonical(measure_function(mref), inf + x)
        @test supports(par) == [1]
        @test !used_by_objective(mref)
        # test nested use
        mref2 = GeneralVariableRef(m, 2, MeasureIndex)
        mref3 = GeneralVariableRef(m, 3, MeasureIndex)
        @test isequal(measure(inf + measure(inf2 + x, data2), data), mref3)
        @test isequal_canonical(measure_function(mref3), inf + mref2)
        @test isequal_canonical(measure_function(mref2), inf2 + x)
        @test supports(par) == [1]
        @test supports(par2) == [1]
        # test DependentParameters
        mref4 = GeneralVariableRef(m, 4, MeasureIndex)
        @test isequal(measure(inf4 + x, data3), mref4)
        @test isequal_canonical(measure_function(mref4), inf4 + x)
        @test supports(pars[1]) == [1]
        @test supports(pars[2]) == [1]
        # test scalar FunctionalDiscreteMeasureData
        mref6 = GeneralVariableRef(m, 5, MeasureIndex)
        @test isequal(measure(inf, data4), mref6)
        @test InfiniteOpt._core_variable_object(mref6).data.label == UniformGrid
        # test multidim FunctionalDiscreteMeasureData
        mref7 = GeneralVariableRef(m, 6, MeasureIndex)
        @test isequal(measure(inf5, data5), mref7)
        @test InfiniteOpt._core_variable_object(mref7).data.label == WeightedSample
        # test analytic evaluation
        @test measure(x, data) isa GeneralVariableRef
        @test is_analytic(measure(x, data))
        @test measure(par2, data) isa GeneralVariableRef
        @test is_analytic(measure(par2, data))
        @test measure(inf4 + measure(inf + x, data3), data) isa GeneralVariableRef
        @test !is_analytic(measure(inf4 + measure(inf + x, data3), data))
        @test is_analytic(measure(x, data6))
        @test is_analytic(measure(par2, data6))
        # test no variables
        @test_throws ErrorException measure(zero(GenericAffExpr{Float64,
                                                     GeneralVariableRef}), data)
    end
end

@testset "Macros" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, pars2[1:2] ~ MvNormal(ones(2), Float64[1 0; 0 1]))
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(pars))
    @variable(m, inf3, Infinite(pars2))
    @variable(m, x)
    coeff_func(x) = 1
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data3 = FunctionalDiscreteMeasureData(par, coeff_func, 10, UniformGrid)
    data4 = FunctionalDiscreteMeasureData(pars2, coeff_func, 10, WeightedSample)
    # test @measure
    @testset "@measure" begin
        @test_macro_throws ErrorException @measure(par, data1, data2)
        @test_macro_throws ErrorException @measure(par, data1, bob = 53)
        # TODO improve these
        @test @measure(par, data) isa GeneralVariableRef
        @test @measure(inf + par, data) isa GeneralVariableRef
        @test @measure(inf2, data2) isa GeneralVariableRef
        @test @measure(inf, data3) isa GeneralVariableRef
        @test @measure(inf3, data4) isa GeneralVariableRef
    end
end

@testset "Naming methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    data = DiscreteMeasureData(par, [1], [1])
    mref = @measure(par, data, name = "test")
    # test name
    @testset "JuMP.name" begin
        @test name(mref) == "test"
    end
    # test set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(mref, "new"), Nothing)
        @test name(mref) == "new"
    end
end

@testset "General Queries" begin
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    data = DiscreteMeasureData(par, [1], [1])
    mref1 = @measure(inf + par, data)
    mref2 = @measure(2 * inf, data)
    # test num_measures
    @testset "num_measures" begin
        @test num_measures(m) == 2
    end
    #test all_measures
    @testset "all_measures" begin
        @test isequal(all_measures(m), [mref1, mref2])
    end
end

@testset "Measure Deletion" begin
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    data = DiscreteMeasureData(par, [1], [1])
    mref1 = @measure(inf + par, data)
    mref2 = @measure(mref1, data)
    mref3 = @measure(mref1 + inf, data)
    con = build_constraint(error, mref1, MOI.GreaterThan(1.0))
    constr1 = add_constraint(m, con)
    constr2 = @constraint(m, 2 * mref1 <= 10)
    obj = @objective(m, Min, mref1)


    m2 = InfiniteModel()
    @infinite_parameter(m2, x in [0, 1])
    @variable(m2, z >= 0)
    data = FunctionalDiscreteMeasureData(x, ones, 0, All, 
               generative_support_info = UniformGenerativeInfo([0.5], All))
    mref4 = @measure(3*x, data)
    obj2 = @objective(m2, Min, z + mref4)

    # test JuMP.delete
    @testset "JuMP.delete" begin
        @test_throws AssertionError delete(m2, mref1)
        @test delete(m, mref1) isa Nothing
        @test isequal(measure_function(mref2), zero(GenericAffExpr{Float64, GeneralVariableRef}))
        @test isequal(measure_function(mref3), 1. * inf)
        @test isequal(constraint_object(constr1).func, zero(GenericAffExpr{Float64, GeneralVariableRef}))
        @test isequal(constraint_object(constr2).func, zero(GenericAffExpr{Float64, GeneralVariableRef}))
        @test isequal(objective_function(m), zero(GenericAffExpr{Float64, GeneralVariableRef}))
        @test delete(m2, mref4) isa Nothing
        @test isequal(objective_function(m2), 1. * z)
    end
end
