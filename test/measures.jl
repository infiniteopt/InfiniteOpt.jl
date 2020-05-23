# Test the core measure data structures
#=
@testset "Core Datatypes" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    # test AbstractMeasureData
    @testset "AbstractMeasureData" begin
        @test AbstractMeasureData isa DataType
    end
    # test DiscreteMeasureData
    @testset "DiscreteMeasureData" begin
        # type and normal usage
        @test DiscreteMeasureData <: AbstractMeasureData
        @test DiscreteMeasureData(par, [1], [1], "", default_weight).name == ""
        @test DiscreteMeasureData(par, [1], [1], "",
                                  default_weight).parameter_ref == par
    end
    # test MultiDiscreteMeasureData
    @testset "MultiDiscreteMeasureData" begin
        # test normal usage and type
        @test MultiDiscreteMeasureData <: AbstractMeasureData
        @test MultiDiscreteMeasureData(pars, [1], zeros(2, 1), "",
                                       default_weight).name == ""
        @test MultiDiscreteMeasureData(pars, [1], zeros(2, 1), "",
                                       default_weight).parameter_refs == pars
    end
    # test Measure
    @testset "Measure" begin
        # prepare data
        data1 = DiscreteMeasureData(par, [1], [1], "", default_weight)
        data2 =  MultiDiscreteMeasureData(pars, [1], zeros(2, 1), "", default_weight)
        # run tests
        @test Measure isa UnionAll
        @test isa(Measure(par, data1), Measure)
        @test Measure(par, data1).func == par
        @test isa(Measure(pars[1] + pars[2], data2),
                  Measure{GenericAffExpr{Float64, ParameterRef},
                          MultiDiscreteMeasureData})
    end
    # test Base.copy for MeasureRefs
    @testset "Base.copy (MeasureRef)" begin
        m2 = InfiniteModel()
        mref = MeasureRef(m, 1)
        @test copy(mref, m2) == MeasureRef(m2, 1)
    end
    # test Base.copy for DiscreteMeasureData
    @testset "Base.copy (DiscreteMeasureData)" begin
        data = DiscreteMeasureData(par, [1], [1], "", default_weight)
        @test copy(data) isa DiscreteMeasureData
    end
    # test Base.copy for MultiDiscreteMeasureData
    @testset "Base.copy (MultiDiscreteMeasureData)" begin
        data = MultiDiscreteMeasureData(pars, [1], zeros(2, 1), "", default_weight)
        @test copy(data) isa MultiDiscreteMeasureData
    end
end
=#

# Test measure data constructors
@testset "Measure Data Constructors" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1]) # Vector
    @infinite_parameter(m, pars2[1:2, 1:3] in [0, 1]) # Matrix
    @infinite_parameter(m, pars3[2:3] in [0, 1]) # DenseAxisArray
    @infinite_parameter(m, pars4[i = 1:2, j = 1:2; i <= j] in [0, 1]) # SparseAxisArray
    @finite_parameter(m, fpar, 42)
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
        @test isa(FunctionalDiscreteMeasureData(par, coeff_func, 10, :label), FunctionalDiscreteMeasureData)
        @test_throws ErrorException FunctionalDiscreteMeasureData(pars[1], coeff_func, 10, :label)
        @test_throws ErrorException FunctionalDiscreteMeasureData(par, coeff_func, -1, :label)
        @test_throws ErrorException FunctionalDiscreteMeasureData(par, coeff_func, 10, MCSample, lower_bound = -1, upper_bound = 2)
    end
    # test _check_multidim_params
    @testset "_check_params" begin
        @test_throws ErrorException InfiniteOpt._check_params(fpar)
        @test_throws ErrorException InfiniteOpt._check_params([par, pars[1]])
        @test_throws ErrorException InfiniteOpt._check_params([pars[1], pars2[1]])
        @test_throws ErrorException InfiniteOpt._check_params([pars4[1, 1], pars4[1, 2]])
        @test InfiniteOpt._check_params(pars) isa Nothing
    end
    # test _convert_param_refs_and_supports (Vector)
    @testset "_convert_param_refs_and_supports (Vector)" begin
        supp = [[0, 0], [1, 1]]
        expected = (pars, [0 1; 0 1])
        @test InfiniteOpt._convert_param_refs_and_supports(pars, supp) == expected
    end
    # test _convert_param_refs_and_supports (Array)
    @testset "_convert_param_refs_and_supports (Array)" begin
        supp = [zeros(2, 3), ones(2, 3)]
        expected = ([pars2...], [zeros(6, 1) ones(6, 1)])
        @test InfiniteOpt._convert_param_refs_and_supports(pars2, supp) == expected
    end
    # test _convert_param_refs_and_supports (DenseAxisArray)
    @testset "_convert_param_refs_and_supports (DenseAxisArray)" begin
        supp = [JuMPC.DenseAxisArray([0, 0], axes(pars3)),
                JuMPC.DenseAxisArray([1, 1], axes(pars3))]
        expected = ([pars3[2], pars3[3]], [0 1; 0 1])
        @test InfiniteOpt._convert_param_refs_and_supports(pars3, supp) == expected
    end
    # test _convert_param_refs_and_supports (SparseAxisArray)
    @testset "_convert_param_refs_and_supports (SparseAxisArray)" begin
        supp = [JuMPC.SparseAxisArray(Dict(k => 0 for k in keys(pars4.data))),
                JuMPC.SparseAxisArray(Dict(k => 1 for k in keys(pars4.data)))]
        expected = ([pars4[1, 1], pars4[1, 2], pars4[2, 2]], [0 1; 0 1; 0 1])
        @test isequal(InfiniteOpt._convert_param_refs_and_supports(pars4, supp), expected)
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
    end
    # test _check_bounds_in_set (IndependentParameter)
    @testset "_check_bounds_in_set (independent parameters)" begin
        disp_par = dispatch_variable_ref(par)
        @test_throws ErrorException InfiniteOpt._check_bounds_in_set([disp_par, disp_par], [-1, -1], [1, 1])
        @test_throws ErrorException InfiniteOpt._check_bounds_in_set([disp_par, disp_par], [0, 0], [2, 2])
        @test InfiniteOpt._check_bounds_in_set([disp_par, disp_par], [0.3, 0.2], [0.7, 0.8]) isa Nothing
    end
    # test _check_bounds_in_set (DependentParameter)
    @testset "_check_bounds_in_set (dependent parameters)" begin
        disp_vars = dispatch_variable_ref.(pars)
        @test_throws ErrorException InfiniteOpt._check_bounds_in_set(disp_vars, [-1, -1], [1, 1])
        @test_throws ErrorException InfiniteOpt._check_bounds_in_set(disp_vars, [0, 0], [2, 2])
        @test InfiniteOpt._check_bounds_in_set(disp_vars, [0.3, 0.2], [0.7, 0.8]) isa Nothing
    end
    # test multidim FunctionalDiscreteMeasureData constructor
    @testset "FunctionalDiscreteMeasureData (multidim)" begin
        @test isa(FunctionalDiscreteMeasureData(pars, coeff_func, 10, :label), FunctionalDiscreteMeasureData)
        @test_throws ErrorException FunctionalDiscreteMeasureData(pars, coeff_func, -1, :label)
    end
end
# Test data access methods
@testset "Data Queries" begin
    # initialize model and references
    m = InfiniteModel()
    f(x) = 1
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data3 = FunctionalDiscreteMeasureData(par, f, 10, MCSample)
    data4 = FunctionalDiscreteMeasureData(pars, f, 10, MCSample)
    # parameter_refs (DiscreteMeasureData)
    @testset "parameter_refs (Single)" begin
        @test parameter_refs(data) == par
    end
    # parameter_refs (MultiDiscreteMeasureData)
    @testset "parameter_refs (Multi)" begin
        @test parameter_refs(data2) == pars
    end
    # parameter_refs (Fallback)
    @testset "parameter_refs (Fallback)" begin
        @test_throws ErrorException parameter_refs(BadData())
    end
    # support_label
    @testset "support_label" begin
        @test isa(support_label(data), Symbol)
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
    # _is_expect
    @testset "_is_expect" begin
        @test !InfiniteOpt._is_expect(data)
        @test !InfiniteOpt._is_expect(BadData())
    end
    # supports (DiscreteMeasureData)
    @testset "supports (Single)" begin
        @test supports(data) == Float64[1]
    end
    # supports (MultiDiscreteMeasureData)
    @testset "supports (Multi)" begin
        @test supports(data2) == ones(Float64, 2, 1)
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
# Test add_measure and helper functions
@testset "Measure Addition" begin
    m = InfiniteModel()
    coeff_func(x) = 1
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= par3 <= 1)
    @infinite_parameter(m, 0 <= pars2[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @hold_variable(m, x)
    data1 = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data3 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample)
    data4 = FunctionalDiscreteMeasureData(pars, coeff_func, 5, MCSample)
#    meas = Measure(par + 2inf - x, data, [1], [1], false)
    # test _add_supports_to_multiple_parameters (independent)
    @testset "_add_supports_to_multiple_parameters (independent)" begin
        prefs = dispatch_variable_ref.([par, par])
        @test isa(InfiniteOpt._add_supports_to_multiple_parameters(prefs, [0.1 0.2;0.3 0.4], :label), Nothing)
        delete_supports(par)
    end
    # test _add_supports_to_multiple_parameters (dependent)
    @testset "_add_supports_to_multiple_parameters (dependent)" begin
        prefs = dispatch_variable_ref.(pars)
        @test isa(InfiniteOpt._add_supports_to_multiple_parameters(prefs, [0. 1.; 0. 1.], :label), Nothing)
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
        @test num_supports(pars) == 5 # wrong with CollectionSet support generation
        # clear supports
        delete_supports(pars)
    end
    # test add_supports_to_parameters (fallbacks)
    @testset "add_supports_to_parameters (fallbacks)" begin
        @test_throws ErrorException add_supports_to_parameters(BadData())
    end
    # test add_measure
    @testset "add_measure" begin
        meas1 = Measure(par + 2inf - x, data1, [1], [1], false)
        mref1 = MeasureRef(m, MeasureIndex(1))
        @test dispatch_variable_ref(add_measure(m, meas1, "measure1")) == mref1
        @test supports(par) == [1.0]
        @test InfiniteOpt._measure_dependencies(par) == [MeasureIndex(1)]
        @test InfiniteOpt._data_object(mref1).name == "measure1"
        @test !InfiniteOpt._data_object(mref1).in_objective
    end
end

# Test measure construction methods
@testset "Measure Construction" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars))
    @hold_variable(m, x)
    coeff_func(x) = 1

    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(par2, [1], [1])
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data4 = FunctionalDiscreteMeasureData(par, coeff_func, 5, MCSample)
    data5 = FunctionalDiscreteMeasureData(pars, coeff_func, 5, MCSample)

    # test measure_data_in_hold_bounds (DiscreteMeasureData)
    @testset "measure_data_in_hold_bounds (Discrete)" begin
        # test empty bounds
        @test measure_data_in_hold_bounds(data, ParameterBounds())
        # test in domain
        @set_parameter_bounds(x, par == 1)
        @test measure_data_in_hold_bounds(data, parameter_bounds(x))
        # test outside domain
        @add_parameter_bounds(x, par == 0)
        @test !measure_data_in_hold_bounds(data, parameter_bounds(x))
        set_parameter_bounds(x, ParameterBounds(), force = true)
        delete_supports(par)
    end
    # test measure_data_in_hold_bounds (MultiDiscreteMeasureData)
    @testset "measure_data_in_hold_bounds (Multi)" begin
        # test empty bounds
        @test measure_data_in_hold_bounds(data3, ParameterBounds())
        # test in domain
        @set_parameter_bounds(x, pars == 1)
        @test measure_data_in_hold_bounds(data3, parameter_bounds(x))
        # test not in the domain
        @add_parameter_bounds(x, pars == 0) # TODO: this macro seems to have problems because of support deletion of single dependent parameter
        @test !measure_data_in_hold_bounds(data3, parameter_bounds(x))
        set_parameter_bounds(x, ParameterBounds(), force = true)
        delete_supports(pars)
    end
    # test measure_data_in_hold_bounds (Fallback)
    @testset "measure_data_in_hold_bounds (Fallback)" begin
        warn = "Unable to check if hold variables bounds are valid in measure " *
               "with measure data type `BadData`. This can be resolved by " *
               "extending `measure_data_in_hold_bounds`."
        @test_logs (:warn, warn) measure_data_in_hold_bounds(BadData(), ParameterBounds())
    end
    # test _check_var_bounds (GeneralVariableRef)
    @testset "_check_var_bounds (General)" begin
        @test isa(InfiniteOpt._check_var_bounds(inf, data), Nothing)
    end

    # test _check_var_bounds (HoldVariableRef with DiscreteMeasureData)
    @testset "_check_var_bounds (Hold with Discrete)" begin
        # test normal
        @set_parameter_bounds(x, par == 1)
        @test isa(InfiniteOpt._check_var_bounds(x, data), Nothing)
        # test error
        @add_parameter_bounds(x, par == 0)
        @test_throws ErrorException InfiniteOpt._check_var_bounds(x, data)
        set_parameter_bounds(x, ParameterBounds(), force = true)
        delete_supports(par)
    end
    # test _check_var_bounds (HoldVariableRef with MultiDiscreteMeasureData)
    @testset "_check_var_bounds (Hold with Multi)" begin
        # test normal
        @set_parameter_bounds(x, pars == 1)
        @test isa(InfiniteOpt._check_var_bounds(x, data3), Nothing)
        # test error
        @add_parameter_bounds(x, pars == 0)
        @test_throws ErrorException InfiniteOpt._check_var_bounds(x, data3)
        set_parameter_bounds(x, ParameterBounds(), force = true)
        delete_supports(pars)
    end
    # test _check_var_bounds (MeasureRef)
    # TODO: maybe push this after Measure is tested
    @testset "_check_var_bounds (Measure)" begin
        # make some measures
        meas1 = Measure(par, data, [1], [1], false)
        mref1 = add_measure(m, meas1)
        meas2 = Measure(x - mref1, data, [1], [1], false)
        mref2 = add_measure(m, meas2)
        # test normal
        @set_parameter_bounds(x, par == 1)
        @test isa(InfiniteOpt._check_var_bounds(mref1, data), Nothing)
        @test isa(InfiniteOpt._check_var_bounds(mref2, data), Nothing)
        set_parameter_bounds(x, ParameterBounds(), force = true)
    end
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

# Test queries and used functions
@testset "Queries and Used" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @hold_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = Measure(par + 2inf - x, data, [1], [1], false)
    mref = add_measure(m, meas)
    # test measure_function
    @testset "measure_function" begin
        @test measure_function(mref) == par + 2inf - x
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
        @test InfiniteOpt._make_param_tuple_element(m, IndependentParameterIndex(1), [1]) == par
    end
    # test _make_param_tuple_element (DependentParameterIndex)
    @testset "_make_param_tuple_element (DependentParameterIndex)" begin
        @test InfiniteOpt._make_param_tuple_element(m, DependentParametersIndex(1), [2, 3]) == pars
        @test InfiniteOpt._make_param_tuple_element(m, DependentParametersIndex(1), [2]) == pars[1]
    end
    # test parameter_refs for measure
    @testset "parameter_refs (measure)" begin
        @test collect(parameter_refs(mref)) == [par]
    end
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(mref)
        # prepare use case
        push!(InfiniteOpt._data_object(mref).constraint_indices, ConstraintIndex(1))
        # test used
        @test used_by_constraint(mref)
        # undo changes
        InfiniteOpt._data_object(mref).constraint_indices = ConstraintIndex[]
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(mref)
        # prepare use case
        push!(InfiniteOpt._data_object(mref).measure_indices, MeasureIndex(1))
        # test used
        @test used_by_measure(mref)
        # undo changes
        InfiniteOpt._data_object(mref).measure_indices = MeasureIndex[]
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
        push!(InfiniteOpt._data_object(mref).constraint_indices, ConstraintIndex(1))
        # test used
        @test is_used(mref)
        # undo changes
        InfiniteOpt._data_object(mref).constraint_indices = ConstraintIndex[]
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
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_parameter(m, pars2[1:2] in MvNormal(ones(2), Float64[1 0; 0 1]))
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars))
    @infinite_variable(m, inf5(pars2))
    @hold_variable(m, x)
    coeff_func(x) = 1
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(par2, [1], [1])
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]])
    data4 = FunctionalDiscreteMeasureData(par, coeff_func, 10, UniformGrid)
    data5 = FunctionalDiscreteMeasureData(pars2, coeff_func, 10, WeightedSample)
    # test measure
    @testset "measure" begin
        # test scalar IndependentParameter
        mref = MeasureRef(m, MeasureIndex(1))
        @test dispatch_variable_ref(measure(inf + x, data)) == mref
#            @test name(mref) == "a(inf(par) + x)"
        @test supports(par) == [1]
        @test !used_by_objective(mref)
        # test nested use
        mref2 = MeasureRef(m, MeasureIndex(2))
        mref3 = MeasureRef(m, MeasureIndex(3))
        @test dispatch_variable_ref(measure(inf + measure(inf2 + x, data2), data)) == mref3
#            @test name(mref3) == "a(inf(par) + b(inf2(par, par2) + x))"
        @test supports(par) == [1]
        @test supports(par2) == [1]
        # test DependentParameters
        mref4 = MeasureRef(m, MeasureIndex(4))
        @test dispatch_variable_ref(measure(inf4 + x, data3)) == mref4
#            @test name(mref4) == "c(inf4(pars) + x)"
        @test supports(pars[1]) == [1]
        @test supports(pars[2]) == [1]
        # test with hold bounds
        @set_parameter_bounds(x, par == 1)
        mref5 = MeasureRef(m, MeasureIndex(5))
        @test dispatch_variable_ref(measure(inf + x, data)) == mref5
#            @test name(mref) == "a(inf(par) + x)"
        @test supports(par) == [1]
        # test scalar FunctionalDiscreteMeasureData
        mref6 = MeasureRef(m, MeasureIndex(6))
        @test dispatch_variable_ref(measure(inf, data4)) == mref6
        @test InfiniteOpt._core_variable_object(mref6).data.label == UniformGrid
        # test multidim FunctionalDiscreteMeasureData
        mref7 = MeasureRef(m, MeasureIndex(7))
        @test dispatch_variable_ref(measure(inf5, data5)) == mref7
        @test InfiniteOpt._core_variable_object(mref7).data.label == WeightedSample
        # test analytic evaluation
        @test dispatch_variable_ref(measure(x, data)) isa MeasureRef
        @test is_analytic(measure(x, data))
        @test dispatch_variable_ref(measure(par2, data)) isa MeasureRef
        @test is_analytic(measure(par2, data))
        @test dispatch_variable_ref(measure(inf4 + measure(inf + x, data3), data)) isa MeasureRef
        @test !is_analytic(measure(inf4 + measure(inf + x, data3), data))
        # test no variables
        @test_throws ErrorException measure(zero(GenericAffExpr{Float64,
                                                     GeneralVariableRef}), data)
    end
end

@testset "Macros" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_parameter(m, pars2[1:2] in MvNormal(ones(2), Float64[1 0; 0 1]))
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(pars))
    @infinite_variable(m, inf3(pars2))
    @hold_variable(m, x)
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
        @test dispatch_variable_ref(@measure(par, data)) isa MeasureRef
        @test dispatch_variable_ref(@measure(inf + par, data)) isa MeasureRef
        @test dispatch_variable_ref(@measure(inf2, data2)) isa MeasureRef
        @test dispatch_variable_ref(@measure(inf, data3)) isa MeasureRef
        @test dispatch_variable_ref(@measure(inf3, data4)) isa MeasureRef
    end
end

@testset "Naming methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
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
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    data = DiscreteMeasureData(par, [1], [1])
    mref1 = @measure(inf + par, data)
    mref2 = @measure(2 * inf, data)
    # test num_measures
    @testset "num_measures" begin
        @test num_measures(m) == 2
    end
    #test all_measures
    @testset "all_measures" begin
        @test all_measures(m) == [mref1, mref2]
    end
end

@testset "Measure Deletion" begin
    # TODO finish this
end
#=
# Test user definition methods
@testset "User Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars))
    @hold_variable(m, x)
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1], name = "a")
    data2 = DiscreteMeasureData(par2, [1], [1], name = "b")
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]], name = "c")
    # test _model_from_expr
    @testset "_model_from_expr" begin
        @test InfiniteOpt._model_from_expr([par]) == m
        @test InfiniteOpt._model_from_expr([inf, par]) == m
        @test InfiniteOpt._model_from_expr([inf]) == m
        @test isa(InfiniteOpt._model_from_expr(GeneralVariableRef[]), Nothing)
    end
    # test measure
    @testset "measure" begin
        # test single use
        index = m.next_meas_index
        mref = MeasureRef(m, index + 1)
        @test measure(inf + x, data) == mref
        @test name(mref) == "a(inf(par) + x)"
        @test supports(par) == [1]
        @test !m.meas_in_objective[index + 1]
        # test nested use
        mref2 = MeasureRef(m, index + 2)
        mref3 = MeasureRef(m, index + 3)
        @test measure(inf + measure(inf2 + x, data2), data) == mref3
        @test name(mref3) == "a(inf(par) + b(inf2(par, par2) + x))"
        @test supports(par) == [1]
        @test supports(par2) == [1]
        # test vector parameter
        mref4 = MeasureRef(m, index + 4)
        @test measure(inf4 + x, data3) == mref4
        @test name(mref4) == "c(inf4(pars) + x)"
        @test supports(pars[1]) == [1]
        @test supports(pars[2]) == [1]
        # test with hold bounds
        @set_parameter_bounds(x, par == 1)
        mref = MeasureRef(m, index + 5)
        @test measure(inf + x, data) == mref
        @test name(mref) == "a(inf(par) + x)"
        @test supports(par) == [1]
        # test no variables
        @test_throws ErrorException measure(zero(GenericAffExpr{Float64,
                                                     GeneralVariableRef}), data)
        # test not dependent
        @test measure(x, data) isa MeasureRef
        @test measure(par2, data) isa MeasureRef
        @test measure(inf4 + measure(inf + x, data3), data) isa MeasureRef
        # test with bad variable bounds
        InfiniteOpt._update_variable_param_bounds(x, ParameterBounds(Dict(par => IntervalSet(0, 0))))
        @test_throws ErrorException measure(inf + x, data)
        # test invalid parameter
        delete!(m.params, JuMP.index(par))
        @test_throws ErrorException measure(inf, data)
    end
end

@testset "User Definition w/o Measure Data" begin
    m = InfiniteModel()
    dist1 = Normal(0., 1.)
    dist2 = MvNormal([0., 0.], [1. 0.;0. 10.])
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= par3 <= 1)
    @infinite_parameter(m, 0 <= pars1[1:2] <= 1, container = SparseAxisArray)
    @infinite_parameter(m, 0 <= pars2[("a", "b")] <= 1)
    @infinite_parameter(m, 0 <= pars3[1:2] <= 1)
    @infinite_parameter(m, rp in dist1)
    @infinite_parameter(m, rp2[1:2] in dist2)
    @infinite_parameter(m, 0 <= par4 <= 1)
    @infinite_parameter(m, 0 <= par5 <= 1)
    @infinite_parameter(m, par6 in [0, Inf])
    @infinite_parameter(m, par7 in [-Inf, Inf])
    sets = [IntervalSet(0,1), DistributionSet(Uniform())]
    @infinite_parameter(m, pars4[i = 1:2], set = sets[i])
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par3))
    @infinite_variable(m, inf4(pars1))
    @infinite_variable(m, inf5(pars1, pars2))
    @infinite_variable(m, inf6(pars2))
    @infinite_variable(m, inf7(pars3))
    @infinite_variable(m, inf8(rp))
    @infinite_variable(m, inf9(rp2))
    @infinite_variable(m, inf10(rp, rp2))
    @infinite_variable(m, inf11(par4))
    @infinite_variable(m, inf12(par5))
    @infinite_variable(m, inf13(par6))
    @infinite_variable(m, inf14(par7))
    @hold_variable(m, x)

    # test integral
    @testset "integral" begin
        meas1 = integral(inf, num_supports = 5, eval_method = gauss_legendre)
        (expected_supps, expected_coeffs) = FGQ.gausslegendre(5)
        expected_supps = expected_supps .* 0.5 .+ 0.5
        expected_coeffs = expected_coeffs .* 0.5
        @test all(measure_data(meas1).supports .== expected_supps)
        @test all(measure_data(meas1).coefficients .== expected_coeffs)

        add_supports(par2, [0.3, 0.7])
        meas2 = integral(inf2, par2, use_existing_supports = true)
        @test measure_data(meas2).parameter_ref == par2
        @test measure_data(meas2).supports == [0.3, 0.7]
        meas2 = integral(inf2, par2, 0.5, 0.9, use_existing_supports = true)
        @test measure_data(meas2).supports == [0.7]
        meas2 = integral(inf2, [par2], [0.5], [0.9], use_existing_supports = true)
        @test measure_data(meas2).supports == [0.7]

        meas3 = integral(inf4, num_supports = 5, eval_method = sampling)
        @test pars1[1] in measure_data(meas3).parameter_refs
        @test pars1[2] in measure_data(meas3).parameter_refs

        meas4 = integral(inf5, pars2, num_supports = 5)
        @test pars2["a"] in measure_data(meas4).parameter_refs
        @test pars2["b"] in measure_data(meas4).parameter_refs

        meas5 = integral(inf6, use_existing_supports = true)
        @test measure_data(meas4).supports == measure_data(meas5).supports

        add_supports(pars3[1], [0.3, 0.7])
        add_supports(pars3[2], [0.3, 0.7])
        meas6 = integral(inf7, pars3, 0.5, 1.0, use_existing_supports = true)
        @test measure_data(meas6).supports == ones(Float64, 2, 1) * 0.7
        meas6 = integral(inf7, pars3, [0.5, 0.5], [1.0, 1.0], use_existing_supports = true)
        @test measure_data(meas6).supports == ones(Float64, 2, 1) * 0.7

        meas7 = integral(inf11, par4, num_supports = 5, eval_method = quadrature)
        (expected_supps, expected_coeffs) = FGQ.gausslegendre(5)
        expected_supps = expected_supps .* 0.5 .+ 0.5
        expected_coeffs = expected_coeffs .* 0.5
        @test all(measure_data(meas7).supports .== expected_supps)
        @test all(measure_data(meas7).coefficients .== expected_coeffs)

        add_supports(par5, [0.3, 0.7])
        warn = "Quadrature method will not be used because use_existing_supports is set as true."
        @test_logs (:warn, warn) integral(inf12, par5,
                               use_existing_supports = true, eval_method = quadrature)

        meas8 = integral(inf13, par6, num_supports = 5, eval_method = quadrature)
        (expected_supps, expected_coeffs) = FGQ.gausslaguerre(5)
        expected_coeffs = expected_coeffs .* exp.(expected_supps)
        @test all(measure_data(meas8).supports .== expected_supps)
        @test all(measure_data(meas8).coefficients .== expected_coeffs)

        meas9 = integral(inf14, par7, num_supports = 5, eval_method = quadrature)
        (expected_supps, expected_coeffs) = FGQ.gausshermite(5)
        expected_coeffs = expected_coeffs .* exp.(expected_supps.^2)
        @test all(measure_data(meas9).supports .== expected_supps)
        @test all(measure_data(meas9).coefficients .== expected_coeffs)

        # test errors
        @test_throws ErrorException integral(x)
        @test_throws ErrorException integral(inf, ParameterRef[])
        @test_throws ErrorException integral(inf2)
        @test_throws ErrorException integral(inf2, [par, par2])
        @test_throws ErrorException integral(inf2, par, 1., 3.)
        @test_throws ErrorException integral(inf2, par, [0., 1.])
        @test_throws ErrorException integral(inf2, par, 0., [1., 1.])
        @test_throws ErrorException integral(inf2, par, 0.5, 0.)
        @test_throws ErrorException integral(inf8, use_existing_supports = true)
        @test_throws ErrorException integral(meas1)
        @test_throws ErrorException integral(inf4, pars1, eval_method = quadrature)
        @test_throws ErrorException integral(inf, par, -1)
        @test_throws ErrorException integral(inf, par, 0, 2)
        @test_throws ErrorException integral(inf4, pars1, -1)
        @test_throws ErrorException integral(inf4, pars1, 0, 2)
        @test_throws ErrorException integral(x, pars4)
    end
    # test support_sum
    @testset "support_sum" begin
        sum1 = support_sum(inf2, par2)
        @test measure_data(sum1).parameter_ref == par2
        @test measure_data(sum1).supports == [0.3, 0.7]
        @test name(sum1) == "sum(inf2(par, par2))"
        sum2 = support_sum(inf7)
        @test pars3[1] in measure_data(sum2).parameter_refs
        @test pars3[2] in measure_data(sum2).parameter_refs
        @test measure_data(sum2).supports == Float64[0.3 0.7; 0.3 0.7]
    end
    # test expectation measure
    @testset "expect" begin
        expect1 = expect(inf8, num_supports = 5)
        expect2 = expect(inf9, num_supports = 5)
        @test_throws ErrorException expect(inf10)
        expect3 = expect(inf10, rp, use_existing_supports = true)
        check1 = expect(inf8, use_existing_supports = true)
        check2 = expect(inf9, use_existing_supports = true)
        @test measure_data(expect1).supports == measure_data(check1).supports
        @test measure_data(expect2).supports == measure_data(check2).supports
        @test measure_data(expect3).supports == measure_data(expect1).supports
    end
    # test set_integral_defaults
    @testset "set_integral_defaults" begin
        set_integral_defaults(m, num_supports = 5, eval_method = quadrature,
                            new_kwarg = true)
        def_vals = integral_defaults(m)
        @test def_vals[:num_supports] == 5
        @test def_vals[:eval_method] == quadrature
        @test def_vals[:name] == "integral"
        @test def_vals[:weight_func] == default_weight
        @test def_vals[:use_existing_supports] == false
        @test def_vals[:new_kwarg] == true
    end
    # test measure with default keyword argument values
    @testset "default measure" begin
        delete!(integral_defaults(m), :new_kwarg)
        meas = integral(inf)
        (expected_supps, expected_coeffs) = FGQ.gausslegendre(5)
        expected_supps = expected_supps .* 0.5 .+ 0.5
        expected_coeffs = expected_coeffs .* 0.5
        @test all(measure_data(meas).supports .== expected_supps)
        @test all(measure_data(meas).coefficients .== expected_coeffs)
    end
end
=#
