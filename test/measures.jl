# Test the core measure data structures
@testset "Core Datatypes" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    # test AbstractMeasureData
    @testset "AbstractMeasureData" begin
        @test AbstractMeasureData isa DataType
    end
    # test DiscreteMeasureData
    @testset "DiscreteMeasureData" begin
        # type and normal usage
        @test DiscreteMeasureData <: AbstractMeasureData
        @test DiscreteMeasureData(par, [1], [1], "", InfiniteOpt._w).name == ""
        @test DiscreteMeasureData(par, [1], [1], "",
                                  InfiniteOpt._w).parameter_ref == par
        # test errors
        @test_throws ErrorException DiscreteMeasureData(par, Int[], [1], "",
                                                        InfiniteOpt._w)
        @test_throws ErrorException DiscreteMeasureData(par, [1], [2], "",
                                                        InfiniteOpt._w)
        @test_throws ErrorException DiscreteMeasureData(par, [1], [-1], "",
                                                        InfiniteOpt._w)
    end
    # test MultiDiscreteMeasureData
    @testset "MultiDiscreteMeasureData" begin
        # test normal usage and type
        @test MultiDiscreteMeasureData <: AbstractMeasureData
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0])
        @test MultiDiscreteMeasureData(pars, [1], [supp], "",
                                       InfiniteOpt._w).name == ""
        @test MultiDiscreteMeasureData(pars, [1], [supp], "",
                                       InfiniteOpt._w).parameter_ref == pars
        # test errors
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1],
                          JuMP.Containers.SparseAxisArray[], "", InfiniteOpt._w)
        @test_throws ErrorException MultiDiscreteMeasureData(pars, Int[],
                          JuMP.Containers.SparseAxisArray[], "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [-1, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [2, 0])
        @test_throws ErrorException MultiDiscreteMeasureData(pars, [1], [supp],
                                                             "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0, 0])
        pars2 = JuMP.Containers.SparseAxisArray(Dict((1,) => pars[1],
                                                  (2,) => pars[2], (3,) => par))
        @test_throws ErrorException MultiDiscreteMeasureData(pars2, [1], [supp],
                                                             "", InfiniteOpt._w)
    end
    # test Measure
    @testset "Measure" begin
        # prepare data
        data1 = DiscreteMeasureData(par, [1], [1], "", InfiniteOpt._w)
        supp = convert(JuMP.Containers.SparseAxisArray, [0, 0])
        data2 =  MultiDiscreteMeasureData(pars, [1], [supp], "", InfiniteOpt._w)
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
end

# Test measure data constructors
@testset "Measure Data Constructors" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @finite_parameter(m, fpar, 42)
    # test default weight function
    @testset "_w" begin
        @test InfiniteOpt._w(42) == 1
        @test InfiniteOpt._w(zeros(10)) == 1
    end
    # test DiscreteMeasureData constructor for scalar parameters
    @testset "DiscreteMeasureData (scalar)" begin
        # test normal usage
        @test isa(DiscreteMeasureData(par, [1], [1]), DiscreteMeasureData)
        @test DiscreteMeasureData(par, [1], [1]).name == "measure"
        @test DiscreteMeasureData(par, [1], [1], name = "test").name == "test"
        @test DiscreteMeasureData(par, [1], [1]).weight_function(42) == 1
        # test errors
        @test_throws ErrorException DiscreteMeasureData(par, [1], Int[])
        @test_throws ErrorException DiscreteMeasureData(par, [1], [2])
        @test_throws ErrorException DiscreteMeasureData(fpar, [1], [1])
    end
    # test MultiDiscreteMeasureData constructor for array parameters
    @testset "DiscreteMeasureData (array)" begin
        # test normal usage
        @test isa(DiscreteMeasureData(pars, [1], [[1,1]]),
                  MultiDiscreteMeasureData)
        @test DiscreteMeasureData(pars, [1], [[1, 1]]).name == "measure"
        @test DiscreteMeasureData(pars, [1], [[1, 1]], name = "t").name == "t"
        @test DiscreteMeasureData(pars, [1], [[1, 1]]).weight_function([0, 1]) == 1
        # test errors
        @test_throws ErrorException DiscreteMeasureData(pars, [1], [[1]])
        @test_throws ErrorException DiscreteMeasureData(pars, [1, 2], [[0, 0]])
        pars2 = JuMP.Containers.SparseAxisArray(Dict((1,) => pars[1],
                                                  (2,) => pars[2], (3,) => par))
        @test_throws ErrorException DiscreteMeasureData(pars2, [1], [[0, 0, 0]])
        @test_throws ErrorException DiscreteMeasureData([fpar, fpar], [1],
                                                        [[1, 1]])
    end
end

# Test basic JuMP extensions
@testset "Basic JuMP Extensions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    data = DiscreteMeasureData(par, [1], [1])
    mref = MeasureRef(m, 1)
    m.meas_to_name[1] = "test"
    m.measures[1] = Measure(par, data)
    # test name
    @testset "JuMP.name" begin
        @test name(mref) == "test"
    end
    # test set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(mref, "new"), Nothing)
        @test name(mref) == "new"
    end
    # test is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, mref)
        @test !is_valid(m, MeasureRef(InfiniteModel(), 1))
        @test !is_valid(m, MeasureRef(m, 2))
    end
end

# Test measure definition
@testset "Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    @infinite_variable(m, inf(par))
    @hold_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1], name = "test")
    data2 = DiscreteMeasureData(pars, [1], [[1, 1]])
    meas = Measure(par + 2inf - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
    # _make_meas_name
    @testset "_make_meas_name" begin
        @test InfiniteOpt._make_meas_name(meas) == "test(par + 2 inf(par) - x)"
    end
    # test _update_var_meas_mapping
    @testset "_update_var_meas_mapping" begin
        # retrieve variables and parameters
        mref = MeasureRef(m, 2)
        vars = [par, inf, x, mref, rv]
        # test initial reference
        @test isa(InfiniteOpt._update_var_meas_mapping(vars, 1), Nothing)
        @test m.var_to_meas[JuMP.index(inf)] == [1]
        @test m.var_to_meas[JuMP.index(x)] == [1]
        @test m.param_to_meas[JuMP.index(par)] == [1]
        @test m.meas_to_meas[JuMP.index(mref)] == [1]
        @test m.reduced_to_meas[JuMP.index(rv)] == [1]
        # test repeated reference
        @test isa(InfiniteOpt._update_var_meas_mapping(vars, 2), Nothing)
        @test m.var_to_meas[JuMP.index(inf)] == [1, 2]
        @test m.var_to_meas[JuMP.index(x)] == [1, 2]
        @test m.param_to_meas[JuMP.index(par)] == [1, 2]
        @test m.meas_to_meas[JuMP.index(mref)] == [1, 2]
        @test m.reduced_to_meas[JuMP.index(rv)] == [1, 2]
        # clear additions
        delete!(m.var_to_meas, JuMP.index(inf))
        delete!(m.var_to_meas, JuMP.index(x))
        delete!(m.param_to_meas, JuMP.index(par))
        delete!(m.meas_to_meas, JuMP.index(mref))
        delete!(m.reduced_to_meas, JuMP.index(rv))
    end
    # test _update_param_data_mapping (scalar)
    @testset "_update_param_data_mapping (scalar)" begin
        # test first addition
        @test isa(InfiniteOpt._update_param_data_mapping(m, data, 1), Nothing)
        @test m.param_to_meas[JuMP.index(par)] == [1]
        # test second addition
        @test isa(InfiniteOpt._update_param_data_mapping(m, data, 2), Nothing)
        @test m.param_to_meas[JuMP.index(par)] == [1, 2]
        # clear additions
        delete!(m.param_to_meas, JuMP.index(par))
    end
    # test _update_param_data_mapping (array)
    @testset "_update_param_data_mapping (array)" begin
        # test first addition
        @test isa(InfiniteOpt._update_param_data_mapping(m, data2, 1), Nothing)
        @test m.param_to_meas[JuMP.index(pars[1])] == [1]
        @test m.param_to_meas[JuMP.index(pars[2])] == [1]
        # test second addition
        @test isa(InfiniteOpt._update_param_data_mapping(m, data2, 2), Nothing)
        @test m.param_to_meas[JuMP.index(pars[1])] == [1, 2]
        @test m.param_to_meas[JuMP.index(pars[2])] == [1, 2]
        # clear additions
        delete!(m.param_to_meas, JuMP.index(pars[1]))
        delete!(m.param_to_meas, JuMP.index(pars[2]))
    end
    # test _update_param_data_mapping (fallback)
    @testset "_update_param_data_mapping (fallback)" begin
        warn = "Unable to map parameter dependence for measure data type " *
               "BadData. Parameter deletion methods should not be used."
        @test_logs (:warn, warn) InfiniteOpt._update_param_data_mapping(m,
                                                                   BadData(), 1)
    end
    # test _add_supports_to_parameters (scalar)
    @testset "_add_supports_to_parameters (scalar)" begin
        # test functionality
        @test isa(InfiniteOpt._add_supports_to_parameters(par, [1]), Nothing)
        @test supports(par) == [1]
        # clear supports
        delete_supports(par)
    end
    # test _add_supports_to_parameters (array)
    @testset "_add_supports_to_parameters (array)" begin
        # test functionality
        supp = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        @test isa(InfiniteOpt._add_supports_to_parameters(pars, [supp]), Nothing)
        @test supports(pars[1]) == [1]
        @test supports(pars[2]) == [1]
        # clear supports
        delete_supports(pars[1])
        delete_supports(pars[2])
    end
    # test add_measure
    @testset "add_measure" begin
        mref = MeasureRef(m, 1)
        @test add_measure(m, meas) == mref
        @test name(mref) == "test(par + 2 inf(par) - x)"
        @test supports(par) == [1]
        @test m.var_to_meas[JuMP.index(inf)] == [1]
        @test m.var_to_meas[JuMP.index(x)] == [1]
        @test m.param_to_meas[JuMP.index(par)] == [1]
        @test !m.meas_in_objective[JuMP.index(mref)]
        # test errors
        @test_throws VariableNotOwned add_measure(m, Measure(@variable(Model()),
                                                             data))
    end
end

# Test queries
@testset "Queries" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @hold_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = Measure(par + 2inf - x, data)
    mref = add_measure(m, meas)
    # test measure_function
    @testset "measure_function" begin
        @test measure_function(mref) == par + 2inf - x
    end
    # test measure_data
    @testset "measure_data" begin
        @test measure_data(mref) == data
    end
end

# Test methods variable search methods
@testset "Parameter Search" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars))
    @hold_variable(m, x)
    m.reduced_info[-1] = ReducedInfiniteInfo(inf2, Dict(2 => 0.5))
    m.infinite_to_reduced[JuMP.index(inf2)] = [-1]
    rv = ReducedInfiniteVariableRef(m, -1)
    # prepare measures
    data = DiscreteMeasureData(par, [1], [1])
    data2 = DiscreteMeasureData(par2, [1], [1])
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]])
    meas1 = Measure(par, data)
    mref1 = add_measure(m, meas1)
    meas2 = Measure(inf3, data2)
    mref2 = add_measure(m, meas2)
    meas3 = Measure(3 - inf2, data2)
    mref3 = add_measure(m, meas3)
    meas4 = Measure(x - mref3, data)
    mref4 = add_measure(m, meas4)
    meas5 = Measure(x + 3mref1, data)
    mref5 = add_measure(m, meas5)
    meas6 = Measure(x + inf4, data3)
    mref6 = add_measure(m, meas6)
    meas7 = Measure(rv, data)
    mref7 = add_measure(m, meas7)
    # test _has_variable
    @testset "_has_variable" begin
        # test in vector
        @test InfiniteOpt._has_variable([par], par)
        @test InfiniteOpt._has_variable([par, x, inf], inf)
        # test not in vector
        @test !InfiniteOpt._has_variable([par], inf)
        @test !InfiniteOpt._has_variable([par, x, inf], inf2)
        # test variable in prior
        @test InfiniteOpt._has_variable([x], inf, prior = [inf, par])
        @test InfiniteOpt._has_variable([x], inf, prior = [x, par, inf])
        # test variable not in prior
        @test !InfiniteOpt._has_variable([x], inf, prior = [x, par])
        # test variable in measure
        @test InfiniteOpt._has_variable([x, mref1], par)
        @test InfiniteOpt._has_variable([x, mref1, inf], par)
        @test InfiniteOpt._has_variable([x, mref2, mref1], par)
        @test InfiniteOpt._has_variable([x, mref2, mref1, x], par)
        # test variable not in measure
        @test !InfiniteOpt._has_variable([x, mref2], par)
        @test !InfiniteOpt._has_variable([x, mref2, inf], par)
        @test !InfiniteOpt._has_variable([x, mref2, mref3, x], par)
        # test deeper measure nesting
        @test InfiniteOpt._has_variable([x, mref2, mref4], inf2)
        @test !InfiniteOpt._has_variable([x, mref2, mref4], par)
    end
    # test _has_parameter
    @testset "_has_parameter" begin
        # test in vector
        @test InfiniteOpt._has_parameter([par], par)
        @test InfiniteOpt._has_parameter([x, par], par)
        @test InfiniteOpt._has_parameter([x, inf], par)
        @test InfiniteOpt._has_parameter([x, rv], par)
        # test not in vector
        @test !InfiniteOpt._has_parameter([par], par2)
        @test !InfiniteOpt._has_parameter([x, par], par2)
        @test !InfiniteOpt._has_parameter([x, inf], par2)
        # test in measure
        @test InfiniteOpt._has_parameter([x, mref1], par)
        @test InfiniteOpt._has_parameter([x, mref1, inf3], par)
        @test InfiniteOpt._has_parameter([x, mref2, mref1], par)
        @test InfiniteOpt._has_parameter([x, mref2, mref3], par)
        @test InfiniteOpt._has_parameter([x, mref2, mref4], par)
        # test not in measure
        @test !InfiniteOpt._has_parameter([x, mref1], par2)
        @test !InfiniteOpt._has_parameter([x, mref1, inf], par2)
        @test !InfiniteOpt._has_parameter([x, mref1, mref1], par2)
        @test !InfiniteOpt._has_parameter([x, mref1, mref5], par2)
    end
    # _check_has_parameter (scalar)
    @testset "_check_has_parameter (scalar)" begin
        # test that is has the parameter
        @test isa(InfiniteOpt._check_has_parameter([par], par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([inf], par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([x, inf], par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([x, mref2, mref4], par),
                  Nothing)
        # test is does not have the parameter
        @test_throws ErrorException InfiniteOpt._check_has_parameter([par], par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter([inf], par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter([x, inf],
                                                                     par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter([x, mref1,
                                                                    mref5], par2)
    end
    # _check_has_parameter (array)
    @testset "_check_has_parameter (array)" begin
        # test that is has the parameter
        @test isa(InfiniteOpt._check_has_parameter([pars[1], pars[2]], pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([inf4], pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([x, inf4], pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter([x, mref6], pars), Nothing)
        # test is does not have the parameter
        @test_throws ErrorException InfiniteOpt._check_has_parameter([pars[1]], pars)
        @test_throws ErrorException InfiniteOpt._check_has_parameter([inf2], pars)
        @test_throws ErrorException InfiniteOpt._check_has_parameter([x, mref1,
                                                                    mref5], pars)
    end
end

# Test user definition methods
@testset "User Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, container = SparseAxisArray)
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
        add_parameter_bound(x, par, 0, 0)
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
        delete_supports.(pars)
    end
    # test _check_var_bounds (HoldVariableRef Fallback)
    @testset "_check_var_bounds (Hold Fallback)" begin
        warn = "Unable to check if hold variables bounds are valid in measure " *
               "with custom measure data type BadData."
        @test_logs (:warn, warn) InfiniteOpt._check_var_bounds(x, BadData())
    end
    # test _check_var_bounds (MeasureRef)
    @testset "_check_var_bounds (Measure)" begin
        # make some measures
        meas1 = Measure(par, data)
        mref1 = add_measure(m, meas1)
        meas2 = Measure(x - mref1, data)
        mref2 = add_measure(m, meas2)
        # test normal
        @set_parameter_bounds(x, par == 1)
        @test isa(InfiniteOpt._check_var_bounds(mref1, data), Nothing)
        @test isa(InfiniteOpt._check_var_bounds(mref2, data), Nothing)
        delete_supports(par)
        set_parameter_bounds(x, ParameterBounds(), force = true)
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
        # test errors
        @test_throws ErrorException measure(x, data)
        @test_throws ErrorException measure(zero(GenericAffExpr{Float64,
                                                     GeneralVariableRef}), data)
        @test_throws ErrorException measure(par2, data)
        @test_throws ErrorException measure(inf4 + measure(inf + x, data3), data)
        # test with bad variable bounds
        InfiniteOpt._update_variable_param_bounds(x, ParameterBounds(Dict(par => IntervalSet(0, 0))))
        @test_throws ErrorException measure(inf + x, data)
    end
end

# Test if used
@testset "Used" begin
    # initialize model and reference
    m = InfiniteModel()
    mref = MeasureRef(m, 1)
    m.meas_in_objective[JuMP.index(mref)] = false
    # used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(mref)
        # prepare use case
        m.meas_to_constrs[JuMP.index(mref)] = [1]
        # test used
        @test used_by_constraint(mref)
        # undo changes
        delete!(m.meas_to_constrs, JuMP.index(mref))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(mref)
        # prepare use case
        m.meas_to_meas[JuMP.index(mref)] = [2]
        # test used
        @test used_by_measure(mref)
        # undo changes
        delete!(m.meas_to_meas, JuMP.index(mref))
    end
    # used_by_objective
    @testset "used_by_objective" begin
        # test not used
        @test !used_by_objective(mref)
        # prepare use case
        m.meas_in_objective[JuMP.index(mref)] = true
        # test used
        @test used_by_objective(mref)
        # undo changes
        m.meas_in_objective[JuMP.index(mref)] = false
    end
    # is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(mref)
        # prepare use case
        m.meas_to_constrs[JuMP.index(mref)] = [1]
        # test used
        @test is_used(mref)
        # undo changes
        delete!(m.meas_to_constrs, JuMP.index(mref))
        # prepare use case
        m.meas_to_meas[JuMP.index(mref)] = [2]
        # test used
        @test is_used(mref)
        # undo changes
        delete!(m.meas_to_meas, JuMP.index(mref))
        # prepare use case
        m.meas_in_objective[JuMP.index(mref)] = true
        # test used
        @test is_used(mref)
        # undo changes
        m.meas_in_objective[JuMP.index(mref)] = false
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
    @hold_variable(m, x)

    # test measure that does not use AbstractMeasureData inputs
    @testset "measure (no AbstractMeasureData)" begin
        meas1 = measure(inf, num_supports = 5, eval_method = Gauss_Legendre)
        (expected_supps, expected_coeffs) = FGQ.gausslegendre(5)
        expected_supps = expected_supps .* 0.5 .+ 0.5
        expected_coeffs = expected_coeffs .* 0.5
        @test all(measure_data(meas1).supports .== expected_supps)
        @test all(measure_data(meas1).coefficients .== expected_coeffs)

        add_supports(par2, [0.3, 0.7])
        meas2 = measure(inf2, par2, use_existing_supports = true)
        @test measure_data(meas2).parameter_ref == par2
        @test measure_data(meas2).supports == [0.3, 0.7]
        meas2 = measure(inf2, par2, 0.5, 0.9, use_existing_supports = true)
        @test measure_data(meas2).supports == [0.7]

        meas3 = measure(inf4, num_supports = 5)
        @test pars1[1] in measure_data(meas3).parameter_ref
        @test pars1[2] in measure_data(meas3).parameter_ref

        meas4 = measure(inf5, pars2, num_supports = 5)
        @test pars2["a"] in measure_data(meas4).parameter_ref
        @test pars2["b"] in measure_data(meas4).parameter_ref

        meas5 = measure(inf6, use_existing_supports = true)
        @test measure_data(meas4).supports == measure_data(meas5).supports

        # test errors
        @test_throws ErrorException measure(x)
        @test_throws ErrorException measure(inf, ParameterRef[])
        @test_throws ErrorException measure(inf2)
        @test_throws ErrorException measure(inf2, par, 1., 3.)
        @test_throws ErrorException measure(inf2, par, [0., 1.])
        @test_throws ErrorException measure(inf2, par, 0., [1., 1.])
        @test_throws ErrorException measure(inf2, par, 0.5, 0.)
        @test_throws ErrorException measure(inf7, use_existing_supports = true)
        @test_throws ErrorException measure(meas1)
    end
    # test support_sum
    @testset "support_sum" begin
        sum1 = support_sum(inf2, par2)
        @test measure_data(sum1).parameter_ref == par2
        @test measure_data(sum1).supports == [0.3, 0.7]
        @test name(sum1) == "sum(inf2(par, par2))"
        add_supports(pars3[1], [0.3, 0.7])
        add_supports(pars3[2], [0.3, 0.7])
        sum2 = support_sum(inf7)
        @test pars3[1] in measure_data(sum2).parameter_ref
        @test pars3[2] in measure_data(sum2).parameter_ref
        supps = [JuMPC.SparseAxisArray(Dict([((1,),0.3), ((2,), 0.3)])),
                 JuMPC.SparseAxisArray(Dict([((1,),0.7), ((2,), 0.7)]))]
        @test measure_data(sum2).supports == supps
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
end
