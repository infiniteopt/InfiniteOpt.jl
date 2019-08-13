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
    @global_variable(m, x)
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
               "other. Parameter deletion methods should not be used."
        @test_logs (:warn, warn) InfiniteOpt._update_param_data_mapping(m,
                                                                   BadData(), 1)
    end
    # test add_measure
    @testset "add_measure" begin
        mref = MeasureRef(m, 1)
        @test add_measure(m, meas) == mref
        @test name(mref) == "test(par + 2 inf(par) - x)"
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
    @global_variable(m, x)
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
    @global_variable(m, x)
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
        @test isa(InfiniteOpt._check_has_parameter(par, par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(inf, par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(x - inf, par), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(x + mref2 + mref4, par),
                  Nothing)
        # test is does not have the parameter
        @test_throws ErrorException InfiniteOpt._check_has_parameter(par, par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter(inf, par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter(x - inf,
                                                                     par2)
        @test_throws ErrorException InfiniteOpt._check_has_parameter(x + mref1 +
                                                                    mref5, par2)
    end
    # _check_has_parameter (array)
    @testset "_check_has_parameter (array)" begin
        # test that is has the parameter
        @test isa(InfiniteOpt._check_has_parameter(pars[1] + pars[2], pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(inf4, pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(x - inf4, pars), Nothing)
        @test isa(InfiniteOpt._check_has_parameter(x + mref6, pars), Nothing)
        # test is does not have the parameter
        @test_throws ErrorException InfiniteOpt._check_has_parameter(pars[1], pars)
        @test_throws ErrorException InfiniteOpt._check_has_parameter(inf2, pars)
        @test_throws ErrorException InfiniteOpt._check_has_parameter(x + mref1 +
                                                                    mref5, pars)
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
    @global_variable(m, x)
    # prepare measure data
    data = DiscreteMeasureData(par, [1], [1], name = "a")
    data2 = DiscreteMeasureData(par2, [1], [1], name = "b")
    data3 = DiscreteMeasureData(pars, [1], [[1, 1]], name = "c")
    # test _model_from_expr
    @testset "_model_from_expr" begin
        @test InfiniteOpt._model_from_expr(par) == m
        @test InfiniteOpt._model_from_expr(inf + par) == m
        @test InfiniteOpt._model_from_expr(inf^2) == m
        @test isa(InfiniteOpt._model_from_expr(zero(GenericAffExpr{Float64,
                                                 GeneralVariableRef})), Nothing)
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
    # test measure
    @testset "measure" begin
        # test single use
        mref = MeasureRef(m, 1)
        @test measure(inf + x, data) == mref
        @test name(mref) == "a(inf(par) + x)"
        @test supports(par) == [1]
        @test !m.meas_in_objective[1]
        # test nested use
        mref2 = MeasureRef(m, 2)
        mref3 = MeasureRef(m, 3)
        @test measure(inf + measure(inf2 + x, data2), data) == mref3
        @test name(mref3) == "a(inf(par) + b(inf2(par, par2) + x))"
        @test supports(par) == [1]
        @test supports(par2) == [1]
        # test vector parameter
        mref4 = MeasureRef(m, 4)
        @test measure(inf4 + x, data3) == mref4
        @test name(mref4) == "c(inf4(pars) + x)"
        @test supports(pars[1]) == [1]
        @test supports(pars[2]) == [1]
        # test errors
        @test_throws ErrorException measure(x, data)
        @test_throws ErrorException measure(zero(GenericAffExpr{Float64,
                                                     GeneralVariableRef}), data)
        @test_throws ErrorException measure(par2, data)
        @test_throws ErrorException measure(inf4 + measure(inf + x, data3), data)
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
