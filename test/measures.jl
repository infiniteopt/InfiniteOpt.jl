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
    @infinite_variable(m, inf(par))
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1], name = "test")
    meas = Measure(par + 2inf - x, data)
    # _make_meas_name
    @testset "_make_meas_name" begin
        @test InfiniteOpt._make_meas_name(meas) == "test(par + 2 inf(par) - x)"
    end
    # test _update_var_meas_mapping
    @testset "_update_var_meas_mapping" begin
        # retrieve variables and parameters
        mref = MeasureRef(m, 2)
        vars = [par, inf, x, mref]
        # test initial reference
        @test isa(InfiniteOpt._update_var_meas_mapping(vars, 1), Nothing)
        @test m.var_to_meas[index(inf)] == [1]
        @test m.var_to_meas[index(x)] == [1]
        @test m.param_to_meas[index(par)] == [1]
        @test m.meas_to_meas[index(mref)] == [1]
        # test repeated reference
        @test isa(InfiniteOpt._update_var_meas_mapping(vars, 2), Nothing)
        @test m.var_to_meas[index(inf)] == [1, 2]
        @test m.var_to_meas[index(x)] == [1, 2]
        @test m.param_to_meas[index(par)] == [1, 2]
        @test m.meas_to_meas[index(mref)] == [1, 2]
        # clear additions
        delete!(m.var_to_meas, index(inf))
        delete!(m.var_to_meas, index(x))
        delete!(m.param_to_meas, index(par))
        delete!(m.meas_to_meas, index(mref))
    end
    # add_measure
    @testset "add_measure" begin
        mref = MeasureRef(m, 1)
        @test add_measure(m, meas) == mref
        @test name(mref) == "test(par + 2 inf(par) - x)"
        @test m.var_to_meas[index(inf)] == [1]
        @test m.var_to_meas[index(x)] == [1]
        @test m.param_to_meas[index(par)] == [1]
        @test !m.meas_in_objective[index(mref)]
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
    m.meas_in_objective[index(mref)] = false
    # used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(mref)
        # prepare use case
        m.meas_to_constrs[index(mref)] = [1]
        # test used
        @test used_by_constraint(mref)
        # undo changes
        delete!(m.meas_to_constrs, index(mref))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(mref)
        # prepare use case
        m.meas_to_meas[index(mref)] = [2]
        # test used
        @test used_by_measure(mref)
        # undo changes
        delete!(m.meas_to_meas, index(mref))
    end
    # used_by_objective
    @testset "used_by_objective" begin
        # test not used
        @test !used_by_objective(mref)
        # prepare use case
        m.meas_in_objective[index(mref)] = true
        # test used
        @test used_by_objective(mref)
        # undo changes
        m.meas_in_objective[index(mref)] = false
    end
    # is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(mref)
        # prepare use case
        m.meas_to_constrs[index(mref)] = [1]
        # test used
        @test is_used(mref)
        # undo changes
        delete!(m.meas_to_constrs, index(mref))
        # prepare use case
        m.meas_to_meas[index(mref)] = [2]
        # test used
        @test is_used(mref)
        # undo changes
        delete!(m.meas_to_meas, index(mref))
        # prepare use case
        m.meas_in_objective[index(mref)] = true
        # test used
        @test is_used(mref)
        # undo changes
        m.meas_in_objective[index(mref)] = false
    end
end

# Test methods for reduced infinite variables
@testset "Reduced Variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    index = m.next_var_index + 1
    rvref = InfiniteOpt._ReducedInfiniteRef(m, index, inf2, Dict(1 => 1))
    rvref2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf2, Dict(2 => 1))
    # test name
    @testset "JuMP.name" begin
        @test name(rvref) == "inf2(1, par2)"
        @test name(rvref2) == "inf2(par1, 1)"
    end
end

# Test measure expansion methods
@testset "Expansion" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_parameter(m, 1 <= pars2[1:2] <= 2)
    @infinite_variable(m, inf1(par1))
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @infinite_variable(m, inf5(pars1, pars2))
    @infinite_variable(m, inf6(pars2))
    @infinite_variable(m, inf7(par1, par2, pars1))
    @global_variable(m, x)
    # prepare measures
    w(t) = 3
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(par2, [0.5, 0.5], [1, 1], weight_function = w)
    data3 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    data4 = DiscreteMeasureData(pars2, [2, 2], [[1, 1], [1, 1]])
    meas1 = measure(inf1, data1)
    meas2 = measure(2 * inf3 * x - par1, data2)
    # prepare point_mapper function and info
    function mapper(m2, pvref, ivref, support)
        set_name(pvref, string(InfiniteOpt._root_name(ivref), support))
        return
    end
    trans_model = Model()
    map_args = (trans_model, mapper)
    # test _make_point_variable
    @testset "_make_point_variable" begin
        index = m.next_var_index + 1
        @test InfiniteOpt._make_point_variable(inf1) == PointVariableRef(m,
                                                                         index)
    end
    # test _make_reduced_variable (from ivref)
    @testset "_make_reduced_variable (from ivref)" begin
        index = m.next_var_index + 1
        rvref = InfiniteOpt._ReducedInfiniteRef(m, index, inf2, Dict(1 => 1))
        @test InfiniteOpt._make_reduced_variable(inf2, 1, 1) == rvref
    end
    # test _make_reduced_variable (from rvref)
    @testset "_make_reduced_variable (from rvref)" begin
        index = m.next_var_index + 1
        rvref = InfiniteOpt._ReducedInfiniteRef(m, index, inf2, Dict(1 => 1))
        @test InfiniteOpt._make_reduced_variable(inf2, Dict(1 => 1)) == rvref
    end
    # test _expand_measure (infinite variable)
    @testset "_expand_measure (Infinite Variable)" begin
        # test single param infinite var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(inf1, data1, map_args...) == expected
        # test single param infinite var without measure param
        expected = 3inf1
        @test InfiniteOpt._expand_measure(inf1, data2, map_args...) == expected
        # test single param infinite var with measure param and others
        index = m.next_var_index + 1
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf2, Dict(1 => 1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf2, Dict(1 => 2))
        expected = 0.5 * (rv1 + rv2)
        @test InfiniteOpt._expand_measure(inf2, data1, map_args...) == expected
        # test array param infinite var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(inf4, data3, map_args...) == expected
        # test array param infinite var without measure param
        expected = 4inf4
        @test InfiniteOpt._expand_measure(inf4, data4, map_args...) == expected
        # test array param infinite var with measure param and others
        index = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf5, Dict(1 => supp1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf5,
                                              Dict(1 => supp2))
        expected = rv1 + rv2
        @test InfiniteOpt._expand_measure(inf5, data3, map_args...) == expected
    end
    # test _expand_measure (reduced infinite variable)
    @testset "_expand_measure (Reduced Variable)" begin
        # test single param reduced var without measure param
        rv = InfiniteOpt._ReducedInfiniteRef(m, -1, inf2, Dict(1 => 1))
        expected = 0.5 * (rv + rv)
        @test InfiniteOpt._expand_measure(rv, data1, map_args...) == expected
        # test single param reduced var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 1.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data2, map_args...) == expected
        # test single param reduced var with measure param and others
        rv = InfiniteOpt._ReducedInfiniteRef(m, -2, inf7, Dict(1 => 1))
        index = m.next_var_index + 1
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf7,
                                              Dict(1 => 1, 2 => 1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf7,
                                              Dict(1 => 2, 2 => 1))
        expected = 1.5 * (rv1 + rv2)
        @test InfiniteOpt._expand_measure(rv, data2, map_args...) == expected
        # test array param reduced var without measure param
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        rv = InfiniteOpt._ReducedInfiniteRef(m, -1, inf5, Dict(1 => supp1))
        expected = rv + rv
        @test InfiniteOpt._expand_measure(rv, data3, map_args...) == expected
        # test array param reduced var with measure param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 2 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(rv, data4, map_args...) == expected
        # test array param reduced var with measure param and others
        rv = InfiniteOpt._ReducedInfiniteRef(m, -2, inf7, Dict(1 => 1))
        index = m.next_var_index + 1
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf7,
                                              Dict(1 => 1, 3 => supp1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf7,
                                              Dict(1 => 2, 3 => supp2))
        expected = rv1 + rv2
        @test InfiniteOpt._expand_measure(rv, data3, map_args...) == expected
    end
    # test _expand_measure (finite variable)
    @testset "_expand_measure (Finite Variable)" begin
        # test with single parameter measure
        expected = 0.5 * (x + x)
        @test InfiniteOpt._expand_measure(x, data1, map_args...) == expected
        # test with multi parameter measure
        expected = x + x
        @test InfiniteOpt._expand_measure(x, data3, map_args...) == expected
    end
    # test _expand_measure (parameter)
    @testset "_expand_measure (Parameter)" begin
        # test with different parameter
        expected = 1.5 * (par1 + par1)
        @test InfiniteOpt._expand_measure(par1, data2, map_args...) == expected
        # test with same parameter
        expected = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        expected.constant = 0.5 * (1 + 2)
        @test InfiniteOpt._expand_measure(par1, data1, map_args...) == expected
    end
    # test _expand_measure (multi measure data with parameter)
    @testset "_expand_measure (Multi Parameter)" begin
        # test with different parameter
        expected = par1 + par1
        @test InfiniteOpt._expand_measure(par1, data3, map_args...) == expected
        # test with same parameter
        expected = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        expected.constant = 1 + 2
        @test InfiniteOpt._expand_measure(pars1[2], data3, map_args...) == expected
    end
    # test _expand_measure (AffExpr)
    @testset "_expand_measure (AffExpr)" begin
        # test single param AffExpr, no measures
        expr = 2inf1 + par1 - x + 3par2 - 3
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (2pts[1] + 2pts[2] + 1 + 2 - x - x + 3par2 + 3par2 - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param AffExpr, with measures
        expr = meas2 - x + par1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (6 * (pts[1] * x + pts[2] * x) - 2x - 6)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test array param AffExpr, no measures
        expr = inf4 + par1 - x + 3pars1[2] - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2] + 2par1 - 2x + 9 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param AffExpr, with measures
        expr = meas2 - x + par1 - inf4
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1),
               PointVariableRef(m, index + 2), PointVariableRef(m, index + 3)]
        expected = 6 * (pts[1] * x + pts[2] * x) - 4par1 - 2x - pts[3] - pts[4]
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
    end
    # test _expand_measure (QuadExpr)
    @testset "_expand_measure (QuadExpr)" begin
        # test single param QuadExpr with both variables integrated or not
        expr = 2 * inf1 * inf2 - inf3 * inf4 + x + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index + 2, inf2, Dict(1 => 1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 3, inf2, Dict(1 => 2))
        expected = 0.5 * (2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf3 *
                          inf4 + 2x + 4)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable not integrated
        expr = 3 * inf3 * inf1 + pars1[1] - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (3 * inf3 * pts[1] + 3 * inf3 * pts[2] + 2pars1[1] - 2)
        @test InfiniteOpt._expand_measure(expr, data1, map_args...) == expected
        # test single param QuadExpr with first variable integrated
        expr = 3 * inf2 * inf1 + pars1[1] + 1
        index = m.next_var_index + 1
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf2, Dict(2 => 1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf2, Dict(2 => 1))
        expected = 1.5 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data2, map_args...) == expected
        # test array param QuadExpr with both variables integrated
        expr = 2 * inf4 * inf5 - inf1 * inf2 + x + 2
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        supp1 = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        supp2 = convert(JuMP.Containers.SparseAxisArray, [2, 2])
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index + 2, inf5,
                                              Dict(1 => supp1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 3, inf5,
                                              Dict(1 => supp2))
        expected = 2 * pts[1] * rv1 + 2 * pts[2] * rv2 - 2 * inf1 * inf2 + 2x + 4
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable not integrated
        expr = 3 * inf1 * inf4 + par1 - 1
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 3 * inf1 * pts[1] + 3 * inf1 * pts[2] + 2par1 - 2
        @test InfiniteOpt._expand_measure(expr, data3, map_args...) == expected
        # test array param QuadExpr with first variable integrated
        expr = 3 * inf5 * inf1 + pars1[1] + 1
        index = m.next_var_index + 1
        supp = convert(JuMP.Containers.SparseAxisArray, [1, 1])
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index, inf5, Dict(2 => supp))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 1, inf5, Dict(2 => supp))
        expected = 2 * (3 * rv1 * inf1 + 3 * rv2 * inf1 + 2pars1[1] + 2)
        @test InfiniteOpt._expand_measure(expr, data4, map_args...) == expected
    end
    # test _expand_measure (measure)
    @testset "_expand_measure (Measure)" begin
        # test single param measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (pts[1] + pts[2])
        @test InfiniteOpt._expand_measure(meas1, data1, map_args...) == expected
        # test another single param
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 0.5 * (6 * pts[1] * x + 6 * pts[2] * x - 3 * (1 + 2))
        @test InfiniteOpt._expand_measure(meas2, data1, map_args...) == expected
        # test array param measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = pts[1] + pts[2]
        @test InfiniteOpt._expand_measure(meas1, data3, map_args...) == expected
    end
    # test _expand_measure (other)
    @testset "_expand_measure (Other)" begin
        # prepare test
        @variable(trans_model, y)
        struct new <: AbstractMeasureData end
        # test it
        @test_throws ErrorException InfiniteOpt._expand_measure(y, new(),
                                                                map_args...)
    end
end

# Test user expansion methods
@testset "User Expansion" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_parameter(m, 1 <= pars1[1:2] <= 2)
    @infinite_variable(m, inf1(par1) >= 1)
    @infinite_variable(m, inf2(par1, par2))
    @infinite_variable(m, inf3(par2))
    @infinite_variable(m, inf4(pars1))
    @global_variable(m, x)
    # prepare measures
    data1 = DiscreteMeasureData(par1, [0.5, 0.5], [1, 2])
    data2 = DiscreteMeasureData(pars1, [1, 1], [[1, 1], [2, 2]])
    meas1 = measure(inf1 + 3x - inf2 + inf3 - 2, data1)
    meas2 = measure(2 * inf4 * x - pars1[2] + inf2, data2)
    # test _add_mapped_point_variable
    @testset "_add_mapped_point_variable" begin
        # prepare the partially made point variable reference
        pref = InfiniteOpt._make_point_variable(inf1)
        # try mapping/adding it
        @test isa(InfiniteOpt._add_mapped_point_variable(m, pref, inf1, (1,)),
                  Nothing)
        @test name(pref) == "inf1(1)"
        @test haskey(m.vars, index(pref))
        @test lower_bound(pref) == 1
    end
    # test expand
    @testset "expand" begin
        # test the first measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        rv1 = InfiniteOpt._ReducedInfiniteRef(m, index + 2, inf2, Dict(1 => 1))
        rv2 = InfiniteOpt._ReducedInfiniteRef(m, index + 3, inf2, Dict(1 => 2))
        expected = 0.5 * (pts[1] + pts[2] + 6x - rv1 - rv2 + 2inf3 - 4)
        @test expand(meas1) == expected
        # test the second measure
        index = m.next_var_index + 1
        pts = [PointVariableRef(m, index), PointVariableRef(m, index + 1)]
        expected = 2 * pts[1] * x + 2 * pts[2] * x - 3 + 2inf2
        @test expand(meas2) == expected
        # test error
        m2 = Model()
        @variable(m2, y)
        struct new <: AbstractMeasureData end
        m.measures[1] = Measure(y, new())
        @test_throws ErrorException expand(meas1)
    end
end
