# Test basic extensions
@testset "Basics" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10])
    con1 = ScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0))
    con2 = BoundedScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0),
                                  ParameterBounds((t => IntervalSet(0, 1),)),
                                  ParameterBounds())
    object1 = ConstraintData(con1, Int[], "c1", MeasureIndex[], false)
    object2 = ConstraintData(con2, [1], "c2", MeasureIndex[], false)
    idx1 = ConstraintIndex(1)
    idx2 = ConstraintIndex(2)
    cref1 = InfOptConstraintRef(m, idx1)
    cref2 = InfOptConstraintRef(m, idx2)
    bad_cref = InfOptConstraintRef(m, ConstraintIndex(-1))
    # owner_model
    @testset "JuMP.owner_model" begin
      @test owner_model(cref1) === m
      @test owner_model(cref2) === m
    end
    # index
    @testset "JuMP.index" begin
      @test JuMP.index(cref1) == ConstraintIndex(1)
      @test JuMP.index(cref2) == ConstraintIndex(2)
    end
    # test Base.:(==) of constraint references
    @testset "Base.:(==) References" begin
        @test !(cref1 == cref2)
        @test cref1 == InfOptConstraintRef(m, idx1)
        @test cref1 != InfOptConstraintRef(m, idx2)
        @test cref1 != InfOptConstraintRef(InfiniteModel(), idx1)
    end
    # test broadcastable
    @testset "Base.broadcastable Reference" begin
        @test isa(Base.broadcastable(cref1), Base.RefValue)
    end
    # test constraint type
    @testset "JuMP.constraint_type" begin
        @test constraint_type(m) == InfOptConstraintRef
    end
    # test shape
    @testset "JuMP.shape" begin
        @test JuMP.shape(con2) == ScalarShape()
    end
    # test jump_function
    @testset "JuMP.jump_function" begin
        @test jump_function(con2) == zero(AffExpr)
    end
    # test moi_set
    @testset "JuMP.moi_set" begin
        @test moi_set(con2) == MOI.EqualTo(0.0)
    end
    # test parameter_bounds
    @testset "parameter_bounds" begin
        @test parameter_bounds(con2) == ParameterBounds((t => IntervalSet(0, 1),))
    end
    # test original_parameter_bounds
    @testset "original_parameter_bounds" begin
        @test original_parameter_bounds(con2) == ParameterBounds()
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object1) == idx1
        @test InfiniteOpt._add_data_object(m, object2) == idx2
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(cref1) === m.constraints
    end
    # is_valid
    @testset "JuMP.is_valid" begin
      @test is_valid(m, cref1)
      @test is_valid(m, cref2)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(cref1) === object1
        @test InfiniteOpt._data_object(cref2) === object2
        @test_throws ErrorException InfiniteOpt._data_object(bad_cref)
    end
    # _core_constraint_object
    @testset "_core_constraint_object" begin
        @test InfiniteOpt._core_constraint_object(cref1) === con1
        @test InfiniteOpt._core_constraint_object(cref2) === con2
    end
    # _set_core_constraint_object
    @testset "_set_core_constraint_object" begin
        @test InfiniteOpt._set_core_constraint_object(cref1, con2) isa Nothing
        @test InfiniteOpt._set_core_constraint_object(cref1, con1) isa Nothing
    end
    # _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(cref1) == []
        @test InfiniteOpt._object_numbers(cref2) == [1]
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(cref1) == MeasureIndex[]
    end
    @testset "_is_info_constraint" begin
        @test !InfiniteOpt._is_info_constraint(cref2)
    end
    # constraint_object
    @testset "JuMP.constraint_object" begin
      @test constraint_object(cref1) == con1
      @test constraint_object(cref2) == con2
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(cref1) == "c1"
        @test name(cref2) == "c2"
        @test name(bad_cref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(cref1, "new"), Nothing)
        @test_throws ErrorException set_name(bad_cref, "")
        @test name(cref1) == "new"
    end
    # test has_parameter_bounds
    @testset "has_parameter_bounds" begin
        @test has_parameter_bounds(cref2)
        @test !has_parameter_bounds(cref1)
    end
    # test parameter_bounds
    @testset "parameter_bounds" begin
        @test_throws ErrorException parameter_bounds(cref1)
        @test parameter_bounds(cref2) == ParameterBounds((t => IntervalSet(0, 1),))
    end
    # _make_constraint_ref
    @testset "_make_constraint_ref" begin
        @test InfiniteOpt._make_constraint_ref(m, idx1) == cref1
        @test InfiniteOpt._make_constraint_ref(m, idx2) == cref2
    end
    # constraint_by_name
    @testset "JuMP.constraint_by_name" begin
        # test normal
        @test constraint_by_name(m, "new") == cref1
        @test isa(constraint_by_name(m, "test2"), Nothing)
        # prepare constraint with duplicate name
        @test isa(set_name(cref2, "new"), Nothing)
        # test for duplciate name error
        @test_throws ErrorException constraint_by_name(m, "new")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(cref1) isa Nothing
        @test length(InfiniteOpt._data_dictionary(cref1)) == 1
        @test !is_valid(m, cref1)
    end
end

# Test definition methods
@testset "Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, pars))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x, parameter_bounds = (par == 1))
    var = build_variable(error, inf2, Dict{Int, Float64}(1 => 0.5), check = false)
    rv = add_variable(m, var)
    data = TestData(par, 1, 1)
    meas = @measure(inf + par - x, data, name = "test")
    mindex = MeasureIndex(1)
    dinf = @deriv(inf, par)
    @testset "JuMP.build_constraint (Single)" begin
        # test bounded constraint
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        @test jump_function(build_constraint(error, inf, MOI.EqualTo(42.0),
                                             parameter_bounds = bounds)) == inf
        @test moi_set(build_constraint(error, inf, MOI.EqualTo(42.0),
                                 parameter_bounds = bounds)) == MOI.EqualTo(42.0)
        @test parameter_bounds(build_constraint(error, inf, MOI.EqualTo(42.0),
                                         parameter_bounds = bounds)) === bounds
        @test original_parameter_bounds(build_constraint(error, inf,
                       MOI.EqualTo(42.0), parameter_bounds = bounds)) !== bounds
        @test original_parameter_bounds(build_constraint(error, inf,
                        MOI.EqualTo(42.0), parameter_bounds = bounds)) == bounds
        # test scalar constraint
        con = ScalarConstraint(inf, MOI.EqualTo(42.0))
        @test jump_function(build_constraint(error, inf, MOI.EqualTo(42.0))) == inf
        @test moi_set(build_constraint(error, inf, MOI.EqualTo(42.0))) == MOI.EqualTo(42.0)
    end
    # test build_constraint (expr)
    @testset "JuMP.build_constraint (Expr)" begin
        # test bounded constraint
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        @test jump_function(build_constraint(error, inf + par, MOI.EqualTo(42.0),
                                             parameter_bounds = bounds)) == inf + par
        @test moi_set(build_constraint(error, inf + par, MOI.EqualTo(42.0),
                                 parameter_bounds = bounds)) == MOI.EqualTo(42.0)
        @test parameter_bounds(build_constraint(error, inf + par, MOI.EqualTo(42.0),
                                         parameter_bounds = bounds)) === bounds
        @test original_parameter_bounds(build_constraint(error, inf + par,
                       MOI.EqualTo(42.0), parameter_bounds = bounds)) !== bounds
        @test original_parameter_bounds(build_constraint(error, inf + par,
                        MOI.EqualTo(42.0), parameter_bounds = bounds)) == bounds
        # test scalar constraint
        con = ScalarConstraint(inf + par, MOI.EqualTo(42.0))
        @test jump_function(build_constraint(error, inf + par, MOI.EqualTo(42.0))) == inf + par
        @test moi_set(build_constraint(error, inf + par, MOI.EqualTo(42.0))) == MOI.EqualTo(42.0)
    end
    # test _check_and_update_bounds (BoundedScalarConstraint)
    @testset "_check_and_update_bounds (Bounded)" begin
        bounds = ParameterBounds()
        c = BoundedScalarConstraint(x + inf, MOI.EqualTo(0.0), bounds, copy(bounds))
        @test InfiniteOpt._check_and_update_bounds(m, c,
                                      [x, inf]).bounds[par] == IntervalSet(1, 1)
        @test original_parameter_bounds(InfiniteOpt._check_and_update_bounds(m, c,
                                                 [x, inf])) == ParameterBounds()
    end
    # test _check_and_update_bounds (ScalarConstraint)
    @testset "_check_and_update_bounds (Scalar)" begin
        c = JuMP.ScalarConstraint(inf + x, MOI.EqualTo(0.0))
        @test InfiniteOpt._check_and_update_bounds(m, c,
                            [x, inf]).bounds[par] == IntervalSet(1, 1)
        set_parameter_bounds(x, ParameterBounds(), force = true)
        m.has_hold_bounds = false
        @test InfiniteOpt._check_and_update_bounds(m, c, [x, inf]) == c
    end
    # test _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # add dumby constraints
        con = ScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0))
        object1 = ConstraintData(con, Int[], "c", MeasureIndex[], false)
        object2 = ConstraintData(con, Int[], "c", MeasureIndex[], false)
        idx1 = ConstraintIndex(1)
        idx2 = ConstraintIndex(2)
        cref1 = InfOptConstraintRef(m, idx1)
        cref2 = InfOptConstraintRef(m, idx2)
        @test InfiniteOpt._add_data_object(m, object1) == idx1
        @test InfiniteOpt._add_data_object(m, object2) == idx2
        # test initial use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], cref1),
                  Nothing)
        @test InfiniteOpt._constraint_dependencies(inf) == [idx1]
        @test InfiniteOpt._constraint_dependencies(par) == [idx1]
        @test InfiniteOpt._constraint_dependencies(meas) == [idx1]
        @test InfiniteOpt._constraint_dependencies(rv) == [idx1]
        @test InfiniteOpt._measure_dependencies(cref1) == [mindex]
        # test secondary use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], cref2),
                  Nothing)
        @test InfiniteOpt._constraint_dependencies(inf) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(par) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(meas) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(rv) == [idx1, idx2]
        @test InfiniteOpt._measure_dependencies(cref2) == [mindex]
        # Undo changes
        empty!(InfiniteOpt._constraint_dependencies(inf))
        empty!(InfiniteOpt._constraint_dependencies(par))
        empty!(InfiniteOpt._constraint_dependencies(meas))
        empty!(InfiniteOpt._constraint_dependencies(rv))
        empty!(InfiniteOpt._measure_dependencies(cref1))
        empty!(InfiniteOpt._measure_dependencies(cref2))
        m.constraints = MOIUC.CleverDict{ConstraintIndex, ConstraintData}()
    end
    # test add_constraint
    @testset "JuMP.add_constraint" begin
        # test reject vector constraint
        con = VectorConstraint([inf + pt], MOI.Zeros(1))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test rejects bad variables
        @hold_variable(InfiniteModel(), z)
        con = ScalarConstraint(z, MOI.EqualTo(42.0))
        @test_throws VariableNotOwned{GeneralVariableRef} add_constraint(m, con, "a")
        # test bounded constraint
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0), bounds, bounds)
        idx = ConstraintIndex(1)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "a") == cref
        # test infinite constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        idx = ConstraintIndex(2)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "b") == cref
        # test finite constraint
        con = ScalarConstraint(x + pt, MOI.EqualTo(42.0))
        idx = ConstraintIndex(3)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "d") == cref
        @test name(cref) == "d"
        @test !InfiniteOpt._is_info_constraint(cref)
        @test !optimizer_model_ready(m)
        @test used_by_constraint(pt)
        # test with bounded hold variables
        @set_parameter_bounds(x, par == 1)
        con = ScalarConstraint(inf + pt + x, MOI.EqualTo(42.0))
        idx = ConstraintIndex(4)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "b") == cref
        @test parameter_bounds(cref)[par] == IntervalSet(1, 1)
        @test InfiniteOpt._core_constraint_object(cref).orig_bounds == ParameterBounds()
        set_parameter_bounds(x, ParameterBounds(), force = true)
    end
    # test macro
    @testset "JuMP.@constraint" begin
        # test scalar constraint
        idx = ConstraintIndex(5)
        cref = InfOptConstraintRef(m, idx)
        @test @constraint(m, e, x + pt -2 <= 2) == cref
        @test isa(InfiniteOpt._core_constraint_object(cref), ScalarConstraint)
        # test bounded scalar constraint
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        idx = ConstraintIndex(6)
        cref = InfOptConstraintRef(m, idx)
        @test @constraint(m, f, inf + meas - dinf <= 2, parameter_bounds = bounds) == cref
        @test InfiniteOpt._object_numbers(cref) == [1]
        @test isa(InfiniteOpt._core_constraint_object(cref), BoundedScalarConstraint)
        @test used_by_constraint(dinf)
    end
end

# Test all the methods for bound constraint macros
@testset "Bounded Constraint Macros" begin
    # test _parse_name_expression
    @testset "_parse_name_expression" begin
        @test InfiniteOpt._parse_name_expression(error, :bob) == :bob
        @test InfiniteOpt._parse_name_expression(error, :([i = 1:2])) == :([i = 1:2])
        @test InfiniteOpt._parse_name_expression(error, :(bob[i = 1:2])) == :(bob[i = 1:2])
        @test_throws ErrorException InfiniteOpt._parse_name_expression(error, :((2, 3)))
    end
    # test _extract_bounds (:call)
    @testset "_extract_bounds (:call)" begin
        # test single anonymous bound
        args = [:in, :t, :([0, 1])]
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(($(dict_arg),)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (nothing, bounds)
        # test with name expression
        args = [:(bob[i = 1:2]), :(t in [0, 1])]
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (args[1], bounds)
    end
    # test @BDconstraint
    @testset "@BDconstraint" begin
        # initialize model and references
        m = InfiniteModel()
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        @infinite_variable(m, inf(par))
        @point_variable(m, inf(0.5), pt)
        @hold_variable(m, x)
        bounds1 = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        m2 = Model()
        # test errors
        @test_macro_throws ErrorException @BDconstraint(m, par = 0, inf == 0,
                                                     parameter_bounds = bounds1)
        @test_macro_throws ErrorException @BDconstraint(m, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con[1:2], inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, [1:2], inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, a.b, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m2, par == 0, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m2, con(par == 0), inf == 0)
        # test anonymous constraint with set
        idx = ConstraintIndex(1)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, par in [0, 1], inf + x == 2) == cref
        @test parameter_bounds(cref) == bounds1
        # test anonymous constraint with comparison
        idx = ConstraintIndex(2)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, 0 <= par <= 1, inf + x == 2) == cref
        @test parameter_bounds(cref) == bounds1
        # test anonymous constraint with equality
        idx = ConstraintIndex(3)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, par == 0, inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test anonymous multiple constraint
        idxs = [ConstraintIndex(4), ConstraintIndex(5)]
        crefs = [InfOptConstraintRef(m, idx) for idx in idxs]
        @test @BDconstraint(m, [1:2](par == 0), inf + x == 2) == crefs
        @test parameter_bounds(crefs[1]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test parameter_bounds(crefs[2]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test anonymous multiple bounds
        idx = ConstraintIndex(6)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, (par == 0, pars[1] in [0, 1]), inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds((par => IntervalSet(0, 0),
                                                  pars[1] => IntervalSet(0, 1)))
        # test anonymous multiple bounds with vector input
        idx = ConstraintIndex(7)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, (par == 0, pars in [0, 1]), inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds((par => IntervalSet(0, 0),
                    pars[1] => IntervalSet(0, 1), pars[2] => IntervalSet(0, 1)))
        # test named constraint
        idx = ConstraintIndex(8)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, c1(par in [0, 1]), inf + x == 2) == cref
        @test parameter_bounds(cref) == bounds1
        # test named constraint with other vector type
        idx = ConstraintIndex(9)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, c2(par in [0, 1], 0 <= pars <= 1), inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds((par => IntervalSet(0, 1),
                    pars[1] => IntervalSet(0, 1), pars[2] => IntervalSet(0, 1)))
        # test named constraints
        idxs = [ConstraintIndex(10), ConstraintIndex(11)]
        crefs = [InfOptConstraintRef(m, idx) for idx in idxs]
        @test @BDconstraint(m, c3[1:2](par == 0), inf + x == 2) == crefs
        @test parameter_bounds(crefs[1]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test parameter_bounds(crefs[2]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test different container type
        idxs = [ConstraintIndex(12), ConstraintIndex(13)]
        crefs = [InfOptConstraintRef(m, idx) for idx in idxs]
        expected = convert(JuMPC.SparseAxisArray, crefs)
        @test @BDconstraint(m, c4[1:2](par == 0), inf + x == 2,
                            container = SparseAxisArray) == expected
        @test parameter_bounds(crefs[1]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test parameter_bounds(crefs[2]) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test with SparseAxisArray input
        idx = ConstraintIndex(14)
        cref = InfOptConstraintRef(m, idx)
        pars = convert(JuMPC.SparseAxisArray, pars)
        @test @BDconstraint(m, c5(pars in [0, 1]), inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds((pars[1] => IntervalSet(0, 1),
                                                  pars[2] => IntervalSet(0, 1)))
        # test with hold variable bounds
        @set_parameter_bounds(x, pars[1] in [0, 2])
        idx = ConstraintIndex(15)
        cref = InfOptConstraintRef(m, idx)
        @test @BDconstraint(m, c10(par in [0, 1]), inf + x == 2) == cref
        @test parameter_bounds(cref) == ParameterBounds(Dict(par => IntervalSet(0, 1),
                                                  pars[1] => IntervalSet(0, 2)))
        @test InfiniteOpt._core_constraint_object(cref).orig_bounds == ParameterBounds(Dict(par => IntervalSet(0, 1)))
    end
end

# Test parameter reference methods
@testset "Parameter References" begin
    m = InfiniteModel()
    @independent_parameter(m, t in [0, 1])
    @independent_parameter(m, y in [0, 1])
    @dependent_parameters(m, x[1:3] in [0, 1])
    @hold_variable(m, z)
    @constraint(m, c1, z >= 0)
    @constraint(m, c2, z + t + x[1] >= 0)
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(c1) == ()
        @test parameter_refs(c2) == (t, x)
    end
end

# Test parameter bound methods
@testset "Parameter Bounds" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 10)
    @infinite_parameter(m, 0 <= pars[1:2] <= 10)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x >= 0, Int)
    @BDconstraint(m, c1(par in [0, 1]), inf + x == 0)
    @constraint(m, c2, x * pt + x == 2)
    # test _update_constr_param_bounds
    @testset "_update_constr_param_bounds" begin
        # test empty bounds
        bounds = ParameterBounds()
        @test isa(InfiniteOpt._update_constr_param_bounds(c1, bounds, bounds),
                  Nothing)
        @test isa(constraint_object(c1), ScalarConstraint)
        # test normalized_coefficient
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 5)))
        @test isa(InfiniteOpt._update_constr_param_bounds(c1, bounds, bounds),
                  Nothing)
        @test parameter_bounds(c1) == bounds
    end
    # test set_parameter_bounds
    @testset "set_parameter_bounds" begin
        # test force error
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        @test_throws ErrorException set_parameter_bounds(c1, bounds)
        # test normal
        @test isa(set_parameter_bounds(c2, bounds), Nothing)
        @test parameter_bounds(c2) == bounds
        @test !optimizer_model_ready(m)
        # test test error with bounds
        bounds = ParameterBounds(Dict(par => IntervalSet(-1, 1)))
        @test_throws ErrorException set_parameter_bounds(c1, bounds, force = true)
    end
    # test add_parameter_bound
    @testset "add_parameter_bounds" begin
        # test already has bounds
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))
        @test isa(add_parameter_bounds(c1, copy(bounds)), Nothing)
        @test parameter_bounds(c1) == bounds
        # test doesn't have bounds
        @constraint(m, c3, inf + pt == 0)
        @test isa(add_parameter_bounds(c3, copy(bounds)), Nothing)
        @test parameter_bounds(c3) == bounds
    end
    # test delete_parameter_bounds
    @testset "delete_parameter_bounds" begin
        # test partial deletion
        bounds = ParameterBounds(Dict(pars[1] => IntervalSet(0, 1)))
        @test isa(add_parameter_bounds(c1, bounds), Nothing)
        @add_parameter_bounds(x, par == 0)
        @test isa(delete_parameter_bounds(c1), Nothing)
        expected = ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test parameter_bounds(c1) == expected
        # test again without hold bounds
        delete_parameter_bounds(x)
        @test isa(constraint_object(c1), ScalarConstraint)
        # test already gone
        @test isa(delete_parameter_bounds(c1), Nothing)
        @test isa(constraint_object(c1), ScalarConstraint)
    end
    # test @set_parameter_bounds
    @testset "@set_parameter_bounds" begin
        # Note errors were already checked with hold variables
        # test force error
        @test_macro_throws ErrorException @set_parameter_bounds(c1, par == 1)
        # test with single
        @set_parameter_bounds(x, pars == 0)
        @test isa(@set_parameter_bounds(c1, par == 1, force = true), Nothing)
        @test parameter_bounds(c1) == ParameterBounds((par => IntervalSet(1, 1),
                                                  pars => IntervalSet(0, 0)))
        set_parameter_bounds(x, ParameterBounds(), force = true)
        # test multiple
        @test isa(@set_parameter_bounds(c1, pars == 0, force = true), Nothing)
        @test parameter_bounds(c1)[pars[2]] == IntervalSet(0, 0)
    end
    # test @add_parameter_bounds
    @testset "@add_parameter_bounds" begin
        # Note errors were already checked with hold variables
        # test bounds error
        @test_macro_throws ErrorException @add_parameter_bounds(c3, par in [-1, 1])
        # test with multiple
        @test isa(@add_parameter_bounds(c2, (pars == 0, par in [0, 1])), Nothing)
        @test parameter_bounds(c2)[pars[2]] == IntervalSet(0, 0)
        @test parameter_bounds(c2)[par] == IntervalSet(0, 1)
    end
end

# Test coefficient queries and modifications
@testset "Coefficient Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x >= 0, Int)
    @BDconstraint(m, c1(par in [0, 1]), inf + x == 0)
    @constraint(m, c2, x * pt + x == 2)
    # test set_normalized_coefficient
    @testset "JuMP.set_normalized_coefficient" begin
        @test isa(set_normalized_coefficient(c1, x, 2), Nothing)
        @test isa(constraint_object(c1), BoundedScalarConstraint)
        @test jump_function(constraint_object(c1)) == inf + 2x
        @test isa(set_normalized_coefficient(c2, inf, 2), Nothing)
        @test isa(constraint_object(c2), ScalarConstraint)
        @test jump_function(constraint_object(c2)) == x * pt + x + 2inf
    end
    # test normalized_coefficient
    @testset "JuMP.normalized_coefficient" begin
        @test normalized_coefficient(c1, x) == 2.
        @test normalized_coefficient(c1, pt) == 0.
        @test normalized_coefficient(c2, inf) == 2.
        @test normalized_coefficient(LowerBoundRef(x), x) == 1.
        @test normalized_coefficient(LowerBoundRef(x), pt) == 0.
    end
    # test _set_set_value
    @testset "_set_set_value" begin
        @test InfiniteOpt._set_set_value(MOI.EqualTo(1.), 2) == MOI.EqualTo(2.)
        @test InfiniteOpt._set_set_value(MOI.GreaterThan(1.), 2) == MOI.GreaterThan(2.)
        @test InfiniteOpt._set_set_value(MOI.LessThan(1.), 2) == MOI.LessThan(2.)
    end
    # test set_normalized_rhs
    @testset "JuMP.set_normalized_rhs" begin
        @test isa(set_normalized_rhs(c1, 42), Nothing)
        @test isa(constraint_object(c1), BoundedScalarConstraint)
        @test constraint_object(c1).set == MOI.EqualTo(42.0)
        @test isa(set_normalized_rhs(c2, 42), Nothing)
        @test isa(constraint_object(c2), ScalarConstraint)
        @test constraint_object(c2).set == MOI.EqualTo(42.0)
        @test isa(set_normalized_rhs(LowerBoundRef(x), 42), Nothing)
        @test isa(constraint_object(LowerBoundRef(x)), ScalarConstraint)
        @test constraint_object(LowerBoundRef(x)).set == MOI.GreaterThan(42.0)
    end
    # test normalized_rhs
    @testset "JuMP.normalized_rhs" begin
        @test normalized_rhs(c1) == 42.
        @test normalized_rhs(c2) == 42.
        @test normalized_rhs(LowerBoundRef(x)) == 42.
    end
    # test add_to_function_constant
    @testset "JuMP.add_to_function_constant" begin
        @test isa(add_to_function_constant(c1, 4), Nothing)
        @test isa(constraint_object(c1), BoundedScalarConstraint)
        @test constraint_object(c1).set == MOI.EqualTo(38.0)
        @test isa(add_to_function_constant(c2, 4), Nothing)
        @test isa(constraint_object(c2), ScalarConstraint)
        @test constraint_object(c2).set == MOI.EqualTo(38.0)
        @test isa(add_to_function_constant(LowerBoundRef(x), 4), Nothing)
        @test isa(constraint_object(LowerBoundRef(x)), ScalarConstraint)
        @test constraint_object(LowerBoundRef(x)).set == MOI.GreaterThan(38.0)
    end
end

# Test num_constraints extensions
@testset "JuMP.num_constraints Extensions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par) >= 0)
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        @test num_constraints(m, GeneralVariableRef, MOI.LessThan) == 1
        @test num_constraints(m, GeneralVariableRef, MOI.GreaterThan) == 3
    end
    # test search function type
    @testset "Function" begin
        @test num_constraints(m, GeneralVariableRef) == 5
    end
    # test search set type
    @testset "Set" begin
        @test num_constraints(m, MOI.LessThan) == 1
        @test num_constraints(m, MOI.GreaterThan) == 3
    end
    # test search total
    @testset "Total" begin
        @test num_constraints(m) == 5
    end
end

# Test all_constraints extensions
@testset "JuMP.all_constraints Extensions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par) >= 0)
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        list = [InfOptConstraintRef(m, ConstraintIndex(4))]
        @test all_constraints(m, GeneralVariableRef, MOI.LessThan) == list
        list = [InfOptConstraintRef(m, ConstraintIndex(1)),
                InfOptConstraintRef(m, ConstraintIndex(2)),
                InfOptConstraintRef(m, ConstraintIndex(3))]
        @test all_constraints(m, GeneralVariableRef, MOI.GreaterThan) == list
    end
    # test search function type
    @testset "Function" begin
        list = [InfOptConstraintRef(m, ConstraintIndex(1)),
                InfOptConstraintRef(m, ConstraintIndex(2)),
                InfOptConstraintRef(m, ConstraintIndex(3)),
                InfOptConstraintRef(m, ConstraintIndex(4)),
                InfOptConstraintRef(m, ConstraintIndex(5))]
        @test all_constraints(m, GeneralVariableRef) == list
    end
    # test search set type
    @testset "Set" begin
        list = [InfOptConstraintRef(m, ConstraintIndex(4))]
        @test all_constraints(m, MOI.LessThan) == list
        list = [InfOptConstraintRef(m, ConstraintIndex(1)),
                InfOptConstraintRef(m, ConstraintIndex(2)),
                InfOptConstraintRef(m, ConstraintIndex(3))]
        @test all_constraints(m, MOI.GreaterThan) == list
    end
    # test search total
    @testset "Total" begin
        list = [InfOptConstraintRef(m, ConstraintIndex(1)),
                InfOptConstraintRef(m, ConstraintIndex(2)),
                InfOptConstraintRef(m, ConstraintIndex(3)),
                InfOptConstraintRef(m, ConstraintIndex(4)),
                InfOptConstraintRef(m, ConstraintIndex(5))]
        @test all_constraints(m) == list
    end
end

# Test list_of_constraint_types
@testset "JuMP.list_of_constraint_types" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par) >= 0)
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, 0 <= x <= 1, Int)
    expected = [(GeneralVariableRef, MOI.GreaterThan{Float64}),
                (GeneralVariableRef, MOI.LessThan{Float64}),
                (GeneralVariableRef, MOI.Integer)]
    @test isempty(setdiff(list_of_constraint_types(m), expected))
end
