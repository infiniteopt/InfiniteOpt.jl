# Test basic extensions
@testset "Basic Extensions" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    cref1 = FiniteConstraintRef(m, 1, ScalarShape())
    cref2 = InfiniteConstraintRef(m, 2, ScalarShape())
    con = BoundedScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0),
                                  ParameterBounds(), ParameterBounds())
    # owner_model
    @testset "JuMP.owner_model" begin
      @test owner_model(cref1) == m
    end
    # index
    @testset "JuMP.index" begin
      @test JuMP.index(cref1) == 1
    end
    # is_valid
    @testset "JuMP.is_valid" begin
      @test !is_valid(m, cref1)
      m.constrs[1] = con
      @test is_valid(m, cref1)
    end
    # constraint_object
    @testset "JuMP.constraint_object" begin
      @test constraint_object(cref1).func == con.func
      @test constraint_object(cref1).set == con.set
      @test constraint_object(cref1).bounds == con.bounds
      @test constraint_object(cref1).orig_bounds == con.orig_bounds
    end
    # test Base.:(==) of constraint references
    @testset "Base.:(==) References" begin
        @test !(cref1 == cref2)
        @test cref1 == FiniteConstraintRef(m, 1, ScalarShape())
        @test cref1 != FiniteConstraintRef(m, 2, ScalarShape())
        @test cref1 != FiniteConstraintRef(InfiniteModel(), 1, ScalarShape())
    end
    # test broadcastable
    @testset "Base.broadcastable Reference" begin
        @test isa(Base.broadcastable(cref1), Base.RefValue)
    end
    # test constraint type
    @testset "JuMP.constraint_type" begin
        @test constraint_type(m) == GeneralConstraintRef
    end
    # test shape
    @testset "shape" begin
        @test JuMP.shape(con) == ScalarShape()
    end
    # test jump_function
    @testset "JuMP.jump_function" begin
        @test jump_function(con) == zero(AffExpr)
    end
    # test moi_set
    @testset "JuMP.moi_set" begin
        @test moi_set(con) == MOI.EqualTo(0.0)
    end
end

# Test name methods
@testset "Name" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    cref = FiniteConstraintRef(m, 1, ScalarShape())
    m.constr_to_name[1] = "test"
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(cref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(cref, "new"), Nothing)
        @test name(cref) == "new"
    end
    # _make_constraint_ref
    @testset "_make_constraint_ref" begin
        # prepare for infinite constraint
        m.constrs[1] = ScalarConstraint(inf + par, MOI.LessThan(1.0))
        @test InfiniteOpt._make_constraint_ref(m, 1) == InfiniteConstraintRef(m,
                                                        1, ScalarShape())
        # prepare for measure constraint
        m.constrs[1] = ScalarConstraint(pt + meas, MOI.LessThan(1.0))
        @test InfiniteOpt._make_constraint_ref(m, 1) == MeasureConstraintRef(m,
                                                        1, ScalarShape())
        # prepare for finite constraint
        m.constrs[1] = ScalarConstraint(pt + x, MOI.LessThan(1.0))
        @test InfiniteOpt._make_constraint_ref(m, 1) == FiniteConstraintRef(m,
                                                        1, ScalarShape())
    end
    # constraint_by_name
    @testset "JuMP.constraint_by_name" begin
        # test normal
        @test constraint_by_name(m, "new") == cref
        @test isa(constraint_by_name(m, "test2"), Nothing)
        # prepare constraint with duplicate name
        m.constrs[2] = ScalarConstraint(inf + par, MOI.LessThan(1.0))
        m.constr_to_name[2] = "new"
        m.name_to_constr = nothing
        # test for duplciate name error
        @test_throws ErrorException constraint_by_name(m, "new")
    end
end

# Test definition methods
@testset "Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @hold_variable(m, x, parameter_bounds = (par == 1))
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
    # test build_constraint (single)
    @testset "JuMP.build_constraint (Single)" begin
        # test bounded constraint
        con = BoundedScalarConstraint(inf, MOI.EqualTo(42.0),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))))
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
             parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).func == con.func
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
               parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).set == con.set
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
         parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).bounds == con.bounds
        # test scalar constraint
        con = ScalarConstraint(inf, MOI.EqualTo(42.0))
        @test build_constraint(error, inf, MOI.EqualTo(42.0)).func == con.func
        @test build_constraint(error, inf, MOI.EqualTo(42.0)).set == con.set
    end
    # test build_constraint (expr)
    @testset "JuMP.build_constraint (Expr)" begin
        # test bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))))
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
             parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).func == con.func
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
               parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).set == con.set
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
         parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))).bounds == con.bounds
        # test scalar constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0)).func == con.func
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0)).set == con.set
    end
    # test _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # test initial use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], 1),
                  Nothing)
        @test m.var_to_constrs[JuMP.index(inf)] == [1]
        @test m.param_to_constrs[JuMP.index(par)] == [1]
        @test m.meas_to_constrs[JuMP.index(meas)] == [1]
        @test m.reduced_to_constrs[JuMP.index(rv)] == [1]
        # test secondary use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], 2),
                  Nothing)
        @test m.var_to_constrs[JuMP.index(inf)] == [1, 2]
        @test m.param_to_constrs[JuMP.index(par)] == [1, 2]
        @test m.meas_to_constrs[JuMP.index(meas)] == [1, 2]
        @test m.reduced_to_constrs[JuMP.index(rv)] == [1, 2]
        # Undo changes
        delete!(m.var_to_constrs, JuMP.index(inf))
        delete!(m.param_to_constrs, JuMP.index(par))
        delete!(m.meas_to_constrs, JuMP.index(meas))
        delete!(m.reduced_to_constrs, JuMP.index(rv))
    end
    # test _check_and_update_bounds (BoundedScalarConstraint with no hold variables)
    @testset "_check_and_update_bounds (Bounded no Hold)" begin
        c = BoundedScalarConstraint(inf, MOI.EqualTo(0.0), ParameterBounds(),
                                    ParameterBounds())
        @test InfiniteOpt._check_and_update_bounds(m, c, [inf]).bounds == c.bounds
    end
    # test _check_and_update_bounds (BoundedScalarConstraint with hold variables)
    @testset "_check_and_update_bounds (Bounded with Hold)" begin
        bounds = ParameterBounds()
        c = BoundedScalarConstraint(x + inf, MOI.EqualTo(0.0), bounds, copy(bounds))
        @test InfiniteOpt._check_and_update_bounds(m, c,
                            [x, inf]).bounds.intervals[par] == IntervalSet(1, 1)
        @test InfiniteOpt._check_and_update_bounds(m, c,
                                      [x, inf]).orig_bounds == ParameterBounds()
    end
    # test _check_and_update_bounds (ScalarConstraint with no hold variables)
    @testset "_check_and_update_bounds (Scalar no Hold)" begin
        c = JuMP.ScalarConstraint(inf, MOI.EqualTo(0.0))
        @test InfiniteOpt._check_and_update_bounds(m, c, [inf]) == c
    end
    # test _check_and_update_bounds (ScalarConstraint with hold variables)
    @testset "_check_and_update_bounds (Scalar with Hold)" begin
        c = JuMP.ScalarConstraint(inf + x, MOI.EqualTo(0.0))
        @test InfiniteOpt._check_and_update_bounds(m, c,
                            [x, inf]).bounds.intervals[par] == IntervalSet(1, 1)
        set_parameter_bounds(x, ParameterBounds(), force = true)
        m.has_hold_bounds = false
    end
    # test add_constraint
    @testset "JuMP.add_constraint" begin
        # test reject vector constraint
        con = VectorConstraint([inf + pt], MOI.Zeros(1))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test rejects parameter only input
        con = ScalarConstraint(par, MOI.EqualTo(42.0))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test rejects bad variables
        @variable(Model(), z)
        con = ScalarConstraint(z, MOI.EqualTo(42.0))
        @test_throws VariableNotOwned add_constraint(m, con, "a")
        @hold_variable(InfiniteModel(), g)
        con = ScalarConstraint(g, MOI.EqualTo(42.0))
        @test_throws VariableNotOwned add_constraint(m, con, "a")
        # test bounded constraint
        index = m.next_constr_index + 1
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))),
                                ParameterBounds(Dict(par => IntervalSet(0, 1))))
        @test add_constraint(m, con, "a") == InfiniteConstraintRef(m, index,
                                                                  ScalarShape())
        # test infinite constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "b") == InfiniteConstraintRef(m, index + 1,
                                                                  ScalarShape())
        # test measure constraint
        con = ScalarConstraint(x + meas, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "c") == MeasureConstraintRef(m, index + 2,
                                                                  ScalarShape())
        # test finite constraint
        con = ScalarConstraint(x + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "d") == FiniteConstraintRef(m, index + 3,
                                                                  ScalarShape())
        @test name(FiniteConstraintRef(m, index + 3, ScalarShape())) == "d"
        @test !m.constr_in_var_info[index + 3]
        @test !optimizer_model_ready(m)
        # test with bounded hold variables
        @set_parameter_bounds(x, par == 1)
        con = ScalarConstraint(inf + pt + x, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "b") == InfiniteConstraintRef(m, index + 4,
                                                                  ScalarShape())
        @test m.constrs[index + 4].bounds.intervals[par] == IntervalSet(1, 1)
        @test m.constrs[index + 4].orig_bounds == ParameterBounds()
        set_parameter_bounds(x, ParameterBounds(), force = true)
    end
    # test macro
    @testset "JuMP.@constraint" begin
        # test scalar constraint
        index = m.next_constr_index + 1
        @test @constraint(m, e, x + pt -2 <= 2) == FiniteConstraintRef(m, index,
                                                                  ScalarShape())
        @test isa(m.constrs[index], ScalarConstraint)
        # test bounded scalar constraint
        @test @constraint(m, f, inf + meas <= 2,
            parameter_bounds = ParameterBounds(Dict(par => IntervalSet(0, 1)))) == InfiniteConstraintRef(m, index + 1, ScalarShape())
        @test isa(m.constrs[index + 1], BoundedScalarConstraint)
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
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg))))
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
        index = m.next_constr_index + 1
        cref = InfiniteConstraintRef(m, index, ScalarShape())
        @test @BDconstraint(m, par in [0, 1], inf + x == 2) == cref
        @test m.constrs[1].bounds == bounds1
        # test anonymous constraint with comparison
        cref = InfiniteConstraintRef(m, index + 1, ScalarShape())
        @test @BDconstraint(m, 0 <= par <= 1, inf + x == 2) == cref
        @test m.constrs[2].bounds == bounds1
        # test anonymous constraint with equality
        cref = InfiniteConstraintRef(m, index + 2, ScalarShape())
        @test @BDconstraint(m, par == 0, inf + x == 2) == cref
        @test m.constrs[3].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test anonymous multiple constraint
        cref1 = InfiniteConstraintRef(m, index + 3, ScalarShape())
        cref2 = InfiniteConstraintRef(m, index + 4, ScalarShape())
        @test @BDconstraint(m, [1:2](par == 0), inf + x == 2) == [cref1, cref2]
        @test m.constrs[4].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test m.constrs[5].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test anonymous multiple bounds
        cref = InfiniteConstraintRef(m, index + 5, ScalarShape())
        @test @BDconstraint(m, (par == 0, pars[1] in [0, 1]), inf + x == 2) == cref
        @test m.constrs[6].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0),
                                                  pars[1] => IntervalSet(0, 1)))
        # test anonymous multiple bounds with vector input
        cref = InfiniteConstraintRef(m, index + 6, ScalarShape())
        @test @BDconstraint(m, (par == 0, pars in [0, 1]), inf + x == 2) == cref
        @test m.constrs[7].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0),
                    pars[1] => IntervalSet(0, 1), pars[2] => IntervalSet(0, 1)))
        # test named constraint
        cref = InfiniteConstraintRef(m, index + 7, ScalarShape())
        @test @BDconstraint(m, c1(par in [0, 1]), inf + x == 2) == cref
        @test m.constrs[8].bounds == bounds1
        # test named constraint with other vector type
        cref = InfiniteConstraintRef(m, index + 8, ScalarShape())
        @test @BDconstraint(m, c2(par in [0, 1], 0 <= pars <= 1), inf + x == 2) == cref
        @test m.constrs[9].bounds == ParameterBounds(Dict(par => IntervalSet(0, 1),
                    pars[1] => IntervalSet(0, 1), pars[2] => IntervalSet(0, 1)))
        # test named constraints
        cref1 = InfiniteConstraintRef(m, index + 9, ScalarShape())
        cref2 = InfiniteConstraintRef(m, index + 10, ScalarShape())
        @test @BDconstraint(m, c3[1:2](par == 0), inf + x == 2) == [cref1, cref2]
        @test m.constrs[10].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test m.constrs[11].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test different container type
        cref1 = InfiniteConstraintRef(m, index + 11, ScalarShape())
        cref2 = InfiniteConstraintRef(m, index + 12, ScalarShape())
        expected = convert(JuMPC.SparseAxisArray, [cref1, cref2])
        @test @BDconstraint(m, c4[1:2](par == 0), inf + x == 2,
                            container = SparseAxisArray) == expected
        @test m.constrs[12].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test m.constrs[13].bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        # test with SparseAxisArray input
        cref = InfiniteConstraintRef(m, index + 13, ScalarShape())
        pars = convert(JuMPC.SparseAxisArray, pars)
        @test @BDconstraint(m, c5(pars in [0, 1]), inf + x == 2) == cref
        @test m.constrs[14].bounds == ParameterBounds(Dict(pars[1] => IntervalSet(0, 1),
                                                  pars[2] => IntervalSet(0, 1)))
        # test with hold variable bounds
        @set_parameter_bounds(x, pars[1] in [0, 2])
        cref = InfiniteConstraintRef(m, index + 14, ScalarShape())
        @test @BDconstraint(m, c10(par in [0, 1]), inf + x == 2) == cref
        @test m.constrs[15].bounds == ParameterBounds(Dict(par => IntervalSet(0, 1),
                                                  pars[1] => IntervalSet(0, 2)))
        @test m.constrs[15].orig_bounds == ParameterBounds(Dict(par => IntervalSet(0, 1)))
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
    # test has_parameter_bounds
    @testset "has_parameter_bounds" begin
        @test has_parameter_bounds(c1)
        @test !has_parameter_bounds(c2)
    end
    # test parameter_bounds
    @testset "parameter_bounds" begin
        @test_throws ErrorException parameter_bounds(c2)
        @test parameter_bounds(c1) == ParameterBounds(Dict(par => IntervalSet(0, 1)))
    end
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
    @testset "add_parameter_bound" begin
        # test already has bounds
        @test isa(add_parameter_bound(c1, par, 0, 1), Nothing)
        @test parameter_bounds(c1) == ParameterBounds(Dict(par => IntervalSet(0, 1)))
        # test doesn't have bounds
        @constraint(m, c3, inf + pt == 0)
        @test isa(add_parameter_bound(c3, par, 0, 1), Nothing)
        @test parameter_bounds(c3) == ParameterBounds(Dict(par => IntervalSet(0, 1)))
    end
    # test delete_parameter_bound
    @testset "delete_parameter_bound" begin
        # test partial deletion
        @test isa(add_parameter_bound(c1, pars[1], 0, 1), Nothing)
        @add_parameter_bounds(x, par == 0)
        @test isa(delete_parameter_bound(c1, par), Nothing)
        expected = ParameterBounds(Dict(par => IntervalSet(0, 0),
                                   pars[1] => IntervalSet(0, 1)))
        @test parameter_bounds(c1) == expected
        # test again without hold bounds
        delete_parameter_bounds(x)
        expected = ParameterBounds(Dict(pars[1] => IntervalSet(0, 1)))
        @test parameter_bounds(c1) == expected
        # test already gone
        @test isa(delete_parameter_bound(c1, par), Nothing)
        @test parameter_bounds(c1) == expected
    end
    # test delete_parameter_bounds
    @testset "delete_parameter_bounds" begin
        # test partial deletion
        @test isa(add_parameter_bound(c1, pars[1], 0, 1), Nothing)
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
    # test
    # test @set_parameter_bounds
    @testset "@set_parameter_bounds" begin
        # Note errors were already checked with hold variables
        # test force error
        @test_macro_throws ErrorException @set_parameter_bounds(c1, par == 1)
        # test with single
        @set_parameter_bounds(x, pars[2] == 0)
        @test isa(@set_parameter_bounds(c1, par == 1, force = true), Nothing)
        @test parameter_bounds(c1) == ParameterBounds(Dict(par => IntervalSet(1, 1),
                                                  pars[2] => IntervalSet(0, 0)))
        set_parameter_bounds(x, ParameterBounds(), force = true)
        # test multiple
        @test isa(@set_parameter_bounds(c1, pars == 0, force = true), Nothing)
        @test parameter_bounds(c1).intervals[pars[2]] == IntervalSet(0, 0)
    end
    # test @add_parameter_bounds
    @testset "@add_parameter_bounds" begin
        # Note errors were already checked with hold variables
        # test bounds error
        @test_macro_throws ErrorException @add_parameter_bounds(c3, par in [-1, 1])
        # test with multiple
        @test isa(@add_parameter_bounds(c2, (pars == 0, par in [0, 1])), Nothing)
        @test parameter_bounds(c2).intervals[pars[2]] == IntervalSet(0, 0)
        @test parameter_bounds(c2).intervals[par] == IntervalSet(0, 1)
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
        @test num_constraints(m, HoldVariableRef, MOI.LessThan) == 1
        @test num_constraints(m, InfiniteVariableRef, MOI.GreaterThan) == 1
    end
    # test search function type
    @testset "Function" begin
        @test num_constraints(m, HoldVariableRef) == 3
        @test num_constraints(m, InfiniteVariableRef) == 1
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
        list = [FiniteConstraintRef(m, 4, ScalarShape())]
        @test all_constraints(m, HoldVariableRef, MOI.LessThan) == list
        list = [InfiniteConstraintRef(m, 1, ScalarShape())]
        @test all_constraints(m, InfiniteVariableRef, MOI.GreaterThan) == list
    end
    # test search function type
    @testset "Function" begin
        list = [FiniteConstraintRef(m, 3, ScalarShape()),
                FiniteConstraintRef(m, 4, ScalarShape()),
                FiniteConstraintRef(m, 5, ScalarShape())]
        @test all_constraints(m, HoldVariableRef) == list
        list = [InfiniteConstraintRef(m, 1, ScalarShape())]
        @test all_constraints(m, InfiniteVariableRef) == list
    end
    # test search set type
    @testset "Set" begin
        list = [FiniteConstraintRef(m, 4, ScalarShape())]
        @test all_constraints(m, MOI.LessThan) == list
        list = [InfiniteConstraintRef(m, 1, ScalarShape()),
                FiniteConstraintRef(m, 2, ScalarShape()),
                FiniteConstraintRef(m, 3, ScalarShape())]
        @test all_constraints(m, MOI.GreaterThan) == list
    end
    # test search total
    @testset "Total" begin
        list = [InfiniteConstraintRef(m, 1, ScalarShape()),
                FiniteConstraintRef(m, 2, ScalarShape()),
                FiniteConstraintRef(m, 3, ScalarShape()),
                FiniteConstraintRef(m, 4, ScalarShape()),
                FiniteConstraintRef(m, 5, ScalarShape())]
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
    expected = [(InfiniteVariableRef, MOI.GreaterThan{Float64}),
                (PointVariableRef, MOI.GreaterThan{Float64}),
                (HoldVariableRef, MOI.GreaterThan{Float64}),
                (HoldVariableRef, MOI.LessThan{Float64}),
                (HoldVariableRef, MOI.Integer)]
    @test list_of_constraint_types(m) == expected
end
