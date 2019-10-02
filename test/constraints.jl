# Test basic extensions
@testset "Basic Extensions" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    cref1 = FiniteConstraintRef(m, 1, ScalarShape())
    cref2 = InfiniteConstraintRef(m, 2, ScalarShape())
    con = BoundedScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0),
                                  Dict{ParameterRef, IntervalSet}())
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
    @global_variable(m, x)
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
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
    # test _expand_parameter_dict(Dict{ParameterRef,IntervalSet}))
    @testset "_expand_parameter_dict (acceptable Form)" begin
        d = Dict(par => IntervalSet(0, 1))
        @test InfiniteOpt._expand_parameter_dict(d) == d
    end
    # test _expand_parameter_dict(Dict{Any,IntervalSet}))
    @testset "_expand_parameter_dict (Array Form)" begin
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(InfiniteOpt._expand_parameter_dict(d),
                  Dict{ParameterRef, IntervalSet})
    end
    # test _expand_parameter_dict(Dict))
    @testset "_expand_parameter_dict (Fallback)" begin
        d = Dict(pars => 1, par => 2)
        @test_throws ErrorException InfiniteOpt._expand_parameter_dict(d)
    end
    # test build_constraint (single)
    @testset "JuMP.build_constraint (Single)" begin
        # test bounded constraint
        con = BoundedScalarConstraint(inf, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 1)))
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
             parameter_bounds = Dict(par => IntervalSet(0, 1))).func == con.func
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
               parameter_bounds = Dict(par => IntervalSet(0, 1))).set == con.set
        @test build_constraint(error, inf, MOI.EqualTo(42.0),
         parameter_bounds = Dict(par => IntervalSet(0, 1))).bounds == con.bounds
        # test bounded constraint with vector bound key
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(build_constraint(error, inf, MOI.EqualTo(42.0),
                  parameter_bounds = d).bounds, Dict{ParameterRef, IntervalSet})
        # test scalar constraint
        con = ScalarConstraint(inf, MOI.EqualTo(42.0))
        @test build_constraint(error, inf, MOI.EqualTo(42.0)).func == con.func
        @test build_constraint(error, inf, MOI.EqualTo(42.0)).set == con.set
    end
    # test build_constraint (expr)
    @testset "JuMP.build_constraint (Expr)" begin
        # test bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 1)))
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
             parameter_bounds = Dict(par => IntervalSet(0, 1))).func == con.func
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
               parameter_bounds = Dict(par => IntervalSet(0, 1))).set == con.set
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0),
         parameter_bounds = Dict(par => IntervalSet(0, 1))).bounds == con.bounds
        # test bounded constraint with vector bound key
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(build_constraint(error, inf + pt, MOI.EqualTo(42.0),
                  parameter_bounds = d).bounds, Dict{ParameterRef, IntervalSet})
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
    # test _check_bounds
    @testset "_check_bounds" begin
        # test normal
        @test isa(InfiniteOpt._check_bounds(m, Dict(par => IntervalSet(0, 1))),
                                                                        Nothing)
        # test errors
        par2 = ParameterRef(InfiniteModel(), 2)
        @test_throws ErrorException InfiniteOpt._check_bounds(m,
                                                Dict(par2 => IntervalSet(0, 1)))
        @test_throws ErrorException InfiniteOpt._check_bounds(m,
                                                Dict(par => IntervalSet(-1, 1)))
        @test_throws ErrorException InfiniteOpt._check_bounds(m,
                                                 Dict(par => IntervalSet(0, 2)))
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
        @global_variable(InfiniteModel(), g)
        con = ScalarConstraint(g, MOI.EqualTo(42.0))
        @test_throws VariableNotOwned add_constraint(m, con, "a")
        # test bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 1)))
        @test add_constraint(m, con, "a") == InfiniteConstraintRef(m, 1,
                                                                  ScalarShape())
        # test bad bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 2)))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test infinite constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "b") == InfiniteConstraintRef(m, 2,
                                                                  ScalarShape())
        # test measure constraint
        con = ScalarConstraint(x + meas, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "c") == MeasureConstraintRef(m, 3,
                                                                  ScalarShape())
        # test finite constraint
        con = ScalarConstraint(x + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "d") == FiniteConstraintRef(m, 4,
                                                                  ScalarShape())
        @test name(FiniteConstraintRef(m, 4, ScalarShape())) == "d"
        @test !m.constr_in_var_info[4]
        @test !optimizer_model_ready(m)
    end
    # test macro
    @testset "JuMP.@constraint" begin
        # test scalar constraint
        @test @constraint(m, e, x + pt -2 <= 2) == FiniteConstraintRef(m, 5,
                                                                  ScalarShape())
        @test isa(m.constrs[5], ScalarConstraint)
        # test bounded scalar constraint
        @test @constraint(m, f, inf + meas <= 2,
                          parameter_bounds = Dict(par => IntervalSet(0, 1))) == InfiniteConstraintRef(m, 6, ScalarShape())
        @test isa(m.constrs[6], BoundedScalarConstraint)
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
    # test _make_interval_set
    @testset "_make_interval_set" begin
        expected = :(IntervalSet(0, 1))
        @test InfiniteOpt._make_interval_set(error, :([0, 1])) == expected
        @test_throws ErrorException InfiniteOpt._make_interval_set(error, :([0]))
    end
    # test _process_parameter
    @testset "_process_parameter" begin
        @test InfiniteOpt._process_parameter(:x) == :(InfiniteOpt._make_vector(x))
    end
    # test _make_bound_pair (in)
    @testset "_make_bound_pair (in)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
    end
    # test _make_bound_pair (==)
    @testset "_make_bound_pair (==)" begin
        expr = :(t == 0)
        set = :(IntervalSet(0, 0))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:(==))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with <=)
    @testset "_make_bound_pair (<=)" begin
        expr = :(0 <= t <= 1)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with >=)
    @testset "_make_bound_pair (>=)" begin
        expr = :(1 >= t >= 0)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(>=)),
                                           Val(:(>=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison fallback)
    @testset "_make_bound_pair (Fallback)" begin
        expr = :(t <= 0)
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                                 Val(:(<=)))
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                         Val(:(<=)), Val(:(>=)))
    end
    # test _make_bound_pair (wrapper)
    @testset "_make_bound_pair (Wrapper)" begin
        # test :call expression
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
        # test :compariosn expression
        expr = :(0 <= t <= 1)
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
        # test error
        expr = :([2, 3])
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr)
    end
    # test _parse_parameter_bounds (Vector)
    @testset "_parse_parameter_bounds (Vector)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg), $(dict_arg)))
        @test InfiniteOpt._parse_parameter_bounds(error, args) == dict
    end
    # test _parse_parameter_bounds (Expr)
    @testset "_parse_parameter_bounds (Expr)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg)))
        @test InfiniteOpt._parse_parameter_bounds(error, expr) == dict
    end
    # test _extract_bounds (:call)
    @testset "_extract_bounds (:call)" begin
        # test single anonymous bound
        args = [:in, :t, :([0, 1])]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(Dict($(dict_arg)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (nothing, bounds)
        # test with name expression
        args = [:(bob[i = 1:2]), :(t in [0, 1])]
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (args[1], bounds)
    end
    # test _extract_bounds (:tuple)
    @testset "_extract_bounds (:tuple)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(Dict($(dict_arg), $(dict_arg)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:tuple)) == (nothing, bounds)
    end
    # test _extract_bounds (:comparison)
    @testset "_extract_bounds (:comparison)" begin
        args = [0, :(<=), :t, :(<=), 1]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(Dict($(dict_arg)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:comparison)) == (nothing, bounds)
    end
    # test @BDconstraint
    @testset "@BDconstraint" begin
        # initialize model and references
        m = InfiniteModel()
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        @infinite_variable(m, inf(par))
        @point_variable(m, inf(0.5), pt)
        @global_variable(m, x)
        bounds1 = Dict(par => IntervalSet(0, 1))
        m2 = Model()
        # test errors
        @test_macro_throws ErrorException @BDconstraint(m, par = 0, inf == 0,
                                                     parameter_bounds = bounds1)
        @test_macro_throws ErrorException @BDconstraint(m, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con[1:2], inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, [1:2], inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, con, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m, par = 2, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m2, par == 0, inf == 0)
        @test_macro_throws ErrorException @BDconstraint(m2, con(par == 0), inf == 0)
        # test anonymous constraint with set
        cref = InfiniteConstraintRef(m, 1, ScalarShape())
        @test @BDconstraint(m, par in [0, 1], inf + x == 2) == cref
        @test m.constrs[1].bounds == bounds1
        # test anonymous constraint with comparison
        cref = InfiniteConstraintRef(m, 2, ScalarShape())
        @test @BDconstraint(m, 0 <= par <= 1, inf + x == 2) == cref
        @test m.constrs[2].bounds == bounds1
        # test anonymous constraint with equality
        cref = InfiniteConstraintRef(m, 3, ScalarShape())
        @test @BDconstraint(m, par == 0, inf + x == 2) == cref
        @test m.constrs[3].bounds == Dict(par => IntervalSet(0, 0))
        # test anonymous multiple constraint
        cref1 = InfiniteConstraintRef(m, 4, ScalarShape())
        cref2 = InfiniteConstraintRef(m, 5, ScalarShape())
        @test @BDconstraint(m, [1:2](par == 0), inf + x == 2) == [cref1, cref2]
        @test m.constrs[4].bounds == Dict(par => IntervalSet(0, 0))
        @test m.constrs[5].bounds == Dict(par => IntervalSet(0, 0))
        # test anonymous multiple bounds
        cref = InfiniteConstraintRef(m, 6, ScalarShape())
        @test @BDconstraint(m, (par == 0, pars[1] in [0, 1]), inf + x == 2) == cref
        @test m.constrs[6].bounds == Dict(par => IntervalSet(0, 0), pars[1] => IntervalSet(0, 1))
        # test anonymous multiple bounds with vector input
        cref = InfiniteConstraintRef(m, 7, ScalarShape())
        @test @BDconstraint(m, (par == 0, pars in [0, 1]), inf + x == 2) == cref
        @test m.constrs[7].bounds == Dict(par => IntervalSet(0, 0), pars[1] => IntervalSet(0, 1),
                                          pars[2] => IntervalSet(0, 1))
        # test named constraint
        cref = InfiniteConstraintRef(m, 8, ScalarShape())
        @test @BDconstraint(m, c1(par in [0, 1]), inf + x == 2) == cref
        @test m.constrs[8].bounds == bounds1
        # test named constraint with other vector type
        cref = InfiniteConstraintRef(m, 9, ScalarShape())
        @test @BDconstraint(m, c2(par in [0, 1], 0 <= pars <= 1), inf + x == 2) == cref
        @test m.constrs[9].bounds == Dict(par => IntervalSet(0, 1), pars[1] => IntervalSet(0, 1),
                                          pars[2] => IntervalSet(0, 1))
        # test named constraints
        cref1 = InfiniteConstraintRef(m, 10, ScalarShape())
        cref2 = InfiniteConstraintRef(m, 11, ScalarShape())
        @test @BDconstraint(m, c3[1:2](par == 0), inf + x == 2) == [cref1, cref2]
        @test m.constrs[10].bounds == Dict(par => IntervalSet(0, 0))
        @test m.constrs[11].bounds == Dict(par => IntervalSet(0, 0))
        # test different container type
        cref1 = InfiniteConstraintRef(m, 12, ScalarShape())
        cref2 = InfiniteConstraintRef(m, 13, ScalarShape())
        expected = convert(JuMPC.SparseAxisArray, [cref1, cref2])
        @test @BDconstraint(m, c4[1:2](par == 0), inf + x == 2,
                            container = SparseAxisArray) == expected
        @test m.constrs[12].bounds == Dict(par => IntervalSet(0, 0))
        @test m.constrs[13].bounds == Dict(par => IntervalSet(0, 0))
        # test with SparseAxisArray input
        cref = InfiniteConstraintRef(m, 14, ScalarShape())
        pars = convert(JuMPC.SparseAxisArray, pars)
        @test @BDconstraint(m, c5(pars in [0, 1]), inf + x == 2) == cref
        @test m.constrs[14].bounds == Dict(pars[1] => IntervalSet(0, 1),
                                           pars[2] => IntervalSet(0, 1))
    end
end

# Test coefficient queries and modifications
@testset "Coefficient Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x >= 0, Int)
    @constraint(m, c1, inf + x == 0,
                parameter_bounds = Dict(par => IntervalSet(0, 1)))
    @constraint(m, c2, x * pt + x == 2)
    # test set_normalized_coefficient
    @testset "JuMP.set_normalized_coefficient" begin
        @test isa(set_normalized_coefficient(c1, x, 2), Nothing)
        @test isa(constraint_object(c1), BoundedScalarConstraint)
        @test jump_function(constraint_object(c1)) == inf + 2x
        @test isa(set_normalized_coefficient(c2, inf, 2), Nothing)
        @test isa(constraint_object(c2),ScalarConstraint)
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
    @global_variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        @test num_constraints(m, GlobalVariableRef, MOI.LessThan) == 1
        @test num_constraints(m, InfiniteVariableRef, MOI.GreaterThan) == 1
    end
    # test search function type
    @testset "Function" begin
        @test num_constraints(m, GlobalVariableRef) == 3
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
    @global_variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        list = [FiniteConstraintRef(m, 4, ScalarShape())]
        @test all_constraints(m, GlobalVariableRef, MOI.LessThan) == list
        list = [InfiniteConstraintRef(m, 1, ScalarShape())]
        @test all_constraints(m, InfiniteVariableRef, MOI.GreaterThan) == list
    end
    # test search function type
    @testset "Function" begin
        list = [FiniteConstraintRef(m, 3, ScalarShape()),
                FiniteConstraintRef(m, 4, ScalarShape()),
                FiniteConstraintRef(m, 5, ScalarShape())]
        @test all_constraints(m, GlobalVariableRef) == list
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
    @global_variable(m, 0 <= x <= 1, Int)
    expected = [(InfiniteVariableRef, MOI.GreaterThan{Float64}),
                (PointVariableRef, MOI.GreaterThan{Float64}),
                (GlobalVariableRef, MOI.GreaterThan{Float64}),
                (GlobalVariableRef, MOI.LessThan{Float64}),
                (GlobalVariableRef, MOI.Integer)]
    @test list_of_constraint_types(m) == expected
end
