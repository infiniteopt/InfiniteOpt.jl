# Test basic queries
@testset "Basic Queries" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    cref = FiniteConstraintRef(m, 1, ScalarShape())
    # owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(cref) == m
    end
    # index
    @testset "JuMP.index" begin
        @test index(cref) == 1
    end
end

# Test basic extensions
@testset "Basic Extensions" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    cref1 = FiniteConstraintRef(m, 1, ScalarShape())
    cref2 = InfiniteConstraintRef(m, 2, ScalarShape())
    con = BoundedScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0),
                                  Dict{ParameterRef, IntervalSet}())
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
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
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
        # test scalar constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0)).func == con.func
        @test build_constraint(error, inf + pt, MOI.EqualTo(42.0)).set == con.set
    end
    # test _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # test initial use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas], 1),
                  Nothing)
        @test m.var_to_constrs[index(inf)] == [1]
        @test m.param_to_constrs[index(par)] == [1]
        @test m.meas_to_constrs[index(meas)] == [1]
        # test secondary use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas], 2),
                  Nothing)
        @test m.var_to_constrs[index(inf)] == [1, 2]
        @test m.param_to_constrs[index(par)] == [1, 2]
        @test m.meas_to_constrs[index(meas)] == [1, 2]
        # Undo changes
        delete!(m.var_to_constrs, index(inf))
        delete!(m.param_to_constrs, index(par))
        delete!(m.meas_to_constrs, index(meas))
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
        @test_throws ErrorException add_constraint(m, con, "a")
        @global_variable(InfiniteModel(), g)
        con = ScalarConstraint(g, MOI.EqualTo(42.0))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 1)))
        @test add_constraint(m, con, "a") == InfiniteConstraintRef(m, 1, ScalarShape())
        # test bad bounded constraint
        con = BoundedScalarConstraint(inf + pt, MOI.EqualTo(42.0),
                                       Dict(par => IntervalSet(0, 2)))
        @test_throws ErrorException add_constraint(m, con, "a")
        # test infinite constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "b") == InfiniteConstraintRef(m, 2, ScalarShape())
        # test measure constraint
        con = ScalarConstraint(x + meas, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "c") == MeasureConstraintRef(m, 3, ScalarShape())
        # test finite constraint
        con = ScalarConstraint(x + pt, MOI.EqualTo(42.0))
        @test add_constraint(m, con, "d") == FiniteConstraintRef(m, 4, ScalarShape())
        @test name(FiniteConstraintRef(m, 4, ScalarShape())) == "d"
        @test !m.constr_in_var_info[4]
        @test !optimizer_model_ready(m)
    end
end
