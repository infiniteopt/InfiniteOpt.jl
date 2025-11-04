# Test variable initializers
@testset "Variable Initializers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], derivative_method = OrthogonalCollocation(2))
    @variable(m, x >= cos, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = (p, ps) -> p + sum(ps))
    x0 = x(0)
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]))
    y1 = y(1, [1, 1])
    unfix(y1)
    y01 = y(0, [1, 1])
    delete_lower_bound(y01)
    delete_upper_bound(y01)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    @parameter_function(m, pf == (par, pars) -> 42)
    point_pf = pf(0, [0, 0])
    semi_pf = pf(0, pars)
    dx = deriv(x, par)
    dy = deriv(y, par)
    dx3 = @deriv(x, par^3)
    set_upper_bound(dx, sin)
    set_start_value(dy, (p, ps) -> p + sum(ps))
    yrv = y(par, [0, 0])
    @variable(m, 0 <= yrv2 <= 3, SemiInfinite(y, par, [1, 1]), start = 1)
    @variable(m, yrv3 <= 3, SemiInfinite(y, 1, pars), start = 1)
    @variable(m, yrv4 == 4, SemiInfinite(y, 0, pars))
    data = DiscreteMeasureData(par, [0.5, 0.5], [0, 1])
    @constraint(m, c1, x + z - 2 <= 0)
    @constraint(m, c2, measure(x + y, data) - w == 0)
    @constraint(m, c3, x0 + y0 == 5)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    # test transcribe_finite_variables!
    @testset "transcribe_finite_variables!" begin
        @test isa(IOTO.transcribe_finite_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).finvar_mappings) == 2
        @test IOTO.transcription_variable(z, tb) isa VariableRef
        @test IOTO.transcription_variable(w, tb) isa VariableRef
        @test name(IOTO.transcription_variable(z, tb)) == name(z)
        @test name(IOTO.transcription_variable(w, tb)) == name(w)
        @test has_lower_bound(IOTO.transcription_variable(z, tb))
        @test has_upper_bound(IOTO.transcription_variable(z, tb))
        @test is_binary(IOTO.transcription_variable(z, tb))
        @test is_fixed(IOTO.transcription_variable(w, tb))
        @test is_integer(IOTO.transcription_variable(w, tb))
        @test start_value(IOTO.transcription_variable(w, tb)) == 1.
    end
    # test _format_infinite_info
    @testset "_format_infinite_info" begin
        info = core_object(y).info
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_lb
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_ub
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_fix
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).fixed_value == 2
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_start
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).start == 1.25
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).binary 
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).integer
        info = core_object(x).info
        @test IOTO._format_infinite_info(info, [0.]).has_lb
        @test IOTO._format_infinite_info(info, [0.]).lower_bound == cos(0)
        @test !IOTO._format_infinite_info(info, [0.]).has_ub
        @test !IOTO._format_infinite_info(info, [0.]).has_fix
        @test !IOTO._format_infinite_info(info, [0.]).has_start
        @test isnan(IOTO._format_infinite_info(info, [0.]).start)
        @test IOTO._format_infinite_info(info, [0.]).integer
        info = core_object(dy).info
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_lb
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_ub
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_fix
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).has_start
        @test IOTO._format_infinite_info(info, [0., 0.5, 0.75]).start == 1.25
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).binary 
        @test !IOTO._format_infinite_info(info, [0., 0.5, 0.75]).integer
        info = core_object(dx).info
        @test IOTO._format_infinite_info(info, [0.]).has_ub
        @test IOTO._format_infinite_info(info, [0.]).upper_bound == sin(0)
        @test !IOTO._format_infinite_info(info, [0.]).has_lb
        @test !IOTO._format_infinite_info(info, [0.]).has_fix
        @test !IOTO._format_infinite_info(info, [0.]).has_start
        @test isnan(IOTO._format_infinite_info(info, [0.]).start)
        @test !IOTO._format_infinite_info(info, [0.]).integer
    end
    # test name formatting
    @testset "_make_var_name" begin
        @test IOTO._make_var_name("x", [1, 2, 3], (0, [1, 1]), (1, 2)) == "x(0, [1, 1])"
        @test IOTO._make_var_name("x", [1, 2, 3, 4, 5], (0, [0, 0, 0, 0]), (1, 2)) == "x[1, 2]"
    end
    # test transcribe_infinite_variables!
    @testset "transcribe_infinite_variables!" begin
        @test isa(IOTO.transcribe_infinite_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 2
        @test IOTO.transcription_variable(x, tb) isa Vector{VariableRef}
        @test IOTO.transcription_variable(y, tb) isa Matrix{VariableRef}
        @test name(IOTO.transcription_variable(x, tb)[1]) == "x(0.0)"
        @test name(IOTO.transcription_variable(y, tb)[2, 1]) == "y(1.0, [0.0, 0.0])"
        @test has_lower_bound(IOTO.transcription_variable(x)[1])
        @test is_binary(IOTO.transcription_variable(y, tb)[2])
        @test is_fixed(IOTO.transcription_variable(y, tb)[4])
        @test is_integer(IOTO.transcription_variable(x, tb)[2])
        @test start_value.(IOTO.transcription_variable(y, tb)) == [0. 2; 1 3]
        @test supports(x) == [(0,), (1,)]
        @test supports(y) == [(0, [0, 0]) (0, [1, 1]); (1, [0, 0]) (1, [1, 1])]
    end
    # test transcribe_parameter_functions!
    @testset "transcribe_parameter_functions!" begin
        @test isa(IOTO.transcribe_parameter_functions!(tb, m), Nothing)
        @test IOTO.transcription_variable(pf, tb) == fill(42, 2, 2)
        @test supports(pf) == [(0, [0, 0]) (0, [1, 1]); (1, [0, 0]) (1, [1, 1])]
    end
    # test transcribe_derivative_variables!
    @testset "transcribe_derivative_variables!" begin
        @test isa(IOTO.transcribe_derivative_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 7
        @test num_derivatives(m) == 4
        @test IOTO.transcription_variable(dx, tb) isa Vector{VariableRef}
        @test IOTO.transcription_variable(dy, tb) isa Matrix{VariableRef}
        @test IOTO.transcription_variable(dx3, tb) isa Vector{VariableRef}
        @test name(IOTO.transcription_variable(dx, tb)[1]) == "d/dpar[x(par)](0.0)"
        @test name(IOTO.transcription_variable(dx3, tb)[1]) == "d/dpar[d/dpar[d/dpar[x(par)]]](0.0)"
        possible = Sys.iswindows() ? "d/dpar[y(par, pars)](1.0, [0.0, 0.0])" : "∂/∂par[y(par, pars)](1.0, [0.0, 0.0])"
        @test name(IOTO.transcription_variable(dy, tb)[2, 1]) == possible
        @test has_upper_bound(IOTO.transcription_variable(dx, tb)[1])
        @test start_value.(IOTO.transcription_variable(dy, tb)) == [0. 2; 1 3]
        @test supports(dx) == [(0,), (1,)]
        @test supports(dy) == [(0, [0, 0]) (0, [1, 1]); (1, [0, 0]) (1, [1, 1])]
    end
    # test transcribe_semi_infinite_variables!
    @testset "transcribe_semi_infinite_variables!" begin 
        @test IOTO.transcribe_semi_infinite_variables!(tb, m) isa Nothing
        @test IOTO.transcription_variable(yrv) isa Vector{VariableRef}
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 12
        @test IOTO.lookup_by_support(y, tb, [1., 0., 0.]) == IOTO.lookup_by_support(yrv, tb, [1.])
        @test IOTO.transcription_variable(semi_pf) == fill(42, 2)
        @test supports(semi_pf) == [([0, 0],), ([1, 1],)]
        @test IOTO.lookup_by_support(semi_pf, tb, [0., 0.]) == 42
        @test upper_bound(IOTO.lookup_by_support(yrv2, tb, [1.])) == 3
        @test start_value(IOTO.lookup_by_support(yrv2, tb, [1.])) == 1
        @test upper_bound(IOTO.lookup_by_support(yrv3, tb, [1., 1])) == 3
        @test start_value(IOTO.lookup_by_support(yrv3, tb, [1., 1])) == 1
        @test fix_value(IOTO.lookup_by_support(yrv4, tb, [0., 0])) == 4
    end
    # test transcribe_point_variables!
    @testset "transcribe_point_variables!" begin
        @test isa(IOTO.transcribe_point_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).finvar_mappings) == 6
        @test IOTO.transcription_variable(x0, tb) == IOTO.lookup_by_support(x, tb, [0.])
        @test IOTO.transcription_variable(y0, tb) == IOTO.lookup_by_support(y, tb, [0., 0., 0.])
        @test name(IOTO.transcription_variable(x0, tb)) == "x(0.0)"
        @test name(IOTO.transcription_variable(y0, tb)) == "y(0.0, [0.0, 0.0])"
        @test lower_bound(IOTO.transcription_variable(x0, tb)) == cos(0)
        @test is_integer(IOTO.transcription_variable(x0, tb))
        @test lower_bound(IOTO.transcription_variable(y0, tb)) == 0
        @test upper_bound(IOTO.transcription_variable(y0, tb)) == 1
        @test is_binary(IOTO.transcription_variable(y0, tb))
        @test start_value(IOTO.transcription_variable(y0, tb)) == 0.
        @test has_lower_bound(IOTO.lookup_by_support(y, tb, [0., 0., 0.]))
        @test length(IOTO.transcription_data(tb).point_pfunc_mappings) == 1
        @test IOTO.transcription_variable(point_pf, tb) == 42
        @test !is_fixed(IOTO.lookup_by_support(y1, tb, [1., 1., 1.]))
        @test !has_lower_bound(IOTO.lookup_by_support(y01, tb, [0., 1., 1.]))
        @test !has_upper_bound(IOTO.lookup_by_support(y01, tb, [0., 1., 1.]))
    end
end

# Test Measure and Objective Transcription 
@testset "Measure/Objective Transcription" begin 
    # Prepare the model 
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par2 in [3, 4], supports = [3, 4])
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @parameter_function(m, pf1 == par -> sin(par))
    @parameter_function(m, pf2 == (par, par2) -> sin(par)*cos(par2))
    @variable(m, x >= 0, Infinite(par))
    @variable(m, y == 2, Infinite(par, pars), start = (p, ps) -> p + sum(ps))
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    meas1 = support_sum(x + 2w, par)
    meas2 = integral(w, par)
    meas3 = integral(y^2, pars)
    # meas4 = integral(pars[1], pars[1])
    meas5 = integral(pf1, par)
    meas6 = integral(pf2, par2)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_finite_variables!(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
    IOTO.transcribe_parameter_functions!(tb, m)
    tx = IOTO.transcription_variable(x)
    ty = IOTO.transcription_variable(y)
    tz = IOTO.transcription_variable(z)
    tw = IOTO.transcription_variable(w)
    # test transcribe_measures!
    @testset "transcribe_measures!" begin 
        @test IOTO.transcribe_measures!(tb, m) isa Nothing 
        @test IOTO.transcription_variable(meas1) == tx[1] + 4tw + tx[2]
        @test IOTO.transcription_variable(meas2) == tw + 0
        @test IOTO.transcription_variable(meas3) isa Vector
        # @test IOTO.transcription_variable(meas4) isa AffExpr
        meas5Eval = 0.5*sum(IOTO.transcription_variable.(pf1))
        @test IOTO.transcription_variable(meas5) == meas5Eval
        tpf2 = IOTO.transcription_variable.(pf2)
        meas6Eval = 0.5.*[tpf2[1] + tpf2[3], tpf2[2] + tpf2[4]]
        @test IOTO.transcription_variable(meas6) == meas6Eval
        @test supports(meas1) == ()
        @test supports(meas2) == ()
        @test supports(meas3) == [(0.,), (1., )]
        @test supports(meas5) == ()
        @test supports(meas6) == [(0.0,), (1.0,)]
    end
    # test transcribe_objective!
    @testset "transcribe_objective!" begin 
        # normal
        @objective(m, Min, 2z^2 - meas1)
        @test IOTO.transcribe_objective!(tb, m) isa Nothing 
        @test objective_sense(tb.model) == MOI.MIN_SENSE
        @test objective_function(tb.model) == 2tz^2 - tx[1] - 4tw - tx[2]
        # nonlinear objective
        @objective(m, Max, z^4)
        @test IOTO.transcribe_objective!(tb, m) isa Nothing 
        @test objective_sense(tb.model) == MOI.MAX_SENSE
        @test isequal(objective_function(tb.model), tz^4)
    end
end

# Test Constraint Methods 
@testset "Constraint Transcription" begin 
    # model setup 
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0, 0.5, 1])
    @variable(m, 0 <= x <= 2, Infinite(par, pars), Int)
    @variable(m, y == 0, Infinite(par))
    @variable(m, x0 == 0, Point(x, 0, [0, 0]))
    @variable(m, yf <= 2, Point(y, 1))
    @variable(m, z, Bin)
    @constraint(m, c1, x^2 + 2y <= 42)
    @constraint(m, c2, y - z^2 == 0, DomainRestriction(iszero, par))
    @constraint(m, c3, 4x - 3 <= 2, DomainRestriction((p, ps) -> isone(p) && all(isone.(ps)), par, pars))
    @constraint(m, c4, y + z == 0, DomainRestriction(p -> 0 <= p <= 0.5, par))
    @constraint(m, c5, 2z^2 == 0)
    @constraint(m, c6, [z, x] in MOI.Zeros(2))
    @constraint(m, c7, sin(z) ^ x == 0)
    @constraint(m, c8, integral(sin(y), par) == 0)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_finite_variables!(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
    IOTO.transcribe_point_variables!(tb, m)
    IOTO.transcribe_measures!(tb, m)
    xt = IOTO.transcription_variable(x)
    yt = IOTO.transcription_variable(y)
    x0t = IOTO.transcription_variable(x0)
    yft = IOTO.transcription_variable(yf)
    zt = IOTO.transcription_variable(z)
    # test _get_info_constr_from_var
    @testset "_get_info_constr_from_var" begin 
        # lower bounds
        set = moi_set(constraint_object(LowerBoundRef(x)))
        expected = LowerBoundRef(IOTO.lookup_by_support(x, tb, [1., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 1.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) isa Nothing
        # upper bounds
        set = moi_set(constraint_object(UpperBoundRef(x)))
        expected = UpperBoundRef(IOTO.lookup_by_support(x, tb, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) isa Nothing
        set = moi_set(constraint_object(UpperBoundRef(yf)))
        expected = UpperBoundRef(yft)
        @test IOTO._get_info_constr_from_var(tb, y, set, [1., 1., 1.]) == expected
        # fix 
        set = moi_set(constraint_object(FixRef(y)))
        expected = FixRef(IOTO.lookup_by_support(y, tb, [0.]))
        @test IOTO._get_info_constr_from_var(tb, y, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, y, set, [0., 0., 1.]) isa Nothing
        set = moi_set(constraint_object(FixRef(x0)))
        expected = FixRef(x0t)
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) == expected
        # binary 
        set = moi_set(constraint_object(BinaryRef(z)))
        @test IOTO._get_info_constr_from_var(tb, z, set, [1., 1., 1.]) == BinaryRef(zt)
        # integer
        set = moi_set(constraint_object(IntegerRef(x)))
        expected = IntegerRef(IOTO.lookup_by_support(x, tb, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 0.]) == expected
    end
    # test _support_in_restrictions
    @testset "_support_in_restrictions" begin 
        @test IOTO._support_in_restrictions(constraint_object(c1), Float64[])
        @test IOTO._support_in_restrictions(constraint_object(c2), [0.])
        @test !IOTO._support_in_restrictions(constraint_object(c2), [1.])
        @test IOTO._support_in_restrictions(constraint_object(c3), [1., 1., 1.])
        @test !IOTO._support_in_restrictions(constraint_object(c3), [0., 1., 1.])
    end
    # test _process_constraint
    @testset "_process_constraint" begin
        # scalar constraint 
        con = constraint_object(c1)
        func = jump_function(con)
        set = moi_set(con)
        @test IOTO._process_constraint(tb, con, func, set, zeros(3), "test1") isa ConstraintRef 
        @test num_constraints(tb.model, typeof(func), typeof(set)) == 1
        cref = constraint_by_name(tb.model, "test1")
        delete(tb.model, cref)
        # nonlinear scalar constraint 
        con = constraint_object(c7)
        func = jump_function(con)
        set = moi_set(con)
        @test IOTO._process_constraint(tb, con, func, set, zeros(3), "test1") isa ConstraintRef
        @test num_constraints(tb.model, typeof(func), typeof(set)) == 1
        cref = constraint_by_name(tb.model, "test1")
        delete(tb.model, cref)
        # vector constraint 
        con = constraint_object(c6)
        func = jump_function(con)
        set = moi_set(con)
        @test IOTO._process_constraint(tb, con, func, set, zeros(3), "test2") isa ConstraintRef 
        @test num_constraints(tb.model, typeof(func), typeof(set)) == 1
        cref = constraint_by_name(tb.model, "test2")
        delete(tb.model, cref)
        # test nonlinear vector constraint 
        con = VectorConstraint([sin(z)], MOI.Zeros(1))
        func = [sin(z)]
        set = MOI.Zeros(1)
        @test IOTO._process_constraint(tb, con, func, set, zeros(3), "test2") isa ConstraintRef
        # fallback
        @test_throws ErrorException IOTO._process_constraint(tb, :bad, func, set, 
                                                             zeros(3), "bad")
    end
    # test transcribe_constraints!
    @testset "transcribe_constraints!" begin 
        @test IOTO.transcribe_constraints!(tb, m) isa Nothing
        # test the info constraint transcriptions
        @test IOTO.transcription_constraint(LowerBoundRef(x)) isa Vector{ConstraintRef}
        @test IOTO.transcription_constraint(UpperBoundRef(x)) isa Vector{ConstraintRef}
        @test IOTO.transcription_constraint(IntegerRef(x)) isa Matrix{ConstraintRef}
        @test length(IOTO.transcription_constraint(LowerBoundRef(x))) == 5
        @test IOTO.transcription_constraint(FixRef(x0)) == FixRef(x0t)
        @test IOTO.transcription_constraint(FixRef(y)) == [FixRef(yt[i]) for i in 1:2]
        @test IOTO.transcription_constraint(UpperBoundRef(yf)) == UpperBoundRef(yft)
        @test IOTO.transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
        # test constraint transcriptions 
        @test IOTO.transcription_constraint(c1) isa Matrix{ConstraintRef}
        @test length(IOTO.transcription_constraint(c1)) == 6
        @test constraint_object(IOTO.transcription_constraint(c2)).func == yt[1] - zt^2
        xf = IOTO.lookup_by_support(x, tb, [1., 1., 1.])
        @test constraint_object(IOTO.transcription_constraint(c3)).func == 4xf 
        @test constraint_object(IOTO.transcription_constraint(c3)).set == MOI.LessThan(5.)
        expected = [yt[1] + zt, yt[2] + zt]
        @test jump_function.(constraint_object.(IOTO.transcription_constraint(c4))) == expected
        @test constraint_object(IOTO.transcription_constraint(c5)).func == 2zt^2 
        @test length(IOTO.transcription_constraint(c6)) == 6
        @test moi_set(constraint_object(first(IOTO.transcription_constraint(c6)))) == MOI.Zeros(2)
        @test length(IOTO.transcription_constraint(c7)) == 6
        @test IOTO.transcription_constraint(c8) isa ConstraintRef
        # test the info constraint supports 
        expected = [([1.0, 1.0], 0.0), ([0.0, 0.0], 0.5), ([1.0, 1.0], 0.5), ([0.0, 0.0], 1.0), ([1.0, 1.0], 1.0)]
        @test supports(LowerBoundRef(x)) == expected
        @test supports(UpperBoundRef(x)) == expected
        expected = [([0.0, 0.0], 0.0) ([0.0, 0.0], 0.5) ([0.0, 0.0], 1.0); 
                    ([1.0, 1.0], 0.0) ([1.0, 1.0], 0.5) ([1.0, 1.0], 1.0)]
        @test supports(IntegerRef(x)) == expected
        @test supports(FixRef(x0)) == ()
        @test supports(UpperBoundRef(yf)) == ()
        @test supports(BinaryRef(z)) == ()
        # test the constraint supports 
        expected = [([0., 0.], 0.) ([0., 0.], 0.5) ([0., 0.], 1.); ([1., 1.], 0.) ([1., 1.], 0.5) ([1., 1.], 1.)]
        @test supports(c1) == expected
        @test supports(c2) == (0.,)
        @test supports(c3) == ([1., 1.], 1.)
        @test supports(c4) == [(0.0,), (0.5,)]
        @test supports(c5) == ()
        @test supports(c6) == expected
        @test supports(c7) == expected
        @test supports(c8) == ()
    end
end

# Test transcribe_derivative_evaluations! (TODO SEE IF THIS CAN BE MADE TO WORK WITH DEPENDENT PARAMETER BASED DERIVATIVES)
@testset "transcribe_derivative_evaluations!" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1], num_supports = 4, 
                        derivative_method = OrthogonalCollocation(4))
    @variable(m, y, Infinite(t))
    d1 = deriv(y, t)
    d2 = @deriv(y, t^2)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
    IOTO.transcribe_derivative_variables!(tb, m)
    # main test 
    @test IOTO.transcribe_derivative_evaluations!(tb, m) isa Nothing 
    @test num_constraints(tb.model, AffExpr, MOI.EqualTo{Float64}) == 18
    @test num_supports(t) == 4
    @test num_supports(t, label = All) == 10
    @test length(supports(d1)) == 4
    @test length(supports(d2, label = All)) == 10
    @test IOTO.has_internal_supports(tb)
    @test has_generative_supports(t)
end

@testset "transcribe_variable_collocation_restictions!" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1], num_supports = 4, 
                        derivative_method = OrthogonalCollocation(3))
    @infinite_parameter(m, x in [0, 1], num_supports = 3, 
                        derivative_method = OrthogonalCollocation(2))
    @variable(m, y, Infinite(t, x))
    constant_over_collocation(y, t)
    constant_over_collocation(y, x)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
    # main test 
    @test IOTO.transcribe_variable_collocation_restictions!(tb, m) isa Nothing
    @test num_constraints(tb.model, count_variable_in_set_constraints = false) == 3 * 3
    yt = IOTO.transcription_variable(y, label = All)
    cons = all_constraints(tb.model, include_variable_in_set_constraints = false)
    @test jump_function(constraint_object(first(cons))) == yt[7, 1] - yt[6, 1] 
    # test assertion error
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1], supports = [0, 0.9], 
                        derivative_method = OrthogonalCollocation(3))
    @variable(m, y, Infinite(t))
    add_supports(t, 1, label = InternalGaussLobatto)
    constant_over_collocation(y, t)
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
    @test_throws AssertionError IOTO.transcribe_variable_collocation_restictions!(tb, m)
end

# Test build_transcription_backend!
@testset "build_transcription_backend!" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], derivative_method = OrthogonalCollocation(2))
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @variable(m, 1 >= x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = 0)
    x0 = x(0)
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]))
    y1 = y(1, [1, 1])
    unfix(y1)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    @finite_parameter(m, fin == 0)
    f = parameter_function(a -> 0, par)
    d1 = deriv(y, par)
    d2 = @deriv(x, par^2)
    set_upper_bound(d1, 2)
    meas1 = support_sum(x - w - d2, par)
    meas2 = support_sum(y, pars)
    @constraint(m, c1, x + par - z + d1 + f == 0, 
                DomainRestriction((p, ps) -> iszero(p) && all(iszero.(ps)), par, pars))
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x + fin <= 1, 
                DomainRestriction(p -> 0.5 <= p <= 1, par))
    @constraint(m, c5, meas2 == 0)
    @constraint(m, x + y == 83)
    @constraint(m, c6, [z, w] in MOI.Zeros(2))
    g(a) = 42
    @operator(m, gr, 1, g)
    @constraint(m, c7, gr(z) == 2)
    @objective(m, Min, x0 + meas1)
    # test basic usage
    tb = m.backend
    @test IOTO.build_transcription_backend!(
        tb,
        m, 
        check_support_dims = false
    ) isa Nothing
    # test finite variables
    zt = IOTO.transcription_variable(z)
    wt = IOTO.transcription_variable(w)
    @test zt isa VariableRef
    @test wt isa VariableRef
    @test name(zt) == name(z)
    @test name(wt) == name(w)
    @test has_lower_bound(zt)
    @test has_upper_bound(zt)
    @test is_binary(zt)
    @test is_fixed(wt)
    @test is_integer(wt)
    @test start_value(wt) == 1.
    # test infinite variables
    @test IOTO.transcription_variable(x) isa Vector{VariableRef}
    @test IOTO.transcription_variable(y) isa Matrix{VariableRef}
    @test name(IOTO.transcription_variable(x)[1]) == "x(0.0)"
    @test name(IOTO.transcription_variable(y)[1, 2]) == "y(0.0, [1.0, 1.0])"
    @test has_lower_bound(IOTO.transcription_variable(x)[1])
    @test is_binary(IOTO.transcription_variable(y)[2])
    @test is_fixed(IOTO.transcription_variable(y)[3])
    @test is_integer(IOTO.transcription_variable(x)[2])
    @test start_value(IOTO.transcription_variable(y)[1]) == 0.
    @test supports(x) == [(0.,), (1.,)]
    @test supports(y) == [(0.0, [0.0, 0.0]) (0.0, [1.0, 1.0]); (1.0, [0.0, 0.0]) (1.0, [1.0, 1.0])]
    # test point variables
    @test IOTO.transcription_variable(x0) isa VariableRef
    @test IOTO.transcription_variable(y0) isa VariableRef
    @test name(IOTO.transcription_variable(x0)) == "x(0.0)"
    @test name(IOTO.transcription_variable(y0)) == "y(0.0, [0.0, 0.0])"
    @test has_lower_bound(IOTO.transcription_variable(x0))
    @test is_integer(IOTO.transcription_variable(x0))
    @test has_lower_bound(IOTO.transcription_variable(y0))
    @test has_upper_bound(IOTO.transcription_variable(y0))
    @test is_binary(IOTO.transcription_variable(y0))
    @test start_value(IOTO.transcription_variable(y0)) == 0.
    @test !is_fixed(IOTO.transcription_variable(y1))
    # test derivatives 
    d1t = IOTO.transcription_variable(d1)
    d2t = IOTO.transcription_variable(d2)
    @test length(d1t) == 4
    @test length(d2t) == 2
    @test upper_bound(d1t[1]) == 2
    @test supports(d2) == [(0.,), (1.,)]
    # test operators
    attr_dict = backend(tb).model_cache.modattr
    @test length(attr_dict) == 1
    @test attr_dict[MOI.UserDefinedFunction(:gr, 1)] == (g,)
    # test objective
    xt = IOTO.transcription_variable(x, tb)
    @test objective_function(tb.model) == 2xt[1] + xt[2] - 2wt - d2t[1] - d2t[2]
    @test objective_sense(tb.model) == MOI.MIN_SENSE
    # test constraints
    yt = IOTO.transcription_variable(y)
    ft = IOTO.transformation_variable(f)
    dt_c1 = IOTO.lookup_by_support(d1, tb, zeros(3))
    @test constraint_object(IOTO.transcription_constraint(c1)).func == -zt + xt[1] + dt_c1 + ft[1]
    @test constraint_object(IOTO.transcription_constraint(c2)).func == zt + xt[1]
    expected = IOTO.transcription_variable(meas2)[2] - 2 * IOTO.transcription_variable(y0) + xt[2] + IOTO.transcription_variable(fin)
    @test constraint_object(IOTO.transcription_constraint(c4)).func == expected
    @test constraint_object(IOTO.transcription_constraint(c3)).func == xt[1] - 2wt + xt[2] + zt - d2t[1] - d2t[2]
    @test constraint_object(IOTO.transcription_constraint(c6)).func == [zt, wt]
    @test IOTO.transcription_constraint(c5) isa Vector{ConstraintRef}
    @test name(IOTO.transcription_constraint(c2)) == "c2"
    @test name(IOTO.transcription_constraint(c1)) == "c1[1, 1]"
    @test supports(c1) == (0., [0., 0.])
    @test IOTO.transcription_constraint(c7) isa ConstraintRef
    @test isequal(constraint_object(IOTO.transcription_constraint(c7)).func, gr(zt) - 2.)
    # test info constraints
    @test IOTO.transcription_constraint(LowerBoundRef(z)) == LowerBoundRef(zt)
    @test IOTO.transcription_constraint(UpperBoundRef(z)) == UpperBoundRef(zt)
    @test IOTO.transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
    @test IOTO.transcription_constraint(FixRef(w)) == FixRef(wt)
    @test IOTO.transcription_constraint(IntegerRef(w)) == IntegerRef(wt)
    @test IOTO.transcription_constraint(LowerBoundRef(x0))[1] == LowerBoundRef(xt[1])
    @test IOTO.transcription_constraint(UpperBoundRef(x0))[1] == UpperBoundRef(xt[1])
    @test IOTO.transcription_constraint(IntegerRef(x0))[1] == IntegerRef(xt[1])
    @test IOTO.transcription_constraint(LowerBoundRef(y0)) isa ConstraintRef
    @test IOTO.transcription_constraint(UpperBoundRef(y0)) isa ConstraintRef
    @test IOTO.transcription_constraint(BinaryRef(y0)) isa Matrix{ConstraintRef}
    @test IOTO.transcription_constraint(LowerBoundRef(x)) == LowerBoundRef.(xt)
    @test IOTO.transcription_constraint(UpperBoundRef(x)) == UpperBoundRef.(xt)
    @test IOTO.transcription_constraint(IntegerRef(x)) == IntegerRef.(xt)
    @test IOTO.transcription_constraint(FixRef(y)) isa Vector{ConstraintRef}
    @test IOTO.transcription_constraint(BinaryRef(y)) isa Matrix{ConstraintRef}
    @test IOTO.transcription_constraint(UpperBoundRef(d1)) == UpperBoundRef.(d1t)

    # test a finite model 
    m = InfiniteModel()
    @variable(m, y >= 0)
    @objective(m, Min, y)
    @constraint(m, y^2 <= 42)
    tb = m.backend
    @test IOTO.build_transcription_backend!(tb, m) isa Nothing 
    @test IOTO.transcription_variable(y) isa VariableRef 
    @test lower_bound(IOTO.transcription_variable(y)) == 0
    @test objective_sense(tb.model) == MOI.MIN_SENSE
end

# Test build_transformation_backend!
@testset "build_transformation_backend!" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    @infinite_parameter(m, par in [0, 1], num_supports = 3)
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @variable(m, 1 >= x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = 0)
    x0 = x(0)
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]))
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    @finite_parameter(m, p == 0)
    @parameter_function(m, pf == sin(par))
    meas1 = support_sum(x - w + pf, par)
    meas2 = integral(y, pars)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == p)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, DomainRestriction(p -> 0.5 <= p <= 1, par))
    @constraint(m, c5, meas2 == 0)
    @constraint(m, deriv(x, par) == 0)
    @constraint(m, sin(w) + integral(x^3, par) == 0)
    @objective(m, Min, x0 + meas1)
    set_silent(m)
    set_time_limit_sec(m, 42.)
    # test normal usage
    @test isa(build_transformation_backend!(m), Nothing)
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 45
    @test time_limit_sec(m.backend.model) == 42
    # test bad keyword
    @test_throws ErrorException build_transformation_backend!(m, bad = 42)
    # test parameter update mode
    new_backend = TranscriptionBackend(update_parameter_functions = true)
    @test set_transformation_backend(m, new_backend) isa Nothing
    @test build_transformation_backend!(m) isa Nothing
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 48
    @test set_parameter_value(p, 42) isa Nothing
    @test parameter_value(transformation_variable(p)) == 42
    @test set_parameter_value(pf, cos) isa Nothing
    @test parameter_value.(transformation_variable(pf)) == [cos(0), cos(0.5), cos(1)]
end