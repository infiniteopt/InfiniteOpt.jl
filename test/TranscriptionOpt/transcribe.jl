# Test variable initializers
@testset "Variable Initializers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], derivative_method = OrthogonalCollocation(2))
    @variable(m, x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = (p, ps) -> p + sum(ps))
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    dx = @deriv(x, par)
    dy = @deriv(y, par)
    dx3 = @deriv(x, par^3)
    set_lower_bound(dx, 0)
    set_start_value_function(dy, (p, ps) -> p + sum(ps))
    var = build_variable(error, y, Dict{Int, Float64}(2 => 0))
    yrv = add_variable(m, var)
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
        var = InfiniteOpt._core_variable_object(y)
        @test !IOTO._format_infinite_info(var, [0., 0.5, 0.75]).has_lb
        @test !IOTO._format_infinite_info(var, [0., 0.5, 0.75]).has_ub
        @test IOTO._format_infinite_info(var, [0., 0.5, 0.75]).has_fix
        @test IOTO._format_infinite_info(var, [0., 0.5, 0.75]).fixed_value == 2
        @test IOTO._format_infinite_info(var, [0., 0.5, 0.75]).has_start
        @test IOTO._format_infinite_info(var, [0., 0.5, 0.75]).start == 1.25
        @test IOTO._format_infinite_info(var, [0., 0.5, 0.75]).binary 
        @test !IOTO._format_infinite_info(var, [0., 0.5, 0.75]).integer
        var = InfiniteOpt._core_variable_object(x)
        @test IOTO._format_infinite_info(var, [0.]).has_lb
        @test IOTO._format_infinite_info(var, [0.]).lower_bound == 0
        @test !IOTO._format_infinite_info(var, [0.]).has_ub
        @test !IOTO._format_infinite_info(var, [0.]).has_fix
        @test !IOTO._format_infinite_info(var, [0.]).has_start
        @test isnan(IOTO._format_infinite_info(var, [0.]).start)
        @test IOTO._format_infinite_info(var, [0.]).integer
    end
    # test transcribe_infinite_variables!
    @testset "transcribe_infinite_variables!" begin
        @test isa(IOTO.transcribe_infinite_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 2
        @test IOTO.transcription_variable(x, tb) isa Vector{VariableRef}
        @test IOTO.transcription_variable(y, tb) isa Vector{VariableRef}
        @test name(IOTO.transcription_variable(x, tb)[1]) == "x(support: 1)"
        @test name(IOTO.transcription_variable(y, tb)[3]) == "y(support: 3)"
        @test has_lower_bound(IOTO.transcription_variable(x)[1])
        @test is_binary(IOTO.transcription_variable(y, tb)[2])
        @test is_fixed(IOTO.transcription_variable(y, tb)[4])
        @test is_integer(IOTO.transcription_variable(x, tb)[2])
        @test sort!(start_value.(IOTO.transcription_variable(y, tb))) == [0., 1, 2, 3]
        @test supports(x) == [(0,), (1,)]
        @test length(supports(y)) == 4
    end
    # test _format_derivative_info
    @testset "_format_derivative_info" begin
        der = InfiniteOpt._core_variable_object(dy)
        @test !IOTO._format_derivative_info(der, [0., 0.5, 0.75]).has_lb
        @test !IOTO._format_derivative_info(der, [0., 0.5, 0.75]).has_ub
        @test !IOTO._format_derivative_info(der, [0., 0.5, 0.75]).has_fix
        @test IOTO._format_derivative_info(der, [0., 0.5, 0.75]).has_start
        @test IOTO._format_derivative_info(der, [0., 0.5, 0.75]).start == 1.25
        @test !IOTO._format_derivative_info(der, [0., 0.5, 0.75]).binary 
        @test !IOTO._format_derivative_info(der, [0., 0.5, 0.75]).integer
        der = InfiniteOpt._core_variable_object(dx)
        @test IOTO._format_derivative_info(der, [0.]).has_lb
        @test IOTO._format_derivative_info(der, [0.]).lower_bound == 0
        @test !IOTO._format_derivative_info(der, [0.]).has_ub
        @test !IOTO._format_derivative_info(der, [0.]).has_fix
        @test !IOTO._format_derivative_info(der, [0.]).has_start
        @test isnan(IOTO._format_derivative_info(der, [0.]).start)
        @test !IOTO._format_derivative_info(der, [0.]).integer
    end
    # test transcribe_derivative_variables!
    @testset "transcribe_derivative_variables!" begin
        @test isa(IOTO.transcribe_derivative_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 6
        @test num_derivatives(m) == 4
        @test IOTO.transcription_variable(dx, tb) isa Vector{VariableRef}
        @test IOTO.transcription_variable(dy, tb) isa Vector{VariableRef}
        @test IOTO.transcription_variable(dx3, tb) isa Vector{VariableRef}
        @test name(IOTO.transcription_variable(dx, tb)[1]) == "d/dpar[x(par)](support: 1)"
        @test name(IOTO.transcription_variable(dx3, tb)[1]) == "d^3/dpar^3[x(par)](support: 1)"
        @test name(IOTO.transcription_variable(deriv(dx, par), tb)[1]) == "d²/dpar²[x(par)](support: 1)"
        @test name(IOTO.transcription_variable(dy, tb)[3]) == (Sys.iswindows() ? "d/dpar[y(par, pars)](support: 3)" : "∂/∂par[y(par, pars)](support: 3)")
        @test has_lower_bound(IOTO.transcription_variable(dx, tb)[1])
        @test sort!(start_value.(IOTO.transcription_variable(dy, tb))) == [0., 1, 2, 3]
        @test supports(dx) == [(0,), (1,)]
        @test length(supports(dy)) == 4
    end
    # test _set_semi_infinite_variable_mapping
    @testset "_set_semi_infinite_variable_mapping" begin 
        var = SemiInfiniteVariable(y, Dict{Int, Float64}(1 => 0), [1, 2], [1])
        vref = GeneralVariableRef(m, -1, SemiInfiniteVariableIndex)
        @test IOTO._set_semi_infinite_variable_mapping(tb, var, vref, SemiInfiniteVariableIndex) isa Nothing 
        @test IOTO.transcription_variable(vref) isa Vector{VariableRef}
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 7
        @test IOTO.lookup_by_support(y, tb, [0., 0, 0]) == IOTO.lookup_by_support(vref, tb, [0., 0])
        @test IOTO._set_semi_infinite_variable_mapping(tb, var, vref, ParameterFunctionIndex) isa Nothing 
    end
    # test transcribe_semi_infinite_variables!
    @testset "transcribe_semi_infinite_variables!" begin 
        @test IOTO.transcribe_semi_infinite_variables!(tb, m) isa Nothing
        @test IOTO.transcription_variable(yrv) isa Vector{VariableRef}
        @test length(IOTO.transcription_data(tb).infvar_mappings) == 8
        @test IOTO.lookup_by_support(y, tb, [1., 0., 0.]) == IOTO.lookup_by_support(yrv, tb, [1., 0])
    end
    # test _update_point_info
    @testset "_update_point_info" begin
        vref = @variable(tb.model)
        # test lower bound update
        fix(vref, 0, force = true)
        set_lower_bound(x0, 0)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test lower_bound(vref) == 0
        # test upper bound update
        fix(vref, 0, force = true)
        set_upper_bound(x0, 0)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test upper_bound(vref) == 0
        # test fix update
        fix(x0, 0, force = true)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test fix_value(vref) == 0
        # test binary update
        unset_integer(x0)
        set_binary(x0)
        set_integer(vref)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test is_binary(vref)
        # test integer update
        unset_binary(x0)
        set_integer(x0)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test is_integer(vref)
        # test start value
        set_start_value(x0, 0.5)
        @test isa(IOTO._update_point_info(x0, vref), Nothing)
        @test start_value(vref) == 0.5
        # Undo changes
        unfix(x0)
        set_lower_bound(x0, 0)
        set_start_value(x0, NaN)
    end
    # test transcribe_point_variables!
    @testset "transcribe_point_variables!" begin
        @test isa(IOTO.transcribe_point_variables!(tb, m), Nothing)
        @test length(IOTO.transcription_data(tb).finvar_mappings) == 4
        @test IOTO.transcription_variable(x0, tb) == IOTO.lookup_by_support(x, tb, [0.])
        @test IOTO.transcription_variable(y0, tb) == IOTO.lookup_by_support(y, tb, [0., 0., 0.])
        @test name(IOTO.transcription_variable(x0, tb)) == "x(support: 1)"
        @test name(IOTO.transcription_variable(y0, tb))[1:end-2] == "y(support: "
        @test lower_bound(IOTO.transcription_variable(x0, tb)) == 0
        @test is_integer(IOTO.transcription_variable(x0, tb))
        @test lower_bound(IOTO.transcription_variable(y0, tb)) == 0
        @test upper_bound(IOTO.transcription_variable(y0, tb)) == 1
        @test is_integer(IOTO.transcription_variable(y0, tb))
        @test start_value(IOTO.transcription_variable(y0, tb)) == 0.
        @test has_lower_bound(IOTO.lookup_by_support(y, tb, [0., 0., 0.]))
    end
end

# Test Measure and Objective Transcription 
@testset "Measure/Objective Transcription" begin 
    # Prepare the model 
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, x >= 0, Infinite(par))
    @variable(m, y == 2, Infinite(par, pars), start = (p, ps) -> p + sum(ps))
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    meas1 = support_sum(x + 2w, par)
    meas2 = integral(w, par)
    meas3 = integral(y^2, pars)
    meas4 = integral(pars[1], pars[1])
    tb = m.backend
    IOTO.set_parameter_supports(tb, m)
    IOTO.transcribe_finite_variables!(tb, m)
    IOTO.transcribe_infinite_variables!(tb, m)
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
        @test IOTO.transcription_variable(meas4) isa AffExpr
        @test supports(meas1) == ()
        @test supports(meas2) == ()
        @test sort!(supports(meas3)) == [(0.,), (1., )]
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
    @variable(m, x0 == 0, Point(x, 0, [0, 0]), Bin)
    @variable(m, yf <= 2, Point(y, 1), Int)
    @variable(m, z, Bin)
    @constraint(m, c1, x^2 + 2y <= 42)
    @constraint(m, c2, y - z^2 == 0, DomainRestrictions(par => 0))
    @constraint(m, c3, 4x - 3 <= 2, DomainRestrictions(pars => 1, par => 1))
    @constraint(m, c4, y + z == 0, DomainRestrictions(par => [0, 0.5]))
    @constraint(m, c5, 2z^2 == 0, DomainRestrictions(par => 1))
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
        set = moi_set(InfiniteOpt._core_constraint_object(LowerBoundRef(x)))
        expected = LowerBoundRef(IOTO.lookup_by_support(x, tb, [1., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 1.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) isa Nothing
        # upper bounds
        set = moi_set(InfiniteOpt._core_constraint_object(UpperBoundRef(x)))
        expected = UpperBoundRef(IOTO.lookup_by_support(x, tb, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(UpperBoundRef(yf)))
        expected = UpperBoundRef(yft)
        @test IOTO._get_info_constr_from_var(tb, y, set, [1., 1., 1.]) == expected
        # fix 
        set = moi_set(InfiniteOpt._core_constraint_object(FixRef(y)))
        expected = FixRef(IOTO.lookup_by_support(y, tb, [0.]))
        @test IOTO._get_info_constr_from_var(tb, y, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, y, set, [0., 0., 1.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(FixRef(x0)))
        expected = FixRef(x0t)
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) == expected
        # binary 
        set = moi_set(InfiniteOpt._core_constraint_object(BinaryRef(x0)))
        expected = BinaryRef(IOTO.lookup_by_support(x, tb, [0., 0., 0.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 1.]) isa Nothing
        @test IOTO._get_info_constr_from_var(tb, z, set, [1., 1., 1.]) == BinaryRef(zt)
        # integer
        set = moi_set(InfiniteOpt._core_constraint_object(IntegerRef(x)))
        expected = IntegerRef(IOTO.lookup_by_support(x, tb, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tb, x, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tb, x, set, [0., 0., 0.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(IntegerRef(yf)))
        expected = IntegerRef(yft)
        @test IOTO._get_info_constr_from_var(tb, yf, set, [1., 1., 1.]) == expected
    end
    # test _support_in_restrictions
    @testset "_support_in_restrictions" begin 
        @test IOTO._support_in_restrictions([0., 0., 0.], [3], [IntervalDomain(0, 0)])
        @test IOTO._support_in_restrictions([0., 0., 0.], Int[], IntervalDomain[])
        @test IOTO._support_in_restrictions([NaN, 0., 0.], [1], [IntervalDomain(0, 1)])
        @test !IOTO._support_in_restrictions([NaN, 0., 0.], [1, 2], [IntervalDomain(1, 1), IntervalDomain(1, 1)])
        @test !IOTO._support_in_restrictions([NaN, 0., 2.], [1, 3], [IntervalDomain(1, 1), IntervalDomain(1, 1)])
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
        @test IOTO.transcription_constraint(IntegerRef(x)) isa Vector{ConstraintRef}
        @test length(IOTO.transcription_constraint(LowerBoundRef(x))) == 5
        @test IOTO.transcription_constraint(FixRef(x0)) == FixRef(x0t)
        @test IOTO.transcription_constraint(BinaryRef(x0)) == BinaryRef(x0t)
        @test IOTO.transcription_constraint(FixRef(y)) == [FixRef(yt[i]) for i in 1:2]
        @test IOTO.transcription_constraint(UpperBoundRef(yf)) == UpperBoundRef(yft)
        @test IOTO.transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
        # test constraint transcriptions 
        @test IOTO.transcription_constraint(c1) isa Vector{ConstraintRef}
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
        expected = [([0., 0.], 0.5), ([0., 0.], 1.), ([1., 1.], 0.), ([1., 1.], 0.5), ([1., 1.], 1.)]
        @test sort(supports(LowerBoundRef(x))) == expected
        @test sort(supports(UpperBoundRef(x))) == expected
        @test sort(supports(IntegerRef(x))) == expected
        @test supports(FixRef(x0)) == ()
        @test supports(UpperBoundRef(yf)) == ()
        @test supports(BinaryRef(z)) == ()
        # test the constraint supports 
        expected = [([0., 0.], 0.), ([0., 0.], 0.5), ([0., 0.], 1.), ([1., 1.], 0.), ([1., 1.], 0.5), ([1., 1.], 1.)]
        @test sort(supports(c1)) == expected
        @test supports(c2) == (0.,)
        @test supports(c3) == ([1., 1.], 1.)
        @test supports(c4) == [(0.0,), (0.5,)]
        @test supports(c5) == ()
        @test sort(supports(c6)) == expected
        @test sort(supports(c7)) == expected
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
    d1 = @deriv(y, t)
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
    yt = IOTO.transcription_variable(y, label = All, ndarray = true)
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
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    @finite_parameter(m, fin == 0)
    f = parameter_function(a -> 0, par)
    d1 = @deriv(y, par)
    d2 = @deriv(x, par^2)
    set_upper_bound(d1, 2)
    meas1 = support_sum(x - w - d2, par)
    meas2 = support_sum(y, pars)
    @constraint(m, c1, x + par - z + d1 + f == 0, 
                DomainRestrictions(pars => 0, par => 0))
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x + fin <= 1, 
                DomainRestrictions(par => [0.5, 1]))
    @constraint(m, c5, meas2 == 0)
    @constraint(m, x + y == 83)
    @constraint(m, c6, [z, w] in MOI.Zeros(2))
    g(a) = 42
    @operator(m, gr, 1, g)
    @constraint(m, c7, gr(z) == 2)
    @objective(m, Min, x0 + meas1)
    # test basic usage
    tb = m.backend
    @test IOTO.build_transcription_backend!(tb, m, 
                                          check_support_dims = false) isa Nothing
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
    @test IOTO.transcription_variable(y) isa Vector{VariableRef}
    @test name(IOTO.transcription_variable(x)[1]) == "x(support: 1)"
    @test name(IOTO.transcription_variable(y)[3]) == "y(support: 3)"
    @test has_lower_bound(IOTO.transcription_variable(x)[1])
    @test is_binary(IOTO.transcription_variable(y)[2])
    @test is_fixed(IOTO.transcription_variable(y)[4])
    @test is_integer(IOTO.transcription_variable(x)[2])
    @test start_value(IOTO.transcription_variable(y)[1]) == 0.
    @test supports(x) == [(0.,), (1.,)]
    @test length(supports(y)) == 4
    # test point variables
    @test IOTO.transcription_variable(x0) isa VariableRef
    @test IOTO.transcription_variable(y0) isa VariableRef
    @test name(IOTO.transcription_variable(x0)) == "x(support: 1)"
    @test name(IOTO.transcription_variable(y0))[1:end-2] == "y(support: "
    @test has_lower_bound(IOTO.transcription_variable(x0))
    @test is_integer(IOTO.transcription_variable(x0))
    @test has_lower_bound(IOTO.transcription_variable(y0))
    @test has_upper_bound(IOTO.transcription_variable(y0))
    @test is_integer(IOTO.transcription_variable(y0))
    @test start_value(IOTO.transcription_variable(y0)) == 0.
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
    dt_c1 = IOTO.lookup_by_support(d1, tb, zeros(3))
    @test constraint_object(IOTO.transcription_constraint(c1)).func == -zt + xt[1] + dt_c1
    @test constraint_object(IOTO.transcription_constraint(c2)).func == zt + xt[1]
    expected = IOTO.transcription_variable(meas2)[2] - 2 * IOTO.transcription_variable(y0) + xt[2]
    @test constraint_object(IOTO.transcription_constraint(c4)).func == expected
    @test constraint_object(IOTO.transcription_constraint(c3)).func == xt[1] - 2wt + xt[2] + zt - d2t[1] - d2t[2]
    @test constraint_object(IOTO.transcription_constraint(c6)).func == [zt, wt]
    @test IOTO.transcription_constraint(c5) isa Vector{ConstraintRef}
    @test name(IOTO.transcription_constraint(c2)) == "c2(support: 1)"
    @test name(IOTO.transcription_constraint(c1)) == "c1(support: 1)"
    @test supports(c1) == (0., [0., 0.])
    @test IOTO.transcription_constraint(c7) isa ConstraintRef
    @test isequal(constraint_object(IOTO.transcription_constraint(c7)).func, gr(zt) - 2.)
    # test info constraints
    @test IOTO.transcription_constraint(LowerBoundRef(z)) == LowerBoundRef(zt)
    @test IOTO.transcription_constraint(UpperBoundRef(z)) == UpperBoundRef(zt)
    @test IOTO.transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
    @test IOTO.transcription_constraint(FixRef(w)) == FixRef(wt)
    @test IOTO.transcription_constraint(IntegerRef(w)) == IntegerRef(wt)
    @test IOTO.transcription_constraint(LowerBoundRef(x0)) == LowerBoundRef(xt[1])
    @test IOTO.transcription_constraint(UpperBoundRef(x0)) == UpperBoundRef(xt[1])
    @test IOTO.transcription_constraint(IntegerRef(x0)) == IntegerRef(xt[1])
    @test IOTO.transcription_constraint(LowerBoundRef(y0)) isa ConstraintRef
    @test IOTO.transcription_constraint(UpperBoundRef(y0)) isa ConstraintRef
    @test IOTO.transcription_constraint(IntegerRef(y0)) isa ConstraintRef
    @test IOTO.transcription_constraint(LowerBoundRef(x)) == LowerBoundRef.(xt)
    @test IOTO.transcription_constraint(UpperBoundRef(x)) == UpperBoundRef.(xt)
    @test IOTO.transcription_constraint(IntegerRef(x)) == IntegerRef.(xt)
    @test IOTO.transcription_constraint(FixRef(y)) isa Vector{ConstraintRef}
    @test IOTO.transcription_constraint(BinaryRef(y)) isa Vector{ConstraintRef}
    @test IOTO.transcription_constraint(UpperBoundRef(d1)) == UpperBoundRef.(d1t)

    # test a finite model 
    m = InfiniteModel()
    @variable(m, y >= 0)
    @objective(m, Min, y)
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
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    meas1 = support_sum(x - w, par)
    meas2 = integral(y, pars)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, DomainRestrictions(par => [0.5, 1]))
    @constraint(m, c5, meas2 == 0)
    @constraint(m, @deriv(x, par) == 0)
    @constraint(m, sin(w) + integral(x^3, par) == 0)
    @objective(m, Min, x0 + meas1)
    set_silent(m)
    set_time_limit_sec(m, 42.)
    # test normal usage
    @test isa(build_transformation_backend!(m, m.backend), Nothing)
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 44
    @test time_limit_sec(m.backend) == 42
end