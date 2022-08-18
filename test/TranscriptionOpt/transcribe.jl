# Test variable initializers
@testset "Variable Initializers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @variable(m, x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = (p, ps) -> p + sum(ps))
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    dx = @deriv(x, par)
    dy = @deriv(y, par)
    set_lower_bound(dx, 0)
    set_start_value_function(dy, (p, ps) -> p + sum(ps))
    var = build_variable(error, y, Dict{Int, Float64}(2 => 0))
    yrv = add_variable(m, var)
    data = DiscreteMeasureData(par, [0.5, 0.5], [0, 1])
    @constraint(m, c1, x + z - 2 <= 0)
    @constraint(m, c2, measure(x + y, data) - w == 0)
    @constraint(m, c3, x0 + y0 == 5)
    tm = optimizer_model(m)
    IOTO.set_parameter_supports(tm, m)
    # test transcribe_finite_variables!
    @testset "transcribe_finite_variables!" begin
        @test isa(IOTO.transcribe_finite_variables!(tm, m), Nothing)
        @test length(transcription_data(tm).finvar_mappings) == 2
        @test transcription_variable(tm, z) isa VariableRef
        @test transcription_variable(tm, w) isa VariableRef
        @test name(transcription_variable(tm, z)) == name(z)
        @test name(transcription_variable(tm, w)) == name(w)
        @test has_lower_bound(transcription_variable(tm, z))
        @test has_upper_bound(transcription_variable(tm, z))
        @test is_binary(transcription_variable(tm, z))
        @test is_fixed(transcription_variable(tm, w))
        @test is_integer(transcription_variable(tm, w))
        @test start_value(transcription_variable(tm, w)) == 1.
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
        @test isa(IOTO.transcribe_infinite_variables!(tm, m), Nothing)
        @test length(transcription_data(tm).infvar_mappings) == 2
        @test transcription_variable(tm, x) isa Vector{VariableRef}
        @test transcription_variable(tm, y) isa Vector{VariableRef}
        @test name(transcription_variable(tm, x)[1]) == "x(support: 1)"
        @test name(transcription_variable(tm, y)[3]) == "y(support: 3)"
        @test has_lower_bound(transcription_variable(tm, x)[1])
        @test is_binary(transcription_variable(tm, y)[2])
        @test is_fixed(transcription_variable(tm, y)[4])
        @test is_integer(transcription_variable(tm, x)[2])
        @test sort!(start_value.(transcription_variable(tm, y))) == [0., 1, 2, 3]
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
        @test isa(IOTO.transcribe_derivative_variables!(tm, m), Nothing)
        @test length(transcription_data(tm).infvar_mappings) == 4
        @test transcription_variable(tm, dx) isa Vector{VariableRef}
        @test transcription_variable(tm, dy) isa Vector{VariableRef}
        @test name(transcription_variable(tm, dx)[1]) == (Sys.iswindows() ? "d/dpar[x(par)](support: 1)" : "∂/∂par[x(par)](support: 1)")
        @test name(transcription_variable(tm, dy)[3]) == (Sys.iswindows() ? "d/dpar[y(par, pars)](support: 3)" : "∂/∂par[y(par, pars)](support: 3)")
        @test has_lower_bound(transcription_variable(tm, dx)[1])
        @test sort!(start_value.(transcription_variable(tm, dy))) == [0., 1, 2, 3]
        @test supports(dx) == [(0,), (1,)]
        @test length(supports(dy)) == 4
    end
    # test _set_semi_infinite_variable_mapping
    @testset "_set_semi_infinite_variable_mapping" begin 
        var = SemiInfiniteVariable(y, Dict{Int, Float64}(1 => 0), [1, 2], [1])
        vref = GeneralVariableRef(m, -1, SemiInfiniteVariableIndex)
        @test IOTO._set_semi_infinite_variable_mapping(tm, var, vref, SemiInfiniteVariableIndex) isa Nothing 
        @test transcription_variable(vref) isa Vector{VariableRef}
        @test length(transcription_data(tm).infvar_mappings) == 5
        @test IOTO.lookup_by_support(tm, y, [0., 0, 0]) == IOTO.lookup_by_support(tm, vref, [0., 0])
        @test IOTO._set_semi_infinite_variable_mapping(tm, var, vref, ParameterFunctionIndex) isa Nothing 
    end
    # test transcribe_semi_infinite_variables!
    @testset "transcribe_semi_infinite_variables!" begin 
        @test IOTO.transcribe_semi_infinite_variables!(tm, m) isa Nothing
        @test transcription_variable(yrv) isa Vector{VariableRef}
        @test length(transcription_data(tm).infvar_mappings) == 6
        @test IOTO.lookup_by_support(tm, y, [1., 0., 0.]) == IOTO.lookup_by_support(tm, yrv, [1., 0])
    end
    # test _update_point_info
    @testset "_update_point_info" begin
        vref = @variable(tm)
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
        @test isa(IOTO.transcribe_point_variables!(tm, m), Nothing)
        @test length(transcription_data(tm).finvar_mappings) == 4
        @test transcription_variable(tm, x0) == IOTO.lookup_by_support(tm, x, [0.])
        @test transcription_variable(tm, y0) == IOTO.lookup_by_support(tm, y, [0., 0., 0.])
        @test name(transcription_variable(tm, x0)) == "x(support: 1)"
        @test name(transcription_variable(tm, y0))[1:end-2] == "y(support: "
        @test lower_bound(transcription_variable(tm, x0)) == 0
        @test is_integer(transcription_variable(tm, x0))
        @test lower_bound(transcription_variable(tm, y0)) == 0
        @test upper_bound(transcription_variable(tm, y0)) == 1
        @test is_integer(transcription_variable(tm, y0))
        @test start_value(transcription_variable(tm, y0)) == 0.
        @test has_lower_bound(IOTO.lookup_by_support(tm, y, [0., 0., 0.]))
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
    tm = transcription_model(m)
    IOTO.set_parameter_supports(tm, m)
    IOTO.transcribe_finite_variables!(tm, m)
    IOTO.transcribe_infinite_variables!(tm, m)
    tx = transcription_variable(x)
    ty = transcription_variable(y)
    tz = transcription_variable(z)
    tw = transcription_variable(w)
    # test transcribe_measures!
    @testset "transcribe_measures!" begin 
        @test IOTO.transcribe_measures!(tm, m) isa Nothing 
        @test transcription_variable(meas1) == tx[1] + 4tw + tx[2]
        @test transcription_variable(meas2) == tw + 0
        @test transcription_variable(meas3) isa Vector
        @test transcription_variable(meas4) isa AffExpr
        @test supports(meas1) == ()
        @test supports(meas2) == ()
        @test sort!(supports(meas3)) == [(0.,), (1., )]
    end
    # test transcribe_objective!
    @testset "transcribe_objective!" begin 
        # normal
        @objective(m, Min, 2z^2 - meas1)
        @test IOTO.transcribe_objective!(tm, m) isa Nothing 
        @test objective_sense(tm) == MOI.MIN_SENSE
        @test objective_function(tm) == 2tz^2 - tx[1] - 4tw - tx[2]
        # nonlinear objective
        @objective(m, Max, z^4)
        @test IOTO.transcribe_objective!(tm, m) isa Nothing 
        @test objective_sense(tm) == MOI.MAX_SENSE
        @test objective_function_string(MIME("text/plain"), tm) == "subexpression[1] + 0.0"
        @test sprint(show, NonlinearExpression(tm, 1)) == "subexpression[1]: z ^ 4.0"
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
    tm = transcription_model(m)
    IOTO.set_parameter_supports(tm, m)
    IOTO.transcribe_finite_variables!(tm, m)
    IOTO.transcribe_infinite_variables!(tm, m)
    IOTO.transcribe_point_variables!(tm, m)
    xt = transcription_variable(x)
    yt = transcription_variable(y)
    x0t = transcription_variable(x0)
    yft = transcription_variable(yf)
    zt = transcription_variable(z)
    # test _get_info_constr_from_var
    @testset "_get_info_constr_from_var" begin 
        # lower bounds
        set = moi_set(InfiniteOpt._core_constraint_object(LowerBoundRef(x)))
        expected = LowerBoundRef(IOTO.lookup_by_support(tm, x, [1., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tm, x, set, [1., 1., 1.]) == expected
        @test IOTO._get_info_constr_from_var(tm, x, set, [0., 0., 0.]) isa Nothing
        # upper bounds
        set = moi_set(InfiniteOpt._core_constraint_object(UpperBoundRef(x)))
        expected = UpperBoundRef(IOTO.lookup_by_support(tm, x, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tm, x, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tm, x, set, [0., 0., 0.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(UpperBoundRef(yf)))
        expected = UpperBoundRef(yft)
        @test IOTO._get_info_constr_from_var(tm, y, set, [1., 1., 1.]) == expected
        # fix 
        set = moi_set(InfiniteOpt._core_constraint_object(FixRef(y)))
        expected = FixRef(IOTO.lookup_by_support(tm, y, [0.]))
        @test IOTO._get_info_constr_from_var(tm, y, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tm, y, set, [0., 0., 1.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(FixRef(x0)))
        expected = FixRef(x0t)
        @test IOTO._get_info_constr_from_var(tm, x, set, [0., 0., 0.]) == expected
        # binary 
        set = moi_set(InfiniteOpt._core_constraint_object(BinaryRef(x0)))
        expected = BinaryRef(IOTO.lookup_by_support(tm, x, [0., 0., 0.]))
        @test IOTO._get_info_constr_from_var(tm, x, set, [0., 0., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tm, x, set, [1., 1., 1.]) isa Nothing
        @test IOTO._get_info_constr_from_var(tm, z, set, [1., 1., 1.]) == BinaryRef(zt)
        # integer
        set = moi_set(InfiniteOpt._core_constraint_object(IntegerRef(x)))
        expected = IntegerRef(IOTO.lookup_by_support(tm, x, [0., 1., 1.]))
        @test IOTO._get_info_constr_from_var(tm, x, set, [1., 1., 0.]) == expected
        @test IOTO._get_info_constr_from_var(tm, x, set, [0., 0., 0.]) isa Nothing
        set = moi_set(InfiniteOpt._core_constraint_object(IntegerRef(yf)))
        expected = IntegerRef(yft)
        @test IOTO._get_info_constr_from_var(tm, yf, set, [1., 1., 1.]) == expected
    end
    # test _support_in_restrictions
    @testset "_support_in_restrictions" begin 
        @test IOTO._support_in_restrictions([0., 0., 0.], [3], [IntervalDomain(0, 0)])
        @test IOTO._support_in_restrictions([0., 0., 0.], Int[], IntervalDomain[])
        @test IOTO._support_in_restrictions([NaN, 0., 0.], [1], [IntervalDomain(0, 1)])
        @test !IOTO._support_in_restrictions([NaN, 0., 0.], [1, 2], [IntervalDomain(1, 1), IntervalDomain(1, 1)])
        @test !IOTO._support_in_restrictions([NaN, 0., 2.], [1, 3], [IntervalDomain(1, 1), IntervalDomain(1, 1)])
    end
    # test _make_constr_ast
    @testset "_make_constr_ast" begin 
        @test IOTO._make_constr_ast(xt, MOI.LessThan(1.0)) == :($xt <= 1.0)
        @test IOTO._make_constr_ast(xt, MOI.GreaterThan(1.0)) == :($xt >= 1.0)
        @test IOTO._make_constr_ast(xt, MOI.EqualTo(1.0)) == :($xt == 1.0)
        @test IOTO._make_constr_ast(xt, MOI.Interval(0.0, 1.0)) == :(0.0 <= $xt <= 1.0)
        @test_throws ErrorException IOTO._make_constr_ast(xt, MOI.Integer())
    end
    # test _process_constraint
    @testset "_process_constraint" begin
        # scalar constraint 
        con = constraint_object(c1)
        func = jump_function(con)
        set = moi_set(con)
        @test IOTO._process_constraint(tm, con, func, set, zeros(3), "test1") isa ConstraintRef 
        @test num_constraints(tm, typeof(func), typeof(set)) == 1
        cref = constraint_by_name(tm, "test1")
        delete(tm, cref)
        # nonlinear scalar constraint 
        con = constraint_object(c7)
        func = jump_function(con)
        set = moi_set(con)
        expected = Sys.iswindows() ? "subexpression[1] - 0.0 == 0" : "subexpression[1] - 0.0 = 0"
        @test sprint(show, IOTO._process_constraint(tm, con, func, set, zeros(3), "test1")) == expected
        expected = ["subexpression[1]: sin(z) ^ x(support: 1) - 0.0", 
                    "subexpression[1]: sin(z) ^ x(support: 2) - 0.0"]
        @test sprint(show, NonlinearExpression(tm, 1)) in expected
        tm.nlp_model = nothing
        # vector constraint 
        con = constraint_object(c6)
        func = jump_function(con)
        set = moi_set(con)
        @test IOTO._process_constraint(tm, con, func, set, zeros(3), "test2") isa ConstraintRef 
        @test num_constraints(tm, typeof(func), typeof(set)) == 1
        cref = constraint_by_name(tm, "test2")
        delete(tm, cref)
        # test nonlinear vector constraint 
        con = VectorConstraint([sin(z)], MOI.Zeros(1))
        func = [sin(z)]
        set = MOI.Zeros(1)
        @test_throws ErrorException IOTO._process_constraint(tm, con, func, set, zeros(3), "test2")
        tm.nlp_model = nothing
        # fallback
        @test_throws ErrorException IOTO._process_constraint(tm, :bad, func, set, 
                                                             zeros(3), "bad")
    end
    # test transcribe_constraints!
    @testset "transcribe_constraints!" begin 
        @test IOTO.transcribe_constraints!(tm, m) isa Nothing
        # test the info constraint transcriptions
        @test transcription_constraint(LowerBoundRef(x)) isa Vector{ConstraintRef}
        @test transcription_constraint(UpperBoundRef(x)) isa Vector{ConstraintRef}
        @test transcription_constraint(IntegerRef(x)) isa Vector{ConstraintRef}
        @test length(transcription_constraint(LowerBoundRef(x))) == 5
        @test transcription_constraint(FixRef(x0)) == FixRef(x0t)
        @test transcription_constraint(BinaryRef(x0)) == BinaryRef(x0t)
        @test transcription_constraint(FixRef(y)) == FixRef.(yt)[1:2]
        @test transcription_constraint(UpperBoundRef(yf)) == UpperBoundRef(yft)
        @test transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
        # test constraint transcriptions 
        @test transcription_constraint(c1) isa Vector{ConstraintRef}
        @test length(transcription_constraint(c1)) == 6
        @test constraint_object(transcription_constraint(c2)).func == yt[1] - zt^2
        xf = IOTO.lookup_by_support(tm, x, [1., 1., 1.])
        @test constraint_object(transcription_constraint(c3)).func == 4xf 
        @test constraint_object(transcription_constraint(c3)).set == MOI.LessThan(5.)
        expected = [yt[1] + zt, yt[2] + zt]
        @test jump_function.(constraint_object.(transcription_constraint(c4))) == expected
        @test constraint_object(transcription_constraint(c5)).func == 2zt^2 
        @test length(transcription_constraint(c6)) == 6
        @test moi_set(constraint_object(first(transcription_constraint(c6)))) == MOI.Zeros(2)
        @test length(transcription_constraint(c7)) == 6
        @test length(keys(tm.nlp_model.constraints)) == 6
        @test length(tm.nlp_model.expressions) == 6
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
    end
end

# Test transcribe_derivative_evaluations! (TODO SEE IF THIS CAN BE MADE TO WORK WITH DEPENDENDENT PARAMETER BASED DERIVATIVES)
@testset "transcribe_derivative_evaluations!" begin 
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1], num_supports = 4, 
                        derivative_method = OrthogonalCollocation(4))
    @variable(m, y, Infinite(t))
    d1 = @deriv(y, t)
    d2 = @deriv(y, t^2)
    tm = transcription_model(m)
    IOTO.set_parameter_supports(tm, m)
    IOTO.transcribe_infinite_variables!(tm, m)
    IOTO.transcribe_derivative_variables!(tm, m)
    # main test 
    @test IOTO.transcribe_derivative_evaluations!(tm, m) isa Nothing 
    @test num_constraints(tm, AffExpr, MOI.EqualTo{Float64}) == 18
    @test num_supports(t) == 4
    @test num_supports(t, label = All) == 10
    @test length(supports(d1)) == 4
    @test length(supports(d2, label = All)) == 10
    @test IOTO.has_internal_supports(tm)
    @test has_generative_supports(t)
end

# Test build_transcription_model!
@testset "build_transcription_model!" begin
    # initialize model
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    m = InfiniteModel(mockoptimizer)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
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
    @register(m, g(a))
    @constraint(m, c7, g(z) == 2)
    @objective(m, Min, x0 + meas1)
    # test basic usage
    tm = optimizer_model(m)
    @test IOTO.build_transcription_model!(tm, m, 
                                          check_support_dims = false) isa Nothing
    # test finite variables
    zt = transcription_variable(z)
    wt = transcription_variable(w)
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
    @test transcription_variable(x) isa Vector{VariableRef}
    @test transcription_variable(y) isa Vector{VariableRef}
    @test name(transcription_variable(x)[1]) == "x(support: 1)"
    @test name(transcription_variable(y)[3]) == "y(support: 3)"
    @test has_lower_bound(transcription_variable(x)[1])
    @test is_binary(transcription_variable(y)[2])
    @test is_fixed(transcription_variable(y)[4])
    @test is_integer(transcription_variable(x)[2])
    @test start_value(transcription_variable(y)[1]) == 0.
    @test supports(x) == [(0.,), (1.,)]
    @test length(supports(y)) == 4
    # test point variables
    @test transcription_variable(x0) isa VariableRef
    @test transcription_variable(y0) isa VariableRef
    @test name(transcription_variable(x0)) == "x(support: 1)"
    @test name(transcription_variable(y0))[1:end-2] == "y(support: "
    @test has_lower_bound(transcription_variable(x0))
    @test is_integer(transcription_variable(x0))
    @test has_lower_bound(transcription_variable(y0))
    @test has_upper_bound(transcription_variable(y0))
    @test is_integer(transcription_variable(y0))
    @test start_value(transcription_variable(y0)) == 0.
    # test derivatives 
    d1t = transcription_variable(d1)
    d2t = transcription_variable(d2)
    @test length(d1t) == 4
    @test length(d2t) == 2
    @test upper_bound(d1t[1]) == 2
    @test supports(d2) == [(0.,), (1.,)]
    # test registration 
    r = tm.nlp_model.operators
    @test length(r.registered_univariate_operators) == 1
    @test r.univariate_operator_f[1].f == g
    # test objective
    xt = transcription_variable(tm, x)
    @test objective_function(tm) == 2xt[1] + xt[2] - 2wt - d2t[1] - d2t[2]
    @test objective_sense(tm) == MOI.MIN_SENSE
    # test constraints
    yt = transcription_variable(y)
    dt_c1 = IOTO.lookup_by_support(tm, d1, zeros(3))
    @test constraint_object(transcription_constraint(c1)).func == -zt + xt[1] + dt_c1
    @test constraint_object(transcription_constraint(c2)).func == zt + xt[1]
    expected = transcription_variable(meas2)[2] - 2 * transcription_variable(y0) + xt[2]
    @test constraint_object(transcription_constraint(c4)).func == expected
    @test constraint_object(transcription_constraint(c3)).func == xt[1] - 2wt + xt[2] + zt - d2t[1] - d2t[2]
    @test constraint_object(transcription_constraint(c6)).func == [zt, wt]
    @test transcription_constraint(c5) isa Vector{ConstraintRef}
    @test name(transcription_constraint(c2)) == "c2(support: 1)"
    @test name(transcription_constraint(c1)) == "c1(support: 1)"
    @test supports(c1) == (0., [0., 0.])
    @test transcription_constraint(c7) isa NonlinearConstraintRef
    # test info constraints
    @test transcription_constraint(LowerBoundRef(z)) == LowerBoundRef(zt)
    @test transcription_constraint(UpperBoundRef(z)) == UpperBoundRef(zt)
    @test transcription_constraint(BinaryRef(z)) == BinaryRef(zt)
    @test transcription_constraint(FixRef(w)) == FixRef(wt)
    @test transcription_constraint(IntegerRef(w)) == IntegerRef(wt)
    @test transcription_constraint(LowerBoundRef(x0)) == LowerBoundRef(xt[1])
    @test transcription_constraint(UpperBoundRef(x0)) == UpperBoundRef(xt[1])
    @test transcription_constraint(IntegerRef(x0)) == IntegerRef(xt[1])
    @test transcription_constraint(LowerBoundRef(y0)) isa ConstraintRef
    @test transcription_constraint(UpperBoundRef(y0)) isa ConstraintRef
    @test transcription_constraint(IntegerRef(y0)) isa ConstraintRef
    @test transcription_constraint(LowerBoundRef(x)) == LowerBoundRef.(xt)
    @test transcription_constraint(UpperBoundRef(x)) == UpperBoundRef.(xt)
    @test transcription_constraint(IntegerRef(x)) == IntegerRef.(xt)
    @test transcription_constraint(FixRef(y)) isa Vector{ConstraintRef}
    @test transcription_constraint(BinaryRef(y)) isa Vector{ConstraintRef}
    @test transcription_constraint(UpperBoundRef(d1)) == UpperBoundRef.(d1t)

    # test a finite model 
    m = InfiniteModel()
    @variable(m, y >= 0)
    @objective(m, Min, y)
    tm = transcription_model(m)
    @test IOTO.build_transcription_model!(tm, m) isa Nothing 
    @test transcription_variable(y) isa VariableRef 
    @test lower_bound(transcription_variable(y)) == 0
    @test objective_sense(tm) == MOI.MIN_SENSE
end
