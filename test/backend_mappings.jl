# Test build_transformation_backend!
@testset "build_transformation_backend!" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], 
                        derivative_method = OrthogonalCollocation(3))
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @variable(m, 1 >= x >= 0, Infinite(par), Int)
    @variable(m, y == 2, Infinite(par, pars), Bin, start = 0)
    @variable(m, x0, Point(x, 0))
    @variable(m, 0 <= y0 <= 1, Point(y, 0, [0, 0]), Int)
    @variable(m, 0 <= z <= 1, Bin)
    @variable(m, w == 1, Int, start = 1)
    @deriv(x, par)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, DomainRestrictions(par => [0.5, 1]))
    @constraint(m, c5, meas2 == 0)
    @objective(m, Min, x0 + meas1)
    # test extra keywords
    @test_throws ErrorException build_transformation_backend!(m, bad = 42)
    # test normal usage
    @test isa(build_transformation_backend!(m), Nothing)
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 14
    # test repeated build
    @test isa(build_transformation_backend!(m), Nothing)
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 14
    # test finite model
    m = InfiniteModel()
    @variable(m, y >= 0)
    @objective(m, Min, y)
    warn = "Finite models (i.e., `InfiniteModel`s with no infinite " * 
           "parameters) should be modeled directly via a `Model` in JuMP.jl."
    @test_logs (:warn, warn) build_transformation_backend!(m)
end

# Test transformation backend querying methods
@testset "Transformation Backend Queries" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], 
                        derivative_method = OrthogonalCollocation(3))
    @variable(m, x, Infinite(par))
    @variable(m, x0, Point(x, 0))
    @variable(m, z)
    d1 = @deriv(x, par)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - z, data1)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    f = parameter_function(sin, par)
    build_transformation_backend!(m)
    tb = m.backend
    tdata = IOTO.transcription_data(tb)
    # Test transformation_variable
    @testset "transformation_variable" begin
        # test normal usage
        @test transformation_variable(x, label = All) == IOTO.transcription_variable(x, label = All)
        @test transformation_variable(x, label = All, ndarray = true) == IOTO.transcription_variable(x, label = All, ndarray = true)
        @test transformation_variable(x0) == IOTO.transcription_variable(x0)
        @test transformation_variable(z) == IOTO.transcription_variable(z)
        @test transformation_variable(d1, label = InternalLabel) == IOTO.transcription_variable(d1, label = InternalLabel)
        @test transformation_variable(f) == [0, sin(1)]
        # test deprecation 
        @test (@test_deprecated optimizer_model_variable(z)) == transformation_variable(z)
        # test fallback
        @test_throws ErrorException transformation_variable(x, TestBackend(), my_key = true)
    end
    # Test variable_supports
    @testset "variable_supports" begin
        # test finite variable fallback
        @test InfiniteOpt.variable_supports(dispatch_variable_ref(z), TestBackend()) == ()
        @test InfiniteOpt.variable_supports(dispatch_variable_ref(x0), TestBackend(), label = All) == ()
        # test fallback
        @test_throws ErrorException InfiniteOpt.variable_supports(x, TestBackend())
    end
    # Test supports (variables)
    @testset "supports (Variables)" begin
        # test normal usage
        @test supports(x) == [(0.,), (1.,)]
        @test supports(x, ndarray = true) == [(0.,), (1.,)]
        @test supports(x, label = All) == [(0.,), (0.5,), (1.,)]
        @test supports(meas1) == () 
        @test supports(d1, label = InternalLabel) == [(0.5,)]
        @test supports(f, label = All) == [(0.,), (0.5,), (1.,)]
    end
    # Test transformation_expression
    @testset "transformation_expression" begin
        # test variable references
        @test transformation_expression(x, label = All) == IOTO.transcription_variable(x, label = All)
        @test transformation_expression(z) == IOTO.transcription_variable(z)
        @test transformation_expression(x0) == IOTO.transcription_variable(x0)
        @test transformation_expression(x0, ndarray = true) == IOTO.transcription_variable(x0)
        # test expression without variables 
        expr = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) + 42
        @test transformation_expression(expr) == zero(AffExpr) + 42
        # test normal expressions 
        xt = IOTO.transcription_variable(x, label = All)
        zt = IOTO.transcription_variable(z)
        @test transformation_expression(x^2 + z) == [xt[1]^2 + zt, xt[3]^2 + zt]
        @test transformation_expression(x^2 + z, ndarray = true) == [xt[1]^2 + zt, xt[3]^2 + zt]
        @test transformation_expression(x^2 + z, label = All) == [xt[1]^2 + zt, xt[2]^2 + zt, xt[3]^2 + zt]
        @test transformation_expression(2z - 3) == 2zt - 3
        @test transformation_expression(2 * f) == [zero(AffExpr), zero(AffExpr) + sin(1) * 2]
        # test deprecation
        @test (@test_deprecated optimizer_model_expression(2z-4)) == 2zt - 4
        # test fallback
        @test_throws ErrorException transformation_expression(c1, TestBackend(), my_key = true)
    end
    # Test expression_supports
    @testset "expression_supports" begin
        # test normal usage
        @test InfiniteOpt.expression_supports(x, tb) == [(0.,), (1.,)]
        @test InfiniteOpt.expression_supports(x^2 - x0, tb) == [(0.,), (1.,)]
        @test InfiniteOpt.expression_supports(x^2 - x0, tb, label = All) == [(0.,), (0.5,), (1.,)]
        # test fallback
        @test_throws ErrorException InfiniteOpt.expression_supports(z^2, TestBackend())
    end
    # Test supports (expressions)
    @testset "supports (Expressions)" begin
        # test normal usage
        @test supports(x) == [(0.,), (1.,)]
        @test supports(2x - z + x0 + 43) == [(0.,), (1.,)]
        @test supports(2x - z + x0 + 43, ndarray = true) == [(0.,), (1.,)]
        expr = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef}) + 42
        @test supports(expr, label = All) == ()
    end
    # Test transformation_constraint
    @testset "transformation_constraint" begin
        # test normal usage
        @test transformation_constraint(c1) == IOTO.transcription_constraint(c1)
        @test transformation_constraint(c2, label = All) == IOTO.transcription_constraint(c2, label = All)
        @test transformation_constraint(c2, label = All, ndarray = true) == IOTO.transcription_constraint(c2, label = All, ndarray = true)
        @test transformation_constraint(c3) == IOTO.transcription_constraint(c3)
        # test deprecation 
        @test (@test_deprecated optimizer_model_constraint(c1)) == transformation_constraint(c1)
        # test fallback
        @test_throws ErrorException transformation_constraint(c1, TestBackend(), my_key = true)
    end
    # Test constraint_supports
    @testset "constraint_supports" begin
        # test normal usage
        @test InfiniteOpt.constraint_supports(c1, tb) == [(0.,), (1.,)]
        @test InfiniteOpt.constraint_supports(c1, tb, label = All) == [(0.,), (0.5,), (1.,)]
        @test InfiniteOpt.constraint_supports(c1, tb, label = All, ndarray = true) == [(0.,), (0.5,), (1.,)]
        # test fallback
        @test_throws ErrorException InfiniteOpt.constraint_supports(c1, TestBackend())
    end
    # Test supports (constraints)
    @testset "supports (Constraints)" begin
        # test normal usage
        @test supports(c1, label = UserDefined) == [(0.,), (1.,)]
        # test fallback
        @test supports(c2) == ()
    end
end

# Test optimize!
@testset "JuMP.optimize!" begin
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
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, DomainRestrictions(par => [0.5, 1]))
    @constraint(m, c5, meas2 == 0)
    @objective(m, Min, x0 + meas1)
    # test normal usage
    @test isa(optimize!(m, check_support_dims = false), Nothing)
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 8
    # test optimize hook
    function myhook(model; n = "", ub = 2, kwargs...)
        if !isempty(n)
            var = variable_by_name(model, n)
            set_upper_bound(var, ub)
        end
        optimize!(model; ignore_optimize_hook = true, kwargs...)
        return
    end
    @test set_optimize_hook(m, myhook) isa Nothing
    @test optimize!(m, n = "x", check_support_dims = false) isa Nothing
    @test transformation_backend_ready(m)
    @test num_variables(m.backend.model) == 8
    @test upper_bound(x) == 2
    @test set_optimize_hook(m, nothing) isa Nothing
    @test isnothing(m.optimize_hook)
    @test_throws ErrorException optimize!(m, n = "x")
    @test optimize!(m) isa Nothing
    @test transformation_backend_ready(m)
    set_transformation_backend(m, TestBackend())
    @test_throws ErrorException optimize!(m)
    @test_throws ErrorException optimize!(TestBackend())
end

# Test that we avoid world age problems with generated parameter functions 
@testset "Generated Parameter Function Build" begin 
    function make_model()
        m =InfiniteModel()
        @infinite_parameter(m, t in [0, 1])
        @variable(m, y, Infinite(t))
        myfunc(ts, a) = ts + a
        @parameter_function(m, d[i = 1:3] == (t) -> myfunc(t, i))
        @constraint(m, [i = 1:3], y >= d[i])
        build_transformation_backend!(m)
        return m
    end
    @test make_model() isa InfiniteModel
end