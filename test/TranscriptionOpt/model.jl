# Test the TranscriptionData datatype
@testset "TranscriptionData" begin
    @test TranscriptionData isa DataType
    @test TranscriptionData().finvar_mappings isa Dict
end

# Test basic definition and queries
@testset "Basic Definition and Queries" begin
    # initialize needed data
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    # test TranscriptionModel (no optimizer)
    @testset "TranscriptionModel (Default)" begin
        @test haskey(TranscriptionModel().ext, :TransData)
    end
    # test TranscriptionModel (with optimizer)
    @testset "TranscriptionModel (Optimizer)" begin
        @test haskey(TranscriptionModel(mockoptimizer).ext, :TransData)
    end
    # test is_transcription_model
    @testset "is_transcription_model" begin
        @test is_transcription_model(TranscriptionModel())
        @test !is_transcription_model(Model())
    end
    # test transcription_data
    @testset "transcription_data" begin
        @test transcription_data(TranscriptionModel()).reduced_vars isa Vector
        @test_throws ErrorException transcription_data(Model())
    end
end

# Test variable mapping queries
@testset "Variable Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @dependent_parameters(m, pars[1:2] in [0, 1])
    @infinite_variable(m, x(par, pars))
    @point_variable(m, x(0, [0, 0]), x0)
    @hold_variable(m, y)
    supps = Dict{Int, Float64}(2 => 1)
    xrv = add_variable(m, build_variable(error, x, supps))
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    # test transcription_variable (finite variables)
    @testset "transcription_variable (Finite)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, y)
        @test_throws ErrorException transcription_variable(tm, x0)
        # test normal
        tm.ext[:TransData].finvar_mappings[y] = a
        @test transcription_variable(tm, y) == a
        tm.ext[:TransData].finvar_mappings[x0] = b
        @test transcription_variable(tm, x0) == b
    end
    # test transcription_variable (Infinite and Reduced)
    @testset "transcription_variable (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, x)
        @test_throws ErrorException transcription_variable(tm, xrv)
        # test normal
        tm.ext[:TransData].infvar_mappings[x] = [b, c, d]
        @test transcription_variable(tm, x) == [b, c, d]
        tm.ext[:TransData].infvar_mappings[xrv] = [c, d]
        @test transcription_variable(tm, xrv) == [c, d]
    end
    # test transcription_variable (Fallback)
    @testset "transcription_variable (Fallback)" begin
        @test_throws ErrorException transcription_variable(tm, par)
    end
    # test transcription_variable (Single argument)
    @testset "transcription_variable (Single)" begin
        @test transcription_variable(y) == a
        @test transcription_variable(x) == [b, c, d]
        @test transcription_variable(x0) == b
    end
    # test optimizer_model_variable extension
    @testset "optimizer_model_variable" begin
        @test optimizer_model_variable(y, Val(:TransData)) == a
        @test optimizer_model_variable(x, Val(:TransData)) == [b, c, d]
        @test optimizer_model_variable(x0, Val(:TransData)) == b
    end
    # test variable_supports for infinite variable with 2 inputs
    @testset "variable_supports (Model, Infinite Variable)" begin
        # test not in mappings error
        dvref = dispatch_variable_ref(x)
        delete!(tm.ext[:TransData].infvar_mappings, x)
        @test_throws ErrorException InfiniteOpt.variable_supports(tm, dvref)
        tm.ext[:TransData].infvar_mappings[x] = [b, c, d]
        # test supports are empty
        lookups = Dict{Vector{Float64}, Int}([0, 0, 0] => 1, [0, 0, 1] => 2, [1, 1, 1] => 3)
        tm.ext[:TransData].infvar_lookup[x] = lookups
        @test InfiniteOpt.variable_supports(tm, dvref) == [(0., [0., 0.]), (0., [0., 1.]), (1., [1., 1.])]
        # test normal
        @test InfiniteOpt.variable_supports(tm, dvref) == [(0., [0., 0.]), (0., [0., 1.]), (1., [1., 1.])]
        # test with reduced variable
        lookups = Dict{Vector{Float64}, Int}([0, 1] => 1, [1, 1] => 2)
        tm.ext[:TransData].infvar_lookup[xrv] = lookups
        dvref = dispatch_variable_ref(xrv)
        @test InfiniteOpt.variable_supports(tm, dvref) == [(0., 1.), (1., 1.)]
    end
    # test supports for infinite variable
    @testset "supports (Infinite)" begin
        @test supports(x) == [(0., [0., 0.]), (0., [0., 1.]), (1., [1., 1.])]
        @test supports(xrv) == [(0., 1.), (1., 1.)]
    end
    # test lookup_by_support (infinite vars)
    @testset "lookup_by_support (Infinite)" begin
        # test errors
        @infinite_variable(m, x2(par))
        @test_throws ErrorException IOTO.lookup_by_support(tm, x2, [0.])
        @test_throws ErrorException IOTO.lookup_by_support(tm, x, [1., 0., 0.])
        @test_throws ErrorException IOTO.lookup_by_support(tm, xrv, [0., 0., 0.])
        # test normal
        @test IOTO.lookup_by_support(tm, x, [0., 0., 0.]) == b
        @test IOTO.lookup_by_support(tm, x, [0., 0., 1.]) == c
        @test IOTO.lookup_by_support(tm, x, [1., 1., 1.]) == d
        @test IOTO.lookup_by_support(tm, xrv, [0., 0., 1.]) == c
        @test IOTO.lookup_by_support(tm, xrv, [1., 1., 1.]) == d
    end
    # test lookup_by_support (finite vars)
    @testset "lookup_by_support (Finite)" begin
        # test errors
        @hold_variable(m, z2)
        @test_throws ErrorException IOTO.lookup_by_support(tm, z2, [0.])
        # test normal
        @test IOTO.lookup_by_support(tm, x0, [0., 0., 0.]) == b
        @test IOTO.lookup_by_support(tm, y, [0., 0., 1.]) == a
    end
    # test internal_reduced_variable
    @testset "internal_reduced_variable" begin
        rv = ReducedVariableRef(m, ReducedVariableIndex(-1))
        # test errors
        @test_throws ErrorException InfiniteOpt.internal_reduced_variable(rv, Val(:TransData))
        # test normal
        var = build_variable(error, x, supps)
        push!(transcription_data(tm).reduced_vars, var)
        @test InfiniteOpt.internal_reduced_variable(rv, Val(:TransData)) == var
        eval_supports(rv) == supps
    end
end

# Test the measure mapping methods
@testset "Measure Queries" begin
    # initialize tbe needed info
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @dependent_parameters(m, pars[1:2] in [0, 1])
    @infinite_variable(m, x(par, pars))
    @point_variable(m, x(0, [0, 0]), x0)
    @hold_variable(m, y)
    meas1 = support_sum(2par -2, par)
    meas2 = support_sum(x^2 - y, pars)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    # test transcription_variable
    @testset "transcription_variable" begin
        # test errors
        @test_throws ErrorException IOTO.transcription_variable(tm, meas1)
        @test_throws ErrorException IOTO.transcription_variable(tm, meas2)
        # test normal
        transcription_data(tm).measure_mappings[meas1] = [-2 * zero(AffExpr)]
        transcription_data(tm).measure_mappings[meas2] = [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
        @test IOTO.transcription_variable(tm, meas1) == -2 * zero(AffExpr)
        @test IOTO.transcription_variable(tm, meas2) == [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
    end
    # test lookup_by_support
    @testset "lookup_by_support" begin
        # test errors
        @test_throws ErrorException IOTO.lookup_by_support(tm, meas1, Float64[])
        @test_throws ErrorException IOTO.lookup_by_support(tm, meas2, [0.])
        # test normal
        transcription_data(tm).measure_lookup[meas1] = Dict(Float64[] => 1)
        transcription_data(tm).measure_lookup[meas2] = Dict{Vector{Float64}, Int}([0] => 1, [1] => 2)
        @test IOTO.lookup_by_support(tm, meas1, Float64[]) == -2 * zero(AffExpr)
        @test IOTO.lookup_by_support(tm, meas2, [1.]) == b^2 + d^2 - 2a
    end
    # test InfiniteOpt.variable_supports
    @testset "InfiniteOpt.variable_supports" begin
        # test not in mappings error
        dvref = dispatch_variable_ref(meas1)
        delete!(tm.ext[:TransData].measure_mappings, meas1)
        @test_throws ErrorException InfiniteOpt.variable_supports(tm, dvref)
        tm.ext[:TransData].measure_mappings[meas1] = [-2 * zero(AffExpr)]
        # test supports are empty
        @test InfiniteOpt.variable_supports(tm, dvref) == ()
        dvref2 = dispatch_variable_ref(meas2)
        @test InfiniteOpt.variable_supports(tm, dvref2) == [(0.,), (1.,)]
        # test normal
        @test InfiniteOpt.variable_supports(tm, dvref) == ()
        @test InfiniteOpt.variable_supports(tm, dvref2) == [(0.,), (1.,)]
    end
end

# Test the support iteration methods 
@testset "Support Iteration Methods" begin 
    # initialize the needed info 
    m = InfiniteModel()
    @dependent_parameters(m, pars[1:2] in [0, 1], num_supports = 2)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    tm = optimizer_model(m)
    # test _temp_parameter_ref
    @testset "_temp_parameter_ref" begin 
        @test IOTO._temp_parameter_ref(m, index(par)) == dispatch_variable_ref(par)
        @test IOTO._temp_parameter_ref(m, index(pars[1]).object_index) == dispatch_variable_ref(pars[1])
    end
    # test _collected_supports
    @testset "_collected_supports" begin 
        # independent parameters
        expected = Float64[0., 1., NaN]
        @test isequal(IOTO._collected_supports(m, index(par)), expected)
        # dependent parameters
        expected = [Float64[0., 0.], Float64[1., 1.], Float64[NaN, NaN]]
        @test isequal(IOTO._collected_supports(m, index(pars[1]).object_index), expected)
    end
    # test set_parameter_supports
    @testset "set_parameter_supports" begin 
        @test IOTO.set_parameter_supports(tm, m) isa Nothing 
        expected = ([[0., 0.], [1., 1.], [NaN, NaN]], [0., 1., NaN])
        @test isequal(transcription_data(tm).supports, expected)
    end
    # test parameter_supports
    @testset "parameter_supports" begin 
        expected = ([[0., 0.], [1., 1.], [NaN, NaN]], [0., 1., NaN])
        @test isequal(IOTO.parameter_supports(tm), expected)
    end
    # test support_index_iterator with 1 argument
    @testset "support_index_iterator (1 Arg)" begin 
        @test IOTO.support_index_iterator(tm) == CartesianIndices((1:2, 1:2))
    end
    # test support_index_iterator with 2 argument2
    @testset "support_index_iterator (2 Args)" begin 
        @test IOTO.support_index_iterator(tm, [1, 2]) == CartesianIndices((1:2, 1:2))
        @test IOTO.support_index_iterator(tm, [1]) == CartesianIndices((1:2, 3:3))
    end
    # test index_to_support
    @testset "index_to_support" begin 
        @test IOTO.index_to_support(tm, first(CartesianIndices((1:2, 1:2)))) == zeros(3)
        expected = [1., 1., NaN]
        @test isequal(IOTO.index_to_support(tm, last(IOTO.support_index_iterator(tm, [1]))), expected)
    end
end

# Test the expression mappings
@testset "Expression Queries" begin
    # initialize tbe needed info
    m = InfiniteModel()
    @dependent_parameters(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par, pars))
    @finite_parameter(m, finpar, 42)
    @point_variable(m, x(0, [0, 0]), x0)
    @hold_variable(m, y)
    meas1 = support_sum(2par -2, par)
    meas2 = support_sum(x^2 - y, pars)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    # transcribe the variables and measures
    data = transcription_data(tm)
    data.finvar_mappings[y] = a
    data.finvar_mappings[x0] = b
    data.infvar_mappings[x] = [b, c, d]
    data.measure_mappings[meas1] = [-2 * zero(AffExpr)]
    data.measure_mappings[meas2] = [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
    lookups = Dict{Vector{Float64}, Int}([0, 0, 0] => 1, [0, 0, 1] => 2, [1, 1, 1] => 3)
    data.infvar_lookup[x] = lookups
    data.measure_lookup[meas1] = Dict(Float64[] => 1)
    data.measure_lookup[meas2] = Dict{Vector{Float64}, Int}([0] => 1, [1] => 2)
    @test IOTO.set_parameter_supports(tm, m) isa Nothing
    # test transcription_expression in accordance with the methods defined in transcribe.jl
    @testset "transcription_expression (Fallback)" begin
        @test_throws ErrorException IOTO.transcription_expression(tm, a, Float64[])
    end
    # test transcription expression for infinite variables with 3 args
    @testset "transcription_expression (Infinite Variable)" begin
        @test IOTO.transcription_expression(tm, x, [0., 1., 0.]) == c
        @test IOTO.transcription_expression(tm, meas1, [0., 0., 1.]) == -2 * zero(AffExpr)
    end
    # test transcription expression for finite variables with 3 args
    @testset "transcription_expression (Finite Variable)" begin
        @test IOTO.transcription_expression(tm, x0, [0., 1., 0.]) == b
        @test IOTO.transcription_expression(tm, y, [0., 0., 1.]) == a
    end
    # test transcription expression for infinite parameters with 3 args
    @testset "transcription_expression (Infinite Parameter)" begin
        @test IOTO.transcription_expression(tm, par, [0.5, 1., 0.]) == 0.
        @test IOTO.transcription_expression(tm, pars[1], [0.5, 0., 1.]) == 0.5
    end
    # test transcription expression for finite parameters with 3 args
    @testset "transcription_expression (Finite Parameter)" begin
        @test IOTO.transcription_expression(tm, finpar, [0.5, 1., 0.]) == 42
    end
    # test transcription expression for AffExprs with 3 args
    @testset "transcription_expression (AffExpr)" begin
        @test IOTO.transcription_expression(tm, x0 - y + 2x - 2.3, [1., 1., 1.]) == b - a + 2d - 2.3
    end
    # test transcription expression for QuadExprs with 3 args
    @testset "transcription_expression (QuadExpr)" begin
        expr = meas2 - 3y^2 - x0 - 2.3
        expected = b^2 + d^2 - 2a - 3a^2 - b - 2.3
        @test IOTO.transcription_expression(tm, expr, [1., 1., 1.]) == expected
    end
    # test transcription expression for Exprs with 2 args
    @testset "transcription_expression (Expr 2 Args)" begin
        # test infinite expr
        expr = meas2 - 3y^2 - x0 - 2.3
        expected = [-2a^2 + c^2 - 2a - b - 2.3, -3a^2 + b^2 + d^2 - 2a - b - 2.3]
        @test IOTO.transcription_expression(tm, expr) == expected
        # test finite expr 
        expr = 2x0 -y 
        expected = 2b- a 
        @test IOTO.transcription_expression(tm, expr) == expected
    end
    # test transcription expression for variables with 2 args
    @testset "transcription_expression (Variable 2 Args)" begin
        @test IOTO.transcription_expression(tm, x0) == b 
        @test IOTO.transcription_expression(tm, x) == [b, c, d] 
    end
    # test transcription expression with 1 argument
    @testset "transcription_expression (1 Arg)" begin
        @test IOTO.transcription_expression(x0) == b 
        @test IOTO.transcription_expression(x0 - y) == b - a 
        @test IOTO.transcription_expression(zero(QuadExpr) + 2) == zero(AffExpr) + 2
    end
    # test optimizer_model_expression
    @testset "optimizer_model_expression" begin
        @test optimizer_model_expression(x0) == b 
        @test optimizer_model_expression(x0 - y) == b - a 
        @test optimizer_model_expression(zero(QuadExpr) + 2) == zero(AffExpr) + 2
    end
    # test expression_supports
    @testset "expression_supports" begin 
        @test supports(x0) == ()
        @test supports(x0 - y) == ()
        @test supports(meas2 + y) == [(0.,), (1.,)]
    end
end


# Test constraint mapping queries
@testset "Constraint Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_variable(m, x(par))
    @point_variable(m, x(0), x0)
    @hold_variable(m, y)
    @constraint(m, c1, x + y - 2 <= 0)
    @constraint(m, c2, support_sum(x, par) == 0)
    @constraint(m, c3, x0 + y == 5)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @constraint(tm, tc1, b + a <= 2)
    @constraint(tm, tc2, c + a <= 2)
    @constraint(tm, tc3, 0.5b + 0.5c == 0)
    @constraint(tm, tc4, b + a == 5)
    # test transcription_constraint (2 Args)
    @testset "transcription_constraint (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_constraint(tm, c1)
        # test normal
        transcription_data(tm).constr_mappings[c1] = [tc1, tc2]
        @test transcription_constraint(tm, c1) == [tc1, tc2]
        # test error
        @test_throws ErrorException transcription_constraint(tm, c2)
        # test normal
        transcription_data(tm).constr_mappings[c2] = [tc3]
        @test transcription_constraint(tm, c2) == tc3
        # test error
        @test_throws ErrorException transcription_constraint(tm, c3)
        # test normal
        transcription_data(tm).constr_mappings[c3] = [tc4]
        @test transcription_constraint(tm, c3) == tc4
    end
    # test transcription_constraint (Single argument)
    @testset "transcription_constraint (1 Arg)" begin
        @test transcription_constraint(c1) == [tc1, tc2]
        @test transcription_constraint(c2) == tc3
        @test transcription_constraint(c3) == tc4
    end
    # test optimizer_model_constraint extension
    @testset "optimizer_model_constraint" begin
        @test optimizer_model_constraint(c1, Val(:TransData)) == [tc1, tc2]
        @test optimizer_model_constraint(c2, Val(:TransData)) == tc3
        @test optimizer_model_constraint(c3, Val(:TransData)) == tc4
    end
    # test constraint_supports 
    @testset "constraint_supports" begin
        # test error
        @test_throws ErrorException InfiniteOpt.constraint_supports(tm, c1)
        # test normal
        transcription_data(tm).constr_supports[c1] = [(0.,), (1.,)]
        @test InfiniteOpt.constraint_supports(tm, c1) == [(0.,), (1.,)]
        # test finite 
        transcription_data(tm).constr_supports[c3] = [()]
        @test InfiniteOpt.constraint_supports(tm, c3) == ()
    end
    # test supports
    @testset "supports" begin
        @test supports(c1) == [(0.,), (1.,)]
        @test supports(c3) == ()
    end
end
