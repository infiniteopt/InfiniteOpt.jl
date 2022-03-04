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
        @test transcription_data(TranscriptionModel()).semi_infinite_vars isa Vector
        @test_throws ErrorException transcription_data(Model())
    end
    # test has_internal_supports
    @testset "has_internal_supports" begin
        m = TranscriptionModel()
        @test !IOTO.has_internal_supports(m)
        transcription_data(m).has_internal_supports = true
        @test IOTO.has_internal_supports(m)
    end
end

# Test query helper functions 
@testset "Query Formatters" begin 
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, x, Infinite(par, pars))
    @variable(m, q, Infinite(pars, par))
    @variable(m, w, Infinite(par))
    @variable(m, x0, Point(x, 0, [0, 0]))
    @variable(m, y)
    add_supports(par, 0.5, label = InternalLabel)
    @variable(m, xrv, SemiInfinite(x, par, [1, pars[2]]))
    tm = optimizer_model(m)
    data = transcription_data(tm)
    data.has_internal_supports = true
    data.supports = ([0., 0.5, 1., NaN], [[0., 0.], [1., 1.], [NaN, NaN]])
    s1 = Set([UserDefined])
    s2 = Set([InternalLabel])
    sd = Set{DataType}()
    data.support_labels = ([s1, s2, s1, sd], [s1, s1, sd])
    # test _get_array_type
    @testset "_get_array_type" begin 
        @test IOTO._get_array_type(ones(Bool, 2)) == Bool
    end
    # test _get_object_numbers
    @testset "_get_object_numbers" begin 
        @test IOTO._get_object_numbers(x) == [1, 2]
        @test IOTO._get_object_numbers(y) == []
        @test IOTO._get_object_numbers(y + x) == [1, 2]
    end
    # test make_ndarray
    @testset "make_ndarray" begin 
        # test finite variable 
        @test IOTO.make_ndarray(tm, y, [1], PublicLabel) == [1]
        # test ordered infinite variable
        @test IOTO.make_ndarray(tm, x, collect(1:6), All) == [1 4; 2 5; 3 6]
        @test IOTO.make_ndarray(tm, x, collect(1:6), PublicLabel) == [1 4; 3 6]
        # test unordered infinite variable 
        @test IOTO.make_ndarray(tm, q, collect(1:6), All) == [1 2 3; 4 5 6]
        @test IOTO.make_ndarray(tm, q, collect(1:6), PublicLabel) == [1 3; 4 6]
        # test infinite variable with single parameter 
        @test IOTO.make_ndarray(tm, w, collect(1:3), All) == [1, 2, 3]
        # test expression 
        @test IOTO.make_ndarray(tm, w + x, collect(1:6), All) == [1 4; 2 5; 3 6]
        # test error 
        @test_throws ErrorException IOTO.make_ndarray(tm, q, collect(1:3), All)
    end
end

# Test variable mapping queries
@testset "Variable Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, x, Infinite(par, pars))
    @variable(m, x0, Point(x, 0, [0, 0]))
    @variable(m, y)
    f1 = parameter_function(sin, par)
    f2 = parameter_function((a, b) -> 1, (par, pars))
    add_supports(par, 0.5, label = InternalLabel)
    supps = Dict{Int, Float64}(2 => 1)
    xrv = add_variable(m, build_variable(error, x, supps))
    tm = optimizer_model(m)
    data = transcription_data(tm)
    data.has_internal_supports = true
    data.supports = ([0., 0.5, 1., NaN], [[0., 0.], [1., 1.], [NaN, NaN]])
    s1 = Set([UserDefined])
    s2 = Set([InternalLabel])
    s3 = union(s1, s2)
    sd = Set{DataType}()
    data.support_labels = ([s1, s2, s1, sd], [s1, s1, sd])
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    @variable(tm, e)
    @variable(tm, f)
    # test _ignore_label
    @testset "_ignore_label" begin 
        @test IOTO._ignore_label(tm, All)
        @test !IOTO._ignore_label(tm, PublicLabel)
        @test !IOTO._ignore_label(tm, UserDefined)
    end
    # test transcription_variable (finite variables)
    @testset "transcription_variable (Finite)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, y)
        @test_throws ErrorException transcription_variable(tm, x0)
        # test normal
        data.finvar_mappings[y] = a
        @test transcription_variable(tm, y) == a
        data.finvar_mappings[x0] = b
        @test transcription_variable(tm, x0, ndarray = true, label = All) == b
    end
    # test transcription_variable (Infinite, semi-infinite, and derivative)
    @testset "transcription_variable (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, x)
        @test_throws ErrorException transcription_variable(tm, xrv)
        # test normal
        data.infvar_mappings[x] = [a, b, c, d, e, f]
        data.infvar_support_labels[x] = [s1, s3, s1, s1, s3, s1]
        @test transcription_variable(tm, x) == [a, b, c, d, e, f]
        @test transcription_variable(tm, x, label = All) == [a, b, c, d, e, f]
        @test transcription_variable(tm, x, label = InternalLabel) == [b, e]
        data.infvar_mappings[xrv] = [d, e, f]
        data.infvar_support_labels[xrv] = [s1, s3, s1]
        @test transcription_variable(tm, xrv) == [d, e, f]
        @test transcription_variable(tm, xrv, label = All) == [d, e, f]
        @test transcription_variable(tm, xrv, label = InternalLabel) == [e]
        # test ndarray 
        @test transcription_variable(tm, x, label = All, ndarray = true) == [a d; b e; c f]
        @test transcription_variable(tm, x, ndarray = true) == [a d; c f]
        @test_throws ErrorException transcription_variable(tm, xrv, ndarray = true)
    end
    # test transcription_variable (Parameter Function)
    @testset "transcription_variable (Parameter Function)" begin
        # test normal
        @test transcription_variable(tm, f1) == [sin(0), sin(1)]
        @test transcription_variable(tm, f1, label = All) == sin.([0, 0.5, 1])
        # test ndarray 
        @test transcription_variable(tm, f1, label = All, ndarray = true) == sin.([0, 0.5, 1])
        @test transcription_variable(tm, f2, ndarray = true) == ones(2, 2)
    end
    # test transcription_variable (Fallback)
    @testset "transcription_variable (Fallback)" begin
        @test_throws ErrorException transcription_variable(tm, par)
    end
    # test transcription_variable (Single argument)
    @testset "transcription_variable (Single)" begin
        @test transcription_variable(y) == a
        @test transcription_variable(x, label = All) == [a, b, c, d, e, f]
        @test transcription_variable(x0) == b
        @test transcription_variable(f2, ndarray = true) == ones(2, 2)
    end
    # test optimizer_model_variable extension
    @testset "optimizer_model_variable" begin
        @test optimizer_model_variable(y, Val(:TransData), label = All) == a
        @test optimizer_model_variable(x, Val(:TransData), label = All) == [a, b, c, d, e, f]
        @test optimizer_model_variable(x, Val(:TransData)) == [a, b, c, d, e, f]
        @test optimizer_model_variable(x0, Val(:TransData)) == b
    end
    # test variable_supports for infinite variable with 2 inputs
    @testset "variable_supports (Model, Infinite)" begin
        # test not in mappings error
        dvref = dispatch_variable_ref(x)
        delete!(data.infvar_mappings, x)
        @test_throws ErrorException InfiniteOpt.variable_supports(tm, dvref)
        data.infvar_mappings[x] = [a, b, c, d, e, f]
        # test supports are empty
        lookups = Dict{Vector{Float64}, Int}([0, 0, 0] => 1, [0.5, 0, 0] => 2, [1, 0, 0] => 3,
                                             [0, 1, 1] => 4, [0.5, 1, 1] => 5, [1, 1, 1] => 6)
        data.infvar_lookup[x] = lookups
        expected = [(0., [0., 0.]), (0.5, [0., 0.]), (1., [0., 0.]), 
                    (0., [1., 1.]), (0.5, [1., 1.]), (1., [1., 1.])]
        @test InfiniteOpt.variable_supports(tm, dvref) == expected
        # test normal
        @test InfiniteOpt.variable_supports(tm, dvref, label = All) == expected
        @test InfiniteOpt.variable_supports(tm, dvref, label = InternalLabel) == expected[[2, 5]]
        # test ndarray 
        expected = permutedims([(0., [0., 0.]) (0.5, [0., 0.]) (1., [0., 0.]); 
                                (0., [1., 1.]) (0.5, [1., 1.]) (1., [1., 1.])], (2, 1))
        @test InfiniteOpt.variable_supports(tm, dvref, ndarray = true) == expected[[1, 3], :]
        @test InfiniteOpt.variable_supports(tm, dvref, ndarray = true, label = All) == expected
        # test with semi-infinite variable
        lookups = Dict{Vector{Float64}, Int}([0, 1] => 1, [0.5, 1] => 2, [1, 1] => 3)
        data.infvar_lookup[xrv] = lookups
        dvref = dispatch_variable_ref(xrv)
        expected = [(0., 1.), (0.5, 1.), (1., 1.)]
        @test InfiniteOpt.variable_supports(tm, dvref) == expected
        @test InfiniteOpt.variable_supports(tm, dvref, label = InternalLabel) == [expected[2]]
    end
    # test variable_supports for infinite parameter functions with 2 inputs
    @testset "variable_supports (Model, Parameter Function)" begin
        # test normal
        df1 = dispatch_variable_ref(f1)
        @test InfiniteOpt.variable_supports(tm, df1) == [(0.,), (1.,)]
        @test InfiniteOpt.variable_supports(tm, df1, label = All) == [(0.,), (0.5,), (1.,)]
        # test ndarray 
        df2 = dispatch_variable_ref(f2)
        @test InfiniteOpt.variable_supports(tm, df1, label = All, ndarray = true) == [(0.,), (0.5,), (1.,)]
        @test InfiniteOpt.variable_supports(tm, df2, ndarray = true) isa Array
    end
    # test supports for infinite variable
    @testset "supports (Infinite)" begin
        @test supports(x, label = InternalLabel) == [(0.5, [0., 0.]), (0.5, [1., 1.])]
        @test supports(xrv) == [(0., 1.), (0.5, 1.), (1., 1.)]
        @test supports(f1, label = All) == [(0.,), (0.5,), (1.,)]
    end
    # test lookup_by_support (infinite vars)
    @testset "lookup_by_support (Infinite)" begin
        # test errors
        @variable(m, x2, Infinite(par))
        @test_throws ErrorException IOTO.lookup_by_support(tm, x2, [0.])
        @test_throws ErrorException IOTO.lookup_by_support(tm, x, [0.9, 0., 0.])
        @test_throws ErrorException IOTO.lookup_by_support(tm, xrv, [0., 0., 0.])
        # test normal
        @test IOTO.lookup_by_support(tm, x, [0., 0., 0.]) == a
        @test IOTO.lookup_by_support(tm, x, [0., 1., 1.]) == d
        @test IOTO.lookup_by_support(tm, x, [1., 1., 1.]) == f
        @test IOTO.lookup_by_support(tm, xrv, [0., 1.]) == d
        @test IOTO.lookup_by_support(tm, xrv, [1., 1.]) == f
    end
    # test lookup_by_support (infinite parameter functions)
    @testset "lookup_by_support (Parameter Function)" begin
        @test IOTO.lookup_by_support(tm, f1, [0.]) == 0
        @test IOTO.lookup_by_support(tm, f2, [0., 0., 1.]) == 1
    end
    # test lookup_by_support (finite vars)
    @testset "lookup_by_support (Finite)" begin
        # test errors
        @variable(m, z2)
        @test_throws ErrorException IOTO.lookup_by_support(tm, z2, [0.])
        # test normal
        @test IOTO.lookup_by_support(tm, x0, [0., 0., 0.]) == b
        @test IOTO.lookup_by_support(tm, y, [0., 0., 1.]) == a
    end
    # test internal_semi_infinite_variable
    @testset "internal_semi_infinite_variable" begin
        rv = SemiInfiniteVariableRef(m, SemiInfiniteVariableIndex(-1))
        # test errors
        @test_throws ErrorException InfiniteOpt.internal_semi_infinite_variable(rv, Val(:TransData))
        # test normal
        var = build_variable(error, x, supps)
        push!(transcription_data(tm).semi_infinite_vars, var)
        @test InfiniteOpt.internal_semi_infinite_variable(rv, Val(:TransData)) == var
        eval_supports(rv) == supps
    end
end

# Test the measure mapping methods
@testset "Measure Queries" begin
    # initialize tbe needed info
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, x, Infinite(par, pars))
    @variable(m, x0, Point(x, 0, [0, 0]))
    @variable(m, y)
    meas1 = support_sum(2par -2, par)
    meas2 = support_sum(x^2 - y, pars)
    tm = optimizer_model(m)
    data = transcription_data(tm)
    data.has_internal_supports = true
    s1 = Set([UserDefined])
    s2 = Set([InternalLabel])
    sd = Set{DataType}()
    data.supports = ([0., 1., NaN], [[0., 0.], [1., 1.], [NaN, NaN]])
    data.support_labels = ([s1, s2, sd], [s1, s1, sd])
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
        data.measure_mappings[meas1] = [-2 * zero(AffExpr)]
        data.measure_support_labels[meas1] = [s1]
        data.measure_mappings[meas2] = [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
        data.measure_support_labels[meas2] = [s1, s2]
        @test IOTO.transcription_variable(tm, meas1) == -2 * zero(AffExpr)
        expected = [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
        @test IOTO.transcription_variable(tm, meas2) == [expected[1]]
        @test IOTO.transcription_variable(tm, meas2, label = All) == expected
        # test ndarray 
        @test IOTO.transcription_variable(tm, meas2, ndarray = true) == [expected[1]]
        @test IOTO.transcription_variable(tm, meas2, label = All, ndarray = true) == expected
    end
    # test lookup_by_support
    @testset "lookup_by_support" begin
        # test errors
        @test_throws ErrorException IOTO.lookup_by_support(tm, meas1, Float64[])
        @test_throws ErrorException IOTO.lookup_by_support(tm, meas2, [0.])
        # test normal
        data.measure_lookup[meas1] = Dict(Float64[] => 1)
        data.measure_lookup[meas2] = Dict{Vector{Float64}, Int}([0] => 1, [1] => 2)
        @test IOTO.lookup_by_support(tm, meas1, Float64[]) == -2 * zero(AffExpr)
        @test IOTO.lookup_by_support(tm, meas2, [1.]) == b^2 + d^2 - 2a
    end
    # test InfiniteOpt.variable_supports
    @testset "InfiniteOpt.variable_supports" begin
        # test not in mappings error
        dvref = dispatch_variable_ref(meas1)
        delete!(data.measure_mappings, meas1)
        @test_throws ErrorException InfiniteOpt.variable_supports(tm, dvref)
        data.measure_mappings[meas1] = [-2 * zero(AffExpr)]
        # test supports are empty
        @test InfiniteOpt.variable_supports(tm, dvref) == ()
        dvref2 = dispatch_variable_ref(meas2)
        @test InfiniteOpt.variable_supports(tm, dvref2, label = All) == [(0.,), (1.,)]
        # test normal
        @test InfiniteOpt.variable_supports(tm, dvref) == ()
        @test InfiniteOpt.variable_supports(tm, dvref2) == [(0.,)]
        # test ndarray 
        @test InfiniteOpt.variable_supports(tm, dvref, ndarray = true) == [()]
        @test InfiniteOpt.variable_supports(tm, dvref2, ndarray = true) == [(0.,)]
        @test InfiniteOpt.variable_supports(tm, dvref2, label = All, ndarray = true) == [(0.,), (1.,)]
    end
end

# Test the support iteration methods 
@testset "Support Iteration Methods" begin 
    # initialize the needed info 
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1], num_supports = 2)
    @infinite_parameter(m, par in [0, 1], supports = [0, 1], 
                        derivative_method = OrthogonalCollocation(3))
    @variable(m, y, Infinite(par))
    d1 = @deriv(y, par)
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
        @test isequal(IOTO._collected_supports(dispatch_variable_ref(par)), expected)
        # dependent parameters
        expected = sort!([Float64[0., 0.], Float64[1., 1.], Float64[NaN, NaN]])
        @test isequal(sort!(IOTO._collected_supports(dispatch_variable_ref(pars[1]))), expected)
    end
    # test _collected_support_labels
    @testset "_collected_support_labels" begin 
        # independent parameters
        supps = Float64[0., 1., NaN]
        expected = [Set([UserDefined]), Set([UserDefined]), Set{DataType}()]
        @test isequal(IOTO._collected_support_labels(dispatch_variable_ref(par), supps), expected)
        # dependent parameters
        supps = [Float64[0., 0.], Float64[1., 1.], Float64[NaN, NaN]]
        expected = [Set([UniformGrid]), Set([UniformGrid]), Set{DataType}()]
        @test isequal(IOTO._collected_support_labels(dispatch_variable_ref(pars[1]), supps), expected)
    end
    # test set_parameter_supports
    @testset "set_parameter_supports" begin 
        add_supports(par, 0.6, label = InternalLabel)
        @test IOTO.set_parameter_supports(tm, m) isa Nothing 
        expected = ([[0., 0.], [1., 1.], [NaN, NaN]], [0., 0.3, 0.6, 0.8, 1., NaN])
        @test isequal(sort.(transcription_data(tm).supports), expected)
        @test IOTO.has_internal_supports(tm)
        expected = ([Set([UniformGrid]), Set([UniformGrid]), Set{DataType}()], 
                    [Set([UserDefined]), Set([InternalGaussLobatto]), 
                     Set([InternalLabel]), Set([InternalGaussLobatto]), 
                     Set([UserDefined]), Set{DataType}()])
        @test isequal(transcription_data(tm).support_labels, expected)
        @test supports(par, label = All) == [0., 0.3, 0.6, 0.8, 1.]
        @test has_generative_supports(par)
    end
    # test parameter_supports
    @testset "parameter_supports" begin 
    expected = ([[0., 0.], [1., 1.], [NaN, NaN]], [0., 0.3, 0.6, 0.8, 1., NaN])
        @test isequal(sort.(IOTO.parameter_supports(tm)), expected)
    end
    # test support_index_iterator with 1 argument
    @testset "support_index_iterator (1 Arg)" begin 
        @test IOTO.support_index_iterator(tm) == CartesianIndices((1:2, 1:5))
    end
    # test support_index_iterator with 2 argument2
    @testset "support_index_iterator (2 Args)" begin 
        @test IOTO.support_index_iterator(tm, [1, 2]) == CartesianIndices((1:2, 1:5))
        @test IOTO.support_index_iterator(tm, [1]) == CartesianIndices((1:2, 6:6))
    end
    # test index_to_support
    @testset "index_to_support" begin 
        @test IOTO.index_to_support(tm, first(CartesianIndices((1:2, 1:5)))) isa Vector
        @test isnan(IOTO.index_to_support(tm, last(IOTO.support_index_iterator(tm, [1])))[3])
    end
    # test index_to_labels
    @testset "index_to_labels" begin 
        idxs = CartesianIndices((1:2, 1:5))
        @test IOTO.index_to_labels(tm, first(idxs)) == Set([UserDefined, UniformGrid])
        @test IOTO.index_to_labels(tm, idxs[3]) == Set([UniformGrid, InternalGaussLobatto])
        @test IOTO.index_to_labels(tm, last(IOTO.support_index_iterator(tm, Int[]))) == Set{DataType}()
    end
end

# Test the expression mappings
@testset "Expression Queries" begin
    # initialize tbe needed info
    m = InfiniteModel()
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @infinite_parameter(m, par in [0, 1], supports = [0])
    @variable(m, x, Infinite(par, pars))
    @finite_parameter(m, finpar == 42)
    @variable(m, x0, Point(x, 0, [0, 0]))
    @variable(m, y)
    f = parameter_function((a,b) -> 1, (par, pars))
    add_supports(par, 1, label = InternalLabel)
    meas1 = support_sum(2par -2, par)
    meas2 = support_sum(x^2 - y, pars)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    # transcribe the variables and measures
    data = transcription_data(tm)
    data.has_internal_supports = true
    s1 = Set([UserDefined])
    s2 = Set([InternalLabel])
    data.finvar_mappings[y] = a
    data.finvar_mappings[x0] = b
    data.infvar_mappings[x] = [b, c, d]
    data.infvar_support_labels[x] = [s1, s2, s1]
    data.measure_mappings[meas1] = [-2 * zero(AffExpr)]
    data.measure_support_labels[meas1] = [s1]
    data.measure_mappings[meas2] = [a^2 + c^2 - 2a, b^2 + d^2 - 2a]
    data.measure_support_labels[meas2] = [s1, s2]
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
        @test IOTO.transcription_expression(tm, f, [0., 1., 0.]) == 1
    end
    # test transcription expression for semi_infinite variables with 3 args
    @testset "transcription_expression (Semi-Infinite Variable)" begin
        # semi_infinite of parameter function 
        rv = add_variable(m, build_variable(error, f, Dict(1=>1.)), 
                          add_support = false)
        @test IOTO.transcription_expression(tm, rv, [0., 1., 0.]) == 1
        # semi_infinite of infinite variable
        rv = add_variable(m, build_variable(error, x, Dict(1=>1.)), 
                          add_support = false)
        data.infvar_mappings[rv] = [b, c]
        lookups = Dict{Vector{Float64}, Int}([0, 0] => 1, [1, 0] => 2)
        data.infvar_lookup[rv] = lookups
        @test IOTO.transcription_expression(tm, rv, [1., 0., 0.]) == c
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
        # test normal
        expr = meas2 - 3y^2 - x0 - 2.3
        expected = b^2 + d^2 - 2a - 3a^2 - b - 2.3
        @test IOTO.transcription_expression(tm, expr, [1., 1., 1.]) == expected
        # test becomes a nonlinear expression
        expr = meas2 * x0
        @test IOTO.transcription_expression(tm, expr, [1., 1., 1.]) isa NonlinearExpression
        @test nl_expr_string(tm, MIME("text/plain"), tm.nlp_data.nlexpr[end]) == "+((+(-2.0 * a) + b * b + d * d) * b)"
    end
    # test transcription expression for NLPExprs with 3 args
    @testset "transcription_expression (NLPExpr)" begin
        expr = sin(y)
        @test IOTO.transcription_expression(tm, expr, [1., 1., 1.]) isa NonlinearExpression
        @test nl_expr_string(tm, MIME("text/plain"), tm.nlp_data.nlexpr[end]) == "sin(a)"
    end
    # test transcription expression for numbers with 3 args
    @testset "transcription_expression (Real)" begin
        expected = zero(AffExpr) + 42
        @test IOTO.transcription_expression(tm, 42, [1., 1., 1.]) == expected
    end
    # test transcription expression for Exprs with 2 args
    @testset "transcription_expression (Expr 2 Args)" begin
        # test infinite expr
        expr = meas2 - 3y^2 - x0 - 2.3
        expected = [-2a^2 + c^2 - 2a - b - 2.3, -3a^2 + b^2 + d^2 - 2a - b - 2.3]
        @test IOTO.transcription_expression(tm, expr, label = All) == expected
        @test IOTO.transcription_expression(tm, expr, label = All, ndarray = true) == expected
        expected = [-2a^2 + c^2 - 2a - b - 2.3]
        @test IOTO.transcription_expression(tm, expr) == expected
        @test IOTO.transcription_expression(tm, expr, ndarray = true) == expected
        # test finite expr 
        expr = 2x0 -y 
        expected = 2b- a 
        @test IOTO.transcription_expression(tm, expr) == expected
        @test IOTO.transcription_expression(tm, expr, ndarray = true) == [expected]
        # test NLPExpr 
        @test IOTO.transcription_expression(tm, sin(x0)) isa NonlinearExpression
        @test nl_expr_string(tm, MIME("text/plain"), tm.nlp_data.nlexpr[end]) == "sin(b)"
    end
    # test transcription expression for variables with 2 args
    @testset "transcription_expression (Variable 2 Args)" begin
        @test IOTO.transcription_expression(tm, x0) == b 
        @test IOTO.transcription_expression(tm, x) == [b, d] 
        @test IOTO.transcription_expression(tm, x, label = All) == [b, c, d] 
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
        @test supports(x0 - y, ndarray = true) == [()]
        @test supports(meas2 + y) == [(0.,)]
        @test supports(meas2 + y, ndarray = true) == [(0.,)]
        @test supports(meas2 + y, label = All) == [(0.,), (1.,)]
        @test supports(meas2 + y, label = All, ndarray = true) == [(0.,), (1.,)]
    end
end

# Test constraint mapping queries
@testset "Constraint Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0])
    add_supports(par, 1, label = InternalLabel)
    @variable(m, x, Infinite(par))
    @variable(m, x0, Point(x, 0))
    @variable(m, y)
    @constraint(m, c1, x + y - 2 <= 0)
    @constraint(m, c2, support_sum(x, par) == 0)
    @constraint(m, c3, x0 + y == 5)
    tm = optimizer_model(m)
    data = transcription_data(tm)
    data.has_internal_supports = true
    s1 = Set([UserDefined])
    s2 = Set([InternalLabel])
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @constraint(tm, tc1, b + a <= 2)
    @constraint(tm, tc2, c + a <= 2)
    @constraint(tm, tc3, 0.5b + 0.5c == 0)
    @constraint(tm, tc4, b + a == 5)
    @test IOTO.set_parameter_supports(tm, m) isa Nothing
    # test transcription_constraint (2 Args)
    @testset "transcription_constraint (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_constraint(tm, c1)
        # test normal
        data.constr_mappings[c1] = [tc1, tc2]
        data.constr_support_labels[c1] = [s1, s2]
        @test transcription_constraint(tm, c1, label = All) == [tc1, tc2]
        @test transcription_constraint(tm, c1) == [tc1]
        @test transcription_constraint(tm, c1, label = All, ndarray = true) == [tc1, tc2]
        @test transcription_constraint(tm, c1, ndarray = true) == [tc1]
        # test error
        @test_throws ErrorException transcription_constraint(tm, c2)
        # test normal
        data.constr_mappings[c2] = [tc3]
        data.constr_support_labels[c2] = [s1]
        @test transcription_constraint(tm, c2) == tc3
        @test transcription_constraint(tm, c2, ndarray = true) == [tc3]
        # test error
        @test_throws ErrorException transcription_constraint(tm, c3)
        # test normal
        data.constr_mappings[c3] = [tc4]
        data.constr_support_labels[c3] = [s1]
        @test transcription_constraint(tm, c3) == tc4
    end
    # test transcription_constraint (Single argument)
    @testset "transcription_constraint (1 Arg)" begin
        @test transcription_constraint(c1) == [tc1]
        @test transcription_constraint(c1, label = All) == [tc1, tc2]
        @test transcription_constraint(c2) == tc3
        @test transcription_constraint(c3, label = All) == tc4
    end
    # test optimizer_model_constraint extension
    @testset "optimizer_model_constraint" begin
        @test optimizer_model_constraint(c1, Val(:TransData), label = All) == [tc1, tc2]
        @test optimizer_model_constraint(c2, Val(:TransData)) == tc3
        @test optimizer_model_constraint(c3, Val(:TransData)) == tc4
    end
    # test constraint_supports 
    @testset "constraint_supports" begin
        # test error
        @test_throws ErrorException InfiniteOpt.constraint_supports(tm, c1)
        # test normal
        data.constr_supports[c1] = [(0.,), (1.,)]
        @test InfiniteOpt.constraint_supports(tm, c1) == [(0.,)]
        @test InfiniteOpt.constraint_supports(tm, c1, label = All) == [(0.,), (1.,)]
        @test InfiniteOpt.constraint_supports(tm, c1, ndarray = true) == [(0.,)]
        @test InfiniteOpt.constraint_supports(tm, c1, label = All, ndarray = true) == [(0.,), (1.,)]
        # test finite 
        data.constr_supports[c3] = [()]
        @test InfiniteOpt.constraint_supports(tm, c3) == ()
        @test InfiniteOpt.constraint_supports(tm, c3, ndarray = true) == [()]
    end
    # test supports
    @testset "supports" begin
        @test supports(c1, label = All) == [(0.,), (1.,)]
        @test supports(c1, label = All, ndarray = true) == [(0.,), (1.,)]
        @test supports(c1, label = InternalLabel) == [(1.,)]
        @test supports(c3) == ()
    end
end
