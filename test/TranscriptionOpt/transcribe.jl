# Test variable initializers
@testset "Variable Initializers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1])
    @infinite_variable(m, x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    data = DiscreteMeasureData(par, [0.5, 0.5], [0, 1])
    @constraint(m, c1, x + z - 2 <= 0)
    @constraint(m, c2, measure(x + y, data) - w == 0)
    @constraint(m, c3, x0 + y0 == 5)
    tm = optimizer_model(m)
    # test _initialize_global_variables
    @testset "_initialize_global_variables" begin
        @global_variable(m)
        @test isa(IOTO._initialize_global_variables(tm, m), Nothing)
        @test length(transcription_data(tm).global_to_var) == 2
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
    # test _list_supports
    @testset "_list_supports" begin
        @test IOTO._list_supports((par, pars)) == [[0, 1], supports(pars)]
    end
    # test _make_supports
    @testset "_make_supports" begin
        expected = sort([(0, 0), (0, 1), (1, 0), (1, 1)])
        @test sort(IOTO._make_supports((par, par))) == expected
        @test IOTO._make_supports((par, pars)) isa Vector
    end
    # test _initialize_infinite_variables
    @testset "_initialize_infinite_variables" begin
        @infinite_variable(m, parameter_refs = (par,))
        @test isa(IOTO._initialize_infinite_variables(tm, m), Nothing)
        @test length(transcription_data(tm).infinite_to_vars) == 2
        @test transcription_variable(tm, x) isa Vector{VariableRef}
        @test transcription_variable(tm, y) isa Vector{VariableRef}
        @test name(transcription_variable(tm, x)[1]) == "x(support: 1)"
        @test name(transcription_variable(tm, y)[3]) == "y(support: 3)"
        @test has_lower_bound(transcription_variable(tm, x)[1])
        @test is_binary(transcription_variable(tm, y)[2])
        @test is_fixed(transcription_variable(tm, y)[4])
        @test is_integer(transcription_variable(tm, x)[2])
        @test start_value(transcription_variable(tm, y)[1]) == 0.
        @test supports(x) == [(0,), (1,)]
        @test length(supports(y)) == 4
    end
    # test _update_point_mapping
    @testset "_update_point_mapping" begin
        # test first point variable
        @test isa(IOTO._update_point_mapping(tm, x0, x, parameter_values(x0)),
                  Nothing)
        @test transcription_data(tm).point_to_var[x0] == transcription_variable(tm, x)[1]
        # test second point variable
        @test isa(IOTO._update_point_mapping(tm, y0, y, parameter_values(y0)),
                  Nothing)
        @test transcription_data(tm).point_to_var[y0] == transcription_variable(tm, y)[1]
        # test error
        @test_throws ErrorException IOTO._update_point_mapping(tm, x0, x, (0.5,))
    end
    # test _update_point_info
    @testset "_update_point_info" begin
        vref = transcription_variable(tm, x0)
        # test lower bound update
        fix(vref, 0, force = true)
        set_lower_bound(x0, 0)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test lower_bound(vref) == 0
        # test upper bound update
        fix(vref, 0, force = true)
        set_upper_bound(x0, 0)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test upper_bound(vref) == 0
        # test fix update
        fix(x0, 0, force = true)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test fix_value(vref) == 0
        # test binary update
        unset_integer(x0)
        set_binary(x0)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test is_binary(vref)
        # test integer update
        unset_binary(x0)
        set_integer(x0)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test is_integer(vref)
        # test start value
        set_start_value(x0, 0.5)
        @test isa(IOTO._update_point_info(tm, x0), Nothing)
        @test start_value(vref) == 0.5
        # Undo changes
        unfix(vref)
        unfix(x0)
        set_lower_bound(vref, 0)
        set_lower_bound(x0, 0)
        set_start_value(vref, NaN)
        set_start_value(x0, NaN)
        delete!(transcription_data(tm).point_to_var, x0)
        delete!(transcription_data(tm).point_to_var, y0)
    end
    # test _map_point_variables
    @testset "_map_point_variables" begin
        @test isa(IOTO._map_point_variables(tm, m), Nothing)
        @test length(transcription_data(tm).infinite_to_vars) == 2
        @test transcription_variable(tm, x0) isa VariableRef
        @test transcription_variable(tm, y0) isa VariableRef
        @test name(transcription_variable(tm, x0)) == "x(support: 1)"
        @test name(transcription_variable(tm, y0)) == "y(support: 1)"
        @test has_lower_bound(transcription_variable(tm, x0))
        @test is_integer(transcription_variable(tm, x0))
        @test has_lower_bound(transcription_variable(tm, y0))
        @test has_upper_bound(transcription_variable(tm, y0))
        @test is_integer(transcription_variable(tm, y0))
        @test start_value(transcription_variable(tm, y0)) == 0.
    end
end

# Test parameter search methods
@testset "Parameter Search Methods" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= x <= 1)
    @infinite_parameter(m, 0 <= y <= 1)
    @infinite_parameter(m, 0 <= z[1:3] <= 1, container = SparseAxisArray)
    supp1 = (1, convert(JuMPC.SparseAxisArray, [2, 3, 4]))
    tup1 = (y, z)
    supp2 = (1, 2, convert(JuMPC.SparseAxisArray, [3, 4, 5]))
    tup2 = (x, y, z)
    supp3 = (1, 2, JuMPC.SparseAxisArray(Dict((1,) => 3, (3,) => 5)))
    tup3 = (x, y, JuMPC.SparseAxisArray(Dict((1,) => z[1], (3,) => z[3])))
    # test _parameter_tuple_index
    @testset "_parameter_tuple_index" begin
        @test IOTO._parameter_tuple_index(x, tup1) == -1
        @test IOTO._parameter_tuple_index(x, tup2) == 1
        @test IOTO._parameter_tuple_index(y, tup2) == 2
        @test IOTO._parameter_tuple_index(z[2], tup2) == 3
        @test IOTO._parameter_tuple_index(z[2], tup3) === 3
        @test IOTO._parameter_tuple_index(z[3], tup3) == 3
    end
    # test _parameter_value
    @testset "_parameter_value" begin
        @test IOTO._parameter_value(x, supp1, tup1) === NaN
        @test IOTO._parameter_value(x, supp2, tup2) == 1
        @test IOTO._parameter_value(y, supp2, tup2) == 2
        @test IOTO._parameter_value(z[2], supp2, tup2) == 4
        @test IOTO._parameter_value(z[2], supp3, tup3) === NaN
        @test IOTO._parameter_value(z[3], supp3, tup3) == 5
    end
end

# Test "_map_to_variable"
@testset "_map_to_variable" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    index = m.next_var_index + 1
    m.reduced_info[index] = ReducedInfiniteInfo(y, Dict(1 => 1))
    yr = ReducedInfiniteVariableRef(m, index)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    @variable(tm, e)
    @variable(tm, f)
    @variable(tm, g)
    @variable(tm, h)
    supp1 = convert(JuMPC.SparseAxisArray, [0, 0])
    supp2 = convert(JuMPC.SparseAxisArray, [1, 1])
    # test for FiniteVariableRefs
    @testset "FiniteVariableRef" begin
        # test cannot find
        @test_throws ErrorException IOTO._map_to_variable(x0, (), (), tm)
        @test_throws ErrorException IOTO._map_to_variable(y0, (), (), tm)
        @test_throws ErrorException IOTO._map_to_variable(z, (), (), tm)
        @test_throws ErrorException IOTO._map_to_variable(w, (), (), tm)
        # prepare transcrition model with mapped variables
        tm.ext[:TransData].point_to_var[x0] = a
        tm.ext[:TransData].point_to_var[y0] = b
        tm.ext[:TransData].global_to_var[z] = c
        tm.ext[:TransData].global_to_var[w] = d
        # test the function
        @test IOTO._map_to_variable(x0, (), (), tm) == a
        @test IOTO._map_to_variable(y0, (), (), tm) == b
        @test IOTO._map_to_variable(z, (), (), tm) == c
        @test IOTO._map_to_variable(w, (), (), tm) == d
    end
    # test for InfiniteVariableRefs
    @testset "InfiniteVariableRef" begin
        # prepare transcrition model with mapped variables
        tm.ext[:TransData].infinite_to_vars[x] = [a, e]
        tm.ext[:TransData].infinite_to_vars[y] = [b, f, g, h]
        tm.ext[:TransData].infvar_to_supports[x] = [(0,), (1,)]
        tm.ext[:TransData].infvar_to_supports[y] = [(0, supp1), (0, supp2),
                                                    (1, supp1), (1, supp2)]
        # test cannot find
        @test_throws ErrorException IOTO._map_to_variable(x, (0.5, supp1),
                                                          (par, pars), tm)
        @test_throws ErrorException IOTO._map_to_variable(y, (0.5, supp1),
                                                          (par, pars), tm)
        # test the function
        @test IOTO._map_to_variable(x, (-0., supp1), (par, pars), tm) == a
        @test IOTO._map_to_variable(x, (1., supp1), (par, pars), tm) == e
        @test IOTO._map_to_variable(y, (-0., supp1), (par, pars), tm) == b
        @test IOTO._map_to_variable(y, (0, supp2), (par, pars), tm) == f
        @test IOTO._map_to_variable(y, (1., supp1), (par, pars), tm) == g
        @test IOTO._map_to_variable(y, (1, supp2), (par, pars), tm) == h
    end
    # test for ReducedInfiniteVariableRefs
    @testset "ReducedInfiniteVariableRef" begin
        # test cannot find
        supp3 = convert(JuMPC.SparseAxisArray, [0.5, 0.5])
        @test_throws ErrorException IOTO._map_to_variable(yr, (0.5, supp3),
                                                          (par, pars), tm)
        # test the function
        @test IOTO._map_to_variable(yr, (0., supp1), (par, pars), tm) == g
        @test IOTO._map_to_variable(yr, (supp1,), (pars,), tm) == g
        @test IOTO._map_to_variable(yr, (0., supp2), (par, pars), tm) == h
        @test IOTO._map_to_variable(yr, (supp2,), (pars,), tm) == h
    end
    # test for ParameterRefs
    @testset "ParameterRef" begin
        # test cannot find
        @test_throws ErrorException IOTO._map_to_variable(par, (supp1,),
                                                          (pars,), tm)
        # test the function
        @test IOTO._map_to_variable(par, (0., supp1), (par, pars), tm) == 0
        @test IOTO._map_to_variable(pars[2], (0., supp1), (par, pars), tm) == 0
    end
end

# Test _support(s)_in_bounds
@testset "Support Bound Checking" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    supp1 = convert(JuMPC.SparseAxisArray, [0, 0])
    supp2 = convert(JuMPC.SparseAxisArray, [1, 1])
    # test _support_in_bounds
    @testset "_support_in_bounds" begin
        # test in bounds
        @test IOTO._support_in_bounds((0, supp1), (par, pars), Dict())
        # test out of bounds
        bounds = Dict(par => IntervalSet(0.5, 1))
        @test !IOTO._support_in_bounds((0, supp1), (par, pars), bounds)
        bounds = Dict(par => IntervalSet(0, 0.5))
        @test !IOTO._support_in_bounds((.7, supp1), (par, pars), bounds)
        # test parameter bounds not relavent
        @test IOTO._support_in_bounds((supp1,), (pars,), bounds)
    end
    @testset "_supports_in_bounds" begin
        supps = [(0, supp1), (0, supp2), (1, supp1), (1, supp2)]
        bounds = Dict(par => IntervalSet(0.5, 1))
        @test IOTO._supports_in_bounds(supps, (par, pars), bounds) == [false,
                                                              false, true, true]
    end
end

# Test _truncate_supports
@testset "_make_transcription_function" begin
    # initialize the models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    m.next_var_index += 1
    m.reduced_info[m.next_var_index] = ReducedInfiniteInfo(y, Dict(1 => 1))
    yr = ReducedInfiniteVariableRef(m, m.next_var_index)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    meas3 = measure(y^2 + z, data1)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    @variable(tm, e)
    @variable(tm, f)
    @variable(tm, g)
    @variable(tm, h)
    supp1 = convert(JuMPC.SparseAxisArray, [0, 0])
    supp2 = convert(JuMPC.SparseAxisArray, [1, 1])
    # transcribe the variables
    tm.ext[:TransData].point_to_var[x0] = a
    tm.ext[:TransData].point_to_var[y0] = b
    tm.ext[:TransData].global_to_var[z] = c
    tm.ext[:TransData].global_to_var[w] = d
    tm.ext[:TransData].infinite_to_vars[x] = [a, e]
    tm.ext[:TransData].infinite_to_vars[y] = [b, f, g, h]
    tm.ext[:TransData].infvar_to_supports[x] = [(0,), (1,)]
    tm.ext[:TransData].infvar_to_supports[y] = [(0, supp1), (0, supp2),
                                                (1, supp1), (1, supp2)]
    # test FiniteVariableRefs
    @testset "FiniteVariableRef" begin
        # test point variables
        @test IOTO._make_transcription_function(x0, tm) == (a, ())
        @test IOTO._make_transcription_function(y0, tm) == (b, ())
        # test global variables
        @test IOTO._make_transcription_function(z, tm) == (c, ())
        @test IOTO._make_transcription_function(w, tm) == (d, ())
    end
    # test InfiniteVariableRefs
    @testset "InfiniteVariableRef" begin
        # test without bounds
        @test IOTO._make_transcription_function(x, tm) == ([a, e], (par,),
                                                           [(0,), (1,)])
        @test IOTO._make_transcription_function(y, tm) == ([b, f, g, h],
                                                           (par, pars),
                                                           [(0, supp1), (0, supp2),
                                                            (1, supp1), (1, supp2)])
        # test with bounds
        bounds = Dict(par => IntervalSet(0.5, 1))
        @test IOTO._make_transcription_function(x, tm, bounds) == ([e], (par,),
                                                                   [(1,)])
        @test IOTO._make_transcription_function(y, tm, bounds) == ([g, h],
                                                                   (par, pars),
                                                                   [(1, supp1),
                                                                   (1, supp2)])
    end
    # test ParameterRefs
    @testset "ParameterRef" begin
        # test without bounds
        @test IOTO._make_transcription_function(par, tm) == ([0, 1], (par,), [0, 1])
        @test IOTO._make_transcription_function(pars[1], tm) == ([0, 1], (pars[1],),
                                                                 [0, 1])
        # test with bounds
        bounds = Dict(par => IntervalSet(0.5, 1))
        @test IOTO._make_transcription_function(par, tm, bounds) == ([1], (par,),
                                                                     [1])
        @test IOTO._make_transcription_function(pars[1], tm, bounds) == ([0, 1],
                                                             (pars[1],), [0, 1])
    end
    # test Finite GenericAffExprs
    @testset "Finite AffExpr" begin
        # prepare expressions
        aff1 = x0 + 2z - 3
        aff2 = 3z - w
        # test expressions
        @test IOTO._make_transcription_function(aff1, tm)[1] == a + 2c - 3
        @test IOTO._make_transcription_function(aff2, tm)[1] == 3c - d
    end
    # test Finite GenericQuadExprs
    @testset "Finite QuadExpr" begin
        # prepare expressions
        quad1 = x0^2 - y0 * w + 2z - 3
        quad2 = 2z * w + 3z - w
        # test expressions
        @test IOTO._make_transcription_function(quad1, tm)[1] == a^2 - b * d + 2c - 3
        @test IOTO._make_transcription_function(quad2, tm)[1] == 2c * d + 3c - d
    end
    # test General GenericAffExprs
    @testset "General AffExpr" begin
        # test only having 1 variable that is finite
        @test IOTO._make_transcription_function(z - 2, tm)[1] == c - 2
        @test IOTO._make_transcription_function(2w - w, tm)[1] == d + 0
        # test having only 1 variable that is infinite
        @test IOTO._make_transcription_function(x - 2, tm) == ([a - 2, e - 2],
                                                           (par,), [(0,), (1,)])
        @test IOTO._make_transcription_function(par - 2, tm) == ([-2., -1.],
                                                                 (par,), [0, 1])
        # test without bounds
        expr = meas1 + 2x - par + z
        expected = ([3a + e - 2d + c, a + 3e - 2d + c - 1], (par,), [(0,), (1,)])
        @test IOTO._make_transcription_function(expr, tm) == expected
        # test with bounds
        bounds = Dict(par => IntervalSet(0.5, 1))
        expr = meas2 + 2x - pars[2] + z
        expected = ([b + g + 2e + c, f + h + 2e + c - 1], (pars, par),
                    [(supp1, 1), (supp2, 1)])
        @test IOTO._make_transcription_function(expr, tm, bounds) == expected
    end
    # test General GenericQuadExprs
    @testset "General QuadExpr" begin
        # test only having 1 variable that is finite
        expr = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        expr.aff = z - 2
        @test IOTO._make_transcription_function(expr, tm)[1] == c - 2
        expr.aff = 2w - w
        @test IOTO._make_transcription_function(expr, tm)[1] == d + 0
        # test having only 1 variable that is infinite
        expr.aff = x - 2
        @test IOTO._make_transcription_function(expr, tm) == ([a - 2, e - 2],
                                                           (par,), [(0,), (1,)])
        expr.aff = par - 2
        @test IOTO._make_transcription_function(expr, tm) == ([-2., -1.], (par,),
                                                              [0, 1])
        # test without bounds
        expr = x0^2 + 2 * z * y0 - par * x + w * par + 2 * par^2 - 2par + meas1 - 3
        expected = ([a^2 + 2 * c * b + a + e - 2d - 3,
                     a^2 + 2 * c * b + a + 0e - d - 3], (par,), [(0,), (1,)])
        @test IOTO._make_transcription_function(expr, tm) == expected
        # test without bounds
        bounds = Dict(par => IntervalSet(0.5, 1))
        expr = x0^2 + 2 * z * y0 - pars[1] * x + w * par + 2 * par^2 - 2par + meas2 - 3
        expected = ([a^2 + 2 * c * b + b + g + d - 3,
                     a^2 + 2 * c * b + f + h - e + d - 3], (par, pars),
                     [(1, supp1), (1, supp2)])
        @test IOTO._make_transcription_function(expr, tm, bounds) == expected
        # test AffExpr dispatch to QuadExpr
        bounds = Dict(par => IntervalSet(0.5, 1))
        expr = 2meas3 - x + 2
        expected = ([2 * (b^2 + g^2 + 2c) - e + 2, 2 * (f^2 + h^2 + 2c) - e + 2],
                    (par, pars), [(1, supp1), (1, supp2)])
        @test IOTO._make_transcription_function(expr, tm, bounds) == expected
    end
    # test MeasureFiniteVariableRef Expressions
    @testset "MeasureFiniteVariableRef Exprs" begin
        # test finite AffExpr
        expr = meas1 + z - 2w
        expected = (a + e - 4d + c, ())
        @test IOTO._make_transcription_function(expr, tm) == expected
        # test infinite AffExpr
        expr = meas1 + x - 2w
        expected = ([2a + e - 4d, a + 2e - 4d], (par,), [(0,), (1,)])
        @test IOTO._make_transcription_function(expr, tm) == expected
        # test finite QuadExpr
        expr = meas1 + z^2 - 2w
        expected = (c^2 + a + e - 4d, ())
        @test IOTO._make_transcription_function(expr, tm) == expected
        # test infinite QuadExpr
        expr = meas1 + x^2 - 2w
        expected = ([a^2 + a + e - 4d, e^2 + a + e - 4d], (par,), [(0,), (1,)])
        @test IOTO._make_transcription_function(expr, tm) == expected
    end
    # test MeasrureRefs
    @testset "MeasureRef" begin
        @test IOTO._make_transcription_function(meas1, tm)[1] == a + e - 2d
        @test IOTO._make_transcription_function(meas2, tm)[1] == [b + g, f + h]
    end
    # test numeric AffExpr
    @testset "Numeric AffExpr" begin
        expr = zero(AffExpr) + 42
        @test IOTO._make_transcription_function(expr, tm) == (expr, ())
    end
    # test the fallback
    @testset "Fallback" begin
        expr = zero(QuadExpr)
        @test_throws ErrorException IOTO._make_transcription_function(expr, tm)
    end
end

# Test _set_objective
@testset "_set_objective" begin
    # initialize the models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    @variable(tm, e)
    @variable(tm, f)
    @variable(tm, g)
    @variable(tm, h)
    supp1 = convert(JuMPC.SparseAxisArray, [0, 0])
    supp2 = convert(JuMPC.SparseAxisArray, [1, 1])
    # transcribe the variables
    tm.ext[:TransData].global_to_var[z] = c
    tm.ext[:TransData].global_to_var[w] = d
    tm.ext[:TransData].infinite_to_vars[x] = [a, e]
    tm.ext[:TransData].infinite_to_vars[y] = [b, f, g, h]
    tm.ext[:TransData].infvar_to_supports[x] = [(0,), (1,)]
    tm.ext[:TransData].infvar_to_supports[y] = [(0, supp1), (0, supp2),
                                                (1, supp1), (1, supp2)]
    # test empty objective
    @test isa(IOTO._set_objective(tm, m), Nothing)
    @test objective_function(tm) == zero(AffExpr)
    # test finite objective
    @objective(m, Min, z - w)
    @test isa(IOTO._set_objective(tm, m), Nothing)
    @test objective_function(tm) == c - d
    @test objective_sense(tm) == MOI.MIN_SENSE
    # test measure objective
    @objective(m, Max, z - 2meas1)
    @test isa(IOTO._set_objective(tm, m), Nothing)
    @test objective_function(tm) == c - 2a - 2e + 4d
    @test objective_sense(tm) == MOI.MAX_SENSE
    # test error with infinite objective
    @objective(m, Max, z - 2meas2)
    @test_throws ErrorException IOTO._set_objective(tm, m)
end

# Test the helper functions for constraint transcription
@testset "Constraint Helpers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par) >= 0, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - z, data1)
    tm = optimizer_model(m)
    @variable(tm, xt[1:2])
    @variable(tm, zt)
    # Setup up the constraints
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(tm, c1t[i = 1:2], xt[i] + supports(par)[i] - zt == 0)
    @constraint(tm, c2t, zt >= -3)
    @constraint(tm, c3t, xt[1] + xt[2] - 3zt == 0)
    # Test _set_mapping (InfiniteConstraintRefs)
    @testset "_set_mapping (Infinite)" begin
        @test isa(IOTO._set_mapping(c1, c1t), Nothing)
        @test tm.ext[:TransData].infinite_to_constrs[c1] == [c1t[1], c1t[2]]
    end
    # Test _set_mapping (MeasureConstraintRefs that are infinite)
    @testset "_set_mapping (Infinite Measures)" begin
        @test isa(IOTO._set_mapping(c3, [c3t]), Nothing)
        @test tm.ext[:TransData].measure_to_constrs[c3] == [c3t]
    end
    # Test _set_mapping (MeasureConstraintRefs that are finite)
    @testset "_set_mapping (Finite Measures)" begin
        @test isa(IOTO._set_mapping(c3, c3t), Nothing)
        @test tm.ext[:TransData].measure_to_constrs[c3] == [c3t]
    end
    # Test _set_mapping (FiniteConstraintRefs)
    @testset "_set_mapping (Finite)" begin
        @test isa(IOTO._set_mapping(c2, c2t), Nothing)
        @test tm.ext[:TransData].finite_to_constr[c2] == c2t
    end
    # Test _set_supports (InfiniteConstraintRefs)
    @testset "_set_supports (Infinite)" begin
        @test isa(IOTO._set_supports(tm, c1, [(0,), (1,)]), Nothing)
        @test tm.ext[:TransData].infconstr_to_supports[c1] == [(0,), (1,)]
    end
    # Test _set_supports (MeasureConstraintRefs)
    @testset "_set_supports (Measure)" begin
        @test isa(IOTO._set_supports(tm, c3, [(0,), (1,)]), Nothing)
        @test tm.ext[:TransData].measconstr_to_supports[c3] == [(0,), (1,)]
    end
    # Test _set_parameter_refs (InfiniteConstraintRefs)
    @testset "_set_parameter_refs (Infinite)" begin
        @test isa(IOTO._set_parameter_refs(tm, c1, (par,)), Nothing)
        @test tm.ext[:TransData].infconstr_to_params[c1] == (par,)
    end
    # Test _set_parameter_refs (MeasureConstraintRefs)
    @testset "_set_parameter_refs (Measure)" begin
        @test isa(IOTO._set_parameter_refs(tm, c3, (par,)), Nothing)
        @test tm.ext[:TransData].measconstr_to_params[c3] == (par,)
    end
    # Test _root_name (constraints)
    @testset "_root_name (Constraints)" begin
        # test simple name return
        @test InfiniteOpt._root_name(c1) == "c1"
        # test anonymous name
        @test InfiniteOpt._root_name(@constraint(m, z == 1)) == "noname"
        # test array of constraints
        @test InfiniteOpt._root_name(@constraint(m, a[1:2], z == 1)[1]) == "a"
    end
end

# Test _set_constraints
@testset "_set_constraints" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @variable(tm, d)
    @variable(tm, e)
    @variable(tm, f)
    @variable(tm, g)
    @variable(tm, h)
    supp1 = convert(JuMPC.SparseAxisArray, [0, 0])
    supp2 = convert(JuMPC.SparseAxisArray, [1, 1])
    # transcribe the variables
    tm.ext[:TransData].point_to_var[x0] = a
    tm.ext[:TransData].point_to_var[y0] = b
    tm.ext[:TransData].global_to_var[z] = c
    tm.ext[:TransData].global_to_var[w] = d
    tm.ext[:TransData].infinite_to_vars[x] = [a, e]
    tm.ext[:TransData].infinite_to_vars[y] = [b, f, g, h]
    tm.ext[:TransData].infvar_to_supports[x] = [(0,), (1,)]
    tm.ext[:TransData].infvar_to_supports[y] = [(0, supp1), (0, supp2),
                                                (1, supp1), (1, supp2)]
    # Setup up the constraints
    bounds = Dict(par => IntervalSet(0.5, 1))
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, parameter_bounds = bounds)
    @constraint(m, c5, meas2 == 0)
    # test the main function
    @test isa(IOTO._set_constraints(tm, m), Nothing)
    @test length(list_of_constraint_types(tm)) == 3
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c1][1]).func == a - c
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c1][2]).func == e - c
    @test constraint_object(tm.ext[:TransData].finite_to_constr[c2]).func == c + a
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c4][1]).func == -b + g + e
    @test name(tm.ext[:TransData].finite_to_constr[c2]) == "c2"
    @test name(tm.ext[:TransData].infinite_to_constrs[c1][1]) == "c1(Support: 1)"
    @test tm.ext[:TransData].infconstr_to_params[c1] == (par, )
    @test tm.ext[:TransData].infconstr_to_params[c4] == (pars, par)
    @test tm.ext[:TransData].infconstr_to_supports[c1] == [(0,), (1,)]
    @test tm.ext[:TransData].measconstr_to_params[c5] == (pars,)
end

# Test variable info mappers
@testset "Variable Info Mappers" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, 1 >= x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    tm = optimizer_model(m)
    IOTO._initialize_global_variables(tm, m)
    IOTO._initialize_infinite_variables(tm, m)
    IOTO._map_point_variables(tm, m)
    xt = transcription_variable(tm, x)
    yt = transcription_variable(tm, y)
    zt = transcription_variable(tm, z)
    wt = transcription_variable(tm, w)
    # test _map_info_constraints (FiniteVariableRef)
    @testset "_map_info_constraints (Finite)" begin
        # test z
        @test isa(IOTO._map_info_constraints(z, tm), Nothing)
        @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(z)] == LowerBoundRef(zt)
        @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(z)] == UpperBoundRef(zt)
        @test tm.ext[:TransData].finite_to_constr[BinaryRef(z)] == BinaryRef(zt)
        # test w
        @test isa(IOTO._map_info_constraints(w, tm), Nothing)
        @test tm.ext[:TransData].finite_to_constr[FixRef(w)] == FixRef(wt)
        @test tm.ext[:TransData].finite_to_constr[IntegerRef(w)] == IntegerRef(wt)
    end
    # test _map_info_constraints (InfiniteVariableRef)
    @testset "_map_info_constraints (Infinite)" begin
        # test x
        @test isa(IOTO._map_info_constraints(x, tm), Nothing)
        @test tm.ext[:TransData].infinite_to_constrs[LowerBoundRef(x)] == LowerBoundRef.(xt)
        @test tm.ext[:TransData].infconstr_to_params[LowerBoundRef(x)] == parameter_refs(x)
        @test tm.ext[:TransData].infconstr_to_supports[LowerBoundRef(x)] == [(0,), (1,)]
        @test tm.ext[:TransData].infinite_to_constrs[UpperBoundRef(x)] == UpperBoundRef.(xt)
        @test tm.ext[:TransData].infconstr_to_params[UpperBoundRef(x)] == parameter_refs(x)
        @test tm.ext[:TransData].infconstr_to_supports[UpperBoundRef(x)] == [(0,), (1,)]
        @test tm.ext[:TransData].infinite_to_constrs[IntegerRef(x)] == IntegerRef.(xt)
        @test tm.ext[:TransData].infconstr_to_params[IntegerRef(x)] == parameter_refs(x)
        @test tm.ext[:TransData].infconstr_to_supports[IntegerRef(x)] == [(0,), (1,)]
        # test y
        @test isa(IOTO._map_info_constraints(y, tm), Nothing)
        @test tm.ext[:TransData].infinite_to_constrs[FixRef(y)] == FixRef.(yt)[2:end]
        @test tm.ext[:TransData].infconstr_to_params[FixRef(y)] == parameter_refs(y)
        @test tm.ext[:TransData].infconstr_to_supports[FixRef(y)] == supports(y)[2:end]
        @test tm.ext[:TransData].infinite_to_constrs[BinaryRef(y)] == [BinaryRef(yt[i]) for i = 2:4]
        @test tm.ext[:TransData].infconstr_to_params[BinaryRef(y)] == parameter_refs(y)
        @test tm.ext[:TransData].infconstr_to_supports[BinaryRef(y)] == supports(y)[2:end]
    end
end

# Test _map_variable_info_constraints
@testset "_map_variable_info_constraints" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, 1 >= x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    tm = optimizer_model(m)
    IOTO._initialize_global_variables(tm, m)
    IOTO._initialize_infinite_variables(tm, m)
    IOTO._map_point_variables(tm, m)
    xt = transcription_variable(tm, x)
    yt = transcription_variable(tm, y)
    zt = transcription_variable(tm, z)
    wt = transcription_variable(tm, w)
    # test running the function
    @test isa(IOTO._map_variable_info_constraints(tm, m), Nothing)
    # test z
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(z)] == LowerBoundRef(zt)
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(z)] == UpperBoundRef(zt)
    @test tm.ext[:TransData].finite_to_constr[BinaryRef(z)] == BinaryRef(zt)
    # test w
    @test tm.ext[:TransData].finite_to_constr[FixRef(w)] == FixRef(wt)
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(w)] == IntegerRef(wt)
    # test x0
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(x0)] == LowerBoundRef(xt[1])
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(x0)] == UpperBoundRef(xt[1])
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(x0)] == IntegerRef(xt[1])
    # test y0
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(y0)] == LowerBoundRef(yt[1])
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(y0)] == UpperBoundRef(yt[1])
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(y0)] == IntegerRef(yt[1])
    # test x
    @test tm.ext[:TransData].infinite_to_constrs[LowerBoundRef(x)] == LowerBoundRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[LowerBoundRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[LowerBoundRef(x)] == [(0,), (1,)]
    @test tm.ext[:TransData].infinite_to_constrs[UpperBoundRef(x)] == UpperBoundRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[UpperBoundRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[UpperBoundRef(x)] == [(0,), (1,)]
    @test tm.ext[:TransData].infinite_to_constrs[IntegerRef(x)] == IntegerRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[IntegerRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[IntegerRef(x)] == [(0,), (1,)]
    # test y
    @test tm.ext[:TransData].infinite_to_constrs[FixRef(y)] == FixRef.(yt)[2:end]
    @test tm.ext[:TransData].infconstr_to_params[FixRef(y)] == parameter_refs(y)
    @test tm.ext[:TransData].infconstr_to_supports[FixRef(y)] == supports(y)[2:end]
    @test tm.ext[:TransData].infinite_to_constrs[BinaryRef(y)] == [BinaryRef(yt[i]) for i = 2:4]
    @test tm.ext[:TransData].infconstr_to_params[BinaryRef(y)] == parameter_refs(y)
    @test tm.ext[:TransData].infconstr_to_supports[BinaryRef(y)] == supports(y)[2:end]
end

# Test TranscriptionModel constructor
@testset "TranscriptionModel (Full Build)" begin
    # initialize model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1],
                        container = SparseAxisArray)
    @infinite_variable(m, 1 >= x(par) >= 0, Int)
    @infinite_variable(m, y(par, pars) == 2, Bin, start = 0)
    @point_variable(m, x(0), x0)
    @point_variable(m, y(0, [0, 0]), 0 <= y0 <= 1, Int)
    @global_variable(m, 0 <= z <= 1, Bin)
    @global_variable(m, w == 1, Int, start = 1)
    data1 = DiscreteMeasureData(par, [1, 1], [0, 1])
    meas1 = measure(x - w, data1)
    meas2 = measure(y, data1)
    bounds = Dict(par => IntervalSet(0.5, 1))
    @constraint(m, c1, x + par - z == 0)
    @constraint(m, c2, z + x0 >= -3)
    @constraint(m, c3, meas1 + z == 0)
    @constraint(m, c4, meas2 - 2y0 + x <= 1, parameter_bounds = bounds)
    @constraint(m, c5, meas2 == 0)
    @objective(m, Min, x0 + meas1)
    # test basic usage
    @test isa(set_optimizer_model(m, TranscriptionModel(m)), Nothing)
    tm = optimizer_model(m)
    # test global variables
    @test length(transcription_data(tm).global_to_var) == 2
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
    # test infinite variables
    @test length(transcription_data(tm).infinite_to_vars) == 2
    @test transcription_variable(tm, x) isa Vector{VariableRef}
    @test transcription_variable(tm, y) isa Vector{VariableRef}
    @test name(transcription_variable(tm, x)[1]) == "x(support: 1)"
    @test name(transcription_variable(tm, y)[3]) == "y(support: 3)"
    @test has_lower_bound(transcription_variable(tm, x)[1])
    @test is_binary(transcription_variable(tm, y)[2])
    @test is_fixed(transcription_variable(tm, y)[4])
    @test is_integer(transcription_variable(tm, x)[2])
    @test start_value(transcription_variable(tm, y)[1]) == 0.
    @test supports(x) == [(0,), (1,)]
    @test length(supports(y)) == 4
    # test point variables
    @test length(transcription_data(tm).infinite_to_vars) == 2
    @test transcription_variable(tm, x0) isa VariableRef
    @test transcription_variable(tm, y0) isa VariableRef
    @test name(transcription_variable(tm, x0)) == "x(support: 1)"
    @test name(transcription_variable(tm, y0)) == "y(support: 1)"
    @test has_lower_bound(transcription_variable(tm, x0))
    @test is_integer(transcription_variable(tm, x0))
    @test has_lower_bound(transcription_variable(tm, y0))
    @test has_upper_bound(transcription_variable(tm, y0))
    @test is_integer(transcription_variable(tm, y0))
    @test start_value(transcription_variable(tm, y0)) == 0.
    # test objective
    xt = transcription_variable(tm, x)
    wt = transcription_variable(tm, w)
    @test objective_function(tm) == 2xt[1] + xt[2] - 2wt
    @test objective_sense(tm) == MOI.MIN_SENSE
    # test constraints
    zt = transcription_variable(tm, z)
    yt = transcription_variable(tm, y)
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c1][1]).func == xt[1] - zt
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c1][2]).func == xt[2] - zt
    @test constraint_object(tm.ext[:TransData].finite_to_constr[c2]).func == zt + xt[1]
    @test constraint_object(tm.ext[:TransData].infinite_to_constrs[c4][1]).func == -yt[1] + yt[2] + xt[2]
    @test name(tm.ext[:TransData].finite_to_constr[c2]) == "c2"
    @test name(tm.ext[:TransData].infinite_to_constrs[c1][1]) == "c1(Support: 1)"
    @test tm.ext[:TransData].infconstr_to_params[c1] == (par, )
    @test tm.ext[:TransData].infconstr_to_params[c4] == (pars, par)
    @test tm.ext[:TransData].infconstr_to_supports[c1] == [(0,), (1,)]
    @test tm.ext[:TransData].measconstr_to_params[c5] == (pars,)
    # test info constraints
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(z)] == LowerBoundRef(zt)
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(z)] == UpperBoundRef(zt)
    @test tm.ext[:TransData].finite_to_constr[BinaryRef(z)] == BinaryRef(zt)
    @test tm.ext[:TransData].finite_to_constr[FixRef(w)] == FixRef(wt)
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(w)] == IntegerRef(wt)
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(x0)] == LowerBoundRef(xt[1])
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(x0)] == UpperBoundRef(xt[1])
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(x0)] == IntegerRef(xt[1])
    @test tm.ext[:TransData].finite_to_constr[LowerBoundRef(y0)] == LowerBoundRef(yt[1])
    @test tm.ext[:TransData].finite_to_constr[UpperBoundRef(y0)] == UpperBoundRef(yt[1])
    @test tm.ext[:TransData].finite_to_constr[IntegerRef(y0)] == IntegerRef(yt[1])
    @test tm.ext[:TransData].infinite_to_constrs[LowerBoundRef(x)] == LowerBoundRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[LowerBoundRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[LowerBoundRef(x)] == [(0,), (1,)]
    @test tm.ext[:TransData].infinite_to_constrs[UpperBoundRef(x)] == UpperBoundRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[UpperBoundRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[UpperBoundRef(x)] == [(0,), (1,)]
    @test tm.ext[:TransData].infinite_to_constrs[IntegerRef(x)] == IntegerRef.(xt)
    @test tm.ext[:TransData].infconstr_to_params[IntegerRef(x)] == parameter_refs(x)
    @test tm.ext[:TransData].infconstr_to_supports[IntegerRef(x)] == [(0,), (1,)]
    @test tm.ext[:TransData].infinite_to_constrs[FixRef(y)] == FixRef.(yt)[2:end]
    @test tm.ext[:TransData].infconstr_to_params[FixRef(y)] == parameter_refs(y)
    @test tm.ext[:TransData].infconstr_to_supports[FixRef(y)] == supports(y)[2:end]
    @test tm.ext[:TransData].infinite_to_constrs[BinaryRef(y)] == [BinaryRef(yt[i]) for i = 2:4]
    @test tm.ext[:TransData].infconstr_to_params[BinaryRef(y)] == parameter_refs(y)
    @test tm.ext[:TransData].infconstr_to_supports[BinaryRef(y)] == supports(y)[2:end]
    # test with solver
    mockoptimizer = with_optimizer(MOIU.MockOptimizer,
                                   MOIU.Model{Float64}(),
                                   eval_objective_value=false)
    @test isa(set_optimizer_model(m, TranscriptionModel(m, mockoptimizer)), Nothing)
end
