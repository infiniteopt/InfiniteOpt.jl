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
        expected = ([b + g + 2e + c, f + h + 2e + c - 1], (par, pars),
                    [(1, supp1), (1, supp2)])
        @test IOTO._make_transcription_function(expr, tm, bounds) == expected
    end
    # test General GenericQuadExprs
    @testset "General QuadExpr" begin
        # TODO include test for affexpr dispatch
    end
end
