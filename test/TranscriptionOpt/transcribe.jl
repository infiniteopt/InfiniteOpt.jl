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

end

# Test _support_in_bounds
@testset "_support_in_bounds" begin

end

# Test _truncate_supports
@testset "_make_transcription_function" begin

end
