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
    # test _make_support_indices
    @testset "_make_support_indices" begin
        expected = sort([(0, 0), (0, 1), (1, 0), (1, 1)])
        @test sort(IOTO._make_support_indices((par, par))) == expected
        @test IOTO._make_support_indices((par, pars)) isa Vector
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
