# test the add_measure_variable extensions
@testset "Expansion Extensions" begin
    # setup the model testing
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_parameter(m, 0 <= pars[1:2] <= 1, supports = [0, 1])
    @infinite_variable(m, x(par))
    @infinite_variable(m, y(par, pars))
    @point_variable(m, x(0), x0)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    data = transcription_data(tm)
    data.infinite_to_vars[x] = [a, b, c]
    data.infvar_to_supports[x] = [(0.,), (0.5,), (1.,)]
    data.point_to_var[x0] = a
    key = Val(:TransData)
    # test add_measure_variable for point variables
    @testset "add_measure_variable (Point)" begin
        var = PointVariable(InfiniteOpt._variable_info(x), x, VectorTuple{Float64}(1))
        idx = data.next_var_index - 1
        vref = PointVariableRef(m, idx)
        set_name(vref, "vref") # for test debugging
        @test add_measure_variable(tm, var, key) == vref
        @test data.point_to_var[vref] == c
    end
    # test add_measure_variable for reduced variables
    @testset "add_measure_variable (Reduced)" begin
        var = ReducedInfiniteInfo(y, Dict{Int, Float64}(1 => 0))
        idx = data.next_var_index - 1
        vref = ReducedInfiniteVariableRef(m, idx)
        @test add_measure_variable(tm, var, key) == vref
        @test data.reduced_info[idx] == var
    end
    # test reduction_info extension
    @testset "reduction_info" begin
        vref = ReducedInfiniteVariableRef(m, data.next_var_index)
        @test reduction_info(vref, key).infinite_variable_ref == y
        @test reduction_info(vref, key).eval_supports == Dict{Int, Float64}(1 => 0)
    end
    # test delete_reduced_variable extension
    @testset "delete_reduced_variable" begin
        vref = ReducedInfiniteVariableRef(m, data.next_var_index)
        @test delete_reduced_variable(tm, vref, key) isa Nothing
        @test_throws KeyError reduction_info(vref, key)
    end
end
