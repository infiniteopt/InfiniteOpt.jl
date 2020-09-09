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
    @variable(tm, d)
    data = transcription_data(tm)
    data.infvar_mappings[x] = [a, b, c]
    data.infvar_supports[x] = [(0.,), (0.5,), (1.,)]
    data.infvar_lookup[x] = Dict{Vector{Float64}, Int}([0] => 1, [0.5] => 2, [1] => 3)
    data.infvar_mappings[y] = [a, b, c, d]
    data.infvar_supports[y] = [(0., [0., 0.]), (0., [1., 1.]), (1., [0., 0.]), (1., [1., 1.])]
    data.infvar_lookup[y] = Dict{Vector{Float64}, Int}([0, 0, 0] => 1, [0, 1, 1] => 2, 
                                                       [1, 0, 0] => 3, [1, 1, 1] => 4)
    data.finvar_mappings[x0] = a
    key = Val(:TransData)
    IOTO.set_parameter_supports(tm, m)
    # test add_measure_variable for point variables
    @testset "add_measure_variable (Point)" begin
        info = JuMP.VariableInfo(true, 0., true, 0., true, 0., true, 0., true, true)
        var = PointVariable(info, x, [1.])
        vref = GeneralVariableRef(m, -1, PointVariableIndex)
        @test InfiniteOpt.add_measure_variable(tm, var, key) == vref
        @test transcription_variable(vref) == c
    end
    # test add_measure_variable for reduced variables
    @testset "add_measure_variable (Reduced)" begin
        var = ReducedVariable(y, Dict{Int, Float64}(1 => 0), [2, 3], [2])
        vref = GeneralVariableRef(m, -1, ReducedVariableIndex)
        @test InfiniteOpt.add_measure_variable(tm, var, key) == vref
        @test data.reduced_vars == [var]
        @test a in transcription_variable(vref)
        @test b in transcription_variable(vref)
        @test sort!(supports(vref)) == [([0., 0.], ), ([1., 1.], )]
    end
    # test delete_reduced_variable extension
    @testset "delete_reduced_variable" begin
        vref = ReducedVariableRef(m, ReducedVariableIndex(-1))
        @test InfiniteOpt.delete_reduced_variable(tm, vref, key) isa Nothing
    end
end
