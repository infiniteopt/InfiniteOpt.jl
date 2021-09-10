# test the add_measure_variable extensions
@testset "Expansion Extensions" begin
    # setup the model testing
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1], supports = [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1], supports = [0, 1])
    @variable(m, x, Infinite(par))
    @variable(m, y, Infinite(par, pars))
    @variable(m, x0, Point(x, 0))
    @variable(m, y0, SemiInfinite(y, 0, pars))
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
    data.infvar_mappings[y0] = [a, b]
    data.infvar_supports[y0] = [(0., [0., 0.]), (0., [1., 1.])]
    data.infvar_lookup[y0] = Dict{Vector{Float64}, Int}([0, 0, 0] => 1, [0, 1, 1] => 2)
    data.finvar_mappings[x0] = a
    key = Val(:TransData)
    IOTO.set_parameter_supports(tm, m)
    # test add_point_variable
    @testset "add_point_variable" begin
        # add one that was already added to the infinite model 
        @test isequal(InfiniteOpt.add_point_variable(tm, x, Float64[0], key), x0)
        @test transcription_variable(x0) == a
        # add one that hasn't been added
        vref = GeneralVariableRef(m, -1, PointVariableIndex)
        @test isequal(InfiniteOpt.add_point_variable(tm, x, Float64[1], key), vref)
        @test transcription_variable(vref) == c
        # add one that has been added internally
        @test isequal(InfiniteOpt.add_point_variable(tm, x, Float64[1], key), vref)
        @test transcription_variable(vref) == c
    end
    # test add_semi_infinite_variable
    @testset "add_semi_infinite_variable" begin
        # add one that was already added to the infinite model
        var = SemiInfiniteVariable(y, Dict{Int, Float64}(1 => 0), [2, 3], [2])
        @test isequal(InfiniteOpt.add_semi_infinite_variable(tm, var, key), y0)
        @test transcription_variable(y0) == [a, b]
        # add a new one
        var = SemiInfiniteVariable(y, Dict{Int, Float64}(1 => 1), [2, 3], [2])
        vref = GeneralVariableRef(m, -1, SemiInfiniteVariableIndex)
        @test isequal(InfiniteOpt.add_semi_infinite_variable(tm, var, key), vref)
        @test isequal(data.semi_infinite_vars, [var])
        @test c in transcription_variable(vref)
        @test d in transcription_variable(vref)
        @test sort!(supports(vref)) == [([0., 0.], ), ([1., 1.], )]
        # add one that has already been added internally
        @test isequal(InfiniteOpt.add_semi_infinite_variable(tm, var, key), vref)
        @test isequal(data.semi_infinite_vars, [var])
        @test c in transcription_variable(vref)
        @test d in transcription_variable(vref)
        @test sort!(supports(vref)) == [([0., 0.], ), ([1., 1.], )]
    end
end
