# test adding new infinite set
@testset "Infinite Sets" begin
    # load in the extension
    include("./extensions/infinite_set.jl")

    # initialize model
    m = InfiniteModel()

    # test definition
    @test MyNewSet(0, 1) isa MyNewSet
    set = MyNewSet(0, 1)
    @test @infinite_parameter(m, set = set) isa ParameterRef
    m.params[1] == set
    @test @infinite_parameter(m, par in set, supports = [0, 1]) isa ParameterRef
    infinite_set(m[:par]) == set
    supports(m[:par]) == [0., 1.]
    @test @infinite_parameter(m, par2 in set, num_supports = 3) isa ParameterRef
    supports(m[:par2]) == [0., 0.5, 1.]
    @test @infinite_parameter(m, [1:2], set = set, num_supports = 3) isa Vector

    # set support methods
    par = m[:par]
    @test_throws ErrorException add_supports(par, 2)

    # test bound methods
    @test has_lower_bound(m[:par2])
    @test lower_bound(m[:par2]) == 0
    @test set_lower_bound(m[:par2], -1) isa Nothing
    @test lower_bound(m[:par2]) == -1
    @test has_upper_bound(m[:par2])
    @test upper_bound(m[:par2]) == 1
    @test set_upper_bound(m[:par2], 0) isa Nothing
    @test upper_bound(m[:par2]) == 0

    # add variables
    @test @infinite_variable(m, x(par) >= 0) isa InfiniteVariableRef
    x = m[:x]

    # test measures
    register_eval_method(m, MyNewSet, default_methods[1])
    @test measure(x^2 + par, par, num_supports = 2, eval_method = gauss_legendre) isa MeasureRef

    # test constraints
    @test @constraint(m, x + par <= 0) isa GeneralConstraintRef

    # transcribe the model
    @test build_optimizer_model!(m) isa Nothing
    @test num_variables(optimizer_model(m)) == 4
end

# Test extensions of measure data
@testset "Measure Data" begin
    # load in the extension
    include("./extensions/measure_data.jl")

    # setup the model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10])
    @infinite_variable(m, x(t) >= 0)
    @hold_variable(m, z, parameter_bounds = (t in [0, 5]))
    data1 = generate_measure_data(t, 2, 0, 5, eval_method = gauss_legendre)
    data2 = generate_measure_data(t, 4, 0, 10)

    # test definition
    @test NewMeasureData("test", data1) isa NewMeasureData
    new_data1 = NewMeasureData("test", data1)
    new_data2 = NewMeasureData("test", data2)

    # test data queries
    @test parameter_refs(new_data1) == t
    @test parameter_refs(new_data2) == t
    @test measure_name(new_data1) == "test"
    @test measure_name(new_data2) == "test"
    @test supports(new_data1) == supports(data1)
    @test supports(new_data2) == supports(data2)
    @test measure_data_in_hold_bounds(data1, ParameterBounds())
    @test measure_data_in_hold_bounds(data1, parameter_bounds(z))
    @test !measure_data_in_hold_bounds(data2, parameter_bounds(z))

    # test measure definition
    @test measure(x + z, new_data1) isa MeasureRef
    @test_throws ErrorException measure(x + z, new_data2)
    @test new_measure(x^2, t, 0, 4, num_supports = 4) isa MeasureRef
    @test_throws ErrorException new_measure(x^2 + z, t, 6, 10)

    # test expansion
    index = m.next_var_index
    pvrefs = [PointVariableRef(m, index + 1), PointVariableRef(m, index + 2)]
    @test expand(measure(x + z, new_data1)) == 2.5 * (pvrefs[1] + pvrefs[2]) + 5z

    # test transcription
    delete_supports(t)
    @test @constraint(m, z == measure(x, new_data1)) isa GeneralConstraintRef
    @test build_optimizer_model!(m) isa Nothing
    @test num_variables(optimizer_model(m)) == 3

    # test deletion
    @test_throws ErrorException delete(m, t)
    indices = collect(keys(m.measures))
    for index in indices
        @test delete(m, MeasureRef(m, index)) isa Nothing
    end
    @test delete(m, t) isa Nothing
end
