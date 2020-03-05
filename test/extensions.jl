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
    @test measure(x^2 + par, par, num_supports = 2, method = gauss_legendre) isa MeasureRef

    # test constraints
    @test @constraint(m, x + par <= 0) isa GeneralConstraintRef

    # transcribe the model
    @test build_optimizer_model!(m) isa Nothing
    @test num_variables(optimizer_model(m)) == 4
end
