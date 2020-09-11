# test adding new infinite set
@testset "Infinite Sets" begin
    # load in the extension
    include("./extensions/infinite_set.jl")

    # initialize model
    m = InfiniteModel()

    # test definition
    @test MyNewSet(0, 1) isa MyNewSet
    set = MyNewSet(0, 1)
    @test @infinite_parameter(m, set = set) isa GeneralVariableRef
    @test infinite_set(first(all_parameters(m))) == set
    @test @infinite_parameter(m, par in set, supports = [0, 1]) isa GeneralVariableRef
    infinite_set(m[:par]) == set
    supports(m[:par]) == [0., 1.]
    @test @infinite_parameter(m, par2 in set, num_supports = 3) isa GeneralVariableRef
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
    @test @infinite_variable(m, x(par) >= 0) isa GeneralVariableRef
    x = m[:x]
    
    # test measures
    @test @integral(x^2 + par, par, num_supports = 2, eval_method = GaussLegendre) isa GeneralVariableRef

    # test constraints
    @test @constraint(m, x + par <= 0) isa InfOptConstraintRef

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
    data1 = DiscreteMeasureData(t, [2.5, 2.5], [2., 4.], lower_bound = 0., upper_bound = 5.)
    data2 = DiscreteMeasureData(t, [1., 1.], [3., 6.], lower_bound = 1., upper_bound = 10.)

    # test definition
    @test NewMeasureData("test", data1) isa NewMeasureData
    new_data1 = NewMeasureData("test", data1)
    new_data2 = NewMeasureData("test", data2)

    # test data queries
    @test parameter_refs(new_data1) == t
    @test parameter_refs(new_data2) == t
    @test supports(new_data1) == supports(data1)
    @test supports(new_data2) == supports(data2)
    @test coefficients(new_data1) == coefficients(data1)
    @test coefficients(new_data2) == coefficients(data2)
    @test weight_function(new_data1) == weight_function(data1)
    @test weight_function(new_data2) == weight_function(data2)
    @test measure_data_in_hold_bounds(data1, ParameterBounds())
    @test measure_data_in_hold_bounds(data1, parameter_bounds(z))
    @test !measure_data_in_hold_bounds(data2, parameter_bounds(z))

    # test measure definition
    @test measure(x + z, new_data1) isa GeneralVariableRef
    @test_throws ErrorException measure(x + z, new_data2)
    @test new_measure(x^2, t, 0, 4, num_supports = 4) isa GeneralVariableRef
    @test_throws ErrorException new_measure(x^2 + z, t, 6, 10)

    # test expansion
    pvrefs = [GeneralVariableRef(m, 1, PointVariableIndex), GeneralVariableRef(m, 2, PointVariableIndex)]
    @test expand(measure(x + z, new_data1)) == 2.5 * (pvrefs[1] + pvrefs[2]) + 5z

    # test transcription
    @test @constraint(m, z == measure(x, new_data1)) isa InfOptConstraintRef
    @test build_optimizer_model!(m) isa Nothing
    @test num_variables(optimizer_model(m)) == 6
    println(optimizer_model(m))
    println(supports(t))

    # test deletion
    @test_throws ErrorException delete(m, t)
    indices = collect(keys(m.measures))
    for index in indices
        @test delete(m, MeasureRef(m, index)) isa Nothing
    end
    @test delete(m, t) isa Nothing
end

#=
# Test extensions of measure evaluation methods
@testset "Measure Evaluation" begin
    # load in the extension
    include("./extensions/measure_eval.jl")

    # set up the model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 5])
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, xi in Normal(0., 1.))
    @infinite_variable(m, y(t) >= 0)
    @infinite_variable(m, f(x))
    mref = integral(y^2 + t, t, 0, 4, num_supports = 5, eval_method = NewEvalMethod)
    @test supports(measure_data(mref)) == Array([0., 1., 2., 3., 4.])
    warn = "The method is implemented for independent multivariate parameters."
    @test_logs (:warn, warn) integral(f, x, num_supports = 3,
                                      eval_method = NewEvalMethod,
                                      independent = false)
    mref2 = integral(f, x, num_supports = 3, eval_method = NewEvalMethod,
                     independent = is_independent(x[1]))
    @test supports(measure_data(mref2)) == Float64[0 0.5 1; 0 0.5 1]
    @test_throws ErrorException integral(xi^2, eval_method = NewEvalMethod)
end

# Test otpimizer model extensions
@testset "Optimizer Model" begin
    # load in the extension
    include("./extensions/optimizer_model.jl")

    # setup the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par))
    @point_variable(m, x(0), x0)
    @hold_variable(m, y)
    data = DiscreteMeasureData(par, [0.5, 0.5], [0, 1])
    @constraint(m, c1, x + y - 2 <= 0)
    @constraint(m, c2, measure(x, data) == 0)
    @constraint(m, c3, x0 + y == 5)
    @objective(m, Min, y)

    # test optimizer model constructor
    @test NewReformModel() isa Model
    @test NewReformModel().ext[:ReformData] isa NewReformData

    # test data extraction
    @test_throws ErrorException reform_data(optimizer_model(m))
    @test reform_data(NewReformModel()) isa NewReformData

    # test replacing the optimizer model
    @test set_optimizer_model(m, NewReformModel()) isa Nothing
    @test haskey(optimizer_model(m).ext, :ReformData)

    # test retrival errors
    @test_throws ErrorException optimizer_model_variable(x)
    @test_throws ErrorException optimizer_model_variable(x0)
    @test_throws ErrorException optimizer_model_variable(y)
    @test_throws ErrorException optimizer_model_constraint(c1)
    @test_throws ErrorException optimizer_model_constraint(c2)
    @test_throws ErrorException optimizer_model_constraint(c3)
    @test_throws ErrorException supports(x)
    @test_throws ErrorException supports(c1)
    @test_throws ErrorException supports(c2)
    @test_throws ErrorException supports(c3)
    @test_throws ErrorException parameter_refs(c1)
    @test_throws ErrorException parameter_refs(c2)
    @test_throws ErrorException parameter_refs(c3)

    # test build_optimizer_model!
    @test build_optimizer_model!(m, my_kwarg = true) isa Nothing
    @test haskey(optimizer_model(m).ext, :ReformData)
    @test num_variables(optimizer_model(m)) == 3

    # test retrivals
    @test optimizer_model_variable(x, my_kwarg = true) isa Vector{VariableRef}
    @test optimizer_model_variable(x0, my_kwarg = true) isa VariableRef
    @test optimizer_model_variable(y, my_kwarg = true) isa VariableRef
    @test optimizer_model_constraint(c1, my_kwarg = true) isa Vector
    @test optimizer_model_constraint(c2, my_kwarg = true) isa Vector
    @test optimizer_model_constraint(c3, my_kwarg = true) isa ConstraintRef
    @test supports(x, my_kwarg = true) == [(0.,), (1.,)]
    @test supports(c1, my_kwarg = true) == [(0.,), (1.,)]
    @test parameter_refs(x) == (par,)
    @test parameter_refs(c1, my_kwarg = true) == (par,)

    # test more retrival errors
    @test_throws ErrorException supports(c2)
    @test_throws ErrorException parameter_refs(c2)
    @test_throws ErrorException supports(c3)
    @test_throws ErrorException parameter_refs(c3)

    # test optimization with rebuild
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    @test set_optimizer_model_ready(m, false) isa Nothing
    @test set_optimizer(m, mockoptimizer) isa Nothing
    @test set_silent(m)
    @test set_time_limit_sec(m, 42) == 42
    @test optimize!(m) isa Nothing
    @test get_optimizer_attribute(m, MOI.Silent())

    # prepare for result queries
    mockoptimizer = JuMP.backend(m).optimizer.model
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ResultCount(), 2)
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.SolveTime(), 0.42)
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(optimizer_model_variable(x)[1]), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(optimizer_model_variable(x)[2]), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(optimizer_model_variable(y)), 1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(optimizer_model_constraint(c1)[1]), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(optimizer_model_constraint(c1)[2]), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(optimizer_model_constraint(c2)[1]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(optimizer_model_constraint(c3)), -1.0)

    # test result queries
    @test termination_status(m) == MOI.OPTIMAL
    @test raw_status(m) == "solver specific string"
    @test result_count(m) == 2
    @test primal_status(m) == MOI.FEASIBLE_POINT
    @test dual_status(m) == MOI.FEASIBLE_POINT
    @test solve_time(m) == 0.42
    @test has_values(m)
    @test has_duals(m)
    @test objective_bound(m) == 2
    @test objective_value(m) == -1
    @test value(x) == [-1, 0]
    @test value(x0) == -1
    @test value(y) == 1
    @test dual(c1) == [-1, -1]
    @test dual(c2) == [0]
    @test dual(c3) == -1
    @test optimizer_index(x) == optimizer_index.(optimizer_model_variable(x))
    @test optimizer_index(x0) == optimizer_index(optimizer_model_variable(x0))
    @test optimizer_index(y) == optimizer_index(optimizer_model_variable(y))
    @test optimizer_index(c1) == optimizer_index.(optimizer_model_constraint(c1))
    @test optimizer_index(c2) == optimizer_index.(optimizer_model_constraint(c2))
    @test optimizer_index(c3) == optimizer_index(optimizer_model_constraint(c3))
    @test shadow_price(c1) == [-1, -1]
    @test shadow_price(c2) == [0]
    @test shadow_price(c3) == -1
end
=#