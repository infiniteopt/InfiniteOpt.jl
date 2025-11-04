# test adding new infinite domain
@testset "Infinite Domains" begin
    # load in the extension
    include("./extensions/infinite_domain.jl")

    # initialize model
    m = InfiniteModel()

    # test definition
    @test MyNewDomain(0, 1) isa MyNewDomain
    domain = MyNewDomain(0, 1)
    @test @infinite_parameter(m, domain = domain) isa GeneralVariableRef
    @test infinite_domain(first(all_parameters(m))) == domain
    @test @infinite_parameter(m, par in domain, supports = [0, 1]) isa GeneralVariableRef
    infinite_domain(m[:par]) == domain
    supports(m[:par]) == [0., 1.]
    @test @infinite_parameter(m, par2 in domain, num_supports = 3) isa GeneralVariableRef
    supports(m[:par2]) == [0., 0.5, 1.]
    @test @infinite_parameter(m, [1:2], domain = domain, num_supports = 3) isa Vector

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
    @test @variable(m, x >= 0, Infinite(par)) isa GeneralVariableRef
    x = m[:x]
    
    # test measures
    @test @integral(x^2 + par, par, num_nodes = 2, eval_method = GaussLegendre()) isa GeneralVariableRef
    @test expect(x, par) isa GeneralVariableRef

    # test constraints
    @test @constraint(m, x + par <= 0) isa InfOptConstraintRef

    # transcribe the model
    @test build_transformation_backend!(m) isa Nothing
    @test num_variables(m.backend.model) == 4
end

# test using the new derivative evaluation method
@testset "Derivative Method" begin
    # load in the extension
    include("./extensions/derivative_method.jl")

    # initialize model and objects
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10], num_supports = 3)
    @variable(m, x, Infinite(t))
    d = deriv(x, t)
    method = MyDerivMethod(0.5)

    # test definition and addition to a parameter
    @test method isa MyDerivMethod
    @test set_derivative_method(t, method) isa Nothing 
    @test derivative_method(t) isa MyDerivMethod

    # test preliminaries to evaluation 
    @test generative_support_info(method) == UniformGenerativeInfo([0.5], InternalLabel)
    @test support_label(method) == InternalLabel

    # test evaluation 
    @test InfiniteOpt.evaluate_derivative(d, x, method, m) isa Vector 
    @test has_generative_supports(t)
    @test num_supports(t, label = All) == 5
    @test evaluate(d) isa Nothing 
    @test length(derivative_constraints(d)) == 4
    @test num_constraints(m) == 4
end

# Test extensions of measure data
@testset "Measure Data" begin
    # load in the extension
    include("./extensions/measure_data.jl")

    # setup the model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10])
    @variable(m, x >= 0, Infinite(t))
    @variable(m, z)
    data1 = DiscreteMeasureData(t, [2.5, 2.5], [2., 4.], lower_bound = 0., upper_bound = 5.)
    data2 = DiscreteMeasureData(t, [1., 1.], [3., 6.], lower_bound = 1., upper_bound = 10.)

    # test definition
    @test NewMeasureData("test", data1) isa NewMeasureData
    new_data1 = NewMeasureData("test", data1)
    new_data2 = NewMeasureData("test", data2)

    # test data queries
    @test isequal(parameter_refs(new_data1), t)
    @test isequal(parameter_refs(new_data2), t)
    @test supports(new_data1) == supports(data1)
    @test supports(new_data2) == supports(data2)
    @test coefficients(new_data1) == coefficients(data1)
    @test coefficients(new_data2) == coefficients(data2)
    @test weight_function(new_data1) == weight_function(data1)
    @test weight_function(new_data2) == weight_function(data2)

    # test measure definition
    @test measure(x + z, new_data1) isa GeneralVariableRef
    @test new_measure(x^2, t, 0, 4, num_supports = 4) isa GeneralVariableRef
    @test_throws ErrorException new_measure(x^2 + z, t, 6, 10)

    # test expansion
    pvrefs = [GeneralVariableRef(m, i, PointVariableIndex) for i in 1:2] 
    @test isequal_canonical(expand(measure(x + z, new_data1)), 2.5 * (pvrefs[1] + pvrefs[2]) + 5z)

    # test transcription
    @test @constraint(m, z == measure(x, new_data1)) isa InfOptConstraintRef
    @test build_transformation_backend!(m) isa Nothing
    @test num_variables(m.backend.model) == 6

    # test deletion
    @test_throws ErrorException delete(m, t)
    indices = collect(keys(m.measures))
    for index in indices
        @test delete(m, MeasureRef(m, index)) isa Nothing
    end
    @test delete(m, x) isa Nothing
    @test delete(m, t) isa Nothing
end

# Test extensions of measure evaluation methods
@testset "Measure Evaluation" begin
    # load in the extension
    include("./extensions/measure_eval.jl")

    # set up the model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 5])
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    @infinite_parameter(m, p[1:2] in [0, 1])
    @infinite_parameter(m, xi ~ Normal(0., 1.))
    @variable(m, y >= 0, Infinite(t))
    @variable(m, f, Infinite(x...))
    mref = integral(y^2 + t, t, 0, 4, num_supports = 5, 
                    eval_method = NewUniEvalMethod())
    @test supports(measure_data(mref)) == Array([0., 1., 2., 3., 4.])
    warn = "The method is implemented for independent multivariate parameters."
    @test_logs (:warn, warn) integral(f, p, num_supports = 3,
                                      eval_method = NewMultiEvalMethod())
    mref2 = integral(f, x, num_supports = 3, eval_method = NewMultiEvalMethod())
    @test supports(measure_data(mref2)) == Float64[0 0.5 1; 0 0.5 1]
end

# Test extensions of generative support info
@testset "Generative Support Info" begin
    # load in the extension
    include("./extensions/generative_info.jl")

    # set up the model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 5], num_supports = 4)
    pref = dispatch_variable_ref(t)
    
    # test creation and basic functions 
    @test MyGenInfo(2) isa AbstractGenerativeInfo
    @test support_label(MyGenInfo(2)) == MyGenLabel 
    @test length(make_generative_supports(MyGenInfo(2), t, supports(t))) == 6

    # test incorporation into infinite parameter 
    @test InfiniteOpt._set_generative_support_info(pref, MyGenInfo(2)) isa Nothing 
    @test add_generative_supports(t) isa Nothing 
    @test num_supports(t, label = All) == 10
end

# Test Backend extensions
@testset "Backend" begin
    # load in the extension
    include("./extensions/backend.jl")

    # setup the infinite model
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, x, Infinite(par))
    @variable(m, x0, Point(x, 0))
    @variable(m, y)
    meas = integral(x, par)
    @parameter_function(m, pf == sin(par))
    @finite_parameter(m, p == 3)
    @constraint(m, c1, x + y + pf - 2 <= 0)
    @constraint(m, c2, meas == 0)
    @constraint(m, c3, x0 + y == 5 + p)
    @objective(m, Min, y)

    # test backend constructor
    @test NewReformBackend() isa NewReformBackend
    @test NewReformBackend().data isa NewReformData

    # test basic extraction
    @test transformation_model(NewReformBackend()) isa Model
    @test transformation_data(NewReformBackend()) isa NewReformData

    # test replacing the backend
    b = NewReformBackend()
    @test set_transformation_backend(m, b) isa Nothing
    @test transformation_backend(m) === b

    # test making InfiniteModel with the new optimizer model
    @test InfiniteModel(NewReformBackend()) isa InfiniteModel 
    @test transformation_backend(InfiniteModel(NewReformBackend())) isa NewReformBackend

    # test retrival errors
    @test_throws ErrorException transformation_variable(x)
    @test_throws ErrorException transformation_variable(x0)
    @test_throws ErrorException transformation_variable(y)
    @test_throws ErrorException transformation_constraint(c1)
    @test_throws ErrorException transformation_constraint(c2)
    @test_throws ErrorException transformation_constraint(c3)
    @test_throws ErrorException supports(x)
    @test_throws ErrorException supports(c1)
    @test_throws ErrorException supports(c2)
    @test_throws ErrorException supports(c3)

    # test build_transformation_backend!
    @test build_transformation_backend!(m, my_kwarg = true) isa Nothing
    @test transformation_backend_ready(m)
    @test num_variables(transformation_model(m)) == 16

    # test retrivals
    @test transformation_variable(x, my_kwarg = true) isa Vector{VariableRef}
    @test transformation_variable(x0, my_kwarg = true) isa VariableRef
    @test transformation_variable(y, my_kwarg = true) isa VariableRef
    @test transformation_variable(pf, my_kwarg = true) isa Vector{VariableRef}
    @test transformation_variable(p, my_kwarg = true) isa VariableRef
    @test transformation_variable(meas, my_kwarg = true) isa Vector{VariableRef}
    @test transformation_constraint(c1, my_kwarg = true) isa Vector{<:ConstraintRef}
    @test transformation_constraint(c2, my_kwarg = true) isa Vector{<:ConstraintRef}
    @test transformation_constraint(c3, my_kwarg = true) isa Vector{<:ConstraintRef}
    @test transformation_expression(x^2) == zero(AffExpr)
    @test supports(x, my_kwarg = true) == [(0.,), (1.,)]
    @test supports(pf) == [(0.,), (1.,)]
    @test supports(y) == ()
    @test supports(meas) == [(-1.,), (-2.,)]
    @test supports(c1, my_kwarg = true) == [(2.,), (3.,)]
    @test supports(x + y) == [(-42.,), (1.,)]
    @test variable_ref_type(m.backend) == VariableRef

    # test parameter updates
    @test set_parameter_value(pf, cos) isa Nothing
    @test set_parameter_value(p, 10) isa Nothing
    @test parameter_value(transformation_variable(p)) == 10

    # test optimization with rebuild
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    @test set_transformation_backend_ready(m, false) isa Nothing
    @test set_optimizer(m, mockoptimizer) isa Nothing
    @test set_silent(m) isa Nothing
    @test set_time_limit_sec(m, 42.) isa Nothing
    @test optimize!(m) isa Nothing
    @test get_attribute(m, MOI.Silent())

    # prepare for result queries
    mockoptimizer = JuMP.backend(m).optimizer.model
    for v in all_variables(m.backend.model)
        MOI.set(mockoptimizer, MOI.VariablePrimal(),
                JuMP.optimizer_index(v), 1.0)
    end
    MOI.set(mockoptimizer, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mockoptimizer, MOI.RawStatusString(), "solver specific string")
    MOI.set(mockoptimizer, MOI.ResultCount(), 2)
    MOI.set(mockoptimizer, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.DualStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mockoptimizer, MOI.SolveTimeSec(), 0.42)
    MOI.set(mockoptimizer, MOI.ObjectiveValue(), -1.0)
    MOI.set(mockoptimizer, MOI.ObjectiveBound(), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(x)[1]), -1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(x)[2]), 0.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(y)), 1.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(x0)), 42.)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(meas)[1]), 2.0)
    MOI.set(mockoptimizer, MOI.VariablePrimal(),
            JuMP.optimizer_index(transformation_variable(meas)[2]), -2.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c1)[1]), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c1)[2]), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c2)[1]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c2)[2]), -1.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c3)[1]), 0.0)
    MOI.set(mockoptimizer, MOI.ConstraintDual(),
            JuMP.optimizer_index(transformation_constraint(c3)[2]), -1.0)

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
    @test value(x0) == 42.
    @test value(y) == 1
    @test value(meas) == [2., -2.]
    @test value(x + y) == 0.
    @test dual(c1) == [-1, -1]
    @test dual(c2) == [0., -1.]
    @test dual(c3) == [0., -1.]
    @test is_solved_and_feasible(m)
    @test warmstart_backend_start_values(m) isa Nothing
    @test start_value(transformation_variable(y)) == 1
#     @test optimizer_index(x) == optimizer_index.(transformation_variable(x))
#     @test optimizer_index(x0) == optimizer_index(transformation_variable(x0))
#     @test optimizer_index(y) == optimizer_index(transformation_variable(y))
#     @test optimizer_index(c1) == optimizer_index.(transformation_constraint(c1))
#     @test optimizer_index(c2) == optimizer_index.(transformation_constraint(c2))
#     @test optimizer_index(c3) == optimizer_index.(transformation_constraint(c3))
#     @test shadow_price(c1) == [-1, -1]
#     @test shadow_price(c2) == [0., -1.]
#     @test shadow_price(c3) == [0., -1.]
end
