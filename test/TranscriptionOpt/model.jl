# Test the TranscriptionData datatype
@testset "TranscriptionData" begin
    @test TranscriptionData isa DataType
    @test TranscriptionData().finite_to_constr isa Dict
end

# Test basic definition and queries
@testset "Basic Definition and Queries" begin
    # initialize needed data
    mockoptimizer = with_optimizer(MOIU.MockOptimizer,
                                   MOIU.Model{Float64}(),
                                   eval_objective_value=false)
    # test TranscriptionModel (no factory)
    @testset "TranscriptionModel (Default)" begin
        @test haskey(TranscriptionModel().ext, :TransData)
    end
    # test TranscriptionModel (with factory)
    @testset "TranscriptionModel (Optimizer)" begin
        @test haskey(TranscriptionModel(mockoptimizer).ext, :TransData)
    end
    # test is_transcription_model
    @testset "is_transcription_model" begin
        @test is_transcription_model(TranscriptionModel())
        @test !is_transcription_model(Model())
    end
    # test transcription_data
    @testset "transcription_data" begin
        @test transcription_data(TranscriptionModel()).point_to_var isa Dict
        @test_throws ErrorException transcription_data(Model())
    end
end

# Test variable mapping queries
@testset "Variable Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par))
    @point_variable(m, x(0), x0)
    @hold_variable(m, y)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    # test transcription_variable (HoldVariableRef)
    @testset "transcription_variable (Hold)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, y)
        # test normal
        tm.ext[:TransData].hold_to_var[y] = a
        @test transcription_variable(tm, y) == a
    end
    # test transcription_variable (InfiniteVariableRef)
    @testset "transcription_variable (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, x)
        # test normal
        tm.ext[:TransData].infinite_to_vars[x] = [b, c]
        @test transcription_variable(tm, x) == [b, c]
    end
    # test transcription_variable (PointVariableRef)
    @testset "transcription_variable (Point)" begin
        # test error
        @test_throws ErrorException transcription_variable(tm, x0)
        # test normal
        tm.ext[:TransData].point_to_var[x0] = b
        @test transcription_variable(tm, x0) == b
    end
    # test supports for infinite variable with 2 inputs
    @testset "supports (Model, Infinite Variable)" begin
        # test error
        @test_throws ErrorException supports(tm, x)
        # test normal
        tm.ext[:TransData].infvar_to_supports[x] = [(0.,), (1.,)]
        @test supports(tm, x) == [(0.,), (1.,)]
        # undo changes
        delete!(tm.ext[:TransData].infvar_to_supports, x)
    end
    # test supports for infinite variable
    @testset "supports (Infinite Variable)" begin
        # test error
        @test_throws ErrorException supports(x)
        # test normal
        tm.ext[:TransData].infvar_to_supports[x] = [(0.,), (1.,)]
        @test supports(x) == [(0.,), (1.,)]
    end
end

# Test constraint mapping queries
@testset "Constraint Mapping Queries" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1, supports = [0, 1])
    @infinite_variable(m, x(par))
    @point_variable(m, x(0), x0)
    @hold_variable(m, y)
    data = DiscreteMeasureData(par, [0.5, 0.5], [0, 1])
    @constraint(m, c1, x + y - 2 <= 0)
    @constraint(m, c2, measure(x, data) == 0)
    @constraint(m, c3, x0 + y == 5)
    tm = optimizer_model(m)
    @variable(tm, a)
    @variable(tm, b)
    @variable(tm, c)
    @constraint(tm, tc1, b + a <= 2)
    @constraint(tm, tc2, c + a <= 2)
    @constraint(tm, tc3, 0.5b + 0.5c == 0)
    @constraint(tm, tc4, b + a == 5)
    # test transcription_constraint (InfiniteConstraintRef)
    @testset "transcription_constraint (Infinite)" begin
        # test error
        @test_throws ErrorException transcription_constraint(tm, c1)
        # test normal
        tm.ext[:TransData].infinite_to_constrs[c1] = [tc1, tc2]
        @test transcription_constraint(tm, c1) == [tc1, tc2]
    end
    # test transcription_constraint (MeasureConstraintRef)
    @testset "transcription_constraint (Measure)" begin
        # test error
        @test_throws ErrorException transcription_constraint(tm, c2)
        # test normal
        tm.ext[:TransData].measure_to_constrs[c2] = [tc3]
        @test transcription_constraint(tm, c2) == [tc3]
    end
    # test transcription_constraint (FiniteConstraintRef)
    @testset "transcription_constraint (Finite)" begin
        # test error
        @test_throws ErrorException transcription_constraint(tm, c3)
        # test normal
        tm.ext[:TransData].finite_to_constr[c3] = tc4
        @test transcription_constraint(tm, c3) == tc4
    end
    # test supports for infinite constraint with 2 inputs
    @testset "supports (Model, Infinite)" begin
        # test error
        @test_throws ErrorException supports(tm, c1)
        # test normal
        tm.ext[:TransData].infconstr_to_supports[c1] = [(0.,), (1.,)]
        @test supports(tm, c1) == [(0.,), (1.,)]
        # undo changes
        delete!(tm.ext[:TransData].infconstr_to_supports, c1)
    end
    # test supports for infinite constraint with 1 input
    @testset "supports (Infinite)" begin
        # test error
        @test_throws ErrorException supports(c1)
        # test normal
        tm.ext[:TransData].infconstr_to_supports[c1] = [(0.,), (1.,)]
        @test supports(c1) == [(0.,), (1.,)]
    end
    # test supports for measure constraint with 2 inputs
    @testset "supports (Model, Measure)" begin
        # test error
        @test_throws ErrorException supports(tm, c2)
        # test normal
        tm.ext[:TransData].measconstr_to_supports[c2] = [(0.,), (1.,)]
        @test supports(tm, c2) == [(0.,), (1.,)]
        # undo changes
        delete!(tm.ext[:TransData].measconstr_to_supports, c2)
    end
    # test supports for measure constraint with 1 input
    @testset "supports (Measure)" begin
        # test error
        @test_throws ErrorException supports(c2)
        # test normal
        tm.ext[:TransData].measconstr_to_supports[c2] = [(0.,), (1.,)]
        @test supports(c2) == [(0.,), (1.,)]
    end
    # test parameter_refs for infinite constraint with 2 inputs
    @testset "parameter_refs (Model, Infinite)" begin
        # test error
        @test_throws ErrorException parameter_refs(tm, c1)
        # test normal
        tm.ext[:TransData].infconstr_to_params[c1] = (par,)
        @test parameter_refs(tm, c1) == (par,)
        # undo changes
        delete!(tm.ext[:TransData].infconstr_to_params, c1)
    end
    # test parameter_refs for infinite constraint with 1 input
    @testset "parameter_refs (Infinite)" begin
        # test error
        @test_throws ErrorException parameter_refs(c1)
        # test normal
        tm.ext[:TransData].infconstr_to_params[c1] = (par,)
        @test parameter_refs(c1) == (par,)
    end
    # test parameter_refs for measure constraint with 2 inputs
    @testset "parameter_refs (Model, Measure)" begin
        # test error
        @test_throws ErrorException parameter_refs(tm, c2)
        # test normal
        tm.ext[:TransData].measconstr_to_params[c2] = (par,)
        @test parameter_refs(tm, c2) == (par,)
        # undo changes
        delete!(tm.ext[:TransData].measconstr_to_params, c2)
    end
    # test parameter_refs for measure constraint with 1 input
    @testset "parameter_refs (Measure)" begin
        # test error
        @test_throws ErrorException parameter_refs(c2)
        # test normal
        tm.ext[:TransData].measconstr_to_params[c2] = (par,)
        @test parameter_refs(c2) == (par,)
    end
end
