function test_infiniteInterpolate()
    Optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), eval_objective_value=false)
    model = InfiniteModel(Optimizer)
    tb = model.backend
    @infinite_parameter(model, t ∈ [1, 2], num_supports = 3)
    @infinite_parameter(model, s ∈ [3, 4], num_supports = 3)
    @infinite_parameter(model, γ ∈ [5, 6], supports = [5, 5.15, 5.25, 6])
    @infinite_parameter(model, δ ∈ [7, 8], num_supports = 3)
    @variable(model, 0 ≤ x ≤ 10, Infinite(t))    # For single variable case
    @variable(model, 0 ≤ y ≤ 10, Infinite(t, s))  # For multi-variable case
    @variable(model, 0.3 ≤ z[1:2] ≤ 8, Infinite(γ))  # For vector, single variable case
    @variable(model, 0.1 ≤ w[1:2] ≤ 8, Infinite(γ, δ))  # For vector, multi-variable case

    # TODO: add semi-infinite variable and point variables for testing

    @constraint(model, 8*∂(x, t)^2 + cos(∂(y,s))^3 ≤ 5*x*y)
    @constraint(model, sin(x)*t ≥ y*s)
    @constraint(model, [i in 1:2], z[i] + w[i] ≤ 1.5)
    @objective(model, Max, ∫(x^2, t) + ∫(∫(y^2 + x*y, t), s) + ∫(∫(z[1] - sin(z[2]) + w[1] - cos(w[2]), γ), δ))
    set_silent(model)
    JuMP.optimize!(model)

    # Set up the results
    mockOptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockOptimizer, MOI.TerminationStatus(), MOI.ALMOST_LOCALLY_SOLVED)

    xVar = transformation_variable(x, label = All)
    yVar = transformation_variable(y, label = All)
    zVar = transformation_variable.(z, label = All)
    wVar = transformation_variable.(w, label = All)

    xVals = [2.858, 2.929, 2.991]

    yVals = [0.093 0.079 0.069;
            0.106 0.090 0.079;
            0.100 0.086 0.075]

    zVals = [[0.575, 0.578, 0.583, 0.581],
            [0.300, 0.300, 0.299, 0.297]]

    wVals = [[0.920 0.925 0.927;
            0.919 0.922 0.928;
            0.913 0.914 0.919;
            0.914 0.917 0.937],
            [1.199 1.200 1.121;
            1.199 1.153 1.435;
            1.199 1.201 1.198;
            1.198 1.200 1.197]]

    for i in eachindex(xVals)
        MOI.set(
            mockOptimizer,
            MOI.VariablePrimal(1),
            JuMP.optimizer_index(xVar[i]),
            xVals[i]
        )
    end

    for i in eachindex(yVals)
        MOI.set(
            mockOptimizer,
            MOI.VariablePrimal(1),
            JuMP.optimizer_index(yVar[i]),
            yVals[i]
        )
    end

    for i in 1:2
        for j in eachindex(zVals[i])
            MOI.set(
                mockOptimizer,
                MOI.VariablePrimal(1),
                JuMP.optimizer_index(zVar[i][j]),
                zVals[i][j]
            )
        end
    end

    for i in 1:2
        for j in eachindex(wVals[i])
            MOI.set(
                mockOptimizer,
                MOI.VariablePrimal(1),
                JuMP.optimizer_index(wVar[i][j]),
                wVals[i][j]
            )
        end
    end

    # Test the interpolation value function
    xFunc = value(x, cubic_spline_interpolation)
    yFunc = value(y, linear_interpolation)
    zFunc = value.(z, linear_interpolation)
    wFunc = value.(w, linear_interpolation)

    # Unit tests
    tol = 1e-06
    @test termination_status(model) == MOI.ALMOST_LOCALLY_SOLVED
    @test value(x) isa Vector{<:Real}
    @test value(y) isa Matrix{<:Real}
    @test value.(z) isa Vector{<:Vector{<:Real}}
    @test value.(w) isa Vector{<:Matrix{<:Real}}

    # Create alternative interpolation method
    quad_interpolation(params, supps) = interpolate(params, supps, Gridded(Quadratic()))

    # Test interpolation method errors
    @test_throws ArgumentError value(x, quad_interpolation)
    @test_throws ArgumentError value.(z, cubic_spline_interpolation)
    @test_throws ArgumentError value.(w, cubic_spline_interpolation)

    # Test the interpolation values
    # println("wFunc[2](5.75,7.35): $(wFunc[2](5.75,7.35))")
    @test isapprox(xFunc(1.55), 2.9355847500000003, atol=tol)
    @test isapprox(yFunc(1.55, 3.6), 0.08739999999999998, atol=tol)
    @test isapprox(zFunc[1](5.4), 0.5825999999999999, atol=tol)
    @test isapprox(zFunc[2](5.75), 0.29766666666666663, atol=tol)
    @test isapprox(wFunc[1](5.75, 7.35), 0.9153, atol=tol)
    @test isapprox(wFunc[2](5.75, 7.35), 1.1997333333333333, atol=tol)
end

@testset "InfiniteInterpolate" begin
    test_infiniteInterpolate()
end