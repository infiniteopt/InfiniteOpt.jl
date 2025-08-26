using Interpolations # TODO: Move to runtests.jl once OffSetArrays type piracy problem is resolved 

function test_infiniteInterpolate()
    # Set up model with mock optimizer
    Optimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()), eval_objective_value=false)
    model = InfiniteModel(Optimizer)
    tb = model.backend

    # Define infinite parameters and variables
    @infinite_parameter(model, t ∈ [1, 2], num_supports = 3)
    @infinite_parameter(model, s ∈ [3, 4], num_supports = 3, derivative_method = OrthogonalCollocation(3))
    @infinite_parameter(model, γ ∈ [5, 6], supports = [5, 5.15, 5.25, 6])
    @infinite_parameter(model, δ ∈ [7, 8], num_supports = 3)
    @variable(model, 0 ≤ x ≤ 10, Infinite(t))                       # single parameter
    @variable(model, 0 ≤ y ≤ 10, Infinite(t, s))                    # multi parameter
    @variable(model, 0.3 ≤ z[1:2] ≤ 8, Infinite(γ))                 # vector, single parameter
    @variable(model, 0.1 ≤ w[1:2] ≤ 8, Infinite(γ, δ))              # vector, multi-parameter
    @variable(model, ySemi, SemiInfinite(y, 1.5, s))       # semi-infinite
    @variable(model, wSemi[i = 1:2], SemiInfinite(w[i], 5.25, δ))  # vector, semi-infinite
    @variable(model, α, Point(y, 1.5, 3.5))                # point variable
    @variable(model, 0 ≤ ω ≤ 10)                                    # finite variable

    @constraint(model, 8*∂(x, t)^2 + cos(∂(y,s))^3 ≤ 5*x*y)
    @constraint(model, sin(x)*t ≥ y*s)
    @constraint(model, [i in 1:2], z[i] + w[i] ≤ 1.5)
    @objective(model, Max, ∫(x^2, t) + ∫(∫(y^2 + x*y, t), s) + ∫(∫(z[1] - sin(z[2]) + w[1] - cos(w[2]), γ), δ))
    set_silent(model)
    JuMP.optimize!(model)

    # Set up the results
    mockOptimizer = JuMP.backend(tb).optimizer.model
    MOI.set(mockOptimizer, MOI.TerminationStatus(), MOI.ALMOST_LOCALLY_SOLVED)

    # Set up the transformation variables
    xVar = transformation_variable(x, label = InfiniteOpt.All)
    yVar = transformation_variable(y, label = InfiniteOpt.All)
    zVar = transformation_variable.(z, label = InfiniteOpt.All)
    wVar = transformation_variable.(w, label = InfiniteOpt.All)
    ySemiVar = transformation_variable(ySemi, label = InfiniteOpt.All)
    wSemiVar = transformation_variable.(wSemi, label = InfiniteOpt.All)
    αVar = transformation_variable(α, label = InfiniteOpt.All)
    ωVar = transformation_variable(ω, label = InfiniteOpt.All)
    dxVar = transformation_variable(∂(x, t), label = InfiniteOpt.All)
    dyVar = transformation_variable(∂(y, s), label = InfiniteOpt.All)

    xVals = [2.858, 2.929, 2.991]

    yVals = [0.093 0.079 0.069;
            0.098 0.086 0.075;
            0.106 0.090 0.079;
            0.115 0.088 0.078;
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

    ySemiVals = [0.193 0.179 0.169 0.165 0.135]
    
    wSemiVals = [[0.620 0.625 0.627],
            [2.199 2.200 2.121]]
    
    dxVals = [1.658, 1.629, 1.691]
    
    dyVals = [2.093 2.079 2.069;
            2.095 2.086 2.076;
            2.106 2.090 2.079;
            2.115 2.096 2.082;
            2.100 2.086 2.075]

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

    for i in eachindex(ySemiVals)
        MOI.set(
            mockOptimizer,
            MOI.VariablePrimal(1),
            JuMP.optimizer_index(ySemiVar[i]),
            ySemiVals[i]
        )
    end

    for i in 1:2
        for j in eachindex(wSemiVals)
            MOI.set(
                mockOptimizer,
                MOI.VariablePrimal(1),
                JuMP.optimizer_index(wSemiVar[i][j]),
                wSemiVals[i][j]
            )
        end
    end

    MOI.set(mockOptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(αVar), 0.090)
    MOI.set(mockOptimizer, MOI.VariablePrimal(1), JuMP.optimizer_index(ωVar), 3.5)

    for i in eachindex(dxVals)
        MOI.set(
            mockOptimizer,
            MOI.VariablePrimal(1),
            JuMP.optimizer_index(dxVar[i]),
            dxVals[i]
        )
    end

    for i in eachindex(dyVals)
        MOI.set(
            mockOptimizer,
            MOI.VariablePrimal(1),
            JuMP.optimizer_index(dyVar[i]),
            dyVals[i]
        )
    end

    # Test the interpolation value function
    xFunc = value(x, cubic_spline_interpolation)
    yFunc = value(y, linear_interpolation)
    zFunc = value.(z, linear_interpolation)
    wFunc = value.(w, linear_interpolation)
    ySemiFunc = value(ySemi, cubic_spline_interpolation)
    wSemiFunc = value.(wSemi, linear_interpolation)
    αValue = value(α, constant_interpolation)
    ωValue = value(ω, constant_interpolation)
    dxFunc = value(∂(x, t), cubic_spline_interpolation)
    dyFunc = value(∂(y, s), linear_interpolation)

    # Unit tests
    tol = 1e-06
    @test termination_status(model) == MOI.ALMOST_LOCALLY_SOLVED
    @test value(x) isa Vector{<:Real}
    @test value(y) isa Matrix{<:Real}
    @test value.(z) isa Vector{<:Vector{<:Real}}
    @test value.(w) isa Vector{<:Matrix{<:Real}}
    @test value(α) isa Real
    @test value(ω) isa Real

    # Create alternative interpolation method
    quad_interpolation(params, supps) = interpolate(params, supps, Gridded(Quadratic()))

    # Test interpolation method errors
    @test_throws ArgumentError value(x, quad_interpolation)
    @test_throws ArgumentError value.(z, cubic_spline_interpolation)
    @test_throws ArgumentError value.(w, cubic_spline_interpolation)

    # Test the interpolation values
    @test isapprox(xFunc(1.55), 2.9355847500000003, atol=tol)
    @test isapprox(yFunc(1.55, 3.6), 0.09764, atol=tol)
    @test isapprox(zFunc[1](5.4), 0.5825999999999999, atol=tol)
    @test isapprox(zFunc[2](5.75), 0.29766666666666663, atol=tol)
    @test isapprox(wFunc[1](5.75, 7.35), 0.8185666666666667, atol=tol)
    @test isapprox(wFunc[2](5.75, 7.35), 1.5328333333333333, atol=tol)
    @test isapprox(ySemiFunc(3.24), 0.12989190399999995, atol=tol)
    @test isapprox(wSemiFunc[1](7.24), 0.6224000000000001, atol=tol)
    @test isapprox(wSemiFunc[2](7.24), 2.19948, atol=tol)
    @test αValue == 0.090
    @test ωValue == 3.5
    @test isapprox(dxFunc(1.35), 1.62957825, atol=tol)
    @test isapprox(dyFunc(1.35, 3.86), 2.083256, atol=tol)
end

@testset "InfiniteInterpolate" begin
    test_infiniteInterpolate()
end