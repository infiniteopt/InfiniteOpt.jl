using InfiniteOpt, Ipopt, Interpolations

function test_infiniteInterpolate()
    model = InfiniteModel(Ipopt.Optimizer)
    @infinite_parameter(model, t ∈ [1, 2], num_supports = 5)
    @infinite_parameter(model, s ∈ [3, 4], num_supports = 5)
    @infinite_parameter(model, γ ∈ [5, 6], supports = [5, 5.15, 5.25, 5.5, 6])
    @infinite_parameter(model, δ ∈ [7, 8], num_supports = 5)
    @variable(model, 0 ≤ x ≤ 10, Infinite(t))    # For single variable case
    @variable(model, 0 ≤ y ≤ 10, Infinite(t, s))  # For multi-variable case
    @variable(model, 0.3 ≤ z[1:2] ≤ 8, Infinite(γ))  # For vector, single variable case
    @variable(model, 0.1 ≤ w[1:2] ≤ 8, Infinite(γ, δ))  # For vector, multi-variable case

    @constraint(model, 8*∂(x, t)^2 + cos(∂(y,s,s))^3 ≤ 5*x*y)
    @constraint(model, sin(x)*t ≥ y*s)
    @constraint(model, [i in 1:2], z[i] + w[i] ≤ 1.5)
    @objective(model, Max, ∫(x^2, t) + ∫(∫(y^2 + x*y, t), s) + ∫(∫(z[1] - sin(z[2]) + w[1] - cos(w[2]), γ), δ))
    set_silent(model)
    optimize!(model)

    # Test the interpolation value function
    xFunc = value(x, cubic_spline_interpolation)
    yFunc = value(y, linear_interpolation)
    zFunc = value.(z, linear_interpolation)
    wFunc = value.(w, linear_interpolation)

    # Unit tests
    tol = 1e-06
    @test termination_status(model) == MOI.ALMOST_LOCALLY_SOLVED
    @test has_values(model)
    @test value.(x) isa Vector{<:Real}
    @test value.(y) isa Matrix{<:Real}
    @test value.(z) isa Vector{<:Vector{<:Real}}
    @test value.(w) isa Vector{<:Matrix{<:Real}}

    # Create alternative interpolation method
    quad_interpolation(params, supps) = interpolate(params, supps, Gridded(Quadratic()))

    # Test interpolation method errors
    @test_throws ArgumentError value(x, quad_interpolation)
    @test_throws ArgumentError value.(z, cubic_spline_interpolation)
    @test_throws ArgumentError value.(w, cubic_spline_interpolation)

    # Test the interpolation values
    @test isapprox(xFunc(1.55), 2.9560057246495184, atol=tol)
    @test isapprox(yFunc(1.55, 3.6), 0.07663947385950569, atol=tol)
    @test isapprox(zFunc[1](5.4), 0.80273946541415354, atol=tol)
    @test isapprox(zFunc[2](5.75), 0.29999999442572844, atol=tol)
    @test isapprox(wFunc[1](5.75, 7.35), 0.7013948501185225, atol=tol)
    @test isapprox(wFunc[2](5.75, 7.35), 1.1999999847259393, atol=tol)
end

@testset "InfiniteInterpolate" begin
    test_infiniteInterpolate()
end