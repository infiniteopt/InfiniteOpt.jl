using InfiniteOpt, Ipopt, Interpolations

<<<<<<< HEAD
function test_infiniteInterpolate()
    model = InfiniteModel(Ipopt.Optimizer)
    @infinite_parameter(model, t ∈ [1, 2], num_supports = 9)
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
    @test xDiscrete isa Vector{<:Real}
    @test yDiscrete isa Matrix{<:Real}
    @test zDiscrete isa Vector{<:Vector{<:Real}}
    @test wDiscrete isa Vector{<:Matrix{<:Real}}
    @test_throws ArgumentError value.(z, cubic_spline_interpolation)
    @test_throws ArgumentError value.(w, cubic_spline_interpolation)

    @test isapprox(xFunc(1.55), 2.95589105143786, atol=tol)
    @test isapprox(yFunc(1.55, 3.6), 0.07670311557146893, atol=tol)
    @test isapprox(zFunc[1](5.4), 0.7451798867405894, atol=tol)
    @test isapprox(zFunc[2](5.75), 0.29999999442572844, atol=tol)
    @test isapprox(wFunc[1](5.75, 7.35), 0.7606595292931279, atol=tol)
    @test isapprox(wFunc[2](5.75, 7.35), 1.1999999847259393, atol=tol)
end

@testset "InfiniteInterpolate" begin
    test_infiniteInterpolate()
end
=======
# TODO: convert code into function to include in runtests.jl
model = InfiniteModel(Ipopt.Optimizer)
@infinite_parameter(model, t ∈ [1, 2], num_supports = 9)
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

# TODO: After finishing development, remove all code related to plotting & preparing arrays for plotting
# Extract the results
xDiscrete = value(x)
yDiscrete = value(y)
zDiscrete = value.(z)
wDiscrete = value.(w)
tDiscrete = value(t)
sDiscrete = value(s)

# Create matrices for graphing y(t, s)
t_points = repeat(tDiscrete, inner=length(sDiscrete))
s_points = repeat(sDiscrete, outer=length(tDiscrete))
y_points = vec(yDiscrete)

# Test the interpolation value function
xFunc = value(x, cubic_spline_interpolation)
yFunc = value(y, linear_interpolation)
zFunc = value.(z, linear_interpolation)
wFunc = value.(w, linear_interpolation)

# Plot the interpolated values (just for dev purposes)
t_vals = collect(LinRange(1, 2, 100))
s_vals = collect(LinRange(3, 4, 100))
x_vals = xFunc(t_vals)
y_vals = yFunc(t_vals, s_vals);
xPlot = plot(t_vals, x_vals, label="x(t)", xlabel="t", ylabel="x(t)", title="Interpolated x(t)")
plot!(tDiscrete, xDiscrete, seriestype=:scatter, label="Data points", markersize=5, color=:blue)
yPlot = surface(t_vals, s_vals, y_vals, label="y(t,s)", xlabel="t", ylabel="s", zlabel="y(t,s)", title="Interpolated y(t,s)")
scatter!(t_points, s_points, y_points, label="Data points", markersize=5, color=:red)
display(xPlot)
display(yPlot)

# Unit tests
using Test
tol = 1e-06
@test termination_status(model) == MOI.ALMOST_LOCALLY_SOLVED
@test has_values(model)
@test xDiscrete isa Vector{<:Real}
@test yDiscrete isa Matrix{<:Real}
@test zDiscrete isa Vector{<:Vector{<:Real}}
@test wDiscrete isa Vector{<:Matrix{<:Real}}
@test_throws ArgumentError value.(z, cubic_spline_interpolation)
@test_throws ArgumentError value.(w, cubic_spline_interpolation)

@test isapprox(xFunc(1.55), 2.95589105143786, atol=tol)
@test isapprox(yFunc(1.55, 3.6), 0.07670311557146893, atol=tol)
@test isapprox(zFunc[1](5.4), 0.7451798867405894, atol=tol)
@test isapprox(zFunc[2](5.75), 0.29999999442572844, atol=tol)
@test isapprox(wFunc[1](5.75, 7.35), 0.7606595292931279, atol=tol)
@test isapprox(wFunc[2](5.75, 7.35), 1.1999999847259393, atol=tol)
>>>>>>> b51f3bf (Added unit tests for catching irregular grid method errors and checking interpolated values)
