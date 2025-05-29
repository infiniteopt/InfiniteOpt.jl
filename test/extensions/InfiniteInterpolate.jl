using InfiniteOpt, Ipopt, Interpolations, Plots

model = InfiniteModel(Ipopt.Optimizer)
@infinite_parameter(model, t ∈ [1, 2], num_supports = 5)
@infinite_parameter(model, s ∈ [3, 4], num_supports = 5)
@variable(model, 0 ≤ x ≤ 10, Infinite(t))    # For single variable case
@variable(model, 0 ≤ y ≤ 10, Infinite(t, s))  # For multi-variable case

# TODO: How do we want to return vector variables?
# @variable(model, 0 ≤ z[1:2] ≤ 8, Infinite(t))  # For vector, single variable case
# @variable(model, 0 ≤ w[1:2] ≤ 8, Infinite(t, s))  # For vector, multi-variable case

# Multi-variable, multi-parameter case?

@constraint(model, ∂(x, t)^2 + cos(∂(y,s,s))^3 ≤ 5*x*y)
@constraint(model, sin(x)*t ≥ y*s)
@objective(model, Max, ∫(x^2, t) + ∫(∫(y^2 + x*y, t), s))
set_silent(model)
optimize!(model)

# Test the regular value function
xDiscrete = value(x, false)
yDiscrete = value(y, false)
tDiscrete = value(t, false)
sDiscrete = value(s, false)

# Create matrices for graphing y(t, s)
t_points = repeat(tDiscrete, inner=length(sDiscrete))
s_points = repeat(sDiscrete, outer=length(tDiscrete))
y_points = vec(yDiscrete)

# Test the interpolation value function
xFunc = value(x, true)
yFunc = value(y, true)
# println("interpolated x = $(xFunc)")
# println("interpolated y = $(yFunc)")

# Plot the interpolated values
t_vals = collect(LinRange(1, 2, 100))
s_vals = collect(LinRange(3, 4, 100))
x_vals = xFunc(t_vals)
y_vals = yFunc(t_vals, s_vals)
xPlot = plot(t_vals, x_vals, label="x(t)", xlabel="t", ylabel="x(t)", title="Interpolated x(t)")
plot!(tDiscrete, xDiscrete, seriestype=:scatter, label="Data points", markersize=5, color=:blue)
yPlot = surface(t_vals, s_vals, y_vals, label="y(t,s)", xlabel="t", ylabel="s", zlabel="y(t,s)", title="Interpolated y(t,s)")
scatter!(t_points, s_points, y_points, label="Data points", markersize=5, color=:red)
display(xPlot)
display(yPlot)