"""
Dynamic path planning problem of hovercraft trying to minimize thrust usage
to hit all the waypoints at the alotted target times.
"""

using InfiniteOpt, Ipopt, PyPlot

# Waypoints
xw = [1 4 6 1; 1 3 0 1] # positions
tw = [0, 25, 50, 60]    # times

# Initialize the model
m = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# Set the parameters and variables
@infinite_parameter(m, t in [0, 60], num_supports = 61)
@infinite_variable(m, x[1:2](t), start = 1) # position
@infinite_variable(m, v[1:2](t), start = 0) # velocity
@infinite_variable(m, u[1:2](t), start = 0) # thruster input

# Specify the objective
@objective(m, Min, integral(u[1]^2 + u[2]^2, t))

# Set the initial conditions
@BDconstraint(m, initial_velocity[i = 1:2](t == 0), v[i] == 0)

# Define the point physics
@constraint(m, [i = 1:2], deriv(x[i], t) == v[i])
@constraint(m, [i = 1:2], deriv(v[i], t) == u[i])

# Hit all the waypoints
@BDconstraint(m, [i = 1:2, j = eachindex(tw)](t == tw[j]), x[i] == xw[i, j])

# Optimize the model
optimize!(m)

# Get the results
if has_values(m)
    x_opt = value.(x)
    v_opt = value.(v)
    u_opt = value.(u)
    opt_supports = supports.(x)
    obj_opt = objective_value(m)
end

# Plot the results
figure()
scatter(xw[1,:], xw[2,:], label = "Waypoints")
plot(x_opt[1], x_opt[2], label = "Trajectory", color = "C1")
xlabel(L"$x_1$")
ylabel(L"$x_2$")
legend()
