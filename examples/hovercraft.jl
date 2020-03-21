"""
Dynamic path planning problem of hovercraft trying to minimize thrust usage
to hit all the waypoints at the alotted target times.
"""

using InfiniteOpt, JuMP, Ipopt, Plots

# Set problem information
max_time = 60
time_points = Vector(0:max_time)
xw = [1 4 6 1; 1 3 0 1] # waypoints
times = [0; 25; 50; 60]

# Initialize the model
m = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# Set the parameters and variables
@infinite_parameter(m, t in [0, max_time])
@infinite_variable(m, x[1:2](t), start = 1) # position
@infinite_variable(m, v[1:2](t), start = 0) # velocity
@infinite_variable(m, u[1:2](t), start = 0) # thruster input

# Set the initial conditions
@BDconstraint(m, initial_velocity[i = 1:2](t == 0), v[i] == 0)

# Manually implement euler scheme for motion equations
@point_variable(m, x[i](time_points[j]), xp[i = 1:2, j = 1:length(time_points)])
@point_variable(m, v[i](time_points[j]), vp[i = 1:2, j = 1:length(time_points)])
@point_variable(m, u[i](time_points[j]), up[i = 1:2, j = 1:length(time_points)])
@constraint(m, [i = 1:2, j = 1:length(time_points)-1], xp[i, j+1] == xp[i, j] + vp[i, j])
@constraint(m, [i = 1:2, j = 1:length(time_points)-1], vp[i, j+1] == vp[i, j] + up[i, j])

# Hit all the waypoints
@BDconstraint(m, [i = 1:2, j = 1:length(times)](t == times[j]), x[i] == xw[i, j])

# Specify the objective
@objective(m, Min, support_sum(u[1]^2 + u[2]^2, t))

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
scatter(xw[1,:], xw[2,:], label = "Waypoints")
plot!(x_opt[1,:], x_opt[2,:], label = "Trajectory")
xlabel!("x1")
ylabel!("x2")
