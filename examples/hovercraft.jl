# include("C:/Users/puls446/.julia/dev/InfOpt/examples/hovercraft.jl")

using InfiniteOpt, JuMP, Ipopt, PyPlot

# Set problem information
max_time = 60
time_points = Vector(0:max_time)
xw = [1 4 6 1; 1 3 0 1]
times = [0; 20; 50; 60]

# Initialize the model
m = InfiniteModel(with_optimizer(Ipopt.Optimizer, print_level = 0))

# Set the parameters and variables
@infinite_parameter(m, 0 <= t <= max_time, supports = time_points)
@infinite_variable(m, x[1:2](t)) # position
@infinite_variable(m, v[1:2](t)) # velocity
@infinite_variable(m, u[1:2](t)) # thruster input

# Specify the objective
time_data = DiscreteMeasureData(t, ones(num_supports(t)), supports(t), name = "sum")
@objective(m, Min, measure(u[1]^2 + u[2]^2, time_data))

# Set the initial conditions
initial_time = Dict(t => IntervalSet(0, 0))
@constraint(m, initial_velocity[i = 1:2], v[i] == 0, parameter_bounds = initial_time)

# Newton equations
for time in 0:max_time-1
    x_curr = @point_variable(m, [i = 1:2], infinite_variable_ref = x[i], parameter_values = time)
    x_next = @point_variable(m, [i = 1:2], infinite_variable_ref = x[i], parameter_values = time + 1)
    v_curr = @point_variable(m, [i = 1:2], infinite_variable_ref = v[i], parameter_values = time)
    v_next = @point_variable(m, [i = 1:2], infinite_variable_ref = v[i], parameter_values = time + 1)
    u_curr = @point_variable(m, [i = 1:2], infinite_variable_ref = u[i], parameter_values = time)
    @constraint(m, [i = 1:2], x_next[i] == x_curr[i] + v_curr[i])
    @constraint(m, [i = 1:2], v_next[i] == v_curr[i] + u_curr[i])
end

# Hit all the waypoints
for i = 1:length(times)
    time_bound = Dict(t => IntervalSet(times[i], times[i]))
    @constraint(m, [j = 1:2], x[j] == xw[j, i], parameter_bounds = time_bound)
end

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
plot(x_opt[1,:], x_opt[2,:], "b.-", markersize=4 )
plot( xw[1,:], xw[2,:], "r.", markersize=12 )
axis("equal");
axis((1.,8.,-.5,3.5));
