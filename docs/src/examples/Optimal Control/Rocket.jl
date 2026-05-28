# # Rocket Launch
# We will solve an optimal control problem that aims to minimize final time.

# # Problem Statement and Model
# In this problem, we seek to minimize the final time ``t_f`` of a rocket launch between 0 and 10 metres.
# Minimizing fuel use is also important for this problem. The force- or thrust- of the rocket can be between
# -1.1 and 1.1, while the velocity of the rocket cannot exceed 1.7. The drag resistance of the rocket is
# proportional to the square of the velocity. The mass of the rocket decreases as fuel is burned over time. 
# Our optimization guidelines are subject to these constraints:
# ```math
# \begin{aligned}
#   &&\min t_f \\
#   &&&\frac{ds}{dt}= v \\
#   &&&\frac{dv}{dt}= (u-0.2*(v^2))/2\\
#   &&&\frac{dm}{dt}= -0.01(u^2)\\
#   &&&0 \leq v(t) \leq 1.7 \\
#   &&&-1.1 \leq u(t) \leq 1.1 \\
#       &&&s(0) = 0 \\
#       &&&v(0) = 0 \\
#       &&&m(0) = 1 \\
#       &&&s(t_f) = 10 \\
#       &&&v(t_f) = 0 \\
# \end{aligned}
# ```
# where ``v(t)`` represents the rocket's velocity, ``s(t)`` represents the rocket's position and
# ``m(t)`` represents the rocket's mass. Additionally, ``u(t)`` represents the rocket's force (or thrust).

# # Modeling in InfiniteOpt
# First we must import ``InfiniteOpt`` and other packages.
using JuMP, InfiniteOpt, Ipopt, Plots;

# # Model Initialization
# We then initialize the infinite model with [`InfiniteModel`](@ref), selecting Ipopt as the desired optimizer to solve this model.
model = InfiniteModel(Ipopt.Optimizer);

# # Infinite Parameter Definition
# We are defining the infinite parameter ``t \in [0, 1]``. The bounds give the optimizer parameters to iterate over.
@infinite_parameter(model, t in [0,1], num_supports = 50, derivative_method = FiniteDifference(Backward()));

# # Infinite Variable Definition
# Moving on to specifying variables, we have v as velocity of the rocket, u as the force of the rocket.
# The position, p, and mass, m, vary over time. The final time, t_f, must be between 0.1 and 100 seconds.
@variable(model, 0 <= v <= 1.7, Infinite(t)) 
@variable(model, -1.1 <= u <= 1.1, Infinite(t)) 
@variable(model, 0.1 <= t_f <= 100) 
@variable(model, p, Infinite(t))
@variable(model, m, Infinite(t));

# # Initial and Final Condition Definition
# Now that the model is defined, we then want to specify both our initial and final conditions.
@constraint(model, v(0) == 0)
@constraint(model, v(1) == 0)
@constraint(model, p(0) == 0)
@constraint(model, p(1) == 10)
@constraint(model, m(0) == 1)
@constraint(model, m(1) >= 0.2);

# # Objective Definition
# We now add the objective, using `@objective`, which is to minimize final time t_f.
@objective(model, Min, t_f);

# # Constraint Definition
# We finally are adding our constraints, which are defining the ODEs for our system model.
@constraint(model, v * t_f == deriv(p, t))
@constraint(model, ((u - (0.2 * (v ^ 2))) * t_f) == (deriv(v, t)) * m) 
@constraint(model, (-0.01 * (u ^ 2)) * t_f == deriv(m, t));

# # Problem Solution
# Now that everything is set up, we can solve the model by using `@optimize!`
optimize!(model)

# # Extract and Plot the Results
# We can now extract the optimized results and plot them to visualize four results graphically:
# Position over time, velocity over time, mass over time and force over time. In addition, we can print the final time value of t_f.
tf_data = value(t_f)
t_data = supports(t) .* tf_data
v_data = value(v)
p_data = value(p)
m_data = value(m)
u_data = value(u);

# Printing the final time value.
println("Minimized Final Time: ", tf_data," seconds");

# Creating the four plots as described above. 
plot1 = plot(t_data, v_data, title="Velocity Over Time", linecolor = :red, xlabel="Time", label = "Velocity", ylabel="Velocity", lw=2);
plot2 = plot(t_data, p_data, title="Position Over Time", linecolor = :orange, xlabel="Time", label = "Position", ylabel="Position", lw=2);
plot3 = plot(t_data, m_data, title="Mass Over Time", linecolor = :yellow, xlabel="Time", label = "Position", ylabel="Position", lw=2);
plot4 = plot(t_data, u_data, title="Force Over Time", linecolor = :limegreen, xlabel="Time", label = "Position", ylabel="Position", lw=2);

# Visualizing all four plots on one display.
display(plot(plot1, plot2, plot3, plot4, layout=(4,1), size = (1000, 1000)))

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
tol = 1E-6
@test termination_status(model) == MOI.LOCALLY_SOLVED
@test has_values(model)