# # Fishing Optimal Control
# We will solve an optimal control problem to maximize profit
# constrained by fish population.

# # Problem Statement and Model
# In this system, our input is the fishing rate ``u``, and the profit
# ``J`` will be defined in the objective to maximize.
# Our profit objective is represented as:
# ```math
# \begin{aligned}
# &&\max_{u(t)} J(t) \\
# &&&&&J = \int_0^{10} \left(E - \frac{c}{x}\right) u U_{max} \, dt \\
# &&\text{s.t.} &&& \frac{dx}{dt}= rx(t)\left(1 - \frac{x(t)}{k}\right) - uU_{max}, t \in [0,10] \\
# &&&&&x(0) = 70 \\
# &&&&&0 \leq u(t) \leq 1 \\
# &&&&&E = 1, \; c = 17.5, \; r = 0.71, \; k = 80.5, \; U_{max} = 20 \\
# &&&&&J(0) = 0 \\
# \end{aligned}
# ```

# # Modeling in InfiniteOpt
# First we must import ``InfiniteOpt`` and other packages.
using InfiniteOpt, Ipopt, Plots;

# Next we specify our initial conditions and problem variables.
x0 = 70
E, c, r, k, Umax = 1, 17.5, 0.71, 80.5, 20;

# # Model Initialization
# We initialize the infinite model with [`InfiniteModel`](@ref) and select Ipopt as our optimizer that will be used to solve it.
m = InfiniteModel(Ipopt.Optimizer);

# # Infinite Parameter Definition
# We now define the infinite parameter ``t \in [0, 10]`` to represent time over a 10 year period. We'll also specify 100 equidistant time points.
@infinite_parameter(m, t in [0, 10], num_supports=100)

# # Infinite Variable Definition
# Now that we have our infinite parameter defined, let's specify our infinite variables:
# - ``1 \leq x(t)`` : fish population at time ``t``
# - ``0 \leq u(t) \leq 1`` : fishing rate
# - ``J(t)`` : profit over time
@variable(m, 1 <= x, Infinite(t))
@variable(m, 0 <= u <= 1, Infinite(t));
@variable(m, J, Infinite(t));

# # Objective Definition
# Now we add the objective using `@objective` to maximize profit ``J`` at the end of the 10 year period:
@objective(m, Max, J(10));

# # Constraint Definition
# The last step is to add our constraints. First, define the ODEs which serve as our system model.
@constraint(m, ∂(J,t) == (E-c/x) * u * Umax)
@constraint(m, ∂(x,t) == r * x *(1 - x/k) - u*Umax);
# We also set our initial conditions for ``x`` and ``J``.
@constraint(m, x(0) == x0)
@constraint(m, J(0) == 0);

# # Problem Solution
# Now we're ready to solve! We can solve the model by invoking `optimize!`:
optimize!(m)

# # Extract and Plot the Results
# Now we can extract the optimal results and plot them to visualize how the fish population compares to the profit. Note that the values of infinite variables are
# returned as arrays corresponding to how the supports were used to discretize our model.
# We can use the `value` function to extract the values of the infinite variables.
ts = value(t)
u_opt = value(u)
x_opt = value(x)
J_opt = value(J);

# Create the plot for the fish population and profit over time.
p1 = plot(ts, [x_opt, J_opt] ,
    label=["Fish Pop" "Profit"], 
    title="State Variables");
# Create the plot for the fishing rate over time.
p2 = plot(ts, u_opt,
    label = "Fishing Rate",
    title = "Fishing Rate vs Time");
# Visualize the two plots on one figure.
plot(p1,p2 ,layout=(2,1), size=(800,600))

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
tol = 1E-6
@test termination_status(m) == MOI.LOCALLY_SOLVED
@test has_values(m)
@test u_opt isa Vector{<:Real}
@test J_opt isa Vector{<:Real}
@test isapprox(u_opt[end], 1.0000000088945395, atol=tol)
@test isapprox(x_opt[end], 31.441105707837544, atol=tol)
@test isapprox(J_opt[end], 106.80870543718251, atol=tol)