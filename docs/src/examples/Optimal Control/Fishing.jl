# # Fishing Optimal Control
# We will solve an optimal control problem to maximize profit
# constrained by fish population
# # Problem Statement and Model
# In this system out input is the rate of fishing ``u``, and the profit
# ``J`` will be defined in the objective to maximize.
# Our profit objective is represented as:
# ```math
# \begin{aligned}
# &&\max_{u(t)} J(t) \\
# &&&&&J = \int_0^{10} \left(E - \frac{c}{x}\right) u U_{max} \, dt \\
# &&\text{s.t.} &&& \frac{dx}{dt}= rx(t)\left(1 - \frac{x(t)}{k}\right) - uU_{max} \\
# &&&&&x(0) = 70 \\
# &&&&&0 \leq u(t) \leq 1 \\
# &&&&&E = 1, \; c = 17.5, \; r = 0.71, \; k = 80.5, \; U_{max} = 20 \\
# &&&&&J(0) = 0 \\
# \end{aligned}
# ```
# # Model Definition
# First we must import ``InfiniteOpt`` and other packages.
using InfiniteOpt, Ipopt, Plots;
# Next we specify an array of initial conditions as well as problem variables.
x0 = 70
E, c, r, k, Umax = 1, 17.5, 0.71, 80.5, 20;
# We initialize the infinite model and opt to use the Ipopt solver
m = InfiniteModel(Ipopt.Optimizer);
# Now let's specify variables. ``u`` is as our fishing rate.
# ``x`` will be used to model the fish population in response to ``u``
# the infinite parameter ``t`` that will span over 10 years.
@infinite_parameter(m, t in [0,10],num_supports=100)
@variable(m, 1 <= x, Infinite(t))
@variable(m, 0 <= u <= 1, Infinite(t));
# ``J`` represents profit over time.
@variable(m, J, Infinite(t));
# Specifying the objective to maximize profit ``J``:
@objective(m, Max, J(10));
# Define the ODEs which serve as our system model.
@constraint(m, ∂(J,t) == (E-c/x) * u * Umax)
@constraint(m, ∂(x,t) == r * x *(1 - x/k) - u*Umax);
# Set our initial conditions.
@constraint(m, x(0) == x0)
@constraint(m, J(0) == 0);
# # Problem Solution
# Optimize the model:
optimize!(m)
# Extract the results.
ts = value(t)
u_opt = value(u)
x_opt = value(x)
J_opt = value(J);

p1 = plot(ts, [x_opt, J_opt] ,
    label=["Fish Pop" "Profit"], 
    title="State Variables")
p2 = plot(ts, u_opt,
    label = "Rate",
    title = "Rate vs Time")
plot(p1,p2 ,layout=(2,1), size=(800,600));
