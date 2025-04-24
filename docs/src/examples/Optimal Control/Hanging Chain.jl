# # Hanging Chain Problem
# We will solve a reformulated version of the Hanging Chain Problem
# from the Constrained Optimization Problem Set.
# # Problem Statement and Model
# In this problem, we seek to minimize the potential energy of a chain of
# length ``L`` suspended between points ``a`` and ``b``.
# The potential energy constrained by length is represented by:
# ```math
# \begin{gathered}
# \int_0^1 x(1 + (x')^2)^{1/2} \, dt
# \end{gathered}
# ```
# Our optimization problem is defined as follows:
# ```math
# \begin{aligned}
# &&\min_{x,u} x_2(t_f)\\
# &&\text{s.t.} &&& \quad \dot{x}_1= u\\
# &&&&&\dot{x}_2 = x_1(1+u^2)^{1/2}\\
# &&&&&\dot{x}_3 = (1+u^2)^{1/2}\\
# &&&&&x(t_0) = (a,0,0)^T\\
# &&&&&x_1(t_f) = b\\
# &&&&&x_3(t_f) = L\\
# &&&&&x(t) \in [0,10]\\
# &&&&&u(t) \in [-10,20]\\
# \end{aligned}
# ```

# # Model Definition
# First we must import ``InfiniteOpt`` and other packages.
using InfiniteOpt, Ipopt, Plots;
# Next we specify an array of initial conditions as well as problem variables.
a, b, L = 1, 3, 4
x0 = [a, 0, 0]
xtf = [b, NaN, L];
# We initialize the infinite model and opt to use the Ipopt solver
m = InfiniteModel(Ipopt.Optimizer);
# t is specified as ``\ t \in [0,1]``. The bounds are arbitrary for this problem:
@infinite_parameter(m, t in [0,1], num_supports = 100);
# Now let's specify variables. ``u`` is our controller variable.
@variable(m, 0 <= x[1:3] <= 10, Infinite(t))
@variable(m, -10 <= u <= 20, Infinite(t));
# Specifying the objective to minimize kinetic energy at the final time:
@objective(m, Min, x[2](1));
# Define the ODEs which serve as our system model.
@constraint(m, ∂(x[1],t) == u)
@constraint(m, ∂(x[2],t) == x[1] * (1 + u^2)^(1/2))
@constraint(m, ∂(x[3],t) == (1 + u^2)^(1/2));
# Set our inital and final conditions.
@constraint(m, [i in 1:3], x[i](0) == x0[i])
@constraint(m, x[1](1) == xtf[1])
@constraint(m, x[3](1) == xtf[3]);
# # Problem Solution
# Optimize the model:
optimize!(m)
# Extract the results.
ts = value(t)
u_opt = value(u)
x1_opt = value(x[1])
x2_opt = value(x[2])
x3_opt = value(x[3])
@show(objective_value(m))
p1 = plot(ts, [x1_opt, x2_opt, x3_opt], 
    label=["x1" "x2" "x3"], 
    title="State Variables")

p2 = plot(ts, u_opt, 
    label="u(t)", 
    title="Input")
plot(p1, p2, layout=(2,1), size=(800,600));

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
@test termination_status(m) == MOI.LOCALLY_SOLVED
@test has_values(m)
@test u_opt isa Vector{<:Real}
@test x1_opt isa Vector{<:Real}