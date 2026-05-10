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

# # Modeling in InfiniteOpt 
# First we must import ``InfiniteOpt`` and other packages.
using InfiniteOpt, Ipopt, Plots;

# Next we specify our initial conditions and problem variables.
a, b, L = 1, 3, 4
x0 = [a, 0, 0]
xtf = [b, NaN, L];

# # Model Definition
# We initialize the infinite model with [`InfiniteModel`](@ref) and select Ipopt as our optimizer that will be used to solve it.
m = InfiniteModel(Ipopt.Optimizer);

# # Infinite Parameter Definition
# We define the infinite parameter t as ``\ t \in [0,1]``. The bounds are arbitrary for this problem:
@infinite_parameter(m, t in [0,1], num_supports = 100);

# # Infinite Variable Definition
# Now let's specify variables. We define ``x`` as a vector containing the state variables for position, energy and chain length. ``u`` is our control variable.
@variable(m, 0 <= x[1:3] <= 10, Infinite(t))
@variable(m, -10 <= u <= 20, Infinite(t));

# # Objective Definition
# We specify the objective using `@objective` to minimize potential energy at the final time:
@objective(m, Min, x[2](1));

# # Constraint Definition
# The last step is to add our constraints. First, define the ODEs which serve as our system model.
@constraint(m, ∂(x[1],t) == u)
@constraint(m, ∂(x[2],t) == x[1] * (1 + u^2)^(1/2))
@constraint(m, ∂(x[3],t) == (1 + u^2)^(1/2));

# We also set our initial and final conditions for ``x``.
@constraint(m, [i in 1:3], x[i](0) == x0[i])
@constraint(m, x[1](1) == xtf[1])
@constraint(m, x[3](1) == xtf[3]);

# # Problem Solution
# Now we can solve the model with `optimize!`:
optimize!(m)

# # Extract and Plot the Results
# Extract the results using `value`. Note that they are returned as arrays corresponding to the supports used to discretize the model.
ts = value(t)
u_opt = value(u)
x1_opt = value(x[1])
x2_opt = value(x[2])
x3_opt = value(x[3]);

# We can also check the objective value:
@show(objective_value(m))

# Create the plot for the state variables and input over time.
p1 = plot(ts, [x1_opt, x2_opt, x3_opt], 
    
    label=["x1" "x2" "x3"], 
    title="State Variables")

p2 = plot(ts, u_opt, 
    label="u(t)", 
    title="Input");

# Visualize the two plots on one figure.
plot(p1, p2, layout=(2,1), size=(800,600))

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
tol = 1E-6
@test termination_status(m) == MOI.LOCALLY_SOLVED
@test has_values(m)
@test u_opt isa Vector{<:Real}
@test x1_opt isa Vector{<:Real}
@test isapprox(objective_value(m), 5.127030122851338, atol=tol)
@test isapprox(u_opt[end], 7.20355021172144, atol=tol)