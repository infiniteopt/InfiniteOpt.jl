# # Minimizing Final Time (Jennings Problem)
# Solve an optimal control problem with a minimal final time. 
# Set up and solve the Jennings optimal control benchmark problem.
# # Problem Statement and Model
# When solving differential equations over a variable time interval ``[0,t_f]``,
# we can apply a time-scaling transformation to normalize the interval to``[0,1]``.
# This is achieved by introducing a final time parameter ``t_f``.
# The Jennings optimal control problem divides derivatives by ``t_f``.
# In practice, ``t_f`` appears on the right hand side to avoid any divisions by 0.
# ```math
# \begin{gathered}
# \frac{\frac{dx}{dt}}{t_f} =  f(x,u) \\
# \frac{dx}{dt} = t_f f(x,u) \\
# \end{gathered}
# ```
# Our specific problem is defined as the following:
# ```math
# \begin{aligned}
# &&\min_{u(t),t_f} t_f \\
# &&\text{s.t.} &&& \frac{dx_1}{dt}= t_f u, && t \in [0,1] \\
# &&&&&\frac{dx_2}{dt} = t_f \cos(x_1(t)), && t \in [0,1] \\
# &&&&&\frac{dx_3}{dt} = t_f \sin(x_1(t)), && t \in [0,1] \\
# &&&&&x(0) = [\pi/2, 4, 0] \\
# &&&&&x_2(t_f) = 0 \\
# &&&&&x_3(t_f) = 0 \\
# &&&&&-2 \leq u(t) \leq 2
# \end{aligned}
# ```

# # Modeling in InfiniteOpt
# First we must import ``InfiniteOpt`` and other packages.
using InfiniteOpt, Ipopt, Plots;

# Next we specify an array of initial conditions.
x0 = [π/2, 4, 0]; # x(0) for x1, x2, x3

# # Model Definition
# We initialize the infinite model with [`InfiniteModel`](@ref), opting to use the Ipopt solver.
m = InfiniteModel(Ipopt.Optimizer);

# # Infinite Parameter Definition
# Recall t is specified as ``\ t \in [0,1]``:
@infinite_parameter(m, t in [0,1],num_supports= 100);

# # Variable Definition
# Now let's specify decision variables. Notice that ``t_f`` is
# not a function of time and is a singular value.
@variable(m, x[1:3], Infinite(t))
@variable(m, -2 <= u <= 2, Infinite(t))
@variable(m, 0.1 <= tf);

# # Objective Definition
# Now we add the objective using `@objective` to minimize final time:
@objective(m, Min, tf);

# # Constraint Definition
# The last step is to add our constraints. First, define the ODEs which serve as our system model.
@constraint(m, ∂(x[1],t) == tf*u)
@constraint(m, ∂(x[2],t) == tf*cos(x[1]))
@constraint(m, ∂(x[3],t) == tf*sin(x[1]));

# Set our inital and final conditions for ``x``.
@constraint(m, [i in 1:3], x[i](0) == x0[i])
@constraint(m, x[2](1) <=0)
@constraint(m, x[3](1) <= 1e-1);

# # Problem Solution
# Now everything is ready for solving! We can solve the model with `@optimize!`:
optimize!(m)

# # Extract and Plot the Results
# We can extract the results as arrays using the `value` function. Notice that we multiply by ``t_f``
# to scale our time.
ts = value(t)*value(tf)
u_opt = value(u)
x1_opt = value(x[1])
x2_opt = value(x[2])
x3_opt = value(x[3]);

# Plot the results
plot(ts, u_opt, label = "u(t)", linecolor = :black, linestyle = :dash)
plot!(ts, x1_opt, linecolor = :blue, linealpha = 0.4, label = "x1")
plot!(ts, x2_opt, linecolor = :green, linealpha = 0.4, label = "x2")
plot!(ts, x3_opt, linecolor = :red, linealpha = 0.4, label = "x3")

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
tol = 1E-6
@test termination_status(m) == MOI.LOCALLY_SOLVED
@test has_values(m)
@test u_opt isa Vector{<:Real}
@test x1_opt isa Vector{<:Real}
@test isapprox(x1_opt[end], 3.2501431326448293, atol=tol)
@test isapprox(objective_value(m), 4.284564834847627, atol=tol)