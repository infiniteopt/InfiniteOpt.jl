# Examples
Here we exemplify the use of `InfiniteOpt` via a few case studies:

## Two-Stage Stochastic Program
First let's consider a standard two-stage stochastic program. Such problems
consider 1st stage variables ``x \in X \subseteq \mathbb{R}^{n_x}`` which denote
upfront (hear-and-now) decisions made before any realization of the random
parameters ``\xi \in \mathbb{R}^{n_\xi}`` is observed, and 2nd stage variables
``y(\xi) \in \mathbb{R}^{n_y}`` which denote recourse (wait-and-see) decisions
that are made in response to realizations of ``\xi``. Moreover, the objective
seeks to optimize 1st stage costs ``f_1(x)`` and second stage costs
``f_2(x, y(\xi))`` which are evaluated over the uncertain domain via a risk
measure ``R_\xi[\cdot]`` (e.g., the expectation ``\mathbb{E}_\xi[\cdot]``).
Putting this together, we obtain the two-stage stochastic program:
```math
\begin{aligned}
	&&\min_{x, y(\xi)} &&& f_1(x) + R_\xi[f_2(x, y(\xi))] \\
	&&\text{s.t.} &&&  g_i(x, y(\xi)) = 0, && i \in I\\
	&&&&& h_j(x, y(\xi)) \leq 0, && j \in J\\
    &&&&&  x \in X\\
\end{aligned}
```
where ``g_i(x, y(\xi)), \ i \in I,`` denote 2nd stage equality constraints,
``h_j(x, y(\xi)), \ j \in J,`` are 2nd stage inequality constraints, and ``X``
denotes the set of feasible 1st stage decisions.

For an example, we consider the classic farmer problem. Here the farmer must
allocate farmland ``x_c`` for each crop ``c \in C`` with random yields per acre
``\xi_c`` such that he minimizes expenses (i.e., maximizes profit) while fulfilling
contractual demand ``d_c``. If needed he can purchase crops from other farmers
to satisfy his contracts. He can also sell extra crop yield that exceeds his
contractual obligations. Thus, here we have 1st stage variables ``x_c`` and
2nd stage variables of crops sold ``w_c(\xi)`` and crops purchased ``y_c(\xi)``.
Putting this together using the expectation ``\mathbb{E}_\xi[\cdot]`` as our
risk measure we obtain:
```math
\begin{aligned}
	&&\underset{x, y(\xi), w(\xi)}{\text{min}} &&& \sum_{c \in C} \alpha_c x_c + \mathbb{E}_{\xi}\left[\sum_{c \in C}\beta_c y_c(\xi) - \lambda_c w_c(\xi)\right] \\
	&&\text{s.t.} &&&  \sum_{c \in C} x_c \leq \bar{x}\\
	&&&&& \xi_c x_c + y_c(\xi) - w_c(\xi) \geq d_c, && c \in C \\
    &&&&& 0 \leq x_c \leq \bar{x}, && c \in C \\
    &&&&& 0 \leq y_c(\xi) \leq \bar{y}_c, && c \in C \\
    &&&&& 0 \leq w_c(\xi) \leq \bar{w}_c, && c \in C \\
    &&&&& \xi_c \in \Xi_c, && c \in C
\end{aligned}
```
where ``\alpha_c`` are production costs, ``\beta_c`` are the purchase prices,
``\lambda_c`` are the selling prices, ``\bar{x}`` is the total acreage,
``\bar{y}_c`` are purchases limits, ``\bar{w}_c`` are selling limits, and
``\Xi_c`` are the underlying distributions.

Now let's implement this first by defining the problem parameters:
```jldoctest 2-stage; output = false
using Distributions

# Model parameters
num_scenarios = 10
C = 1:3

# Data
α = [150, 230, 260] # land cost
β = [238, 210, 0]   # purchasing cost
λ = [170, 150, 36]  # selling price
d = [200, 240, 0]   # contract demand
xbar = 500          # total land
wbar3 = 6000        # no upper bound on the other crops
ybar3 = 0           # no upper bound on the other crops

# Define the distributions
Ξ = [Uniform(0, 5), Uniform(0, 5), Uniform(10, 30)]

# output
3-element Array{Uniform{Float64},1}:
 Uniform{Float64}(a=0.0, b=5.0)
 Uniform{Float64}(a=0.0, b=5.0)
 Uniform{Float64}(a=10.0, b=30.0)
```
Great now we can formulate and solve the problem using `InfiniteOpt`:
```jldoctest 2-stage; output = false
using InfiniteOpt, JuMP, Ipopt

# Initialize the model
model = InfiniteModel(Ipopt.Optimizer, seed = true) # seed to test output
set_optimizer_attribute(model, "print_level", 0)

# Define the random parameters
@infinite_parameter(model, ξ[c in C] in Ξ[c], num_supports = num_scenarios)

# Define the variables and bounds
@hold_variable(model, 0 <= x[C] <= xbar)
@infinite_variable(model, 0 <= y[C](ξ))
@infinite_variable(model, 0 <= w[C](ξ))

# Define the objective
@objective(model, Min, sum(α[c] * x[c] for c in C) +
           expect(sum(β[c] * y[c] - λ[c] * w[c] for c in C), ξ,
                  use_existing_supports = true))

# Define the constraints
@constraint(model, capacity, sum(x[c] for c in C) <= xbar)
@constraint(model, balance[c in C], ξ[c] * x[c] + y[c] - w[c] >= d[c])
@constraint(model, w[3] <= wbar3)
@constraint(model, y[3] <= ybar3)

# Optimize and get the results
optimize!(model)
x_opt = value.(x)
profit = -objective_value(model)

# Print the results
println("Land Allocations: ", [round(x_opt[k], digits = 2) for k in keys(x_opt)])
println("Expected Profit: \$", round(profit, digits = 2))

# output

Land Allocations: [48.56, 214.77, 236.67]
Expected Profit: $57100.21
```
```
Land Allocations: [48.56, 214.77, 236.67]
Expected Profit: $57100.21
```

We did it! An interesting modification would be to use a ``CVaR`` risk measure
instead of an expectation. This also can be readily achieved via `InfiniteOpt`.
The ``CVaR`` measure is defined:
```math
CVaR_\epsilon(X) = \underset{t \in \mathbb{R}}{\text{inf}}\left\{t + \frac{1}{1-\epsilon} \mathbb{E}[\text{max}(0, X - t)] \right\}
```
where ``\epsilon`` is the confidence level. Inserting this into the formulation,
we now obtain:
```math
\begin{aligned}
	&&\underset{x, y(\xi), w(\xi), t, q(\xi)}{\text{min}} &&& \sum_{c \in C} \alpha_c x_c + t + \frac{1}{1-\epsilon} \mathbb{E}_{\xi}[q(\xi)] \\
	&&\text{s.t.} &&& \sum_{c \in C} x_c \leq \bar{x}\\
	&&&&& \xi_c x_c + y_c(\xi) - w_c(\xi) \geq d_c, && c \in C \\
    &&&&& 0 \leq x_c \leq \bar{x}, && c \in C \\
    &&&&& 0 \leq y_c(\xi) \leq \bar{y}_c, && c \in C \\
    &&&&& 0 \leq w_c(\xi) \leq \bar{w}_c, && c \in C \\
    &&&&& \xi_c \in \Xi_c, && c \in C \\
    &&&&& q(\xi) \geq \sum_{c \in C}\beta_c y_c(\xi) - \lambda_c w_c(\xi) - t \\
    &&&&& q(\xi) \geq 0
\end{aligned}
```
where ``q(\xi)`` is introduced to handle the max operator. Let's update and
resolve our `InfiniteOpt` model using ``\epsilon = 0.95``:
```jldoctest 2-stage; output = false
# Define the additional variables
@hold_variable(model, t)
@infinite_variable(model, q(ξ) >= 0)

# Redefine the objective
@objective(model, Min, sum(α[c] * x[c] for c in C) + t + 1 \ (1 - 0.95) *
           expect(q, ξ, use_existing_supports = true))

# Add the max constraint
@constraint(model, max, q >= sum(β[c] * y[c] - λ[c] * w[c] for c in C) - t)

# Optimize and get the results
optimize!(model)
x_opt = value.(x)
y_opt = value.(y)
w_opt = value.(w)
profit = -sum(α[c] * x_opt[c] for c in C) - 1 / num_scenarios *
         sum(β[c] * y_opt[c][k] - λ[c] * w_opt[c][k] for c in C, k in 1:num_scenarios)

# Print the results
println("Land Allocations: ", [round(x_opt[k], digits = 2) for k in keys(x_opt)])
println("Expected Profit: \$", round(profit, digits = 2))

# output

Land Allocations: [224.58, 83.63, 191.88]
Expected Profit: $47840.01
```
```
Land Allocations: [224.58, 83.63, 191.88]
Expected Profit: $47840.01
```

## Optimal Control
In this case study, we seek to determine an optimal control policy for the
trajectory of a hovercraft that travels to a set of dynamic waypoints while
trying to minimize the thrust input. The corresponding dynamic optimization
problem is expressed:
```math
\begin{aligned}
	&&\underset{x(t), v(t), u(t)}{\text{min}} &&& \int_{t \in T} |u(t)|_2^2 dt  \\
	&&\text{s.t.} &&& v(0) = v0\\
	&&&&& \frac{dx}{dt} = v(t), && t \in T\\
    &&&&& \frac{dv}{dt} = u(t), && t \in T\\
    &&&&& x(t_i) = xw_i, && i \in I
\end{aligned}
```
where ``x(t)`` is the Cartesian position, ``v(t)`` is the velocity, ``u(t)`` is
the thrust input, ``xw_i, \ i \in I,`` are the waypoints, and ``T`` is the time
horizon. This contains two ordinary differential equations which aren't currently
supported by `InfiniteOpt`, so we'll need to reformulate them. This can be
achieved via a simple Euler integration. Thus, we define a set of time points
``T_m \subseteq T`` with an equal time step ``\Delta t`` and apply Euler's
method to obtain:
```math
\begin{aligned}
	&&\underset{x(t), v(t), u(t)}{\text{min}} &&& \int_{t \in T} |u(t)|_2^2 dt  \\
	&&\text{s.t.} &&& v(0) = v0\\
	&&&&& x(t_{j+1}) = v(t_{j}) \Delta t + x(t_{j}), && j \in T_m \setminus T_{mf}\\
    &&&&& v(t_{j+1}) = u(t_{j}) \Delta t + v(t_{j}), && j \in T_m \setminus T_{mf}\\
    &&&&& x(t_i) = xw_i, && i \in I
\end{aligned}
```
where ``T_{mf}`` is the final time point in ``T_m``.

Let's implement this in `InfiniteOpt` by first defining our problem parameters:
```jldoctest hovercraft; output = false
# Initial condition
v0 = [0, 0]

# Time horizon
max_time = 60

# Euler parameters
Δt = 1
time_points = Vector(0:Δt:max_time)

# Waypoints
xw = [1 4 6 1; 1 3 0 1] # positions
tw = [0, 25, 50, 60]    # times

# output
4-element Array{Int64,1}:
  0
 25
 50
 60
```

With our parameters set, let's now construct the `InfiniteModel` and solve it:
```jldoctest hovercraft; output = false
using InfiniteOpt, JuMP, Ipopt

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
@BDconstraint(m, [i = 1:2, j = 1:length(tw)](t == tw[j]), x[i] == xw[i, j])

# Specify the objective
@objective(m, Min, integral(u[1]^2 + u[2]^2, t, use_existing_supports = true))

# Optimize the model
optimize!(m)

# output

```
This thus demonstrates how point variables can be used to enable functionality
that is not built in `InfinitOpt`. Finally, let's extract the solution and plot
it to see the finalized trajectory:
```julia
# Get the results
if has_values(m)
    x_opt = value.(x)
end

# Plot the results
scatter(xw[1,:], xw[2,:], label = "Waypoints")
plot!(x_opt[1,:], x_opt[2,:], label = "Trajectory")
xlabel!("x1")
ylabel!("x2")
```
![answer](./assets/hovercraft.png)

## Optimal Control of a Distillation Column
The following example demonstrates the use of `InfiniteOpt` to model and optimize
the control of a binary distillation in a 30 stage distillation column. In
defining the objective function, our goal will be to minimize the deviation
between the instantaneous distillate purity (``y_{A,1}``) and the desired purity (``\bar{y}``). Along with this, we will add an economic consideration by
minimizing the deviation of the instantaneous reflux ratio (``u``) from the
steady-state reflux ratio (``\bar{u}``). To simplify the analysis, the assumption
of constant molal overflow will be implemented, creating constant flowrates
between stages. Similarly, we will simplify equilibrium constraints within the
column by assuming that the relative volatility (``\alpha_{A,B}``) of the two
species is constant at each stage. With these assumptions, material and
equilibrium constraints in the condenser, reboiler, feed stage, and the stripping
and rectifying sections may be translated into a mathematical model for
optimization as described below.
```math
\begin{aligned}
&&\text{min} &&& \int_{t_{0}}^{t_{f}} (y_{A,1}-\bar{y})^2+(u-\bar{u})^2 dt \\
&&\text{s.t.} &&& \frac{d  x_{A, 0} }{d t} = \frac{1}{A_{c}}V(y_{A, 1}-x_{A, 0})\\
	&&&&&  \frac{d  x_{A, i} }{d t} = \frac{1}{A_{t}}(L_{1}(y_{A, i-1}-x_{A, i})-V(y_{A, i}-y_{A, i+1})), &&  i \in {1, ..., FT-1}\\
	&&&&&  \frac{d  x_{A, FT} }{d t} = \frac{1}{A_{t}}(Fx_{A, f}+L_{1}x_{A,FT-1}-L_{2}x_{A,FT}-V(y_{A, FT}-y_{A, FT+1}))\\
	&&&&&  \frac{d  x_{A, i} }{d t} = \frac{1}{A_{t}}(L_{2}(y_{A, i-1}-x_{A, i})-V(y_{A, i}-y_{A, i+1})), && i \in {FT+1, ..., NT}\\
	&&&&&  \frac{d  x_{A, NT+1} }{d t} = \frac{1}{A_{r}}(L_{2}x_{A,NT}-(F-D)x_{A,NT+1}-Vy_{A,NT+1})\\
	&&&&& x_{A,i}=x_{0_{A,i}}, && i \in {0, ...,NT+1}\\
	&&&&& V=L_{1}+D,\hspace{1cm} L_{2}=L_{1}+F,\hspace{1cm} u=\frac{L_{1}}{D}\\
	&&&&& \alpha_{A,B}=\frac{y_{A}(1-x_{A})}{x_{A}(1-y_{A})}
	\end{aligned}
```
Within these equations, ``x_{A,i}`` and ``y_{A,i}`` are the liquid and vapor mole fractions of component A in the ``i^{th}`` equilibrium stage. ``F`` and ``D`` are
the flowrates of the feed and distillate. ``V``is the vapor flow through the
column.  ``L_{1}`` and ``L_{2}`` are the liquid flows in the rectifying and
stripping sections. Finally, ``A_{c}``, `` A_{t}``, and ``A_{r}`` are defined as
the total molar hold up in the condenser, trays, and reboiler. By choosing the
feed stage (``FT``) to be 17, we now have a mathematical model that we can
optimize. To do so, we’ll use `infiniteopt`!

We will begin by defining our constants and inputing the data within the script.
Since the objective function is quadratic, we will use `Ipopt` as our solver.
```jldoctest hovercraft; output = false
# Optimal Control of 30 Tray Distillation Column

using InfiniteOpt, JuMP, Ipopt, Plots

# Define General Parameters

F = 0.4                      # Feed flow rate (kmol/min)
x_A_Feed = 0.5               # Mole fraction of component A in the feed
NT = 30                      # Number of trays in distillation column
FT = 17                      # Feed tray for the distillation column
alpha_AB = 1.6               # Relative volatility of component A and B (assumed to be constant)
t_initial = 0                # Initial time (min)
t_final = 40                 # Final time (min)
Ac = 0.5                     # Total molar hold up in the condensor (kmol)
At = 0.25                    # Total molar hold up in the trays (kmol)
Ar = 1                       # Total molar hold up in the reboiler (kmol)
D = 0.2                      # Distillate flow rate (kmol/min)
u_bar = 2                    # Objective reflux ratio
y_bar = 0.8958               # Objective purity of distillate
init_cond=ones(32,1)*0.6     # Initial Condition
num=15;
# output

```
Next, we need to initialize the model. For this case, I have defined it as `Control`.
```jldoctest hovercraft; output = false
# Initialize the Model

Control = InfiniteModel(Ipopt.Optimizer)
# output
Control = InfiniteModel(Ipopt.Optimizer)
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```
We need to define our infinite parameter and our infinite variables. For this
problem, our mole fractions, rectifying and stripping liquid flows, vapor flow,
and reflux ratio will be time dependent infinite variables, while time will be
our infinite parameter.
```jldoctest hovercraft; output = false
# Define time as the infinite parameter

@infinite_parameter(Control, t in [0, 40])

# Define infinite variables as functions of time
# Note: 0 = Condensor; 31 = Reboiler

@infinite_variable(Control, 0 <= x_A[i=0:NT+1](t) <= 1, start = 0.5)
@infinite_variable(Control, 0 <= y_A[i=0:NT+1](t) <= 1, start = 0.5)
@infinite_variable(Control, V(t) >= 0, start = 0.3)
@infinite_variable(Control, L1(t) >= 0, start = 0.1)
@infinite_variable(Control, u(t) >= 0, start = 2)
@infinite_variable(Control, L2(t) >= 0, start = 0.5)
# output
L2(t)
```
With our infinite parameter and variables defined, we need to address the issue
of the differential equations within the model. To do this, we will reformulate
the expressions by separating and integrating. This requires us to define initial
and final point variables for each differential equation.
```jldoctest hovercraft; output = false
# Define point variables

@point_variable(Control, x_A[i](t_initial), x_A_0[i = 0:NT+1])
@point_variable(Control, x_A[i](t_final), x_A_F[i = 0:NT+1])
# output
1-dimensional DenseAxisArray{PointVariableRef,1,...} with index sets:
    Dimension 1, 0:31
And data, a 32-element Array{PointVariableRef,1}:
 x_A_F[0]
 x_A_F[1]
 x_A_F[2]
 x_A_F[3]
 x_A_F[4]
 ⋮        
 x_A_F[28]
 x_A_F[29]
 x_A_F[30]
 x_A_F[31]
```
Now that we have defined our point variables, we can define our constraints. We
will evaluate the differential constraints using the integral function. Along
with this, we will provide an initial condition.
```jldoctest hovercraft; output = false
# Set Integral Default for Integral Functions

set_integral_defaults(Control, eval_method = gauss_legendre, num_supports = num)

# Define Constraints

# (20b) Condenser Constraint
@constraint(Control, x_A_F[0] - x_A_0[0] == integral(1/Ac * V * (y_A[1]-x_A[0]), t))
# (20c) Rectifying Section Constraint
@constraint(Control, [i=1:FT-1], x_A_F[i] - x_A_0[i] == integral(1/At * (L1 * (y_A[i-1] - x_A[i]) - V * (y_A[i]-y_A[i+1])), t))
# (20d) Feed Tray Constraints
@constraint(Control, x_A_F[FT] - x_A_0[FT] == integral(1/At * (F * x_A_Feed + L1 * x_A[FT-1] - L2 * x_A[FT] - V * (y_A[FT] - y_A[FT+1])), t))
# (20e) Stripping Section Constraints
@constraint(Control,[i=FT+1:NT], x_A_F[i] - x_A_0[i] == integral(1/At * (L2 * (y_A[i-1] - x_A[i]) - V * (y_A[i] - y_A[i+1])),  t))
# (20f) Reboiler Constraint
@constraint(Control, x_A_F[NT+1] - x_A_0[NT+1] == integral(1/Ar * (L2 * x_A[NT] - (F-D) * x_A[NT+1] - V * y_A[NT+1]),t))
# (20g) Initial Condition
@BDconstraint(Control, initial[i in 0:NT+1](t==0), x_A[i] == init_cond[i+1])
# (20h) Material Balance Constraints
@constraint(Control, V == L1 + D)
@constraint(Control, L2 == L1 + F)
@constraint(Control, u * D == L1)
# (20i) Equilibrium Constraint
@constraint(Control,[i=0:NT+1], alpha_AB * (x_A[i] * (1 - y_A[i])) == y_A[i] * (1 - x_A[i]))

# output
1-dimensional DenseAxisArray{InfiniteConstraintRef{ScalarShape},1,...} with index sets:
    Dimension 1, 0:31
And data, a 32-element Array{InfiniteConstraintRef{ScalarShape},1}:
 -0.6000000000000001 y_A[0](t)*x_A[0](t) + 1.6 x_A[0](t) - y_A[0](t) = 0.0    
 -0.6000000000000001 y_A[1](t)*x_A[1](t) + 1.6 x_A[1](t) - y_A[1](t) = 0.0    
 -0.6000000000000001 y_A[2](t)*x_A[2](t) + 1.6 x_A[2](t) - y_A[2](t) = 0.0    
 -0.6000000000000001 y_A[3](t)*x_A[3](t) + 1.6 x_A[3](t) - y_A[3](t) = 0.0    
 -0.6000000000000001 y_A[4](t)*x_A[4](t) + 1.6 x_A[4](t) - y_A[4](t) = 0.0    
 ⋮                                                                            
 -0.6000000000000001 y_A[28](t)*x_A[28](t) + 1.6 x_A[28](t) - y_A[28](t) = 0.0
 -0.6000000000000001 y_A[29](t)*x_A[29](t) + 1.6 x_A[29](t) - y_A[29](t) = 0.0
 -0.6000000000000001 y_A[30](t)*x_A[30](t) + 1.6 x_A[30](t) - y_A[30](t) = 0.0
 -0.6000000000000001 y_A[31](t)*x_A[31](t) + 1.6 x_A[31](t) - y_A[31](t) = 0.0
```
With our constraints defined, we need to create our objective function.
```jldoctest hovercraft; output = false
# Define the objective function

@objective(Control, Min, integral((y_A[1] - y_bar) ^ 2 + (u - u_bar) ^ 2, t))
# output
integral(y_A[1](t)² + u(t)² - 1.7916 y_A[1](t) - 4 u(t) + 4.80245764)
```
With that, the entire model is specified. All we need to do is optimize!
```jldoctest hovercraft; output = false
# Optimize the model

optimize!(Control)
# output
4
```
