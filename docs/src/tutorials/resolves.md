```@meta
DocTestFilters = [r"≤|<=", r" == | = ", r" ∈ | in ", r" for all | ∀ ", r"d|∂", 
                  r"E|𝔼", r"integral|∫"]
```

# Re-solve Tutorial
There are settings where one may want to repeatedly solve problems defined in `InfiniteOpt`, such as in model predictive control. To support this, `InfiniteOpt` provides a re-solve framework that enables model updates between solves without rebuilding the backend or reinitializing the solver. This is facilitated through the following APIs:

- [`set_parameter_value`](@ref) for updating parameters and parameter functions
- [`warmstart_backend_start_values`](@ref) for warmstarting the backend using a previous solution

As a result, we get a persistent backend and solver instance across subsequent solves, reducing initialization overhead and improving performance.

Below, we demonstrate how to do efficient re-solves via model predictive control of a continuous stirred tank reactor (CSTR). For re-solving on GPUs, see the [InfiniteExaModels](@ref examodels_solves) section.

## Preliminaries 
### Software Setup
First, we need to make sure everything is installed. This will include:
- installing Julia 
- installing `InfiniteOpt.jl`
- installing `DifferentialEquations.jl` to simulate our system via an ODE problem
- installing wanted optimizers e.g., `Ipopt.jl`
See [Installation](@ref) for more information.

### Problem Formulation
Consider a model predictive control (MPC) problem over a time domain ``D_{mpc} = [t_{0}, t_{f}]``, consisting of ``N_{mpc}`` control steps with a control interval of ``\Delta t``. At each control step ``t_k`` for ``k \in \{1, 2, ..., N_{mpc}\}``, an optimal control problem is solved to predict system behavior and determine the optimal input for the next control step. As such, this sets up the stage for optimal control re-solves. See below for a diagram of the MPC workflow:

![mpc](../assets/mpc_diagram.png)

We define our optimal control problem over the prediction horizon ``\mathcal{D_p} = [0, t_p]``, where the objective is set to track a setpoint ``T_{sp}(t)``. We also define the control horizon ``\mathcal{D_c} = [0, t_c]``, where ``t_c \leq t_p``, representing the subset of ``\mathcal{D_p}`` over which the control input is allowed to vary. Here is the full optimal control formulation:
```math
\begin{aligned}
	&&\underset{C_A(t), T(t), T_c(t)}{\text{min}} &&& \int_{t \in \mathcal{D_p}} (T(t) - T_{sp}(t))^2 dt\\
	&&\text{s.t.} &&& C_a(0) = C_{A}(t_k)\\
    &&&&& T(0) = T(t_k)\\
    &&&&& k = k_0 \exp\left(\frac{-E_R}{T(t)}\right), & \forall t \in \mathcal{D_p}\\
    &&&&& r_A = kC_A(t),& \forall t \in \mathcal{D_p}\\
	&&&&& \frac{\partial C_A(t)}{\partial t} = \frac{F(C_f - C_A(t)) - Vr_A}{V}, & \forall t \in \mathcal{D_p}\\
    &&&&& \frac{\partial T(t)}{\partial t} = \frac{F\rho c_p(T_f - T(t)) + r_AH_RV + U_A(T_c(t) - T(t))}{\rho c_p V}, & \forall t \in \mathcal{D_p}\\
    &&&&& T_c(t) = T_c(t_c), & \forall t \in \mathcal{D_p} \setminus \mathcal{D_c}\\
    &&&&& C_A(t) \geq 0\\
    &&&&& 273.15 \leq T(t) \leq 400\\
    &&&&& 250 \leq T_c(t) \leq 350\\
\end{aligned}
```
Here, our state variables are concentration ``C_A(t)`` and reactor temperature ``T(t)``, while the control input is the jacket temperature ``T_c(t)``. There is also a constraint to keep input ``T_c(t)`` constant beyond the control horizon ``D_c``. Overall, this results in a dynamic model based on time ``t``.

## Problem setup
### Parameter specification
Before moving on, we'll need to define the necessary constants and problem 
parameters. We'll define the following in our 
Julia session (these could also be put into a script as shown later on):
```jldoctest quick
julia> Δt = 0.1; t0 = 0; tf = 2; tp = 3; tc = 2.5; # set MPC simulation parameters

julia> Dmpc = t0:Δt:tf; Dp = 0:Δt:tp; # set problem domains

julia> states = [0.9, 305]; # set the initial conditions for state

julia> F, V, rho, cp, Hr, E, k₀, Ua, cf, Tf = [100, 100, 1000, 0.239, 5e4, 8750, 7.2e10, 5e4, 1.0, 350]; # set the problem parameters
```

We'll also create a function for our temperature setpoint ``T_{sp}(t)``, which is defined below:
```math
T_{sp}(t) = 
\left\{
\begin{aligned}
&& 310, & \quad t_0 \leq t < 0.7 \\
&& 323, & \quad 0.7 \leq t < 1.3 \\
&& 318, & \quad 1.3 \leq t \leq t_f \\
\end{aligned}
\right.
```
```jldoctest quick
julia> function setpoint(t, offset)
        t += offset
        if t < 0.7
            return 310
        elseif t < 1.3
            return 323
        else
            return 318
        end
        end;
```
Note that an `offset` argument is added to ensure the setpoint function returns accurate values depending on what control step ``t_k`` the problem is posed at. 

### Problem model 
From here, we'll initialize our `InfiniteModel`, using Ipopt as the optimizer:
```jldoctest quick
julia> using InfiniteOpt, Ipopt;

julia> model = InfiniteModel(TranscriptionBackend(Ipopt.Optimizer, update_parameter_functions = true))
An InfiniteOpt Model
Feasibility problem with:
  Finite parameters: 0
  Infinite parameters: 0
  Variables: 0
  Derivatives: 0
  Measures: 0
Transformation backend information:
  Backend type: TranscriptionBackend
  Solver: Ipopt
  Transformation built and up-to-date: false
```
Note that we also set the keyword argument `update_parameter_functions = true`. This is only necessary if you plan on updating parameter functions while working with `TranscriptionBackend`s. In this case, we will need it to update our setpoint function ``T_{sp}(t)``. Learn more about `TranscriptionBackend`s on our 
[Model Transcription](@ref transcription_docs) page.

Continuing on, we define our InfiniteOpt problem as normal:
```jldoctest quick
julia> @finite_parameter(model, Ca0 == states[1]);

julia> @finite_parameter(model, T0 == states[2]);

julia> @infinite_parameter(model, t in [0, tp], supports = collect(Dp), 
                           derivative_method = OrthogonalCollocation(3));

julia> @variable(model, 0 ≤ Ca, Infinite(t), start = states[1]);

julia> @variable(model, 273.15 ≤ T ≤ 400, Infinite(t), start = states[2]);

julia> @variable(model, 250 ≤ Tc ≤ 350, Infinite(t), start = 300);

julia> @parameter_function(model, Tsp == t -> setpoint(t, tk))
Tsp(t)

julia> @objective(model, Min, integral((T - Tsp)^2, t))
∫{t ∈ [0, 3]}[T(t)² - 2 Tsp(t)*T(t) + Tsp(t)²]

julia> @constraint(model, Ca(0) == Ca0)
Ca(0) - Ca0 = 0

julia> @constraint(model, T(0) == T0)
T(0) - T0 = 0

julia> @expression(model, k, k₀ * exp(-E/T))
7.2e10 * exp(-8750.0 / T(t))

julia> @expression(model, rate, k*Ca)
7.2e10 * Ca(t) * exp(-8750.0 / T(t))

julia> @constraint(model, deriv(Ca, t) == (F*(cf - Ca) - V*rate)/V);

julia> @constraint(model, deriv(T, t) == (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T)));
```

To add a constraint for the control horizon, we can use [`DomainRestriction`](@ref).
 ```jldoctest quick
julia> control_func(t_c) = (tc ≤ t_c ≤ tp); # Function for domain values

julia> control_domain = DomainRestriction(control_func, t);   # Create domain restriction for infinite parameter t

julia> @constraint(model, Tc == Tc(tc), control_domain);  # Add constraint defined on restricted domain
```

We also call [`constant_over_collocation`](@ref constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)) 
on our control variable `Tc` to address any unwanted degrees of freedom introduced by internal collocation 
nodes from [`OrthogonalCollocation`](@ref):
```jldoctest quick
julia> constant_over_collocation.(Tc, t);
``` 
That's it for defining our `InfiniteOpt` problem! Now, let's solve via [`optimize!`](@ref):
```jldoctest quick; setup = :(set_attribute(model, "print_level", 0))
julia> tk = 0;   # Offset for first control step

julia> optimize!(model)
```
```
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:      969
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:      427

Total number of variables............................:      305
                     variables with only lower bounds:       61
                variables with lower and upper bounds:      122
                     variables with only upper bounds:        0
Total number of equality constraints.................:      274
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  5.0280000e+02 3.92e+01 6.25e-01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  2.7236558e+01 5.22e+00 5.37e+00  -1.0 3.92e+01    -  7.04e-01 1.00e+00f  1
   2  9.0788361e+00 8.53e-01 5.57e-01  -1.0 2.14e+01    -  8.38e-01 1.00e+00f  1
   3  4.2566639e+00 3.38e-01 1.78e-01  -1.0 3.33e+01    -  9.02e-01 1.00e+00f  1
   4  3.1557133e+00 7.86e-02 2.07e-02  -1.0 2.14e+01    -  1.00e+00 1.00e+00h  1
   5  2.8579262e+00 2.42e-02 2.30e-03  -1.7 1.46e+01    -  1.00e+00 1.00e+00h  1
   6  2.7627030e+00 1.23e-02 1.58e-04  -2.5 8.72e+00    -  1.00e+00 1.00e+00h  1
   7  2.7341074e+00 3.13e-03 4.70e-05  -3.8 4.18e+00    -  1.00e+00 1.00e+00h  1
   8  2.7262156e+00 4.54e-04 8.51e-06  -3.8 1.55e+00    -  1.00e+00 1.00e+00h  1
   9  2.7245921e+00 2.43e-05 4.28e-07  -5.7 3.59e-01    -  1.00e+00 1.00e+00h  1
iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
  10  2.7245135e+00 6.10e-08 1.14e-09  -5.7 1.78e-02    -  1.00e+00 1.00e+00h  1
  11  2.7245096e+00 1.38e-10 1.74e-12  -8.6 8.88e-04    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 11

                                   (scaled)                 (unscaled)
Objective...............:   2.7245096440910856e+00    2.7245096440910856e+00
Dual infeasibility......:   1.7412247422617024e-12    1.7412247422617024e-12
Constraint violation....:   1.3787371244688984e-10    1.3787371244688984e-10
Variable bound violation:   2.8948174417564587e-06    2.8948174417564587e-06
Complementarity.........:   2.7463356154477961e-09    2.7463356154477961e-09
Overall NLP error.......:   2.7463356154477961e-09    2.7463356154477961e-09


Number of objective function evaluations             = 12
Number of objective gradient evaluations             = 12
Number of equality constraint evaluations            = 12
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 12
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 11
Total seconds in IPOPT                               = 0.009

EXIT: Optimal Solution Found.
```
As shown in the solver output, we've converged in 11 iterations!

### Simulating the system
Now we'll extract the input for the next control step via
[`value`](@ref JuMP.value(::GeneralVariableRef)), which will be used to simulate the system forward:
```jldoctest quick
julia> Tc_opt = value.(Tc)[2];
```
To simulate our system, we choose to set up an ODE problem via `DifferentialEquations.jl`. Start by creating a function for the CSTR dynamics:
```jldoctest quick
julia> using DifferentialEquations

julia> function cstr_dynamics!(dx, x, p, t)
        Ca, T, Tc = x
        k = k₀ * exp(-E/T)
        rate = k * Ca
        dx[1] = F/V * (cf - Ca) - V * rate
        dx[2] = (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T))
        end;
```
Then create a function that solves the ODE problem and returns the updated state:
```jldoctest quick
julia> function cstr_sim(x0, tspan)
        prob = ODEProblem(cstr_dynamics!, x0, tspan)
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)
        xsol = sol.u[end]
        return xsol
        end;
```
Now we simulate the system forward to get our updated states. These will become the initial conditions for the next optimal control solve.
```jldoctest quick
julia> tspan = (tk, tk + Δt); # Timespan to integrate over

julia> xk = [states..., Tc_opt];  # Current state + input

julia> xsol = cstr_sim(xk, tspan);    # Integrate the ODE problem

julia> states[1:2] = xsol[1:2];       # Obtain updated states
```

## Re-solves
### Updating the model
Now we need to update our model for the next solve. First, we update our control step accordingly:
```jldoctest quick
julia> tk += Δt;
```
We then update the start values via [`warmstart_backend_start_values`](@ref). This will use the previous solution to warmstart all possible variables (e.g., primal, dual, etc...).
```jldoctest quick
julia> warmstart_backend_start_values(model);
```
!!! note
    This will only update the start values in the backend.
    The start values in the `InfiniteModel` will remain the same.

We'll also update the initial condition parameters using [`set_parameter_value`](@ref).
```jldoctest quick
julia> set_parameter_value(Ca0, states[1]);

julia> set_parameter_value(T0, states[2]);
```
Lastly, we'll need to update the setpoint parameter function, which can also be done via [`set_parameter_value`](@ref). We'll create a new function `Tsp_new` with an updated offset value for the new control step. Then we can update `Tsp` with `Tsp_new`.
```jldoctest quick
julia> Tsp_new = (t) -> setpoint(t, tk);

julia> set_parameter_value(Tsp, Tsp_new);
```
!!! note
    The new function must have the same infinite parameter format as the original parameter function.

Now we can call `optimize!` again to solve our updated model! 
```jldoctest quick
julia> optimize!(model)
```
```
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:      989
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:      427

Total number of variables............................:      305
                     variables with only lower bounds:       61
                variables with lower and upper bounds:      122
                     variables with only upper bounds:        0
Total number of equality constraints.................:      285
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  5.0280000e+02 3.92e+01 1.80e+00  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  3.3245030e+01 4.66e+00 1.79e+00  -1.0 3.92e+01    -  7.14e-01 1.00e+00f  1
   2  9.1455568e+00 1.05e+00 7.82e-01  -1.0 2.24e+01    -  8.22e-01 1.00e+00f  1
   3  4.2533979e+00 3.66e-01 2.05e-01  -1.0 3.32e+01    -  9.02e-01 1.00e+00f  1
   4  3.1557492e+00 7.79e-02 2.06e-02  -1.0 2.13e+01    -  1.00e+00 1.00e+00h  1
   5  2.8576232e+00 2.40e-02 2.29e-03  -1.7 1.46e+01    -  1.00e+00 1.00e+00h  1
   6  2.7626143e+00 1.22e-02 1.58e-04  -2.5 8.71e+00    -  1.00e+00 1.00e+00h  1
   7  2.7340812e+00 3.12e-03 4.68e-05  -3.8 4.17e+00    -  1.00e+00 1.00e+00h  1
   8  2.7262098e+00 4.52e-04 8.47e-06  -3.8 1.54e+00    -  1.00e+00 1.00e+00h  1
   9  2.7245920e+00 2.41e-05 4.25e-07  -5.7 3.58e-01    -  1.00e+00 1.00e+00h  1
iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
  10  2.7245139e+00 6.02e-08 1.13e-09  -5.7 1.77e-02    -  1.00e+00 1.00e+00h  1
  11  2.7245100e+00 1.38e-10 1.74e-12  -8.6 8.87e-04    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 11

                                   (scaled)                 (unscaled)
Objective...............:   2.7245100343448030e+00    2.7245100343448030e+00
Dual infeasibility......:   1.7391275767757273e-12    1.7391275767757273e-12
Constraint violation....:   1.3756107364315540e-10    1.3756107364315540e-10
Variable bound violation:   2.8948873023182387e-06    2.8948873023182387e-06
Complementarity.........:   2.7460183036212907e-09    2.7460183036212907e-09
Overall NLP error.......:   2.7460183036212907e-09    2.7460183036212907e-09


Number of objective function evaluations             = 12
Number of objective gradient evaluations             = 12
Number of equality constraint evaluations            = 12
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 12
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 11
Total seconds in IPOPT                               = 0.009

EXIT: Optimal Solution Found.
```
!!! warning
    The re-solve framework's efficiency is based on the idea that the problem structure remains the same between solves (AKA no new variables or constraints are added).
    If changing the structure or using APIs other than `warmstart_backend_start_values` or `set_parameter_value` for updating, then a new backend and solver instance will be made from scratch.

### [Performance tips](@id resolve_performance)
Given a warmstart, we can also reduce the number of iterations by decreasing certain solver options via [`set_optimizer_attribute`](https://jump.dev/JuMP.jl/stable/api/JuMP/#set_optimizer_attribute). Here, `bound_push` and `bound_frac` are solver options that determine how much the initial point might have to be adjusted to be sufficiently inside the problem bounds. Since our warmstart is presumably close to the actual solution, we can decrease both values for better performance. We can also decrease the initial barrier parameter `mu_init` for the same reason.
```jldoctest quick
julia> set_optimizer_attribute(model, "bound_push", 1e-8); # Desired minimum distance from intial point to bounds

julia> set_optimizer_attribute(model, "bound_frac", 1e-8); # Desired minimum relative distance from initial point to bounds

julia> set_optimizer_attribute(model, "mu_init", 1e-11); # Initial barrier parameter value
```
Note that the chosen values here are not always good, as optimal values are problem dependent. Moving on, solving with these new solver options reduces the number of iterations to 8:
```jldoctest quick
julia> optimize!(model)
```
```
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:      989
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:      427

Total number of variables............................:      305
                     variables with only lower bounds:       61
                variables with lower and upper bounds:      122
                     variables with only upper bounds:        0
Total number of equality constraints.................:      285
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  1.7342046e+01 5.78e+00 1.45e+00 -11.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  1.5717415e+01 2.93e-01 6.72e-01 -11.0 2.10e+01    -  8.46e-01 1.00e+00f  1
   2  3.8283977e+00 2.28e+00 5.69e-01 -11.0 9.30e+01    -  3.94e-01 1.00e+00f  1
   3  3.8665197e+00 6.46e-03 5.07e-02 -11.0 6.13e+00    -  6.15e-01 1.00e+00h  1
   4  3.6493980e+00 4.97e-03 9.58e-02 -11.0 4.58e+00    -  9.45e-07 1.00e+00f  1
   5  1.9537673e+00 9.30e-01 4.58e-02 -11.0 7.46e+01    -  6.12e-02 8.36e-01f  1
   6  1.9164754e+00 2.85e-02 2.87e-03 -11.0 8.29e+00    -  9.28e-01 1.00e+00h  1
   7  1.9165162e+00 2.77e-06 2.29e-06 -11.0 1.21e-01    -  9.99e-01 1.00e+00h  1
   8  1.9165162e+00 3.04e-12 8.34e-14 -11.0 1.24e-04    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 8

                                   (scaled)                 (unscaled)
Objective...............:   1.9165161630626244e+00    1.9165161630626244e+00
Dual infeasibility......:   8.3353452221015503e-14    8.3353452221015503e-14
Constraint violation....:   3.0420110874729289e-12    3.0420110874729289e-12
Variable bound violation:   3.4981688941115863e-06    3.4981688941115863e-06
Complementarity.........:   1.3620439999832149e-11    1.3620439999832149e-11
Overall NLP error.......:   1.3620439999832149e-11    1.3620439999832149e-11


Number of objective function evaluations             = 9
Number of objective gradient evaluations             = 9
Number of equality constraint evaluations            = 9
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 9
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 8
Total seconds in IPOPT                               = 0.005

EXIT: Optimal Solution Found.
```

## Model Predictive Control Script 
The steps outlined in the sections above can be captured in an MPC loop. This is summarized in the script below:
```julia
using InfiniteOpt, Ipopt, DifferentialEquations

# DEFINE THE PROBLEM CONSTANTS
F = 100       # m³/s
V = 100       # m³
rho = 1000    # kg/m³
cp = 0.239    # J/kg K
Hr = 5e4      # Heat of reaction J/mol
E = 8750      # E/R K
k₀ = 7.2e10   # Pre-exponential factor 1/s
Ua = 5e4      # Heat transfer coefficient J/s K
cf = 1.0      # Feed concentration mol/m³
Tf = 350      # Feed temperature K
states = [0.9, 305]  # Initial conditions

# DEFINE MPC PARAMETERS
t0 = 0            # Initial simulation time
tf = 2            # Final simulation time
Δt = 0.1          # Control interval
tp = 3            # Prediction endpoint
tc = 2.5          # Control endpoint
Dmpc = t0:Δt:tf   # MPC Simulation domain
Dp = 0:Δt:tp      # Prediction horizon
Dc = 0:Δt:tc      # Control horizon

# INITIALIZE RELEVANT FUNCTIONS
function setpoint(t, offset)
    t += offset
    if t < 0.7
        return 310
    elseif t < 1.3
        return 323
    else
        return 318
    end
end

function cstr_dynamics!(dx, x, p, t)
    Ca, T, Tc = x
    k = k₀ * exp(-E/T)
    rate = k * Ca
    dx[1] = F/V * (cf - Ca) - V * rate
    dx[2] = (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T))
end

function cstr_sim(x0, tspan)
    prob = ODEProblem(cstr_dynamics!, x0, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)
    xsol = sol.u[end]
    return xsol
end

# INITIALIZE THE MODEL
model = InfiniteModel(TranscriptionBackend(Ipopt.Optimizer, update_parameter_functions = true))

# INITIALIZE PARAMETERS & PARAMETER FUNCTIONS
@infinite_parameter(model, t ∈ [0, tp],
                supports = collect(Dp),
                derivative_method = OrthogonalCollocation(3))
@finite_parameter(model, Ca0 == states[1])
@finite_parameter(model, T0 == states[2])
@parameter_function(model, Tsp == t -> setpoint(t, 0))

# INITIALIZE THE VARIABLES
@variables(model, begin
    0 ≤ Ca, Infinite(t), (start = states[1])
    273.15 ≤ T ≤ 400, Infinite(t), (start = states[2])
    250 ≤ Tc ≤ 350, Infinite(t), (start = 300)
end)

# SET THE OBJECTIVE
@objective(model, Min, ∫((T - Tsp)^2, t))

# SET THE INITIAL CONDITIONS
@constraint(model, initialCa, Ca(0) == Ca0)
@constraint(model, initialT, T(0) == T0)

# SET THE PROBLEM CONSTRAINTS
@expression(model, k, k₀ * exp(-E/T))
@expression(model, rate, k * Ca)
@constraint(model, ∂(Ca, t) == (F * (cf - Ca) - V * rate)/V)
@constraint(model, ∂(T, t) == (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T)))

# SET THE CONTROL HORIZON CONSTRAINT
control_func(t_c) = (tc ≤ t_c ≤ tp)
control_domain = DomainRestriction(control_func, t)
@constraint(model, Tc == Tc(tc), control_domain)

# ADJUST DEGREES OF FREEDOM FOR CONTROL VARIABLES
constant_over_collocation.(Tc, t)

# MPC LOOP
for tk in Dmpc
    # SOLVE THE MODEL
    optimize!(model)

    # GET THE OPTIMAL INPUT
    Tc_opt = value.(Tc)[2]
    
    # SIMULATE SYSTEM FORWARD
    tspan = (tk, tk + Δt)
    xk = [states..., Tc_opt]
    xsol = cstr_sim(xk, tspan)
    states[1:2] = xsol[1:2]

    # WARMSTART MODEL FOR NEXT SOLVE
    warmstart_backend_start_values(model)

    # UPDATE PARAMETERS
    set_parameter_value(Ca0, states[1])
    set_parameter_value(T0, states[2])
    Tsp_new = (t) -> setpoint(t, tk + Δt)
    set_parameter_value(Tsp, Tsp_new)

    # ADJUST SOLVER PARAMETERS
    set_optimizer_attribute(model, "bound_push", 1e-8)
    set_optimizer_attribute(model, "bound_frac", 1e-8)
    set_optimizer_attribute(model, "mu_init", 1e-11)
end
```

## [Re-solves with InfiniteExaModels.jl](@id examodels_solves)
### Software Setup
Alternatively, [`InfiniteExaModels.jl`](https://github.com/infiniteopt/InfiniteExaModels.jl) can be used to improve performance. Instead of a `JuMP` model, we can work with an `ExaModel` that exploits recurrent problem structure, cutting down on both model and automatic differentiation (AD) time. First, we'll need to ensure the following packages are installed:
- `InfiniteExaModels.jl`
- `NLPModelsIpopt.jl` (if solving on CPU)

If solving on GPU, we'll also need to install:
- `MadNLP.jl`
- `CUDA.jl`
- `CUDSS.jl`

!!! note
    Currently, the GPU workflow is only available on NVIDIA GPUs that support CUDA.

### Solving on CPU
When initializing the model, we specify an `ExaTranscriptionBackend` which will transcribe the `InfiniteOpt` problem into an `ExaModel`. For solving on CPU specifically, we use `NLPModelsIpopt` as the solver.

```julia
julia> using InfiniteExaModels, NLPModelsIpopt

julia> model = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver));
```

### Solving on GPU
For GPU, we'll want to define an `ExaTranscriptionBackend` that employs a `CUDABackend`. We'll also use the `MadNLPGPU` module from `MadNLP.jl` for solving.
```julia
julia> using InfiniteExaModels, MadNLPGPU, CUDA

julia> model = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver, backend = CUDABackend()));
```
Since we'll be indexing `Tc` for the optimal control input, we'll also need to set `allowscalar` to true.
```julia
julia> CUDA.allowscalar(true);
```
As mentioned in the [`Performance tips`](@ref resolve_performance) section earlier, we can be more performant by adjusting `MadNLP` solver options when given a warmstart. Note that the numerical value chosen here is different than what we originally chose when warmstarting with `Ipopt`. This highlights that the optimal value is not only problem-dependent, but also depends on the chosen solver.
```
julia> set_optimizer_attribute(model, "mu_init", 2E-2)
```
That's all we need to change for the GPU setup! From here, we follow the same steps as the sections above to define the rest of our problem & solve.

!!! note
    Although `MadNLP.jl` does support `bound_push` and `bound_fac` as solver options, they currently do not have an effect in re-solves.
