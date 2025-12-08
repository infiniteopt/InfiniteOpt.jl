```@meta
DocTestFilters = [r"≤|<=", r" == | = ", r" ∈ | in ", r" for all | ∀ ", r"d|∂", 
                  r"E|𝔼", r"integral|∫"]
```

# Resolve Guide
Below we show how to do efficient resolves of problems defined in `InfiniteOpt`, such as in model predictive control. Please refer to the 
Guide on our subsequent pages for more complete information on modeling a stand-alone problem.

## Preliminaries 
### Software Setup
First, we need to make sure everything is installed. This will include:
- installing Julia 
- installing the `InfiniteOpt.jl`, `JuMP.jl`, and `DifferentialEquations.jl` packages
- installing wanted optimizers e.g., `Ipopt.jl`
See [Installation](@ref) for more information.

### Optimal Control Formulation
Now we need to formulate the problem we want to solve mathematically. For example, 
let's define a simple optimal control model for a continuous stirred tank reactor (CSTR):
```math
\begin{aligned}
	&&\underset{C_A(t), T(t), T_c(t)}{\text{min}} &&& \int_{t \in \mathcal{D_t}} (T(t) - T_{sp}(t))^2 dt\\
	&&\text{s.t.} &&& C_a(0) = C_{A, 0},\\
    &&&&& T(0) = T_{0},\\
    &&&&& k = k_0 \exp\left(\frac{-E_R}{T(t)}\right),\\
    &&&&& r_A = kC_A(t),\\
	&&&&& \frac{\partial C_A(t)}{\partial t} = \frac{F(C_f - C_A(t)) - Vr_A}{V},\\
    &&&&& \frac{\partial T(t)}{\partial t} = \frac{F\rho c_p(T_f - T(t)) + r_AH_RV + U_A(T_c(t) - T(t))}{\rho c_p V},\\
    &&&&& C_A(t) \geq 0,\\
    &&&&& 273.15 \leq T(t) \leq 400,\\
    &&&&& 250 \leq T_c(t) \leq 350,\\
\end{aligned}
```
Here, our state variables are concentration ``C_A(t)`` and reactor temperature ``T(t)``, while the control input is the jacket temperature ``T_c(t)``. This results in a dynamic model based on time ``t`` over the prediction horizon ``D_t``.

## Parameter Specification
Before moving on, we'll need to define the necessary constants and problem 
parameters. We'll define the following in our 
Julia session (these could also be put into a script as shown later on):
```jldoctest quick
julia> Δt = 0.1; Dt = 3; # set time domain parameters 

julia> init = [0.9, 305, 300]; # set the initial conditions

julia> p = [100, 100, 1000, 0.239, 5e4, 8750, 7.2e10, 5e4, 1.0, 350]; # set the problem parameters

julia> F, V, rho, cp, Hr, E, k₀, Ua, cf, Tf = p; # assign variables to each parameter
```

We'll also define a function for the temperature setpoint:
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
        end
setpoint (generic function with 1 method)
```
Note that an `offset` argument is added to ensure the setpoint function returns different values depending on what control step the problem is posed at. 

## Optimal Control Problem
### Model Initialization 
First, we'll initialize our `InfiniteModel` and assign an 
appropriate optimizer that will be used to solve its transcripted variant. Here, let's choose to use Ipopt:
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
Note that `update_parameter_functions` is set to true to enable setpoint function updates for ``T_{sp}(t)`` later on. Learn more about `InfiniteModel`s and optimizers on our 
[Infinite Models](@ref infinite_model_docs) page.

Before moving on, let's make finite parameters via [`@finite_parameter`](@ref) 
to represent initial conditions for ``C_A`` and ``T`` since we'll want to update them between optimal control solves: 
```jldoctest quick
julia> @finite_parameter(model, Ca0 == init[1])
Ca0

julia> @finite_parameter(model, T0 == init[2])
T0
```
Learn more about finite parameters on our [Finite Parameters](@ref finite_param_docs) 
page.

### Infinite Parameters 
The next thing we need to do is identify the infinite domains our problem contains 
and define an infinite parameter(s) for each one via [`@infinite_parameter`]. For 
this problem, we have the time domain ``t \in \mathcal{D}_t``:
```jldoctest quick
julia> @infinite_parameter(model, t in [0, Dt], supports = collect(0:Δt:Dt), 
                           derivative_method = OrthogonalCollocation(3))
t
```
We specify the domain the parameter depends on via `in`.
Here we also directly provide the supports
that will be used to reformulate and solve the problem (i.e., discretize). 
We also specify the derivative evaluation method associated with `t` that will be 
used to evaluate the derivatives numerically. See more information about parameters 
on our [Infinite Parameters] (@ref inf_par_docs) page. Also learn more about 
derivative methods on our [Derivative Operators](@ref deriv_docs) page.

### Variables 
Now that we have an `InfiniteModel` with infinite parameters, let's define our 
decision variables. First, infinite variables (ones that depend on infinite 
parameters) are defined via 
[`@variable`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@variable) 
with the addition of the [`Infinite`](@ref) variable type argument to specify the 
infinite parameters it depends on:
 ```jldoctest quick
julia> @variable(model, 0 ≤ Ca, Infinite(t), start = init[1])
Ca(t)

julia> @variable(model, 273.15 ≤ T ≤ 400, Infinite(t), start = init[2])
T(t)

julia> @variable(model, 250 ≤ Tc ≤ 350, Infinite(t), start = init[3])
Tc(t)
```
Notice that we specify the initial guess for all of them via `start`, which is the same as the initial conditions in this case. We also 
can symbolically define variable bounds accordingly.

We'll also define the setpoint parameter function via [`@parameter_function`](@ref):
```jldoctest quick
julia> @parameter_function(model, Tsp == t -> setpoint(t, tk))
Tsp(t)
```
For more information, please see the [Parameter Function](@ref par_func_docs) page.

That's it for this example, but other problems might also employ the following:
- Finite variables: variables that do not depend on infinite parameters 
  (defined using `@variable`)
- Semi-infinite variables: infinite variables where 1 or more parameters are 
  set a particular point (defined via [Restricted Variables](@ref)).
- Point variables: infinite variables at a particular point (defined via [Restricted Variables](@ref)).

### Objective & Constraints 
Now that the variables and parameters are ready to go, let's define our problem. 
We can define the objective using 
[`@objective`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@objective):
 ```jldoctest quick
julia> @objective(model, Min, integral((T - Tsp)^2, t))
∫{t ∈ [0, 3]}[T(t)² - 2 Tsp(t)*T(t) + Tsp(t)²]
```
Here, we employ [`integral`](@ref) to define the integral. Note that 
objectives must evaluate over all included infinite domains. 

Now let's define the initial conditions using 
[`@constraint`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@constraint) 
in combination with [Restricted Variables](@ref) which will restrict the domain 
of the variables to only be enforced at the initial time:
 ```jldoctest quick
julia> @constraint(model, Ca(0) == Ca0)
Ca(0) - Ca0 = 0

julia> @constraint(model, T(0) == T0)
T(0) - T0 = 0
```
Note that it is important that we include appropriate boundary conditions when using 
derivatives in our model. For more information please see 
[Derivative Operators](@ref deriv_docs).

Next, we add our model constraints that have derivatives using 
[`@constraint`](https://jump.dev/JuMP.jl/v1/api/JuMP/#JuMP.@constraint) 
and [`deriv`](@ref):
 ```jldoctest quick
julia> @expression(model, k, k₀ * exp(-E/T))
7.2e10 * exp(-8750.0 / T(t))

julia> @expression(model, rate, k*Ca)
7.2e10 * Ca(t) * exp(-8750.0 / T(t))

julia> @constraint(model, deriv(Ca, t) == (F*(cf - Ca) - V*rate)/V);

julia> @constraint(model, deriv(T, t) == (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T)));
```
Finally, to address any unwanted degrees of freedom introduced by internal collocation 
nodes with [`OrthogonalCollocation`](@ref). We should call [`constant_over_collocation`](@ref constant_over_collocation(::InfiniteVariableRef, ::GeneralVariableRef)) 
on any control variables:
```jldoctest quick
julia> constant_over_collocation.(Tc, t);
``` 
That's it, now we have our problem defined in `InfiniteOpt`!

## Solution & Queries
### Optimize 
Now that our model is defined, let's optimize it via [`optimize!`](@ref):
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

### Simulate the system
Now that we've solved the optimal control problem, we'll extract the input for the next control step. This can be done via
[`value`](@ref JuMP.value(::GeneralVariableRef)):
```jldoctest quick
julia> Tc_opt = value.(Tc)[2];
```
Next, the optimal input is used to simulate the system forward. Here, we choose to set up an ODE problem via `DifferentialEquations.jl`. Start by creating a function for the CSTR dynamics:
```jldoctest quick
julia> using DifferentialEquations

julia> function cstrDynamics!(dx, x, p, t)
        Ca, T, Tc = x
        k = k₀ * exp(-E/T)
        rate = k * Ca
        dx[1] = F/V * (cf - Ca) - V * rate
        dx[2] = (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T))
        end
```
Then create a function that solves the ODE problem and returns the updated state:
```jldoctest quick
julia> function cstrSim(x0, tspan)
        prob = ODEProblem(cstrDynamics!, x0, tspan)
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)
        xsol = sol.u[end]
        return xsol
        end
```
Now let's simulate the system forward!
```jldoctest quick
julia> tspan = (tk, tk + Δt); # Timespan to integrate over

julia> xk = [init[1], init[2], Tc_opt];  # Current state + input

julia> xsol = cstrSim(xk, tspan);    # Integrate the ODE problem

julia> Ca_opt, T_opt = xsol[1:2];    # Obtain updated states
```
Note that the updated states `Ca_opt` and `T_opt` will become the initial conditions for the next optimal control solve.

## Resolves
### Updating the model
Now we need to update our model for the next solve. First, we update some parameters accordingly:
```jldoctest quick
julia> tk += Δt;     # Updated offset for the next control step
```
We'll then update the initial guesses in the model using the previous solution via `warmstart_backend_start_values`.
```jldoctest quick
julia> warmstart_backend_start_values(model);
```
We'll also update the initial condition parameters using `set_parameter_value`.
```jldoctest quick
julia> set_parameter_value(Ca0, Ca_opt);

julia> set_parameter_value(T0, T_opt);
```
Lastly, we'll need to update the setpoint parameter function, which can also be done via `set_parameter_value`. In this case, we'll create a new function with an updated offset value.
```jldoctest quick
julia> newTsp = (t) -> setpoint(t, tk);

julia> set_parameter_value(Tsp, newTsp);
```
Now we can solve our updated model! 
```jldoctest quick
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
   0  1.7342045e+01 5.78e+00 5.16e-02  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  1.5692725e+01 2.93e-01 5.98e-01  -1.0 2.10e+01    -  8.40e-01 1.00e+00f  1
   2  8.5431043e+00 3.45e-01 1.02e-01  -1.0 2.80e+01    -  6.84e-01 1.00e+00f  1
   3  4.8015451e+00 2.55e-01 5.17e-02  -1.0 3.46e+01    -  5.42e-01 1.00e+00f  1
   4  2.7068147e+00 3.71e-01 8.82e-02  -1.0 4.90e+01    -  9.70e-01 1.00e+00f  1
   5  2.1446416e+00 3.26e-02 5.37e-03  -1.7 1.83e+01    -  9.79e-01 1.00e+00h  1
   6  1.9800640e+00 1.98e-02 3.83e-04  -2.5 1.15e+01    -  1.00e+00 1.00e+00h  1
   7  1.9351403e+00 4.17e-03 9.26e-05  -2.5 5.15e+00    -  1.00e+00 1.00e+00h  1
   8  1.9200825e+00 8.60e-04 1.66e-05  -3.8 2.35e+00    -  1.00e+00 1.00e+00h  1
   9  1.9169920e+00 5.12e-05 1.22e-06  -3.8 5.59e-01    -  1.00e+00 1.00e+00h  1
iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
  10  1.9165133e+00 1.25e-06 2.11e-08  -5.7 9.03e-02    -  1.00e+00 1.00e+00h  1
  11  1.9165046e+00 4.24e-10 8.01e-12  -8.6 1.65e-03    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 11

                                   (scaled)                 (unscaled)
Objective...............:   1.9165045692830063e+00    1.9165045692830063e+00
Dual infeasibility......:   8.0068109934202508e-12    8.0068109934202508e-12
Constraint violation....:   4.2356873564131092e-10    4.2356873564131092e-10
Variable bound violation:   2.8904434543619573e-06    2.8904434543619573e-06
Complementarity.........:   3.3290100135960107e-09    3.3290100135960107e-09
Overall NLP error.......:   3.3290100135960107e-09    3.3290100135960107e-09


Number of objective function evaluations             = 12
Number of objective gradient evaluations             = 12
Number of equality constraint evaluations            = 12
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 12
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 11
Total seconds in IPOPT                               = 0.008

EXIT: Optimal Solution Found.
```
!!! note
    This tutorial assumes that the problem structure remains the same between consecutive solves (AKA no new variables or constraints are added to the problem).

### Performance tips
Given a warmstart, we can also reduce the number of iterations by decreasing certain solver parameters via `set_optimizer_attribute`.
```jldoctest quick
julia> set_optimizer_attribute(model, "bound_push", 1e-8); # Desired minimum distance from intial point to bounds

julia> set_optimizer_attribute(model, "bound_frac", 1e-8); # Desired minimum relative distance from initial point to bounds

julia> set_optimizer_attribute(model, "mu_init", 1e-11); # Initial barrier parameter value
```
Solving with these new solver parameters reduces the number of iterations to 7:
```jldoctest quick
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
   0  1.7342045e+01 5.78e+00 5.16e-02 -11.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  1.5715639e+01 2.93e-01 8.03e-01 -11.0 2.10e+01    -  8.46e-01 1.00e+00f  1
   2  3.8284002e+00 2.28e+00 5.74e-01 -11.0 9.30e+01    -  3.94e-01 1.00e+00f  1
   3  3.8665095e+00 6.50e-03 5.07e-02 -11.0 6.11e+00    -  6.16e-01 1.00e+00h  1
   4  1.9535736e+00 1.08e+00 5.00e-02 -11.0 7.90e+01    -  8.03e-08 8.47e-01f  1
   5  1.9164656e+00 2.81e-02 2.96e-03 -11.0 8.23e+00    -  9.29e-01 1.00e+00h  1
   6  1.9165046e+00 3.14e-06 2.47e-06 -11.0 1.30e-01    -  9.99e-01 1.00e+00h  1
   7  1.9165046e+00 3.98e-12 1.05e-13 -11.0 1.37e-04    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 7

                                   (scaled)                 (unscaled)
Objective...............:   1.9165045626427855e+00    1.9165045626427855e+00
Dual infeasibility......:   1.0451206186354420e-13    1.0451206186354420e-13
Constraint violation....:   3.9843683907747618e-12    3.9843683907747618e-12
Variable bound violation:   3.4981688941115863e-06    3.4981688941115863e-06
Complementarity.........:   1.4336253488765137e-11    1.4336253488765137e-11
Overall NLP error.......:   1.4336253488765137e-11    1.4336253488765137e-11


Number of objective function evaluations             = 8
Number of objective gradient evaluations             = 8
Number of equality constraint evaluations            = 8
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 8
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 7
Total seconds in IPOPT                               = 0.004

EXIT: Optimal Solution Found
```

## Model Predictive Control Script 
The steps outlined in the sections above can be captured in a loop. This is summarized in the model predictive control script below :
```julia
using InfiniteOpt, Ipopt, DifferentialEquations

# DEFINE THE PROBLEM CONSTANTS
F = 100       # m³/s
V = 100       # m³
rho = 1000      # kg/m³
cp = 0.239    # J/kg K
Hr = 5e4       # Heat of reaction J/mol
E = 8750     # E/R K
k₀ = 7.2e10   # Pre-exponential factor 1/s
Ua = 5e4      # Heat transfer coefficient J/s K
cf = 1.0     # Feed concentration mol/m³
Tf = 350     # Feed temperature K
init = [0.9, 305, 300]  # Initial conditions
Ca_k, T_k = init[1:2]

# DEFINE MPC PARAMETERS
Δt = 0.1      # Control interval
Dt = 3        # Prediction horizon
t0 = 0        # Initial simulation time
tf = 2       # Final simulation time
t_vals = collect(t0:Δt:tf)

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

function cstrDynamics!(dx, x, p, t)
    Ca, T, Tc = x
    k = k₀ * exp(-E/T)
    rate = k * Ca
    dx[1] = F/V * (cf - Ca) - V * rate
    dx[2] = (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T))
end

function cstrSim(x0, tspan)
    prob = ODEProblem(cstrDynamics!, x0, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)
    xsol = sol.u[end]
    return xsol
end

# INITIALIZE THE MODEL
model = InfiniteModel(TranscriptionBackend(Ipopt.Optimizer, update_parameter_functions = true))

# INITIALIZE THE PARAMETERS
@infinite_parameter(model, t ∈ [0, Dt],
                supports = collect(0:Δt:Dt),
                derivative_method = OrthogonalCollocation(3))
@finite_parameter(model, Ca0 == init[1])
@finite_parameter(model, T0 == init[2])
@parameter_function(model, Tsp == t -> setpoint(t, 0))

# INITIALIZE THE VARIABLES
@variables(model, begin
    0 ≤ Ca, Infinite(t), (start = init[1])
    273.15 ≤ T ≤ 400, Infinite(t), (start = init[2])
    250 ≤ Tc ≤ 350, Infinite(t), (start = init[3])
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

# ADJUST DEGREES OF FREEDOM FOR CONTROL VARIABLES
constant_over_collocation.(Tc, t)

# MPC LOOP
for i in eachindex(t_vals)
    # SOLVE THE MODEL
    optimize!(model)

    # GET THE OPTIMAL INPUT
    Tc_opt = value.(Tc)[2]
    
    # SIMULATE SYSTEM FORWARD
    tk = t_vals[i]  # Offset
    tspan = (tk, tk + Δt)
    xk = [Ca_k, T_k, Tc_opt]
    xsol = cstrSim(xk, tspan)
    global Ca_k, T_k = xsol[1:2]

    # WARMSTART MODEL FOR NEXT SOLVE
    warmstart_backend_start_values(model)

    # UPDATE PARAMETERS
    set_parameter_value(Ca0, Ca_k)
    set_parameter_value(T0, T_k)
    newTsp = (t) -> setpoint(t, tk + Δt)
    set_parameter_value(Tsp, newTsp)

    # ADJUST SOLVER PARAMETERS
    set_optimizer_attribute(model, "bound_push", 1e-8)
    set_optimizer_attribute(model, "bound_frac", 1e-8)
    set_optimizer_attribute(model, "mu_init", 1e-11)
end
```

## GPU-Accelerated Resolves
### Software Setup
We can also facilitate resolves on GPUs. First, we'll need to ensure the following packages are installed:
- `InfiniteExaModels.jl`
- `MadNLP.jl`
- `CUDA.jl`

!!! note
    Currently, this workflow is only available on NVIDIA GPUs that support CUDA.

### Problem setup
Now, we'll need to initialize our model with an `ExaTranscriptionBackend` that employs a `CUDABackend`. This will transcribe the InfiniteOpt problem into an ExaModel that is GPU compatible. We'll also use a GPU-based solver like `MadNLP.jl`.
```julia
julia> using InfiniteExaModels, InfiniteExaModels, MadNLPGPU, CUDA

julia> model = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver, backend = CUDABackend()));
```

Since we'll be indexing `Tc` for the optimal control input, we'll also need to set `allowscalar` to true.
```julia
julia> CUDA.allowscalar(true);
```

That's it for the GPU setup! From here, we can follow the same steps as above to formulate our problem.

### Model Predictive Control Script
The GPU version of the MPC script is given below:
```julia
using InfiniteOpt, InfiniteExaModels, MadNLPGPU, CUDA, DifferentialEquations

# DEFINE THE PROBLEM CONSTANTS
F = 100       # m³/s
V = 100       # m³
rho = 1000      # kg/m³
cp = 0.239    # J/kg K
Hr = 5e4       # Heat of reaction J/mol
E = 8750     # E/R K
k₀ = 7.2e10   # Pre-exponential factor 1/s
Ua = 5e4      # Heat transfer coefficient J/s K
cf = 1.0     # Feed concentration mol/m³
Tf = 350     # Feed temperature K
init = [0.9, 305, 300]  # Initial conditions
Ca_k, T_k = init[1:2]

# DEFINE MPC PARAMETERS
Δt = 0.1      # Control interval
Dt = 3        # Prediction horizon
t0 = 0        # Initial simulation time
tf = 2       # Final simulation time
t_vals = collect(t0:Δt:tf)

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

function cstrDynamics!(dx, x, p, t)
    Ca, T, Tc = x
    k = k₀ * exp(-E/T)
    rate = k * Ca
    dx[1] = F/V * (cf - Ca) - V * rate
    dx[2] = (1/(rho * cp * V)) * (F * rho * cp * (Tf - T) + V * Hr * rate + Ua * (Tc - T))
end

function cstrSim(x0, tspan)
    prob = ODEProblem(cstrDynamics!, x0, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)
    xsol = sol.u[end]
    return xsol
end

# INITIALIZE THE MODEL
model = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver, backend = CUDABackend()))
CUDA.allowscalar(true)

# INITIALIZE THE PARAMETERS
@infinite_parameter(model, t ∈ [0, Dt],
                supports = collect(0:Δt:Dt),
                derivative_method = OrthogonalCollocation(3))
@finite_parameter(model, Ca0 == init[1])
@finite_parameter(model, T0 == init[2])
@parameter_function(model, Tsp == t -> setpoint(t, 0))

# INITIALIZE THE VARIABLES
@variables(model, begin
    0 ≤ Ca, Infinite(t), (start = init[1])
    273.15 ≤ T ≤ 400, Infinite(t), (start = init[2])
    250 ≤ Tc ≤ 350, Infinite(t), (start = init[3])
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

# ADJUST DEGREES OF FREEDOM FOR CONTROL VARIABLES
constant_over_collocation.(Tc, t)

# MPC LOOP
for i in eachindex(t_vals)
    # SOLVE THE MODEL
    optimize!(model)

    # GET THE OPTIMAL INPUT
    Tc_opt = value.(Tc)[2]
    
    # SIMULATE SYSTEM FORWARD
    tk = t_vals[i]  # Offset
    tspan = (tk, tk + Δt)
    xk = [Ca_k, T_k, Tc_opt]
    xsol = cstrSim(xk, tspan)
    global Ca_k, T_k = xsol[1:2]

    # WARMSTART MODEL FOR NEXT SOLVE
    warmstart_backend_start_values(model)

    # UPDATE PARAMETERS
    set_parameter_value(Ca0, Ca_k)
    set_parameter_value(T0, T_k)
    newTsp = (t) -> setpoint(t, tk + Δt)
    set_parameter_value(Tsp, newTsp)

    # ADJUST SOLVER PARAMETERS
    set_optimizer_attribute(model, "mu_init", 2e-2)
end
```
!!! note
    Although `MadNLP.jl` does support `bound_push` and `bound_fac` as options, they currently do not have an effect in resolves.
