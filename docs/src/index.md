![Logo](assets/full_logo.png)

A `JuMP` extension for expressing and solving infinite dimensional optimization
problems.

## Overview
`InfiniteOpt.jl` provides a mathematical interface to express and solve
optimization problems that entail an infinite dimensional decision space. Such
problems stem from areas such as dynamic programming, state-space models, and
stochastic programming. `InfiniteOpt` is meant to facilitate intuitive model
definition, automatic transcription into solvable models, permit a wide range
of user-defined extensions/behavior, and more. Currently, its capabilities
include:
- `JuMP`-like symbolic macro interface
- Infinite parameter support and parameterization of variables/constraints
- Finite parameters support(similar to `ParameterJuMP`)
- Direct support of infinite, point, and hold variables
- Symbolic measure (integral) expression
- Infinite/finite constraint definition
- Ordinary differential equation support (coming soon)
- Automated model transcription/reformulation and solution
- Readily extendable to accommodate user defined abstractions and solution techniques

!!! note
    Currently, `InfiniteOpt` only accepts linear and quadratic expressions.
    Development is underway to allow for general nonlinear constraints.  

## Installation
`InfiniteOpt.jl` is still under initial development but can be
installed by entering the following in the package manager.

```julia
(v1.3) pkg> add https://github.com/pulsipher/InfiniteOpt.jl
```

## Quick Start
Below is a brief example of the high-level API.

```julia
using InfiniteOpt, JuMP, Ipopt, Distributions

# Set the problem information
θ_nom, covar = [0.; 60.; 10.], [80. 0 0; 0 80. 0; 0 0 120.]
n_z, n_θ, n_d = 3, 3, 3

# Initialize the model
m = InfiniteModel(with_optimizer(Ipopt.Optimizer))

# Set the uncertainty parameters
dist = MvNormal(θ_nom, covar)
@infinite_parameter(m, θ[i = 1:n_θ] in dist, num_supports = 100)
@infinite_parameter(m, t in [0, 10])

# Initialize the variables
@infinite_variable(m, z[1:n_z](θ, t))
@infinite_variable(m, 0 <= y(θ) <= 100)
@hold_variable(m, d[1:n_d] >= 0)

# Set objective function
@objective(m, Min, expect(1 - y, θ))

# Set first stage constraints
@constraint(m, max_cost, sum(1 / 3 * d[i] for i = 1:n_d) <= 5)

# Set the second stage constraints
@constraint(m, f1, -z[1] - 35 - d[1] + y <= 0)
@constraint(m, f2, z[1] - 35 - d[1] + y <= 0)
@constraint(m, f3, -z[2] - 50 - d[2] + y <= 0)
@constraint(m, f4, z[1] - 50 - d[2] + y <= 0)
@constraint(m, h1, z[1] - θ[1] == 0)
@constraint(m, h2, -z[1] -z[2] + z[3] - θ[2] == 0)
@constraint(m, h3, z[2] - θ[3] == 0)

# Solve and and obtain results
optimize!(m)
if has_values(m)
    opt_y = value(y)
    opt_d = value.(d)
    opt_obj = objective_value(m)
end
```

## Acknowledgements
We acknowledge our support from the Department of Energy under grant
DE-SC0014114.
