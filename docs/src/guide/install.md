# Installation Guide
A complete guide to installing all that is necessary to use `InfiniteOpt.jl`.

## Install Julia
`InfiniteOpt` is a [Julia](https://julialang.org/) package. Thus, we first need
have an installation of Julia to get started. An appropriate download can be
found [here](https://julialang.org/downloads/). Alternatively, users with a
limited programming background might find it easier to install and use
[JuliaPro](https://juliacomputing.com/products/juliapro.html).

!!! note
    This version of `InfiniteOpt` requires that Julia 1.0 or newer be used.

## Install Packages
Now that Julia has been installed we can add the needed packages. Open up a
Julia terminal and enter the package manager:
```julia
julia> ]

(v1.2) pkg>
```
Use the `add` command in the package to manager to add the following packages:

- `JuMP`
- `Distributions`

For example, to install `JuMP` we would enter:
```julia
(v1.2) pkg> add JuMP
```

Now let's install `InfiniteOpt`, because it is not yet a registered Julia
package we have to specify the GitHub repository address:
```julia
(v1.2) pkg> add https://github.com/pulsipher/InfiniteOpt.jl
```

## Install Optimization Solvers
`InfiniteOpt` relies on solvers to solve optimization problems. Many solvers are
not native to Julia and might require commercial licenses. A list of currently
supported solvers and their corresponding Julia packages is provided in
[Supported Optimizers](@ref).

For example, we can install Ipopt which is an open-source nonlinear solver:
```julia
(v1.2) pkg> add Ipopt
```
Now Ipopt can be used as the optimizer (solver) for an infinite model by running:
```julia
julia> using InfiniteOpt, JuMP, Ipopt

julia> model = InfiniteModel(with_optimizer(Ipopt.Optimizer))
```
Most solver packages follow the `ModuleName.Optimizer` naming convention, but
this may not always be the case. See [Infinite Models](@ref) for more
information on defining infinite models and specifying solvers.
