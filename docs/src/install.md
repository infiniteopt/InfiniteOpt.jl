# Installation Guide
A complete guide to installing all that is necessary to use `InfiniteOpt.jl`.

## Install Julia
`InfiniteOpt` is a [Julia](https://julialang.org/) package. Thus, we first need
have an installation of Julia to get started. An appropriate download can be
found [here](https://julialang.org/downloads/). We recommend using 
[VSCode](https://www.julia-vscode.org/) to edit and run Julia scripts. 
Alternatively, users with a limited programming background might find it easier 
to install and use [JuliaPro](https://juliacomputing.com/products/juliapro/) 
(however this may introduce performance degradations and compatibility issues).

!!! note
    This version of `InfiniteOpt` requires that Julia 1.0 or newer be used.

## Install Packages
Now that Julia has been installed we can add the needed packages. Open up a
Julia terminal and enter the package manager:
```julia-repl
julia> ]

(v1.6) pkg>
```

!!! tip
    We recommend you create a Pkg _environment_ for each project you use `InfiniteOpt`
    for, instead of adding lots of packages to the global environment. The
    [Pkg manager documentation](https://pkgdocs.julialang.org/v1/environments/)
    has more information on this topic.

Use the `add` command in the package to manager to add the following packages:

- `Distributions` (required for stochastic programming)

For example, to install `Distributions` we would enter:
```julia-repl
(v1.6) pkg> add Distributions
```

Now let's install `InfiniteOpt`:
```julia-repl
(v1.6) pkg> add InfiniteOpt
```

!!! info
    Installation troubles? Check the [Common Installation Problems](@ref) section
    below.

Alternatively, we can install the current experimental version of 
`InfiniteOpt` via:
```julia-repl
(v1.6) pkg> add https://github.com/pulsipher/InfiniteOpt.jl
```

## Install Optimization Solvers
`InfiniteOpt` relies on solvers to solve optimization problems. Many solvers are
not native to Julia and might require commercial licenses. A list of currently
supported solvers and their corresponding Julia packages is provided in
[Supported Optimizers](@ref).

For example, we can install Ipopt which is an open-source nonlinear solver:
```julia-repl
(v1.6) pkg> add Ipopt
```
Now Ipopt can be used as the optimizer (solver) for an infinite model by running:
```julia-repl
julia> using InfiniteOpt, Ipopt

julia> model = InfiniteModel(Ipopt.Optimizer)
```
Most solver packages follow the `ModuleName.Optimizer` naming convention, but
this may not always be the case. See [Infinite Models](@ref infinite_model_docs)
for more information on defining infinite models and specifying solvers.

## Common Installation Problems
!!! tip
    When in doubt, run `import Pkg; Pkg.update()` to see if updating your
    packages fixes the issue. Remember you will need to exit Julia and start a
    new session for the changes to take effect.

### Check the version of your packages
Each package is versioned with a [three-part number](https://semver.org) of the
form `vX.Y.Z`. You can check which versions you have installed with:
```julia-repl
julia> ]

(v1.6) pkg> status
```
This should almost always be the most-recent release. You can check the releases
of a package by going to the relevant Github page, and navigating to the
"releases" page. For example, the list of JuMP releases is available at:
[https://github.com/pulsipher/InfiniteOpt.jl/releases](https://github.com/pulsipher/InfiniteOpt.jl/releases).

If you need to ask question for help, please include the output of `status`!

### Unsatisfiable requirements detected
Did you get an error like 
`Unsatisfiable requirements detected for package InfiniteOpt`? The Pkg 
documentation has a 
[section on how to understand and manage these conflicts](https://pkgdocs.julialang.org/v1/managing-packages/).

Typically, these conflicts can be resolved by using 
[package environments](https://pkgdocs.julialang.org/v1/environments/).

### Installing new packages can make InfiniteOpt downgrade to an earlier version
Another common issue is that after adding a new package, code that previously 
worked no longer runs.

This usually happens because the new package is not compatible with the latest
version of `InfiniteOpt`. Therefore, the package manager downgrades `InfiniteOpt` 
to an earlier version!

Thus, please Pay careful attention to the output of the package manager when 
adding new packages, especially when you see a package being downgraded!