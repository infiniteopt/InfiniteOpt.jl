# Infinite Models
A guide and manual for initializing infinite models and understanding how
they work.

## Overview
Infinite models are expressed via the [`InfiniteModel`](@ref) datatype which is at the
core of `InfiniteOpt`. These model objects are designed to emulate the behavior
of [`Model`](@ref) objects in `JuMP`. These data objects store the parameters, variables,
measures, objective, constraints, and all other data used in `InfiniteOpt`. This
differs from `JuMP` models which store such information in a `MathOptInterface`
model backend.

## Basic Usage
Infinite models can be initialized with no arguments by default:
```julia
julia> using InfiniteOpt, JuMP

julia> model = InfiniteModel()
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```

The optimizer that will be used to solve the model can also be specified at
model definition using `JuMP`'s [`with_optimizer`](@ref) constructor:
```julia
julia> using InfiniteOpt, JuMP, Gurobi

julia> model = InfiniteModel(with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Gurobi
```
Note that any optimizer currently supported by `JuMP v0.19.0` or newer is
supported for use in `InfiniteOpt`. For completeness, the table of currently
supported optimizers is provided below in [Supported Optimizers](@ref).

Now you have an initialized `InfiniteModel` that is ready for your mathematical
model to be defined and optimized!

## Advanced Definition Information
As you may have noticed in the above examples, `InfiniteModel`s contain an
optimizer model backend which simply corresponds to a `JuMP` [`Model`](@ref) that
will be used to store and optimize the reformulation of the infinite mathematical
model stored in `InfiniteModel`. It also will contain a mapping between its
optimization model and that of the `InfiniteModel` (e.g., a mapping between the
variables and constraints). By default, `InfiniteModel`s use a
[`TranscriptionModel`](@ref) optimizer model backend which will store a
transcribed (discretized) version of the infinite model. More information on
the internal use of `TranscriptionModel`s is provided in
[Model Transcription](@ref).

All the arguments used with the `InfiniteModel` constructor (e.g., the optimizer)
are simply passed on and stored in the optimizer model backend. Thus, any
argument supported by [`JuMP.Model`](@ref) can be passed on to the optimizer
model by including it in the `InfiniteModel` constructor. For example, we can
specify the `caching_mode` keyword argument in the `InfiniteModel` call to use
in the definition of the optimizer model:
```julia
julia> using InfiniteOpt, JuMP, Gurobi, MathOptInterface

julia> const MOIU = MathOptInterface.Utilities

julia> model = InfiniteModel(with_optimizer(Gurobi.Optimizer),
                             caching_mode = MOIU.MANUAL)
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: MANUAL
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Gurobi
```
Notice that the model mode of the optimizer model is now `MANUAL`.

## Supported Optimizers
`InfiniteOpt` can use any optimizer that is supported by `JuMP v0.19.0` or newer.
In spirit of providing complete documentation, the table of optimizers currently
supported by [`JuMP.jl`](https://github.com/JuliaOpt/JuMP.jl) is provided below.
This information comes directly from their documentation pages. Please refer to
`JuMP` documentation for additional optimizer information on installation and use.

| Solver                                                                         | Julia Package                                                                    | License | Supports                           |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------- | ---------------------------------- |
| [Artelys Knitro](https://www.artelys.com/knitro)                               | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP, NLP, MINLP |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Cbc.jl)                                     | EPL     | MILP                               |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              | GPL     | LP, SOCP, SDP                      |
| [Clp](https://projects.coin-or.org/Clp)                                        | [Clp.jl](https://github.com/JuliaOpt/Clp.jl)                                     | EPL     | LP                                 |
| [COSMO](https://github.com/oxfordcontrol/COSMO.jl)                             | [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                            | Apache  | LP, QP, SOCP, SDP                  |
| [CPLEX](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)                                 | Comm.   | LP, MILP, SOCP, MISOCP             |
| [CSDP](https://projects.coin-or.org/Csdp/)                                     | [CSDP.jl](https://github.com/JuliaOpt/CSDP.jl)                                   | EPL     | LP, SDP                            |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl)                                   | GPL     | LP, SOCP                           |
| [FICO Xpress](http://www.fico.com/en/products/fico-xpress-optimization-suite)  | [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP             |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl)                                   | GPL     | LP, MILP                           |
| [Gurobi](http://gurobi.com)                                                    | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)                               | Comm.   | LP, MILP, SOCP, MISOCP             |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)                                 | EPL     | LP, QP, NLP                        |
| [Juniper](https://github.com/lanl-ansi/Juniper.jl)                             | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                            | MIT     | MISOCP, MINLP                      |
| [MOSEK](http://www.mosek.com/)                                                 | [MosekTools.jl](https://github.com/JuliaOpt/MosekTools.jl)                       | Comm.   | LP, MILP, SOCP, MISOCP, SDP        |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              | Apache  | LP, QP                             |
| [ProxSDP](https://github.com/mariohsouto/ProxSDP.jl)                           | [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                          | MIT     | LP, SOCP, SDP                      |
| [SCIP](https://scip.zib.de/)                                                   | [SCIP.jl](https://github.com/SCIP-Interfaces/SCIP.jl)                            | ZIB     | MILP, MINLP                        |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/JuliaOpt/SCS.jl)                                     | MIT     | LP, SOCP, SDP                      |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl)                                   | GPL     | LP, SDP                            |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/JuliaOpt/SeDuMi.jl)                               | GPL     | LP, SOCP, SDP                      |


Where:

-   LP = Linear programming
-   QP = Quadratic programming
-   SOCP = Second-order conic programming (including problems with convex
    quadratic constraints and/or objective)
-   MILP = Mixed-integer linear programming
-   NLP = Nonlinear programming
-   MINLP = Mixed-integer nonlinear programming
-   SDP = Semidefinite programming
-   MISDP = Mixed-integer semidefinite programming

You may also use [AmplNLWriter](https://github.com/JuliaOpt/AmplNLWriter.jl) to
access solvers that support the [nl format](https://en.wikipedia.org/wiki/Nl_(format)).
Such solvers include [Bonmin](https://projects.coin-or.org/Bonmin) and
[Couenne](https://projects.coin-or.org/Couenne). See a more complete list
[here](https://ampl.com/products/solvers/all-solvers-for-ampl/).

## Datatypes
```@docs
InfiniteModel
```

## Methods
```@docs
InfiniteModel()
```
