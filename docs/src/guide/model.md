# [Infinite Models] (@id infinite_model_docs)
A guide and manual for initializing infinite models and understanding how
they work. The Datatypes and Methods sections at the end comprise the manual,
and the above sections comprise the guide.  

## Overview
Infinite models are expressed via the [`InfiniteModel`](@ref) datatype which is at the
core of `InfiniteOpt`. These model objects are designed to emulate the behavior
of [`Model`](@ref) objects in `JuMP`. These data objects store the parameters, variables,
measures, objective, constraints, and all other data used in `InfiniteOpt`. This
differs from `JuMP` models which store such information in a `MathOptInterface`
model backend.

## Basic Usage
Infinite models can be initialized with no arguments by default:
```jldoctest
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
model definition:
```jldoctest
julia> using InfiniteOpt, JuMP, Ipopt

julia> model = InfiniteModel(Ipopt.Optimizer)
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```
Note that any optimizer currently supported by `JuMP v0.19.0` or newer is
supported for use in `InfiniteOpt`. For completeness, the table of currently
supported optimizers is provided below in [Supported Optimizers](@ref).

We can also specify optimizer attributes via
[`optimizer_with_attributes`](@ref JuMP.optimizer_with_attributes(::Any, ::Pair))
which allows us to append as many attributes as we like, for example:
```jldoctest
julia> using InfiniteOpt, JuMP, Ipopt

julia> model = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer,
                                                       "output_level" => 0))
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```

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
[Model Transcription](@ref transcription_docs).

All the arguments used with the `InfiniteModel` constructor (e.g., the optimizer)
are simply passed on and stored in the optimizer model backend. Thus, any
argument supported by [`JuMP.Model`](@ref) can be passed on to the optimizer
model by including it in the `InfiniteModel` constructor. For example, we can
specify the `caching_mode` keyword argument in the `InfiniteModel` call to use
in the definition of the optimizer model:
```jldoctest
julia> using InfiniteOpt, JuMP, Ipopt, MathOptInterface

julia> const MOIU = MathOptInterface.Utilities;

julia> model = InfiniteModel(Ipopt.Optimizer,
                             caching_mode = MOIU.MANUAL)
An InfiniteOpt Model
Feasibility problem with:
Variables: 0
Optimizer model backend information:
Model mode: MANUAL
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Ipopt
```
Notice that the model mode of the optimizer model is now `MANUAL`.

## Supported Optimizers
`InfiniteOpt` can use any optimizer that is supported by `JuMP v0.19.0` or newer.
In spirit of providing complete documentation, the table of optimizers currently
supported by [`JuMP.jl`](https://github.com/JuliaOpt/JuMP.jl) is provided below.
This information comes directly from their documentation pages. Please refer to
`JuMP` documentation for additional optimizer information on installation and use.

| Solver                                                                         | Julia Package                                                                    | License  | Supports                           |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | -------- | ---------------------------------- |
| [Artelys Knitro](https://www.artelys.com:443/solvers/knitro/)                  | [KNITRO.jl](https://github.com/JuliaOpt/KNITRO.jl)                               | Comm.    | LP, MILP, SOCP, MISOCP, NLP, MINLP |
| [Cbc](https://projects.coin-or.org/Cbc)                                        | [Cbc.jl](https://github.com/JuliaOpt/Cbc.jl)                                     | EPL      | MILP                               |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              | GPL      | LP, SOCP, SDP                      |
| [CDD](https://github.com/cddlib/cddlib)                                        | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)                         | GPL      | LP                                 |
| [Clp](https://projects.coin-or.org/Clp)                                        | [Clp.jl](https://github.com/JuliaOpt/Clp.jl)                                     | EPL      | LP                                 |
| [COSMO](https://github.com/oxfordcontrol/COSMO.jl)                             | [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                            | Apache   | LP, QP, SOCP, SDP                  |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)                         | [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl)                                 | Comm.    | LP, MILP, SOCP, MISOCP             |
| [CSDP](https://github.com/coin-or/csdp/)                                       | [CSDP.jl](https://github.com/JuliaOpt/CSDP.jl)                                   | EPL      | LP, SDP                            |
| [ECOS](https://github.com/embotech/ecos)                                       | [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl)                                   | GPL      | LP, SOCP                           |
| [FICO Xpress](https://www.fico.com/en/products/fico-xpress-optimization)       | [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl)                               | Comm.    | LP, MILP, SOCP, MISOCP             |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl)                                   | GPL      | LP, MILP                           |
| [Gurobi](https://www.gurobi.com/)                                              | [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl)                               | Comm.    | LP, MILP, SOCP, MISOCP             |
| [Ipopt](https://projects.coin-or.org/Ipopt)                                    | [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl)                                 | EPL      | LP, QP, NLP                        |
| [Juniper](https://github.com/lanl-ansi/Juniper.jl)                             | [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                            | MIT      | MISOCP, MINLP                      |
| [MOSEK](https://www.mosek.com/)                                                | [MosekTools.jl](https://github.com/JuliaOpt/MosekTools.jl)                       | Comm.    | LP, MILP, SOCP, MISOCP, SDP        |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              | Apache   | LP, QP                             |
| [ProxSDP](https://github.com/mariohsouto/ProxSDP.jl)                           | [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                          | MIT      | LP, SOCP, SDP                      |
| [SCIP](https://scip.zib.de/)                                                   | [SCIP.jl](https://github.com/SCIP-Interfaces/SCIP.jl)                            | ZIB      | MILP, MINLP                        |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/JuliaOpt/SCS.jl)                                     | MIT      | LP, SOCP, SDP                      |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/JuliaOpt/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl)                                   | GPL      | LP, SDP                            |
| [SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)               | [SDPNAL.jl](https://github.com/JuliaOpt/SDPNAL.jl)                               | CC BY-SA | LP, SDP                            |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/JuliaOpt/SDPT3.jl)                                 | GPL      | LP, SOCP, SDP                      |
| SeDuMi                                                                         | [SeDuMi.jl](https://github.com/JuliaOpt/SeDuMi.jl)                               | GPL      | LP, SOCP, SDP                      |
| [Tulip](https://github.com/ds4dm/Tulip.jl)                                     | [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                    | MPL-2    | LP                                 |


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
