![Logo](assets/full_logo.png)

A `JuMP` extension for expressing and solving infinite-dimensional optimization 
problems.

## Overview
`InfiniteOpt.jl` provides a mathematical interface to express and solve 
optimization problems that entail an infinite-dimensional decision space. Such 
problems stem from areas such as space-time programming and 
stochastic programming. `InfiniteOpt` is meant to facilitate intuitive model 
definition, automatic transcription into solvable models, permit a wide range 
of user-defined extensions/behavior, and more. Currently, its capabilities 
include:
- `JuMP`-like symbolic macro interface
- Infinite domain abstractions for parameterization of variables/constraints
- Finite parameters support and use (similar to `ParameterJuMP`)
- Direct support of infinite, semi-infinite, point, and finite variables
- Straightforward measure operator definition (e.g., integrals, risk measures)
- Infinite/finite constraint definition
- Event constraint definition (e.g., chance constraints)
- Compact ordinary/partial differential operator expression
- Efficient automated model transcription/reformulation and solution
- Compatible with all [JuMP-supported solvers](https://www.juliaopt.org/JuMP.jl/dev/installation/#Getting-Solvers-1)
- Readily extendable to accommodate user defined abstractions and solution techniques.

!!! note 
    Nonlinear objects as defined by `JuMP.@NL[macro_name]` are not currently 
    supported by `InfiniteOpt`. See [Nonlinear Expressions](@ref) for more 
    information and possible workarounds. 

### Modeling Infinite-Dimensional Problems with InfiniteOpt.jl
See our YouTube overview of infinite-dimensional programming and InfiniteOpt.jl's 
capabilities: 
[![youtube](assets/youtube.PNG)](https://www.youtube.com/watch?v=q5ETFLZbxiU "Modeling Infinite-Dimensional Problems with InfiniteOpt.jl")

## Installation
`InfiniteOpt.jl` is a registered `Julia` package and can be added simply by 
inputting the following in the package manager:
```julia
(v1.6) pkg> add InfiniteOpt
```
Please visit our [Installation Guide](@ref) for more details and information
on how to get started.

Moreover, `InfiniteOpt` is under constant develop with new features being added 
often. Thus, the latest pre-release experimental version can be obtained via the 
following command:
```julia
(v1.6) pkg> add https://github.com/pulsipher/InfiniteOpt.jl
```

## How to Use the Documentation
`InfiniteOpt` is intended to serve both as a high-level interface for 
infinite-dimensional optimization and as a highly customizable/extendable 
platform for implementing advanced techniques. With this in mind, we provide the 
`User Guide` sections to walk through the ins and outs of `InfiniteOpt`. Each 
page in the `User Guide` typically contains the following:
- An `Overview` section describing the purpose of the page (at the top)
- A `Basic Usage` section to guide using `InfiniteOpt` at a high level 
  (near the top)
- `Methods` and/or `DataTypes` sections serving as a technical manual for all the 
  public methods and datatypes (at the bottom)
- Other sections offering more in-depth information/guidance beyond basic usage 
  (in the middle)

Details, instructions, templates, and tutorials on how to write user-defined 
extensions in `InfiniteOpt` are provided on the [Extensions](@ref) page.

Finally, case study examples are provided on the [Examples](@ref) page.

## Contribution
`InfiniteOpt` is a powerful tool with a broad scope lending to a large realm of 
possible feature additions and enhancements. So, we are thrilled to support 
anyone who would like to contribute to this project in any way big or small.

For small documentation fixes (such as typos or wording clarifications) please 
do the following:
1. Click on `Edit on GitHub` at the top of the documentation page
2. Make the desired changes
3. Submit a pull request

For other contributions, please visit our [Developers Guide](@ref) for step by 
step instructions and to review our style guide.

## Acknowledgements
We acknowledge our support from the Department of Energy under grant 
DE-SC0014114.
