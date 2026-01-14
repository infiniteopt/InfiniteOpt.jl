```@raw html
<img class="display-light-only" src="assets/full_logo.png" alt="InfiniteOpt logo"/>
<img class="display-dark-only" src="assets/full_logo-dark.png" alt="InfiniteOpt logo"/>
```

A `JuMP` extension for modeling and efficiently solving infinite-dimensional optimization 
problems.

!!! note
    `InfiniteOpt v0.6` introduces quite a few new features and makes some 
    breaking changes to the underlying data structures. The core API remains 
    largely the same.

## What is InfiniteOpt?
`InfiniteOpt.jl` provides a general modeling environment to express and efficiently solve 
infinite-dimensional optimization problems (e.g., dynamic, stochastic, and 
PDE-constrained optimization). `InfiniteOpt` is meant to facilitate intuitive model 
definition that is compatible with a large suite of solution techniques 
(including automatic transcription) to tackle large-scale real world problems. 

It builds upon `JuMP` to add support for many complex modeling objects which 
include:
- Infinite parameters (e.g., time, space, uncertainty, etc)
- Infinite variables (e.g., ``y(t, x)``)
- Derivatives (e.g., ``\frac{\partial y(t, x)}{\partial t}``)
- Measures (e.g., ``\int_{t \in \mathcal{D}_t}y(t,x) dt``, ``\mathbb{E}[y(\xi)]``)
- Chance constraints
- More

Moreover, `InfiniteOpt` modeling interface decouples infinite-dimensional 
formulations from the solution/transformations used to solve them. This enables rapid 
diverse switching between techniques with state-of-the-art implementations. 
By default, we employ direct transcription (i.e., discretization) transformations 
whose features include:
- Highly efficient automatic transformation that scales linearly 
- Large suite of approximation schemes (e.g., quadratures, orthogonal collocation)
- Compatible with all [JuMP-supported solvers](https://jump.dev/JuMP.jl/v1/installation/#Supported-solvers)
- Compatible with GPU solvers via [InfiniteExaModels](https://github.com/infiniteopt/InfiniteExaModels.jl) that make InfiniteOpt **orders-of-magnitude faster than state-of-the-art** alternatives.

Accepted infinite/finite problem forms currently include:
- Variables
    - Continuous and semi-continuous
    - Binary
    - Integer and semi-integer
    - Semi-definite
- Derivatives
    - Ordinary derivative operators (of any order)
    - Partial derivative operators (of any order)
- Measures
    - Univariate and multivariate integrals 
    - Univariate and multivariate expectations 
    - Arbitrary measure operators (via general measure API)
- Objectives
    - Linear
    - Quadratic (convex and non-convex)
    - General nonlinear
- Constraints
    - Linear
    - Quadratic (convex and non-convex)
    - General nonlinear
    - Conic
    - Semi-definite
    - Indicator
    - Anything else supported by JuMP
- Machine learning surrogate models (via [MathOptAI](@ref mathoptai_guide))
- Generalized disjunctive programming (via [DisjunctiveProgramming.jl](https://github.com/infiniteopt/DisjunctiveProgramming.jl))

### Infinite-Dimensional Optimization with InfiniteOpt.jl
See our YouTube overview of infinite-dimensional programming and InfiniteOpt.jl's 
capabilities from our JuMP-dev 2025 talk: 
[![youtube](assets/youtube.jpg)](https://www.youtube.com/watch?v=oug8ttdzc2g&t=0s "The State of InfiniteOpt.jl")

## Installation
`InfiniteOpt.jl` is a registered `Julia` package and can be added simply by 
inputting the following in the package manager:
```julia
(v1.12) pkg> add InfiniteOpt
```
Please visit our [Installation Guide](@ref) for more details and information
on how to get started.

## First Steps 
`InfiniteOpt` is extension of `JuMP`, thus new users should start by 
familiarizing themselves with how to use `JuMP` first. See 
[`JuMP`'s documentation](https://jump.dev/JuMP.jl/v1/) to learn more.

Once you're familiar with `JuMP` check out our [Quick Start Guide](@ref) to get 
started. From there we provide extensive documentation to help you along, please 
see [How to Use the Documentation](@ref) for more info.

## How to Use the Documentation
`InfiniteOpt` is intended to serve both as a high-level interface for 
infinite-dimensional optimization and as a highly customizable/extendable 
platform for implementing advanced techniques. With this in mind, we provide the 
`User Guide` sections to walk through the ins and outs of `InfiniteOpt`. Each 
page in the `User Guide` typically contains the following:
- An `Overview` section describing the purpose of the page.
- A `Basic Usage` section to guide using `InfiniteOpt` at a high level.
- Other sections offering more in-depth information/guidance beyond basic usage.

We also provide a technical manual in the `API Manual` sections which are 
comprised of the docstrings of all the public methods and types that comprise 
`InfiniteOpt`. These detail the technical aspects of each function, method, macro, 
and composite type.

Details, instructions, templates, and tutorials on how to write user-defined 
extensions in `InfiniteOpt` are provided on the [Extensions](@ref) page.

Finally, case study examples are provided in the `Examples` sections.

## Questions
For additional help please visit and post in our 
[Discussion Forum](https://github.com/infiniteopt/InfiniteOpt.jl/discussions).

## Contribution
`InfiniteOpt` is a powerful tool with a broad scope lending to a large realm of 
possible feature additions and enhancements. So, we are thrilled to support 
anyone who would like to contribute to this project in any way big or small.

For small documentation fixes (such as typos or wording clarifications) please 
do the following:
1. Click on `Edit on GitHub` at the top of the documentation page
2. Make the desired changes
3. Submit a pull request

For other contributions, please visit our 
[Developers Guide](@ref contribute_guide) for step-by-step instructions and to 
review our style guide.

## Citing
[![DOI](https://img.shields.io/badge/Elsevier-CompChemEng%3A107567-yellow.svg)](https://doi.org/10.1016/j.compchemeng.2021.107567) 
[![DOI](https://img.shields.io/badge/math.OC-arXiv%3A2106.12689-B31B1B.svg)](https://arxiv.org/abs/2106.12689)

If you use InfiniteOpt.jl in your research, we would greatly appreciate your 
citing it.
```latex
@article{pulsipher2022unifying,
      title = {A unifying modeling abstraction for infinite-dimensional optimization},
      journal = {Computers & Chemical Engineering},
      volume = {156},
      year = {2022},
      issn = {0098-1354},
      doi = {https://doi.org/10.1016/j.compchemeng.2021.107567},
      url = {https://www.sciencedirect.com/science/article/pii/S0098135421003458},
      author = {Joshua L. Pulsipher and Weiqi Zhang and Tyler J. Hongisto and Victor M. Zavala},
}
```
A pre-print version is freely available though [arXiv](https://arxiv.org/abs/2106.12689).

## Acknowledgements
We acknowledge previous support from the Department of Energy under grant 
DE-SC0014114.
