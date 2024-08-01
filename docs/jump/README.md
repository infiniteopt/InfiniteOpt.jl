![Logo](https://raw.githubusercontent.com/infiniteopt/InfiniteOpt.jl/master/full_logo.png)

[`InfiniteOpt.jl`](https://github.com/infiniteopt/InfiniteOpt.jl) is a `JuMP` extension for expressing and solving infinite-dimensional optimization
problems. Such areas include [stochastic programming](https://en.wikipedia.org/wiki/Stochastic_programming),
[dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming),
space-time optimization, and more. `InfiniteOpt` serves as an easy-to-use modeling
interface for these advanced problem types that can be used by those with little
background in these areas. The package also contains a wealth of capabilities
making `InfiniteOpt` a powerful and convenient tool for advanced users.  

| **Current Version**                     | **Documentation**                                                               | **Build Status**                                                                                | **Citation** |
|:---------------------------------------:|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:--------------------------------------:|
| [![](https://docs.juliahub.com/InfiniteOpt/version.svg)](https://juliahub.com/ui/Packages/InfiniteOpt/p3GvY) | [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://infiniteopt.github.io/InfiniteOpt.jl/stable) | [![Build Status](https://github.com/infiniteopt/InfiniteOpt.jl/workflows/CI/badge.svg?branch=master)](https://github.com/infiniteopt/InfiniteOpt.jl/actions?query=workflow%3ACI) [![codecov.io](https://codecov.io/github/infiniteopt/InfiniteOpt.jl/coverage.svg?branch=master)](https://codecov.io/github/infiniteopt/InfiniteOpt.jl?branch=master) | [![DOI](https://img.shields.io/badge/Elsevier-CompChemEng%3A107567-yellow.svg)](https://doi.org/10.1016/j.compchemeng.2021.107567) |

`InfiniteOpt` builds upon `JuMP` to add support for many complex modeling objects which 
include:
- Infinite parameters (for example, time, space, and/or uncertainty)
- Finite parameters (similar to `ParameterJuMP`)
- Infinite variables (decision functions) (for example, ``y(t, x)``)
- Derivatives (for example, ``\frac{\partial y(t, x)}{\partial t}``)
- Measures (for example, ``\int_{t \in T}y(t, x)dt`` and ``\mathbb{E}[y(\xi)]``)

The unifying modeling abstraction behind `InfiniteOpt` captures a wide spectrum 
of disciplines which include dynamic, PDE, stochastic, and semi-infinite 
optimization. Moreover, we facilitate transferring techniques between these 
to synthesize new optimization paradigms.

## License
`InfiniteOpt` is licensed under the [MIT "Expat" license](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/LICENSE).

## Installation
`InfiniteOpt.jl` is a registered [Julia](https://julialang.org/) package and 
can be installed by entering the following in the REPL.

```julia
julia> import Pkg; Pkg.add("InfiniteOpt")
```

## Documentation
Please visit our [documentation pages](https://infiniteopt.github.io/InfiniteOpt.jl/stable) 
to learn more. These pages are quite extensive and feature overviews, guides,
manuals, tutorials, examples, and more.

## Questions
For additional help please visit and post in our 
[discussion forum](https://github.com/infiniteopt/InfiniteOpt.jl/discussions).

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
