![Logo](full_logo.png)
---

A `JuMP` extension for expressing and solving infinite-dimensional optimization
problems. Such areas include [stochastic programming](https://en.wikipedia.org/wiki/Stochastic_programming),
[dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming),
space-time optimization, and more. `InfiniteOpt` serves as an easy to use modeling
interface for these advanced problem types that can be used by those with little
to no background in these areas. It also it contains a wealth of capabilities
making it a powerful and convenient tool for advanced users.  

| **Documentation**                                                               | **Build Status**                                                                                | **Citation** |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:--------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/stable) | [![Build Status](https://api.travis-ci.com/pulsipher/InfiniteOpt.jl.svg?branch=v0.3.2)](https://travis-ci.com/pulsipher/InfiniteOpt.jl) [![Build Status2](https://ci.appveyor.com/api/projects/status/p3srfp3uuvchfg3j/branch/v0.3.2?svg=true)](https://ci.appveyor.com/project/pulsipher/InfiniteOpt-jl) [![codecov.io](https://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=release-0.3)](https://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=release-0.3) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291106.svg)](https://doi.org/10.5281/zenodo.4291106) |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev) | [![Build Status](https://github.com/pulsipher/InfiniteOpt.jl/workflows/CI/badge.svg?branch=master)](https://github.com/pulsipher/InfiniteOpt.jl/actions?query=workflow%3ACI) [![codecov.io](https://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=master)](https://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=master) | |

It builds upon `JuMP` to add support for many complex modeling objects which 
include:
- Infinite parameters (e.g., time, space, uncertainty, etc.)
- Finite parameters (similar to `ParameterJuMP`)
- Infinite variables (e.g., ``y(t, x)``)
- Derivatives (e.g., ``\frac{\partial y(t, x)}{\partial t}``)
- Measures (e.g., ``\int_{t \in \mathcal{D}_t}y(t,x) dt``, ``\mathbb{E}[y(\xi)]``)
- More

Moreover, `InfiniteOpt`'s default direct transcription (i.e., discretization) 
features include:
- Efficient implementations that scale **linearly**!
- Diverse integral approximations (e.g., quadratures, sampling)
- Diverse derivative approximations (e.g., finite difference, orthogonal 
  collocation)
- Sophisticated support point management system
- Compatible with all [JuMP-supported solvers](https://jump.dev/JuMP.jl/v0.21.8/installation/#Supported-solvers)

Currently, the following infinite and finite problem types are accepted:
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
    - Higher-order powers (via place holder variables)
- Constraints
    - Linear
    - Quadratic (convex and non-convex)
    - Higher-order powers (via place holder variables)
    - Conic
    - Semi-definite
    - Indicator

Comments, suggestions and improvements are welcome and appreciated.

## License
`InfiniteOpt` is licensed under the [MIT "Expat" license](./LICENSE).

## Installation
`InfiniteOpt.jl` is a registered package and can be installed by entering the
following in the package manager.

```julia
(v1.6) pkg> add InfiniteOpt
```

## Modeling Infinite-Dimensional Problems with InfiniteOpt.jl
See our YouTube overview of infinite-dimensional programming and InfiniteOpt.jl's 
capabilities (Note the syntax in the video is now deprecated):
[![youtube](docs/src/assets/youtube.PNG)](http://www.youtube.com/watch?v=q5ETFLZbxiU "Modeling Infinite-Dimensional Problems with InfiniteOpt.jl")

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/stable)

Please visit our [documentation pages](https://pulsipher.github.io/InfiniteOpt.jl/stable) 
to learn more. These pages are quite extensive and feature overviews, guides,
manuals, tutorials, examples, and more!

## Citing
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291106.svg)](https://doi.org/10.5281/zenodo.4291106)

If you use InfiniteOpt.jl in your research, we would greatly appreciate your 
citing it.
```latex
@misc{pulsipher20,
  author       = {J. Pulsipher and W. Zhang and V. Zavala},
  title        = {InfiniteOpt.jl -- A Julia package for modeling infinite-dimensional optimization problems},
  year         = 2020,
  doi          = {10.5281/zenodo.4291106},
  url          = {https://doi.org/10.5281/zenodo.4291106}
}
```

## Project Status
The package is tested against Julia `1.0` and `1.6` on Linux and Windows.

## Contributing
`InfiniteOpt` is being actively developed and suggestions or other forms of contribution are encouraged.
There are many ways to contribute to this package. For more information please
visit [CONTRIBUTING](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/CONTRIBUTING.md).
