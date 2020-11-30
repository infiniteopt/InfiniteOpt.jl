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
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev) | [![Build Status](https://travis-ci.com/pulsipher/InfiniteOpt.jl.svg?branch=master)](https://travis-ci.com/pulsipher/InfiniteOpt.jl) [![Build Status2](https://ci.appveyor.com/api/projects/status/github/pulsipher/InfiniteOpt.jl?branch=master&svg=true)](https://ci.appveyor.com/project/pulsipher/InfiniteOpt-jl) [![codecov.io](https://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=master)](https://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=master) | |

Its capabilities include:
- `JuMP`-like symbolic macro interface
- Infinite set abstractions for parameterization of variables/constraints
- Finite parameters support and use (similar to `ParameterJuMP`)
- Direct support of infinite, point, and hold variables
- Straightforward measure operator definition (e.g., integrals, risk measures)
- Infinite/finite constraint definition
- Event constraint definition (e.g., chance constraints)
- Compact ordinary/partial differential operator expression
- Efficient automated model transcription/reformulation and solution
- Compatible with all [JuMP-supported solvers](https://www.juliaopt.org/JuMP.jl/dev/installation/#Getting-Solvers-1)
- Readily extendable to accommodate user defined abstractions and solution techniques.

Currently, the following infinite and finite problem types are accepted:
- Variables
    - Continuous and semi-continuous
    - Binary
    - Integer and semi-integer
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

Comments, suggestions and improvements are welcome and appreciated.

## License
`InfiniteOpt` is licensed under the [MIT "Expat" license](./LICENSE).

## Installation
`InfiniteOpt.jl` is a registered package and can be installed by entering the
following in the package manager.

```julia
(v1.5) pkg> add InfiniteOpt
```

## Modeling Infinite-Dimensional Problems with InfiniteOpt.jl
See our YouTube overview of infinite-dimensional programming and InfiniteOpt.jl's 
capabilities:
[![youtube](docs/src/assets/youtube.PNG)](http://www.youtube.com/watch?v=q5ETFLZbxiU "Modeling Infinite-Dimensional Problems with InfiniteOpt.jl")

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/stable)

Please visit our [documentation pages](https://pulsipher.github.io/InfiniteOpt.jl/stable) to learn more. These pages are quite extensive and feature overviews, guides, manuals,
tutorials, examples, and more!

## Citing
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4291106.svg)](https://doi.org/10.5281/zenodo.4291106)

If you use InfiniteOpt.jl in your research, we would greatly appreciate your citing it.
```latex
@misc{pulsipher20,
  author       = {J. Pulsipher and W. Zhang and V. Zavala},
  title        = {InfiniteOpt.jl -- A Julia package for modeling infinite-dimensional optimizataion problems},
  year         = 2020,
  doi          = {10.5281/zenodo.4291106},
  url          = {https://doi.org/10.5281/zenodo.4291106}
}
```

## Project Status
The package is tested against Julia `1.0` and `1.5` on Linux and Windows.

## Contributing
`InfiniteOpt` is being actively developed and suggestions or other forms of contribution are encouraged.
There are many ways to contribute to this package. For more information please
visit [CONTRIBUTING](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/CONTRIBUTING.md).
