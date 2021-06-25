![Logo](full_logo.png)
---

A `JuMP` extension for expressing and solving infinite-dimensional optimization
problems. Such areas include [stochastic programming](https://en.wikipedia.org/wiki/Stochastic_programming),
[dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming),
space-time optimization, and more. `InfiniteOpt` serves as an easy to use modeling
interface for these advanced problem types that can be used by those with little
to no background in these areas. It also it contains a wealth of capabilities
making it a powerful and convenient tool for advanced users.  

:warning: **`v0.4` introduced breaking changes**: See the documentation for details.

| **Documentation**                                                               | **Build Status**                                                                                | **Citation** |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:--------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/stable) | [![Build Status](https://github.com/pulsipher/InfiniteOpt.jl/workflows/CI/badge.svg?branch=release-0.4)](https://github.com/pulsipher/InfiniteOpt.jl/actions?query=workflow%3ACI) [![codecov.io](https://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=release-0.4)](https://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=release-0.4) | [![DOI](https://img.shields.io/badge/math.OC-arXiv%3A2106.12689-B31B1B.svg)](https://arxiv.org/abs/2106.12689) |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev) | [![Build Status](https://github.com/pulsipher/InfiniteOpt.jl/workflows/CI/badge.svg?branch=master)](https://github.com/pulsipher/InfiniteOpt.jl/actions?query=workflow%3ACI) [![codecov.io](https://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=master)](https://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=master) | |

It builds upon `JuMP` to add support for many complex modeling objects which 
include:
- Infinite parameters (e.g., time, space, uncertainty, etc.)
- Finite parameters (similar to `ParameterJuMP`)
- Infinite variables (e.g., `y(t, x)`)
- Derivatives (e.g., `∂y(t, x)/∂t`)
- Measures (e.g., `∫y(t,x)dt`, `𝔼[y(ξ)]`)
- More

The unifying modeling abstraction behind `InfiniteOpt` captures a wide spectrum 
of disciplines which include dynamic, PDE, stochastic, and semi-infinite 
optimization. Moreover, we facilitate transferring techniques between these 
to synthesize new optimization paradigms!

![abstract](abstraction.gif)

Comments, suggestions and improvements are welcome and appreciated.

## License
`InfiniteOpt` is licensed under the [MIT "Expat" license](./LICENSE).

## Installation
`InfiniteOpt.jl` is a registered package and can be installed by entering the
following in the package manager.

```julia
(v1.6) pkg> add InfiniteOpt
```

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/stable)

Please visit our [documentation pages](https://pulsipher.github.io/InfiniteOpt.jl/stable) 
to learn more. These pages are quite extensive and feature overviews, guides,
manuals, tutorials, examples, and more!

## Citing
[![DOI](https://img.shields.io/badge/math.OC-arXiv%3A2106.12689-B31B1B.svg)](https://arxiv.org/abs/2106.12689)

If you use InfiniteOpt.jl in your research, we would greatly appreciate your 
citing it.
```latex
@misc{pulsipher2021unifying,
      title={A Unifying Modeling Abstraction for Infinite-Dimensional Optimization}, 
      author={Joshua L. Pulsipher and Weiqi Zhang and Tyler J. Hongisto and Victor M. Zavala},
      year={2021},
      eprint={2106.12689},
      archivePrefix={arXiv},
      primaryClass={math.OC}
}
```

## Project Status
The package is tested against Julia `1.0` and `1.6` on Linux and Windows.

## Contributing
`InfiniteOpt` is being actively developed and suggestions or other forms of contribution are encouraged.
There are many ways to contribute to this package. For more information please
visit [CONTRIBUTING](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/CONTRIBUTING.md).
