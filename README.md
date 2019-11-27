![Logo](full_logo.png)
---

A `JuMP` extension for expressing and solving infinite dimensional optimization
problems.

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev) | [![Build Status](https://travis-ci.com/pulsipher/InfiniteOpt.jl.svg?branch=master)](https://travis-ci.com/pulsipher/InfiniteOpt.jl) [![codecov.io](http://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=master)](http://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=master) |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev) | [![Build Status](https://travis-ci.com/pulsipher/InfiniteOpt.jl.svg?branch=master)](https://travis-ci.com/pulsipher/InfiniteOpt.jl) [![codecov.io](http://codecov.io/github/pulsipher/InfiniteOpt.jl/coverage.svg?branch=master)](http://codecov.io/github/pulsipher/InfiniteOpt.jl?branch=master) |

Its capabilities include:
- `JuMP`-like symbolic macro interface
- Infinite parameter support and parameterization of variables/constraints
- Finite parameters support and use (similar to `ParameterJuMP`)
- Direct support of infinite, point, and hold variables
- Symbolic measure (integral) expression
- Infinite/finite constraint definition
- Ordinary differential equation support (coming soon)
- Automated model transcription/reformulation and solution
- Readily extendable to accommodate user defined abstractions and solution techniques

Comments, suggestions and improvements are welcome and appreciated.

## License
`InfiniteOpt` is licensed under the [MIT "Expat" license](./LICENSE).

## Installation
`InfiniteOpt.jl` is still in its early stages of development but can be
installed by entering the following in the package manager.

```julia
(v1.3) pkg> add https://github.com/pulsipher/InfiniteOpt.jl
```

## Documentation
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://pulsipher.github.io/InfiniteOpt.jl/dev)

Please visit our [documentation pages](https://pulsipher.github.io/InfiniteOpt.jl/dev) to learn more. Note these are still under development and will be much more extensive soon.

## Project Status

The package is tested against Julia `1.0`, `1.1`, `1.2`, `1.3`, and nightly on Linux, macOS, and Windows.

## Contributing
`InfiniteOpt` is being actively developed and suggestions or other forms of contribution are encouraged.
There are many ways to contribute to this package:

- Suggest new/improved functionality
- Report an issue if you encounter some odd behavior, or if you have suggestions to improve the package.
- Contribute with code addressing some open issues, that add new functionality or that improve the performance.
- When contributing with code, add docstrings and comments, so others may understand the methods implemented.
- Contribute by updating and improving the documentation.
