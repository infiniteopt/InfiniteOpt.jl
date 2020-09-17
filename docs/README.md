# Documentation Builds
InfiniteOpt.jl's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). 
To install it, run the following command in a Julia session:

```julia
julia> ] 

(v1.5) pkg> add Documenter
```

Once you have Documenter installed you can build the documentation via:
```julia
julia --project=. --color=yes make.jl
```
The compiled documents can then be viewed at `build/index.html`.