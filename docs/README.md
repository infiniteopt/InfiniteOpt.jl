# Documentation Builds
InfiniteOpt.jl's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). 
To install it, run the following command in a Julia session:

```julia
julia> ] 

(v1.6) pkg> add Documenter
```

You'll also need to make sure you have JuMP.jl, Distributions.jl, and Ipopt.jl 
installed for use in the guide's examples.

Once you have these packages installed you can build the documentation via:
```julia
julia --project=. --color=yes make.jl
```

in a Bash terminal. Or you can just run the `make.jl` file in a Julia REPL:
```julia-repl
julia> include("PUT/YOUR/PATH/HERE/InfiniteOpt/docs/make.jl") 
```

The compiled documents can then be viewed at `build/index.html` via a web browser 
of your choice (e.g., Chrome).
