using Documenter, InfiniteOpt, JuMP

makedocs(;
    modules=[InfiniteOpt],
    format = Documenter.HTML(assets = ["assets/custom.css"]),
    pages=[
        "Home" => "index.md",
        "Library" => "library.md",
    ],
    repo="https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename="InfiniteOpt.jl",
    authors="Joshua Pulsipher"
)

deploydocs(;
    repo="github.com/pulsipher/InfiniteOpt.jl",
)
