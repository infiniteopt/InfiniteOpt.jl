using Documenter, InfiniteOpt

makedocs(;
    modules=[InfiniteOpt],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename="InfiniteOpt.jl",
    authors="Joshua Pulsipher",
    assets=String[],
)

deploydocs(;
    repo="github.com/pulsipher/InfiniteOpt.jl",
)
