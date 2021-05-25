using Documenter, InfiniteOpt, JuMP, Distributions

makedocs(;
    pages = ["Home" => "index.md",
             "Quick Start Guide" => "quick_start.md",
             "User Guide" => ["Installation" => "guide/install.md",
                              "Infinite Models" => "guide/model.md",
                              "Infinite Domains" => "guide/domains.md",
                              "Infinite Parameters" => "guide/parameter.md",
                              "Finite Parameters" => "guide/finite_parameter.md",
                              "Variables" => "guide/variable.md",
                              "Derivatives" => "guide/derivative.md",
                              "Expressions" => "guide/expression.md",
                              "Measures" => "guide/measure.md",
                              "Objectives" => "guide/objective.md",
                              "Constraints" => "guide/constraint.md",
                              "Model Transcription" => "guide/transcribe.md",
                              "Optimization" => "guide/optimize.md",
                              "Results" => "guide/result.md"
                              ]
                              ,
             "Examples" => "examples.md",
             "Extensions" => "extensions.md",
             "Development" => "develop.md",
             "Library" => "library.md",
             ],
    repo = "https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename = "InfiniteOpt.jl",
    authors = "Joshua Pulsipher and Weiqi Zhang",
    doctest = true,
    linkcheck = true,
    strict = false,
    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "UA-178297470-1",
        collapselevel = 1,
    )
)

deploydocs(;
    repo = "github.com/pulsipher/InfiniteOpt.jl",
    push_preview = true
)

