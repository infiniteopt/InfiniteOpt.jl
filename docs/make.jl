using Documenter, InfiniteOpt, JuMP, Distributions

makedocs(;
    pages = ["Home" => "index.md",
            "User Guide" => ["Installation" => "guide/install.md",
                             "Infinite Models" => "guide/model.md",
                             "Infinite Parameters" => "guide/parameter.md",
                             "Finite Parameters" => "guide/finite_parameter.md",
                             "Variables" => "guide/variable.md",
                             "Expressions" => "guide/expression.md",
                             "Measures" => "guide/measure.md",
                             "Objectives" => "guide/objective.md",
                             "Constraints" => "guide/constraint.md",
                             "Model Transcription" => "guide/transcribe.md",
                             "Optimization" => "guide/optimize.md",
                             "Results" => "guide/result.md"],
            "Examples" => "examples.md",
            "Extensions" => "extensions.md",
            "Development" => "develop.md",
            "Library" => "library.md",
            hide("JuMP Docs" => "JuMP.md")],
    repo = "https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename = "InfiniteOpt.jl",
    authors = "Joshua Pulsipher and Weiqi Zhang",
    doctest = false,
    linkcheck = true
)

deploydocs(;
    repo = "github.com/pulsipher/InfiniteOpt.jl",
)
