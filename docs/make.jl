using Documenter, InfiniteOpt, JuMP

makedocs(;
    modules = [InfiniteOpt, InfiniteOpt.TranscriptionOpt, JuMP],
    format = Documenter.HTML(assets = ["assets/custom.css"]),
    pages = ["Home" => "index.md",
            "User Guide" => ["Installation" => "guide/install.md",
                             "Infinite Models" => "guide/model.md",
                             "Infinite Parameters" => "guide/parameter.md",
                             "Variables" => "guide/variable.md",
                             "Expressions" => "guide/expression.md",
                             "Measures" => "guide/measure.md",
                             "Objectives" => "guide/objective.md",
                             "Constraints" => "guide/constraint.md",
                             "Model Transcription" => "guide/transcribe.md",
                             "Optimization" => "guide/optimize.md",
                             "Results" => "guide/result.md"],
            "Examples" => "examples.md",
            "Library" => "library.md", # TODO breakdown into subpages
            "Extensions" => "extensions.md",
            "Development" => "develop.md"],
    repo = "https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename = "InfiniteOpt.jl",
    authors = "Joshua Pulsipher",
    doctest = false
)

deploydocs(;
    repo = "github.com/pulsipher/InfiniteOpt.jl",
)
