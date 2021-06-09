using Documenter, InfiniteOpt, Distributions, Test, Literate

if !@isdefined(EXAMPLE_DIR)
    const EXAMPLE_DIR = joinpath(@__DIR__, "src", "examples")
    const EXAMPLE_SUBDIR = ["Optimal Control", "Stochastic Optimization"]
end

function link_example(content)
    edit_url = match(r"EditURL = \"(.+?)\"", content)[1]
    footer = match(r"^(---\n\n\*This page was generated using)"m, content)[1]
    content = replace(
        content, footer => "[View this file on Github]($(edit_url)).\n\n" * footer
    )
    return content
end

function file_list(full_dir, relative_dir, extension)
    return map(
        file -> joinpath(relative_dir, file),
        filter(file -> endswith(file, extension), sort(readdir(full_dir))),
    )
end

# Include the `filename` in a temporary module that acts as a sandbox. (Ensuring
# no constants or functions leak into other files.)
function include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

function literate_directory(dir)
    rm.(file_list(dir, dir, ".md"))
    for filename in file_list(dir, dir, ".jl")
        # `include` the file to test it before `#src` lines are removed. It is
        # in a testset to isolate local variables between files.
        @testset "$(filename)" begin
            include_sandbox(filename)
        end
        Literate.markdown(
            filename,
            dir;
            documenter = true,
            postprocess = link_example,
        )
    end
    return
end

literate_directory.(joinpath.(EXAMPLE_DIR, EXAMPLE_SUBDIR))

makedocs(;
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Tutorials" => ["Quick Start" => "tutorials/quick_start.md"], # TODO add more tutorials
        "Examples" =>  map(
            subdir -> subdir => map(
                file -> joinpath("examples", subdir, file),
                filter(
                    file -> endswith(file, ".md"),
                    sort(readdir(joinpath(EXAMPLE_DIR, subdir))),
                ),
            ),
            EXAMPLE_SUBDIR,
        ),
        # "Background" => [], # TODO add the background
        "Guide" => [
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
            ],
        "API Manual" => [
            "Infinite Models" => "manual/model.md",
            "Infinite Domains" => "manual/domains.md",
            "Infinite Parameters" => "manual/parameter.md",
            "Finite Parameters" => "manual/finite_parameter.md",
            "Variables" => "manual/variable.md",
            "Derivatives" => "manual/derivative.md",
            "Expressions" => "manual/expression.md",
            "Measures" => "manual/measure.md",
            "Objectives" => "manual/objective.md",
            "Constraints" => "manual/constraint.md",
            "Model Transcription" => "manual/transcribe.md",
            "Optimization" => "manual/optimize.md",
            "Results" => "manual/result.md"
            ],
        "Development" => [
            "Extensions" => "develop/extensions.md",
            "Getting Started" => "develop/start_guide.md",
            "Style Guide" => "develop/style.md"
            ],
        ],
    repo = "https://github.com/pulsipher/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename = "InfiniteOpt.jl",
    authors = "Joshua Pulsipher and Weiqi Zhang",
    doctest = false,
    # checkdocs = :exports,
    linkcheck = false,
    linkcheck_ignore = [r"https://www.youtube.com/.*"],
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

