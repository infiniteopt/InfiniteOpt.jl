using Documenter, InfiniteOpt, Distributions, Literate, Random

if !@isdefined(EXAMPLE_DIR)
    const EXAMPLE_DIR = joinpath(@__DIR__, "src", "examples")
end

# Sort and collect all the files of a certain extension in a folder
sorted_files(dir, ext) = sort!(filter!(f -> endswith(f, ext), readdir(dir)))

# Remove any files generated by Literate that are still in a subdirectory
function remove_literate_files(dir)
    for file in append!(sorted_files(dir, ".md"), sorted_files(dir, ".png"))
        rm(joinpath(dir, file))
    end
    return
end

# Post process corrections
function format_fix(content)
    return replace(content, "nothing #hide" => "")
end

# Generate all the markdown files via Literate for a particular doc section
# Expect the directory to only have sub-directories that contain files
function generate_pages(dir)
    for subdir in sort!(readdir(dir))
        current_dir = joinpath(dir, subdir)
        remove_literate_files(current_dir)
        for file in sorted_files(current_dir, ".jl")
            jl_filename = joinpath(current_dir, file)
            Random.seed!(0)
            Literate.markdown(jl_filename, current_dir; documenter = true, 
                              execute = true, postprocess = format_fix)
        end
    end
    return
end

# Generate the list of subsection pages with relative paths
function subsection_entry(dir, subdir)
    files = sorted_files(joinpath(dir, subdir), ".md")
    return map(file -> joinpath(basename(dir), subdir, file), files)
end

# Generate the list of subsections and their pages with relative paths
function section_entry(dir)
    return map(subdir -> subdir => subsection_entry(dir, subdir), readdir(dir))
end

# Run and convert the examples
generate_pages(EXAMPLE_DIR)

# Make the documentation via Documenter.jl
makedocs(;
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Tutorials" => ["Quick Start" => "tutorials/quick_start.md"], # TODO add more tutorials
        "Examples" =>  section_entry(EXAMPLE_DIR),
        # "Background" => [], # TODO add the background
        "User Guide" => [
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
    repo = "https://github.com/infiniteopt/InfiniteOpt.jl/blob/{commit}{path}#L{line}",
    sitename = "InfiniteOpt.jl",
    authors = "Joshua Pulsipher and Weiqi Zhang",
    doctest = true,
    checkdocs = :exports,
    linkcheck = true,
    linkcheck_ignore = [r"https://www.youtube.com/.*"],
    strict = true,
    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "UA-178297470-1",
        collapselevel = 1,
        assets = ["assets/extra_styles.css"],

    )
)

# Push the built HTML files to the correct GitHub branch if building via CI
deploydocs(;
    repo = "github.com/infiniteopt/InfiniteOpt.jl",
    push_preview = true
)
