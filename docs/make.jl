using Documenter
using CanonicalFlows

makedocs(
    sitename = "CanonicalFlows.jl",
    authors  = "Lachlan Whyborn and contributors",
    modules  = [CanonicalFlows],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://lachlanwhyborn.github.io/CanonicalFlows.jl",
    ),
    pages = [
        "Introduction" => "index.md",
        "Flows" => [
            "flows/boundary_layers.md",
        ],
        "Fluid Properties" => "fluid_properties.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)

deploydocs(
    repo   = "github.com/lachlanwhyborn/CanonicalFlows.jl",
    branch = "gh-pages",
)
