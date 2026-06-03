using Documenter
using DiversityIndices

DocMeta.setdocmeta!(
    DiversityIndices,
    :DocTestSetup,
    :(using DiversityIndices);
    recursive=true,
)

makedocs(;
    modules=[DiversityIndices],
    sitename="DiversityIndices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mroughan.github.io/DiversityIndices.jl",
        repolink="https://github.com/mroughan/DiversityIndices.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Examples" => "examples.md",
        "From R vegan" => "vegan-migration.md",
        "Framework" => "framework.md",
        "Workflow Case Study" => "case-study.md",
        "Indices" => [
            "Overview" => "indices.md",
            "Diversity Indices" => "diversity-indices.md",
            "Choosing An Estimator" => "choosing-an-estimator.md",
            "Similarity And Dissimilarity Indices" => "similarity-indices.md",
            "Similarity And Dissimilarity Catalog" => "similarity-index-catalog.md",
            "Scaling And Performance" => "scaling.md",
            "Benchmarks" => "benchmarks.md",
            "Reference Examples" => "reference-examples.md",
            "Release Readiness" => "release-readiness.md",
        ],
        "Availability Checklists" => [
            "Diversity Indices" => "index-checklist.md",
            "Similarity And Dissimilarity" => "similarity-checklist.md",
        ],
        "API Reference" => "api.md",
    ],
    checkdocs=:none,
)

if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/mroughan/DiversityIndices.jl",
    )
end
