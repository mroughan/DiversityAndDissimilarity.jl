using Documenter
using DiversityAndDissimilarity

DocMeta.setdocmeta!(
    DiversityAndDissimilarity,
    :DocTestSetup,
    :(using DiversityAndDissimilarity);
    recursive=true,
)

makedocs(;
    modules=[DiversityAndDissimilarity],
    sitename="DiversityAndDissimilarity.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mroughan.github.io/DiversityAndDissimilarity.jl",
        repolink="https://github.com/mroughan/DiversityAndDissimilarity.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Examples" => "examples.md",
        "From R vegan" => "vegan-migration.md",
        "Indices" => [
            "Indices" => "indices.md",
            "Diversity Indices" => "diversity-indices.md",
            "Choosing An Estimator" => "choosing-an-estimator.md",
            "Similarity And Dissimilarity Indices" => "similarity-indices.md",
            "Similarity And Dissimilarity Catalog" => "similarity-index-catalog.md",
            "Scaling And Performance" => "scaling.md",
            "Benchmarks" => "benchmarks.md",
        ],
        "Framework" => [
            "Framework" => "framework.md",
            "Type Structure" => "type-structure.md",
            "Workflow Case Study" => "case-study.md",
            "Reference Examples" => "reference-examples.md",
            "Release Readiness" => "release-readiness.md",
        ],
        "Availability Checklists" => [
            "Availability Checklists" => "checklists.md",
            "Diversity Indices" => "index-checklist.md",
            "Similarity And Dissimilarity" => "similarity-checklist.md",
        ],
        "API Reference" => "api.md",
    ],
    checkdocs=:none,
)

github_ref = get(ENV, "GITHUB_REF", "")
should_deploy_docs =
    get(ENV, "CI", "false") == "true" &&
    get(ENV, "GITHUB_EVENT_NAME", "") == "push" &&
    (github_ref in ("refs/heads/main", "refs/heads/master") ||
     startswith(github_ref, "refs/tags/"))

if should_deploy_docs
    deploydocs(;
        repo="github.com/mroughan/DiversityAndDissimilarity.jl",
    )
end
