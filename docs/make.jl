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
        "Introduction" => "index.md",
        "Data Input Formats" => "data-input.md",
        "Diversity Indices" => "diversity-indices.md",
        "Dissimilarity Indices" => "similarity-indices.md",
        "Framework" => "framework.md",
        "Additional Information" => "additional-information.md",
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
