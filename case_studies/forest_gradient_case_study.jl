using DiversityAndDissimilarity
using Random
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DATA_FILE = joinpath(ROOT, "validation", "data", "community_counts.csv")
const REPORT_FILE = joinpath(@__DIR__, "forest_gradient_case_study.md")
const TEX_FILE = joinpath(ROOT, "notes", "case-study-results.tex")

function read_community_counts(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    taxon_start = findfirst(==("taxon_01"), header)
    labels = String[]
    habitats = String[]
    gradients = String[]
    rows = Vector{Vector{Int}}()
    for line in lines[2:end]
        parts = split(line, ",")
        push!(labels, parts[1])
        push!(habitats, parts[2])
        push!(gradients, parts[3])
        push!(rows, parse.(Int, parts[taxon_start:end]))
    end
    matrix = reduce(vcat, permutedims.(rows))
    taxa = header[taxon_start:end]
    return (; labels, habitats, gradients, taxa, matrix)
end

function fmt(x; digits=3)
    return string(round(float(x); digits))
end

function interval_text(result)
    return "$(fmt(result.estimate)) [$(fmt(result.lower)), $(fmt(result.upper))]"
end

function quantity_result(report, quantity)
    return only(filter(result -> result.quantity == quantity, report.estimates))
end

function write_markdown(path, data, audit, uncertainty)
    open(path, "w") do io
        println(io, "# Forest Gradient Workflow Case Study")
        println(io)
        println(io, "This simulated case study uses `validation/data/community_counts.csv` as a")
        println(io, "transparent stand-in for a community monitoring workflow. The point is to")
        println(io, "exercise convention-sensitive biodiversity calculations, not to infer from")
        println(io, "real ecological observations.")
        println(io)
        println(io, "## Workflow")
        println(io)
        println(io, "1. Read a sample-by-taxon matrix with site labels and habitat metadata.")
        println(io, "2. Run `diversity_audit` for alpha summaries, estimator diagnostics, and a labeled Bray-Curtis matrix.")
        println(io, "3. Run `uncertainty_audit` to bootstrap Shannon entropy and effective diversity for each site.")
        println(io, "4. Inspect low-coverage warnings before using the pairwise matrix in downstream analyses.")
        println(io)
        println(io, "## Dataset")
        println(io)
        println(io, "- Samples: $(length(data.labels))")
        println(io, "- Taxa: $(length(data.taxa))")
        println(io, "- Habitats: $(join(unique(data.habitats), ", "))")
        println(io, "- Total abundance range: $(minimum(audit.row_totals)) to $(maximum(audit.row_totals))")
        println(io)
        println(io, "## Alpha And Uncertainty Summary")
        println(io)
        println(io, "| Site | Habitat | Gradient | Richness | Coverage | Shannon H bootstrap | Effective diversity bootstrap |")
        println(io, "|---|---|---|---:|---:|---:|---:|")
        for i in eachindex(data.labels)
            alpha = audit.alpha[i]
            report = uncertainty.reports[i]
            entropy = quantity_result(report, :entropy)
            diversity = quantity_result(report, :diversity)
            println(io, "| $(data.labels[i]) | $(data.habitats[i]) | $(data.gradients[i]) | $(alpha.richness) | $(fmt(alpha.sample_coverage)) | $(interval_text(entropy)) | $(interval_text(diversity)) |")
        end
        println(io)
        println(io, "## Audit Warnings")
        println(io)
        if isempty(audit.warnings) && isempty(uncertainty.warnings)
            println(io, "No audit warnings were produced.")
        else
            for warning in audit.warnings
                println(io, "- diversity audit: $warning")
            end
            for warning in uncertainty.warnings
                println(io, "- uncertainty audit: $warning")
            end
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "The workflow keeps alpha summaries, coverage diagnostics, bootstrap")
        println(io, "intervals, and the Bray-Curtis distance matrix tied to the same labels and")
        println(io, "input matrix. This is the core reproducibility benefit: the report exposes")
        println(io, "which samples are convention-sensitive or coverage-sensitive before the")
        println(io, "distance matrix is reused for ordination, clustering, or modelling.")
    end
end

function write_latex(path, data, audit, uncertainty)
    open(path, "w") do io
        println(io, "\\begin{table}[t]")
        println(io, "\\centering")
        println(io, "\\caption{Selected workflow case-study audit results. Shannon values are bootstrap estimates with percentile intervals from \\texttt{uncertainty\\_audit}.}")
        println(io, "\\label{tab:case-study}")
        println(io, "\\small")
        println(io, "\\resizebox{\\textwidth}{!}{%")
        println(io, "\\begin{tabular}{llrlll}")
        println(io, "\\toprule")
        println(io, "Site & Habitat & Richness & Coverage & Shannon H & Effective D\\\\")
        println(io, "\\midrule")
        for i in eachindex(data.labels)
            alpha = audit.alpha[i]
            report = uncertainty.reports[i]
            entropy = quantity_result(report, :entropy)
            diversity = quantity_result(report, :diversity)
            println(io, "$(replace(data.labels[i], "_" => "\\_")) & $(data.habitats[i]) & $(alpha.richness) & $(fmt(alpha.sample_coverage)) & $(interval_text(entropy)) & $(interval_text(diversity))\\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "}")
        println(io, "\\end{table}")
    end
end

function main()
    data = read_community_counts(DATA_FILE)
    audit = diversity_audit(data.matrix; labels=data.labels, pairwise_index=BrayCurtis())
    uncertainty = uncertainty_audit(
        data.matrix;
        labels=data.labels,
        nboot=400,
        level=0.95,
        rng=MersenneTwister(20260524),
    )
    write_markdown(REPORT_FILE, data, audit, uncertainty)
    write_latex(TEX_FILE, data, audit, uncertainty)
    println("Wrote $REPORT_FILE")
    println("Wrote $TEX_FILE")
end

main()
