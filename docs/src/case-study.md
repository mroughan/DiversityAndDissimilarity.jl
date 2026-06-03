# Workflow Case Study

This case study uses the simulated validation matrix in
`validation/data/community_counts.csv`. The dataset is not intended to
represent a real survey; it is a transparent workflow fixture with uneven
abundances, singletons, zero-heavy taxa, and habitat labels.

The script is in `case_studies/forest_gradient_case_study.jl` and produces:

- `case_studies/forest_gradient_case_study.md`
- `notes/case-study-results.tex`

The workflow demonstrates how the package keeps conventions, labels, coverage,
pairwise matrices, and uncertainty diagnostics together.

```julia
using DiversityIndices
using Random

community = [
    18 7 0 3 1 0 0 12 2 0 0 4
    12 5 1 0 0 2 0 9 1 1 0 3
    5 8 4 2 0 0 1 3 0 2 1 0
    0 3 9 6 2 1 0 0 0 4 0 1
]

labels = ["plot_A", "plot_B", "plot_C", "plot_D"]

audit = diversity_audit(community; labels, pairwise_index=BrayCurtis())

uncertainty = uncertainty_audit(
    community;
    labels,
    nboot=400,
    rng=MersenneTwister(20260524),
)
```

The audit report contains alpha-diversity summaries, estimator reports, a
labeled Bray-Curtis matrix, and warnings:

```julia
audit.alpha
audit.estimator_report
audit.pairwise
audit.warnings
```

The uncertainty audit reports Shannon entropy and effective Shannon diversity
intervals per sample:

```julia
uncertainty.reports
```

This is useful as an early workflow check: low coverage, unstable estimates, or
unexpected labels can be addressed before the pairwise matrix is reused in
ordination, clustering, or modelling.
