# DiversityAndDissimilarity.jl

`DiversityAndDissimilarity` provides ecology-style diversity indices for species
abundance data, observation vectors, and pairwise assemblage comparisons.

The package is built around small dispatchable index types and generic
operations:

- [`entropy`](@ref) for Shannon, Renyi, and Tsallis entropy values.
- [`diversity`](@ref) for alpha diversity and effective diversity.
- [`effective_diversity`](@ref) for effective species numbers.
- [`similarity`](@ref), [`dissimilarity`](@ref), and [`distance`](@ref) for
  pairwise comparisons.
- [`counts`](@ref) and [`proportions`](@ref) for preparing abundance data.
- [`alpha_diversity`](@ref) for a compact exploratory alpha-diversity summary.
- [`index_metadata`](@ref), [`index_bounds`](@ref), descriptor helpers such as
  [`is_metric`](@ref) and [`is_similarity`](@ref), [`reference_cases`](@ref),
  and [`diversity_audit`](@ref) for convention-aware and validation-oriented
  workflows.

If you are coming from R, start with the
[vegan migration guide](vegan-migration.md). If you are choosing between
packages, the [diversity](index-checklist.md) and
[similarity/dissimilarity](similarity-checklist.md) checklists summarize where
`DiversityAndDissimilarity.jl` overlaps with common ecology and biodiversity tools.

```@contents
Pages = [
    "getting-started.md",
    "examples.md",
    "vegan-migration.md",
    "indices.md",
    "diversity-indices.md",
    "choosing-an-estimator.md",
    "similarity-indices.md",
    "similarity-index-catalog.md",
    "scaling.md",
    "benchmarks.md",
    "framework.md",
    "case-study.md",
    "reference-examples.md",
    "release-readiness.md",
    "checklists.md",
    "index-checklist.md",
    "similarity-checklist.md",
    "api.md",
]
Depth = 2
```

## Supported Julia Versions

`DiversityAndDissimilarity` supports Julia 1.10 and newer 1.x releases.

## Quality Checks

The standard test suite runs across the supported Julia compatibility matrix.
A separate GitHub Actions quality workflow runs Aqua.jl and JET.jl checks from
the `test/quality` environment.

## Index Types

Alpha diversity:

- [`Richness`](@ref)
- [`Shannon`](@ref), which defaults to `base=2` for entropy in bits
- [`Renyi`](@ref)
- [`Tsallis`](@ref)
- [`Simpson`](@ref)
- [`GiniSimpson`](@ref)
- [`InverseSimpson`](@ref)
- [`Hill`](@ref)
- [`Chao1`](@ref)
- [`ACE`](@ref)
- [`SampleCoverage`](@ref)
- [`PielouEvenness`](@ref)
- [`FisherAlpha`](@ref)

Pairwise comparisons:

- [`Jaccard`](@ref)
- [`SorensenDice`](@ref)
- [`BrayCurtis`](@ref)
