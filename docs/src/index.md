# Introduction

`DiversityAndDissimilarity.jl` provides lightweight Julia tools for ecological
diversity, entropy, richness estimation, and pairwise assemblage comparison.
It is designed for count or abundance data where rows are usually samples,
sites, plots, languages, or other assemblages, and columns are taxa or
categories.

The package is built around small index objects:

```julia
Richness()
Shannon()
BrayCurtis()
JensenShannon()
```

Those objects are passed to generic operations:

```julia
diversity(Richness(), assemblage)
entropy(Shannon(), assemblage)
dissimilarity(BrayCurtis(), left, right)
distance(JensenShannon(), left, right)
```

This makes the calculation explicit while keeping workflows compact. It also
lets the package attach metadata to every index: whether it is a similarity or
dissimilarity, whether it is bounded, whether it is a metric, what its bounds
mean, and which conventions differ from other packages.

## Why Use This Package?

Use `DiversityAndDissimilarity.jl` when you want:

- alpha-diversity summaries for one assemblage or every row of a community
  matrix;
- pairwise similarity, dissimilarity, distance, and divergence calculations;
- entropy estimators and low-count corrections with explicit estimator
  objects;
- convention-aware metadata for reports, validation, and generic workflows;
- small, inspectable implementations that can be checked against published and
  cross-package reference values.

The package is intentionally modest in scope. It aims to make core diversity
and dissimilarity calculations reliable and transparent rather than becoming a
full ecological modelling environment.

## Quickstart

```jldoctest intro
julia> using DiversityAndDissimilarity

julia> assemblage = Dict(:oak => 12, :ash => 5, :elm => 3);

julia> richness(assemblage)
3

julia> shannon_entropy(assemblage)
1.3527241956246545

julia> shannon_diversity(assemblage)
2.5539392274300625

julia> simpson_index(assemblage)
0.445

julia> gini_simpson_index(assemblage)
0.5549999999999999
```

For a community matrix, rows are samples and columns are taxa:

```jldoctest intro
julia> community = [
           1 1 2 0 5
           3 0 1 1 0
       ];

julia> richness(community)
2-element Vector{Int64}:
 4
 3

julia> bray_curtis_distance(community)
2×2 Matrix{Float64}:
 0.0       0.714286
 0.714286  0.0
```

For pairwise comparisons, use an explicit index object:

```jldoctest intro
julia> left = Dict(:oak => 12, :ash => 5);

julia> right = Dict(:ash => 4, :elm => 7);

julia> similarity(Jaccard(), left, right)
0.3333333333333333

julia> dissimilarity(BrayCurtis(), left, right)
0.7142857142857143
```

The same object can be inspected before using it:

```jldoctest intro
julia> is_metric(JensenShannon())
true

julia> is_symmetric(KullbackLeibler())
false

julia> index_bounds(BrayCurtis()).lower_meaning
"minimal dissimilarity; identical or indistinguishable inputs"
```

## Goals

- Keep common diversity and dissimilarity calculations easy to call and easy
  to inspect.
- Make convention choices explicit, especially for names such as Simpson,
  Shannon diversity, Jaccard, Bray-Curtis, and Jensen-Shannon.
- Support common data shapes: vectors, dictionaries, observation vectors,
  matrices, and Tables.jl-compatible tables.
- Expose metadata that software can use, not only prose explanations.
- Maintain reference examples and validation checks against other ecosystems
  such as vegan, scikit-bio, SciPy, and iNEXT.

## Non-Goals

- The package is not an ordination, regression, rarefaction-curve, or null-model
  framework.
- It does not try to reimplement every index from every biodiversity package.
- It avoids a string-based `method="..."` API in favor of dispatchable index
  objects.
- It avoids heavy runtime dependencies for plotting, external package bridges,
  documentation generation, and manuscript tooling.

## Where To Go Next

- [Diversity Indices](diversity-indices.md) covers alpha-diversity, entropy,
  estimators, uncertainty, and the diversity availability checklist.
- [Dissimilarity Indices](similarity-indices.md) covers pairwise incidence,
  abundance, probability, and information-theoretic comparisons.
- [Framework](framework.md) describes the type structure, data semantics,
  metadata traits, architecture, goals, and non-goals in more detail.
- [Additional Information](additional-information.md) collects estimator
  guidance, validation examples, scaling notes, migration notes, and generated
  documentation assets.
- [API Reference](api.md) is the compact function and type index.

```@contents
Pages = [
    "diversity-indices.md",
    "similarity-indices.md",
    "framework.md",
    "additional-information.md",
    "api.md",
]
Depth = 2
```
