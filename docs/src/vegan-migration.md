# From R vegan

This guide maps common `vegan` diversity workflows to `DiversityAndDissimilarity.jl`.
The goal is not to duplicate every `vegan` feature, but to make the common
single-community and pairwise calculations easy to translate.

## Setup

In R:

```r
library(vegan)
```

In Julia:

```julia
using DiversityAndDissimilarity
```

Both packages commonly use a community matrix with samples/sites in rows and
taxa/species in columns.

## Alpha Diversity

| vegan | DiversityAndDissimilarity.jl | Notes |
|---|---|---|
| `specnumber(x)` | `richness(x)` | Observed positive-abundance taxa. |
| `diversity(x, index = "shannon")` | `entropy(Shannon(; base=ℯ), x)` | `vegan` defaults to natural logs; this package defaults to `base=2`. |
| `diversity(x, index = "shannon", base = 2)` | `shannon_entropy(x; base=2)` | Entropy in bits. |
| `exp(diversity(x, "shannon"))` | `effective_diversity(Shannon(; base=ℯ), x)` | Shannon effective diversity. |
| `diversity(x, index = "simpson")` | `gini_simpson_index(x)` | `vegan`'s Simpson output is `1 - sum(p_i^2)`. |
| `diversity(x, index = "invsimpson")` | `inverse_simpson_index(x)` | Inverse Simpson effective diversity. |
| `renyi(x, scales = q)` | `renyi_entropy(x, q; base=ℯ)` | `effective_diversity(Renyi(q; base=ℯ), x)` gives the Hill-number form. |
| `estimateR(x)` | `chao1(x)` / `ace(x)` | This package exposes Chao1 and ACE directly. |

```jldoctest
julia> using DiversityAndDissimilarity

julia> x = [1, 1, 2];

julia> entropy(Shannon(; base=ℯ), x)
1.0397207708399179

julia> gini_simpson_index(x)
0.625

julia> inverse_simpson_index(x)
2.6666666666666665
```

For exploratory work, [`alpha_diversity`](@ref) returns a compact named tuple
with observed richness, Shannon entropy and effective diversity, Simpson-family
indices, Chao1, ACE, and sample coverage.

```julia
alpha_diversity(x; base=ℯ)
```

For a community matrix, the same function returns one summary per row.

## Pairwise Distances

`vegan::vegdist` returns a compact R `dist` object. `DiversityAndDissimilarity.jl`
returns a dense symmetric matrix for community-matrix input.

| vegan | DiversityAndDissimilarity.jl | Notes |
|---|---|---|
| `vegdist(x, method = "bray")` | `bray_curtis_distance(x)` | Bray-Curtis dissimilarity. |
| `vegdist(x, method = "jaccard", binary = TRUE)` | `jaccard_distance(x)` | Presence/absence Jaccard. |
| `vegdist(x, method = "jaccard")` | `ruzicka_distance(x)` | Quantitative Jaccard/Ruzicka for abundance data. |
| `vegdist(decostand(x, "hellinger"), method = "euclidean")` | `hellinger_distance(x)` | Hellinger distance on normalized abundances. |
| `vegdist(decostand(x, "normalize"), method = "euclidean")` | `chord_distance(x)` | Chord-style distance. |
| `vegdist(x, method = "canberra")` | `canberra_distance(x)` | This package uses the averaged positive-denominator form documented in its API. |

```jldoctest
julia> using DiversityAndDissimilarity

julia> community = [
           1 1 2 0 5
           3 0 1 1 0
       ];

julia> bray_curtis_distance(community)
2×2 Matrix{Float64}:
 0.0       0.714286
 0.714286  0.0
```

Use [`similarity`](@ref), [`dissimilarity`](@ref), and [`distance`](@ref) with
index objects when you prefer explicit dispatch:

```julia
distance(BrayCurtis(), community)
similarity(Jaccard(), community)
```

## Table Inputs

`vegan` users often keep site metadata beside species columns. In Julia,
Tables.jl-compatible inputs can be converted with [`community_matrix`](@ref),
and most high-level functions accept a `species` column list.

```julia
community_matrix(table; species=[:oak, :ash, :elm])
alpha_diversity(table; species=[:oak, :ash, :elm])
bray_curtis_distance(table; species=[:oak, :ash, :elm])
```

Passing `species` explicitly is recommended when a table contains numeric
metadata such as site IDs, blocks, or coordinates.

## Convention Differences

- Shannon entropy defaults to `base=2` in `DiversityAndDissimilarity.jl`; use `base=ℯ`
  to match `vegan`'s default natural-log convention.
- [`Simpson`](@ref) in this package is the concentration
  ``\sum_i p_i^2``. Use [`GiniSimpson`](@ref) or
  [`gini_simpson_index`](@ref) to match `vegan`'s `index = "simpson"`.
- Pairwise community-matrix functions return dense matrices. For very large
  numbers of sites, memory grows quadratically.
- Broader ordination, constrained analysis, rarefaction curves, and many
  ecological modeling tools remain `vegan` strengths. `DiversityAndDissimilarity.jl`
  focuses on lightweight diversity, entropy, richness, and pairwise
  similarity/dissimilarity calculations.

See also the [reference examples](reference-examples.md), which pin several
values to `vegan` and scikit-bio conventions.
