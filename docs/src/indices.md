# Index Overview

`DiversityAndDissimilarity` separates its index types into two related groups:

- Diversity indices summarize one assemblage at a time.
- Similarity and dissimilarity indices compare two assemblages.

Both groups share the same small-dispatch style: choose an index object, then
pass it to a generic operation.

```julia
entropy(Shannon(), assemblage)
diversity(Shannon(), assemblage)
similarity(Jaccard(), left, right)
dissimilarity(BrayCurtis(), left, right)
```

## Relationship Between Operations

Alpha-diversity functions answer "how diverse is this assemblage?"

```math
\\mathrm{data} \\rightarrow \\mathrm{diversity\\ value}
```

Pairwise comparison functions answer "how similar or different are these two
assemblages?"

```math
(\\mathrm{left}, \\mathrm{right}) \\rightarrow \\mathrm{similarity\\ or\\ dissimilarity\\ value}
```

Some concepts line up across the two groups. Incidence-based pairwise indices,
such as [`Jaccard`](@ref) and [`SorensenDice`](@ref), only use species presence
or absence. Abundance-based indices, such as [`BrayCurtis`](@ref), use the
amount associated with each species. The alpha-diversity indices similarly
depend on proportions, ``p_i``, derived from positive abundances.

The same proportion-based ideas also appear outside ecology. For example,
[`LinguisticDiversityIndex`](@ref) and [`GreenbergDiversityIndex`](@ref) are
linguistic-demography names for the same quantity as [`GiniSimpson`](@ref):
the probability that two randomly selected individuals belong to different
categories, interpreted as different mother tongues.

## Pages

- [Diversity Indices](diversity-indices.md) covers [`entropy`](@ref),
  [`diversity`](@ref), [`effective_diversity`](@ref), entropy families,
  Simpson-family indices, and convenience functions such as
  [`shannon_entropy`](@ref).
- [Similarity And Dissimilarity Indices](similarity-indices.md) covers [`similarity`](@ref),
  [`dissimilarity`](@ref), incidence comparisons, abundance comparisons, and
  convenience functions such as [`jaccard_index`](@ref).
