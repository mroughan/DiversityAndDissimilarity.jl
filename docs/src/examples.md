# Examples

These examples show common workflows without introducing extra package
dependencies in the core library.

## Alpha Summary

```julia
using DiversityIndices

assemblage = Dict(:oak => 12, :ash => 5, :elm => 3)

alpha_diversity(assemblage)
pielou_evenness(assemblage)
fisher_alpha(assemblage)
```

For a community matrix, summaries are returned row-wise:

```julia
community = [
    1 1 2 0 5
    3 0 1 1 0
]

alpha_diversity(community)
```

## Labeled Pairwise Matrices

When rows have sample/site names, pass labels explicitly:

```julia
sites = ["plot-a", "plot-b"]
labeled_distance(BrayCurtis(), community; labels=sites)
```

For Tables.jl-compatible inputs, pass the label column and the species columns:

```julia
labeled_distance(
    BrayCurtis(),
    table;
    label=:site,
    species=[:oak, :ash, :elm],
)
```

The result is a named tuple:

```julia
result = labeled_distance(BrayCurtis(), community; labels=sites)
result.labels
result.matrix
```

## Vegan-Style Calls

Use natural logarithms to match `vegan`'s default Shannon convention:

```julia
entropy(Shannon(; base=ℯ), assemblage)
```

Use Gini-Simpson to match `vegan`'s `index = "simpson"` output:

```julia
gini_simpson_index(assemblage)
```

Use Bray-Curtis for the common `vegdist(method = "bray")` workflow:

```julia
bray_curtis_distance(community)
```
