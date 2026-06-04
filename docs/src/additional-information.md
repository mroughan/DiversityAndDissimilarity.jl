# Additional Information

This page collects supporting material that does not belong only to diversity
indices, dissimilarity indices, the framework, or the API reference.

## Choosing Shannon Estimators

Shannon entropy estimation is sensitive to finite sample size, unseen support,
and whether the category universe is known. The package uses estimator objects
so the choice is visible in the call.

![Flowchart for Shannon entropy estimator choice](assets/flow_entropy_estimation.svg)

A practical starting point:

- Use [`Plugin`](@ref) when sample size and coverage are high.
- Compare [`MillerMadow`](@ref) when an observed-support first-order bias
  correction is appropriate.
- Use [`HausserStrimmer`](@ref), [`Basharin`](@ref), or [`AddGamma`](@ref) when
  the finite support is known.
- Use [`ChaoShen`](@ref) when support is unknown and unseen categories are
  plausible.

```jldoctest additional
julia> using DiversityAndDissimilarity

julia> assemblage = Dict(:oak => 12, :ash => 5, :elm => 3);

julia> entropy(Shannon(; estimator=Plugin()), assemblage)
1.3527241956246545

julia> entropy(Shannon(; estimator=MillerMadow()), assemblage)
1.4248589476691027

julia> entropy(Shannon(; estimator=AddGamma(1)), assemblage; support=5)
1.7792365361682794
```

Use [`estimator_report`](@ref) to compare estimators and coverage diagnostics
on one assemblage:

```jldoctest additional
julia> report = estimator_report([1, 1, 2, 0, 5]);

julia> report.observed_richness
4

julia> report.sample_coverage
0.7777777777777778
```

## Reference Values

Reference examples are part of the test suite. They pin conventions to simple
published or cross-package values.

```jldoctest additional
julia> x = [1, 1, 2];

julia> entropy(Shannon(; base=ℯ), x)
1.0397207708399179

julia> gini_simpson_index(x)
0.625

julia> inverse_simpson_index(x)
2.6666666666666665
```

For pairwise comparisons:

```jldoctest additional
julia> bray_curtis_distance([1, 2, 3], [2, 2, 0])
0.4

julia> jaccard_similarity([1, 1, 0, 1], [1, 0, 1, 1])
0.5
```

The public validation API exposes a small reference corpus:

```jldoctest additional
julia> results = validate_reference_cases();

julia> all(r -> r.passed, results)
true

julia> first(results).name
"vegan_shannon_natural_log"
```

The `validation` directory also records candidate external datasets such as
the scikit-bio tutorial matrix, vegan's `dune`, `varespec`, `BCI`, and `mite`,
and iNEXT's `spider` data.

## Scaling Notes

Most single-assemblage calculations are linear in the number of inspected
categories. The most important scaling variables are:

- `S`: observed support, the number of positive taxa/categories;
- `K`: known finite support supplied through `support`;
- `M`: number of samples/sites, usually rows in a community matrix;
- `P`: number of taxa/categories columns;
- `B`: number of bootstrap replicates.

Rules of thumb:

- alpha-diversity on one abundance vector is usually `O(S)` or `O(P)`;
- row-wise alpha diversity on a dense community matrix is `O(MP)`;
- pairwise community matrices are usually `O(M^2 P)` time and `O(M^2)` memory;
- finite-support estimators can scale with `K`, not only observed `S`;
- bootstrap cost grows roughly linearly with `nboot`.

For large pairwise workflows, prefer bounded measures with interpretable ranges
when that suits the scientific question:

```jldoctest additional
julia> is_bounded(BrayCurtis())
true

julia> is_bounded(KullbackLeibler())
false
```

## Translating Common vegan Calls

Both vegan and this package commonly use community matrices with samples in
rows and species in columns. Important convention differences:

| vegan | DiversityAndDissimilarity.jl | Notes |
|---|---|---|
| `specnumber(x)` | `richness(x)` | Observed positive-abundance taxa. |
| `diversity(x, "shannon")` | `entropy(Shannon(; base=ℯ), x)` | vegan defaults to natural logs; this package defaults to bits. |
| `diversity(x, "simpson")` | `gini_simpson_index(x)` | vegan's Simpson output is `1 - sum(p_i^2)`. |
| `diversity(x, "invsimpson")` | `inverse_simpson_index(x)` | Inverse Simpson diversity. |
| `estimateR(x)` | `chao1(x)` / `ace(x)` | Richness estimators. |
| `vegdist(x, method="bray")` | `bray_curtis_distance(x)` | Bray-Curtis dissimilarity. |
| `vegdist(x, method="jaccard", binary=TRUE)` | `jaccard_distance(x)` | Presence/absence Jaccard. |
| `vegdist(decostand(x, "hellinger"), method="euclidean")` | `hellinger_distance(x)` | Hellinger distance. |

Broader ordination, constrained analysis, rarefaction curves, and ecological
modelling remain vegan strengths. This package focuses on lightweight,
convention-aware calculations.

## Generated Documentation Assets

The type tree is generated from exported package types:

```bash
julia --project=docs docs/make_type_trees.jl
```

This writes:

- `docs/src/assets/type-tree.dot`
- `docs/src/assets/type-tree.svg`
- `docs/src/assets/type-tree.pdf`

Graphviz's `dot` command is required for SVG and PDF rendering. The DOT file is
still written if Graphviz is unavailable.

## Local Documentation Build

Build the manual from the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'
```

The documentation uses doctests for many examples, so examples on these pages
are checked as part of the build.
