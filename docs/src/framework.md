# Framework

`DiversityAndDissimilarity.jl` exposes index metadata and workflow diagnostics in
addition to numeric index values. These helpers make conventions explicit and
make validation reproducible.

For the package's abstract and concrete index hierarchy, see
[Type Structure](type-structure.md).

## Index Metadata And Traits

Use [`index_metadata`](@ref) to inspect how an index is represented:

```jldoctest
julia> using DiversityAndDissimilarity

julia> m = index_metadata(BrayCurtis());

julia> m.family
:abundance

julia> m.output_mode
:dissimilarity

julia> m.is_metric
false

julia> m.is_semimetric
true

julia> m.bounds.lower_meaning
"minimal dissimilarity; identical or indistinguishable inputs"
```

The lower-level trait helpers are useful when writing generic workflows.
Properties that are not encoded for an index return `:unknown` where a Boolean
answer would overstate what the package knows.

```jldoctest
julia> using DiversityAndDissimilarity

julia> index_family(BrayCurtis())
:abundance

julia> input_mode(BrayCurtis())
:pairwise

julia> output_mode(BrayCurtis())
:dissimilarity

julia> is_metric(Jaccard())
true

julia> is_triangular(BrayCurtis())
false

julia> is_triangular(Overlap())
:unknown

julia> is_nonnegative(KullbackLeibler())
true

julia> is_bounded(JensenShannon())
true

julia> is_semimetric(BrayCurtis())
true

julia> is_similarity(Jaccard())
true

julia> is_dissimilarity(BrayCurtis())
true

julia> is_symmetric(KullbackLeibler())
false

julia> index_range(PielouEvenness())
(lower = 0.0, upper = 1.0)

julia> requires_probabilities(Hellinger())
true

julia> supports_matrix_kernel(Jaccard())
true
```

## Index Bounds

[`index_bounds`](@ref) adds semantic meaning to the numeric range. For
similarities, the lower bound conventionally means no overlap; for
dissimilarities and distances, it means identical inputs.

```jldoctest
julia> using DiversityAndDissimilarity

julia> b = index_bounds(Jaccard());

julia> b.lower, b.upper
(0.0, 1.0)

julia> b.lower_meaning
"minimal similarity; conventionally complete dissimilarity or no overlap"

julia> b.upper_meaning
"maximal similarity; conventionally identical or complete overlap"
```

```jldoctest
julia> using DiversityAndDissimilarity

julia> b = index_bounds(KullbackLeibler());

julia> b.upper
Inf

julia> b.upper_meaning
"unbounded dissimilarity; larger values mean greater separation"
```

## Metric Hierarchy

[`is_metric`](@ref) follows the usual
[metric space](https://en.wikipedia.org/wiki/Metric_space) convention:
nonnegative, zero only for identical inputs, symmetric, and satisfying the
[triangle inequality](https://en.wikipedia.org/wiki/Triangle_inequality).
The related helpers expose progressively weaker properties:

| Trait | Wikipedia | Requires | Notes |
|---|---|---|---|
| [`is_metric`](@ref) | [Metric space](https://en.wikipedia.org/wiki/Metric_space) | nonneg + d=0 iff same + symmetric + triangular | Full metric |
| [`is_pseudometric`](@ref) | [Pseudometric space](https://en.wikipedia.org/wiki/Pseudometric_space) | nonneg + symmetric + d=0 for same + triangular | Distinct inputs may have zero distance |
| [`is_quasimetric`](@ref) | [Quasimetric](https://en.wikipedia.org/wiki/Quasimetric_space) | nonneg + d=0 iff same + triangular | No symmetry required |
| [`is_metametric`](@ref) | [Generalizations](https://en.wikipedia.org/wiki/Metric_space#Generalizations) | nonneg + symmetric + d=0 for same | No triangle inequality; no identity of indiscernibles |
| [`is_semimetric`](@ref) | [Semimetric space](https://en.wikipedia.org/wiki/Semimetric_space) | nonneg + symmetric + d=0 iff same | No triangle inequality |
| [`is_premetric`](@ref) | [Premetric space](https://en.wikipedia.org/wiki/Premetric_space) | nonneg + d=0 for same | Most permissive |
| [`is_supermetric`](@ref) | [Ultrametric space](https://en.wikipedia.org/wiki/Ultrametric_space) | reverse-triangle condition | Uncommon; most indices return `false` or `:unknown` |

```jldoctest
julia> using DiversityAndDissimilarity

julia> is_metric(Jaccard())
true

julia> is_pseudometric(ShannonDifference())
true

julia> is_metametric(BrayCurtis())
true

julia> is_semimetric(BrayCurtis())
true

julia> is_premetric(KullbackLeibler())
true

julia> is_supermetric(Jaccard())
false
```

The full descriptor set with typical results:

```jldoctest
julia> using DiversityAndDissimilarity

julia> is_finite(KullbackLeibler())
false

julia> is_symmetric(KullbackLeibler())
false

julia> is_nonnegative(KullbackLeibler())
true

julia> is_bounded(KullbackLeibler())
false

julia> is_metric(KullbackLeibler())
false

julia> is_triangular(KullbackLeibler())
false

julia> is_pseudometric(KullbackLeibler())
false

julia> is_quasimetric(KullbackLeibler())
false

julia> is_metametric(KullbackLeibler())
false

julia> is_premetric(KullbackLeibler())
true

julia> is_supermetric(KullbackLeibler())
false

julia> is_similarity(KullbackLeibler())
false

julia> is_dissimilarity(KullbackLeibler())
true
```

## Reference Validation

[`reference_cases`](@ref) returns curated examples that pin package conventions
to formulas and external package behavior. [`validate_reference_cases`](@ref)
evaluates them:

```jldoctest
julia> using DiversityAndDissimilarity

julia> results = validate_reference_cases();

julia> all(r -> r.passed, results)
true

julia> results[1].name
"vegan_shannon_natural_log"

julia> results[1].observed ≈ 1.0397207708399179
true
```

## Estimator Reports

[`estimator_report`](@ref) compares Shannon estimators and reports basic
coverage diagnostics:

```jldoctest
julia> using DiversityAndDissimilarity

julia> r = estimator_report([1, 1, 2, 0, 5]);

julia> r.observed_richness
4

julia> r.singletons
2

julia> r.sample_coverage
0.7777777777777778

julia> r.estimates[1].name
:plugin
```

The report includes observed richness, singleton and doubleton counts, sample
coverage, estimator outputs, and warnings about common support/coverage issues.

## Diversity Audit

[`diversity_audit`](@ref) combines input validation, row-wise alpha summaries,
estimator diagnostics, and a labeled pairwise distance matrix:

```jldoctest
julia> using DiversityAndDissimilarity

julia> community = [1 1 2 0 5; 3 0 1 1 0];

julia> audit = diversity_audit(community; labels=["plot-a", "plot-b"]);

julia> audit.n_samples
2

julia> audit.n_taxa
5

julia> audit.pairwise.labels
2-element Vector{String}:
 "plot-a"
 "plot-b"

julia> audit.alpha[1].richness
4
```

For Tables.jl-compatible inputs, use `species` and `label` to separate species
columns from metadata:

```julia
diversity_audit(table; species=[:oak, :ash, :elm], label=:site)
```

## Uncertainty Audit

[`uncertainty_audit`](@ref) wraps the bootstrap tools into a compact report for
workflow checks. It reports Shannon entropy and effective-diversity intervals
for each sample together with labels, sample coverage, richness, and warnings:

```jldoctest
julia> using DiversityAndDissimilarity

julia> community = [1 1 2 0 5; 3 0 1 1 0];

julia> ua = uncertainty_audit(community; labels=["plot-a", "plot-b"], nboot=50);

julia> length(ua.reports)
2

julia> ua.reports[1].label
"plot-a"

julia> ua.reports[1].richness
4
```

Use a modest `nboot` while developing a workflow and a larger value for final
reports.

```@docs
index_metadata
index_family
input_mode
output_mode
is_finite
is_metric
is_triangular
is_nonnegative
is_bounded
is_pseudometric
is_quasimetric
is_metametric
is_semimetric
is_premetric
is_supermetric
is_similarity
is_dissimilarity
is_dissimiliarty
is_symmetric
index_range
index_bounds
requires_probabilities
supports_matrix_kernel
reference_cases
validate_reference_cases
estimator_report
diversity_audit
uncertainty_audit
```
