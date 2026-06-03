# Framework

`DiversityAndDissimilarity.jl` exposes index metadata and workflow diagnostics in
addition to numeric index values. These helpers make conventions explicit and
make validation reproducible.

## Index Metadata And Traits

Use [`index_metadata`](@ref) to inspect how an index is represented:

```julia
index_metadata(BrayCurtis())
index_metadata(GiniSimpson())
```

The lower-level trait helpers are useful when writing generic workflows:

```julia
index_family(BrayCurtis())
input_mode(BrayCurtis())
output_mode(BrayCurtis())
is_metric(Jaccard())
is_triangular(Jaccard())
is_nonnegative(KullbackLeibler())
is_bounded(JensenShannon())
is_semimetric(BrayCurtis())
is_similarity(Jaccard())
is_dissimilarity(BrayCurtis())
is_symmetric(KullbackLeibler())
index_range(PielouEvenness())
index_bounds(BrayCurtis())
requires_probabilities(Hellinger())
supports_matrix_kernel(Jaccard())
```

These traits are intentionally simple symbols and named tuples, so users can
branch on them without depending on internal implementation details.
Properties that are not encoded for an index return `:unknown` where a Boolean
answer would overstate what the package knows.

`is_metric` follows the usual metric convention: nonnegative, zero only for
identical inputs, symmetric, and triangular. The related helpers expose weaker
properties:

- `is_triangular` checks the triangle inequality only.
- `is_pseudometric` allows distinct inputs to have zero distance.
- `is_quasimetric` allows asymmetry while retaining the triangle inequality.
- `is_metametric` means nonnegative, symmetric, and zero for identical inputs,
  without requiring identity of indiscernibles or the triangle inequality.
- `is_semimetric` means nonnegative, symmetric, and zero only for identical
  inputs, without requiring the triangle inequality.
- `is_premetric` means nonnegative and zero for identical inputs.
- `is_supermetric` is reserved for reverse-triangle or supermetric-style
  properties; most implemented indices return `false` or `:unknown`.

Use [`index_bounds`](@ref) when the interpretation of a range matters. For a
similarity, the lower bound commonly means no overlap or complete dissimilarity;
for a dissimilarity or distance, the lower bound commonly means identical or
indistinguishable inputs.

The full descriptor set is:

```julia
is_finite(index)
is_symmetric(index)
is_nonnegative(index)
is_bounded(index)
is_metric(index)
is_triangular(index)
is_pseudometric(index)
is_quasimetric(index)
is_metametric(index)
is_semimetric(index)
is_premetric(index)
is_supermetric(index)
is_similarity(index)
is_dissimilarity(index)
```

`index_metadata(index)` includes these values, the basic
[`index_range`](@ref), and the interpreted [`index_bounds`](@ref) tuple.

## Reference Validation

[`reference_cases`](@ref) returns curated examples that pin package conventions
to formulas and external package behavior. [`validate_reference_cases`](@ref)
evaluates them:

```julia
validate_reference_cases()
```

This is a lightweight validation corpus rather than a replacement for the full
test suite. It is meant to make convention choices visible and reproducible.

## Estimator Reports

[`estimator_report`](@ref) compares Shannon estimators and reports basic
coverage diagnostics:

```julia
estimator_report([1, 1, 2, 0, 5]; support=6)
```

The report includes observed richness, singleton and doubleton counts, sample
coverage, estimator outputs, and warnings about common support/coverage issues.

## Diversity Audit

[`diversity_audit`](@ref) combines input validation, row-wise alpha summaries,
estimator diagnostics, and a labeled pairwise distance matrix:

```julia
audit = diversity_audit(
    community;
    labels=["plot-a", "plot-b"],
    pairwise_index=BrayCurtis(),
)

audit.alpha
audit.estimator_report
audit.pairwise
audit.warnings
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

```julia
uncertainty_audit(
    community;
    labels=["plot-a", "plot-b"],
    nboot=500,
)
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
