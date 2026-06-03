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
index_range(PielouEvenness())
requires_probabilities(Hellinger())
supports_matrix_kernel(Jaccard())
```

These traits are intentionally simple symbols and named tuples, so users can
branch on them without depending on internal implementation details.

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
is_metric
index_range
requires_probabilities
supports_matrix_kernel
reference_cases
validate_reference_cases
estimator_report
diversity_audit
uncertainty_audit
```
