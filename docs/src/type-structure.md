# Type Structure

`DiversityAndDissimilarity.jl` is organized around small dispatchable types.
The main abstract roots are:

- [`DiversityIndex`](@ref), with [`AlphaDiversityIndex`](@ref) and
  [`PairwiseIndex`](@ref) beneath it.
- [`ShannonEstimator`](@ref), used by entropy and divergence estimators.

Concrete index values such as [`Shannon`](@ref), [`Richness`](@ref),
[`Jaccard`](@ref), and [`JensenShannon`](@ref) are lightweight objects passed to
generic operations such as [`entropy`](@ref), [`diversity`](@ref),
[`similarity`](@ref), and [`dissimilarity`](@ref). This keeps formulas,
metadata, validation, and matrix methods attached to the same dispatch surface.

## Type Tree

The diagram below is generated from the package's exported type definitions.
Abstract types are shown as boxes and concrete types as ellipses.

![DiversityAndDissimilarity type tree](assets/type-tree.svg)

Download the generated [DOT source](assets/type-tree.dot) or
[PDF version](assets/type-tree.pdf).

## Regenerating The Diagram

Run the generator from the repository root:

```bash
julia --project=docs docs/make_type_trees.jl
```

The script inspects the exported types defined by `DiversityAndDissimilarity`,
writes `docs/src/assets/type-tree.dot`, and then uses Graphviz `dot` to create
`docs/src/assets/type-tree.svg` and `docs/src/assets/type-tree.pdf`. If
Graphviz is not installed, the script still writes the DOT file and reports
that rendering was skipped.

## Design Notes

The type tree is intentionally shallow. Most behavior lives in methods and
metadata traits rather than deep inheritance. For example, [`index_metadata`](@ref)
and helpers such as [`is_metric`](@ref), [`is_similarity`](@ref), and
[`index_bounds`](@ref) describe conventions without forcing separate abstract
types for every possible property.
